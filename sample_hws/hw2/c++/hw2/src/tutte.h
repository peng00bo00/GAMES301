#pragma once

#include <deque>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/flipped_triangles.h>

using triplet = Eigen::Triplet<double>;

enum class weight_method
{
    uniform,
    floater,
    harmonic
};

enum class border_type
{
    circle,
    square
};

inline auto map_vertices_to_square(const Eigen::MatrixXd& V, const Eigen::VectorXi& bnd, Eigen::MatrixXd& UV)
{
    // Get sorted list of boundary vertices
    std::vector<int> interior, map_ij;
    map_ij.resize(V.rows());

    std::deque is_on_bnd(V.rows(), false);
    for (int i = 0; i < bnd.size(); i++)
    {
        is_on_bnd[bnd[i]] = true;
        map_ij[bnd[i]] = i;
    }

    for (int i = 0; i < is_on_bnd.size(); i++)
    {
        if (!is_on_bnd[i])
        {
            map_ij[i] = interior.size();
            interior.push_back(i);
        }
    }

    // Map boundary to unit square
    std::vector<double> len(bnd.size());
    len[0] = 0.;

    for (int i = 1; i < bnd.size(); i++)
    {
        len[i] = len[i - 1] + (V.row(bnd[i - 1]) - V.row(bnd[i])).norm();
    }
    const double first_len = len[len.size() / 4];
    const double second_len = len[len.size() / 2];
    const double third_len = len[3 * len.size() / 4];
    const double total_len = len[len.size() - 1] + (V.row(bnd[0]) - V.row(bnd[bnd.size() - 1])).norm();

    UV.resize(bnd.size(), 2);
    for (int i = 0; i < bnd.size(); i++)
    {
        if (i < bnd.size() / 4)
        {
            const double frac = len[i] / first_len;
            UV.row(map_ij[bnd[i]]) << 1, 2 * frac - 1;
        }
        else if (i < bnd.size() / 2)
        {
            const double frac = (len[i] - first_len) / (second_len - first_len);
            UV.row(map_ij[bnd[i]]) << 1 - 2 * frac, 1;
        }
        else if (i < 3 * bnd.size() / 4)
        {
            const double frac = (len[i] - second_len) / (third_len - second_len);
            UV.row(map_ij[bnd[i]]) << -1, 1 - 2 * frac;
        }
        else
        {
            const double frac = (len[i] - third_len) / (total_len - third_len);
            UV.row(map_ij[bnd[i]]) << 2 * frac - 1, -1;
        }
    }
}

template <weight_method WeightMethod, border_type BorderType>
auto tutte(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV)
{
    UV = Eigen::MatrixXd::Zero(V.rows(), 2);
    const auto& ones = Eigen::MatrixXd::Ones(V.rows(), 2);
    Eigen::VectorXi bnd;
    igl::boundary_loop(F, bnd);

    if (bnd.size() == 0)
    {
        std::cerr << "ERROR: no boundary found!\n";
        return false;
    }

    Eigen::MatrixXd bnd_uv;

    if constexpr (BorderType == border_type::circle)
    {
        igl::map_vertices_to_circle(V, bnd, bnd_uv);
    }
    else if constexpr (BorderType == border_type::square)
    {
        if (bnd.size() < 4)
        {
            std::cerr << "ERROR: square border parametrization requires a longer border\n";
            return false;
        }
        map_vertices_to_square(V, bnd, bnd_uv);
    }

    if constexpr (WeightMethod == weight_method::harmonic)
    {
        igl::harmonic(V, F, bnd, bnd_uv, 1, UV);
        if (igl::flipped_triangles(UV, F).size() != 0)
        {
            igl::harmonic(F, bnd, bnd_uv, 1, UV);
        }
        return true;
    }

    std::deque is_boundary(V.rows(), false);
    for (int i = 0; i < bnd.rows(); i++)
    {
        is_boundary[bnd[i]] = true;
    }

    std::vector<triplet> triple_list;

    for (int f = 0; f < F.rows(); f++)
    {
        for (int i = 0; i < 3; i++)
        {
            int idx = F(f, i);

            if (is_boundary[idx])
            {
                continue;
            }

            const auto j = (i + 1) % 3;
            const auto k = (i + 2) % 3;

            double w_ij, w_ik;

            if constexpr (WeightMethod == weight_method::uniform)
            {
                w_ij = w_ik = 1;
            }
            else if constexpr (WeightMethod == weight_method::floater)
            {
                const auto& ij = V.row(F(f, j)) - V.row(idx);
                const auto& ik = V.row(F(f, k)) - V.row(idx);

                const auto theta = std::acos(ij.dot(ik) / (ij.norm() * ik.norm()));
                const auto tan = std::tan(theta / 2);
                w_ij = tan / ij.norm();
                w_ik = tan / ik.norm();
            }

            triple_list.emplace_back(idx, F(f, j), w_ij);
            triple_list.emplace_back(idx, F(f, k), w_ik);

            triple_list.emplace_back(idx, idx, -w_ij);
            triple_list.emplace_back(idx, idx, -w_ik);
        }
    }

    for (int i = 0; i < bnd.rows(); i++)
    {
        triple_list.emplace_back(bnd[i], bnd[i], 1);
    }

    Eigen::SparseMatrix<double> L(V.rows(), V.rows());
    L.setFromTriplets(triple_list.begin(), triple_list.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;

    Eigen::VectorXd b_u = Eigen::VectorXd::Zero(V.rows());
    Eigen::VectorXd b_v = Eigen::VectorXd::Zero(V.rows());

    for (int i = 0; i < bnd.rows(); i++)
    {
        b_u[bnd[i]] = bnd_uv(i, 0);
        b_v[bnd[i]] = bnd_uv(i, 1);
    }

    solver.compute(L);

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "ERROR: decomposition failed\n";
        return false;
    }

    const auto& u = solver.solve(b_u);

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "ERROR: solving u failed\n";
        return false;
    }

    const auto& v = solver.solve(b_v);

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "ERROR: solving v failed\n";
        return false;
    }

    UV << u, v;
    UV += ones;
    UV /= 2;
    return true;
}

// Auto Generation by Enum Reflection.
inline auto tutte(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV, int method, int border)
{
    if (static_cast<weight_method>(method) == weight_method::uniform && static_cast<border_type>(border) ==
        border_type::circle)
        return tutte<weight_method::uniform, border_type::circle>(V, F, UV);
    if (static_cast<weight_method>(method) == weight_method::uniform && static_cast<border_type>(border) ==
        border_type::square)
        return tutte<weight_method::uniform, border_type::square>(V, F, UV);
    if (static_cast<weight_method>(method) == weight_method::floater && static_cast<border_type>(border) ==
        border_type::circle)
        return tutte<weight_method::floater, border_type::circle>(V, F, UV);
    if (static_cast<weight_method>(method) == weight_method::floater && static_cast<border_type>(border) ==
        border_type::square)
        return tutte<weight_method::floater, border_type::square>(V, F, UV);
    if (static_cast<weight_method>(method) == weight_method::harmonic && static_cast<border_type>(border) ==
        border_type::circle)
        return tutte<weight_method::harmonic, border_type::circle>(V, F, UV);
    if (static_cast<weight_method>(method) == weight_method::harmonic && static_cast<border_type>(border) ==
        border_type::square)
        return tutte<weight_method::harmonic, border_type::square>(V, F, UV);
    return tutte<weight_method::harmonic, border_type::circle>(V, F, UV);
}
