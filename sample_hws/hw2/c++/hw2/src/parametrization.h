#pragma once

#include "tutte.h"
#include "utils.h"
#include "arap.h"
#include "symmetric_dirichlet.h"
#include "mips.h"
#include "eigen_system.h"
#include "line_search.h"

#include <omp.h>
#include <igl/doublearea.h>
#include <igl/adjacency_matrix.h>

inline auto projected_newton_precompute(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    bool unused;

    ParamData data;

    data.v = V;
    data.f = F;

    data.DmInv.resize(F.rows());
    data.PFPx.resize(F.rows());
    data.areas.resize(F.rows());

    for (int idx = 0; idx < F.rows(); ++idx)
    {
        const auto& f = F.row(idx);

        const auto& i = V.row(f[0]);
        const auto& j = V.row(f[1]);
        const auto& k = V.row(f[2]);

        const Eigen::Vector3d& eij = j - i;
        const Eigen::Vector3d& eik = k - i;

        Matrix2x2 x;

        x << Eigen::Vector2d(eij.norm(), 0),
            Eigen::Vector2d(eij.dot(eik), (eij.cross(eik)).norm()) / eij.norm();

        x.computeInverseWithCheck(data.DmInv[idx], unused);

        data.PFPx[idx] = compute_PFPx(data.DmInv[idx]);
    }

    igl::doublearea(V, F, data.areas);
    data.areas /= 2;

    data.mesh_area = data.areas.sum();

    Eigen::SparseMatrix<int> A;
    igl::adjacency_matrix(F, A);

    Eigen::SparseVector<int> Asum;
    igl::sum(A, 1, Asum);

    data.max_reserve = Asum.coeffs().maxCoeff();

    return data;
}

template <typename T>
auto deformation_gradient_and_hessian(const Eigen::Vector2d& i, const Eigen::Vector2d& j, const Eigen::Vector2d& k,
                                      const ParamData& data, const int idx)
{
    Matrix2x2 Ds;

    Ds << j - i, k - i;

    const auto& F = Ds * data.DmInv[idx];
    const auto& PFPx = data.PFPx[idx];

    const auto& eig = T::eval_eigen_system(F);

    Eigen::Matrix4d P2PhiP2f = 0.1 * Eigen::Matrix4d::Identity();

    for (int ei = 0; ei < 4; ++ei)
    {
        P2PhiP2f += std::max(eig.lambda[ei], 0.) * eig.eigenvector[ei] * eig.eigenvector[ei].transpose();
    }

    return std::tuple<Vector6, Matrix6x6>
    {
        data.areas[idx] * PFPx.transpose() * T::compute_PPhiPf(F, eig),
        data.areas[idx] * PFPx.transpose() * P2PhiP2f * PFPx
    };
}

template <typename T>
auto energy(const Eigen::MatrixXd& x, const ParamData& data)
{
    double energy = 0.0;

    for (int idx = 0; idx < data.f.rows(); ++idx)
    {
        const auto& f = data.f.row(idx);
        const auto i = f[0], j = f[1], k = f[2];

        const auto& vi = x.middleRows<2>(2 * i);
        const auto& vj = x.middleRows<2>(2 * j);
        const auto& vk = x.middleRows<2>(2 * k);

        Matrix2x2 Ds;
        Ds << vj - vi, vk - vi;
        const auto& F = Ds * data.DmInv[idx];

        energy += data.areas[idx] * T::energy(F);
    }
    return energy;
}

template <typename T>
auto project_gradient_and_hessian(const Eigen::MatrixXd& x, const ParamData& data)
{
    std::vector<triplet> triplet_list(36 * data.f.rows());
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(x.size());
    std::vector<std::vector<double>> gradients(x.size());
    Eigen::SparseMatrix<double> hessian(x.size(), x.size());

#pragma omp parallel for num_threads(4) schedule(static)
    for (int i = 0; i < gradients.size(); ++i)
    {
        gradients[i].reserve(data.max_reserve);
    }

#pragma omp parallel for num_threads(4) schedule(dynamic)
    for (int idx = 0; idx < data.f.rows(); ++idx)
    {
        const auto& f = data.f.row(idx);
        const auto i = f[0], j = f[1], k = f[2];

        const auto& globals = std::array{2 * i, 2 * i + 1, 2 * j, 2 * j + 1, 2 * k, 2 * k + 1};

        const auto& vi = x.middleRows<2>(2 * i);
        const auto& vj = x.middleRows<2>(2 * j);
        const auto& vk = x.middleRows<2>(2 * k);

        const auto& [g, h] = deformation_gradient_and_hessian<T>(vi, vj, vk, data, idx);

        for (int m = 0; m < 6; ++m)
        {
            for (int n = 0; n < 6; ++n)
            {
                triplet_list[idx * 36 + m * 6 + n] = {globals[m], globals[n], h(m, n)};
            }

#pragma omp critical
            {
                gradients[globals[m]].push_back(g[m]);
            }
        }
    }

    kahan_accumulation init = { 0, 0 };
#pragma omp parallel for num_threads(4) schedule(static)
    for (int i = 0; i < gradients.size(); ++i)
    {
        gradient[i] = std::accumulate(gradients[i].begin(), gradients[i].end(), init, kahan_sum).sum;
    }

    hessian.setFromTriplets(triplet_list.begin(), triplet_list.end());

    return std::tuple{gradient, hessian};
}

template <typename T>
auto projected_newton_solver(
    Eigen::MatrixXd& UV,
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const ParamData& data, const int max_iter = 3000,
    const double esp = 1e-4)
{
    Eigen::MatrixXd x = UV.reshaped<Eigen::RowMajor>();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    bool is_analyzed = false;

    std::vector<double> energy_vector(max_iter + 1, -1);

    auto cur_energy = energy<T>(x, data);
    energy_vector[0] = cur_energy / data.mesh_area;

    for (int i = 0; i < max_iter; ++i)
    {
        const auto& [bi, Hi] = project_gradient_and_hessian<T>(x, data);

        if (bi.template lpNorm<Eigen::Infinity>() < esp)
        {
            break;
        }

        if (!is_analyzed)
        {
            solver.analyzePattern(Hi);
            is_analyzed = true;
        }

        solver.factorize(Hi);

        if (solver.info() != Eigen::Success)
        {
            break;
        }

        const Eigen::VectorXd& di = solver.solve(-bi);

        if (solver.info() != Eigen::Success)
        {
            break;
        }

        cur_energy = flip_avoiding_line_search(x, di, bi, data, energy<T>, cur_energy);
        energy_vector[i + 1] = cur_energy / data.mesh_area;
    }

    UV = x.reshaped<Eigen::RowMajor>(V.rows(), 2);
    return energy_vector;
}

inline auto projected_newton_solver(
    Eigen::MatrixXd& UV,
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const ParamData& data, const int type, const int max_iter = 200,
    const double esp = 1e-4)
{
    if (static_cast<energy_type>(type) == energy_type::arap)
        return projected_newton_solver<arap>(
            UV, V, F, data, max_iter, esp);
    if (static_cast<energy_type>(type) == energy_type::sd)
        return projected_newton_solver<symmetric_dirichlet>(
            UV, V, F, data, max_iter, esp);
    if (static_cast<energy_type>(type) == energy_type::mips)
        return projected_newton_solver<mips>(
            UV, V, F, data, max_iter, esp);
    return projected_newton_solver<symmetric_dirichlet>(UV, V, F, data, max_iter, esp);
}

template <typename T>
auto projected_newton_solver_animation(
    Eigen::MatrixXd& UV,
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const ParamData& data,
    const double esp = 1e-4)
{
    Eigen::MatrixXd x = UV.reshaped<Eigen::RowMajor>();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    const auto old_energy = energy<T>(x, data);
    double new_energy = -1;

    const auto& [bi, Hi] = project_gradient_and_hessian<T>(x, data);

    if (bi.template lpNorm<Eigen::Infinity>() < esp)
    {
        return std::tuple{old_energy, new_energy};
    }

    solver.compute(Hi);

    if (solver.info() != Eigen::Success)
    {
        return std::tuple{old_energy, new_energy};
    }

    const Eigen::VectorXd& di = solver.solve(-bi);

    if (solver.info() != Eigen::Success)
    {
        return std::tuple{old_energy, new_energy};
    }

    new_energy = flip_avoiding_line_search(x, di, bi, data, energy<T>, old_energy);

    UV = x.reshaped<Eigen::RowMajor>(V.rows(), 2);
    return std::tuple{old_energy, new_energy};
}

inline auto projected_newton_solver_animation(
    Eigen::MatrixXd& UV,
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const ParamData& data, const int type,
    const double esp = 1e-4)
{
    if (static_cast<energy_type>(type) == energy_type::arap)
        return projected_newton_solver_animation<arap>(
            UV, V, F, data, esp);
    if (static_cast<energy_type>(type) == energy_type::sd)
        return projected_newton_solver_animation<symmetric_dirichlet>(
            UV, V, F, data, esp);
    if (static_cast<energy_type>(type) == energy_type::mips)
        return projected_newton_solver_animation<mips>(
            UV, V, F, data, esp);
    return projected_newton_solver_animation<symmetric_dirichlet>(UV, V, F, data, esp);
}
