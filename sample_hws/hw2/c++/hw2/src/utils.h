#pragma once

#include <Eigen/Core>
#include <igl/polar_svd.h>
#include <omp.h>

using Matrix2x2 = Eigen::Matrix2d;
using Matrix4x6 = Eigen::Matrix<double, 4, 6>;

using Vector4 = Eigen::Matrix<double, 4, 1>;
using Vector6 = Eigen::Matrix<double, 6, 1>;
using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

struct ParamData
{
    Eigen::MatrixXd v;
    Eigen::MatrixXi f;

    std::vector<Matrix2x2> DmInv;
    std::vector<Matrix4x6> PFPx;

    Eigen::VectorXd areas;

    double mesh_area;
    int max_reserve;
};

namespace EigenSystem
{
    struct eigen_system;
}

enum class energy_type
{
    sd,
    arap,
    mips
};

struct kahan_accumulation
{
    double sum;
    double correction;
};

inline auto kahan_sum(const kahan_accumulation& accumulation, const double value)
{
    const auto y = value - accumulation.correction;
    const auto t = accumulation.sum + y;

    return kahan_accumulation{t, (t - accumulation.sum) - y };
}

template <typename Scalar, int Rows, int Cols>
auto vec(const Eigen::Matrix<Scalar, Rows, Cols>& f)
{
    return f.reshaped();
}

inline auto vec(const Eigen::MatrixXd& f)
{
    return f.reshaped();
}

inline auto compute_PFPx(const Matrix2x2& DmInv)
{
    const auto m = DmInv(0, 0);
    const auto n = DmInv(0, 1);
    const auto p = DmInv(1, 0);
    const auto q = DmInv(1, 1);

    const auto t1 = -m - p;
    const auto t2 = -n - q;

    Matrix4x6 PFPx = Matrix4x6::Zero();
    PFPx(0, 0) = t1;
    PFPx(0, 2) = m;
    PFPx(0, 4) = p;
    PFPx(1, 1) = t1;
    PFPx(1, 3) = m;
    PFPx(1, 5) = p;
    PFPx(2, 0) = t2;
    PFPx(2, 2) = n;
    PFPx(2, 4) = q;
    PFPx(3, 1) = t2;
    PFPx(3, 3) = n;
    PFPx(3, 5) = q;

    return PFPx;
}

namespace flip_avoiding
{
    inline auto get_min_pos_root(const Eigen::MatrixXd& x, const Eigen::VectorXd& d, const ParamData& data,
                                 const int idx)
    {
        const auto& f = data.f.row(idx);
        const auto i = f[0], j = f[1], k = f[2];

        const auto& vi = x.middleRows<2>(2 * i);
        const auto& vj = x.middleRows<2>(2 * j);
        const auto& vk = x.middleRows<2>(2 * k);

        const auto& di = d.segment<2>(2 * i);
        const auto& dj = d.segment<2>(2 * j);
        const auto& dk = d.segment<2>(2 * k);

        Matrix2x2 Ds, D;
        Ds << vj - vi, vk - vi;
        D << dj - di, dk - di;
        const auto& step_matrix = D * Ds.inverse();

        Eigen::Vector2d Sigma;
        Matrix2x2 U, V, R, S;

        igl::polar_svd(step_matrix, R, S, U, Sigma, V);

        const auto t1 = std::max(Sigma[0], Sigma[1]);

        if (t1 <= 0)
        {
            return std::numeric_limits<double>::max();
        }
        return 1 / t1;
    }

    inline double compute_max_step_from_singularities(const Eigen::MatrixXd& x, const Eigen::VectorXd& d,
                                                      const ParamData& data)
    {
        double max_step = std::numeric_limits<double>::max();
        std::vector<double> min_roots(data.f.rows());

#pragma omp parallel for num_threads(4) schedule(static)
        for (int idx = 0; idx < data.f.rows(); ++idx)
        {
            min_roots[idx] = get_min_pos_root(x, d, data, idx);
        }

        for (const auto& root: min_roots)
        {
            max_step = std::min(max_step, root);
        }

        return max_step;
    }
}
