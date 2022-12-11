#include "mips.h"
#include "eigen_system.h"

#include <igl/polar_svd.h>

EigenSystem::eigen_system mips::eval_eigen_system(const Matrix2x2& F)
{
    Eigen::Vector2d Sigma;
    Matrix2x2 U, V, R, S;

    igl::polar_svd(F, R, S, U, Sigma, V);

    EigenSystem::eigen_system eigens;

    eigens.R = R;
    eigens.S = S;

    eigens.I[1] = S.trace();
    eigens.I[2] = S.squaredNorm();
    eigens.I[3] = S.determinant();

    const auto alpha = std::sqrt(std::pow(eigens.I[2], 2) - 3 * std::pow(eigens.I[3], 3));
    const auto beta = eigens.I[2] / (Sigma[1] * Sigma[1] - Sigma[0] * Sigma[0] + alpha);
    const auto gamma = std::sqrt(1 + beta * beta);

    eigens.lambda = {
        2 / eigens.I[3] - eigens.I[2] / std::pow(eigens.I[3], 2),
        2 / eigens.I[3] + eigens.I[2] / std::pow(eigens.I[3], 2),
        eigens.I[2] * (eigens.I[2] - alpha) - 2 * std::pow(eigens.I[3], 2),
        eigens.I[2] * (eigens.I[2] + alpha) - 2 * std::pow(eigens.I[3], 2)
    };

    const auto& t = vec(1 / std::sqrt(2) * U * EigenSystem::twist * V.transpose());
    const auto& l = vec(1 / std::sqrt(2) * U * EigenSystem::flip * V.transpose());
    const auto& d1 = vec(U * Eigen::Vector2d(1, 0).asDiagonal() * V.transpose());
    const auto& d2 = vec(U * Eigen::Vector2d(0, 1).asDiagonal() * V.transpose());

    eigens.eigenvector = {
        t, l, (beta * d1 + d2) / gamma, (d1 + beta * d2) / gamma
    };

    return eigens;
}

double mips::energy(const Matrix2x2& F)
{
    const auto& eig = eval_eigen_system(F);

    return eig.I[2] / eig.I[3];
}

Vector4 mips::compute_PPhiPf(const Matrix2x2& F, const EigenSystem::eigen_system& eig)
{
    const auto& f_vec = vec(F);

    const auto& G = EigenSystem::twist * F * EigenSystem::twist.transpose();

    const auto& g_vec = vec(G);

    return f_vec / eig.I[3] - eig.I[2] / std::pow(eig.I[3], 2) * g_vec;
}
