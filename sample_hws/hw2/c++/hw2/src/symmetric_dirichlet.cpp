#include "symmetric_dirichlet.h"
#include "eigen_system.h"

#include <igl/polar_svd.h>

EigenSystem::eigen_system symmetric_dirichlet::eval_eigen_system(const Matrix2x2& F)
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

    eigens.lambda = {
        1 + 3 / std::pow(Sigma[0], 4),
        1 + 3 / std::pow(Sigma[1], 4),
        1 + 1 / std::pow(eigens.I[3], 2) - eigens.I[2] / std::pow(eigens.I[3], 3),
        1 + 1 / std::pow(eigens.I[3], 2) + eigens.I[2] / std::pow(eigens.I[3], 3)
    };

    eigens.eigenvector = {
        vec(U * Eigen::Vector2d(1, 0).asDiagonal() * V.transpose()),
        vec(U * Eigen::Vector2d(0, 1).asDiagonal() * V.transpose()),
        vec(1 / std::sqrt(2) * U * EigenSystem::flip * V.transpose()),
        vec(1 / std::sqrt(2) * U * EigenSystem::twist * V.transpose())
    };

    return eigens;
}

double symmetric_dirichlet::energy(const Matrix2x2& F)
{
    return F.squaredNorm() + F.inverse().squaredNorm();
}

Vector4 symmetric_dirichlet::compute_PPhiPf(const Matrix2x2& F, const EigenSystem::eigen_system& eig)
{
    const auto& f_vec = vec(F);

    const auto& G = EigenSystem::twist * F * EigenSystem::twist.transpose();

    const auto& g_vec = vec(G);

    return (1 + 1 / (eig.I[3] * eig.I[3])) * f_vec - eig.I[2] / (eig.I[3] * eig.I[3] * eig.I[3]) * g_vec;
}
