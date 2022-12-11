#include "arap.h"
#include "eigen_system.h"

#include <igl/polar_svd.h>

EigenSystem::eigen_system arap::eval_eigen_system(const Matrix2x2& F)
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
        2 - 4 / (Sigma[0] + Sigma[1]),
        2,
        2,
        2
    };

    eigens.eigenvector = {
        vec(1 / std::sqrt(2) * U * EigenSystem::twist * V.transpose()),
        vec(1 / std::sqrt(2) * U * EigenSystem::flip * V.transpose()),
        vec(U * Eigen::Vector2d(1, 0).asDiagonal() * V.transpose()),
        vec(U * Eigen::Vector2d(0, 1).asDiagonal() * V.transpose())
    };

    return eigens;
}

double arap::energy(const Matrix2x2& F)
{
    const auto& eig = eval_eigen_system(F);

    return eig.I[2] - 2 * eig.I[1] + 2;
}

Vector4 arap::compute_PPhiPf(const Matrix2x2& F, const EigenSystem::eigen_system& eig)
{
    const auto& f_vec = vec(F);

    const auto& r_vec = vec(eig.R);

    return f_vec - r_vec;
}
