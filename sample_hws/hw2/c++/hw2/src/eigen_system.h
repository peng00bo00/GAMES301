#pragma once

#include "utils.h"

namespace EigenSystem
{
    inline static Matrix2x2 twist = (Matrix2x2() << 0, -1, 1, 0).finished();
    inline static Matrix2x2 flip = (Matrix2x2() << 0, 1, 1, 0).finished();
    inline static Matrix2x2 pinch = (Matrix2x2() << 1, 0, 0, -1).finished();

    struct eigen_system
    {
        std::array<double, 4> lambda;
        std::array<Eigen::Vector4d, 4> eigenvector;
        std::array<double, 4> I;
        Matrix2x2 R, S;
    };
}
