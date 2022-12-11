#pragma once

#include "utils.h"

struct mips
{
    static EigenSystem::eigen_system eval_eigen_system(const Matrix2x2&);
    static double energy(const Matrix2x2&);
    static Vector4 compute_PPhiPf(const Matrix2x2&, const EigenSystem::eigen_system&);
};
