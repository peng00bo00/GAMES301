#pragma once

#include "utils.h"

template <typename EnergyFun>
auto line_search(Eigen::MatrixXd& x, const Eigen::VectorXd& d, const Eigen::VectorXd& g, const ParamData& data,
                 EnergyFun fun, const double cur_energy = -1, const double max_step_size = 1.0,
                 const double shrink = 0.8, const int max_iters = 64,
                 const double armijo_esp = 1e-4)
{
    double old_energy;
    if (cur_energy > 0)
    {
        old_energy = cur_energy;
    }
    else
    {
        old_energy = fun(x, data);
    }

    auto step_size = max_step_size;

    for (int i = 0; i < max_iters; ++i)
    {
        Eigen::MatrixXd new_x = x + step_size * d;

        const auto new_energy = fun(new_x, data);

        if (new_energy <= old_energy + armijo_esp * step_size * d.dot(g))
        {
            x = new_x;
            return new_energy;
        }

        step_size *= shrink;
    }
    return old_energy;
}

template <typename EnergyFun>
auto flip_avoiding_line_search(Eigen::MatrixXd& x, const Eigen::VectorXd& d, const Eigen::VectorXd& g,
                               const ParamData& data, EnergyFun fun, const double cur_energy = -1, double shrink = 0.8,
                               int max_iters = 64,
                               double armijo_esp = 1e-4)
{
    const auto min_step_to_singularity = flip_avoiding::compute_max_step_from_singularities(x, d, data);
    const auto max_step = std::min(1.0, 0.99 * min_step_to_singularity);

    return line_search(x, d, g, data, fun, cur_energy, max_step, shrink, max_iters, armijo_esp);
}
