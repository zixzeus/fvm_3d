#pragma once

#include "riemann_solver.hpp"

namespace fvm3d::spatial {

/**
 * Lax-Friedrichs Riemann solver.
 * Simple and robust but more diffusive than HLL/HLLC.
 *
 * Flux: F_LF = 0.5 * (F_L + F_R) - 0.5 * lambda * (U_R - U_L)
 *
 * where lambda = max(|u| + a) is the maximum wave speed.
 *
 * Advantages: Simple, always convergent
 * Disadvantages: High numerical dissipation
 */
class LaxFriedrichsSolver : public RiemannSolver {
public:
    Eigen::VectorXd solve(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const override;

    double max_wave_speed(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const override;

    std::string name() const override { return "Lax-Friedrichs"; }

private:
    /**
     * Compute flux in a given direction.
     */
    Eigen::VectorXd compute_flux(
        const Eigen::VectorXd& U,
        int direction
    ) const;
};

} // namespace fvm3d::spatial
