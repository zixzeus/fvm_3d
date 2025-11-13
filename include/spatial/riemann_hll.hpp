#pragma once

#include "riemann_solver.hpp"

namespace fvm3d::spatial {

/**
 * HLL (Harten-Lax-van Leer) Riemann solver.
 * More accurate than Lax-Friedrichs, especially for shear flows.
 */
class HLLSolver : public RiemannSolver {
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

    std::string name() const override { return "HLL"; }

private:
    double estimate_s_left(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const;

    double estimate_s_right(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const;

    Eigen::VectorXd compute_flux(
        const Eigen::VectorXd& U,
        int direction
    ) const;
};

} // namespace fvm3d::spatial
