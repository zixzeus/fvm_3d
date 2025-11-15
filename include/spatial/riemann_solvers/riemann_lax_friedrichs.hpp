#pragma once

#include "riemann_solver.hpp"
#include "physics/physics_base.hpp"
#include <memory>

namespace fvm3d::spatial {

/**
 * Lax-Friedrichs Riemann solver.
 *
 * Simple and robust but more diffusive than HLL/HLLC.
 *
 * Flux: F_LF = 0.5 * (F_L + F_R) - 0.5 * lambda * (U_R - U_L)
 *
 * where lambda = max(|u| + a) is the maximum wave speed.
 */
class LaxFriedrichsSolver : public RiemannSolver {
public:
    /**
     * Constructor with physics object.
     * @param physics: Physics equation system (Euler, MHD, etc.)
     */
    explicit LaxFriedrichsSolver(const std::shared_ptr<physics::PhysicsBase>& physics);

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
    std::shared_ptr<physics::PhysicsBase> physics_;
};

} // namespace fvm3d::spatial
