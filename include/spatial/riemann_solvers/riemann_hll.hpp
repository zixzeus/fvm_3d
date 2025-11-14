#pragma once

#include "riemann_solver.hpp"
#include "physics/physics_base.hpp"
#include <memory>

namespace fvm3d::spatial {

/**
 * HLL (Harten-Lax-van Leer) Riemann solver.
 * More accurate than Lax-Friedrichs, especially for shear flows.
 *
 * Now uses PhysicsBase interface for generality - works with any physics
 * (Euler, MHD, etc.) by delegating flux calculation to the physics object.
 */
class HLLSolver : public RiemannSolver {
public:
    /**
     * Constructor with physics object.
     * @param physics: Physics equation system (Euler, MHD, etc.)
     */
    explicit HLLSolver(const std::shared_ptr<physics::PhysicsBase>& physics);

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
    std::shared_ptr<physics::PhysicsBase> physics_;  // Physics object for flux calculation

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
};

} // namespace fvm3d::spatial
