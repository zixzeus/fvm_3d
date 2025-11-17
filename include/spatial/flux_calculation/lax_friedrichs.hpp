#pragma once

#include "spatial/flux_calculation/flux_calculator_base.hpp"
#include "physics/physics_base.hpp"
#include <memory>

namespace fvm3d::spatial {

/**
 * Lax-Friedrichs flux calculator.
 * Simple and robust but more diffusive than HLL/HLLC.
 *
 * Flux: F_LF = 0.5 * (F_L + F_R) - 0.5 * lambda * (U_R - U_L)
 *
 * where lambda = max(|u| + a) is the maximum wave speed.
 *
 * Uses PhysicsBase interface for generality - works with any physics
 * (Euler, MHD, etc.) by delegating flux calculation to the physics object.
 *
 * Advantages: Simple, always convergent, physics-agnostic
 * Disadvantages: High numerical dissipation
 */
class LaxFriedrichsFlux : public FluxCalculator {
public:
    /**
     * Constructor with physics object.
     * @param physics: Physics equation system (Euler, MHD, etc.)
     */
    explicit LaxFriedrichsFlux(std::shared_ptr<physics::PhysicsBase> physics);

    Eigen::VectorXd compute_flux(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        const physics::PhysicsBase& physics,
        int direction
    ) const override;

    double compute_max_wave_speed(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        const physics::PhysicsBase& physics,
        int direction
    ) const override;

    bool is_physics_agnostic() const override { return true; }
    bool needs_wave_speed() const override { return true; }

private:
    std::shared_ptr<physics::PhysicsBase> physics_;  // Internal physics for compatibility
};

} // namespace fvm3d::spatial
