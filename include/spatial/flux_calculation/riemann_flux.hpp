#pragma once

#include "spatial/flux_calculation/flux_calculator_base.hpp"
#include "spatial/riemann_solvers/riemann_solver.hpp"
#include "physics/physics_base.hpp"
#include <memory>
#include <string>

namespace fvm3d::spatial {

/**
 * Generic Riemann solver-based numerical flux calculator.
 *
 * Uses approximate Riemann solvers (HLL, HLLC, HLLD, etc.) to compute
 * numerical fluxes from left and right interface states.
 *
 * This class wraps different Riemann solver implementations and provides
 * a unified FluxCalculator interface.
 */
class RiemannFlux : public FluxCalculator {
public:
    /**
     * Constructor.
     *
     * @param solver_type: Type of Riemann solver ("hll", "hllc", "hlld")
     * @param physics: Physics equation system (Euler, MHD, etc.)
     */
    RiemannFlux(const std::string& solver_type, std::shared_ptr<physics::PhysicsBase> physics);

    /**
     * Compute numerical flux using Riemann solver.
     *
     * @param U_L: Left interface state
     * @param U_R: Right interface state
     * @param physics: Physics equation object (ignored, uses internal physics)
     * @param direction: 0=x, 1=y, 2=z
     * @return: Numerical flux vector
     */
    Eigen::VectorXd compute_flux(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        const physics::PhysicsBase& physics,
        int direction
    ) const override;

    /**
     * Compute maximum wave speed for CFL condition.
     *
     * @param U_L: Left state
     * @param U_R: Right state
     * @param physics: Physics equation object (ignored, uses internal physics)
     * @param direction: Direction (0=x, 1=y, 2=z)
     * @return: Maximum wave speed
     */
    double compute_max_wave_speed(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        const physics::PhysicsBase& physics,
        int direction
    ) const override;

    /**
     * Check if this flux calculator is physics-agnostic.
     */
    bool is_physics_agnostic() const override;

    /**
     * Check if this flux calculator supports the given physics type.
     */
    bool supports_physics(const std::string& physics_type) const override;

    /**
     * Riemann solvers compute wave speeds internally.
     */
    bool needs_wave_speed() const override { return true; }

    /**
     * Get the type of Riemann solver.
     */
    std::string get_solver_type() const { return solver_type_; }

private:
    std::unique_ptr<RiemannSolver> riemann_solver_;  ///< Underlying Riemann solver
    std::shared_ptr<physics::PhysicsBase> physics_;  ///< Physics equation system
    std::string solver_type_;                        ///< Type of solver (hll, hllc, hlld)
};

} // namespace fvm3d::spatial
