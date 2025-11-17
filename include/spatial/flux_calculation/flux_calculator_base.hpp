#pragma once

#include "physics/physics_base.hpp"
#include <Eigen/Dense>
#include <memory>
#include <string>

namespace fvm3d::spatial {

/**
 * Abstract base class for all numerical flux calculation methods.
 *
 * Flux calculators compute numerical fluxes from left and right interface states
 * that have been reconstructed by spatial reconstruction methods.
 *
 * Design Philosophy:
 * - Unified interface for all flux calculation methods (Riemann solvers, central schemes, etc.)
 * - Physics-agnostic: works with any PhysicsBase implementation
 * - Supports both single-interface and vectorized computations
 *
 * Concrete implementations include:
 * - Riemann solver-based: Lax-Friedrichs, HLL, HLLC, HLLD
 * - Central schemes: Central differencing, spectral methods
 * - Hybrid methods: AUSM, Roe-type schemes
 */
class FluxCalculator {
public:
    /**
     * Constructor.
     * @param name: Name of the flux calculator (e.g., "HLL", "HLLC")
     * @param flux_type: Type category ("riemann", "central", "hybrid")
     */
    explicit FluxCalculator(const std::string& name, const std::string& flux_type = "generic")
        : name_(name), flux_type_(flux_type) {}

    virtual ~FluxCalculator() = default;

    /**
     * Compute numerical flux from left and right interface states.
     *
     * This is the core method that all flux calculators must implement.
     *
     * @param U_L: Left interface state (conservative variables)
     * @param U_R: Right interface state (conservative variables)
     * @param physics: Physics equation object
     * @param direction: 0=x, 1=y, 2=z
     * @return: Numerical flux vector at the interface
     */
    virtual Eigen::VectorXd compute_flux(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        const physics::PhysicsBase& physics,
        int direction
    ) const = 0;

    /**
     * Compute maximum wave speed for stability analysis.
     *
     * Used for CFL condition calculation.
     *
     * @param U_L: Left state
     * @param U_R: Right state
     * @param physics: Physics equation object
     * @param direction: Direction (0=x, 1=y, 2=z)
     * @return: Maximum wave speed
     */
    virtual double compute_max_wave_speed(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        const physics::PhysicsBase& physics,
        int direction
    ) const = 0;

    /**
     * Get name of this flux calculator.
     */
    std::string name() const { return name_; }

    /**
     * Get flux type category.
     */
    std::string flux_type() const { return flux_type_; }

    /**
     * Check if this flux calculator is physics-agnostic.
     *
     * Physics-agnostic methods work with any PhysicsBase implementation.
     * Non-agnostic methods (e.g., HLLD) require specific physics types (e.g., MHD).
     *
     * @return: True if physics-agnostic
     */
    virtual bool is_physics_agnostic() const { return true; }

    /**
     * Check if this flux calculator supports the given physics type.
     *
     * @param physics_type: Type string ("euler", "mhd_advanced", etc.)
     * @return: True if supported
     */
    virtual bool supports_physics(const std::string& physics_type) const {
        return is_physics_agnostic();
    }

    /**
     * Check if this flux calculator needs explicit wave speed computation.
     *
     * Some methods (like Lax-Friedrichs) explicitly use wave speeds,
     * while others (like exact Riemann solvers) compute them implicitly.
     *
     * @return: True if explicit wave speed is needed
     */
    virtual bool needs_wave_speed() const { return false; }

protected:
    std::string name_;       ///< Name of the flux calculator
    std::string flux_type_;  ///< Type category
};

} // namespace fvm3d::spatial
