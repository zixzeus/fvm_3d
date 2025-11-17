#pragma once

#include "flux_calculator_base.hpp"
#include "physics/physics_base.hpp"
#include <memory>
#include <string>
#include <vector>

namespace fvm3d::spatial {

/**
 * Factory for creating flux calculators.
 *
 * Provides unified interface to create different types of flux calculators:
 * - Riemann solver-based (HLL, HLLC, HLLD, Lax-Friedrichs)
 * - Central schemes (future)
 * - Hybrid methods (future)
 *
 * Usage (preferred):
 *   auto physics = PhysicsFactory::create("euler");
 *   auto flux_calc = FluxCalculatorFactory::create("hllc", physics);
 *   Eigen::VectorXd flux = flux_calc->compute_flux(U_L, U_R, *physics, 0);
 */
class FluxCalculatorFactory {
public:
    /**
     * Create a flux calculator with existing physics object (PREFERRED).
     *
     * This method avoids duplicate physics object creation and allows
     * the flux calculator to use the same physics configuration as the solver.
     *
     * Supported methods:
     * - Riemann solvers: "laxfriedrichs", "lf", "hll", "hllc", "hlld"
     * - Central schemes: "central" (future)
     * - Hybrid: "ausm" (future)
     *
     * @param name: Name of the flux calculator
     * @param physics: Shared pointer to physics object
     * @return: Unique pointer to FluxCalculator instance
     * @throws std::invalid_argument if name is unknown
     */
    static std::unique_ptr<FluxCalculator> create(
        const std::string& name,
        const std::shared_ptr<physics::PhysicsBase>& physics
    );

    /**
     * Create a flux calculator by name (DEPRECATED).
     *
     * This method creates a new physics object internally, which may have
     * different parameters than the solver's physics object. Use the
     * new create(name, physics) method instead.
     *
     * @param name: Name of the flux calculator
     * @param physics_type: Type of physics ("euler", "mhd_advanced", etc.)
     * @param num_vars: Number of variables (5 for Euler, 9 for MHD)
     * @return: Unique pointer to FluxCalculator instance
     * @throws std::invalid_argument if name or physics_type is unknown
     * @deprecated Use create(name, physics) instead
     */
    static std::unique_ptr<FluxCalculator> create(
        const std::string& name,
        const std::string& physics_type,
        int num_vars = 5
    );

    /**
     * Get list of supported flux calculators.
     */
    static std::vector<std::string> supported_calculators();

    /**
     * Check if a flux calculator is available.
     *
     * @param name: Flux calculator name
     * @return: True if available
     */
    static bool is_available(const std::string& name);

    /**
     * Print available flux calculators to stdout.
     */
    static void print_available();
};

} // namespace fvm3d::spatial
