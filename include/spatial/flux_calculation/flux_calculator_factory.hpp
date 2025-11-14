#pragma once

#include "flux_calculator_base.hpp"
#include "riemann_flux_adapter.hpp"
#include "spatial/riemann_solver_factory.hpp"
#include <memory>
#include <string>
#include <vector>

namespace fvm3d::spatial {

/**
 * Factory for creating flux calculators.
 *
 * Provides unified interface to create different types of flux calculators:
 * - Riemann solver-based (via RiemannSolverFactory + adapter)
 * - Central schemes
 * - Hybrid methods
 *
 * Usage:
 *   auto flux_calc = FluxCalculatorFactory::create("hllc", "mhd_advanced");
 *   Eigen::VectorXd flux = flux_calc->compute_flux(U_L, U_R, physics, 0);
 */
class FluxCalculatorFactory {
public:
    /**
     * Create a flux calculator by name.
     *
     * Supported methods:
     * - Riemann solvers: "laxfriedrichs", "lf", "hll", "hllc", "hlld"
     * - Central schemes: "central" (future)
     * - Hybrid: "ausm" (future)
     *
     * @param name: Name of the flux calculator
     * @param physics_type: Type of physics ("euler", "mhd_advanced", etc.)
     * @param num_vars: Number of variables (5 for Euler, 9 for MHD)
     * @return: Unique pointer to FluxCalculator instance
     * @throws std::invalid_argument if name or physics_type is unknown
     */
    static std::unique_ptr<FluxCalculator> create(
        const std::string& name,
        const std::string& physics_type = "euler",
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
