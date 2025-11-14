#pragma once

#include "riemann_solver.hpp"
#include <memory>
#include <string>
#include <vector>

namespace fvm3d::spatial {

/**
 * Factory class for creating Riemann solvers by name.
 * Allows dynamic selection of solver algorithms at runtime.
 */
class RiemannSolverFactory {
public:
    /**
     * Create a Riemann solver by name.
     * Supported solvers: "laxfriedrichs" (or "lf"), "hll", "hllc", "hlld"
     *
     * @param name: name of the solver
     * @param physics_type: type of physics ("euler" or "mhd_advanced")
     * @param num_vars: number of variables (5 for Euler, 9 for MHD)
     * @return: unique_ptr to a RiemannSolver instance
     * @throws std::invalid_argument if name or physics_type is unknown
     */
    static std::unique_ptr<RiemannSolver> create(
        const std::string& name,
        const std::string& physics_type = "euler",
        int num_vars = 5
    );

    /**
     * Get list of supported Riemann solvers.
     */
    static std::vector<std::string> supported_solvers();

    /**
     * Print available solvers to stdout.
     */
    static void print_available_solvers();
};

} // namespace fvm3d::spatial
