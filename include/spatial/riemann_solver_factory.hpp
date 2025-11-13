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
     * Supported solvers: "laxfriedrichs" (or "lf"), "hll", "hllc"
     *
     * @param name: name of the solver
     * @return: unique_ptr to a RiemannSolver instance
     * @throws std::invalid_argument if name is unknown
     */
    static std::unique_ptr<RiemannSolver> create(const std::string& name);

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
