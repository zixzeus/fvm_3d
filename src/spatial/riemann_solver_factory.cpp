#include "spatial/riemann_solver_factory.hpp"
#include "spatial/riemann_laxfriedrichs.hpp"
#include "spatial/riemann_hll.hpp"
#include "spatial/riemann_hllc.hpp"
#include <stdexcept>
#include <iostream>
#include <algorithm>

namespace fvm3d::spatial {

std::unique_ptr<RiemannSolver> RiemannSolverFactory::create(const std::string& name) {
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (name_lower == "laxfriedrichs" || name_lower == "lf") {
        return std::make_unique<LaxFriedrichsSolver>();
    } else if (name_lower == "hll") {
        return std::make_unique<HLLSolver>();
    } else if (name_lower == "hllc") {
        return std::make_unique<HLLCSolver>();
    } else {
        throw std::invalid_argument("Unknown Riemann solver: " + name);
    }
}

std::vector<std::string> RiemannSolverFactory::supported_solvers() {
    return {
        "laxfriedrichs (lf)",
        "hll",
        "hllc"
    };
}

void RiemannSolverFactory::print_available_solvers() {
    std::cout << "Available Riemann solvers:\n";
    for (const auto& solver : supported_solvers()) {
        std::cout << "  - " << solver << "\n";
    }
}

} // namespace fvm3d::spatial
