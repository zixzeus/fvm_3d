#include "spatial/riemann_solver_factory.hpp"
#include "spatial/riemann_laxfriedrichs.hpp"
#include "spatial/riemann_hll.hpp"
#include "spatial/riemann_hllc.hpp"
#include "spatial/riemann_hlld.hpp"
#include "physics/resistive_mhd3d_advanced.hpp"
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <memory>

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
    } else if (name_lower == "hlld") {
        // HLLD solver with AdvancedResistiveMHD3D physics (gamma=5/3, GLM ch=0.2, cr=0.2)
        physics::AdvancedResistiveMHD3D::ResistivityModel resistivity;
        resistivity.eta0 = 1e-3;
        resistivity.eta1 = 0.01667;
        resistivity.localization_scale = 1.0;

        physics::AdvancedResistiveMHD3D::GLMParameters glm(0.2, 0.2);

        auto mhd_physics = std::make_shared<physics::AdvancedResistiveMHD3D>(resistivity, glm);
        return std::make_unique<HLLDSolver>(mhd_physics);
    } else {
        throw std::invalid_argument("Unknown Riemann solver: " + name);
    }
}

std::vector<std::string> RiemannSolverFactory::supported_solvers() {
    return {
        "laxfriedrichs (lf)",
        "hll",
        "hllc",
        "hlld"
    };
}

void RiemannSolverFactory::print_available_solvers() {
    std::cout << "Available Riemann solvers:\n";
    for (const auto& solver : supported_solvers()) {
        std::cout << "  - " << solver << "\n";
    }
}

} // namespace fvm3d::spatial
