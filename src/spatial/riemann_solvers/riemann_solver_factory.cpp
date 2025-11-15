#include "spatial/riemann_solvers/riemann_solver_factory.hpp"
#include "spatial/riemann_solvers/riemann_hll.hpp"
#include "spatial/riemann_solvers/riemann_hllc.hpp"
#include "spatial/riemann_solvers/riemann_hlld.hpp"
#include "physics/euler3d.hpp"
#include "physics/resistive_mhd3d_advanced.hpp"
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <memory>

namespace fvm3d::spatial {

std::unique_ptr<RiemannSolver> RiemannSolverFactory::create(
    const std::string& name,
    const std::string& physics_type,
    int num_vars
) {
    // Convert solver name to lowercase
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    // Create physics object based on physics_type
    std::shared_ptr<physics::PhysicsBase> physics;

    if (physics_type == "euler") {
        physics = std::make_shared<physics::EulerEquations3D>();
    } else if (physics_type == "mhd_advanced" || physics_type == "mhd") {
        // Advanced resistive MHD with GLM divergence cleaning
        physics::AdvancedResistiveMHD3D::ResistivityModel resistivity;
        resistivity.eta0 = 1e-3;               // Background resistivity
        resistivity.eta1 = 0.01667;            // Enhanced resistivity
        resistivity.localization_scale = 1.0;  // Localization width

        physics::AdvancedResistiveMHD3D::GLMParameters glm(0.2, 0.2);  // ch=0.2, cr=0.2

        physics = std::make_shared<physics::AdvancedResistiveMHD3D>(resistivity, glm);
    } else {
        throw std::invalid_argument("Unknown physics type: " + physics_type);
    }

    // Create Riemann solver with physics object
    if (name_lower == "hll") {
        return std::make_unique<HLLSolver>(physics);
    } else if (name_lower == "hllc") {
        return std::make_unique<HLLCSolver>(physics);
    } else if (name_lower == "hlld") {
        // HLLD only works with MHD physics
        if (physics_type != "mhd_advanced" && physics_type != "mhd") {
            throw std::invalid_argument("HLLD solver requires MHD physics, got: " + physics_type);
        }
        auto mhd_physics = std::dynamic_pointer_cast<physics::AdvancedResistiveMHD3D>(physics);
        return std::make_unique<HLLDSolver>(mhd_physics);
    } else {
        throw std::invalid_argument("Unknown Riemann solver: " + name);
    }
}

std::vector<std::string> RiemannSolverFactory::supported_solvers() {
    return {
        "hll",
        "hllc",
        "hlld"
    };
}

void RiemannSolverFactory::print_available_solvers() {
    std::cout << "Available Riemann solvers via RiemannSolverFactory:\n";
    std::cout << "  - hll\n";
    std::cout << "  - hllc\n";
    std::cout << "  - hlld (MHD only)\n";
    std::cout << "\nNote: For unified flux calculation interface, use FluxCalculatorFactory.\n";
}

} // namespace fvm3d::spatial
