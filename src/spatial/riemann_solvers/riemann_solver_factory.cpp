#include "spatial/riemann_solvers/riemann_solver_factory.hpp"
#include "spatial/riemann_solvers/riemann_hll.hpp"
#include "spatial/riemann_solvers/riemann_hllc.hpp"
#include "spatial/riemann_solvers/riemann_hlld.hpp"
#include "physics/euler3d.hpp"
#include "physics/resistive_mhd3d_advanced.hpp"
#include "physics/physics_factory.hpp"
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <functional>

namespace fvm3d::spatial {

// Preferred method: Create Riemann solver with provided physics object
std::unique_ptr<RiemannSolver> RiemannSolverFactory::create(
    const std::string& name,
    const std::shared_ptr<physics::PhysicsBase>& physics
) {
    using SolverCreator = std::function<std::unique_ptr<RiemannSolver>(std::shared_ptr<physics::PhysicsBase>)>;

    // Registry for Riemann solver creation
    static const std::unordered_map<std::string, SolverCreator> solver_registry = {
        {"hll", [](auto phys) {
            return std::make_unique<HLLSolver>(phys);
        }},
        {"hllc", [](auto phys) {
            return std::make_unique<HLLCSolver>(phys);
        }},
        {"hlld", [](auto phys) {
            // HLLD requires MHD physics
            auto mhd_phys = std::dynamic_pointer_cast<physics::AdvancedResistiveMHD3D>(phys);
            if (!mhd_phys) {
                throw std::invalid_argument("HLLD solver requires MHD physics");
            }
            return std::make_unique<HLLDSolver>(mhd_phys);
        }}
    };

    // Convert solver name to lowercase
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    // Create Riemann solver
    auto solver_it = solver_registry.find(name_lower);
    if (solver_it == solver_registry.end()) {
        throw std::invalid_argument("Unknown Riemann solver: " + name);
    }

    return solver_it->second(physics);
}

// Deprecated method: Create Riemann solver with internal physics object creation
std::unique_ptr<RiemannSolver> RiemannSolverFactory::create(
    const std::string& name,
    const std::string& physics_type,
    int num_vars
) {
    // Create physics object using PhysicsFactory
    // Note: This creates a new physics object each time - prefer using create(name, physics) to avoid duplication
    auto physics = physics::PhysicsFactory::create(physics_type);

    // Delegate to the preferred method
    return create(name, physics);
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
