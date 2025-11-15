#include "spatial/flux_calculation/flux_calculator_factory.hpp"
#include "spatial/flux_calculation/lax_friedrichs.hpp"
#include "spatial/flux_calculation/riemann_flux.hpp"
#include "physics/euler3d.hpp"
#include "physics/resistive_mhd3d_advanced.hpp"
#include <algorithm>
#include <iostream>

namespace fvm3d::spatial {

std::unique_ptr<FluxCalculator> FluxCalculatorFactory::create(
    const std::string& name,
    const std::string& physics_type,
    int num_vars
) {
    // Convert name to lowercase for case-insensitive matching
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

    // Lax-Friedrichs flux calculator
    if (name_lower == "laxfriedrichs" || name_lower == "lf") {
        return std::make_unique<LaxFriedrichsFlux>(physics);
    }

    // Riemann solver-based flux calculators
    if (name_lower == "hll" || name_lower == "hllc" || name_lower == "hlld") {
        return std::make_unique<RiemannFlux>(name_lower, physics);
    }

    // Unknown flux calculator
    throw std::invalid_argument("Unknown flux calculator: " + name);

    // Future: add non-Riemann methods here
    // if (name_lower == "central") {
    //     return std::make_unique<CentralFluxScheme>(physics);
    // } else if (name_lower == "ausm") {
    //     return std::make_unique<AUSMFluxCalculator>(physics);
    // }
}

std::vector<std::string> FluxCalculatorFactory::supported_calculators() {
    std::vector<std::string> calculators;

    // Riemann solver-based flux calculators
    calculators.push_back("laxfriedrichs (lf)");
    calculators.push_back("hll");
    calculators.push_back("hllc");
    calculators.push_back("hlld (MHD only)");

    // Future: add additional non-Riemann methods
    // calculators.push_back("central");
    // calculators.push_back("ausm");

    return calculators;
}

bool FluxCalculatorFactory::is_available(const std::string& name) {
    auto supported = supported_calculators();
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    for (const auto& calc : supported) {
        std::string calc_lower = calc;
        std::transform(calc_lower.begin(), calc_lower.end(), calc_lower.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        if (calc_lower.find(name_lower) != std::string::npos) {
            return true;
        }
    }
    return false;
}

void FluxCalculatorFactory::print_available() {
    std::cout << "Available flux calculators:\n";
    std::cout << "\nRiemann Solver-Based Methods:\n";
    for (const auto& calc : supported_calculators()) {
        std::cout << "  - " << calc << "\n";
    }

    std::cout << "\nFuture Methods (not yet implemented):\n";
    std::cout << "  - central: Central differencing scheme\n";
    std::cout << "  - ausm: Advection Upstream Splitting Method\n";
    std::cout << "  - kfvs: Kinetic Flux Vector Splitting\n";
}

} // namespace fvm3d::spatial
