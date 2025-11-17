#include "spatial/flux_calculation/flux_calculator_factory.hpp"
#include "spatial/flux_calculation/lax_friedrichs.hpp"
#include "spatial/flux_calculation/riemann_flux.hpp"
#include "physics/physics_factory.hpp"
#include <algorithm>
#include <iostream>

namespace fvm3d::spatial {

// Preferred method: accepts existing physics object
std::unique_ptr<FluxCalculator> FluxCalculatorFactory::create(
    const std::string& name,
    const std::shared_ptr<physics::PhysicsBase>& physics
) {
    // Convert name to lowercase for case-insensitive matching
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

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

// Deprecated method: creates physics object internally
std::unique_ptr<FluxCalculator> FluxCalculatorFactory::create(
    const std::string& name,
    const std::string& physics_type,
    int num_vars
) {
    // Create physics object using PhysicsFactory
    // Note: Uses default parameters, may not match solver's physics configuration
    auto physics = physics::PhysicsFactory::create(physics_type, 0.0);

    // Delegate to new create() method
    return create(name, physics);
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
