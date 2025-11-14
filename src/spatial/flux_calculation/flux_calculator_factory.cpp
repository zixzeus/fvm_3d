#include "spatial/flux_calculation/flux_calculator_factory.hpp"
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

    // For now, all flux calculators are Riemann solver-based
    // Create Riemann solver using existing factory
    auto riemann_solver = RiemannSolverFactory::create(name, physics_type, num_vars);

    // Wrap in adapter to provide FluxCalculator interface
    return std::make_unique<RiemannFluxAdapter>(std::move(riemann_solver));

    // Future: add non-Riemann methods here
    // if (name_lower == "central") {
    //     return std::make_unique<CentralFluxScheme>(physics_type);
    // } else if (name_lower == "ausm") {
    //     return std::make_unique<AUSMFluxCalculator>(physics_type);
    // }
}

std::vector<std::string> FluxCalculatorFactory::supported_calculators() {
    // Currently delegates to Riemann solver factory
    // Future: add additional non-Riemann methods
    return RiemannSolverFactory::supported_solvers();
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
