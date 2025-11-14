#include "spatial/reconstruction/reconstruction_factory.hpp"
#include "spatial/reconstruction/constant_reconstruction.hpp"
#include "spatial/reconstruction/muscl_reconstruction.hpp"
#include <iostream>
#include <algorithm>
#include <stdexcept>

namespace fvm3d::spatial {

std::unique_ptr<ReconstructionMethod> ReconstructionFactory::create(
    const std::string& name,
    int num_vars,
    const std::string& limiter,
    double kappa
) {
    // Convert name to lowercase for case-insensitive matching
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), ::tolower);

    if (name_lower == "constant" || name_lower == "1st" || name_lower == "first") {
        return std::make_unique<ConstantReconstruction>(num_vars);
    }
    else if (name_lower == "muscl" || name_lower == "2nd" || name_lower == "second") {
        return std::make_unique<MUSCLReconstruction>(num_vars, limiter, kappa);
    }
    else {
        throw std::invalid_argument(
            "ReconstructionFactory: Unknown reconstruction method '" + name + "'. "
            "Supported methods: constant, muscl"
        );
    }
}

std::vector<std::string> ReconstructionFactory::supported_methods() {
    return {
        "constant",  // First-order piecewise constant
        "muscl"      // Second-order MUSCL with limiters
        // Future: "weno3", "weno5", "eno", "ppm"
    };
}

bool ReconstructionFactory::is_available(const std::string& name) {
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), ::tolower);

    auto methods = supported_methods();
    return std::find(methods.begin(), methods.end(), name_lower) != methods.end() ||
           name_lower == "1st" || name_lower == "first" ||
           name_lower == "2nd" || name_lower == "second";
}

void ReconstructionFactory::print_available() {
    std::cout << "Available reconstruction methods:\n";
    std::cout << "  - constant (1st order, piecewise constant)\n";
    std::cout << "      Aliases: 1st, first\n";
    std::cout << "      Properties: TVD, no limiters, robust\n";
    std::cout << "      Ghost cells: 1\n";
    std::cout << "\n";
    std::cout << "  - muscl (2nd order, MUSCL with limiters)\n";
    std::cout << "      Aliases: 2nd, second\n";
    std::cout << "      Properties: TVD with limiters, higher accuracy\n";
    std::cout << "      Limiters: minmod, van_leer, superbee, mc\n";
    std::cout << "      Ghost cells: 2\n";
    std::cout << "      Parameters: kappa (default 1/3)\n";
    std::cout << "\n";
    std::cout << "Future methods (planned):\n";
    std::cout << "  - weno3, weno5 (Weighted ENO, 3rd/5th order)\n";
    std::cout << "  - eno (Essentially Non-Oscillatory)\n";
    std::cout << "  - ppm (Piecewise Parabolic Method)\n";
}

std::vector<std::string> ReconstructionFactory::supported_limiters() {
    return {
        "minmod",     // Most dissipative, most robust
        "van_leer",   // Good balance
        "superbee",   // Least dissipative, can be less stable
        "mc"          // Monotonized Central, good compromise
    };
}

} // namespace fvm3d::spatial
