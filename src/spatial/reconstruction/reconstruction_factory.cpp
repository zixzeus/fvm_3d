#include "spatial/reconstruction/reconstruction_factory.hpp"
#include "spatial/reconstruction/constant_reconstruction.hpp"
#include "spatial/reconstruction/muscl_reconstruction.hpp"
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <functional>

namespace fvm3d::spatial {

std::unique_ptr<ReconstructionMethod> ReconstructionFactory::create(
    const std::string& name,
    int num_vars,
    const std::string& limiter,
    double kappa
) {
    using CreatorFunc = std::function<std::unique_ptr<ReconstructionMethod>(int, const std::string&, double)>;

    // Registry with lambda creators for each reconstruction method
    static const std::unordered_map<std::string, CreatorFunc> registry = {
        // Constant reconstruction and aliases
        {"constant", [](int n, const std::string&, double) {
            return std::make_unique<ConstantReconstruction>(n);
        }},
        {"1st", [](int n, const std::string&, double) {
            return std::make_unique<ConstantReconstruction>(n);
        }},
        {"first", [](int n, const std::string&, double) {
            return std::make_unique<ConstantReconstruction>(n);
        }},

        // MUSCL reconstruction and aliases
        {"muscl", [](int n, const std::string& lim, double k) {
            return std::make_unique<MUSCLReconstruction>(n, lim, k);
        }},
        {"2nd", [](int n, const std::string& lim, double k) {
            return std::make_unique<MUSCLReconstruction>(n, lim, k);
        }},
        {"second", [](int n, const std::string& lim, double k) {
            return std::make_unique<MUSCLReconstruction>(n, lim, k);
        }}
    };

    // Convert name to lowercase for case-insensitive matching
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), ::tolower);

    auto it = registry.find(name_lower);
    if (it != registry.end()) {
        return it->second(num_vars, limiter, kappa);
    }

    throw std::invalid_argument(
        "ReconstructionFactory: Unknown reconstruction method '" + name + "'. "
        "Supported methods: constant, muscl"
    );
}

std::vector<std::string> ReconstructionFactory::supported_methods() {
    return {
        "constant",  // First-order piecewise constant
        "muscl"      // Second-order MUSCL with limiters
        // Future: "weno3", "weno5", "eno", "ppm"
    };
}

bool ReconstructionFactory::is_available(const std::string& name) {
    static const std::unordered_map<std::string, bool> valid_names = {
        {"constant", true}, {"1st", true}, {"first", true},
        {"muscl", true}, {"2nd", true}, {"second", true}
    };

    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), ::tolower);

    return valid_names.count(name_lower) > 0;
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
