#include "temporal/time_integrator_factory.hpp"
#include "temporal/forward_euler.hpp"
#include "temporal/rk_integrators.hpp"
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <functional>

namespace fvm3d::temporal {

std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create(const std::string& name) {
    using CreatorFunc = std::function<std::unique_ptr<TimeIntegrator>()>;

    // Registry with lambda creators for each time integrator
    static const std::unordered_map<std::string, CreatorFunc> registry = {
        // Forward Euler and aliases
        {"euler", []() {
            return std::make_unique<ForwardEuler>();
        }},
        {"forward_euler", []() {
            return std::make_unique<ForwardEuler>();
        }},

        // Runge-Kutta schemes
        {"rk2", []() {
            return std::make_unique<RK2>();
        }},
        {"rk3", []() {
            return std::make_unique<RK3>();
        }}
    };

    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    auto it = registry.find(name_lower);
    if (it != registry.end()) {
        return it->second();
    }

    throw std::invalid_argument("Unknown time integrator: " + name);
}

std::vector<std::string> TimeIntegratorFactory::supported_schemes() {
    return {
        "euler",
        "forward_euler",
        "rk2",
        "rk3"
    };
}

} // namespace fvm3d::temporal
