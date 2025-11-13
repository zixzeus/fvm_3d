#include "temporal/time_integrator_factory.hpp"
#include "temporal/forward_euler.hpp"
#include "temporal/rk_integrators.hpp"
#include <algorithm>
#include <stdexcept>

namespace fvm3d::temporal {

std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create(const std::string& name) {
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (name_lower == "euler" || name_lower == "forward_euler") {
        return std::make_unique<ForwardEuler>();
    } else if (name_lower == "rk2") {
        return std::make_unique<RK2>();
    } else if (name_lower == "rk3") {
        return std::make_unique<RK3>();
    } else {
        throw std::invalid_argument("Unknown time integrator: " + name);
    }
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
