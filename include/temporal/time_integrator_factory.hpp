#pragma once

#include "time_integrator.hpp"
#include <string>
#include <memory>
#include <vector>

namespace fvm3d::temporal {

/**
 * Factory for creating time integrator schemes.
 */
class TimeIntegratorFactory {
public:
    /**
     * Create a time integrator by name.
     * Supported names: "euler", "forward_euler", "rk2", "rk3"
     */
    static std::unique_ptr<TimeIntegrator> create(const std::string& name);

    /**
     * List all supported schemes.
     */
    static std::vector<std::string> supported_schemes();
};

} // namespace fvm3d::temporal
