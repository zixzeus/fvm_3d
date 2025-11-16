#pragma once

#include "physics/physics_base.hpp"
#include "physics/euler3d.hpp"
#include "physics/resistive_mhd3d.hpp"
#include "physics/resistive_mhd3d_advanced.hpp"
#include <memory>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <functional>
#include <algorithm>

namespace fvm3d::physics {

/**
 * Factory for creating physics objects using registry pattern.
 *
 * Uses a map-based approach instead of if-else chains for better
 * extensibility and maintainability.
 */
class PhysicsFactory {
public:
    using CreatorFunc = std::function<std::shared_ptr<PhysicsBase>(double)>;

    /**
     * Create a physics object by type name.
     *
     * @param physics_type Type of physics: "euler", "mhd", "mhd_advanced"
     * @param resistivity Resistivity coefficient for MHD (default: 0.0)
     * @return Shared pointer to physics object
     *
     * @throws std::invalid_argument if physics_type is unknown
     */
    static std::shared_ptr<PhysicsBase> create(
        const std::string& physics_type,
        double resistivity = 0.0
    ) {
        // Convert to lowercase for case-insensitive matching
        std::string type_lower = physics_type;
        std::transform(type_lower.begin(), type_lower.end(),
                      type_lower.begin(), ::tolower);

        auto& registry = get_registry();
        auto it = registry.find(type_lower);

        if (it != registry.end()) {
            return it->second(resistivity);
        }

        throw std::invalid_argument(
            "Unknown physics type: '" + physics_type + "'. " +
            "Supported types: " + get_supported_types()
        );
    }

    /**
     * Create advanced MHD physics with custom resistivity model.
     *
     * @param resistivity_model Custom resistivity model configuration
     * @param glm_params GLM divergence cleaning parameters
     * @return Shared pointer to AdvancedResistiveMHD3D
     */
    static std::shared_ptr<AdvancedResistiveMHD3D> create_advanced_mhd(
        const AdvancedResistiveMHD3D::ResistivityModel& resistivity_model,
        const AdvancedResistiveMHD3D::GLMParameters& glm_params =
            AdvancedResistiveMHD3D::GLMParameters()
    ) {
        return std::make_shared<AdvancedResistiveMHD3D>(resistivity_model, glm_params);
    }

    /**
     * Check if a physics type is supported.
     *
     * @param physics_type Type to check (case-insensitive)
     * @return True if supported
     */
    static bool is_supported(const std::string& physics_type) {
        std::string type_lower = physics_type;
        std::transform(type_lower.begin(), type_lower.end(),
                      type_lower.begin(), ::tolower);
        return get_registry().count(type_lower) > 0;
    }

    /**
     * Get comma-separated list of supported physics types.
     */
    static std::string get_supported_types() {
        std::string types;
        for (const auto& pair : get_registry()) {
            if (!types.empty()) types += ", ";
            types += pair.first;
        }
        return types;
    }

private:
    /**
     * Get the physics creator registry (singleton pattern).
     * Registry is initialized on first access.
     */
    static std::unordered_map<std::string, CreatorFunc>& get_registry() {
        static std::unordered_map<std::string, CreatorFunc> registry = {
            {"euler", [](double) {
                return std::make_shared<EulerEquations3D>();
            }},
            {"mhd", [](double resistivity) {
                return std::make_shared<ResistiveMHD3D>(resistivity);
            }},
            {"mhd_advanced", [](double) {
                return std::make_shared<AdvancedResistiveMHD3D>();
            }}
        };
        return registry;
    }
};

} // namespace fvm3d::physics
