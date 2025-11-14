#pragma once

#include "physics/physics_base.hpp"
#include <Eigen/Dense>

namespace fvm3d::physics {

/**
 * Base class for MagnetoHydroDynamic (MHD) physics.
 *
 * Provides common functionality for MHD equations including:
 * - Conservative to primitive variable conversion
 * - Primitive to conservative variable conversion
 * - Sound speed and Alfvén speed calculations
 * - Fast magnetosonic speed calculations
 *
 * This class bridges the gap between the generic PhysicsBase interface
 * and the specific directional flux methods (flux_x, flux_y, flux_z)
 * commonly used in MHD implementations.
 */
class MHDBase : public ConservationLaw {
public:
    /**
     * Constructor.
     * @param name: Name of the MHD variant
     * @param num_variables: Number of variables (8 for pure MHD, 9 for GLM-MHD)
     */
    MHDBase(const std::string& name, int num_variables)
        : ConservationLaw(name, num_variables, 3) {}

    virtual ~MHDBase() = default;

    // ========== MHD-Specific Pure Virtual Methods ==========

    /**
     * Convert conservative to primitive variables.
     * MHD primitive: [rho, u, v, w, p, Bx, By, Bz]
     */
    virtual void conservative_to_primitive(
        const Eigen::VectorXd& U,
        double& rho, double& u, double& v, double& w, double& p,
        double& Bx, double& By, double& Bz
    ) const = 0;

    /**
     * Convert primitive to conservative variables.
     * MHD conservative: [rho, rho*u, rho*v, rho*w, E, Bx, By, Bz, (psi)]
     */
    virtual void primitive_to_conservative(
        double rho, double u, double v, double w, double p,
        double Bx, double By, double Bz,
        Eigen::VectorXd& U
    ) const = 0;

    /**
     * Compute flux in X direction.
     */
    virtual Eigen::VectorXd flux_x(const Eigen::VectorXd& U) const = 0;

    /**
     * Compute flux in Y direction.
     */
    virtual Eigen::VectorXd flux_y(const Eigen::VectorXd& U) const = 0;

    /**
     * Compute flux in Z direction.
     */
    virtual Eigen::VectorXd flux_z(const Eigen::VectorXd& U) const = 0;

    /**
     * Compute sound speed (gas pressure waves).
     */
    virtual double sound_speed(double rho, double p) const = 0;

    /**
     * Compute Alfvén speed (magnetic waves).
     */
    virtual double alfven_speed(double B_component, double rho) const = 0;

    /**
     * Compute fast magnetosonic speed (fastest MHD wave).
     */
    virtual double fast_speed(
        double rho, double p,
        double Bx, double By, double Bz
    ) const = 0;

    // ========== PhysicsBase Interface Implementation ==========

    /**
     * Unified flux computation interface.
     * Routes to appropriate directional flux method.
     */
    Eigen::VectorXd compute_flux(
        const Eigen::VectorXd& U,
        int direction
    ) const override {
        switch (direction) {
            case 0: return flux_x(U);
            case 1: return flux_y(U);
            case 2: return flux_z(U);
            default:
                throw std::invalid_argument(
                    "Invalid direction: " + std::to_string(direction) +
                    " (must be 0, 1, or 2)"
                );
        }
    }

    /**
     * Validate MHD state.
     * Checks for positive density and pressure.
     */
    bool validate_state(const Eigen::VectorXd& U) const override {
        // Basic checks from base class
        if (!PhysicsBase::validate_state(U)) {
            return false;
        }

        // MHD-specific checks
        double rho, u, v, w, p, Bx, By, Bz;
        try {
            conservative_to_primitive(U, rho, u, v, w, p, Bx, By, Bz);
        } catch (...) {
            return false;
        }

        // Density must be positive
        if (rho <= 0.0) {
            return false;
        }

        // Pressure must be positive (or at least non-negative)
        if (p < 0.0) {
            return false;
        }

        return true;
    }

    /**
     * MHD variable names.
     */
    std::vector<std::string> get_variable_names() const override {
        std::vector<std::string> names = {
            "rho", "rho_u", "rho_v", "rho_w", "E",
            "Bx", "By", "Bz"
        };

        // Add GLM field if 9 variables
        if (num_variables_ == 9) {
            names.push_back("psi");
        }

        return names;
    }

    /**
     * MHD primitive variable names.
     */
    std::vector<std::string> get_primitive_names() const override {
        std::vector<std::string> names = {
            "rho", "u", "v", "w", "p",
            "Bx", "By", "Bz"
        };

        if (num_variables_ == 9) {
            names.push_back("psi");
        }

        return names;
    }

    /**
     * MHD has pressure.
     */
    bool has_pressure() const override {
        return true;
    }

    /**
     * MHD can compute temperature if equation of state is known.
     */
    bool has_temperature() const override {
        return true;
    }
};

} // namespace fvm3d::physics
