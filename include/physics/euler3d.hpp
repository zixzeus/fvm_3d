#pragma once

#include <Eigen/Dense>
#include <cmath>
#include "physics_base.hpp"

namespace fvm3d::physics {

/**
 * Compressible Euler equations for ideal gas.
 * Conservative variables: U = [rho, rho_u, rho_v, rho_w, E]
 * Primitive variables: V = [rho, u, v, w, p]
 */
class EulerEquations3D : public ConservationLaw {
public:
    static constexpr int nvars = 5;
    static constexpr double GAMMA = 1.4;  // Adiabatic index for air
    static constexpr double RHO_FLOOR = 1e-10;
    static constexpr double P_FLOOR = 1e-11;

    // ========== Constructor ==========

    explicit EulerEquations3D()
        : ConservationLaw("Euler3D", nvars, 3) {}

    // ========== PhysicsBase Interface Implementation ==========

    /**
     * Convert conservative to primitive variables (PhysicsBase interface).
     * @param U: Conservative variables [rho, rho_u, rho_v, rho_w, E]
     * @return V: Primitive variables [rho, u, v, w, p]
     */
    Eigen::VectorXd conservative_to_primitive(const Eigen::VectorXd& U) const override {
        double rho, u, v, w, p;
        conservative_to_primitive(U, rho, u, v, w, p);

        Eigen::VectorXd V(nvars);
        V << rho, u, v, w, p;
        return V;
    }

    /**
     * Convert primitive to conservative variables (PhysicsBase interface).
     * @param V: Primitive variables [rho, u, v, w, p]
     * @return U: Conservative variables [rho, rho_u, rho_v, rho_w, E]
     */
    Eigen::VectorXd primitive_to_conservative(const Eigen::VectorXd& V) const override {
        Eigen::VectorXd U(nvars);
        primitive_to_conservative(V(0), V(1), V(2), V(3), V(4), U);
        return U;
    }

    /**
     * Compute flux in given direction (PhysicsBase interface).
     */
    Eigen::VectorXd compute_flux(const Eigen::VectorXd& U, int direction) const override {
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

    // ========== Legacy Interface (for backwards compatibility) ==========

    /**
     * Convert conservative variables to primitive variables.
     * @param U: conservative state [rho, rho_u, rho_v, rho_w, E]
     * @param rho, u, v, w, p: primitive variables (output)
     */
    void conservative_to_primitive(
        const Eigen::VectorXd& U,
        double& rho, double& u, double& v, double& w, double& p
    ) const {
        rho = std::max(U(0), RHO_FLOOR);
        u = U(1) / rho;
        v = U(2) / rho;
        w = U(3) / rho;

        // Direct kinetic energy calculation to avoid division issues
        double kinetic_energy = 0.5 * (U(1)*U(1) + U(2)*U(2) + U(3)*U(3)) / rho;
        double internal_energy = U(4) / rho - kinetic_energy;
        p = std::max((GAMMA - 1.0) * rho * internal_energy, P_FLOOR);
    }

    /**
     * Convert primitive variables to conservative variables.
     */
    void primitive_to_conservative(
        double rho, double u, double v, double w, double p,
        Eigen::VectorXd& U
    ) const {
        rho = std::max(rho, RHO_FLOOR);
        p = std::max(p, P_FLOOR);

        U(0) = rho;
        U(1) = rho * u;
        U(2) = rho * v;
        U(3) = rho * w;

        double kinetic_energy = 0.5 * rho * (u*u + v*v + w*w);
        double internal_energy = p / ((GAMMA - 1.0) * rho);
        U(4) = rho * internal_energy + kinetic_energy;
    }

    /**
     * Compute sound speed: a = sqrt(gamma * p / rho)
     */
    double sound_speed(double rho, double p) const {
        rho = std::max(rho, RHO_FLOOR);
        p = std::max(p, P_FLOOR);
        return std::sqrt(GAMMA * p / rho);
    }

    /**
     * Compute maximum wave speed in a given direction.
     * Used for CFL condition.
     */
    double max_wave_speed(const Eigen::VectorXd& U, int direction) const override {
        double rho, u, v, w, p;
        conservative_to_primitive(U, rho, u, v, w, p);

        double u_normal = (direction == 0) ? u : (direction == 1) ? v : w;
        double a = sound_speed(rho, p);

        return std::abs(u_normal) + a;
    }

    /**
     * Compute flux in X direction.
     * F = [rho*u, rho*u*u+p, rho*u*v, rho*u*w, (E+p)*u]
     */
    Eigen::VectorXd flux_x(const Eigen::VectorXd& U) const;

    /**
     * Compute flux in Y direction.
     * G = [rho*v, rho*v*u, rho*v*v+p, rho*v*w, (E+p)*v]
     */
    Eigen::VectorXd flux_y(const Eigen::VectorXd& U) const;

    /**
     * Compute flux in Z direction.
     * H = [rho*w, rho*w*u, rho*w*v, rho*w*w+p, (E+p)*w]
     */
    Eigen::VectorXd flux_z(const Eigen::VectorXd& U) const;
};

} // namespace fvm3d::physics
