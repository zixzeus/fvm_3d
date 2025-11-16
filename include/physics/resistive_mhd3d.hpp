#pragma once

#include "physics/physics_base.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>

namespace fvm3d::physics {

/**
 * Resistive MagnetoHydroDynamic (MHD) equations for 3D.
 *
 * Conservative variables (8):
 *   U = [rho, rho_u, rho_v, rho_w, E, Bx, By, Bz]
 *
 * Primitive variables (8):
 *   V = [rho, u, v, w, p, Bx, By, Bz]
 *
 * Where:
 *   rho  : density
 *   u,v,w: velocity components
 *   p    : gas pressure (not total pressure)
 *   B    : magnetic field components
 *   E    : total energy density (kinetic + internal + magnetic)
 *
 * Total energy:
 *   E = rho*e + 0.5*rho*|u|^2 + 0.5*|B|^2
 *   where e is specific internal energy
 *
 * Total pressure (for wave speeds):
 *   p_total = p + 0.5*|B|^2  (magnetic + gas pressure)
 *
 * Physical constants (SI units):
 *   gamma      : adiabatic index (typically 1.4 for air)
 *   mu0        : permeability (set to 1.0 in normalized units)
 *   eta        : resistivity coefficient (0 for ideal MHD)
 */
class ResistiveMHD3D : public ConservationLaw {
public:
    static constexpr int nvars = 8;
    static constexpr double GAMMA = 1.4;        // Adiabatic index
    static constexpr double MU0 = 1.0;          // Permeability (normalized)
    static constexpr double RHO_FLOOR = 1e-10;  // Density floor
    static constexpr double P_FLOOR = 1e-11;    // Pressure floor
    static constexpr double B_FLOOR = 1e-12;    // Magnetic field floor

    /**
     * Constructor with optional resistivity.
     * @param resistivity: Ohmic resistivity eta (0 for ideal MHD)
     */
    explicit ResistiveMHD3D(double resistivity = 0.0)
        : ConservationLaw("ResistiveMHD3D", nvars, 3), eta_(resistivity) {}

    /**
     * Convert conservative to primitive variables (PhysicsBase interface).
     * @param U: Conservative variables [rho, rho_u, rho_v, rho_w, E, Bx, By, Bz]
     * @return V: Primitive variables [rho, u, v, w, p, Bx, By, Bz]
     */
    Eigen::VectorXd conservative_to_primitive(const Eigen::VectorXd& U) const override {
        // Density with floor
        double rho = std::max(U(0), RHO_FLOOR);

        // Velocity
        double u = U(1) / rho;
        double v = U(2) / rho;
        double w = U(3) / rho;

        // Magnetic field (directly from conservative)
        double Bx = U(5);
        double By = U(6);
        double Bz = U(7);

        // Compute kinetic and magnetic energy
        double u_sq = u*u + v*v + w*w;
        double B_sq = Bx*Bx + By*By + Bz*Bz;
        double kinetic_energy = 0.5 * rho * u_sq;
        double magnetic_energy = 0.5 * B_sq / MU0;

        // Internal energy
        double internal_energy = U(4) / rho - 0.5 * u_sq - (B_sq / (2.0 * MU0 * rho));
        double p = std::max((GAMMA - 1.0) * rho * internal_energy, P_FLOOR);

        Eigen::VectorXd V(nvars);
        V << rho, u, v, w, p, Bx, By, Bz;
        return V;
    }

    /**
     * Convert primitive to conservative variables (PhysicsBase interface).
     * @param V: Primitive variables [rho, u, v, w, p, Bx, By, Bz]
     * @return U: Conservative variables [rho, rho_u, rho_v, rho_w, E, Bx, By, Bz]
     */
    Eigen::VectorXd primitive_to_conservative(const Eigen::VectorXd& V) const override {
        double rho = std::max(V(0), RHO_FLOOR);
        double u = V(1);
        double v = V(2);
        double w = V(3);
        double p = std::max(V(4), P_FLOOR);
        double Bx = V(5);
        double By = V(6);
        double Bz = V(7);

        Eigen::VectorXd U(nvars);
        U(0) = rho;
        U(1) = rho * u;
        U(2) = rho * v;
        U(3) = rho * w;
        U(5) = Bx;
        U(6) = By;
        U(7) = Bz;

        // Total energy density
        double u_sq = u*u + v*v + w*w;
        double B_sq = Bx*Bx + By*By + Bz*Bz;
        double kinetic_energy = 0.5 * rho * u_sq;
        double internal_energy = p / ((GAMMA - 1.0) * rho);
        double magnetic_energy = 0.5 * B_sq / MU0;

        U(4) = rho * internal_energy + kinetic_energy + magnetic_energy;
        return U;
    }

    /**
     * Compute sound speed (gas sound speed, not Alfvén).
     * a = sqrt(gamma * p / rho)
     */
    double sound_speed(double rho, double p) const {
        rho = std::max(rho, RHO_FLOOR);
        p = std::max(p, P_FLOOR);
        return std::sqrt(GAMMA * p / rho);
    }

    /**
     * Compute Alfvén speed in a given direction.
     * va = |B_dir| / sqrt(mu0 * rho)
     */
    double alfven_speed(double B_component, double rho) const {
        rho = std::max(rho, RHO_FLOOR);
        return std::abs(B_component) / std::sqrt(MU0 * rho);
    }

    /**
     * Compute fast magnetosonic speed.
     * cf = sqrt(a^2 + (B_x^2 + B_y^2 + B_z^2)/(mu0*rho))
     */
    double fast_speed(double rho, double p, double Bx, double By, double Bz) const {
        rho = std::max(rho, RHO_FLOOR);
        p = std::max(p, P_FLOOR);
        double a = sound_speed(rho, p);
        double B_sq = Bx*Bx + By*By + Bz*Bz;
        double B_sq_norm = B_sq / (MU0 * rho);
        return std::sqrt(a*a + B_sq_norm);
    }

    /**
     * Compute maximum wave speed in a given direction.
     * For MHD, this includes acoustic waves and fast magnetosonic waves.
     * @param direction: 0=X, 1=Y, 2=Z
     */
    double max_wave_speed(const Eigen::VectorXd& U, int direction) const override {
        Eigen::VectorXd V = conservative_to_primitive(U);
        double rho = V(0);
        double u = V(1);
        double v = V(2);
        double w = V(3);
        double p = V(4);
        double Bx = V(5);
        double By = V(6);
        double Bz = V(7);

        double u_normal = (direction == 0) ? u : (direction == 1) ? v : w;
        double a = sound_speed(rho, p);
        double cf = fast_speed(rho, p, Bx, By, Bz);

        // Maximum wave speed = |u_normal| + fast_speed
        return std::abs(u_normal) + cf;
    }

    /**
     * Compute flux in X direction.
     * F = [rho*u, rho*u*u+p_tot-Bx^2/mu0, rho*u*v-Bx*By/mu0, rho*u*w-Bx*Bz/mu0,
     *      (E+p_tot)*u - Bx*(u*Bx+v*By+w*Bz)/mu0, 0, v*Bx-u*By, w*Bx-u*Bz]
     */
    Eigen::VectorXd flux_x(const Eigen::VectorXd& U) const;

    /**
     * Compute flux in Y direction.
     */
    Eigen::VectorXd flux_y(const Eigen::VectorXd& U) const;

    /**
     * Compute flux in Z direction.
     */
    Eigen::VectorXd flux_z(const Eigen::VectorXd& U) const;

    // ========== PhysicsBase Interface Implementation ==========

    /**
     * Unified flux computation (required by PhysicsBase).
     * Routes to appropriate directional flux method.
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

    /**
     * Compute resistive dissipation source term.
     * For resistive MHD: S_B = eta * J where J = curl(B)
     * Requires calculation of current density from magnetic field gradients.
     */
    Eigen::VectorXd resistive_source(
        const Eigen::VectorXd& U,
        double dx, double dy, double dz,
        // Neighboring states for gradient computation
        const Eigen::VectorXd& U_xm, const Eigen::VectorXd& U_xp,
        const Eigen::VectorXd& U_ym, const Eigen::VectorXd& U_yp,
        const Eigen::VectorXd& U_zm, const Eigen::VectorXd& U_zp
    ) const;

    /**
     * Get resistivity coefficient.
     */
    double resistivity() const { return eta_; }

    /**
     * Set resistivity coefficient.
     */
    void set_resistivity(double eta) { eta_ = eta; }

private:
    double eta_;  // Resistivity coefficient
};

} // namespace fvm3d::physics
