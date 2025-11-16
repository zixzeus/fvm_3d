#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <functional>
#include "physics_base.hpp"

namespace fvm3d::physics {

/**
 * Advanced Resistive MagnetoHydroDynamic (MHD) equations for 3D magnetic reconnection.
 *
 * Based on OpenMHD 3D reconnection physics with:
 * - Position-dependent resistivity η(x,y,z)
 * - GLM divergence cleaning (∇·B = 0 constraint)
 * - Accurate Ohmic dissipation
 * - Harris sheet magnetic field configuration
 *
 * Conservative variables (9):
 *   U = [ρ, ρu, ρv, ρw, E, Bx, By, Bz, ψ]
 *   where ψ is the GLM divergence cleaning scalar
 *
 * Primitive variables (8):
 *   V = [ρ, u, v, w, p, Bx, By, Bz]
 *   Note: ψ is not a primitive variable (it's a numerical artifact)
 *
 * Key differences from basic MHD:
 * 1. Resistivity varies with position: η = η(x,y,z)
 * 2. GLM field ψ maintains ∇·B constraint
 * 3. Ohmic dissipation: η·J² appears in energy equation
 * 4. Magnetic field evolution includes resistive terms
 * 5. Harris equilibrium support for initial conditions
 */
class AdvancedResistiveMHD3D : public ConservationLaw {
public:
    static constexpr int nvars = 9;          // Including GLM field ψ
    static constexpr double GAMMA = 5.0/3.0; // Adiabatic index (monatomic gas)
    static constexpr double MU0 = 1.0;       // Permeability (normalized)
    static constexpr double RHO_FLOOR = 1e-10;
    static constexpr double P_FLOOR = 1e-11;
    static constexpr double B_FLOOR = 1e-12;

    /**
     * Resistivity model: two-level with gaussian localization
     * Used for magnetic reconnection:
     * η = η₀ + (η₁ - η₀) × sech²(√(x²+y²))
     *
     * Physical meaning:
     * - η₀: Background resistivity (weak, Rm₀ = 1/η₀ ~ 1000)
     * - η₁: Enhanced resistivity at origin (strong, Rm₁ = 1/η₁ ~ 60)
     * - Transition region: ~2-3 cell widths
     */
    struct ResistivityModel {
        double eta0;              // Background: 1e-3
        double eta1;              // Enhanced: ~1.67e-2
        double localization_scale; // Width of enhanced region (typically 1.0)

        ResistivityModel()
            : eta0(1.0 / 1000.0),
              eta1(1.0 / 60.0),
              localization_scale(1.0) {}

        /**
         * Compute resistivity at position (x, y, z).
         * For 2D reconnection, z-dependence is trivial.
         */
        double operator()(double x, double y, double z = 0.0) const {
            double r_sq = x*x + y*y;
            double r = std::sqrt(r_sq);
            // Ensure reasonable bounds
            r = std::min(r, 25.0);  // Cap to prevent overflow

            double sech_sq = 1.0 / std::cosh(r / localization_scale);
            sech_sq *= sech_sq;

            return eta0 + (eta1 - eta0) * sech_sq;
        }
    };

    /**
     * GLM divergence cleaning parameters.
     * Maintains ∇·B = 0 through:
     * ∂ψ/∂t = -ch·(∇·B) - (cr/ch)·ψ
     *
     * Parameters from Dedner et al. (2002):
     * - ch: Wave speed for hyperbolic part
     * - cr: Ratio of parabolic to hyperbolic dissipation
     */
    struct GLMParameters {
        double ch;  // Divergence wave speed (typically 0.2 × max wave speed)
        double cr;  // Parabolic dissipation ratio (typically 0.2)

        GLMParameters(double ch_val = 0.2, double cr_val = 0.2)
            : ch(ch_val), cr(cr_val) {}
    };

    /**
     * Harris sheet equilibrium configuration.
     * Initial condition for magnetic reconnection:
     * - Antiparallel magnetic field: Bx = tanh(y)
     * - Density variation: ρ = ρ₀ × sech²(y)
     * - Temperature balance: p = p₀ × (1 - (tanh²(y) - sech²(y))/β)
     */
    struct HarrisSheetConfig {
        double L_sheet;     // Thickness of current sheet
        double n0;          // Reference density
        double p0;          // Reference pressure
        double B0;          // Reference magnetic field
        double beta;        // Plasma beta = 2μ₀p₀/B₀²
        double perturbation_amplitude;  // Reconnection triggering (typically 0.03)

        HarrisSheetConfig()
            : L_sheet(1.0), n0(1.0), p0(0.1), B0(1.0),
              beta(0.2), perturbation_amplitude(0.03) {}
    };

    // ========== Constructors ==========

    /**
     * Constructor with resistivity and GLM models.
     */
    explicit AdvancedResistiveMHD3D(
        const ResistivityModel& res_model = ResistivityModel(),
        const GLMParameters& glm = GLMParameters()
    ) : ConservationLaw("AdvancedResistiveMHD3D", nvars, 3),
        resistivity_model_(res_model), glm_params_(glm) {}

    // ========== PhysicsBase Interface Implementation ==========

    /**
     * Convert conservative to primitive variables (PhysicsBase interface).
     * @param U: Conservative variables [ρ, ρu, ρv, ρw, E, Bx, By, Bz, ψ]
     * @return V: Primitive variables [ρ, u, v, w, p, Bx, By, Bz]
     * Note: ψ is dropped as it's not a primitive variable
     */
    Eigen::VectorXd conservative_to_primitive(const Eigen::VectorXd& U) const override {
        // Density with floor
        double rho = std::max(U(0), RHO_FLOOR);

        // Velocity components
        double u = U(1) / rho;
        double v = U(2) / rho;
        double w = U(3) / rho;

        // Magnetic field (directly from conservative)
        double Bx = U(5);
        double By = U(6);
        double Bz = U(7);

        // Compute energies
        double u_sq = u*u + v*v + w*w;
        double B_sq = Bx*Bx + By*By + Bz*Bz;
        double kinetic_energy = 0.5 * rho * u_sq;
        double magnetic_energy = 0.5 * B_sq / MU0;

        // Extract pressure from total energy
        double internal_energy = U(4) / rho - 0.5 * u_sq - (B_sq / (2.0 * MU0 * rho));
        double p = std::max((GAMMA - 1.0) * rho * internal_energy, P_FLOOR);

        Eigen::VectorXd V(8);
        V << rho, u, v, w, p, Bx, By, Bz;
        return V;
    }

    /**
     * Convert primitive to conservative variables (PhysicsBase interface).
     * @param V: Primitive variables [ρ, u, v, w, p, Bx, By, Bz]
     * @return U: Conservative variables [ρ, ρu, ρv, ρw, E, Bx, By, Bz, ψ]
     * Note: ψ is set to 0.0 by default
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
        double psi = 0.0;  // Default GLM field value

        Eigen::VectorXd U(nvars);
        U(0) = rho;
        U(1) = rho * u;
        U(2) = rho * v;
        U(3) = rho * w;
        U(5) = Bx;
        U(6) = By;
        U(7) = Bz;
        U(8) = psi;  // GLM field

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
     * Compute flux in given direction (PhysicsBase interface).
     * Note: Advanced MHD flux depends on position for resistivity.
     * This method uses position (0,0,0) by default.
     * For position-dependent flux, use flux_x/y/z directly.
     */
    Eigen::VectorXd compute_flux(const Eigen::VectorXd& U, int direction) const override {
        // Use default position (0,0,0) for interface compatibility
        switch (direction) {
            case 0: return flux_x(U, 0.0, 0.0, 0.0);
            case 1: return flux_y(U, 0.0, 0.0, 0.0);
            case 2: return flux_z(U, 0.0, 0.0, 0.0);
            default:
                throw std::invalid_argument(
                    "Invalid direction: " + std::to_string(direction) +
                    " (must be 0, 1, or 2)"
                );
        }
    }

    // ========== Wave Speed Calculations ==========

    /**
     * Sound speed (gas pressure wave).
     */
    double sound_speed(double rho, double p) const;

    /**
     * Alfvén speed along a direction.
     */
    double alfven_speed(double B_component, double rho) const;

    /**
     * Fast magnetosonic speed.
     */
    double fast_speed(double rho, double p, double Bx, double By, double Bz) const;

    /**
     * Maximum wave speed including GLM wave.
     */
    double max_wave_speed(const Eigen::VectorXd& U, int direction) const override;

    // ========== Flux Functions ==========

    /**
     * Compute flux in X direction (including GLM term).
     */
    Eigen::VectorXd flux_x(const Eigen::VectorXd& U, double x, double y, double z = 0.0) const;

    /**
     * Compute flux in Y direction.
     */
    Eigen::VectorXd flux_y(const Eigen::VectorXd& U, double x, double y, double z = 0.0) const;

    /**
     * Compute flux in Z direction.
     */
    Eigen::VectorXd flux_z(const Eigen::VectorXd& U, double x, double y, double z = 0.0) const;

    // ========== Source Terms ==========

    /**
     * Resistive dissipation source term.
     * Includes:
     * - Ohmic heating: η·J² in energy equation
     * - Magnetic field diffusion: ∇²B terms
     * - GLM constraint maintenance: -ch·(∇·B) - (cr/ch)·ψ
     *
     * @param dx, dy, dz: Grid spacing
     * @param U_center: State at cell center
     * @param neighbors: States at neighboring cell centers [x±, y±, z±]
     */
    Eigen::VectorXd resistive_source(
        const Eigen::VectorXd& U_center,
        double x, double y, double z,  // Position for η(x,y,z)
        double dx, double dy, double dz,
        const Eigen::VectorXd& U_xm, const Eigen::VectorXd& U_xp,
        const Eigen::VectorXd& U_ym, const Eigen::VectorXd& U_yp,
        const Eigen::VectorXd& U_zm, const Eigen::VectorXd& U_zp
    ) const;

    /**
     * GLM constraint source term.
     * dψ/dt = -ch·(∇·B) - (cr/ch)·ψ
     */
    double glm_source(
        const Eigen::VectorXd& U_center,
        double dx, double dy, double dz,
        const Eigen::VectorXd& U_xm, const Eigen::VectorXd& U_xp,
        const Eigen::VectorXd& U_ym, const Eigen::VectorXd& U_yp,
        const Eigen::VectorXd& U_zm, const Eigen::VectorXd& U_zp
    ) const;

    // ========== Initial Conditions ==========

    /**
     * Initialize Harris sheet equilibrium with perturbations.
     * Suitable for magnetic reconnection simulations.
     *
     * @param x, y, z: Position coordinates
     * @param harris: Harris sheet configuration
     * @return U: Initial conservative state
     */
    Eigen::VectorXd harris_sheet_initial(
        double x, double y, double z,
        const HarrisSheetConfig& harris = HarrisSheetConfig()
    ) const;

    /**
     * Initialize simple test: uniform field with small perturbation.
     */
    Eigen::VectorXd uniform_field_initial(
        double rho0, double p0, double Bx0, double By0, double Bz0,
        double perturbation = 0.01
    ) const;

    // ========== Accessors ==========

    const ResistivityModel& resistivity_model() const { return resistivity_model_; }
    void set_resistivity_model(const ResistivityModel& res) { resistivity_model_ = res; }

    const GLMParameters& glm_parameters() const { return glm_params_; }
    void set_glm_parameters(const GLMParameters& glm) { glm_params_ = glm; }

    // ========== Stability/Debug ==========

    /**
     * Check for unphysical states.
     * @return: true if state is valid (ρ > 0, p > 0)
     */
    bool is_valid_state(const Eigen::VectorXd& U) const;

    /**
     * Compute diagnostic quantities.
     * @return: [kinetic_energy, magnetic_energy, internal_energy, total_pressure]
     */
    Eigen::Vector4d compute_energies(const Eigen::VectorXd& U) const;

private:
    ResistivityModel resistivity_model_;
    GLMParameters glm_params_;

    /**
     * Compute current density J = ∇×B using finite differences.
     * @return: [Jx, Jy, Jz]
     */
    Eigen::Vector3d compute_current_density(
        double dx, double dy, double dz,
        const Eigen::VectorXd& U_xm, const Eigen::VectorXd& U_xp,
        const Eigen::VectorXd& U_ym, const Eigen::VectorXd& U_yp,
        const Eigen::VectorXd& U_zm, const Eigen::VectorXd& U_zp
    ) const;

    /**
     * Compute divergence of B using finite differences.
     */
    double compute_div_B(
        double dx, double dy, double dz,
        const Eigen::VectorXd& U_xm, const Eigen::VectorXd& U_xp,
        const Eigen::VectorXd& U_ym, const Eigen::VectorXd& U_yp,
        const Eigen::VectorXd& U_zm, const Eigen::VectorXd& U_zp
    ) const;

    /**
     * Compute laplacian of B components for resistive diffusion.
     */
    Eigen::Vector3d compute_laplacian_B(
        double dx, double dy, double dz,
        const Eigen::VectorXd& U_xm, const Eigen::VectorXd& U_xp,
        const Eigen::VectorXd& U_ym, const Eigen::VectorXd& U_yp,
        const Eigen::VectorXd& U_zm, const Eigen::VectorXd& U_zp
    ) const;
};

} // namespace fvm3d::physics
