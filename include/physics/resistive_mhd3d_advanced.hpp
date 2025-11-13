#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <functional>

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
 *
 * Key differences from basic MHD:
 * 1. Resistivity varies with position: η = η(x,y,z)
 * 2. GLM field ψ maintains ∇·B constraint
 * 3. Ohmic dissipation: η·J² appears in energy equation
 * 4. Magnetic field evolution includes resistive terms
 * 5. Harris equilibrium support for initial conditions
 */
class AdvancedResistiveMHD3D {
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
    ) : resistivity_model_(res_model), glm_params_(glm) {}

    // ========== Variable Conversion ==========

    /**
     * Convert conservative variables to primitive variables.
     * @param U: conservative state [ρ, ρu, ρv, ρw, E, Bx, By, Bz, ψ]
     * @param rho, u, v, w, p, Bx, By, Bz: primitive variables (output)
     */
    void conservative_to_primitive(
        const Eigen::VectorXd& U,
        double& rho, double& u, double& v, double& w, double& p,
        double& Bx, double& By, double& Bz
    ) const;

    /**
     * Convert primitive variables to conservative variables.
     * @param psi: GLM divergence cleaning scalar
     */
    void primitive_to_conservative(
        double rho, double u, double v, double w, double p,
        double Bx, double By, double Bz, double psi,
        Eigen::VectorXd& U
    ) const;

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
    double max_wave_speed(const Eigen::VectorXd& U, int direction) const;

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
