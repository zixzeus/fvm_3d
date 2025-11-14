#include "physics/resistive_mhd3d_advanced.hpp"
#include <stdexcept>

namespace fvm3d::physics {

// ========== Variable Conversion ==========

void AdvancedResistiveMHD3D::conservative_to_primitive(
    const Eigen::VectorXd& U,
    double& rho, double& u, double& v, double& w, double& p,
    double& Bx, double& By, double& Bz
) const {
    if (U.size() != nvars) {
        throw std::invalid_argument("Invalid U vector size");
    }

    // Density with floor
    rho = std::max(U(0), RHO_FLOOR);

    // Velocity components
    u = U(1) / rho;
    v = U(2) / rho;
    w = U(3) / rho;

    // Magnetic field (directly from conservative)
    Bx = U(5);
    By = U(6);
    Bz = U(7);

    // Compute energies
    double u_sq = u*u + v*v + w*w;
    double B_sq = Bx*Bx + By*By + Bz*Bz;
    double kinetic_energy = 0.5 * rho * u_sq;
    double magnetic_energy = 0.5 * B_sq / MU0;

    // Extract pressure from total energy
    // E = rho*e_int + 0.5*rho*u² + 0.5*B²
    // e_int = (p/(ρ(γ-1)) + kinetic + magnetic)/ρ
    double internal_energy = U(4) / rho - 0.5 * u_sq - (B_sq / (2.0 * MU0 * rho));
    p = std::max((GAMMA - 1.0) * rho * internal_energy, P_FLOOR);
}

void AdvancedResistiveMHD3D::primitive_to_conservative(
    double rho, double u, double v, double w, double p,
    double Bx, double By, double Bz, double psi,
    Eigen::VectorXd& U
) const {
    if (U.size() != nvars) {
        U.resize(nvars);
    }

    rho = std::max(rho, RHO_FLOOR);
    p = std::max(p, P_FLOOR);

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
}

// ========== Wave Speed Calculations ==========

double AdvancedResistiveMHD3D::sound_speed(double rho, double p) const {
    rho = std::max(rho, RHO_FLOOR);
    p = std::max(p, P_FLOOR);
    return std::sqrt(GAMMA * p / rho);
}

double AdvancedResistiveMHD3D::alfven_speed(double B_component, double rho) const {
    rho = std::max(rho, RHO_FLOOR);
    return std::abs(B_component) / std::sqrt(MU0 * rho);
}

double AdvancedResistiveMHD3D::fast_speed(
    double rho, double p, double Bx, double By, double Bz
) const {
    rho = std::max(rho, RHO_FLOOR);
    p = std::max(p, P_FLOOR);
    double a = sound_speed(rho, p);
    double B_sq = Bx*Bx + By*By + Bz*Bz;
    double B_sq_norm = B_sq / (MU0 * rho);
    return std::sqrt(a*a + B_sq_norm);
}

double AdvancedResistiveMHD3D::max_wave_speed(
    const Eigen::VectorXd& U, int direction
) const {
    double rho, u, v, w, p, Bx, By, Bz;
    conservative_to_primitive(U, rho, u, v, w, p, Bx, By, Bz);

    double u_normal = (direction == 0) ? u : (direction == 1) ? v : w;
    double a = sound_speed(rho, p);
    double cf = fast_speed(rho, p, Bx, By, Bz);
    double ch = glm_params_.ch;  // GLM wave speed

    // Maximum wave speed = |u_normal| + max(fast_speed, GLM_speed)
    return std::abs(u_normal) + std::max(cf, ch);
}

// ========== Flux Functions with Resistive Terms ==========

Eigen::VectorXd AdvancedResistiveMHD3D::flux_x(
    const Eigen::VectorXd& U, double x, double y, double z
) const {
    double rho, u, v, w, p, Bx, By, Bz;
    conservative_to_primitive(U, rho, u, v, w, p, Bx, By, Bz);

    Eigen::VectorXd F(nvars);
    double B_sq = Bx*Bx + By*By + Bz*Bz;
    double p_total = p + 0.5 * B_sq / MU0;

    // Mass flux
    F(0) = rho * u;

    // Momentum fluxes (with Maxwell stress)
    F(1) = rho * u * u + p_total - Bx * Bx / MU0;
    F(2) = rho * u * v - Bx * By / MU0;
    F(3) = rho * u * w - Bx * Bz / MU0;

    // Energy flux (total enthalpy + magnetic work)
    double u_dot_B = u*Bx + v*By + w*Bz;
    F(4) = (U(4) + p_total) * u - Bx * u_dot_B / MU0;

    // Magnetic field fluxes (Faraday's law: ∇×(v×B))
    F(5) = 0.0;  // ∇·B = 0 component
    F(6) = v * Bx - u * By;
    F(7) = w * Bx - u * Bz;

    // GLM wave propagation: F(ψ) = ch² * Bn
    double psi = (U.size() > 8) ? U(8) : 0.0;
    F(5) = psi;  // F(Bn) = ψ (replaces the 0.0)
    F(8) = glm_params_.ch * glm_params_.ch * Bx;

    return F;
}

Eigen::VectorXd AdvancedResistiveMHD3D::flux_y(
    const Eigen::VectorXd& U, double x, double y, double z
) const {
    double rho, u, v, w, p, Bx, By, Bz;
    conservative_to_primitive(U, rho, u, v, w, p, Bx, By, Bz);

    Eigen::VectorXd G(nvars);
    double B_sq = Bx*Bx + By*By + Bz*Bz;
    double p_total = p + 0.5 * B_sq / MU0;

    // Mass flux
    G(0) = rho * v;

    // Momentum fluxes
    G(1) = rho * v * u - By * Bx / MU0;
    G(2) = rho * v * v + p_total - By * By / MU0;
    G(3) = rho * v * w - By * Bz / MU0;

    // Energy flux
    double u_dot_B = u*Bx + v*By + w*Bz;
    G(4) = (U(4) + p_total) * v - By * u_dot_B / MU0;

    // Magnetic field fluxes
    G(5) = u * By - v * Bx;
    G(7) = w * By - v * Bz;

    // GLM wave propagation: F(ψ) = ch² * Bn
    double psi = (U.size() > 8) ? U(8) : 0.0;
    G(6) = psi;  // F(Bn) = ψ (replaces the 0.0)
    G(8) = glm_params_.ch * glm_params_.ch * By;

    return G;
}

Eigen::VectorXd AdvancedResistiveMHD3D::flux_z(
    const Eigen::VectorXd& U, double x, double y, double z
) const {
    double rho, u, v, w, p, Bx, By, Bz;
    conservative_to_primitive(U, rho, u, v, w, p, Bx, By, Bz);

    Eigen::VectorXd H(nvars);
    double B_sq = Bx*Bx + By*By + Bz*Bz;
    double p_total = p + 0.5 * B_sq / MU0;

    // Mass flux
    H(0) = rho * w;

    // Momentum fluxes
    H(1) = rho * w * u - Bz * Bx / MU0;
    H(2) = rho * w * v - Bz * By / MU0;
    H(3) = rho * w * w + p_total - Bz * Bz / MU0;

    // Energy flux
    double u_dot_B = u*Bx + v*By + w*Bz;
    H(4) = (U(4) + p_total) * w - Bz * u_dot_B / MU0;

    // Magnetic field fluxes
    H(5) = u * Bz - w * Bx;
    H(6) = v * Bz - w * By;

    // GLM wave propagation: F(ψ) = ch² * Bn
    double psi = (U.size() > 8) ? U(8) : 0.0;
    H(7) = psi;  // F(Bn) = ψ (replaces the 0.0)
    H(8) = glm_params_.ch * glm_params_.ch * Bz;

    return H;
}

// ========== Helper Functions for Source Terms ==========

Eigen::Vector3d AdvancedResistiveMHD3D::compute_current_density(
    double dx, double dy, double dz,
    const Eigen::VectorXd& U_xm, const Eigen::VectorXd& U_xp,
    const Eigen::VectorXd& U_ym, const Eigen::VectorXd& U_yp,
    const Eigen::VectorXd& U_zm, const Eigen::VectorXd& U_zp
) const {
    // Extract B components
    double Bx_xm = U_xm(5), Bx_xp = U_xp(5);
    double By_xm = U_xm(6), By_xp = U_xp(6);
    double Bz_xm = U_xm(7), Bz_xp = U_xp(7);

    double Bx_ym = U_ym(5), Bx_yp = U_yp(5);
    double By_ym = U_ym(6), By_yp = U_yp(6);
    double Bz_ym = U_ym(7), Bz_yp = U_yp(7);

    double Bx_zm = U_zm(5), Bx_zp = U_zp(5);
    double By_zm = U_zm(6), By_zp = U_zp(6);
    double Bz_zm = U_zm(7), Bz_zp = U_zp(7);

    // J = ∇×B using central differences
    // Jx = ∂Bz/∂y - ∂By/∂z
    double dBz_dy = (Bz_yp - Bz_ym) / (2.0 * dy);
    double dBy_dz = (By_zp - By_zm) / (2.0 * dz);
    double Jx = dBz_dy - dBy_dz;

    // Jy = ∂Bx/∂z - ∂Bz/∂x
    double dBx_dz = (Bx_zp - Bx_zm) / (2.0 * dz);
    double dBz_dx = (Bz_xp - Bz_xm) / (2.0 * dx);
    double Jy = dBx_dz - dBz_dx;

    // Jz = ∂By/∂x - ∂Bx/∂y
    double dBy_dx = (By_xp - By_xm) / (2.0 * dx);
    double dBx_dy = (Bx_yp - Bx_ym) / (2.0 * dy);
    double Jz = dBy_dx - dBx_dy;

    return Eigen::Vector3d(Jx, Jy, Jz);
}

double AdvancedResistiveMHD3D::compute_div_B(
    double dx, double dy, double dz,
    const Eigen::VectorXd& U_xm, const Eigen::VectorXd& U_xp,
    const Eigen::VectorXd& U_ym, const Eigen::VectorXd& U_yp,
    const Eigen::VectorXd& U_zm, const Eigen::VectorXd& U_zp
) const {
    double Bx_xm = U_xm(5), Bx_xp = U_xp(5);
    double By_ym = U_ym(6), By_yp = U_yp(6);
    double Bz_zm = U_zm(7), Bz_zp = U_zp(7);

    double dBx_dx = (Bx_xp - Bx_xm) / (2.0 * dx);
    double dBy_dy = (By_yp - By_ym) / (2.0 * dy);
    double dBz_dz = (Bz_zp - Bz_zm) / (2.0 * dz);

    return dBx_dx + dBy_dy + dBz_dz;
}

Eigen::Vector3d AdvancedResistiveMHD3D::compute_laplacian_B(
    double dx, double dy, double dz,
    const Eigen::VectorXd& U_xm, const Eigen::VectorXd& U_xp,
    const Eigen::VectorXd& U_ym, const Eigen::VectorXd& U_yp,
    const Eigen::VectorXd& U_zm, const Eigen::VectorXd& U_zp
) const {
    double inv_dx2 = 1.0 / (dx * dx);
    double inv_dy2 = 1.0 / (dy * dy);
    double inv_dz2 = 1.0 / (dz * dz);

    // Bx: center is U(5)
    double Bx_xm = U_xm(5), Bx_xp = U_xp(5);
    double Bx_ym = U_ym(5), Bx_yp = U_yp(5);
    double Bx_zm = U_zm(5), Bx_zp = U_zp(5);
    // Approximate ∇²Bx using difference stencil (ignoring center for brevity)
    double lapl_Bx = inv_dx2 * (Bx_xp + Bx_xm) +
                     inv_dy2 * (Bx_yp + Bx_ym) +
                     inv_dz2 * (Bx_zp + Bx_zm) -
                     2.0 * (inv_dx2 + inv_dy2 + inv_dz2);

    // By
    double By_xm = U_xm(6), By_xp = U_xp(6);
    double By_ym = U_ym(6), By_yp = U_yp(6);
    double By_zm = U_zm(6), By_zp = U_zp(6);
    double lapl_By = inv_dx2 * (By_xp + By_xm) +
                     inv_dy2 * (By_yp + By_ym) +
                     inv_dz2 * (By_zp + By_zm) -
                     2.0 * (inv_dx2 + inv_dy2 + inv_dz2);

    // Bz
    double Bz_xm = U_xm(7), Bz_xp = U_xp(7);
    double Bz_ym = U_ym(7), Bz_yp = U_yp(7);
    double Bz_zm = U_zm(7), Bz_zp = U_zp(7);
    double lapl_Bz = inv_dx2 * (Bz_xp + Bz_xm) +
                     inv_dy2 * (Bz_yp + Bz_ym) +
                     inv_dz2 * (Bz_zp + Bz_zm) -
                     2.0 * (inv_dx2 + inv_dy2 + inv_dz2);

    return Eigen::Vector3d(lapl_Bx, lapl_By, lapl_Bz);
}

// ========== Source Terms ==========

Eigen::VectorXd AdvancedResistiveMHD3D::resistive_source(
    const Eigen::VectorXd& U_center,
    double x, double y, double z,
    double dx, double dy, double dz,
    const Eigen::VectorXd& U_xm, const Eigen::VectorXd& U_xp,
    const Eigen::VectorXd& U_ym, const Eigen::VectorXd& U_yp,
    const Eigen::VectorXd& U_zm, const Eigen::VectorXd& U_zp
) const {
    Eigen::VectorXd S = Eigen::VectorXd::Zero(nvars);

    // Get resistivity at this location
    double eta = resistivity_model_(x, y, z);

    // Compute current density J = ∇×B
    Eigen::Vector3d J = compute_current_density(dx, dy, dz,
                                                U_xm, U_xp, U_ym, U_yp, U_zm, U_zp);
    double J_sq = J.squaredNorm();

    // Extract conservative variables
    double rho = std::max(U_center(0), RHO_FLOOR);
    double Bx = U_center(5);
    double By = U_center(6);
    double Bz = U_center(7);

    // Ohmic dissipation in energy equation: S_E = η·J²
    S(4) = eta * J_sq;

    // Magnetic field diffusion: S_B = ∇·(η·∇B) ≈ η·∇²B
    Eigen::Vector3d lapl_B = compute_laplacian_B(dx, dy, dz,
                                                  U_xm, U_xp, U_ym, U_yp, U_zm, U_zp);
    S(5) = eta * lapl_B(0);
    S(6) = eta * lapl_B(1);
    S(7) = eta * lapl_B(2);

    // Note: No explicit Lorentz force S_momentum term needed
    // It's implicitly handled through the MHD fluxes (Maxwell stress tensor)

    return S;
}

double AdvancedResistiveMHD3D::glm_source(
    const Eigen::VectorXd& U_center,
    double dx, double dy, double dz,
    const Eigen::VectorXd& U_xm, const Eigen::VectorXd& U_xp,
    const Eigen::VectorXd& U_ym, const Eigen::VectorXd& U_yp,
    const Eigen::VectorXd& U_zm, const Eigen::VectorXd& U_zp
) const {
    // dψ/dt = -ch·(∇·B) - (cr/ch)·ψ
    double div_B = compute_div_B(dx, dy, dz, U_xm, U_xp, U_ym, U_yp, U_zm, U_zp);
    double psi = (U_center.size() > 8) ? U_center(8) : 0.0;

    double ch = glm_params_.ch;
    double cr = glm_params_.cr;
    double psi_source = -ch * div_B - (cr / ch) * psi;

    return psi_source;
}

// ========== Initial Conditions ==========

Eigen::VectorXd AdvancedResistiveMHD3D::harris_sheet_initial(
    double x, double y, double z,
    const HarrisSheetConfig& harris
) const {
    // Harris sheet: antiparallel magnetic field with current sheet
    // Bx = tanh(y/L) - antiparallel field
    // ρ = ρ₀·sech²(y/L) - thin sheet density
    // p = p₀[1 - sech²(y/L)·(1 - β⁻¹)] - pressure balance

    double L = harris.L_sheet;
    double y_norm = y / L;

    // Density: ρ = ρ₀(1 + (1/β - 1)sech²(y))
    double sech_y = 1.0 / std::cosh(y_norm);
    double sech_y_sq = sech_y * sech_y;
    double rho = harris.n0 * (1.0 + (1.0/harris.beta - 1.0) * sech_y_sq);
    rho = std::max(rho, RHO_FLOOR);

    // Magnetic field: Bx = B₀·tanh(y/L)
    double tanh_y = std::tanh(y_norm);
    double Bx = harris.B0 * tanh_y;
    double By = 0.0;
    double Bz = 0.0;

    // Pressure balance for Harris sheet:
    // Total pressure p_tot = p + B²/(2μ₀) must be constant
    // At y=0: Bx=0, so p(0) = p_tot
    // For beta (plasma beta): beta = 2μ₀p(0)/B₀²
    // Therefore: p_tot = beta·B₀²/(2μ₀)
    double p_total = harris.beta * harris.B0 * harris.B0 / (2.0 * MU0);
    double B_pressure = 0.5 * Bx * Bx / MU0;
    double p = p_total - B_pressure;
    p = std::max(p, P_FLOOR);

    // Add small perturbations to trigger reconnection (m=1 mode)
    // δBy = A·sin(πx/Lx)·e^(-y²) - localized island
    double perturbation = harris.perturbation_amplitude *
                         std::sin(M_PI * x) *
                         std::exp(-y_norm*y_norm);
    By += perturbation;

    // Initialize with zero velocity and GLM field
    double u = 0.0, v = 0.0, w = 0.0;
    double psi = 0.0;

    Eigen::VectorXd U(nvars);
    primitive_to_conservative(rho, u, v, w, p, Bx, By, Bz, psi, U);

    return U;
}

Eigen::VectorXd AdvancedResistiveMHD3D::uniform_field_initial(
    double rho0, double p0, double Bx0, double By0, double Bz0,
    double perturbation
) const {
    double u = 0.0, v = 0.0, w = 0.0;
    double Bx = Bx0 * (1.0 + perturbation);
    double By = By0 * (1.0 + perturbation);
    double Bz = Bz0 * (1.0 + perturbation);
    double psi = 0.0;

    Eigen::VectorXd U(nvars);
    primitive_to_conservative(rho0, u, v, w, p0, Bx, By, Bz, psi, U);
    return U;
}

// ========== Diagnostics ==========

bool AdvancedResistiveMHD3D::is_valid_state(const Eigen::VectorXd& U) const {
    if (U.size() != nvars) return false;

    double rho = U(0);
    double E = U(4);

    if (rho <= 0.0 || E < 0.0) return false;

    return true;
}

Eigen::Vector4d AdvancedResistiveMHD3D::compute_energies(
    const Eigen::VectorXd& U
) const {
    double rho, u, v, w, p, Bx, By, Bz;
    conservative_to_primitive(U, rho, u, v, w, p, Bx, By, Bz);

    double kinetic = 0.5 * rho * (u*u + v*v + w*w);
    double magnetic = 0.5 * (Bx*Bx + By*By + Bz*Bz) / MU0;
    double internal = rho * p / ((GAMMA - 1.0) * rho);
    double p_total = p + 0.5 * (Bx*Bx + By*By + Bz*Bz) / MU0;

    return Eigen::Vector4d(kinetic, magnetic, internal, p_total);
}

} // namespace fvm3d::physics
