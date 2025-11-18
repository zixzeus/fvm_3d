#include "spatial/riemann_solvers/riemann_hlld.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

HLLDSolver::HLLDSolver(const std::shared_ptr<physics::AdvancedResistiveMHD3D>& mhd)
    : RiemannSolver(mhd), mhd_(mhd) {
    if (!mhd_) {
        throw std::invalid_argument("MHD physics object cannot be null");
    }
}

Eigen::VectorXd HLLDSolver::solve(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    int nvars = U_L.size();

    // Rotate to normal direction (align with x-axis)
    Eigen::VectorXd UL_rot = rotate_to_normal(U_L, direction);
    Eigen::VectorXd UR_rot = rotate_to_normal(U_R, direction);

    // Extract rotated primitive variables
    Eigen::VectorXd VL_rot = mhd_->conservative_to_primitive(UL_rot);
    Eigen::VectorXd VR_rot = mhd_->conservative_to_primitive(UR_rot);

    double rho_L = VL_rot(0);
    double u_L = VL_rot(1);
    double v_L = VL_rot(2);
    double w_L = VL_rot(3);
    double p_L = VL_rot(4);
    double B_x = VL_rot(5);  // Normal magnetic field (same for left and right)
    double By_L = VL_rot(6);
    double Bz_L = VL_rot(7);

    double rho_R = VR_rot(0);
    double u_R = VR_rot(1);
    double v_R = VR_rot(2);
    double w_R = VR_rot(3);
    double p_R = VR_rot(4);
    double By_R = VR_rot(6);
    double Bz_R = VR_rot(7);

    // Compute sound speeds and fast magnetosonic speeds
    double a_L = mhd_->sound_speed(rho_L, p_L);
    double a_R = mhd_->sound_speed(rho_R, p_R);
    double cf_L = mhd_->fast_speed(rho_L, p_L, B_x, By_L, Bz_L);
    double cf_R = mhd_->fast_speed(rho_R, p_R, B_x, By_R, Bz_R);

    // Step 1: Estimate outer wave speeds (HLL-type)
    double S_L = std::min(u_L - cf_L, u_R - cf_R);
    double S_R = std::max(u_L + cf_L, u_R + cf_R);

    // Safety check
    if (S_R <= S_L) {
        Eigen::VectorXd F_L = mhd_->flux_x(UL_rot, 0.0, 0.0, 0.0);
        Eigen::VectorXd F_R = mhd_->flux_x(UR_rot, 0.0, 0.0, 0.0);
        return rotate_from_normal(0.5 * (F_L + F_R), direction);
    }

    // Compute original fluxes
    Eigen::VectorXd F_L = mhd_->flux_x(UL_rot, 0.0, 0.0, 0.0);
    Eigen::VectorXd F_R = mhd_->flux_x(UR_rot, 0.0, 0.0, 0.0);

    // Step 2: Estimate intermediate pressure
    double p_m = estimate_p_middle(VL_rot, VR_rot);

    // Step 3: Estimate contact discontinuity speed
    double S_M = estimate_s_contact(rho_L, u_L, p_L, S_L, rho_R, u_R, p_R, S_R);

    // Step 4: Compute intermediate states U_*L and U_*R
    Eigen::VectorXd U_Lstar = compute_state_L(UL_rot, VL_rot, S_L, S_M, p_m, B_x);
    Eigen::VectorXd U_Rstar = compute_state_R(UR_rot, VR_rot, S_R, S_M, p_m, B_x);

    // Step 5: Select flux based on location relative to wave structure
    Eigen::VectorXd F_hlld(nvars);

    if (0.0 <= S_L) {
        // Left of left fast wave
        F_hlld = F_L;
    } else if (0.0 >= S_R) {
        // Right of right fast wave
        F_hlld = F_R;
    } else if (0.0 <= S_M) {
        // Between left fast wave and contact discontinuity
        Eigen::VectorXd F_Lstar = F_L + S_L * (U_Lstar - UL_rot);
        F_hlld = F_Lstar;
    } else {
        // Between contact discontinuity and right fast wave
        Eigen::VectorXd F_Rstar = F_R + S_R * (U_Rstar - UR_rot);
        F_hlld = F_Rstar;
    }

    // Rotate flux back to original direction
    return rotate_from_normal(F_hlld, direction);
}

Eigen::VectorXd HLLDSolver::rotate_to_normal(const Eigen::VectorXd& U, int direction) const {
    Eigen::VectorXd U_rot = U;

    if (direction == 0) {
        // X is normal: no rotation
        return U_rot;
    } else if (direction == 1) {
        // Y is normal: rotate x->y, y->x
        std::swap(U_rot(1), U_rot(2));  // u <-> v
        std::swap(U_rot(5), U_rot(6));  // Bx <-> By
    } else if (direction == 2) {
        // Z is normal: rotate x->z, z->x
        std::swap(U_rot(1), U_rot(3));  // u <-> w
        std::swap(U_rot(5), U_rot(7));  // Bx <-> Bz
    }

    return U_rot;
}

Eigen::VectorXd HLLDSolver::rotate_from_normal(const Eigen::VectorXd& U, int direction) const {
    Eigen::VectorXd U_orig = U;

    if (direction == 0) {
        // X was normal: no rotation
        return U_orig;
    } else if (direction == 1) {
        // Y was normal: rotate back x->y, y->x
        std::swap(U_orig(1), U_orig(2));
        std::swap(U_orig(5), U_orig(6));
    } else if (direction == 2) {
        // Z was normal: rotate back x->z, z->x
        std::swap(U_orig(1), U_orig(3));
        std::swap(U_orig(5), U_orig(7));
    }

    return U_orig;
}

double HLLDSolver::estimate_s_contact(
    double rho_L, double u_L, double p_L, double S_L,
    double rho_R, double u_R, double p_R, double S_R
) const {
    // Contact discontinuity speed using weighted average
    // S_M = (ρ_L * u_L * (S_L - u_L) + ρ_R * u_R * (S_R - u_R) + (p_L - p_R))
    //       / (ρ_L * (S_L - u_L) + ρ_R * (S_R - u_R))

    double denom = rho_L * (S_L - u_L) + rho_R * (S_R - u_R);
    if (std::abs(denom) < 1.0e-12) {
        // Fallback to Roe average
        double sqrt_rho_L = std::sqrt(rho_L);
        double sqrt_rho_R = std::sqrt(rho_R);
        return (sqrt_rho_L * u_L + sqrt_rho_R * u_R) / (sqrt_rho_L + sqrt_rho_R);
    }

    double numer = rho_L * u_L * (S_L - u_L) + rho_R * u_R * (S_R - u_R) + (p_L - p_R);
    return numer / denom;
}

double HLLDSolver::estimate_p_middle(
    const Eigen::VectorXd& V_L, const Eigen::VectorXd& V_R
) const {
    // Roe-type averaging for middle pressure
    // p_m = (a_R * p_L + a_L * p_R + a_L * a_R * (u_L - u_R)) / (a_L + a_R)

    double rho_L = V_L(0);
    double u_L = V_L(1);
    double p_L = V_L(4);

    double rho_R = V_R(0);
    double u_R = V_R(1);
    double p_R = V_R(4);

    double a_L = mhd_->sound_speed(rho_L, p_L);
    double a_R = mhd_->sound_speed(rho_R, p_R);

    double denom = a_L + a_R;
    if (std::abs(denom) < 1.0e-12) {
        return 0.5 * (p_L + p_R);
    }

    double numer = a_R * p_L + a_L * p_R + a_L * a_R * (u_L - u_R);
    double p_m = numer / denom;

    // Apply floor
    return std::max(p_m, physics_->p_floor());
}

Eigen::VectorXd HLLDSolver::compute_state_L(
    const Eigen::VectorXd& U_L, const Eigen::VectorXd& V_L,
    double S_L, double S_M, double p_m, double B_x
) const {
    int nvars = U_L.size();
    Eigen::VectorXd U_star = U_L;

    double rho_L = V_L(0);
    double u_L = V_L(1);
    double v_L = V_L(2);
    double w_L = V_L(3);
    double p_L = V_L(4);
    double By_L = V_L(6);
    double Bz_L = V_L(7);

    // Density in intermediate state
    double rho_Lstar = rho_L * (S_L - u_L) / (S_L - S_M);

    // Normal velocity becomes contact speed
    double u_Lstar = S_M;

    // Tangential velocities modified by magnetic pressure gradient
    double factor = rho_L * (S_L - u_L) / (S_L - S_M);
    double By_Lstar, Bz_Lstar;

    // Magnetic field evolution (complex: involves Alfvén wave transmission)
    // Simplified version: use sign of B_x and pressure jump
    if (std::abs(B_x) > 1.0e-12) {
        double c_A = mhd_->alfven_speed(B_x, rho_Lstar);  // Alfvén speed
        double px_jump = p_m - p_L;
        double delta_p_fact = px_jump / (rho_L * (S_L - u_L) * (S_L - S_M));

        double v_Lstar = v_L - B_x * By_L / (rho_L * (S_L - u_L)) * delta_p_fact;
        double w_Lstar = w_L - B_x * Bz_L / (rho_L * (S_L - u_L)) * delta_p_fact;

        By_Lstar = By_L * (S_L - u_L) / (S_L - S_M);
        Bz_Lstar = Bz_L * (S_L - u_L) / (S_L - S_M);

        // Set velocities
        U_star(1) = rho_Lstar * u_Lstar;
        U_star(2) = rho_Lstar * v_Lstar;
        U_star(3) = rho_Lstar * w_Lstar;
    } else {
        // No normal magnetic field: simple transformation
        U_star(1) = rho_Lstar * u_Lstar;
        U_star(2) = rho_Lstar * v_L;
        U_star(3) = rho_Lstar * w_L;
        By_Lstar = By_L;
        Bz_Lstar = Bz_L;
    }

    // Update state vector
    U_star(0) = rho_Lstar;
    // Components 5 (B_x) stays same, set 6, 7 to new values
    U_star(6) = By_Lstar;
    U_star(7) = Bz_Lstar;

    // Energy update
    double p_L_scaled = p_L + 0.5 * (By_L * By_L + Bz_L * Bz_L);  // Total pressure including magnetic
    double p_m_scaled = p_m + 0.5 * (By_Lstar * By_Lstar + Bz_Lstar * Bz_Lstar);
    U_star(4) = U_L(4) + (p_m_scaled - p_L_scaled) * (S_L - u_L) / (S_L - S_M);

    // Preserve GLM variable if present
    if (nvars > 8) {
        U_star(8) = U_L(8);  // ψ from GLM
    }

    return U_star;
}

Eigen::VectorXd HLLDSolver::compute_state_R(
    const Eigen::VectorXd& U_R, const Eigen::VectorXd& V_R,
    double S_R, double S_M, double p_m, double B_x
) const {
    int nvars = U_R.size();
    Eigen::VectorXd U_star = U_R;

    double rho_R = V_R(0);
    double u_R = V_R(1);
    double v_R = V_R(2);
    double w_R = V_R(3);
    double p_R = V_R(4);
    double By_R = V_R(6);
    double Bz_R = V_R(7);

    // Density in intermediate state
    double rho_Rstar = rho_R * (S_R - u_R) / (S_R - S_M);

    // Normal velocity becomes contact speed
    double u_Rstar = S_M;

    // Tangential velocities modified by magnetic pressure gradient
    double By_Rstar, Bz_Rstar;

    // Magnetic field evolution
    if (std::abs(B_x) > 1.0e-12) {
        double c_A = mhd_->alfven_speed(B_x, rho_Rstar);  // Alfvén speed
        double px_jump = p_m - p_R;
        double delta_p_fact = px_jump / (rho_R * (S_R - u_R) * (S_R - S_M));

        double v_Rstar = v_R - B_x * By_R / (rho_R * (S_R - u_R)) * delta_p_fact;
        double w_Rstar = w_R - B_x * Bz_R / (rho_R * (S_R - u_R)) * delta_p_fact;

        By_Rstar = By_R * (S_R - u_R) / (S_R - S_M);
        Bz_Rstar = Bz_R * (S_R - u_R) / (S_R - S_M);

        // Set velocities
        U_star(1) = rho_Rstar * u_Rstar;
        U_star(2) = rho_Rstar * v_Rstar;
        U_star(3) = rho_Rstar * w_Rstar;
    } else {
        // No normal magnetic field: simple transformation
        U_star(1) = rho_Rstar * u_Rstar;
        U_star(2) = rho_Rstar * v_R;
        U_star(3) = rho_Rstar * w_R;
        By_Rstar = By_R;
        Bz_Rstar = Bz_R;
    }

    // Update state vector
    U_star(0) = rho_Rstar;
    // Component 5 (B_x) stays same, set 6, 7 to new values
    U_star(6) = By_Rstar;
    U_star(7) = Bz_Rstar;

    // Energy update
    double p_R_scaled = p_R + 0.5 * (By_R * By_R + Bz_R * Bz_R);  // Total pressure including magnetic
    double p_m_scaled = p_m + 0.5 * (By_Rstar * By_Rstar + Bz_Rstar * Bz_Rstar);
    U_star(4) = U_R(4) + (p_m_scaled - p_R_scaled) * (S_R - u_R) / (S_R - S_M);

    // Preserve GLM variable if present
    if (nvars > 8) {
        U_star(8) = U_R(8);  // ψ from GLM
    }

    return U_star;
}

} // namespace fvm3d::spatial
