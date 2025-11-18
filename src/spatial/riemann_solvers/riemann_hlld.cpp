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

    // Compute fast magnetosonic speeds (vf^2)
    double B2_L = By_L * By_L + Bz_L * Bz_L + B_x * B_x;
    double B2_R = By_R * By_R + Bz_R * Bz_R + B_x * B_x;

    double f1_L = physics_->gamma() * p_L;
    double f2_L = 4.0 * B_x * B_x;
    double vf2_L = ((f1_L + B2_L) + std::sqrt(std::max((f1_L + B2_L) * (f1_L + B2_L) - f1_L * f2_L, 0.0))) / (2.0 * rho_L);

    double f1_R = physics_->gamma() * p_R;
    double f2_R = 4.0 * B_x * B_x;
    double vf2_R = ((f1_R + B2_R) + std::sqrt(std::max((f1_R + B2_R) * (f1_R + B2_R) - f1_R * f2_R, 0.0))) / (2.0 * rho_R);

    double vf = std::sqrt(std::max(vf2_L, vf2_R));

    // Step 1: Estimate outer wave speeds
    double S_L = std::min(std::min(u_L, u_R) - vf, 0.0);
    double S_R = std::max(std::max(u_L, u_R) + vf, 0.0);

    // Safety check
    if (S_R <= S_L) {
        Eigen::VectorXd F_L = mhd_->flux_x(UL_rot, 0.0, 0.0, 0.0);
        Eigen::VectorXd F_R = mhd_->flux_x(UR_rot, 0.0, 0.0, 0.0);
        return rotate_from_normal(0.5 * (F_L + F_R), direction);
    }

    // Compute original fluxes
    Eigen::VectorXd F_L = mhd_->flux_x(UL_rot, 0.0, 0.0, 0.0);
    Eigen::VectorXd F_R = mhd_->flux_x(UR_rot, 0.0, 0.0, 0.0);

    // Handle simple cases first
    if (0.0 <= S_L) {
        return rotate_from_normal(F_L, direction);
    } else if (0.0 >= S_R) {
        return rotate_from_normal(F_R, direction);
    }

    // Step 2: Compute HLL intermediate state
    double inv_denom = 1.0 / (S_R - S_L);
    Eigen::VectorXd U_hll(nvars);
    Eigen::VectorXd V_hll(nvars);

    U_hll = inv_denom * (S_R * UL_rot - S_L * UR_rot - (F_L - F_R)) / 1.0;
    // Note: For conservative form: U_hll = (S_R * U_R - S_L * U_L - F_R + F_L) / (S_R - S_L)
    // But we need to be careful with the formula

    // Better formulation
    for (int i = 0; i < nvars; i++) {
        if (i == 0) {  // density uses primitive
            U_hll(i) = inv_denom * (S_R * VR_rot(i) - S_L * VL_rot(i) - F_R(i) + F_L(i));
        } else {
            U_hll(i) = inv_denom * (S_R * UL_rot(i) - S_L * UR_rot(i) - (F_L(i) - F_R(i)));
        }
    }
    V_hll = mhd_->conservative_to_primitive(U_hll);

    // Step 3: Contact velocity from HLL state
    double S_M = U_hll(1) / U_hll(0);  // u_hll = (rho*u)_hll / rho_hll

    // Step 4: Middle pressure (using total pressure with magnetic contribution)
    double pt_L = p_L + 0.5 * B2_L;
    double pt_R = p_R + 0.5 * B2_R;
    double p_m = pt_L + rho_L * (S_L - u_L) * (S_M - u_L);

    // Step 5: Compute densities in intermediate states
    double rho_Lstar = rho_L * (S_L - u_L) / (S_L - S_M);
    double rho_Rstar = rho_R * (S_R - u_R) / (S_R - S_M);

    // Step 6: Compute full intermediate states U_*L and U_*R
    Eigen::VectorXd V_Lstar = mhd_->conservative_to_primitive(U_hll);
    V_Lstar(0) = rho_Lstar;
    V_Lstar(1) = S_M;
    Eigen::VectorXd U_Lstar = compute_state_L(UL_rot, VL_rot, U_hll, V_hll, S_L, S_M, p_m, B_x, rho_Lstar);

    Eigen::VectorXd V_Rstar = mhd_->conservative_to_primitive(U_hll);
    V_Rstar(0) = rho_Rstar;
    V_Rstar(1) = S_M;
    Eigen::VectorXd U_Rstar = compute_state_R(UR_rot, VR_rot, U_hll, V_hll, S_R, S_M, p_m, B_x, rho_Rstar);

    // Step 7: Compute Alfvén wave speeds
    double aL1, aR1;
    compute_alfven_speeds(B_x, rho_Lstar, rho_Rstar, aL1, aR1);

    // Step 8: Select flux based on wave structure
    Eigen::VectorXd F_hlld(nvars);
    const double hllg_factor = 1.001;  // Threshold for HLLC-G fallback

    // HLLC-G fallback: revert to HLLC when Alfvén waves degenerate
    if ((S_L >= (S_M - hllg_factor * std::abs(aL1))) ||
        ((S_M + hllg_factor * std::abs(aR1)) >= S_R)) {

        // Use HLLC-G (simplified HLLD without Alfvén wave splitting)
        if (0.0 <= S_M) {
            F_hlld = F_L + S_L * (U_Lstar - UL_rot);
        } else {
            F_hlld = F_R + S_R * (U_Rstar - UR_rot);
        }
    } else {
        // Full HLLD with Alfvén waves
        if (0.0 <= aL1) {
            // Left of left Alfvén wave
            F_hlld = F_L + S_L * (U_Lstar - UL_rot);
        } else if (0.0 >= aR1) {
            // Right of right Alfvén wave
            F_hlld = F_R + S_R * (U_Rstar - UR_rot);
        } else {
            // Between Alfvén waves: compute central state
            Eigen::VectorXd U_central = compute_central_state(U_Lstar, U_Rstar, V_Lstar, V_Rstar,
                                                             S_M, B_x, direction);

            if (0.0 <= S_M) {
                // Central left region
                F_hlld = F_L + S_L * (U_Lstar - UL_rot) + aL1 * (U_central - U_Lstar);
            } else {
                // Central right region
                F_hlld = F_R + S_R * (U_Rstar - UR_rot) + aR1 * (U_central - U_Rstar);
            }
        }
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

void HLLDSolver::compute_alfven_speeds(
    double B_x, double rho_Lstar, double rho_Rstar,
    double& aL1, double& aR1
) const {
    // Alfvén wave speeds: S_M ± c_A
    // c_A = |B_x| / sqrt(rho)
    double c_A_L = std::abs(B_x) / std::sqrt(std::max(rho_Lstar, 1.0e-12));
    double c_A_R = std::abs(B_x) / std::sqrt(std::max(rho_Rstar, 1.0e-12));

    aL1 = -c_A_L;  // S_M - c_A (left Alfvén, negative going)
    aR1 = c_A_R;   // S_M + c_A (right Alfvén, positive going)
}

Eigen::VectorXd HLLDSolver::compute_state_L(
    const Eigen::VectorXd& U_L, const Eigen::VectorXd& V_L,
    const Eigen::VectorXd& U_hll, const Eigen::VectorXd& V_hll,
    double S_L, double S_M, double p_m, double B_x,
    double rho_Lstar
) const {
    int nvars = U_L.size();
    Eigen::VectorXd U_star(nvars);

    double rho_L = V_L(0);
    double u_L = V_L(1);
    double v_L = V_L(2);
    double w_L = V_L(3);
    double p_L = V_L(4);
    double By_L = V_L(6);
    double Bz_L = V_L(7);

    double B_n = U_hll(5);  // B_x from HLL state
    double By_hll = V_hll(6);
    double Bz_hll = V_hll(7);

    // Standard HLLD formula for intermediate state L*
    // HLLD denominator: rho*(S_L - u_L)*(S_L - S_M) - B_x^2
    double denom = rho_L * (S_L - u_L) * (S_L - S_M) - B_n * B_n;

    if (std::abs(denom) < 1.0e-12) {
        // Degenerate case: return U_L
        return U_L;
    }

    double inv_denom = 1.0 / denom;

    // Tangential magnetic field in intermediate state
    double By_Lstar = By_L * inv_denom * (rho_L * (S_L - u_L) * (S_L - u_L) - B_n * B_n);
    double Bz_Lstar = Bz_L * inv_denom * (rho_L * (S_L - u_L) * (S_L - u_L) - B_n * B_n);

    // Tangential velocity components
    double v_Lstar = v_L - inv_denom * B_n * By_L * (S_M - u_L);
    double w_Lstar = w_L - inv_denom * B_n * Bz_L * (S_M - u_L);

    // Assemble intermediate state
    U_star(0) = rho_Lstar;  // density
    U_star(1) = rho_Lstar * S_M;  // momentum (normal)
    U_star(2) = rho_Lstar * v_Lstar;  // momentum (tangential 1)
    U_star(3) = rho_Lstar * w_Lstar;  // momentum (tangential 2)
    U_star(5) = B_n;  // normal magnetic field (unchanged)
    U_star(6) = By_Lstar;  // tangential magnetic field 1
    U_star(7) = Bz_Lstar;  // tangential magnetic field 2

    // Total pressure including magnetic contribution
    double B2_L = By_L * By_L + Bz_L * Bz_L + B_n * B_n;
    double B2_Lstar = By_Lstar * By_Lstar + Bz_Lstar * Bz_Lstar + B_n * B_n;
    double pt_L = p_L + 0.5 * B2_L;
    double pt_Lstar = p_m + 0.5 * B2_Lstar;

    // Energy update with Poynting flux (Miyoshi & Kusano 2005, Eq. 31)
    // Compute v·B terms for Poynting flux
    double vdotB_L = u_L * B_n + v_L * By_L + w_L * Bz_L;
    double vdotB_Lstar = S_M * B_n + v_Lstar * By_Lstar + w_Lstar * Bz_Lstar;

    // Full HLLD energy formula including Poynting flux
    double E_Lstar = ((S_L - u_L) * U_L(4) - pt_L * u_L + pt_Lstar * S_M +
                      B_n * (vdotB_L - vdotB_Lstar)) / (S_L - S_M);

    U_star(4) = E_Lstar;

    // Preserve GLM variable if present
    if (nvars > 8) {
        U_star(8) = U_L(8);  // ψ from GLM
    }

    return U_star;
}

Eigen::VectorXd HLLDSolver::compute_state_R(
    const Eigen::VectorXd& U_R, const Eigen::VectorXd& V_R,
    const Eigen::VectorXd& U_hll, const Eigen::VectorXd& V_hll,
    double S_R, double S_M, double p_m, double B_x,
    double rho_Rstar
) const {
    int nvars = U_R.size();
    Eigen::VectorXd U_star(nvars);

    double rho_R = V_R(0);
    double u_R = V_R(1);
    double v_R = V_R(2);
    double w_R = V_R(3);
    double p_R = V_R(4);
    double By_R = V_R(6);
    double Bz_R = V_R(7);

    double B_n = U_hll(5);  // B_x from HLL state
    double By_hll = V_hll(6);
    double Bz_hll = V_hll(7);

    // Standard HLLD formula for intermediate state R*
    // HLLD denominator: rho*(S_R - u_R)*(S_R - S_M) - B_x^2
    double denom = rho_R * (S_R - u_R) * (S_R - S_M) - B_n * B_n;

    if (std::abs(denom) < 1.0e-12) {
        // Degenerate case: return U_R
        return U_R;
    }

    double inv_denom = 1.0 / denom;

    // Tangential magnetic field in intermediate state
    double By_Rstar = By_R * inv_denom * (rho_R * (S_R - u_R) * (S_R - u_R) - B_n * B_n);
    double Bz_Rstar = Bz_R * inv_denom * (rho_R * (S_R - u_R) * (S_R - u_R) - B_n * B_n);

    // Tangential velocity components
    double v_Rstar = v_R - inv_denom * B_n * By_R * (S_M - u_R);
    double w_Rstar = w_R - inv_denom * B_n * Bz_R * (S_M - u_R);

    // Assemble intermediate state
    U_star(0) = rho_Rstar;  // density
    U_star(1) = rho_Rstar * S_M;  // momentum (normal)
    U_star(2) = rho_Rstar * v_Rstar;  // momentum (tangential 1)
    U_star(3) = rho_Rstar * w_Rstar;  // momentum (tangential 2)
    U_star(5) = B_n;  // normal magnetic field (unchanged)
    U_star(6) = By_Rstar;  // tangential magnetic field 1
    U_star(7) = Bz_Rstar;  // tangential magnetic field 2

    // Total pressure including magnetic contribution
    double B2_R = By_R * By_R + Bz_R * Bz_R + B_n * B_n;
    double B2_Rstar = By_Rstar * By_Rstar + Bz_Rstar * Bz_Rstar + B_n * B_n;
    double pt_R = p_R + 0.5 * B2_R;
    double pt_Rstar = p_m + 0.5 * B2_Rstar;

    // Energy update with Poynting flux (Miyoshi & Kusano 2005, Eq. 31)
    // Compute v·B terms for Poynting flux
    double vdotB_R = u_R * B_n + v_R * By_R + w_R * Bz_R;
    double vdotB_Rstar = S_M * B_n + v_Rstar * By_Rstar + w_Rstar * Bz_Rstar;

    // Full HLLD energy formula including Poynting flux
    double E_Rstar = ((S_R - u_R) * U_R(4) - pt_R * u_R + pt_Rstar * S_M +
                      B_n * (vdotB_R - vdotB_Rstar)) / (S_R - S_M);

    U_star(4) = E_Rstar;

    // Preserve GLM variable if present
    if (nvars > 8) {
        U_star(8) = U_R(8);  // ψ from GLM
    }

    return U_star;
}

Eigen::VectorXd HLLDSolver::compute_central_state(
    const Eigen::VectorXd& U_Lstar, const Eigen::VectorXd& U_Rstar,
    const Eigen::VectorXd& V_Lstar, const Eigen::VectorXd& V_Rstar,
    double S_M, double B_x, int direction
) const {
    int nvars = U_Lstar.size();
    Eigen::VectorXd U_central(nvars);

    double rho_Lstar = V_Lstar(0);
    double rho_Rstar = V_Rstar(0);
    double v_Lstar_t1 = V_Lstar(2);  // tangential 1
    double v_Lstar_t2 = V_Lstar(3);  // tangential 2
    double v_Rstar_t1 = V_Rstar(2);
    double v_Rstar_t2 = V_Rstar(3);
    double By_Lstar = V_Lstar(6);
    double By_Rstar = V_Rstar(6);
    double Bz_Lstar = V_Lstar(7);
    double Bz_Rstar = V_Rstar(7);

    // Roe averaging for central state
    double sqrt_rho_L = std::sqrt(rho_Lstar);
    double sqrt_rho_R = std::sqrt(rho_Rstar);
    double wsum = sqrt_rho_L + sqrt_rho_R;
    double w1 = sqrt_rho_L / wsum;
    double w2 = sqrt_rho_R / wsum;

    // Averaged density
    double rho_central = sqrt_rho_L * sqrt_rho_R;

    // Averaged tangential velocities (Miyoshi & Kusano 2005, Eq. 44)
    double sign_B = (B_x >= 0.0) ? 1.0 : -1.0;
    double inv_wsum = 1.0 / wsum;

    double v_central_t1 = w1 * v_Lstar_t1 + w2 * v_Rstar_t1 +
                         sign_B * inv_wsum * (By_Rstar - By_Lstar);
    double v_central_t2 = w1 * v_Lstar_t2 + w2 * v_Rstar_t2 +
                         sign_B * inv_wsum * (Bz_Rstar - Bz_Lstar);

    // Averaged tangential magnetic fields - CROSS averaging (Miyoshi & Kusano 2005, Eq. 43)
    // Note: Left weight multiplies RIGHT magnetic field (and vice versa)
    double By_central = w1 * By_Rstar + w2 * By_Lstar +
                       sign_B * sqrt_rho_L * sqrt_rho_R * inv_wsum * (v_Rstar_t1 - v_Lstar_t1);
    double Bz_central = w1 * Bz_Rstar + w2 * Bz_Lstar +
                       sign_B * sqrt_rho_L * sqrt_rho_R * inv_wsum * (v_Rstar_t2 - v_Lstar_t2);

    // Assemble central state
    U_central(0) = rho_central;
    U_central(1) = rho_central * S_M;
    U_central(2) = rho_central * v_central_t1;
    U_central(3) = rho_central * v_central_t2;
    U_central(5) = B_x;
    U_central(6) = By_central;
    U_central(7) = Bz_central;

    // Energy: internal + kinetic + magnetic
    // Extract pressure from intermediate states and average
    // For U*L and U*R, we can reconstruct pressure from total energy
    // Simplified approach: use averaged energies from L* and R* states
    double E_Lstar = U_Lstar(4);
    double E_Rstar = U_Rstar(4);

    // Roe-averaged energy as starting point
    double E_avg = w1 * E_Lstar + w2 * E_Rstar;

    // Apply Poynting flux correction for central state (Miyoshi Eq. 45)
    // Correction term accounts for change in v·B across Alfvén waves
    double vdotB_Lstar = S_M * B_x + v_Lstar_t1 * By_Lstar + v_Lstar_t2 * Bz_Lstar;
    double vdotB_central = S_M * B_x + v_central_t1 * By_central + v_central_t2 * Bz_central;
    double vdotB_Rstar = S_M * B_x + v_Rstar_t1 * By_Rstar + v_Rstar_t2 * Bz_Rstar;

    double E_correction = sign_B * sqrt_rho_L * inv_wsum * (vdotB_Lstar - vdotB_central);
    double E_central = E_avg - E_correction;

    U_central(4) = E_central;

    // Preserve GLM variable if present
    if (nvars > 8) {
        U_central(8) = 0.0;
    }

    return U_central;
}

} // namespace fvm3d::spatial
