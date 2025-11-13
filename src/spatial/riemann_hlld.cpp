#include "spatial/riemann_hlld.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

HLLDSolver::HLLDSolver(const std::shared_ptr<physics::ResistiveMHD3D>& mhd)
    : mhd_(mhd) {
    if (!mhd_) {
        throw std::invalid_argument("ResistiveMHD3D instance required for HLLD solver");
    }
}

Eigen::VectorXd HLLDSolver::solve(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    // Extract primitive variables
    double rho_L, u_L, v_L, w_L, p_L, Bx_L, By_L, Bz_L;
    mhd_->conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L, Bx_L, By_L, Bz_L);

    double rho_R, u_R, v_R, w_R, p_R, Bx_R, By_R, Bz_R;
    mhd_->conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R, Bx_R, By_R, Bz_R);

    // Rotate states to align with normal direction
    Eigen::VectorXd U_L_rot = rotate_to_normal(U_L, direction);
    Eigen::VectorXd U_R_rot = rotate_to_normal(U_R, direction);

    // Extract rotated primitive variables
    double rho_L_r, u_L_r, v_L_r, w_L_r, p_L_r, Bx_r, By_L_r, Bz_L_r;
    mhd_->conservative_to_primitive(U_L_rot, rho_L_r, u_L_r, v_L_r, w_L_r, p_L_r, Bx_r, By_L_r, Bz_L_r);

    double rho_R_r, u_R_r, v_R_r, w_R_r, p_R_r, By_R_r, Bz_R_r;
    mhd_->conservative_to_primitive(U_R_rot, rho_R_r, u_R_r, v_R_r, w_R_r, p_R_r, Bx_r, By_R_r, Bz_R_r);

    // Normal magnetic field component (same in all directions)
    double Bx = Bx_r;

    // Estimate wave speeds
    WaveSpeeds speeds = estimate_wave_speeds(
        rho_L_r, u_L_r, v_L_r, w_L_r, p_L_r, Bx, By_L_r, Bz_L_r,
        rho_R_r, u_R_r, v_R_r, w_R_r, p_R_r, By_R_r, Bz_R_r
    );

    // Compute HLLD intermediate states
    HLLDStates states = compute_hlld_states(U_L_rot, U_R_rot, speeds, direction);

    // Select appropriate state and compute flux
    Eigen::VectorXd U_star(8);
    if (0.0 <= speeds.S_fL) {
        U_star = states.U_L;
    } else if (speeds.S_fL < 0.0 && 0.0 <= speeds.S_aL) {
        U_star = states.U_fL_star;
    } else if (speeds.S_aL < 0.0 && 0.0 <= speeds.S_c) {
        U_star = states.U_aL_star;
    } else if (speeds.S_c < 0.0 && 0.0 <= speeds.S_aR) {
        U_star = states.U_aR_star;
    } else if (speeds.S_aR < 0.0 && 0.0 <= speeds.S_fR) {
        U_star = states.U_fR_star;
    } else {
        U_star = states.U_R;
    }

    // Compute flux from selected state (rotated back)
    Eigen::VectorXd U_star_orig = rotate_from_normal(U_star, direction);

    // Return flux in original (non-rotated) direction
    if (direction == 0) {
        return mhd_->flux_x(U_star_orig);
    } else if (direction == 1) {
        return mhd_->flux_y(U_star_orig);
    } else {
        return mhd_->flux_z(U_star_orig);
    }
}

double HLLDSolver::max_wave_speed(const Eigen::VectorXd& U_L, const Eigen::VectorXd& U_R, int direction) const {
    // Use fast magnetosonic speed as upper bound in specified direction
    double speed_L = mhd_->max_wave_speed(U_L, direction);
    double speed_R = mhd_->max_wave_speed(U_R, direction);
    return std::max(speed_L, speed_R);
}

HLLDSolver::WaveSpeeds HLLDSolver::estimate_wave_speeds(
    double rho_L, double u_L, double v_L, double w_L, double p_L,
    double Bx, double By_L, double Bz_L,
    double rho_R, double u_R, double v_R, double w_R, double p_R,
    double By_R, double Bz_R
) const {
    // Roe average
    double sqrt_rho_L = std::sqrt(rho_L);
    double sqrt_rho_R = std::sqrt(rho_R);
    double w_sum = sqrt_rho_L + sqrt_rho_R;

    double u_roe = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) / w_sum;
    double v_roe = (sqrt_rho_L * v_L + sqrt_rho_R * v_R) / w_sum;
    double w_roe = (sqrt_rho_L * w_L + sqrt_rho_R * w_R) / w_sum;

    // Compute sound speeds
    double a_L = mhd_->sound_speed(rho_L, p_L);
    double a_R = mhd_->sound_speed(rho_R, p_R);
    double a_roe = 0.5 * (a_L + a_R);

    // Alfvén speeds
    double va_L = mhd_->alfven_speed(Bx, rho_L);
    double va_R = mhd_->alfven_speed(Bx, rho_R);
    double va_roe = 0.5 * (va_L + va_R);

    // Fast magnetosonic speed
    double B_sq_L = By_L*By_L + Bz_L*Bz_L + Bx*Bx;
    double B_sq_R = By_R*By_R + Bz_R*Bz_R + Bx*Bx;
    double B_sq_norm_L = B_sq_L / (1.0 * rho_L);  // MU0 = 1.0
    double B_sq_norm_R = B_sq_R / (1.0 * rho_R);

    double cf_L = std::sqrt(a_L*a_L + B_sq_norm_L);
    double cf_R = std::sqrt(a_R*a_R + B_sq_norm_R);
    double cf_roe = std::max(cf_L, cf_R);

    // Estimate contact speed (approximated by Roe average velocity)
    double S_c = u_roe;

    // Wave speed estimates
    WaveSpeeds speeds;
    speeds.S_fL = std::min(u_L, u_roe) - cf_roe;
    speeds.S_aL = u_roe - va_roe;
    speeds.S_c = u_roe;
    speeds.S_aR = u_roe + va_roe;
    speeds.S_fR = std::max(u_R, u_roe) + cf_roe;

    return speeds;
}

HLLDSolver::HLLDStates HLLDSolver::compute_hlld_states(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    const WaveSpeeds& speeds,
    int direction
) const {
    HLLDStates states;
    states.U_L = U_L;
    states.U_R = U_R;

    // For simplicity, use HLL-type averaging for intermediate states
    // A full HLLD implementation would compute exact intermediate states
    // This is a simplified version that still captures main wave structure

    // Left intermediate state (approximate)
    if (speeds.S_fL < 0.0 && speeds.S_c >= 0.0) {
        double weight = (speeds.S_c - speeds.S_fL) / (0.0 - speeds.S_fL + 1e-16);
        states.U_fL_star = weight * U_L + (1.0 - weight) * U_R;
        states.U_aL_star = states.U_fL_star;
    } else {
        states.U_fL_star = U_L;
        states.U_aL_star = U_L;
    }

    // Right intermediate state (approximate)
    if (speeds.S_fR > 0.0 && speeds.S_c <= 0.0) {
        double weight = (0.0 - speeds.S_fR) / (speeds.S_c - speeds.S_fR + 1e-16);
        states.U_fR_star = weight * U_L + (1.0 - weight) * U_R;
        states.U_aR_star = states.U_fR_star;
    } else {
        states.U_fR_star = U_R;
        states.U_aR_star = U_R;
    }

    return states;
}

Eigen::VectorXd HLLDSolver::rotate_to_normal(const Eigen::VectorXd& U, int direction) const {
    Eigen::VectorXd U_rot = U;

    if (direction == 0) {
        // X is normal: no rotation
        return U_rot;
    } else if (direction == 1) {
        // Y is normal: rotate x→y, y→x
        std::swap(U_rot(1), U_rot(2));  // u ↔ v
        std::swap(U_rot(5), U_rot(6));  // Bx ↔ By
    } else if (direction == 2) {
        // Z is normal: rotate x→z, y→z, z→x
        double temp = U_rot(1);
        U_rot(1) = U_rot(3);
        U_rot(3) = temp;
        temp = U_rot(5);
        U_rot(5) = U_rot(7);
        U_rot(7) = temp;
    }

    return U_rot;
}

Eigen::VectorXd HLLDSolver::rotate_from_normal(const Eigen::VectorXd& U, int direction) const {
    Eigen::VectorXd U_orig = U;

    if (direction == 0) {
        // X is normal: no rotation
        return U_orig;
    } else if (direction == 1) {
        // Y was normal: rotate back x→y, y→x
        std::swap(U_orig(1), U_orig(2));
        std::swap(U_orig(5), U_orig(6));
    } else if (direction == 2) {
        // Z was normal: rotate back
        double temp = U_orig(1);
        U_orig(1) = U_orig(3);
        U_orig(3) = temp;
        temp = U_orig(5);
        U_orig(5) = U_orig(7);
        U_orig(7) = temp;
    }

    return U_orig;
}

} // namespace fvm3d::spatial
