#include "spatial/riemann_hlld.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

HLLDSolver::HLLDSolver(const std::shared_ptr<physics::AdvancedResistiveMHD3D>& mhd)
    : mhd_(mhd) {
    if (!mhd_) {
        throw std::invalid_argument("MHD physics object cannot be null");
    }
}

Eigen::VectorXd HLLDSolver::solve(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    constexpr double MU0 = 1.0;
    int nvars = U_L.size();

    // Convert to primitive variables
    double rho_L, u_L, v_L, w_L, p_L, Bx_L, By_L, Bz_L;
    double rho_R, u_R, v_R, w_R, p_R, Bx_R, By_R, Bz_R;
    mhd_->conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L, Bx_L, By_L, Bz_L);
    mhd_->conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R, Bx_R, By_R, Bz_R);

    // Rotate to normal direction
    Eigen::VectorXd UL_rot = rotate_to_normal(U_L, direction);
    Eigen::VectorXd UR_rot = rotate_to_normal(U_R, direction);

    // Extract rotated primitive variables
    double rho_L_r, u_L_r, v_L_r, w_L_r, p_L_r, Bx_r, By_L_r, Bz_L_r;
    double rho_R_r, u_R_r, v_R_r, w_R_r, p_R_r, By_R_r, Bz_R_r;
    mhd_->conservative_to_primitive(UL_rot, rho_L_r, u_L_r, v_L_r, w_L_r, p_L_r, Bx_r, By_L_r, Bz_L_r);
    mhd_->conservative_to_primitive(UR_rot, rho_R_r, u_R_r, v_R_r, w_R_r, p_R_r, Bx_r, By_R_r, Bz_R_r);

    // Compute sound speeds
    double a_L = mhd_->sound_speed(rho_L_r, p_L_r);
    double a_R = mhd_->sound_speed(rho_R_r, p_R_r);

    // Compute fast magnetosonic speeds
    double cf_L = mhd_->fast_speed(rho_L_r, p_L_r, Bx_r, By_L_r, Bz_L_r);
    double cf_R = mhd_->fast_speed(rho_R_r, p_R_r, Bx_r, By_R_r, Bz_R_r);

    // Estimate wave speeds (simple HLL-type)
    double S_L = std::min(u_L_r - cf_L, u_R_r - cf_R);
    double S_R = std::max(u_L_r + cf_L, u_R_r + cf_R);

    // Compute fluxes using physics object (position doesn't matter, use 0,0,0)
    Eigen::VectorXd F_L = mhd_->flux_x(UL_rot, 0.0, 0.0, 0.0);
    Eigen::VectorXd F_R = mhd_->flux_x(UR_rot, 0.0, 0.0, 0.0);

    // HLL flux
    Eigen::VectorXd F_hll(nvars);
    if (S_L >= 0.0) {
        F_hll = F_L;
    } else if (S_R <= 0.0) {
        F_hll = F_R;
    } else {
        // HLL intermediate state and flux
        F_hll = (S_R * F_L - S_L * F_R + S_L * S_R * (UR_rot - UL_rot)) / (S_R - S_L);
    }

    // Rotate flux back to original direction
    return rotate_from_normal(F_hll, direction);
}

double HLLDSolver::max_wave_speed(const Eigen::VectorXd& U_L, const Eigen::VectorXd& U_R, int direction) const {
    // Convert to primitive
    double rho_L, u_L, v_L, w_L, p_L, Bx_L, By_L, Bz_L;
    double rho_R, u_R, v_R, w_R, p_R, Bx_R, By_R, Bz_R;
    mhd_->conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L, Bx_L, By_L, Bz_L);
    mhd_->conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R, Bx_R, By_R, Bz_R);

    // Select normal velocity
    double un_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double un_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    // Compute fast magnetosonic speeds
    double cf_L = mhd_->fast_speed(rho_L, p_L, Bx_L, By_L, Bz_L);
    double cf_R = mhd_->fast_speed(rho_R, p_R, Bx_R, By_R, Bz_R);

    // Maximum wave speed
    return std::max(std::abs(un_L) + cf_L, std::abs(un_R) + cf_R);
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

} // namespace fvm3d::spatial
