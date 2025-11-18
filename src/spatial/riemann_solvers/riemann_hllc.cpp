#include "spatial/riemann_solvers/riemann_hllc.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

HLLCSolver::HLLCSolver(const std::shared_ptr<physics::PhysicsBase>& physics)
    : RiemannSolver(physics) {
    if (!physics_) {
        throw std::invalid_argument("HLLCSolver: physics object cannot be null");
    }
}

Eigen::VectorXd HLLCSolver::solve(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    auto [S_L, S_R] = estimate_wave_speeds(U_L, U_R, direction);

    // Ensure proper ordering
    if (S_R <= S_L) {
        return 0.5 * (physics_->compute_flux(U_L, direction) + physics_->compute_flux(U_R, direction));
    }

    double p_m = estimate_p_middle(U_L, U_R, direction);
    double S_M = estimate_s_contact(U_L, U_R, direction, S_L, S_R);

    Eigen::VectorXd F_L = physics_->compute_flux(U_L, direction);
    Eigen::VectorXd F_R = physics_->compute_flux(U_R, direction);

    if (0.0 <= S_L) {
        return F_L;
    } else if (0.0 >= S_R) {
        return F_R;
    } else if (0.0 <= S_M) {
        // In left or middle region
        Eigen::VectorXd U_Lm = compute_u_left_hllc(U_L, S_L, S_M, p_m, direction);
        Eigen::VectorXd F_Lm = F_L + S_L * (U_Lm - U_L);
        return F_Lm;
    } else {
        // In right region
        Eigen::VectorXd U_Rm = compute_u_right_hllc(U_R, S_R, S_M, p_m, direction);
        Eigen::VectorXd F_Rm = F_R + S_R * (U_Rm - U_R);
        return F_Rm;
    }
}

double HLLCSolver::estimate_p_middle(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    // Convert to primitive variables using physics object
    Eigen::VectorXd V_L = physics_->conservative_to_primitive(U_L);
    Eigen::VectorXd V_R = physics_->conservative_to_primitive(U_R);

    double rho_L = V_L(0);
    double u_L = V_L(1);
    double v_L = V_L(2);
    double w_L = V_L(3);
    double p_L = V_L(4);

    double rho_R = V_R(0);
    double u_R = V_R(1);
    double v_R = V_R(2);
    double w_R = V_R(3);
    double p_R = V_R(4);

    double a_L = sound_speed(rho_L, p_L);
    double a_R = sound_speed(rho_R, p_R);

    double u_normal_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double u_normal_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    // Roe's formula for middle pressure
    double numerator = a_R * p_L + a_L * p_R + a_L * a_R * (u_normal_L - u_normal_R);
    double denominator = a_L + a_R;

    return std::max(numerator / denominator, physics_->p_floor());
}

double HLLCSolver::estimate_s_contact(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction,
    double S_L,
    double S_R
) const {
    // Convert to primitive variables using physics object
    Eigen::VectorXd V_L = physics_->conservative_to_primitive(U_L);
    Eigen::VectorXd V_R = physics_->conservative_to_primitive(U_R);

    double rho_L = V_L(0);
    double u_L = V_L(1);
    double v_L = V_L(2);
    double w_L = V_L(3);
    double p_L = V_L(4);

    double rho_R = V_R(0);
    double u_R = V_R(1);
    double v_R = V_R(2);
    double w_R = V_R(3);
    double p_R = V_R(4);

    double u_normal_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double u_normal_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    // Contact speed using Roe averaging
    double sqrt_rho_L = std::sqrt(rho_L);
    double sqrt_rho_R = std::sqrt(rho_R);
    double S_M = (sqrt_rho_L * u_normal_L + sqrt_rho_R * u_normal_R) /
                 (sqrt_rho_L + sqrt_rho_R);

    // Alternative: simple average (less accurate but more stable)
    // double S_M = 0.5 * (u_normal_L + u_normal_R);

    return S_M;
}

Eigen::VectorXd HLLCSolver::compute_u_left_hllc(
    const Eigen::VectorXd& U_L,
    double S_L,
    double S_M,
    double p_m,
    int direction
) const {
    // Convert to primitive variables using physics object
    Eigen::VectorXd V_L = physics_->conservative_to_primitive(U_L);

    double rho_L = V_L(0);
    double u_L = V_L(1);
    double v_L = V_L(2);
    double w_L = V_L(3);
    double p_L = V_L(4);

    // Get velocity in direction
    double velocities[3] = {u_L, v_L, w_L};
    double u_normal_L = velocities[direction];

    Eigen::VectorXd U_Lm(5);
    double coeff = rho_L * (S_L - u_normal_L) / (S_L - S_M);

    U_Lm(0) = coeff;
    U_Lm(direction + 1) = coeff * S_M;

    // Tangential momentum components unchanged (SIMD vectorization)
    #pragma omp simd
    for (int d = 1; d <= 3; d++) {
        if (d != direction + 1) {
            U_Lm(d) = coeff * velocities[d - 1];
        }
    }

    // Energy with updated pressure
    U_Lm(4) = U_L(4) + (p_m - p_L) * (S_L - u_normal_L) / (S_L - S_M);

    return U_Lm;
}

Eigen::VectorXd HLLCSolver::compute_u_right_hllc(
    const Eigen::VectorXd& U_R,
    double S_R,
    double S_M,
    double p_m,
    int direction
) const {
    // Convert to primitive variables using physics object
    Eigen::VectorXd V_R = physics_->conservative_to_primitive(U_R);

    double rho_R = V_R(0);
    double u_R = V_R(1);
    double v_R = V_R(2);
    double w_R = V_R(3);
    double p_R = V_R(4);

    // Get velocity in direction
    double velocities[3] = {u_R, v_R, w_R};
    double u_normal_R = velocities[direction];

    Eigen::VectorXd U_Rm(5);
    double coeff = rho_R * (S_R - u_normal_R) / (S_R - S_M);

    U_Rm(0) = coeff;
    U_Rm(direction + 1) = coeff * S_M;

    // Tangential momentum components unchanged (SIMD vectorization)
    #pragma omp simd
    for (int d = 1; d <= 3; d++) {
        if (d != direction + 1) {
            U_Rm(d) = coeff * velocities[d - 1];
        }
    }

    // Energy with updated pressure
    U_Rm(4) = U_R(4) + (p_m - p_R) * (S_R - u_normal_R) / (S_R - S_M);

    return U_Rm;
}

} // namespace fvm3d::spatial
