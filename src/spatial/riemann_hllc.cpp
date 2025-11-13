#include "spatial/riemann_hllc.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

Eigen::VectorXd HLLCSolver::solve(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double S_L = estimate_s_left(U_L, U_R, direction);
    double S_R = estimate_s_right(U_L, U_R, direction);

    // Ensure proper ordering
    if (S_R <= S_L) {
        return 0.5 * (compute_flux(U_L, direction) + compute_flux(U_R, direction));
    }

    double p_m = estimate_p_middle(U_L, U_R, direction);
    double S_M = estimate_s_contact(U_L, U_R, direction, S_L, S_R);

    Eigen::VectorXd F_L = compute_flux(U_L, direction);
    Eigen::VectorXd F_R = compute_flux(U_R, direction);

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

double HLLCSolver::max_wave_speed(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double S_L = estimate_s_left(U_L, U_R, direction);
    double S_R = estimate_s_right(U_L, U_R, direction);
    return std::max(std::abs(S_L), std::abs(S_R));
}

double HLLCSolver::estimate_s_left(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double rho_L, u_L, v_L, w_L, p_L;
    double rho_R, u_R, v_R, w_R, p_R;

    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

    double a_L = sound_speed(rho_L, p_L);
    double a_R = sound_speed(rho_R, p_R);

    double u_normal_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double u_normal_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    return std::min(u_normal_L - a_L, u_normal_R - a_R);
}

double HLLCSolver::estimate_s_right(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double rho_L, u_L, v_L, w_L, p_L;
    double rho_R, u_R, v_R, w_R, p_R;

    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

    double a_L = sound_speed(rho_L, p_L);
    double a_R = sound_speed(rho_R, p_R);

    double u_normal_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double u_normal_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    return std::max(u_normal_L + a_L, u_normal_R + a_R);
}

double HLLCSolver::estimate_p_middle(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double rho_L, u_L, v_L, w_L, p_L;
    double rho_R, u_R, v_R, w_R, p_R;

    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

    double a_L = sound_speed(rho_L, p_L);
    double a_R = sound_speed(rho_R, p_R);

    double u_normal_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double u_normal_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    // Roe's formula for middle pressure
    double numerator = a_R * p_L + a_L * p_R + a_L * a_R * (u_normal_L - u_normal_R);
    double denominator = a_L + a_R;

    return std::max(numerator / denominator, P_FLOOR);
}

double HLLCSolver::estimate_s_contact(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction,
    double S_L,
    double S_R
) const {
    double rho_L, u_L, v_L, w_L, p_L;
    double rho_R, u_R, v_R, w_R, p_R;

    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

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

Eigen::VectorXd HLLCSolver::compute_flux(
    const Eigen::VectorXd& U,
    int direction
) const {
    double rho, u, v, w, p;
    conservative_to_primitive(U, rho, u, v, w, p);

    Eigen::VectorXd flux(5);

    if (direction == 0) {
        flux(0) = rho * u;
        flux(1) = rho * u * u + p;
        flux(2) = rho * u * v;
        flux(3) = rho * u * w;
        flux(4) = (U(4) + p) * u;
    } else if (direction == 1) {
        flux(0) = rho * v;
        flux(1) = rho * v * u;
        flux(2) = rho * v * v + p;
        flux(3) = rho * v * w;
        flux(4) = (U(4) + p) * v;
    } else {
        flux(0) = rho * w;
        flux(1) = rho * w * u;
        flux(2) = rho * w * v;
        flux(3) = rho * w * w + p;
        flux(4) = (U(4) + p) * w;
    }

    return flux;
}

Eigen::VectorXd HLLCSolver::compute_u_left_hllc(
    const Eigen::VectorXd& U_L,
    double S_L,
    double S_M,
    double p_m,
    int direction
) const {
    double rho_L, u_L, v_L, w_L, p_L;
    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);

    Eigen::VectorXd U_Lm(5);
    double coeff = rho_L * (S_L - velocity_in_direction(U_L, direction)) /
                   (S_L - S_M);

    U_Lm(0) = coeff;
    U_Lm(direction + 1) = coeff * S_M;

    // Tangential momentum components unchanged
    for (int d = 1; d <= 3; d++) {
        if (d != direction + 1) {
            U_Lm(d) = coeff * velocity_in_direction(U_L, d - 1);
        }
    }

    // Energy with updated pressure
    U_Lm(4) = U_L(4) + (p_m - p_L) * (S_L - velocity_in_direction(U_L, direction)) /
                        (S_L - S_M);

    return U_Lm;
}

Eigen::VectorXd HLLCSolver::compute_u_right_hllc(
    const Eigen::VectorXd& U_R,
    double S_R,
    double S_M,
    double p_m,
    int direction
) const {
    double rho_R, u_R, v_R, w_R, p_R;
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

    Eigen::VectorXd U_Rm(5);
    double coeff = rho_R * (S_R - velocity_in_direction(U_R, direction)) /
                   (S_R - S_M);

    U_Rm(0) = coeff;
    U_Rm(direction + 1) = coeff * S_M;

    // Tangential momentum components unchanged
    for (int d = 1; d <= 3; d++) {
        if (d != direction + 1) {
            U_Rm(d) = coeff * velocity_in_direction(U_R, d - 1);
        }
    }

    // Energy with updated pressure
    U_Rm(4) = U_R(4) + (p_m - p_R) * (S_R - velocity_in_direction(U_R, direction)) /
                        (S_R - S_M);

    return U_Rm;
}

} // namespace fvm3d::spatial
