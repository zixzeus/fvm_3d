#include "spatial/riemann_hll.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

Eigen::VectorXd HLLSolver::solve(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double S_L = estimate_s_left(U_L, U_R, direction);
    double S_R = estimate_s_right(U_L, U_R, direction);

    Eigen::VectorXd F_L = compute_flux(U_L, direction);
    Eigen::VectorXd F_R = compute_flux(U_R, direction);

    if (S_R <= S_L) {
        return 0.5 * (F_L + F_R);
    }

    if (0.0 <= S_L) {
        return F_L;
    } else if (0.0 >= S_R) {
        return F_R;
    } else {
        // HLL flux in the middle region
        Eigen::VectorXd flux = (S_R * F_L - S_L * F_R + S_L * S_R * (U_R - U_L)) / (S_R - S_L);
        return flux;
    }
}

double HLLSolver::max_wave_speed(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double S_L = estimate_s_left(U_L, U_R, direction);
    double S_R = estimate_s_right(U_L, U_R, direction);
    return std::max(std::abs(S_L), std::abs(S_R));
}

double HLLSolver::estimate_s_left(
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

double HLLSolver::estimate_s_right(
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

Eigen::VectorXd HLLSolver::compute_flux(
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

} // namespace fvm3d::spatial
