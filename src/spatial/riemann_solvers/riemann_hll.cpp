#include "spatial/riemann_solvers/riemann_hll.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

HLLSolver::HLLSolver(const std::shared_ptr<physics::PhysicsBase>& physics)
    : RiemannSolver(physics) {
    if (!physics_) {
        throw std::invalid_argument("HLLSolver: physics object cannot be null");
    }
}

Eigen::VectorXd HLLSolver::solve(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    double S_L = estimate_s_left(U_L, U_R, direction);
    double S_R = estimate_s_right(U_L, U_R, direction);

    // Compute fluxes using physics object
    Eigen::VectorXd F_L = physics_->compute_flux(U_L, direction);
    Eigen::VectorXd F_R = physics_->compute_flux(U_R, direction);

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
    // Convert to primitive variables using physics object
    Eigen::VectorXd V_L = physics_->conservative_to_primitive(U_L);
    Eigen::VectorXd V_R = physics_->conservative_to_primitive(U_R);

    // Extract primitive variables (common to both Euler and MHD)
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

    return std::min(u_normal_L - a_L, u_normal_R - a_R);
}

double HLLSolver::estimate_s_right(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    // Convert to primitive variables using physics object
    Eigen::VectorXd V_L = physics_->conservative_to_primitive(U_L);
    Eigen::VectorXd V_R = physics_->conservative_to_primitive(U_R);

    // Extract primitive variables (common to both Euler and MHD)
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

    return std::max(u_normal_L + a_L, u_normal_R + a_R);
}

} // namespace fvm3d::spatial
