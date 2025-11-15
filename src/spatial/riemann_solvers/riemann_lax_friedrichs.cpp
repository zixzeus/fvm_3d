#include "spatial/riemann_solvers/riemann_lax_friedrichs.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

LaxFriedrichsSolver::LaxFriedrichsSolver(const std::shared_ptr<physics::PhysicsBase>& physics)
    : physics_(physics) {
    if (!physics_) {
        throw std::invalid_argument("LaxFriedrichsSolver: physics object cannot be null");
    }
}

Eigen::VectorXd LaxFriedrichsSolver::solve(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    // Compute fluxes using physics object
    Eigen::VectorXd F_L = physics_->compute_flux(U_L, direction);
    Eigen::VectorXd F_R = physics_->compute_flux(U_R, direction);

    // Compute maximum wave speed
    double lambda = max_wave_speed(U_L, U_R, direction);

    // Lax-Friedrichs flux: 0.5 * (F_L + F_R) - 0.5 * lambda * (U_R - U_L)
    Eigen::VectorXd flux = 0.5 * (F_L + F_R) - 0.5 * lambda * (U_R - U_L);

    return flux;
}

double LaxFriedrichsSolver::max_wave_speed(
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

    // Compute sound speeds
    double a_L = sound_speed(rho_L, p_L);
    double a_R = sound_speed(rho_R, p_R);

    // Select normal velocity based on direction
    double u_normal_L = (direction == 0) ? u_L : (direction == 1) ? v_L : w_L;
    double u_normal_R = (direction == 0) ? u_R : (direction == 1) ? v_R : w_R;

    // Maximum wave speed: max(|u| + a)
    double lambda_L = std::abs(u_normal_L) + a_L;
    double lambda_R = std::abs(u_normal_R) + a_R;

    return std::max(lambda_L, lambda_R);
}

} // namespace fvm3d::spatial
