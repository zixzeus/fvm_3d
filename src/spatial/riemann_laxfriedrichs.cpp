#include "spatial/riemann_laxfriedrichs.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

Eigen::VectorXd LaxFriedrichsSolver::solve(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction
) const {
    // Compute fluxes
    Eigen::VectorXd F_L = compute_flux(U_L, direction);
    Eigen::VectorXd F_R = compute_flux(U_R, direction);

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
    double rho_L, u_L, v_L, w_L, p_L;
    double rho_R, u_R, v_R, w_R, p_R;

    conservative_to_primitive(U_L, rho_L, u_L, v_L, w_L, p_L);
    conservative_to_primitive(U_R, rho_R, u_R, v_R, w_R, p_R);

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

Eigen::VectorXd LaxFriedrichsSolver::compute_flux(
    const Eigen::VectorXd& U,
    int direction
) const {
    double rho, u, v, w, p;
    conservative_to_primitive(U, rho, u, v, w, p);

    Eigen::VectorXd flux(5);

    if (direction == 0) {
        // X-direction flux
        flux(0) = rho * u;
        flux(1) = rho * u * u + p;
        flux(2) = rho * u * v;
        flux(3) = rho * u * w;
        flux(4) = (U(4) + p) * u;
    } else if (direction == 1) {
        // Y-direction flux
        flux(0) = rho * v;
        flux(1) = rho * v * u;
        flux(2) = rho * v * v + p;
        flux(3) = rho * v * w;
        flux(4) = (U(4) + p) * v;
    } else {
        // Z-direction flux
        flux(0) = rho * w;
        flux(1) = rho * w * u;
        flux(2) = rho * w * v;
        flux(3) = rho * w * w + p;
        flux(4) = (U(4) + p) * w;
    }

    return flux;
}

} // namespace fvm3d::spatial
