#include "spatial/flux_calculation/lax_friedrichs.hpp"
#include <algorithm>
#include <cmath>

namespace fvm3d::spatial {

LaxFriedrichsFlux::LaxFriedrichsFlux(std::shared_ptr<physics::PhysicsBase> physics)
    : FluxCalculator("Lax-Friedrichs", "riemann"), physics_(physics) {
    if (!physics_) {
        throw std::invalid_argument("LaxFriedrichsFlux: physics object cannot be null");
    }
}

Eigen::VectorXd LaxFriedrichsFlux::compute_flux(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    const physics::PhysicsBase& physics,
    int direction
) const {
    // Compute fluxes using provided physics object
    Eigen::VectorXd F_L = physics.compute_flux(U_L, direction);
    Eigen::VectorXd F_R = physics.compute_flux(U_R, direction);

    // Compute maximum wave speed
    double lambda = compute_max_wave_speed(U_L, U_R, physics, direction);

    // Lax-Friedrichs flux: 0.5 * (F_L + F_R) - 0.5 * lambda * (U_R - U_L)
    Eigen::VectorXd flux = 0.5 * (F_L + F_R) - 0.5 * lambda * (U_R - U_L);

    return flux;
}

double LaxFriedrichsFlux::compute_max_wave_speed(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    const physics::PhysicsBase& physics,
    int direction
) const {
    // Delegate to physics object - it knows the correct wave speeds for the equation type
    // (e.g., acoustic waves for Euler, magnetosonic/Alfvén waves for MHD)
    double lambda_L = physics.max_wave_speed(U_L, direction);
    double lambda_R = physics.max_wave_speed(U_R, direction);

    return std::max(lambda_L, lambda_R);
}

} // namespace fvm3d::spatial
