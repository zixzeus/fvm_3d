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
    auto [S_L, S_R] = estimate_wave_speeds(U_L, U_R, direction);

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

} // namespace fvm3d::spatial
