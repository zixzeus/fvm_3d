#include "spatial/flux_calculation/riemann_flux.hpp"
#include "spatial/riemann_solvers/riemann_solver_factory.hpp"
#include <algorithm>

namespace fvm3d::spatial {

RiemannFlux::RiemannFlux(const std::string& solver_type, std::shared_ptr<physics::PhysicsBase> physics)
    : FluxCalculator("Riemann_" + solver_type, "riemann"),
      physics_(physics),
      solver_type_(solver_type) {

    if (!physics_) {
        throw std::invalid_argument("RiemannFlux: physics object cannot be null");
    }

    // Convert solver type to lowercase for case-insensitive matching
    std::string solver_lower = solver_type;
    std::transform(solver_lower.begin(), solver_lower.end(), solver_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    // Determine physics type for Riemann solver factory
    std::string physics_type = "euler";  // default
    int num_vars = 5;  // default for Euler

    // Try to determine physics type from physics object
    // This is a simple heuristic - could be improved
    if (physics_->num_variables() == 9) {
        physics_type = "mhd_advanced";
        num_vars = 9;
    }

    // Create the appropriate Riemann solver
    try {
        riemann_solver_ = RiemannSolverFactory::create(solver_lower, physics_type, num_vars);
    } catch (const std::exception& e) {
        throw std::invalid_argument("RiemannFlux: Failed to create Riemann solver '" +
                                    solver_type + "': " + e.what());
    }
}

Eigen::VectorXd RiemannFlux::compute_flux(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    const physics::PhysicsBase& physics,
    int direction
) const {
    // Use internal Riemann solver to compute flux
    // Note: We ignore the 'physics' parameter and use our internal physics_
    return riemann_solver_->solve(U_L, U_R, direction);
}

double RiemannFlux::compute_max_wave_speed(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    const physics::PhysicsBase& physics,
    int direction
) const {
    // Use internal Riemann solver to estimate wave speeds
    auto [s_left, s_right] = riemann_solver_->estimate_wave_speeds(U_L, U_R, direction);
    return std::max(std::abs(s_left), std::abs(s_right));
}

bool RiemannFlux::is_physics_agnostic() const {
    // HLLD is MHD-specific, others are more general
    std::string solver_lower = solver_type_;
    std::transform(solver_lower.begin(), solver_lower.end(), solver_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    return solver_lower != "hlld";
}

bool RiemannFlux::supports_physics(const std::string& physics_type) const {
    std::string physics_lower = physics_type;
    std::transform(physics_lower.begin(), physics_lower.end(), physics_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    std::string solver_lower = solver_type_;
    std::transform(solver_lower.begin(), solver_lower.end(), solver_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    // HLLD is specialized for MHD
    if (solver_lower == "hlld") {
        return physics_lower == "mhd" || physics_lower == "mhd_advanced" ||
               physics_lower == "magnetohydrodynamics";
    }

    // Other Riemann solvers support most hyperbolic systems
    return physics_lower == "euler" || physics_lower == "mhd" ||
           physics_lower == "mhd_advanced" || physics_lower == "compressible_flow" ||
           physics_lower == "gas_dynamics";
}

} // namespace fvm3d::spatial
