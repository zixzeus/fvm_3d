#pragma once

#include "flux_calculator_base.hpp"
#include "spatial/riemann_solver.hpp"
#include <memory>

namespace fvm3d::spatial {

/**
 * Adapter to use existing RiemannSolver implementations as FluxCalculator.
 *
 * This adapter wraps the existing Riemann solver interface (which uses
 * physics objects internally) into the new FluxCalculator interface.
 *
 * Design pattern: Adapter
 * Purpose: Maintain backwards compatibility while transitioning to new architecture
 */
class RiemannFluxAdapter : public RiemannSolverFlux {
public:
    /**
     * Constructor.
     * @param riemann_solver: Existing Riemann solver implementation
     */
    explicit RiemannFluxAdapter(std::unique_ptr<RiemannSolver> riemann_solver)
        : RiemannSolverFlux(riemann_solver->name()),
          riemann_solver_(std::move(riemann_solver)) {}

    /**
     * Compute flux using wrapped Riemann solver.
     */
    Eigen::VectorXd compute_flux(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        const physics::PhysicsBase& physics,
        int direction
    ) const override {
        // Riemann solver already has physics object, just call solve()
        return riemann_solver_->solve(U_L, U_R, direction);
    }

    /**
     * Compute maximum wave speed using wrapped Riemann solver.
     */
    double compute_max_wave_speed(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        const physics::PhysicsBase& physics,
        int direction
    ) const override {
        return riemann_solver_->max_wave_speed(U_L, U_R, direction);
    }

    /**
     * Access underlying Riemann solver.
     */
    const RiemannSolver& riemann_solver() const { return *riemann_solver_; }

private:
    std::unique_ptr<RiemannSolver> riemann_solver_;
};

} // namespace fvm3d::spatial
