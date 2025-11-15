#pragma once

#include "riemann_solver.hpp"
#include "spatial/flux_calculation/flux_calculator_base.hpp"
#include <memory>

namespace fvm3d::spatial {

/**
 * Adapter to use FluxCalculator as RiemannSolver.
 *
 * This is the reverse of RiemannFluxAdapter - it wraps a FluxCalculator
 * (which has physics object in its compute_flux signature) to work with
 * the legacy RiemannSolver interface (which has physics internally).
 *
 * Design pattern: Adapter
 * Purpose: Allow central schemes (like LaxFriedrichs) to be used via
 *          RiemannSolver interface for backward compatibility
 */
class FluxCalculatorAdapter : public RiemannSolver {
public:
    /**
     * Constructor.
     * @param flux_calculator: FluxCalculator implementation
     * @param physics: Physics object to use for flux computation
     */
    explicit FluxCalculatorAdapter(
        std::unique_ptr<FluxCalculator> flux_calculator,
        std::shared_ptr<physics::PhysicsBase> physics
    ) : flux_calculator_(std::move(flux_calculator)),
        physics_(physics) {}

    /**
     * Solve Riemann problem using wrapped FluxCalculator.
     */
    Eigen::VectorXd solve(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const override {
        return flux_calculator_->compute_flux(U_L, U_R, *physics_, direction);
    }

    /**
     * Get maximum wave speed using wrapped FluxCalculator.
     */
    double max_wave_speed(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction
    ) const override {
        return flux_calculator_->compute_max_wave_speed(U_L, U_R, *physics_, direction);
    }

    /**
     * Get name from wrapped FluxCalculator.
     */
    std::string name() const override {
        return flux_calculator_->name();
    }

    /**
     * Access underlying FluxCalculator.
     */
    const FluxCalculator& flux_calculator() const { return *flux_calculator_; }

private:
    std::unique_ptr<FluxCalculator> flux_calculator_;
    std::shared_ptr<physics::PhysicsBase> physics_;
};

} // namespace fvm3d::spatial
