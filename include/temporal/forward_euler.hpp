#pragma once

#include "time_integrator.hpp"

namespace fvm3d::temporal {

/**
 * Forward Euler time integrator.
 * U^(n+1) = U^n + dt * RHS(U^n)
 * Order: 1st order, CFL <= 1
 */
class ForwardEuler : public TimeIntegrator {
public:
    void step(
        StateField3D& U_current,
        double dt,
        const RHSFunction& rhs
    ) override;

    int order() const override { return 1; }
    std::string name() const override { return "Forward Euler"; }
};

} // namespace fvm3d::temporal
