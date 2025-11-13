#pragma once

#include "time_integrator.hpp"

namespace fvm3d::temporal {

/**
 * SSP-RK2 (Strong Stability Preserving RK2)
 * U^(1) = U^n + dt * RHS(U^n)
 * U^(n+1) = 0.5 * U^n + 0.5 * U^(1) + 0.5 * dt * RHS(U^(1))
 * Order: 2nd order, CFL <= 1, TVD property
 */
class RK2 : public TimeIntegrator {
public:
    void step(
        StateField3D& U_current,
        double dt,
        const RHSFunction& rhs
    ) override;

    int order() const override { return 2; }
    std::string name() const override { return "SSP-RK2"; }
};

/**
 * SSP-RK3 (Strong Stability Preserving RK3)
 * U^(1) = U^n + dt * RHS(U^n)
 * U^(2) = 0.75 * U^n + 0.25 * U^(1) + 0.25 * dt * RHS(U^(1))
 * U^(n+1) = (1/3) * U^n + (2/3) * U^(2) + (2/3) * dt * RHS(U^(2))
 * Order: 3rd order, CFL <= 1, TVD property
 */
class RK3 : public TimeIntegrator {
public:
    void step(
        StateField3D& U_current,
        double dt,
        const RHSFunction& rhs
    ) override;

    int order() const override { return 3; }
    std::string name() const override { return "SSP-RK3"; }
};

} // namespace fvm3d::temporal
