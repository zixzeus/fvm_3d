#pragma once

#include "core/field3d.hpp"
#include <functional>
#include <memory>
#include <string>

namespace fvm3d::temporal {

using StateField3D = fvm3d::core::StateField3D;

/**
 * Abstract base class for time integration schemes.
 * Using method-of-lines: dU/dt = RHS(U)
 */
class TimeIntegrator {
public:
    using RHSFunction = std::function<void(const StateField3D&, StateField3D&)>;

    virtual ~TimeIntegrator() = default;

    /**
     * Perform one time step.
     */
    virtual void step(
        StateField3D& U_current,
        double dt,
        const RHSFunction& rhs
    ) = 0;

    /**
     * Get the order of accuracy in time.
     */
    virtual int order() const = 0;

    /**
     * Get the name of the integrator.
     */
    virtual std::string name() const = 0;

protected:
    mutable std::unique_ptr<StateField3D> temp_stage1_;
    mutable std::unique_ptr<StateField3D> temp_stage2_;
    mutable std::unique_ptr<StateField3D> temp_rhs_;

    /**
     * Allocate temporary arrays if needed.
     */
    void allocate_temporaries(const StateField3D& U) const {
        if (!temp_rhs_) {
            temp_rhs_ = std::make_unique<StateField3D>(U.nvars(), U.nx(), U.ny(), U.nz());
            temp_stage1_ = std::make_unique<StateField3D>(U.nvars(), U.nx(), U.ny(), U.nz());
            temp_stage2_ = std::make_unique<StateField3D>(U.nvars(), U.nx(), U.ny(), U.nz());
        }
    }
};

} // namespace fvm3d::temporal
