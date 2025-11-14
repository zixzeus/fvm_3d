#include "temporal/forward_euler.hpp"

namespace fvm3d::temporal {

void ForwardEuler::step(
    StateField3D& U_current,
    double dt,
    const RHSFunction& rhs
) {
    allocate_temporaries(U_current);

    // Compute RHS(U^n)
    rhs(U_current, *temp_rhs_);

    // U^(n+1) = U^n + dt * RHS(U^n)
    // Vectorized in-place AXPY: U_current = U_current + dt * temp_rhs
    U_current.add_scaled(dt, *temp_rhs_);
}

} // namespace fvm3d::temporal
