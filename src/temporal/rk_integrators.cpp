#include "temporal/rk_integrators.hpp"

namespace fvm3d::temporal {

void RK2::step(
    StateField3D& U_current,
    double dt,
    const RHSFunction& rhs
) {
    allocate_temporaries(U_current);

    // Stage 1: U^(1) = U^n + dt * RHS(U^n)
    // Vectorized: temp_stage1 = U_current + dt * temp_rhs
    rhs(U_current, *temp_rhs_);
    temp_stage1_->axpy(dt, *temp_rhs_, U_current);

    // Stage 2: RHS(U^(1))
    rhs(*temp_stage1_, *temp_rhs_);

    // U^(n+1) = 0.5 * U^n + 0.5 * U^(1) + 0.5 * dt * RHS(U^(1))
    // Vectorized: U_current = 0.5 * U_current + 0.5 * temp_stage1 + 0.5*dt * temp_rhs
    U_current.linear_combination_3(
        0.5, U_current,
        0.5, *temp_stage1_,
        0.5 * dt, *temp_rhs_
    );
}

void RK3::step(
    StateField3D& U_current,
    double dt,
    const RHSFunction& rhs
) {
    allocate_temporaries(U_current);

    // Stage 1: U^(1) = U^n + dt * RHS(U^n)
    // Vectorized: temp_stage1 = U_current + dt * temp_rhs
    rhs(U_current, *temp_rhs_);
    temp_stage1_->axpy(dt, *temp_rhs_, U_current);

    // Stage 2: U^(2) = 0.75 * U^n + 0.25 * U^(1) + 0.25 * dt * RHS(U^(1))
    // Vectorized: temp_stage2 = 0.75 * U_current + 0.25 * temp_stage1 + 0.25*dt * temp_rhs
    rhs(*temp_stage1_, *temp_rhs_);
    temp_stage2_->linear_combination_3(
        0.75, U_current,
        0.25, *temp_stage1_,
        0.25 * dt, *temp_rhs_
    );

    // Stage 3: U^(n+1) = (1/3) * U^n + (2/3) * U^(2) + (2/3) * dt * RHS(U^(2))
    // Vectorized: U_current = (1/3) * U_current + (2/3) * temp_stage2 + (2/3)*dt * temp_rhs
    rhs(*temp_stage2_, *temp_rhs_);
    U_current.linear_combination_3(
        1.0/3.0, U_current,
        2.0/3.0, *temp_stage2_,
        2.0/3.0 * dt, *temp_rhs_
    );
}

} // namespace fvm3d::temporal
