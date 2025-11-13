#include "temporal/rk_integrators.hpp"

namespace fvm3d::temporal {

void RK2::step(
    StateField3D& U_current,
    double dt,
    const RHSFunction& rhs
) {
    allocate_temporaries(U_current);

    int nvars = U_current.nvars();
    int nx = U_current.nx();
    int ny = U_current.ny();
    int nz = U_current.nz();

    // Stage 1: U^(1) = U^n + dt * RHS(U^n)
    rhs(U_current, *temp_rhs_);
    temp_stage1_->assign(U_current);
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    (*temp_stage1_)(v, i, j, k) += dt * (*temp_rhs_)(v, i, j, k);
                }
            }
        }
    }

    // Stage 2: RHS(U^(1))
    rhs(*temp_stage1_, *temp_rhs_);

    // U^(n+1) = 0.5 * U^n + 0.5 * U^(1) + 0.5 * dt * RHS(U^(1))
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    U_current(v, i, j, k) = 0.5 * U_current(v, i, j, k) +
                                           0.5 * (*temp_stage1_)(v, i, j, k) +
                                           0.5 * dt * (*temp_rhs_)(v, i, j, k);
                }
            }
        }
    }
}

void RK3::step(
    StateField3D& U_current,
    double dt,
    const RHSFunction& rhs
) {
    allocate_temporaries(U_current);

    int nvars = U_current.nvars();
    int nx = U_current.nx();
    int ny = U_current.ny();
    int nz = U_current.nz();

    // Stage 1: U^(1) = U^n + dt * RHS(U^n)
    rhs(U_current, *temp_rhs_);
    temp_stage1_->assign(U_current);
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    (*temp_stage1_)(v, i, j, k) += dt * (*temp_rhs_)(v, i, j, k);
                }
            }
        }
    }

    // Stage 2: U^(2) = 0.75 * U^n + 0.25 * U^(1) + 0.25 * dt * RHS(U^(1))
    rhs(*temp_stage1_, *temp_rhs_);
    temp_stage2_->assign(U_current);
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    (*temp_stage2_)(v, i, j, k) = 0.75 * U_current(v, i, j, k) +
                                                 0.25 * (*temp_stage1_)(v, i, j, k) +
                                                 0.25 * dt * (*temp_rhs_)(v, i, j, k);
                }
            }
        }
    }

    // Stage 3: U^(n+1) = (1/3) * U^n + (2/3) * U^(2) + (2/3) * dt * RHS(U^(2))
    rhs(*temp_stage2_, *temp_rhs_);
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    U_current(v, i, j, k) = (1.0/3.0) * U_current(v, i, j, k) +
                                           (2.0/3.0) * (*temp_stage2_)(v, i, j, k) +
                                           (2.0/3.0) * dt * (*temp_rhs_)(v, i, j, k);
                }
            }
        }
    }
}

} // namespace fvm3d::temporal
