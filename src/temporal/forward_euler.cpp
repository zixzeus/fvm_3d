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
    int nvars = U_current.nvars();
    int nx = U_current.nx();
    int ny = U_current.ny();
    int nz = U_current.nz();

    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    U_current(v, i, j, k) += dt * (*temp_rhs_)(v, i, j, k);
                }
            }
        }
    }
}

} // namespace fvm3d::temporal
