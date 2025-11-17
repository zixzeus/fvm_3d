#include "boundary/reflective_bc.hpp"

namespace fvm3d::boundary {

void ReflectiveBC::apply(StateField3D& state, const Grid3D& grid) {
    if (reflect_x_) apply_reflective_x(state, grid);
    if (reflect_y_) apply_reflective_y(state, grid);
    if (reflect_z_) apply_reflective_z(state, grid);
}

void ReflectiveBC::apply_reflective_x(StateField3D& state, const Grid3D& grid) {
    const int nvars = state.nvars();
    const int nx_total = state.nx();
    const int ny = state.ny();
    const int nz = state.nz();
    const int nghost = grid.nghost();
    const int plane_size = ny * nz;

    // Vectorized approach: process each variable separately
    for (int v = 0; v < nvars; v++) {
        const double sign = (v == 1) ? -1.0 : 1.0;  // Flip X-momentum (v=1)

        // Left ghost: mirror from right interior with u -> -u
        for (int g = 0; g < nghost; g++) {
            const int src_i = nghost + (nghost - 1 - g);
            const int dst_i = nghost - 1 - g;
            #pragma omp simd
            for (int idx = 0; idx < plane_size; idx++) {
                int j = idx / nz;
                int k = idx % nz;
                state(v, dst_i, j, k) = sign * state(v, src_i, j, k);
            }
        }

        // Right ghost: mirror from left interior with u -> -u
        for (int g = 0; g < nghost; g++) {
            const int src_i = nx_total - nghost - 1 - (nghost - 1 - g);
            const int dst_i = nx_total - nghost + g;
            #pragma omp simd
            for (int idx = 0; idx < plane_size; idx++) {
                int j = idx / nz;
                int k = idx % nz;
                state(v, dst_i, j, k) = sign * state(v, src_i, j, k);
            }
        }
    }
}

void ReflectiveBC::apply_reflective_y(StateField3D& state, const Grid3D& grid) {
    const int nvars = state.nvars();
    const int nx = state.nx();
    const int ny_total = state.ny();
    const int nz = state.nz();
    const int nghost = grid.nghost();

    // Vectorized approach: process each variable and x-slice separately
    for (int v = 0; v < nvars; v++) {
        const double sign = (v == 2) ? -1.0 : 1.0;  // Flip Y-momentum (v=2)

        for (int i = 0; i < nx; i++) {
            // Bottom ghost: mirror with v -> -v
            for (int g = 0; g < nghost; g++) {
                const int src_j = nghost + (nghost - 1 - g);
                const int dst_j = nghost - 1 - g;
                #pragma omp simd
                for (int k = 0; k < nz; k++) {
                    state(v, i, dst_j, k) = sign * state(v, i, src_j, k);
                }
            }

            // Top ghost: mirror with v -> -v
            for (int g = 0; g < nghost; g++) {
                const int src_j = ny_total - nghost - 1 - (nghost - 1 - g);
                const int dst_j = ny_total - nghost + g;
                #pragma omp simd
                for (int k = 0; k < nz; k++) {
                    state(v, i, dst_j, k) = sign * state(v, i, src_j, k);
                }
            }
        }
    }
}

void ReflectiveBC::apply_reflective_z(StateField3D& state, const Grid3D& grid) {
    const int nvars = state.nvars();
    const int nx = state.nx();
    const int ny = state.ny();
    const int nz_total = state.nz();
    const int nghost = grid.nghost();

    // Vectorized approach: process each variable and (i,j) pair separately
    // Z-direction is innermost in memory, so these are naturally contiguous
    for (int v = 0; v < nvars; v++) {
        const double sign = (v == 3) ? -1.0 : 1.0;  // Flip Z-momentum (v=3)

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                // Bottom ghost: mirror with w -> -w
                for (int g = 0; g < nghost; g++) {
                    const int src_k = nghost + (nghost - 1 - g);
                    const int dst_k = nghost - 1 - g;
                    state(v, i, j, dst_k) = sign * state(v, i, j, src_k);
                }

                // Top ghost: mirror with w -> -w
                for (int g = 0; g < nghost; g++) {
                    const int src_k = nz_total - nghost - 1 - (nghost - 1 - g);
                    const int dst_k = nz_total - nghost + g;
                    state(v, i, j, dst_k) = sign * state(v, i, j, src_k);
                }
            }
        }
    }
}

} // namespace fvm3d::boundary
