#include "boundary/periodic_bc.hpp"

namespace fvm3d::boundary {

void PeriodicBC::apply(StateField3D& state, const Grid3D& grid) {
    if (periodic_x_) apply_periodic_x(state, grid);
    if (periodic_y_) apply_periodic_y(state, grid);
    if (periodic_z_) apply_periodic_z(state, grid);
}

void PeriodicBC::apply_periodic_x(StateField3D& state, const Grid3D& grid) {
    const int nvars = state.nvars();
    const int nx_total = state.nx();
    const int ny = state.ny();
    const int nz = state.nz();
    const int nghost = grid.nghost();
    const int plane_size = ny * nz;

    // Vectorized approach: process each variable separately
    for (int v = 0; v < nvars; v++) {
        // Left ghost: copy from right interior
        for (int g = 0; g < nghost; g++) {
            const int src_i = nx_total - nghost - 1 - (nghost - 1 - g);
            const int dst_i = g;
            #pragma omp simd
            for (int idx = 0; idx < plane_size; idx++) {
                int j = idx / nz;
                int k = idx % nz;
                state(v, dst_i, j, k) = state(v, src_i, j, k);
            }
        }

        // Right ghost: copy from left interior
        for (int g = 0; g < nghost; g++) {
            const int src_i = nghost + g;
            const int dst_i = nx_total - nghost + g;
            #pragma omp simd
            for (int idx = 0; idx < plane_size; idx++) {
                int j = idx / nz;
                int k = idx % nz;
                state(v, dst_i, j, k) = state(v, src_i, j, k);
            }
        }
    }
}

void PeriodicBC::apply_periodic_y(StateField3D& state, const Grid3D& grid) {
    const int nvars = state.nvars();
    const int nx = state.nx();
    const int ny_total = state.ny();
    const int nz = state.nz();
    const int nghost = grid.nghost();

    // Vectorized approach: process each variable and x-slice separately
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            // Bottom ghost: copy from top interior
            for (int g = 0; g < nghost; g++) {
                const int src_j = ny_total - nghost - 1 - (nghost - 1 - g);
                const int dst_j = g;
                #pragma omp simd
                for (int k = 0; k < nz; k++) {
                    state(v, i, dst_j, k) = state(v, i, src_j, k);
                }
            }

            // Top ghost: copy from bottom interior
            for (int g = 0; g < nghost; g++) {
                const int src_j = nghost + g;
                const int dst_j = ny_total - nghost + g;
                #pragma omp simd
                for (int k = 0; k < nz; k++) {
                    state(v, i, dst_j, k) = state(v, i, src_j, k);
                }
            }
        }
    }
}

void PeriodicBC::apply_periodic_z(StateField3D& state, const Grid3D& grid) {
    const int nvars = state.nvars();
    const int nx = state.nx();
    const int ny = state.ny();
    const int nz_total = state.nz();
    const int nghost = grid.nghost();

    // Vectorized approach: process each variable and (i,j) pair separately
    // Z-direction is innermost in memory, so these are naturally contiguous
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                // Bottom ghost: copy from top interior
                #pragma omp simd
                for (int g = 0; g < nghost; g++) {
                    const int src_k = nz_total - nghost - 1 - (nghost - 1 - g);
                    state(v, i, j, g) = state(v, i, j, src_k);
                }

                // Top ghost: copy from bottom interior
                #pragma omp simd
                for (int g = 0; g < nghost; g++) {
                    const int src_k = nghost + g;
                    const int dst_k = nz_total - nghost + g;
                    state(v, i, j, dst_k) = state(v, i, j, src_k);
                }
            }
        }
    }
}

} // namespace fvm3d::boundary
