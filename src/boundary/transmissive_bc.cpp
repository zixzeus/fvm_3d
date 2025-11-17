#include "boundary/transmissive_bc.hpp"
#include <cstring>

namespace fvm3d::boundary {

void TransmissiveBC::apply(StateField3D& state, const Grid3D& grid) {
    if (transmissive_x_) apply_transmissive_x(state, grid);
    if (transmissive_y_) apply_transmissive_y(state, grid);
    if (transmissive_z_) apply_transmissive_z(state, grid);
}

void TransmissiveBC::apply_transmissive_x(StateField3D& state, const Grid3D& grid) {
    const int nvars = state.nvars();
    const int nx_total = state.nx();
    const int ny = state.ny();
    const int nz = state.nz();
    const int nghost = grid.nghost();
    const int plane_size = ny * nz;  // Size of one y-z plane

    // Vectorized approach: process each variable separately
    // For each variable, ghost cells form contiguous memory blocks in SoA layout

    for (int v = 0; v < nvars; v++) {
        // Left ghost layers: copy from first interior plane (i = nghost)
        // Each ghost layer g gets the same values from i = nghost
        for (int g = 0; g < nghost; g++) {
            #pragma omp simd
            for (int idx = 0; idx < plane_size; idx++) {
                int j = idx / nz;
                int k = idx % nz;
                state(v, g, j, k) = state(v, nghost, j, k);
            }
        }

        // Right ghost layers: copy from last interior plane
        const int src_i = nx_total - nghost - 1;
        for (int g = 0; g < nghost; g++) {
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

void TransmissiveBC::apply_transmissive_y(StateField3D& state, const Grid3D& grid) {
    const int nvars = state.nvars();
    const int nx = state.nx();
    const int ny_total = state.ny();
    const int nz = state.nz();
    const int nghost = grid.nghost();

    // Vectorized approach: process each variable and x-slice separately
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            // Bottom ghost layers: copy from first interior row (j = nghost)
            for (int g = 0; g < nghost; g++) {
                #pragma omp simd
                for (int k = 0; k < nz; k++) {
                    state(v, i, g, k) = state(v, i, nghost, k);
                }
            }

            // Top ghost layers: copy from last interior row
            const int src_j = ny_total - nghost - 1;
            for (int g = 0; g < nghost; g++) {
                const int dst_j = ny_total - nghost + g;
                #pragma omp simd
                for (int k = 0; k < nz; k++) {
                    state(v, i, dst_j, k) = state(v, i, src_j, k);
                }
            }
        }
    }
}

void TransmissiveBC::apply_transmissive_z(StateField3D& state, const Grid3D& grid) {
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
                // Bottom ghost layers: copy from first interior layer (k = nghost)
                const double src_value = state(v, i, j, nghost);
                #pragma omp simd
                for (int g = 0; g < nghost; g++) {
                    state(v, i, j, g) = src_value;
                }

                // Top ghost layers: copy from last interior layer
                const double src_value_top = state(v, i, j, nz_total - nghost - 1);
                #pragma omp simd
                for (int g = 0; g < nghost; g++) {
                    state(v, i, j, nz_total - nghost + g) = src_value_top;
                }
            }
        }
    }
}

} // namespace fvm3d::boundary
