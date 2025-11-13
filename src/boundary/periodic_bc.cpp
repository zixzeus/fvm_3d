#include "boundary/periodic_bc.hpp"

namespace fvm3d::boundary {

void PeriodicBC::apply(StateField3D& state, const Grid3D& grid) {
    if (periodic_x_) apply_periodic_x(state, grid);
    if (periodic_y_) apply_periodic_y(state, grid);
    if (periodic_z_) apply_periodic_z(state, grid);
}

void PeriodicBC::apply_periodic_x(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx_total = state.nx();
    int ny = state.ny();
    int nz = state.nz();
    int nghost = grid.nghost();

    // Left ghost: copy from right interior
    for (int v = 0; v < nvars; v++) {
        for (int g = 0; g < nghost; g++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int src_i = nx_total - nghost - 1 - (nghost - 1 - g);
                    state(v, g, j, k) = state(v, src_i, j, k);
                }
            }
        }
    }

    // Right ghost: copy from left interior
    for (int v = 0; v < nvars; v++) {
        for (int g = 0; g < nghost; g++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int src_i = nghost + g;
                    int dst_i = nx_total - nghost + g;
                    state(v, dst_i, j, k) = state(v, src_i, j, k);
                }
            }
        }
    }
}

void PeriodicBC::apply_periodic_y(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx = state.nx();
    int ny_total = state.ny();
    int nz = state.nz();
    int nghost = grid.nghost();

    // Bottom ghost: copy from top interior
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int g = 0; g < nghost; g++) {
                for (int k = 0; k < nz; k++) {
                    int src_j = ny_total - nghost - 1 - (nghost - 1 - g);
                    state(v, i, g, k) = state(v, i, src_j, k);
                }
            }
        }
    }

    // Top ghost: copy from bottom interior
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int g = 0; g < nghost; g++) {
                for (int k = 0; k < nz; k++) {
                    int src_j = nghost + g;
                    int dst_j = ny_total - nghost + g;
                    state(v, i, dst_j, k) = state(v, i, src_j, k);
                }
            }
        }
    }
}

void PeriodicBC::apply_periodic_z(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx = state.nx();
    int ny = state.ny();
    int nz_total = state.nz();
    int nghost = grid.nghost();

    // Bottom ghost: copy from top interior
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int g = 0; g < nghost; g++) {
                    int src_k = nz_total - nghost - 1 - (nghost - 1 - g);
                    state(v, i, j, g) = state(v, i, j, src_k);
                }
            }
        }
    }

    // Top ghost: copy from bottom interior
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int g = 0; g < nghost; g++) {
                    int src_k = nghost + g;
                    int dst_k = nz_total - nghost + g;
                    state(v, i, j, dst_k) = state(v, i, j, src_k);
                }
            }
        }
    }
}

} // namespace fvm3d::boundary
