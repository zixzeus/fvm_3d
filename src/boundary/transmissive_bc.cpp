#include "boundary/transmissive_bc.hpp"

namespace fvm3d::boundary {

void TransmissiveBC::apply(StateField3D& state, const Grid3D& grid) {
    if (transmissive_x_) apply_transmissive_x(state, grid);
    if (transmissive_y_) apply_transmissive_y(state, grid);
    if (transmissive_z_) apply_transmissive_z(state, grid);
}

void TransmissiveBC::apply_transmissive_x(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx_total = state.nx();
    int ny = state.ny();
    int nz = state.nz();
    int nghost = grid.nghost();

    // Left ghost: copy from first interior cell
    for (int v = 0; v < nvars; v++) {
        for (int g = 0; g < nghost; g++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    state(v, g, j, k) = state(v, nghost, j, k);
                }
            }
        }
    }

    // Right ghost: copy from last interior cell
    for (int v = 0; v < nvars; v++) {
        for (int g = 0; g < nghost; g++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int src_i = nx_total - nghost - 1;
                    int dst_i = nx_total - nghost + g;
                    state(v, dst_i, j, k) = state(v, src_i, j, k);
                }
            }
        }
    }
}

void TransmissiveBC::apply_transmissive_y(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx = state.nx();
    int ny_total = state.ny();
    int nz = state.nz();
    int nghost = grid.nghost();

    // Bottom ghost: copy from first interior row
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int g = 0; g < nghost; g++) {
                for (int k = 0; k < nz; k++) {
                    state(v, i, g, k) = state(v, i, nghost, k);
                }
            }
        }
    }

    // Top ghost: copy from last interior row
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int g = 0; g < nghost; g++) {
                for (int k = 0; k < nz; k++) {
                    int src_j = ny_total - nghost - 1;
                    int dst_j = ny_total - nghost + g;
                    state(v, i, dst_j, k) = state(v, i, src_j, k);
                }
            }
        }
    }
}

void TransmissiveBC::apply_transmissive_z(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx = state.nx();
    int ny = state.ny();
    int nz_total = state.nz();
    int nghost = grid.nghost();

    // Bottom ghost: copy from first interior layer
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int g = 0; g < nghost; g++) {
                    state(v, i, j, g) = state(v, i, j, nghost);
                }
            }
        }
    }

    // Top ghost: copy from last interior layer
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int g = 0; g < nghost; g++) {
                    int src_k = nz_total - nghost - 1;
                    int dst_k = nz_total - nghost + g;
                    state(v, i, j, dst_k) = state(v, i, j, src_k);
                }
            }
        }
    }
}

} // namespace fvm3d::boundary
