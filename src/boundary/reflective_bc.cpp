#include "boundary/reflective_bc.hpp"

namespace fvm3d::boundary {

void ReflectiveBC::apply(StateField3D& state, const Grid3D& grid) {
    if (reflect_x_) apply_reflective_x(state, grid);
    if (reflect_y_) apply_reflective_y(state, grid);
    if (reflect_z_) apply_reflective_z(state, grid);
}

void ReflectiveBC::apply_reflective_x(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx_total = state.nx();
    int ny = state.ny();
    int nz = state.nz();
    int nghost = grid.nghost();

    // Left ghost: mirror from right interior with u -> -u
    for (int v = 0; v < nvars; v++) {
        for (int g = 0; g < nghost; g++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int src_i = nghost + (nghost - 1 - g);
                    int dst_i = nghost - 1 - g;

                    if (v == 1) {  // rho_u component (X-momentum)
                        state(v, dst_i, j, k) = -state(v, src_i, j, k);
                    } else {
                        state(v, dst_i, j, k) = state(v, src_i, j, k);
                    }
                }
            }
        }
    }

    // Right ghost: mirror from left interior with u -> -u
    for (int v = 0; v < nvars; v++) {
        for (int g = 0; g < nghost; g++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int src_i = nx_total - nghost - 1 - (nghost - 1 - g);
                    int dst_i = nx_total - nghost + g;

                    if (v == 1) {  // rho_u component (X-momentum)
                        state(v, dst_i, j, k) = -state(v, src_i, j, k);
                    } else {
                        state(v, dst_i, j, k) = state(v, src_i, j, k);
                    }
                }
            }
        }
    }
}

void ReflectiveBC::apply_reflective_y(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx = state.nx();
    int ny_total = state.ny();
    int nz = state.nz();
    int nghost = grid.nghost();

    // Bottom ghost: mirror with v -> -v
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int g = 0; g < nghost; g++) {
                for (int k = 0; k < nz; k++) {
                    int src_j = nghost + (nghost - 1 - g);
                    int dst_j = nghost - 1 - g;

                    if (v == 2) {  // rho_v component (Y-momentum)
                        state(v, i, dst_j, k) = -state(v, i, src_j, k);
                    } else {
                        state(v, i, dst_j, k) = state(v, i, src_j, k);
                    }
                }
            }
        }
    }

    // Top ghost: mirror with v -> -v
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int g = 0; g < nghost; g++) {
                for (int k = 0; k < nz; k++) {
                    int src_j = ny_total - nghost - 1 - (nghost - 1 - g);
                    int dst_j = ny_total - nghost + g;

                    if (v == 2) {  // rho_v component (Y-momentum)
                        state(v, i, dst_j, k) = -state(v, i, src_j, k);
                    } else {
                        state(v, i, dst_j, k) = state(v, i, src_j, k);
                    }
                }
            }
        }
    }
}

void ReflectiveBC::apply_reflective_z(StateField3D& state, const Grid3D& grid) {
    int nvars = state.nvars();
    int nx = state.nx();
    int ny = state.ny();
    int nz_total = state.nz();
    int nghost = grid.nghost();

    // Bottom ghost: mirror with w -> -w
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int g = 0; g < nghost; g++) {
                    int src_k = nghost + (nghost - 1 - g);
                    int dst_k = nghost - 1 - g;

                    if (v == 3) {  // rho_w component (Z-momentum)
                        state(v, i, j, dst_k) = -state(v, i, j, src_k);
                    } else {
                        state(v, i, j, dst_k) = state(v, i, j, src_k);
                    }
                }
            }
        }
    }

    // Top ghost: mirror with w -> -w
    for (int v = 0; v < nvars; v++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int g = 0; g < nghost; g++) {
                    int src_k = nz_total - nghost - 1 - (nghost - 1 - g);
                    int dst_k = nz_total - nghost + g;

                    if (v == 3) {  // rho_w component (Z-momentum)
                        state(v, i, j, dst_k) = -state(v, i, j, src_k);
                    } else {
                        state(v, i, j, dst_k) = state(v, i, j, src_k);
                    }
                }
            }
        }
    }
}

} // namespace fvm3d::boundary
