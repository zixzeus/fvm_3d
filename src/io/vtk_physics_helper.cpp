#include "io/vtk_physics_helper.hpp"
#include <cmath>

namespace fvm3d::io {

void VTKPhysicsHelper::compute_euler_primitives(
    const core::StateField3D& state,
    const core::Grid3D& grid,
    std::vector<double>& pressure,
    std::vector<double>& vel_x,
    std::vector<double>& vel_y,
    std::vector<double>& vel_z,
    std::vector<double>* temperature,
    double gamma,
    double R
) {
    int nx = grid.nx_local();
    int ny = grid.ny_local();
    int nz = grid.nz_local();
    int ng = grid.nghost();

    // Resize output vectors
    int total_cells = nx * ny * nz;
    pressure.resize(total_cells);
    vel_x.resize(total_cells);
    vel_y.resize(total_cells);
    vel_z.resize(total_cells);
    if (temperature) {
        temperature->resize(total_cells);
    }

    // Compute primitive variables from conservative state
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int idx = k + nz * (j + ny * i);

                double rho = state(0, i + ng, j + ng, k + ng);
                double rho_u = state(1, i + ng, j + ng, k + ng);
                double rho_v = state(2, i + ng, j + ng, k + ng);
                double rho_w = state(3, i + ng, j + ng, k + ng);
                double E = state(4, i + ng, j + ng, k + ng);

                if (rho > MIN_DENSITY) {
                    vel_x[idx] = rho_u / rho;
                    vel_y[idx] = rho_v / rho;
                    vel_z[idx] = rho_w / rho;

                    double ke = 0.5 * (rho_u * rho_u + rho_v * rho_v + rho_w * rho_w) / rho;
                    pressure[idx] = (gamma - 1.0) * (E - ke);

                    if (temperature) {
                        (*temperature)[idx] = pressure[idx] / (rho * R);
                    }
                } else {
                    vel_x[idx] = 0.0;
                    vel_y[idx] = 0.0;
                    vel_z[idx] = 0.0;
                    pressure[idx] = 0.0;
                    if (temperature) {
                        (*temperature)[idx] = 0.0;
                    }
                }
            }
        }
    }
}

void VTKPhysicsHelper::compute_mhd_primitives(
    const core::StateField3D& state,
    const core::Grid3D& grid,
    std::vector<double>& pressure,
    std::vector<double>& vel_x,
    std::vector<double>& vel_y,
    std::vector<double>& vel_z,
    std::vector<double>& B_magnitude,
    double gamma
) {
    int nx = grid.nx_local();
    int ny = grid.ny_local();
    int nz = grid.nz_local();
    int ng = grid.nghost();

    // Resize output vectors
    int total_cells = nx * ny * nz;
    pressure.resize(total_cells);
    vel_x.resize(total_cells);
    vel_y.resize(total_cells);
    vel_z.resize(total_cells);
    B_magnitude.resize(total_cells);

    // Compute primitive variables from conservative MHD state
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int idx = k + nz * (j + ny * i);

                double rho = state(0, i + ng, j + ng, k + ng);
                double rho_u = state(1, i + ng, j + ng, k + ng);
                double rho_v = state(2, i + ng, j + ng, k + ng);
                double rho_w = state(3, i + ng, j + ng, k + ng);
                double E = state(4, i + ng, j + ng, k + ng);
                double Bx = state(5, i + ng, j + ng, k + ng);
                double By = state(6, i + ng, j + ng, k + ng);
                double Bz = state(7, i + ng, j + ng, k + ng);

                double B2 = Bx * Bx + By * By + Bz * Bz;
                B_magnitude[idx] = std::sqrt(B2);

                if (rho > MIN_DENSITY) {
                    vel_x[idx] = rho_u / rho;
                    vel_y[idx] = rho_v / rho;
                    vel_z[idx] = rho_w / rho;

                    double ke = 0.5 * (rho_u * rho_u + rho_v * rho_v + rho_w * rho_w) / rho;
                    pressure[idx] = (gamma - 1.0) * (E - ke - 0.5 * B2);
                } else {
                    vel_x[idx] = 0.0;
                    vel_y[idx] = 0.0;
                    vel_z[idx] = 0.0;
                    pressure[idx] = 0.0;
                }
            }
        }
    }
}

void VTKPhysicsHelper::extract_magnetic_field(
    const core::StateField3D& state,
    const core::Grid3D& grid,
    std::vector<double>& Bx,
    std::vector<double>& By,
    std::vector<double>& Bz
) {
    int nx = grid.nx_local();
    int ny = grid.ny_local();
    int nz = grid.nz_local();
    int ng = grid.nghost();

    // Resize output vectors
    int total_cells = nx * ny * nz;
    Bx.resize(total_cells);
    By.resize(total_cells);
    Bz.resize(total_cells);

    // Extract magnetic field components
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int idx = k + nz * (j + ny * i);
                Bx[idx] = state(5, i + ng, j + ng, k + ng);
                By[idx] = state(6, i + ng, j + ng, k + ng);
                Bz[idx] = state(7, i + ng, j + ng, k + ng);
            }
        }
    }
}

} // namespace fvm3d::io
