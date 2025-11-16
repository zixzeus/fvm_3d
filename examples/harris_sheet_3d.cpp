#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
#include "physics/resistive_mhd3d_advanced.hpp"

using namespace fvm3d::physics;

/**
 * 3D Harris Sheet Magnetic Reconnection Simulation
 *
 * Single-process (serial) implementation for testing and validation.
 * This demonstrates how to:
 * - Initialize the Harris sheet equilibrium
 * - Compute fluxes and sources
 * - Perform time stepping
 * - Monitor energy and magnetic divergence
 */

// Grid parameters
struct GridConfig {
    int nx = 64;      // Resolution in X
    int ny = 32;      // Resolution in Y
    int nz = 32;      // Resolution in Z
    double Lx = 10.0; // Domain size X
    double Ly = 5.0;  // Domain size Y
    double Lz = 4.0;  // Domain size Z
    double dx() const { return Lx / nx; }
    double dy() const { return Ly / ny; }
    double dz() const { return Lz / nz; }
};

// Simulation parameters
struct SimConfig {
    int nsteps = 100;        // Number of time steps
    double CFL = 0.35;       // CFL number
    int diag_freq = 10;      // Diagnostic output frequency
    double dt_resistive_factor = 0.25;  // Parabolic stability factor
};

// 3D state array: state[variable, x, y, z]
class State3D {
public:
    State3D(int nx, int ny, int nz)
        : nx_(nx), ny_(ny), nz_(nz) {
        // Allocate for 9 variables
        data_.resize(9, Eigen::MatrixXd::Zero(nx*ny, nz));
    }

    double& operator()(int v, int i, int j, int k) {
        return data_[v](i*ny_ + j, k);
    }

    const double& operator()(int v, int i, int j, int k) const {
        return data_[v](i*ny_ + j, k);
    }

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }

private:
    int nx_, ny_, nz_;
    std::vector<Eigen::MatrixXd> data_;  // [9 variables]
};

/**
 * Initialize Harris sheet with perturbation
 */
void initialize_harris_sheet(
    State3D& U,
    const GridConfig& grid,
    const AdvancedResistiveMHD3D& mhd,
    const AdvancedResistiveMHD3D::HarrisSheetConfig& harris
) {
    double dx = grid.dx();
    double dy = grid.dy();
    double dz = grid.dz();

    for(int i = 0; i < grid.nx; i++) {
        for(int j = 0; j < grid.ny; j++) {
            for(int k = 0; k < grid.nz; k++) {
                // Map to physical coordinates (centered domain)
                double x = (i + 0.5) * dx - grid.Lx/2.0;
                double y = (j + 0.5) * dy - grid.Ly/2.0;
                double z = (k + 0.5) * dz - grid.Lz/2.0;

                // Get Harris sheet state
                Eigen::VectorXd U_ij = mhd.harris_sheet_initial(x, y, z, harris);

                // Copy to 3D array
                for(int v = 0; v < 9; v++) {
                    U(v, i, j, k) = U_ij(v);
                }
            }
        }
    }
}

/**
 * Compute maximum wave speed for CFL condition
 */
double compute_max_wave_speed(
    const State3D& U,
    const GridConfig& grid,
    const AdvancedResistiveMHD3D& mhd
) {
    double max_speed = 0.0;
    Eigen::VectorXd U_vec(9);

    for(int i = 0; i < grid.nx; i++) {
        for(int j = 0; j < grid.ny; j++) {
            for(int k = 0; k < grid.nz; k++) {
                // Pack state vector
                for(int v = 0; v < 9; v++) {
                    U_vec(v) = U(v, i, j, k);
                }

                // Compute max wave speed in each direction
                double s_x = mhd.max_wave_speed(U_vec, 0);
                double s_y = mhd.max_wave_speed(U_vec, 1);
                double s_z = mhd.max_wave_speed(U_vec, 2);

                max_speed = std::max({max_speed, s_x, s_y, s_z});
            }
        }
    }

    return max_speed;
}

/**
 * Compute time step with both hyperbolic CFL and parabolic stability
 */
double compute_timestep(
    const State3D& U,
    const GridConfig& grid,
    const SimConfig& sim,
    const AdvancedResistiveMHD3D& mhd
) {
    // Hyperbolic CFL
    double max_speed = compute_max_wave_speed(U, grid, mhd);
    double min_dx = std::min({grid.dx(), grid.dy(), grid.dz()});
    double dt_hyperbolic = sim.CFL * min_dx / max_speed;

    // Parabolic (resistive diffusion) stability
    // dt < 0.25 * min(dx²,dy²,dz²) / eta_max
    double eta_max = 0.01667;  // eta1
    double min_dx2 = std::min({grid.dx()*grid.dx(),
                               grid.dy()*grid.dy(),
                               grid.dz()*grid.dz()});
    double dt_parabolic = sim.dt_resistive_factor * min_dx2 / eta_max;

    // Use most restrictive
    double dt = std::min(dt_hyperbolic, dt_parabolic);

    if(dt <= 0) {
        std::cerr << "ERROR: Invalid timestep!" << std::endl;
        dt = 1e-6;
    }

    return dt;
}

/**
 * Compute diagnostics
 */
struct Diagnostics {
    double kinetic_energy = 0.0;
    double magnetic_energy = 0.0;
    double internal_energy = 0.0;
    double total_energy = 0.0;
    double max_div_B = 0.0;
    double max_rho = 0.0;
    double min_rho = 1e10;
    double max_By = 0.0;  // Reconnected field
};

Diagnostics compute_diagnostics(
    const State3D& U,
    const GridConfig& grid,
    const AdvancedResistiveMHD3D& mhd
) {
    Diagnostics diag;
    Eigen::VectorXd U_vec(9);
    double dx = grid.dx();
    double dy = grid.dy();
    double dz = grid.dz();

    for(int i = 0; i < grid.nx; i++) {
        for(int j = 0; j < grid.ny; j++) {
            for(int k = 0; k < grid.nz; k++) {
                // Get state
                for(int v = 0; v < 9; v++) {
                    U_vec(v) = U(v, i, j, k);
                }

                Eigen::VectorXd V_vec = mhd.conservative_to_primitive(U_vec);
                double rho = V_vec(0);
                double u = V_vec(1);
                double v = V_vec(2);
                double w = V_vec(3);
                double p = V_vec(4);
                double Bx = V_vec(5);
                double By = V_vec(6);
                double Bz = V_vec(7);

                // Accumulate energies
                auto E = mhd.compute_energies(U_vec);
                diag.kinetic_energy += E(0) * dx * dy * dz;
                diag.magnetic_energy += E(1) * dx * dy * dz;
                diag.internal_energy += E(2) * dx * dy * dz;

                // Track max/min
                diag.max_rho = std::max(diag.max_rho, rho);
                diag.min_rho = std::min(diag.min_rho, rho);
                diag.max_By = std::max(diag.max_By, std::abs(By));
            }
        }
    }

    diag.total_energy = diag.kinetic_energy + diag.magnetic_energy + diag.internal_energy;

    // Estimate div B using central differences (sample)
    for(int i = 1; i < grid.nx-1; i++) {
        for(int j = 1; j < grid.ny-1; j++) {
            for(int k = 1; k < grid.nz-1; k++) {
                double dBx = (U(5, i+1, j, k) - U(5, i-1, j, k)) / (2.0 * dx);
                double dBy = (U(6, i, j+1, k) - U(6, i, j-1, k)) / (2.0 * dy);
                double dBz = (U(7, i, j, k+1) - U(7, i, j, k-1)) / (2.0 * dz);
                double div_B = std::abs(dBx + dBy + dBz);
                diag.max_div_B = std::max(diag.max_div_B, div_B);
            }
        }
    }

    return diag;
}

/**
 * Simple Euler time step (for demonstration)
 * Note: This is explicit and first-order for simplicity.
 * Production code would use TVD RK2 or higher.
 */
void euler_step(
    State3D& U_new,
    const State3D& U_old,
    double dt,
    const GridConfig& grid,
    const AdvancedResistiveMHD3D& mhd
) {
    double dx = grid.dx();
    double dy = grid.dy();
    double dz = grid.dz();
    Eigen::VectorXd U_vec(9), F(9), G(9), H(9), S(9);

    // For simplicity, only compute interior cells (skip boundaries)
    for(int i = 1; i < grid.nx-1; i++) {
        for(int j = 1; j < grid.ny-1; j++) {
            for(int k = 1; k < grid.nz-1; k++) {
                double x = (i + 0.5) * dx - grid.Lx/2.0;
                double y = (j + 0.5) * dy - grid.Ly/2.0;
                double z = (k + 0.5) * dz - grid.Lz/2.0;

                // Pack current state
                for(int v = 0; v < 9; v++) {
                    U_vec(v) = U_old(v, i, j, k);
                }

                // Compute fluxes
                F = mhd.flux_x(U_vec, x, y, z);
                G = mhd.flux_y(U_vec, x, y, z);
                H = mhd.flux_z(U_vec, x, y, z);

                // Get neighbor states
                Eigen::VectorXd U_xm(9), U_xp(9), U_ym(9), U_yp(9), U_zm(9), U_zp(9);
                for(int v = 0; v < 9; v++) {
                    U_xm(v) = U_old(v, i-1, j, k);
                    U_xp(v) = U_old(v, i+1, j, k);
                    U_ym(v) = U_old(v, i, j-1, k);
                    U_yp(v) = U_old(v, i, j+1, k);
                    U_zm(v) = U_old(v, i, j, k-1);
                    U_zp(v) = U_old(v, i, j, k+1);
                }

                // Compute fluxes at neighbors (simplified: use central differences)
                Eigen::VectorXd F_xm = mhd.flux_x(U_xm, x - dx, y, z);
                Eigen::VectorXd F_xp = mhd.flux_x(U_xp, x + dx, y, z);
                Eigen::VectorXd G_ym = mhd.flux_y(U_ym, x, y - dy, z);
                Eigen::VectorXd G_yp = mhd.flux_y(U_yp, x, y + dy, z);
                Eigen::VectorXd H_zm = mhd.flux_z(U_zm, x, y, z - dz);
                Eigen::VectorXd H_zp = mhd.flux_z(U_zp, x, y, z + dz);

                // Divergence of fluxes
                Eigen::VectorXd dF = (F_xp - F_xm) / (2.0 * dx);
                Eigen::VectorXd dG = (G_yp - G_ym) / (2.0 * dy);
                Eigen::VectorXd dH = (H_zp - H_zm) / (2.0 * dz);

                // Compute source term
                S = mhd.resistive_source(U_vec, x, y, z, dx, dy, dz,
                                         U_xm, U_xp, U_ym, U_yp, U_zm, U_zp);

                // Update state
                for(int v = 0; v < 9; v++) {
                    U_new(v, i, j, k) = U_old(v, i, j, k) -
                                       dt * (dF(v) + dG(v) + dH(v)) +
                                       dt * S(v);
                }
            }
        }
    }

    // Copy boundaries (periodic or reflective)
    // For simplicity, just copy old state at boundaries
    for(int i = 0; i < grid.nx; i++) {
        for(int j = 0; j < grid.ny; j++) {
            for(int k = 0; k < grid.nz; k++) {
                if(i == 0 || i == grid.nx-1 || j == 0 || j == grid.ny-1 ||
                   k == 0 || k == grid.nz-1) {
                    for(int v = 0; v < 9; v++) {
                        U_new(v, i, j, k) = U_old(v, i, j, k);
                    }
                }
            }
        }
    }
}

/**
 * Main simulation
 */
int main(int argc, char* argv[]) {
    std::cout << "════════════════════════════════════════════════════" << std::endl;
    std::cout << "  3D Harris Sheet Magnetic Reconnection Simulation" << std::endl;
    std::cout << "  (Serial Version - Single Process)" << std::endl;
    std::cout << "════════════════════════════════════════════════════" << std::endl << std::endl;

    // Configuration
    GridConfig grid;
    if(argc > 1) grid.nx = std::atoi(argv[1]);
    if(argc > 2) grid.ny = std::atoi(argv[2]);
    if(argc > 3) grid.nz = std::atoi(argv[3]);

    SimConfig sim;

    std::cout << "Grid Configuration:" << std::endl;
    std::cout << "  Resolution:  " << grid.nx << " × " << grid.ny << " × " << grid.nz << std::endl;
    std::cout << "  Domain:      " << grid.Lx << " × " << grid.Ly << " × " << grid.Lz << std::endl;
    std::cout << "  Cell size:   " << grid.dx() << " × " << grid.dy() << " × " << grid.dz() << std::endl;
    std::cout << "  Total cells: " << (grid.nx * grid.ny * grid.nz) << std::endl << std::endl;

    std::cout << "Physics Configuration:" << std::endl;
    std::cout << "  CFL number:           " << sim.CFL << std::endl;
    std::cout << "  Resistivity model:    Position-dependent" << std::endl;
    std::cout << "    η₀ (background):    " << 1e-3 << " (Rm₀ = 1000)" << std::endl;
    std::cout << "    η₁ (enhanced):      " << 0.01667 << " (Rm₁ = 60)" << std::endl;
    std::cout << "  GLM parameters:       ch=0.2, cr=0.2" << std::endl << std::endl;

    // Create physics
    AdvancedResistiveMHD3D::ResistivityModel resistivity;
    resistivity.eta0 = 1e-3;
    resistivity.eta1 = 0.01667;
    resistivity.localization_scale = 1.0;

    AdvancedResistiveMHD3D::GLMParameters glm(0.2, 0.2);
    AdvancedResistiveMHD3D mhd(resistivity, glm);

    // Harris sheet configuration
    AdvancedResistiveMHD3D::HarrisSheetConfig harris;
    harris.L_sheet = 1.0;
    harris.n0 = 1.0;
    harris.p0 = 0.1;
    harris.B0 = 1.0;
    harris.beta = 0.2;
    harris.perturbation_amplitude = 0.03;

    std::cout << "Harris Sheet Equilibrium:" << std::endl;
    std::cout << "  Sheet thickness:      " << harris.L_sheet << std::endl;
    std::cout << "  Plasma beta:          " << harris.beta << std::endl;
    std::cout << "  Perturbation (m=1):   " << (harris.perturbation_amplitude*100) << "%" << std::endl << std::endl;

    // Initialize
    std::cout << "Initializing Harris sheet equilibrium..." << std::endl;
    State3D U_old(grid.nx, grid.ny, grid.nz);
    State3D U_new(grid.nx, grid.ny, grid.nz);
    initialize_harris_sheet(U_old, grid, mhd, harris);
    std::cout << "✓ Initialization complete" << std::endl << std::endl;

    // Time stepping
    std::cout << "Starting time integration..." << std::endl;
    std::cout << "Step | Time    | dt      | KE      | BE      | IE      | max|B_y| | div(B)max" << std::endl;
    std::cout << "─────┼─────────┼─────────┼─────────┼─────────┼─────────┼──────────┼──────────" << std::endl;

    double time = 0.0;
    for(int n = 0; n < sim.nsteps; n++) {
        // Compute timestep
        double dt = compute_timestep(U_old, grid, sim, mhd);
        time += dt;

        // Step
        euler_step(U_new, U_old, dt, grid, mhd);

        // Swap states
        std::swap(U_old, U_new);

        // Diagnostics every diag_freq steps
        if(n % sim.diag_freq == 0) {
            Diagnostics diag = compute_diagnostics(U_old, grid, mhd);

            std::cout << std::setw(4) << n << " | "
                      << std::scientific << std::setprecision(2)
                      << std::setw(7) << time << " | "
                      << std::setw(7) << dt << " | "
                      << std::setw(7) << diag.kinetic_energy << " | "
                      << std::setw(7) << diag.magnetic_energy << " | "
                      << std::setw(7) << diag.internal_energy << " | "
                      << std::setw(8) << diag.max_By << " | "
                      << std::setw(8) << diag.max_div_B
                      << std::endl;
        }
    }

    std::cout << "═════╧═════════╧═════════╧═════════╧═════════╧═════════╧══════════╧══════════" << std::endl;
    std::cout << std::endl << "✓ Simulation complete!" << std::endl;

    // Final diagnostics
    Diagnostics final_diag = compute_diagnostics(U_old, grid, mhd);
    std::cout << std::endl << "Final Diagnostics:" << std::endl;
    std::cout << "  Kinetic energy:       " << std::scientific << std::setprecision(4)
              << final_diag.kinetic_energy << std::endl;
    std::cout << "  Magnetic energy:      " << final_diag.magnetic_energy << std::endl;
    std::cout << "  Internal energy:      " << final_diag.internal_energy << std::endl;
    std::cout << "  Total energy:         " << final_diag.total_energy << std::endl;
    std::cout << "  Max |B_y|:            " << final_diag.max_By << std::endl;
    std::cout << "  Max |∇·B|:            " << final_diag.max_div_B << std::endl;

    return 0;
}
