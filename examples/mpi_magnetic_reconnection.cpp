/**
 * MPI-Parallel 3D Magnetic Reconnection Simulation
 *
 * Demonstrates:
 * - MPIFVMSolver3D with MHD physics
 * - Advanced resistive MHD equations
 * - Harris sheet initial condition
 * - Parallel execution and statistics
 * - HLLD Riemann solver for MHD
 *
 * Run with: mpirun -np 8 ./mpi_magnetic_reconnection
 */

#include "core/mpi_fvm_solver3d.hpp"
#include "physics/resistive_mhd3d_advanced.hpp"
#include "parallel/mpi_utils.hpp"
#include "io/parallel_vtk_writer.hpp"
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace fvm3d;

int main(int argc, char** argv) {
    // Initialize MPI
    parallel::MPIGuard mpi_guard(argc, argv);

    int rank = parallel::MPIUtils::rank();
    int size = parallel::MPIUtils::size();

    if (rank == 0) {
        std::cout << "=== MPI-Parallel 3D Magnetic Reconnection ===\n";
        std::cout << "Running on " << size << " MPI processes\n\n";
    }

    // Configuration
    core::MPIFVMSolverConfig config;

    // Global grid (Harris sheet domain)
    config.xmin = -2.0;
    config.ymin = -1.0;
    config.zmin = -1.0;
    config.Lx = 4.0;
    config.Ly = 2.0;
    config.Lz = 2.0;
    config.nx = 64;  // Moderate resolution for faster demo
    config.ny = 32;
    config.nz = 32;
    config.nghost = 2;

    // MPI decomposition (auto-determine)
    config.px = 0;  // Will be auto-determined
    config.py = 0;
    config.pz = 0;

    // Physics (Advanced resistive MHD with GLM)
    config.physics_type = "mhd_advanced";
    config.num_vars = 9;  // [rho, rho_u, rho_v, rho_w, E, Bx, By, Bz, psi]

    // Numerical schemes
    config.riemann_solver = "hlld";  // HLLD solver for MHD
    config.reconstruction = "constant";  // Use constant for stability
    config.reconstruction_limiter = "minmod";
    config.time_integrator = "rk2";
    config.boundary_condition = "transmissive";

    // Boundary conditions (periodic in X, transmissive in Y, Z)
    config.bc_x = true;   // Periodic in reconnection direction
    config.bc_y = false;  // Transmissive in current sheet direction
    config.bc_z = false;  // Transmissive in out-of-plane direction

    // Time stepping (more restrictive for MHD)
    config.cfl = 0.3;
    config.t_final = 2.0;  // Short run for quick visualization
    config.num_steps = 500;  // Limit steps for quick demo
    config.output_interval = 50;
    config.checkpoint_interval = 0;  // Disable checkpoints for demo

    // Verbosity
    config.verbose = 1;

    try {
        // Create MPI-parallel solver
        core::MPIFVMSolver3D solver(config);

        if (rank == 0) {
            std::cout << "\nInitializing Harris sheet with reconnection trigger...\n";
        }

        // Harris sheet initial condition
        // This is a force-balanced magnetic configuration with current sheet
        physics::AdvancedResistiveMHD3D::HarrisSheetConfig harris_config;
        harris_config.B0 = 1.0;           // Asymptotic field strength
        harris_config.n0 = 1.0;           // Background density
        harris_config.p0 = 0.1;           // Reference pressure
        harris_config.L_sheet = 0.5;      // Current sheet thickness
        harris_config.beta = 1.0;         // Plasma beta

        // Perturbation to trigger reconnection
        harris_config.perturbation_amplitude = 0.03;  // 3% perturbation

        // Resistivity model (position-dependent)
        physics::AdvancedResistiveMHD3D::ResistivityModel resistivity;
        resistivity.eta0 = 1e-3;               // Background resistivity (Rm ~ 1000)
        resistivity.eta1 = 0.01667;            // Enhanced resistivity (Rm ~ 60)
        resistivity.localization_scale = 1.0;  // Width of enhanced region

        // GLM parameters for divergence cleaning
        physics::AdvancedResistiveMHD3D::GLMParameters glm;
        glm.ch = 0.2;  // Wave speed for divergence propagation
        glm.cr = 0.2;  // Damping rate

        // Create MHD physics object
        physics::AdvancedResistiveMHD3D mhd_physics(resistivity, glm);

        // Initial condition function
        auto harris_init = [&](double x, double y, double z, Eigen::VectorXd& U) {
            U = mhd_physics.harris_sheet_initial(x, y, z, harris_config);
        };

        // Initialize
        solver.initialize(harris_init);

        if (rank == 0) {
            std::cout << "\nHarris sheet configuration:\n";
            std::cout << "  B0 = " << harris_config.B0 << "\n";
            std::cout << "  Current sheet thickness L = " << harris_config.L_sheet << "\n";
            std::cout << "  Plasma beta = " << harris_config.beta << "\n";
            std::cout << "  Background Rm = " << 1.0 / resistivity.eta0 << "\n";
            std::cout << "  Enhanced Rm = " << 1.0 / resistivity.eta1 << "\n";
            std::cout << "  Perturbation amplitude = " << harris_config.perturbation_amplitude << "\n\n";
            std::cout << "Running simulation with VTK output...\n\n";
        }

        // Export initial condition to VTK
        io::export_parallel_state_to_vtk(
            "reconnection_0000",
            solver.state(),
            solver.grid(),
            solver.decomposer(),
            9,  // num_vars (MHD with GLM)
            "mhd",
            0.0,
            MPI_COMM_WORLD
        );

        if (rank == 0) {
            std::cout << "Initial condition exported to reconnection_0000.pvti\n\n";
        }

        // Run simulation with periodic VTK output
        int vtk_output_interval = 50;  // Export every 50 steps for quick demo
        int vtk_count = 1;

        while (solver.step_count() < config.num_steps && solver.time() < config.t_final) {
            solver.step();

            // Export to VTK periodically
            if (solver.step_count() % vtk_output_interval == 0) {
                char base_filename[256];
                snprintf(base_filename, sizeof(base_filename), "reconnection_%04d", vtk_count);

                io::export_parallel_state_to_vtk(
                    base_filename,
                    solver.state(),
                    solver.grid(),
                    solver.decomposer(),
                    9,
                    "mhd",
                    solver.time(),
                    MPI_COMM_WORLD
                );

                if (rank == 0) {
                    std::cout << "  -> Exported VTK: " << base_filename << ".pvti (t=" << solver.time() << ")\n";
                }
                vtk_count++;
            }
        }

        // Export final state
        char final_base_filename[256];
        snprintf(final_base_filename, sizeof(final_base_filename), "reconnection_%04d", vtk_count);
        io::export_parallel_state_to_vtk(
            final_base_filename,
            solver.state(),
            solver.grid(),
            solver.decomposer(),
            9,
            "mhd",
            solver.time(),
            MPI_COMM_WORLD
        );

        // Print final statistics
        if (rank == 0) {
            std::cout << "\n=== Simulation Complete ===\n";
            std::cout << "Final time: " << solver.time() << "\n";
            std::cout << "Total steps: " << solver.step_count() << "\n";
            std::cout << "VTK files: reconnection_*.pvti (" << (vtk_count + 1) << " timesteps)\n\n";

            std::cout << "Expected physics:\n";
            std::cout << "  - Magnetic reconnection at X-point (x=0, y=0)\n";
            std::cout << "  - Formation of magnetic islands\n";
            std::cout << "  - Out-of-plane current concentration\n";
            std::cout << "  - Plasma heating and acceleration\n";
            std::cout << "  - Energy conversion: magnetic -> kinetic + thermal\n\n";

            std::cout << "To visualize:\n";
            std::cout << "  1. Open ParaView\n";
            std::cout << "  2. Load reconnection_*.pvti files\n";
            std::cout << "  3. Visualize magnetic field lines, density, pressure\n";
            std::cout << "  4. Look for X-point and magnetic island structures\n";
        }

        // TODO: Add energy diagnostics
        // - Magnetic energy
        // - Kinetic energy
        // - Thermal energy
        // - Reconnection rate
        // - Maximum out-of-plane field (By_max indicator)

    } catch (const std::exception& e) {
        if (rank == 0) {
            std::cerr << "Error: " << e.what() << "\n";
        }
        return 1;
    }

    return 0;
}
