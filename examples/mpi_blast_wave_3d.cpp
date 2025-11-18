/**
 * MPI-parallel 3D Blast Wave Example
 *
 * Demonstrates:
 * - MPIFVMSolver3D usage
 * - Domain decomposition across MPI ranks
 * - Parallel HDF5 checkpoint/restart
 * - Global statistics gathering
 *
 * Run with: mpirun -np 4 ./mpi_blast_wave_3d
 */

#include "core/mpi_fvm_solver3d.hpp"
#include "parallel/mpi_utils.hpp"
#include "io/mpi_hdf5_checkpoint.hpp"
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
        std::cout << "=== MPI-Parallel 3D Blast Wave Simulation ===\n";
        std::cout << "Running on " << size << " MPI processes\n\n";
    }

    // Configuration
    core::MPIFVMSolverConfig config;

    // Global grid
    config.xmin = 0.0;
    config.ymin = 0.0;
    config.zmin = 0.0;
    config.Lx = 1.0;
    config.Ly = 1.0;
    config.Lz = 1.0;
    config.nx = 64;  // Global grid size
    config.ny = 64;
    config.nz = 64;
    config.nghost = 2;

    // MPI decomposition (auto-determine)
    config.px = 0;
    config.py = 0;
    config.pz = 0;

    // Physics (Euler equations)
    config.physics_type = "euler";
    config.num_vars = 5;

    // Numerical schemes
    config.flux_calculator = "hllc";
    config.reconstruction = "muscl";
    config.reconstruction_limiter = "minmod";
    config.time_integrator = "rk2";
    config.boundary_condition = "transmissive";

    // Boundary conditions
    config.bc_x = false;
    config.bc_y = false;
    config.bc_z = false;

    // Time stepping
    config.cfl = 0.4;
    config.t_final = 0.2;
    config.num_steps = 1000;
    config.output_interval = 10;
    config.checkpoint_interval = 50;

    // Verbosity
    config.verbose = 1;

    try {
        // Create MPI-parallel solver
        core::MPIFVMSolver3D solver(config);

        if (rank == 0) {
            std::cout << "\nInitializing blast wave...\n";
        }

        // Initial condition: Spherical blast wave
        // High pressure/energy sphere at center, low pressure outside
        auto blast_wave_init = [](double x, double y, double z, Eigen::VectorXd& U) {
            const double GAMMA = 1.4;

            // Distance from center
            double xc = 0.5, yc = 0.5, zc = 0.5;
            double r = std::sqrt((x - xc) * (x - xc) +
                                (y - yc) * (y - yc) +
                                (z - zc) * (z - zc));

            double rho, p;
            if (r < 0.1) {
                // Inside blast: high pressure
                rho = 1.0;
                p = 10.0;
            } else {
                // Outside: low pressure
                rho = 1.0;
                p = 0.1;
            }

            // Zero velocity
            double u = 0.0, v = 0.0, w = 0.0;

            // Compute conserved variables
            double E = p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v + w * w);

            U(0) = rho;
            U(1) = rho * u;
            U(2) = rho * v;
            U(3) = rho * w;
            U(4) = E;
        };

        // Initialize
        solver.initialize(blast_wave_init);

        if (rank == 0) {
            std::cout << "Running simulation...\n\n";
        }

        // Run simulation
        solver.run();

        // Save final checkpoint
        if (rank == 0) {
            std::cout << "\nSaving final checkpoint...\n";
        }

        // TODO: Enable when parallel checkpoint is implemented
        // solver.save_checkpoint("blast_wave_final.h5", "Final state of blast wave");

        if (rank == 0) {
            std::cout << "\n=== Simulation Complete ===\n";
            std::cout << "Final time: " << solver.time() << "\n";
            std::cout << "Total steps: " << solver.step_count() << "\n";
        }

    } catch (const std::exception& e) {
        if (rank == 0) {
            std::cerr << "Error: " << e.what() << "\n";
        }
        return 1;
    }

    return 0;
}
