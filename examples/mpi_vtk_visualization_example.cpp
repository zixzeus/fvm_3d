/**
 * MPI Parallel VTK Visualization Example
 *
 * Demonstrates:
 * - Running MPI-parallel 3D Euler solver
 * - Exporting parallel snapshots to PVTK format
 * - Visualizing parallel results with ParaView
 *
 * Usage:
 *   mpirun -np 8 ./mpi_vtk_visualization_example
 *   # Then open output_*.pvti files in ParaView
 */

#include "core/mpi_fvm_solver3d.hpp"
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
        std::cout << "=== MPI Parallel VTK Visualization Example ===\n";
        std::cout << "Running on " << size << " MPI processes\n\n";
    }

    // Configuration
    core::MPIFVMSolverConfig config;

    // Global grid parameters
    config.xmin = -1.0;
    config.ymin = -1.0;
    config.zmin = -1.0;
    config.Lx = 2.0;
    config.Ly = 2.0;
    config.Lz = 2.0;
    config.nx = 64;
    config.ny = 64;
    config.nz = 64;
    config.nghost = 2;

    // MPI decomposition (auto-determine)
    config.px = 0;
    config.py = 0;
    config.pz = 0;

    // Physics
    config.physics_type = "euler";
    config.num_vars = 5;

    // Numerical schemes
    config.riemann_solver = "hllc";
    config.reconstruction = "constant";
    config.reconstruction_limiter = "minmod";
    config.time_integrator = "rk2";
    config.boundary_condition = "transmissive";

    config.bc_x = true;
    config.bc_y = true;
    config.bc_z = true;

    // Time stepping
    config.cfl = 0.4;
    config.t_final = 0.25;
    config.num_steps = 2000;
    config.output_interval = 50;

    config.verbose = 1;

    // Create MPI solver
    core::MPIFVMSolver3D solver(config);

    // 3D Blast wave initial condition
    auto blast_wave_init = [](double x, double y, double z, Eigen::VectorXd& U) {
        double r = std::sqrt(x*x + y*y + z*z);
        double gamma = 1.4;

        double rho, p;
        if (r < 0.4) {
            // High pressure region (blast)
            rho = 1.0;
            p = 10.0;
        } else {
            // Low pressure ambient
            rho = 1.0;
            p = 0.1;
        }

        // Zero velocity
        double u = 0.0;
        double v = 0.0;
        double w = 0.0;

        // Conservative variables
        U.resize(5);
        U(0) = rho;
        U(1) = rho * u;
        U(2) = rho * v;
        U(3) = rho * w;
        U(4) = p / (gamma - 1.0) + 0.5 * rho * (u*u + v*v + w*w);
    };

    // Initialize
    if (rank == 0) {
        std::cout << "Initializing 3D blast wave on all ranks...\n";
    }
    solver.initialize(blast_wave_init);

    // Export initial condition
    if (rank == 0) {
        std::cout << "\nExporting initial condition to parallel VTK...\n";
    }

    io::export_parallel_state_to_vtk(
        "parallel_output_0000",
        solver.state(),
        solver.grid(),
        solver.decomposer(),
        5,
        "euler",
        0.0,
        MPI_COMM_WORLD
    );

    // Run simulation with periodic parallel VTK output
    if (rank == 0) {
        std::cout << "\nRunning parallel simulation with PVTK output...\n\n";
    }

    int vtk_output_interval = 100;  // Export every 100 steps
    int vtk_count = 1;

    while (solver.step_count() < config.num_steps && solver.time() < config.t_final) {
        solver.step();

        // Export to parallel VTK periodically
        if (solver.step_count() % vtk_output_interval == 0) {
            char base_filename[256];
            snprintf(base_filename, sizeof(base_filename), "parallel_output_%04d", vtk_count);

            io::export_parallel_state_to_vtk(
                base_filename,
                solver.state(),
                solver.grid(),
                solver.decomposer(),
                5,
                "euler",
                solver.time(),
                MPI_COMM_WORLD
            );

            if (rank == 0) {
                std::cout << "  -> Exported PVTK: " << base_filename << ".pvti\n";
            }
            vtk_count++;
        }
    }

    // Export final state
    if (rank == 0) {
        std::cout << "\nExporting final state to parallel VTK...\n";
    }

    char final_base_filename[256];
    snprintf(final_base_filename, sizeof(final_base_filename), "parallel_output_%04d", vtk_count);
    io::export_parallel_state_to_vtk(
        final_base_filename,
        solver.state(),
        solver.grid(),
        solver.decomposer(),
        5,
        "euler",
        solver.time(),
        MPI_COMM_WORLD
    );

    if (rank == 0) {
        std::cout << "\n=== Parallel Simulation Complete ===\n";
        std::cout << "  Final time: " << solver.time() << "\n";
        std::cout << "  Total steps: " << solver.step_count() << "\n";
        std::cout << "  PVTK master files: parallel_output_*.pvti (" << (vtk_count + 1) << " files)\n";
        std::cout << "  PVTK piece files: " << (size * (vtk_count + 1)) << " .vti files\n\n";

        std::cout << "To visualize parallel data:\n";
        std::cout << "  1. Open ParaView\n";
        std::cout << "  2. File > Open > parallel_output_*.pvti (master file)\n";
        std::cout << "  3. Click 'Apply' to load all subdomains\n";
        std::cout << "  4. ParaView will automatically assemble the parallel data\n";
        std::cout << "  5. Select field to visualize (density, pressure, velocity, etc.)\n\n";
    }

    return 0;
}
