/**
 * VTK Visualization Example
 *
 * Demonstrates:
 * - Running a 3D Euler solver
 * - Exporting snapshots to VTK format
 * - Visualizing results with ParaView
 *
 * Usage:
 *   ./vtk_visualization_example
 *   # Then open output_*.vti files in ParaView
 */

#include "core/fvm_solver3d.hpp"
#include "io/vtk_writer.hpp"
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace fvm3d;

int main() {
    std::cout << "=== VTK Visualization Example ===\n\n";

    // Configuration
    core::FVMSolverConfig config;

    // Grid parameters
    config.xmin = -1.0;
    config.ymin = -1.0;
    config.zmin = -1.0;
    config.Lx = 2.0;
    config.Ly = 2.0;
    config.Lz = 2.0;
    config.nx = 32;
    config.ny = 32;
    config.nz = 32;
    config.nghost = 2;

    // Numerical schemes
    config.flux_calculator = "hllc";
    config.reconstruction = "constant";
    config.reconstruction_limiter = "minmod";
    config.time_integrator = "rk2";
    config.boundary_condition = "transmissive";

    // Boundary conditions
    config.bc_x = true;
    config.bc_y = true;
    config.bc_z = true;

    // Time stepping
    config.cfl = 0.4;
    config.t_final = 0.2;
    config.num_steps = 1000;
    config.output_interval = 20;  // Print every 20 steps

    config.verbose = 1;

    // Create solver
    core::FVMSolver3D solver(config);

    // 3D Blast wave initial condition
    auto blast_wave_init = [](double x, double y, double z, Eigen::VectorXd& U) {
        double r = std::sqrt(x*x + y*y + z*z);
        double gamma = 1.4;

        double rho, p;
        if (r < 0.3) {
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

        // Convert to conservative variables
        U.resize(5);
        U(0) = rho;
        U(1) = rho * u;
        U(2) = rho * v;
        U(3) = rho * w;
        U(4) = p / (gamma - 1.0) + 0.5 * rho * (u*u + v*v + w*w);
    };

    // Initialize
    std::cout << "Initializing 3D blast wave...\n";
    solver.initialize(blast_wave_init);

    // Export initial condition to VTK
    std::cout << "\nExporting initial condition to VTK...\n";
    io::export_state_to_vtk(
        "output_0000.vti",
        solver.state(),
        solver.grid(),
        5,  // num_vars
        "euler",
        0.0,  // time
        true  // binary
    );

    // Run simulation with periodic VTK output
    std::cout << "\nRunning simulation with VTK output...\n\n";

    int vtk_output_interval = 50;  // Export every 50 steps
    int vtk_count = 1;

    while (solver.step_count() < config.num_steps && solver.time() < config.t_final) {
        solver.step();

        // Print progress
        if (config.output_interval > 0 && solver.step_count() % config.output_interval == 0) {
            std::cout << "Step " << solver.step_count()
                      << " | t = " << solver.time()
                      << "\n";
        }

        // Export to VTK periodically
        if (solver.step_count() % vtk_output_interval == 0) {
            char filename[256];
            snprintf(filename, sizeof(filename), "output_%04d.vti", vtk_count);

            io::export_state_to_vtk(
                filename,
                solver.state(),
                solver.grid(),
                5,
                "euler",
                solver.time(),
                true
            );

            std::cout << "  -> Exported VTK: " << filename << "\n";
            vtk_count++;
        }
    }

    // Export final state
    std::cout << "\nExporting final state to VTK...\n";
    char final_filename[256];
    snprintf(final_filename, sizeof(final_filename), "output_%04d.vti", vtk_count);
    io::export_state_to_vtk(
        final_filename,
        solver.state(),
        solver.grid(),
        5,
        "euler",
        solver.time(),
        true
    );

    std::cout << "\n=== Simulation Complete ===\n";
    std::cout << "  Final time: " << solver.time() << "\n";
    std::cout << "  Total steps: " << solver.step_count() << "\n";
    std::cout << "  VTK files: output_*.vti (" << (vtk_count + 1) << " files)\n\n";

    std::cout << "To visualize:\n";
    std::cout << "  1. Open ParaView\n";
    std::cout << "  2. File > Open > output_..vti\n";
    std::cout << "  3. Click 'Apply' to load data\n";
    std::cout << "  4. Select field to visualize (density, pressure, velocity, etc.)\n";
    std::cout << "  5. Use animation controls to step through time\n\n";

    return 0;
}
