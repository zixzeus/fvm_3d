#include "core/fvm_solver3d.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace fvm3d;

/**
 * 3D Sod Shock Tube / Blast Wave Test Case
 *
 * Classic shock tube test for validating FVM solvers.
 * Left state: high pressure (P=1.0, rho=1.0)
 * Right state: low pressure (P=0.1, rho=0.125)
 * Center position: x = 0.5
 */
void blast_wave_3d(double x, double y, double z, Eigen::VectorXd& U) {
    const double gamma = 1.4;
    double rho, u, v, w, p;

    // Use x-coordinate for shock tube test
    if (x < 0.5) {
        // Left state: high pressure
        rho = 1.0;
        p = 1.0;
    } else {
        // Right state: low pressure
        rho = 0.125;
        p = 0.1;
    }

    u = 0.0;
    v = 0.0;
    w = 0.0;

    // Convert to conservative variables
    double rho_u = rho * u;
    double rho_v = rho * v;
    double rho_w = rho * w;
    double E = p / (gamma - 1.0) + 0.5 * rho * (u*u + v*v + w*w);

    U.resize(5);
    U(0) = rho;
    U(1) = rho_u;
    U(2) = rho_v;
    U(3) = rho_w;
    U(4) = E;
}

int main() {
    std::cout << "\n=== FVM3D Blast Wave Simulation with Full Solver ===" << std::endl;

    // Configure solver
    core::FVMSolverConfig config;
    config.xmin = 0.0;
    config.ymin = 0.0;
    config.zmin = 0.0;
    config.Lx = 1.0;
    config.Ly = 0.2;  // Thin domain in Y
    config.Lz = 0.2;  // Thin domain in Z
    config.nx = 64;
    config.ny = 8;
    config.nz = 8;
    config.nghost = 2;

    config.riemann_solver = "hllc";
    config.reconstruction = "muscl";
    config.reconstruction_limiter = "van_leer";
    config.time_integrator = "rk2";
    config.boundary_condition = "transmissive";
    config.bc_x = true;  // Transmissive in X
    config.bc_y = true;  // Transmissive in Y
    config.bc_z = true;  // Transmissive in Z

    config.cfl = 0.4;
    config.t_final = 0.2;
    config.num_steps = 10000;
    config.output_interval = 100;
    config.verbose = 1;

    try {
        // Create and initialize solver
        core::FVMSolver3D solver(config);
        solver.initialize(blast_wave_3d);

        // Print initial state
        std::cout << "\nInitial state at domain center:\n";
        const auto& state = solver.state();
        int i_center = config.nx / 2 + config.nghost;
        int j_center = config.ny / 2 + config.nghost;
        int k_center = config.nz / 2 + config.nghost;
        std::cout << "  U = [" << state(0, i_center, j_center, k_center);
        for (int v = 1; v < 5; v++) {
            std::cout << ", " << state(v, i_center, j_center, k_center);
        }
        std::cout << "]\n\n";

        // Run simulation
        std::cout << "Starting simulation..." << std::endl;
        solver.run();

        // Print final state
        std::cout << "\nFinal state at domain center:\n";
        std::cout << "  U = [" << state(0, i_center, j_center, k_center);
        for (int v = 1; v < 5; v++) {
            std::cout << ", " << state(v, i_center, j_center, k_center);
        }
        std::cout << "]\n";

        std::cout << "\n=== Simulation Completed Successfully ===" << std::endl;
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
