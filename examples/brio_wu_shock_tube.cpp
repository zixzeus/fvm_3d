#include "core/fvm_solver3d.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace fvm3d;

/**
 * Brio-Wu Magnetic Shock Tube Test Case
 *
 * Classic 1D MHD Riemann problem for validating HLLD Riemann solver.
 * Features magnetic field flip across discontinuity.
 *
 * Initial conditions (Brio & Wu 1988):
 * Left state:  rho=1.0, u=0, v=0, w=0, p=1.0, Bx=0.75, By=1.0, Bz=0
 * Right state: rho=0.125, u=0, v=0, w=0, p=0.1, Bx=0.75, By=-1.0, Bz=0
 *
 * This test case:
 * - Has 5 MHD waves: fast MS, slow MS, Alfvén, contact, slow MS, fast MS
 * - Magnetic field has opposite sign (By flips)
 * - Good test for Alfvén wave handling in HLLD solver
 * - Published exact solution available for validation
 */
void brio_wu_initial_condition(double x, double y, double z, Eigen::VectorXd& U) {
    const double gamma = 5.0/3.0;  // Adiabatic index for MHD
    double rho, u, v, w, p, Bx, By, Bz;

    // Use x-coordinate for shock tube test
    if (x < 0.5) {
        // Left state
        rho = 1.0;
        u = 0.0;
        v = 0.0;
        w = 0.0;
        p = 1.0;
        Bx = 0.75;
        By = 1.0;
        Bz = 0.0;
    } else {
        // Right state
        rho = 0.125;
        u = 0.0;
        v = 0.0;
        w = 0.0;
        p = 0.1;
        Bx = 0.75;
        By = -1.0;
        Bz = 0.0;
    }

    // Convert to conservative variables (9-var GLM-MHD system)
    double rho_u = rho * u;
    double rho_v = rho * v;
    double rho_w = rho * w;

    // Magnetic energy
    double B2 = Bx * Bx + By * By + Bz * Bz;

    // Total energy
    double E = p / (gamma - 1.0) + 0.5 * rho * (u*u + v*v + w*w) + 0.5 * B2;

    // GLM divergence cleaning variable
    double psi = 0.0;

    U.resize(9);
    U(0) = rho;
    U(1) = rho_u;
    U(2) = rho_v;
    U(3) = rho_w;
    U(4) = E;
    U(5) = Bx;
    U(6) = By;
    U(7) = Bz;
    U(8) = psi;
}

int main() {
    std::cout << "\n╔════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  FVM3D Brio-Wu MHD Shock Tube Benchmark Test                ║" << std::endl;
    std::cout << "║  HLLD Riemann Solver Validation                             ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════╝\n" << std::endl;

    // Configure solver
    core::FVMSolverConfig config;
    config.xmin = 0.0;
    config.ymin = 0.0;
    config.zmin = 0.0;
    config.Lx = 1.0;
    config.Ly = 0.01;   // Very thin domain in Y
    config.Lz = 0.01;   // Very thin domain in Z
    config.nx = 128;    // 1D shock tube with refined grid
    config.ny = 2;      // Minimal cells in Y
    config.nz = 2;      // Minimal cells in Z
    config.nghost = 2;

    // Physics configuration - Advanced MHD with GLM divergence cleaning
    config.physics_type = "mhd_advanced";
    config.num_vars = 9;  // GLM-MHD: [rho, rho_u, rho_v, rho_w, E, Bx, By, Bz, psi]

    // HLLD Riemann solver for MHD
    config.flux_calculator = "hlld";
    config.reconstruction = "muscl";
    config.reconstruction_limiter = "minmod";
    config.time_integrator = "rk2";
    config.boundary_condition = "transmissive";
    config.bc_x = true;  // Transmissive in X
    config.bc_y = true;  // Transmissive in Y
    config.bc_z = true;  // Transmissive in Z

    config.cfl = 0.4;
    config.t_final = 0.2;      // Run to t=0.2 for solution structure
    config.num_steps = 5000;
    config.output_interval = 50;
    config.verbose = 1;

    try {
        // Create and initialize solver
        core::FVMSolver3D solver(config);
        solver.initialize(brio_wu_initial_condition);

        // Print initial state info
        std::cout << "\n📋 Initial Conditions:" << std::endl;
        std::cout << "  Left state (x < 0.5):" << std::endl;
        std::cout << "    ρ=1.0, u=0, p=1.0, Bx=0.75, By=1.0, Bz=0" << std::endl;
        std::cout << "  Right state (x > 0.5):" << std::endl;
        std::cout << "    ρ=0.125, u=0, p=0.1, Bx=0.75, By=-1.0, Bz=0" << std::endl;
        std::cout << "\n  Key feature: Magnetic field flip (By changes sign)" << std::endl;
        std::cout << "  Expected waves: Fast MS, Slow MS, Alfvén, Contact, Slow MS, Fast MS\n" << std::endl;

        // Print solver configuration
        std::cout << "⚙️  Solver Configuration:" << std::endl;
        std::cout << "  Riemann Solver: HLLD (5-wave for MHD)" << std::endl;
        std::cout << "  Reconstruction: MUSCL with minmod limiter" << std::endl;
        std::cout << "  Time Integration: RK2" << std::endl;
        std::cout << "  CFL: " << config.cfl << std::endl;
        std::cout << "  Grid: " << config.nx << " × " << config.ny << " × " << config.nz << std::endl;
        std::cout << "  Final time: " << config.t_final << "\n" << std::endl;

        // Run simulation
        std::cout << "▶️  Running simulation..." << std::endl;
        solver.run();

        // Print final statistics
        std::cout << "\n✅ Simulation completed successfully!" << std::endl;
        std::cout << "  Output files written to: ./output/" << std::endl;
        std::cout << "\n📊 Expected Results (from Brio & Wu 1988):" << std::endl;
        std::cout << "  - Density range: [0.125, 1.0]" << std::endl;
        std::cout << "  - Pressure range: [0.1, 1.0]" << std::endl;
        std::cout << "  - Clear Alfvén wave discontinuities in By" << std::endl;
        std::cout << "  - Smooth density and pressure profiles between waves" << std::endl;
        std::cout << "\n🔍 Analysis:" << std::endl;
        std::cout << "  1. Check density profile for 5-wave structure" << std::endl;
        std::cout << "  2. Verify Alfvén waves by examining By discontinuities" << std::endl;
        std::cout << "  3. Validate against published exact solution" << std::endl;
        std::cout << "  4. Check divergence cleaning (psi should remain ~0)" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\n❌ Error during simulation:" << std::endl;
        std::cerr << "  " << e.what() << std::endl;
        return 1;
    }
}
