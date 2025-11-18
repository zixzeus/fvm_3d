#include "core/fvm_solver3d.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace fvm3d;

/**
 * Orszag-Tang Vortex MHD Test Case
 *
 * Smooth 2D MHD problem with periodic boundary conditions.
 * Tests the solver's ability to maintain smooth structures and
 * handle complex magnetic field evolution.
 *
 * Initial conditions (Orszag & Tang 1979):
 * Domain: [0, 2π] × [0, 2π] (periodic in both directions)
 *
 * Velocity: u = -sin(y), v = sin(x), w = 0
 * Magnetic: Bx = -sin(y), By = sin(2x), Bz = 0
 * Pressure: p = (γ-1)/(γ*π²) * (Bx² + By²)/2 + ρ*γ*p/(γ-1)
 *
 * where:
 *   ρ = γ² (to get p₀ = 1 when γ=5/3)
 *   p₀ = 1 (ambient pressure)
 *
 * Features:
 * - Smooth initial conditions with vortex structures
 * - Magnetic field creates shear and reconnection
 * - Good test for stability and accuracy with smooth flows
 * - Can develop small-scale structures from smooth initial conditions
 * - Published numerical solutions available for validation
 */
void orszag_tang_initial_condition(double x, double y, double z, Eigen::VectorXd& U) {
    const double gamma = 5.0/3.0;
    const double M_PI_local = 3.141592653589793;

    // Normalize coordinates to [0, 2π]
    double x_norm = x;
    double y_norm = y;

    // Density (constant)
    double rho = gamma * gamma;

    // Velocity components
    double u = -std::sin(y_norm);
    double v = std::sin(x_norm);
    double w = 0.0;

    // Magnetic field components
    double Bx = -std::sin(y_norm);
    double By = std::sin(2.0 * x_norm);
    double Bz = 0.0;

    // Pressure
    // p = (γ-1)/(γ*π²) * (Bx² + By²)/2 + p₀
    // where p₀ = 1
    double B2 = Bx * Bx + By * By;
    double p = (gamma - 1.0) / (gamma * M_PI_local * M_PI_local) * (B2 / 2.0) + 1.0;

    // Convert to conservative variables
    double rho_u = rho * u;
    double rho_v = rho * v;
    double rho_w = rho * w;

    // Total energy
    double E = p / (gamma - 1.0) + 0.5 * rho * (u*u + v*v + w*w) + 0.5 * B2;

    // GLM divergence cleaning
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
    std::cout << "║  FVM3D Orszag-Tang Vortex MHD Benchmark Test                ║" << std::endl;
    std::cout << "║  Smooth 2D MHD Problem with Complex Structures              ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════╝\n" << std::endl;

    const double M_PI_local = 3.141592653589793;

    // Configure solver
    core::FVMSolverConfig config;
    config.xmin = 0.0;
    config.ymin = 0.0;
    config.zmin = 0.0;
    config.Lx = 2.0 * M_PI_local;
    config.Ly = 2.0 * M_PI_local;
    config.Lz = 0.1;        // Thin domain in Z for 2D problem
    config.nx = 96;         // 2D resolution: 96×96
    config.ny = 96;
    config.nz = 2;          // Minimal cells in Z
    config.nghost = 2;

    // Physics configuration - Advanced MHD with GLM divergence cleaning
    config.physics_type = "mhd_advanced";
    config.num_vars = 9;    // GLM-MHD: [rho, rho_u, rho_v, rho_w, E, Bx, By, Bz, psi]

    // HLLD Riemann solver for MHD
    config.flux_calculator = "hlld";
    config.reconstruction = "muscl";
    config.reconstruction_limiter = "minmod";
    config.time_integrator = "rk2";
    config.boundary_condition = "periodic";   // Periodic BC for vortex
    config.bc_x = false;    // Periodic in X
    config.bc_y = false;    // Periodic in Y
    config.bc_z = true;     // Transmissive in Z

    config.cfl = 0.3;       // Conservative for 2D
    config.t_final = 5.0;   // Run to t=5.0 to see development
    config.num_steps = 10000;
    config.output_interval = 100;
    config.verbose = 1;

    try {
        // Create and initialize solver
        core::FVMSolver3D solver(config);
        solver.initialize(orszag_tang_initial_condition);

        // Print initial state info
        std::cout << "\n📋 Initial Conditions (Orszag & Tang 1979):" << std::endl;
        std::cout << "  Domain: [0, 2π] × [0, 2π] with periodic boundaries" << std::endl;
        std::cout << "\n  Velocity field:" << std::endl;
        std::cout << "    u(x,y) = -sin(y)" << std::endl;
        std::cout << "    v(x,y) = sin(x)" << std::endl;
        std::cout << "    w = 0" << std::endl;
        std::cout << "\n  Magnetic field:" << std::endl;
        std::cout << "    Bx(x,y) = -sin(y)" << std::endl;
        std::cout << "    By(x,y) = sin(2x)" << std::endl;
        std::cout << "    Bz = 0" << std::endl;
        std::cout << "\n  Density: ρ = γ² = " << (5.0/3.0 * 5.0/3.0) << std::endl;
        std::cout << "  Pressure: p = 1.0 + magnetic pressure contribution" << std::endl;
        std::cout << "\n  Key features:" << std::endl;
        std::cout << "  - Smooth initial conditions" << std::endl;
        std::cout << "  - Complex magnetic field with shear layers" << std::endl;
        std::cout << "  - Can develop turbulence-like structures" << std::endl;
        std::cout << "  - Tests solver accuracy on smooth problems\n" << std::endl;

        // Print solver configuration
        std::cout << "⚙️  Solver Configuration:" << std::endl;
        std::cout << "  Riemann Solver: HLLD (5-wave for MHD)" << std::endl;
        std::cout << "  Reconstruction: MUSCL with minmod limiter" << std::endl;
        std::cout << "  Time Integration: RK2" << std::endl;
        std::cout << "  Boundary Conditions: Periodic in X-Y, Transmissive in Z" << std::endl;
        std::cout << "  CFL: " << config.cfl << std::endl;
        std::cout << "  Grid: " << config.nx << " × " << config.ny << " × " << config.nz << std::endl;
        std::cout << "  Final time: " << config.t_final << "\n" << std::endl;

        // Run simulation
        std::cout << "▶️  Running simulation..." << std::endl;
        solver.run();

        // Print final statistics
        std::cout << "\n✅ Simulation completed successfully!" << std::endl;
        std::cout << "  Output files written to: ./output/" << std::endl;
        std::cout << "\n📊 Expected Results:" << std::endl;
        std::cout << "  - Smooth initial structures remain stable" << std::endl;
        std::cout << "  - Small-scale structures develop from smooth initial conditions" << std::endl;
        std::cout << "  - Magnetic field lines remain well-represented" << std::endl;
        std::cout << "  - Density remains positive throughout" << std::endl;
        std::cout << "  - Energy is approximately conserved" << std::endl;
        std::cout << "\n🔍 Analysis:" << std::endl;
        std::cout << "  1. Check divergence cleaning (∇·B, psi evolution)" << std::endl;
        std::cout << "  2. Examine density and pressure positivity preservation" << std::endl;
        std::cout << "  3. Verify energy evolution and dissipation" << std::endl;
        std::cout << "  4. Compare with published benchmark results" << std::endl;
        std::cout << "  5. Look for expected magnetic reconnection patterns" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\n❌ Error during simulation:" << std::endl;
        std::cerr << "  " << e.what() << std::endl;
        return 1;
    }
}
