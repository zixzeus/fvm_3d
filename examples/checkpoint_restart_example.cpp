#include "core/fvm_solver3d.hpp"
#include "io/hdf5_checkpoint.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace fvm3d;

/**
 * Simple initialization: Gaussian pulse
 */
void gaussian_pulse_init(double x, double y, double z, Eigen::VectorXd& U) {
    const double gamma = 1.4;
    const double x0 = 0.5, y0 = 0.5, z0 = 0.5;
    const double sigma = 0.1;

    // Base state
    double rho_base = 1.0;
    double p_base = 1.0;

    // Add Gaussian perturbation to pressure
    double r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0) + (z - z0) * (z - z0);
    double p_pert = 0.1 * std::exp(-r2 / (sigma * sigma));

    double rho = rho_base;
    double p = p_base + p_pert;
    double u = 0.0;
    double v = 0.0;
    double w = 0.0;

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
    std::cout << "\n=== FVM3D Checkpoint/Restart Example ===" << std::endl;

    // Configure solver
    core::FVMSolverConfig config;
    config.xmin = 0.0;
    config.ymin = 0.0;
    config.zmin = 0.0;
    config.Lx = 1.0;
    config.Ly = 1.0;
    config.Lz = 1.0;
    config.nx = 32;
    config.ny = 32;
    config.nz = 32;
    config.nghost = 2;

    config.riemann_solver = "hllc";
    config.reconstruction = "muscl";
    config.reconstruction_limiter = "van_leer";
    config.time_integrator = "rk2";
    config.boundary_condition = "transmissive";
    config.bc_x = true;
    config.bc_y = true;
    config.bc_z = true;

    config.cfl = 0.3;
    config.t_final = 0.1;
    config.num_steps = 1000;
    config.output_interval = 50;
    config.verbose = 1;

    std::string checkpoint_file = "simulation_checkpoint.h5";

    try {
        // Phase 1: Run simulation and save checkpoint
        std::cout << "\n=== Phase 1: Run Simulation and Save Checkpoint ===" << std::endl;

        core::FVMSolver3D solver(config);
        solver.initialize(gaussian_pulse_init);

        std::cout << "\nRunning for 500 steps (t=0 to ~t=0.05)...\n";
        for (int step = 0; step < 500; step++) {
            solver.step();
            if (solver.step_count() % config.output_interval == 0) {
                std::cout << "Step " << std::setw(5) << solver.step_count()
                         << " | Time: " << std::scientific << std::setprecision(6) << solver.time() << "\n";
            }
        }

        double time_checkpoint = solver.time();
        int step_checkpoint = solver.step_count();

        std::cout << "\nCheckpoint at:\n"
                 << "  Time: " << std::scientific << std::setprecision(6) << time_checkpoint << "\n"
                 << "  Step: " << step_checkpoint << "\n";

        // Save checkpoint
        std::cout << "\nSaving checkpoint to " << checkpoint_file << "...\n";
        solver.save_checkpoint(checkpoint_file, "Gaussian pulse at mid-simulation");

        // Get state value at center
        const auto& state_phase1 = solver.state();
        int i_center = config.nx / 2 + config.nghost;
        int j_center = config.ny / 2 + config.nghost;
        int k_center = config.nz / 2 + config.nghost;
        double rho_checkpoint = state_phase1(0, i_center, j_center, k_center);

        std::cout << "State value at center (rho): " << std::scientific << std::setprecision(6)
                 << rho_checkpoint << "\n";

        // Phase 2: Load checkpoint and continue
        std::cout << "\n=== Phase 2: Load Checkpoint and Continue ===" << std::endl;

        core::FVMSolver3D solver2(config);
        solver2.initialize(gaussian_pulse_init);

        // Load checkpoint (this overwrites current state with saved state)
        std::cout << "Loading checkpoint from " << checkpoint_file << "...\n";
        if (solver2.load_checkpoint(checkpoint_file)) {
            std::cout << "Checkpoint loaded successfully!\n"
                     << "  Time: " << std::scientific << std::setprecision(6) << solver2.time() << "\n"
                     << "  Step: " << solver2.step_count() << "\n";

            // Verify state matches
            const auto& state_phase2 = solver2.state();
            double rho_loaded = state_phase2(0, i_center, j_center, k_center);
            std::cout << "State value at center (rho): " << std::scientific << std::setprecision(6)
                     << rho_loaded << "\n";

            double state_diff = std::abs(rho_loaded - rho_checkpoint);
            std::cout << "State difference: " << state_diff << "\n";

            if (state_diff < 1e-10) {
                std::cout << "✓ State matches saved checkpoint!\n";
            } else {
                std::cout << "⚠ State difference detected (may be normal due to FP precision)\n";
            }

            // Continue simulation for another 500 steps
            std::cout << "\nContinuing simulation for another 500 steps (t~0.05 to t~0.1)...\n";
            int target_steps = step_checkpoint + 500;
            while (solver2.step_count() < target_steps && solver2.time() < config.t_final) {
                solver2.step();
                if (solver2.step_count() % config.output_interval == 0) {
                    std::cout << "Step " << std::setw(5) << solver2.step_count()
                             << " | Time: " << std::scientific << std::setprecision(6) << solver2.time() << "\n";
                }
            }

            std::cout << "\n✓ Restart completed successfully!\n"
                     << "  Final Time: " << std::scientific << std::setprecision(6) << solver2.time() << "\n"
                     << "  Final Step: " << solver2.step_count() << "\n";
        } else {
            std::cerr << "Failed to load checkpoint\n";
            return 1;
        }

        // Phase 3: Demonstrate metadata reading
        std::cout << "\n=== Phase 3: Read Checkpoint Metadata ===" << std::endl;

        double metadata_time;
        int metadata_step;
        std::string metadata_desc;

        if (io::HDF5Checkpoint::read_metadata(checkpoint_file, metadata_time, metadata_step, metadata_desc)) {
            std::cout << "Checkpoint metadata:\n"
                     << "  Time: " << std::scientific << std::setprecision(6) << metadata_time << "\n"
                     << "  Step: " << metadata_step << "\n"
                     << "  Description: " << metadata_desc << "\n";
            std::cout << "✓ Metadata read successfully!\n";
        }

        std::cout << "\n=== Example Completed Successfully ===" << std::endl;
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
