#include "core/fvm_solver3d.hpp"
#include "io/hdf5_checkpoint.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace fvm3d;

void simple_init(double x, double y, double z, Eigen::VectorXd& U) {
    const double gamma = 1.4;
    U.resize(5);
    U(0) = 1.0;           // rho
    U(1) = 0.0;           // rho_u
    U(2) = 0.0;           // rho_v
    U(3) = 0.0;           // rho_w
    U(4) = 1.0 / (gamma - 1.0);  // E
}

int main() {
    std::cout << "=== HDF5 Checkpoint Test ===" << std::endl;

    try {
        // Create minimal solver
        core::FVMSolverConfig config;
        config.xmin = 0.0;
        config.ymin = 0.0;
        config.zmin = 0.0;
        config.Lx = 1.0;
        config.Ly = 1.0;
        config.Lz = 1.0;
        config.nx = 8;     // Small grid for fast test
        config.ny = 8;
        config.nz = 8;
        config.nghost = 2;

        config.flux_calculator = "hll";
        config.reconstruction = "muscl";
        config.reconstruction_limiter = "van_leer";
        config.time_integrator = "rk2";
        config.boundary_condition = "periodic";
        config.bc_x = true;
        config.bc_y = true;
        config.bc_z = true;

        config.cfl = 0.4;
        config.t_final = 1.0;
        config.num_steps = 100;
        config.output_interval = 10;
        config.verbose = 1;

        std::cout << "\n1. Creating solver and initializing...\n";
        core::FVMSolver3D solver(config);
        solver.initialize(simple_init);

        std::cout << "2. Running 10 steps...\n";
        for (int i = 0; i < 10; i++) {
            solver.step();
        }

        std::cout << "   Current time: " << std::scientific << std::setprecision(6) << solver.time() << "\n"
                 << "   Current step: " << solver.step_count() << "\n";

        std::string checkpoint_file = "test_checkpoint.h5";

        std::cout << "\n3. Saving checkpoint to " << checkpoint_file << "...\n";
        solver.save_checkpoint(checkpoint_file, "Test checkpoint at step 10");

        std::cout << "4. Checking file size...\n";
        std::system(("ls -lh " + checkpoint_file).c_str());

        std::cout << "\n5. Reading checkpoint metadata...\n";
        double time_read;
        int step_read;
        std::string desc_read;
        if (io::HDF5Checkpoint::read_metadata(checkpoint_file, time_read, step_read, desc_read)) {
            std::cout << "   Metadata read successfully!\n"
                     << "   Time: " << std::scientific << std::setprecision(6) << time_read << "\n"
                     << "   Step: " << step_read << "\n"
                     << "   Description: " << desc_read << "\n";
        } else {
            std::cout << "   Failed to read metadata\n";
        }

        std::cout << "\n6. Creating new solver and loading checkpoint...\n";
        core::FVMSolver3D solver2(config);
        solver2.initialize(simple_init);

        if (solver2.load_checkpoint(checkpoint_file)) {
            std::cout << "   Checkpoint loaded successfully!\n"
                     << "   Time: " << std::scientific << std::setprecision(6) << solver2.time() << "\n"
                     << "   Step: " << solver2.step_count() << "\n";

            std::cout << "\n✓ HDF5 Checkpoint Test Passed!\n";
        } else {
            std::cout << "   Failed to load checkpoint\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
