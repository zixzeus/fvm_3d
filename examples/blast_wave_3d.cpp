#include "core/grid3d.hpp"
#include "core/field3d.hpp"
#include "physics/euler3d.hpp"
#include "spatial/flux_calculation/flux_calculator_factory.hpp"
#include "temporal/forward_euler.hpp"
#include "boundary/periodic_bc.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace fvm3d;

/**
 * 3D Spherical Blast Wave Initial Condition.
 * High pressure sphere at center, ambient pressure elsewhere.
 */
void initialize_blast_wave(
    double x, double y, double z,
    double xmin, double ymin, double zmin,
    Eigen::VectorXd& U
) {
    // Domain center
    double xc = xmin + 0.5;
    double yc = ymin + 0.5;
    double zc = zmin + 0.5;

    // Distance from center
    double r = std::sqrt((x - xc)*(x - xc) +
                        (y - yc)*(y - yc) +
                        (z - zc)*(z - zc));

    double gamma = 1.4;
    double rho = 1.0;
    double p;

    // Blast radius
    double blast_radius = 0.1;

    if (r < blast_radius) {
        p = 10.0;  // High pressure in blast region
    } else {
        p = 0.1;   // Ambient pressure
    }

    U(0) = rho;
    U(1) = 0.0;  // rho_u
    U(2) = 0.0;  // rho_v
    U(3) = 0.0;  // rho_w
    U(4) = p / (gamma - 1.0);  // E (internal energy only, no kinetic energy)
}

int main(int argc, char** argv) {
    std::cout << "=== FVM3D Blast Wave Example ===\n\n";

    // Create grid
    core::GridGeometry3D geom(0.0, 0.0, 0.0,  // xmin, ymin, zmin
                              1.0, 1.0, 1.0,  // Lx, Ly, Lz
                              32, 32, 32);     // nx, ny, nz
    core::Grid3D grid(geom, 2);

    std::cout << "Grid created:\n"
              << "  Domain: [" << geom.xmin << ", " << geom.xmax() << "] x "
              << "[" << geom.ymin << ", " << geom.ymax() << "] x "
              << "[" << geom.zmin << ", " << geom.zmax() << "]\n"
              << "  Cells: " << geom.nx << " x " << geom.ny << " x " << geom.nz << "\n"
              << "  dx, dy, dz = " << geom.dx << ", " << geom.dy << ", " << geom.dz << "\n"
              << "  Ghost cells: " << grid.nghost() << "\n\n";

    // Create physics object
    auto physics = std::make_shared<physics::EulerEquations3D>();

    // Create state field (conservative variables)
    core::StateField3D state(5, grid.nx_total(), grid.ny_total(), grid.nz_total());
    state.fill(0.0);

    // Initialize with blast wave
    std::cout << "Initializing blast wave...\n";
    for (int i = 0; i < grid.nx_total(); i++) {
        for (int j = 0; j < grid.ny_total(); j++) {
            for (int k = 0; k < grid.nz_total(); k++) {
                double x = grid.cell_center_x(i);
                double y = grid.cell_center_y(j);
                double z = grid.cell_center_z(k);

                Eigen::VectorXd U(5);
                initialize_blast_wave(x, y, z, geom.xmin, geom.ymin, geom.zmin, U);

                for (int v = 0; v < 5; v++) {
                    state(v, i, j, k) = U(v);
                }
            }
        }
    }

    // Apply boundary conditions
    boundary::PeriodicBC bc(physics, true, true, true);
    bc.apply(state, grid);

    std::cout << "Initial condition set.\n\n";

    // Check initial values at domain center
    int i_center = grid.i_begin() + geom.nx / 2;
    int j_center = grid.j_begin() + geom.ny / 2;
    int k_center = grid.k_begin() + geom.nz / 2;

    Eigen::VectorXd U_center(5);
    for (int v = 0; v < 5; v++) {
        U_center(v) = state(v, i_center, j_center, k_center);
    }
    
    Eigen::VectorXd dummy(5);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Values at domain center:\n";
    std::cout << "  rho = " << U_center(0) << "\n";
    std::cout << "  rho_u = " << U_center(1) << "\n";
    std::cout << "  rho_v = " << U_center(2) << "\n";
    std::cout << "  rho_w = " << U_center(3) << "\n";
    std::cout << "  E = " << U_center(4) << "\n";

    // Test all flux calculators
    std::cout << "\nTesting Flux Calculators:\n";

    Eigen::VectorXd U_L(5), U_R(5);
    U_L << 1.0, 0.0, 0.0, 0.0, 2.5;      // Sod shock tube left state
    U_R << 0.125, 0.0, 0.0, 0.0, 0.25;   // Sod shock tube right state

    // Test using factory pattern (unified flux calculator interface)
    std::vector<std::string> flux_calc_names = {"laxfriedrichs", "hll", "hllc"};

    for (const auto& calc_name : flux_calc_names) {
        auto flux_calc = spatial::FluxCalculatorFactory::create(calc_name, physics);
        Eigen::VectorXd flux = flux_calc->compute_flux(U_L, U_R, *physics, 0);
        double wave_speed = flux_calc->compute_max_wave_speed(U_L, U_R, *physics, 0);

        std::cout << "  " << flux_calc->name() << ":\n";
        std::cout << "    Flux = [" << flux.transpose() << "]\n";
        std::cout << "    Max wave speed = " << wave_speed << "\n";
    }

    // Print available flux calculators
    std::cout << "\n  Available flux calculators:\n";
    spatial::FluxCalculatorFactory::print_available();

    // Test time integrators
    std::cout << "\nTesting Time Integrators:\n";

    temporal::ForwardEuler euler_integrator;
    std::cout << "  Forward Euler: order = " << euler_integrator.order() << "\n";

    // Create test field for ODE: du/dt = -u
    core::StateField3D test_field(1, 3, 3, 3);
    test_field.fill(1.0);

    auto rhs_func = [](const core::StateField3D& U, core::StateField3D& dUdt) {
        // Simple exponential decay: du/dt = -u
        dUdt.fill(0.0);
        for (int i = 0; i < U.nx(); i++) {
            for (int j = 0; j < U.ny(); j++) {
                for (int k = 0; k < U.nz(); k++) {
                    dUdt(0, i, j, k) = -U(0, i, j, k);
                }
            }
        }
    };

    double dt = 0.01;
    euler_integrator.step(test_field, dt, rhs_func);

    std::cout << "  After one step with dt=" << dt << ":\n";
    std::cout << "  test_field(0,1,1,1) = " << test_field(0, 1, 1, 1) << "\n";
    std::cout << "  Expected (exp(-dt)) = " << std::exp(-dt) << "\n";

    std::cout << "\n=== Test Completed Successfully ===\n";

    return 0;
}
