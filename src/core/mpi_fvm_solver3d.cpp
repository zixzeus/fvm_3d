#include "core/mpi_fvm_solver3d.hpp"
#include "spatial/riemann_solver_factory.hpp"
#include "spatial/reconstruction.hpp"
#include "temporal/time_integrator_factory.hpp"
#include "boundary/periodic_bc.hpp"
#include "boundary/reflective_bc.hpp"
#include "boundary/transmissive_bc.hpp"
#include "io/mpi_hdf5_checkpoint.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

namespace fvm3d::core {

MPIFVMSolver3D::MPIFVMSolver3D(const MPIFVMSolverConfig& config)
    : config_(config),
      t_current_(0.0),
      step_count_(0),
      dt_(0.0)
{
    // Get MPI info
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    // Create domain decomposer
    decomposer_ = std::make_unique<parallel::MPIDomainDecomposer>(
        config.nx, config.ny, config.nz,
        config.px, config.py, config.pz
    );

    // Get local dimensions from decomposer
    auto local_cells = decomposer_->local_interior_cells();
    int local_nx = local_cells[0];
    int local_ny = local_cells[1];
    int local_nz = local_cells[2];

    // Get global index range for this process
    auto idx_range = decomposer_->global_index_range();

    // Compute local domain geometry
    double dx = config.Lx / config.nx;
    double dy = config.Ly / config.ny;
    double dz = config.Lz / config.nz;

    double local_xmin = config.xmin + idx_range.i_min * dx;
    double local_ymin = config.ymin + idx_range.j_min * dy;
    double local_zmin = config.zmin + idx_range.k_min * dz;

    double local_Lx = local_nx * dx;
    double local_Ly = local_ny * dy;
    double local_Lz = local_nz * dz;

    // Create local grid
    GridGeometry3D local_geom(
        local_xmin, local_ymin, local_zmin,
        local_Lx, local_Ly, local_Lz,
        local_nx, local_ny, local_nz
    );
    local_grid_ = std::make_unique<Grid3D>(local_geom, config.nghost);

    // Get total dimensions (including ghost cells)
    int nx_total = local_grid_->nx_total();
    int ny_total = local_grid_->ny_total();
    int nz_total = local_grid_->nz_total();

    // Create state fields
    state_ = std::make_unique<StateField3D>(config.num_vars, nx_total, ny_total, nz_total);
    rhs_ = std::make_unique<StateField3D>(config.num_vars, nx_total, ny_total, nz_total);
    u_temp_ = std::make_unique<StateField3D>(config.num_vars, nx_total, ny_total, nz_total);
    flux_x_ = std::make_unique<StateField3D>(config.num_vars, nx_total, ny_total, nz_total);
    flux_y_ = std::make_unique<StateField3D>(config.num_vars, nx_total, ny_total, nz_total);
    flux_z_ = std::make_unique<StateField3D>(config.num_vars, nx_total, ny_total, nz_total);

    // Create MPI communication objects
    halo_exchange_ = std::make_unique<parallel::MPIHaloExchange>(
        std::make_shared<parallel::MPIDomainDecomposer>(*decomposer_),
        config.num_vars
    );

    global_reduction_ = std::make_unique<parallel::MPIGlobalReduction>();

    // Initialize Riemann solver
    riemann_solver_ = spatial::RiemannSolverFactory::create(config.riemann_solver);

    // Initialize reconstruction scheme
    reconstruction_ = spatial::ReconstructionFactory::create(
        config.reconstruction,
        config.reconstruction_limiter
    );

    // Initialize time integrator
    time_integrator_ = temporal::TimeIntegratorFactory::create(config.time_integrator);

    // Initialize boundary conditions
    if (config.boundary_condition == "periodic") {
        boundary_condition_ = std::make_unique<boundary::PeriodicBC>(
            config.bc_x, config.bc_y, config.bc_z
        );
    } else if (config.boundary_condition == "reflective") {
        boundary_condition_ = std::make_unique<boundary::ReflectiveBC>(
            config.bc_x, config.bc_y, config.bc_z
        );
    } else if (config.boundary_condition == "transmissive") {
        boundary_condition_ = std::make_unique<boundary::TransmissiveBC>(
            config.bc_x, config.bc_y, config.bc_z
        );
    } else {
        throw std::invalid_argument("Unknown boundary condition: " + config.boundary_condition);
    }

    // Print info (rank 0 only)
    if (rank_ == 0 && config.verbose > 0) {
        auto proc_grid = decomposer_->process_grid();
        std::cout << "MPI-Parallel FVM3D Solver initialized:\n"
                  << "  MPI processes: " << size_ << " ("
                  << proc_grid[0] << " x " << proc_grid[1] << " x " << proc_grid[2] << ")\n"
                  << "  Global grid: " << config.nx << " x " << config.ny << " x " << config.nz << "\n"
                  << "  Local grid (per process): " << local_nx << " x " << local_ny << " x " << local_nz << "\n"
                  << "  Domain: [" << config.xmin << "," << (config.xmin + config.Lx) << "] x "
                  << "[" << config.ymin << "," << (config.ymin + config.Ly) << "] x "
                  << "[" << config.zmin << "," << (config.zmin + config.Lz) << "]\n"
                  << "  Variables: " << config.num_vars << " (" << config.physics_type << ")\n"
                  << "  Riemann Solver: " << riemann_solver_->name() << "\n"
                  << "  Reconstruction: " << reconstruction_->name() << "\n"
                  << "  Time Integrator: " << time_integrator_->name() << "\n"
                  << "  Boundary Condition: " << boundary_condition_->name() << "\n"
                  << "  CFL: " << config.cfl << "\n";
    }
}

void MPIFVMSolver3D::initialize(const InitFunction& init_func) {
    int nx_total = local_grid_->nx_total();
    int ny_total = local_grid_->ny_total();
    int nz_total = local_grid_->nz_total();

    std::cout << "Rank " << rank_ << ": Starting initialization with grid "
              << nx_total << "x" << ny_total << "x" << nz_total << std::endl;

    // Set initial conditions at all local grid points (including ghosts)
    for (int i = 0; i < nx_total; i++) {
        for (int j = 0; j < ny_total; j++) {
            for (int k = 0; k < nz_total; k++) {
                double x = local_grid_->cell_center_x(i);
                double y = local_grid_->cell_center_y(j);
                double z = local_grid_->cell_center_z(k);

                Eigen::VectorXd U(config_.num_vars);
                init_func(x, y, z, U);

                for (int v = 0; v < config_.num_vars; v++) {
                    (*state_)(v, i, j, k) = U(v);
                }
            }
        }
    }

    std::cout << "Rank " << rank_ << ": Finished setting values, starting halo exchange" << std::endl;

    // Exchange halo data to fill ghost cells
    halo_exchange_->exchange(*state_);

    std::cout << "Rank " << rank_ << ": Halo exchange complete" << std::endl;

    if (rank_ == 0 && config_.verbose > 0) {
        std::cout << "Initial conditions set on all ranks\n";
    }
}

void MPIFVMSolver3D::run() {
    if (rank_ == 0 && config_.verbose > 0) {
        std::cout << "\nStarting time integration...\n\n";
    }

    while (step_count_ < config_.num_steps && t_current_ < config_.t_final) {
        step();

        // Output
        if (config_.output_interval > 0 && step_count_ % config_.output_interval == 0) {
            compute_statistics();
            if (rank_ == 0) {
                print_progress();
            }
        }

        // Checkpoint
        if (config_.checkpoint_interval > 0 && step_count_ % config_.checkpoint_interval == 0) {
            if (rank_ == 0 && config_.verbose > 0) {
                std::cout << "Saving checkpoint at step " << step_count_ << "...\n";
            }
            save_checkpoint("checkpoint_step_" + std::to_string(step_count_) + ".h5");
        }
    }

    if (rank_ == 0 && config_.verbose > 0) {
        std::cout << "\nSimulation complete!\n";
        std::cout << "  Final time: " << t_current_ << "\n";
        std::cout << "  Total steps: " << step_count_ << "\n";
    }
}

void MPIFVMSolver3D::step() {
    // Compute time step (global CFL constraint)
    dt_ = compute_dt();

    // Ensure we don't overshoot final time
    if (t_current_ + dt_ > config_.t_final) {
        dt_ = config_.t_final - t_current_;
    }

    // Define RHS function for time integrator
    auto rhs_function = [this](const StateField3D& U_in, StateField3D& rhs_out) {
        // Copy state
        for (int v = 0; v < config_.num_vars; v++) {
            for (int i = 0; i < local_grid_->nx_total(); i++) {
                for (int j = 0; j < local_grid_->ny_total(); j++) {
                    for (int k = 0; k < local_grid_->nz_total(); k++) {
                        (*state_)(v, i, j, k) = U_in(v, i, j, k);
                    }
                }
            }
        }

        // Exchange halos
        halo_exchange_->exchange(*state_);

        // Apply boundary conditions
        apply_boundary_conditions();

        // Compute RHS
        compute_rhs();

        // Copy rhs
        for (int v = 0; v < config_.num_vars; v++) {
            for (int i = 0; i < local_grid_->nx_total(); i++) {
                for (int j = 0; j < local_grid_->ny_total(); j++) {
                    for (int k = 0; k < local_grid_->nz_total(); k++) {
                        rhs_out(v, i, j, k) = (*rhs_)(v, i, j, k);
                    }
                }
            }
        }
    };

    // Time integration step
    time_integrator_->step(*state_, dt_, rhs_function);

    // Update time and step counter
    t_current_ += dt_;
    step_count_++;
}

void MPIFVMSolver3D::compute_rhs() {
    // Zero out RHS
    rhs_->fill(0.0);

    // Compute fluxes in each direction and accumulate flux divergence
    compute_fluxes(0, *flux_x_);
    add_flux_divergence(*flux_x_, 0);

    compute_fluxes(1, *flux_y_);
    add_flux_divergence(*flux_y_, 1);

    compute_fluxes(2, *flux_z_);
    add_flux_divergence(*flux_z_, 2);
}

void MPIFVMSolver3D::compute_fluxes(int direction, StateField3D& flux_out) {
    int nx = local_grid_->nx_local();
    int ny = local_grid_->ny_local();
    int nz = local_grid_->nz_local();
    int ng = local_grid_->nghost();

    // Loop over interior cell interfaces (need ng+1 interfaces for ng interior cells)
    int imax = (direction == 0) ? nx + 1 : nx;
    int jmax = (direction == 1) ? ny + 1 : ny;
    int kmax = (direction == 2) ? nz + 1 : nz;

    for (int i = 0; i <= imax; i++) {
        for (int j = 0; j <= jmax; j++) {
            for (int k = 0; k <= kmax; k++) {
                // Interface indices (in total grid including ghosts)
                int ii = i + ng;
                int jj = j + ng;
                int kk = k + ng;

                // Reconstruct left and right states
                Eigen::VectorXd U_L(config_.num_vars);
                Eigen::VectorXd U_R(config_.num_vars);

                reconstruct_1d(*state_, direction, ii, jj, kk, U_L, U_R);

                // Compute Riemann flux
                Eigen::VectorXd F = riemann_solver_->solve(U_L, U_R, direction);

                // Store flux
                for (int v = 0; v < config_.num_vars; v++) {
                    flux_out(v, ii, jj, kk) = F(v);
                }
            }
        }
    }
}

void MPIFVMSolver3D::add_flux_divergence(const StateField3D& flux, int direction) {
    int nx = local_grid_->nx_local();
    int ny = local_grid_->ny_local();
    int nz = local_grid_->nz_local();
    int ng = local_grid_->nghost();

    const auto& geom = local_grid_->geometry();
    double dx = geom.dx;
    double dy = geom.dy;
    double dz = geom.dz;
    double inv_dx = (direction == 0) ? -1.0 / dx : ((direction == 1) ? -1.0 / dy : -1.0 / dz);

    // Loop over interior cells
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int ii = i + ng;
                int jj = j + ng;
                int kk = k + ng;

                // Flux difference
                for (int v = 0; v < config_.num_vars; v++) {
                    double flux_diff = 0.0;
                    if (direction == 0) {
                        flux_diff = flux(v, ii+1, jj, kk) - flux(v, ii, jj, kk);
                    } else if (direction == 1) {
                        flux_diff = flux(v, ii, jj+1, kk) - flux(v, ii, jj, kk);
                    } else {
                        flux_diff = flux(v, ii, jj, kk+1) - flux(v, ii, jj, kk);
                    }

                    (*rhs_)(v, ii, jj, kk) += inv_dx * flux_diff;
                }
            }
        }
    }
}

void MPIFVMSolver3D::reconstruct_1d(
    const StateField3D& state,
    int direction,
    int i, int j, int k,
    Eigen::VectorXd& U_L,
    Eigen::VectorXd& U_R
) {
    // Extract stencil values for each variable
    for (int v = 0; v < config_.num_vars; v++) {
        double Um2, Um1, U0, Up1, Up2;

        if (direction == 0) {
            Um2 = state(v, i-2, j, k);
            Um1 = state(v, i-1, j, k);
            U0 = state(v, i, j, k);
            Up1 = state(v, i+1, j, k);
            Up2 = state(v, i+2, j, k);
        } else if (direction == 1) {
            Um2 = state(v, i, j-2, k);
            Um1 = state(v, i, j-1, k);
            U0 = state(v, i, j, k);
            Up1 = state(v, i, j+1, k);
            Up2 = state(v, i, j+2, k);
        } else {
            Um2 = state(v, i, j, k-2);
            Um1 = state(v, i, j, k-1);
            U0 = state(v, i, j, k);
            Up1 = state(v, i, j, k+1);
            Up2 = state(v, i, j, k+2);
        }

        // Reconstruct using scalar interface
        double U_L_scalar, U_R_scalar;
        reconstruction_->reconstruct(Um2, Um1, U0, Up1, Up2, U_L_scalar, U_R_scalar);

        U_L(v) = U_L_scalar;
        U_R(v) = U_R_scalar;
    }
}

void MPIFVMSolver3D::apply_boundary_conditions() {
    boundary_condition_->apply(*state_, *local_grid_);
}

double MPIFVMSolver3D::compute_dt() {
    int nx = local_grid_->nx_local();
    int ny = local_grid_->ny_local();
    int nz = local_grid_->nz_local();
    int ng = local_grid_->nghost();

    const auto& geom = local_grid_->geometry();
    double dx = geom.dx;
    double dy = geom.dy;
    double dz = geom.dz;

    double min_dt = 1e10;

    // Compute local minimum timestep
    for (int i = ng; i < nx + ng; i++) {
        for (int j = ng; j < ny + ng; j++) {
            for (int k = ng; k < nz + ng; k++) {
                // Extract state
                Eigen::VectorXd U(config_.num_vars);
                for (int v = 0; v < config_.num_vars; v++) {
                    U(v) = (*state_)(v, i, j, k);
                }

                // Compute maximum wave speed (this is physics-dependent)
                // For now, use a simple estimate based on density and energy
                double rho = U(0);
                double max_speed = 1.0;  // Placeholder

                if (rho > 1e-10) {
                    // Simple estimate for compressible flow
                    double rho_u = U(1);
                    double rho_v = U(2);
                    double rho_w = U(3);
                    double E = U(4);

                    double u = rho_u / rho;
                    double v = rho_v / rho;
                    double w = rho_w / rho;
                    double ke = 0.5 * (u*u + v*v + w*w);
                    double p = (1.4 - 1.0) * (E - rho * ke);

                    if (p > 0) {
                        double cs = std::sqrt(1.4 * p / rho);
                        double vel_mag = std::sqrt(u*u + v*v + w*w);
                        max_speed = vel_mag + cs;
                    }
                }

                double dt_cell = std::min({dx, dy, dz}) / (max_speed + 1e-10);
                min_dt = std::min(min_dt, dt_cell);
            }
        }
    }

    // Global reduction to get minimum across all ranks
    double global_min_dt = global_reduction_->reduce_cfl_timestep(min_dt);

    return config_.cfl * global_min_dt;
}

void MPIFVMSolver3D::compute_statistics() {
    int nx = local_grid_->nx_local();
    int ny = local_grid_->ny_local();
    int nz = local_grid_->nz_local();
    int ng = local_grid_->nghost();

    // Compute local min/max
    double local_min_rho = 1e10;
    double local_max_rho = -1e10;

    for (int i = ng; i < nx + ng; i++) {
        for (int j = ng; j < ny + ng; j++) {
            for (int k = ng; k < nz + ng; k++) {
                double rho = (*state_)(0, i, j, k);
                local_min_rho = std::min(local_min_rho, rho);
                local_max_rho = std::max(local_max_rho, rho);
            }
        }
    }

    // Global reductions
    stats_.min_rho = global_reduction_->global_min(local_min_rho);
    stats_.max_rho = global_reduction_->global_max(local_max_rho);
}

void MPIFVMSolver3D::print_progress() {
    std::cout << "Step " << std::setw(6) << step_count_
              << " | t = " << std::scientific << std::setprecision(4) << t_current_
              << " | dt = " << dt_
              << " | rho: [" << stats_.min_rho << ", " << stats_.max_rho << "]"
              << "\n";
}

void MPIFVMSolver3D::save_checkpoint(const std::string& filename, const std::string& description) {
    // MPIHDFCheckpoint has static methods
    io::MPIHDFCheckpoint::save(
        filename, *state_, *local_grid_, *decomposer_,
        t_current_, step_count_, description, decomposer_->cartesian_comm()
    );
}

bool MPIFVMSolver3D::load_checkpoint(const std::string& filename) {
    // MPIHDFCheckpoint has static methods
    return io::MPIHDFCheckpoint::load(
        filename, *state_, *local_grid_, *decomposer_,
        t_current_, step_count_, decomposer_->cartesian_comm()
    );
}

} // namespace fvm3d::core
