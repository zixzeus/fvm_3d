#include "core/mpi_fvm_solver3d.hpp"
#include "physics/resistive_mhd3d_advanced.hpp"
#include "io/mpi_hdf5_checkpoint.hpp"
#include "parallel/mpi_utils.hpp"
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

    // Initialize physics using base class method
    initialize_physics(
        config.physics_type,
        config.mhd_resistivity,
        config.mhd_eta0,
        config.mhd_eta1,
        config.mhd_localization_scale,
        config.mhd_glm_ch,
        config.mhd_glm_cr
    );

    // Initialize flux calculator using base class method
    initialize_flux_calculator(config.flux_calculator);

    // Initialize reconstruction scheme using base class method
    initialize_reconstruction(
        config.reconstruction,
        config.num_vars,
        config.reconstruction_limiter
    );

    // Initialize time integrator using base class method
    initialize_time_integrator(config.time_integrator);

    // Initialize boundary conditions using base class method
    initialize_boundary_conditions(
        config.boundary_condition,
        config.bc_x, config.bc_y, config.bc_z
    );

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
                  << "  Riemann Solver: " << flux_calculator_->name() << "\n"
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

    // ================================================================
    // STRANG SPLITTING: D^(dt/2) ∘ L^(dt) ∘ D^(dt/2)
    // ================================================================
    // Reference: OpenMHD glm_ss2.f90, Dedner et al. (2002)
    // ================================================================

    // ========== GLM DAMPING: First half-step ==========
    auto mhd = std::dynamic_pointer_cast<physics::AdvancedResistiveMHD3D>(physics_);
    if (mhd) {
        const int nx = local_grid_->nx_total();
        const int ny = local_grid_->ny_total();
        const int nz = local_grid_->nz_total();
        const int nvars = config_.num_vars;

        #pragma omp parallel
        {
            Eigen::VectorXd U(nvars);

            #pragma omp for collapse(3)
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    for (int k = 0; k < nz; k++) {
                        // Extract state vector at (i,j,k)
                        for (int v = 0; v < nvars; v++) {
                            U(v) = (*state_)(v, i, j, k);
                        }

                        // Apply GLM damping
                        mhd->apply_glm_damping(U, dt_, 0.5);

                        // Write back
                        for (int v = 0; v < nvars; v++) {
                            (*state_)(v, i, j, k) = U(v);
                        }
                    }
                }
            }
        }
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

    // ========== MHD EVOLUTION: Full time step (RK2/RK3) ==========
    time_integrator_->step(*state_, dt_, rhs_function);

    // ========== GLM DAMPING: Second half-step ==========
    if (mhd) {
        const int nx = local_grid_->nx_total();
        const int ny = local_grid_->ny_total();
        const int nz = local_grid_->nz_total();
        const int nvars = config_.num_vars;

        #pragma omp parallel
        {
            Eigen::VectorXd U(nvars);

            #pragma omp for collapse(3)
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    for (int k = 0; k < nz; k++) {
                        // Extract state vector at (i,j,k)
                        for (int v = 0; v < nvars; v++) {
                            U(v) = (*state_)(v, i, j, k);
                        }

                        // Apply GLM damping
                        mhd->apply_glm_damping(U, dt_, 0.5);

                        // Write back
                        for (int v = 0; v < nvars; v++) {
                            (*state_)(v, i, j, k) = U(v);
                        }
                    }
                }
            }
        }
    }

    // ========== POSITIVITY LIMITER ==========
    // Ensure ρ > 0 and p > 0 after time evolution
    apply_positivity_limiter();

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

    // Add source terms: S(U, x, y, z)
    const auto& geom = local_grid_->geometry();
    const int i_begin = local_grid_->i_begin();
    const int i_end = local_grid_->i_end();
    const int j_begin = local_grid_->j_begin();
    const int j_end = local_grid_->j_end();
    const int k_begin = local_grid_->k_begin();
    const int k_end = local_grid_->k_end();
    const int nvars = state_->nvars();

    #pragma omp parallel for collapse(2) if((i_end - i_begin) * (j_end - j_begin) > 100)
    for (int i = i_begin; i < i_end; i++) {
        for (int j = j_begin; j < j_end; j++) {
            for (int k = k_begin; k < k_end; k++) {
                // Get cell position
                double x = local_grid_->cell_center_x(i);
                double y = local_grid_->cell_center_y(j);
                double z = local_grid_->cell_center_z(k);

                // Get state at this cell
                Eigen::VectorXd U(nvars);
                for (int v = 0; v < nvars; v++) {
                    U(v) = (*state_)(v, i, j, k);
                }

                // Compute source term
                Eigen::VectorXd S = physics_->compute_source(
                    U, x, y, z, geom.dx, geom.dy, geom.dz
                );

                // Add source to RHS
                for (int v = 0; v < nvars; v++) {
                    (*rhs_)(v, i, j, k) += S(v);
                }
            }
        }
    }
}

void MPIFVMSolver3D::compute_fluxes(int direction, StateField3D& flux_out) {
    const int nx = local_grid_->nx_local();
    const int ny = local_grid_->ny_local();
    const int nz = local_grid_->nz_local();
    const int ng = local_grid_->nghost();
    const int nvars = config_.num_vars;

    // Optimization: reuse Eigen vectors to reduce allocations
    Eigen::VectorXd U_L(nvars), U_R(nvars), F(nvars);

    // Loop over interior cells to compute fluxes at their interfaces
    // MUSCL needs stencil [i-1, i, i+1, i+2], so we can safely compute
    // fluxes from i=ng to i=ng+nx-1 (which accesses up to ng+nx+1)
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                // Interface indices (in total grid including ghosts)
                const int ii = i + ng;
                const int jj = j + ng;
                const int kk = k + ng;

                // Reconstruct left and right states
                reconstruct_1d(*state_, direction, ii, jj, kk, U_L, U_R);

                // Compute numerical flux
                F = flux_calculator_->compute_flux(U_L, U_R, *physics_, direction);

                // Store flux with SIMD vectorization
                #pragma omp simd
                for (int v = 0; v < nvars; v++) {
                    flux_out(v, ii, jj, kk) = F(v);
                }
            }
        }
    }
}

void MPIFVMSolver3D::add_flux_divergence(const StateField3D& flux, int direction) {
    const int nx = local_grid_->nx_local();
    const int ny = local_grid_->ny_local();
    const int nz = local_grid_->nz_local();
    const int ng = local_grid_->nghost();
    const int nvars = config_.num_vars;

    const auto& geom = local_grid_->geometry();
    const double dx = geom.dx;
    const double dy = geom.dy;
    const double dz = geom.dz;
    const double inv_dx = (direction == 0) ? -1.0 / dx : ((direction == 1) ? -1.0 / dy : -1.0 / dz);

    // Vectorized approach: reorder loops to exploit SoA memory layout
    // Process each variable separately to enable SIMD vectorization

    if (direction == 0) {
        // X-direction flux divergence
        for (int v = 0; v < nvars; v++) {
            for (int i = 0; i < nx; i++) {
                const int ii = i + ng;
                for (int j = 0; j < ny; j++) {
                    const int jj = j + ng;
                    // Vectorize innermost loop (k direction - contiguous in memory)
                    #pragma omp simd
                    for (int k = 0; k < nz; k++) {
                        const int kk = k + ng;
                        const double flux_diff = flux(v, ii+1, jj, kk) - flux(v, ii, jj, kk);
                        (*rhs_)(v, ii, jj, kk) += inv_dx * flux_diff;
                    }
                }
            }
        }
    } else if (direction == 1) {
        // Y-direction flux divergence
        for (int v = 0; v < nvars; v++) {
            for (int i = 0; i < nx; i++) {
                const int ii = i + ng;
                for (int j = 0; j < ny; j++) {
                    const int jj = j + ng;
                    #pragma omp simd
                    for (int k = 0; k < nz; k++) {
                        const int kk = k + ng;
                        const double flux_diff = flux(v, ii, jj+1, kk) - flux(v, ii, jj, kk);
                        (*rhs_)(v, ii, jj, kk) += inv_dx * flux_diff;
                    }
                }
            }
        }
    } else {
        // Z-direction flux divergence
        for (int v = 0; v < nvars; v++) {
            for (int i = 0; i < nx; i++) {
                const int ii = i + ng;
                for (int j = 0; j < ny; j++) {
                    const int jj = j + ng;
                    #pragma omp simd
                    for (int k = 0; k < nz; k++) {
                        const int kk = k + ng;
                        const double flux_diff = flux(v, ii, jj, kk+1) - flux(v, ii, jj, kk);
                        (*rhs_)(v, ii, jj, kk) += inv_dx * flux_diff;
                    }
                }
            }
        }
    }
}

// Note: reconstruct_1d() is now inherited from FVMSolverBase

void MPIFVMSolver3D::apply_boundary_conditions() {
    boundary_condition_->apply(*state_, *local_grid_);
}

double MPIFVMSolver3D::compute_dt() {
    const int nx = local_grid_->nx_local();
    const int ny = local_grid_->ny_local();
    const int nz = local_grid_->nz_local();
    const int ng = local_grid_->nghost();
    const int nvars = config_.num_vars;

    const auto& geom = local_grid_->geometry();
    const double dx = geom.dx;
    const double dy = geom.dy;
    const double dz = geom.dz;
    const double min_dx = std::min({dx, dy, dz});

    double min_dt = 1e10;
    double max_wave_speed_global = 0.0;
    int debug_count = 0;

    // Optimization: separate paths for MHD and Euler to reduce branching
    const bool is_mhd = (nvars >= 8);

    // Compute local minimum timestep
    if (is_mhd) {
        // MHD path - optimized for magnetic fields
        constexpr double gamma_mhd = 5.0/3.0;
        constexpr double gamma_m1 = gamma_mhd - 1.0;
        constexpr double p_floor = 1e-10;
        constexpr double rho_floor = 1e-10;

        for (int i = ng; i < nx + ng; i++) {
            for (int j = ng; j < ny + ng; j++) {
                for (int k = ng; k < nz + ng; k++) {
                    // Direct access to state - avoid Eigen overhead
                    const double rho = std::max((*state_)(0, i, j, k), rho_floor);
                    const double rho_u = (*state_)(1, i, j, k);
                    const double rho_v = (*state_)(2, i, j, k);
                    const double rho_w = (*state_)(3, i, j, k);
                    const double E = (*state_)(4, i, j, k);
                    const double Bx = (*state_)(5, i, j, k);
                    const double By = (*state_)(6, i, j, k);
                    const double Bz = (*state_)(7, i, j, k);

                    const double inv_rho = 1.0 / rho;
                    const double u = rho_u * inv_rho;
                    const double v = rho_v * inv_rho;
                    const double w = rho_w * inv_rho;

                    const double B_sq = Bx*Bx + By*By + Bz*Bz;
                    const double ke = 0.5 * rho * (u*u + v*v + w*w);
                    const double me = 0.5 * B_sq;
                    const double internal_energy = E - ke - me;
                    double p = gamma_m1 * internal_energy;

                    if (p <= 0) {
                        static int warning_count = 0;
                        if (warning_count < 5 && parallel::MPIUtils::is_root()) {
                            std::cout << "WARNING: Negative pressure at cell[" << i << "," << j << "," << k << "]: "
                                      << "p=" << p << " -> floor " << p_floor << std::endl;
                            warning_count++;
                        }
                        p = p_floor;
                    }

                    // Fast magnetosonic speed
                    const double a_sq = gamma_mhd * p * inv_rho;
                    const double va_sq = B_sq * inv_rho;
                    const double cf = std::sqrt(a_sq + va_sq);
                    const double vel_mag = std::sqrt(u*u + v*v + w*w);
                    const double max_speed = vel_mag + cf;

                    max_wave_speed_global = std::max(max_wave_speed_global, max_speed);
                    const double dt_cell = min_dx / (max_speed + 1e-10);
                    min_dt = std::min(min_dt, dt_cell);
                    debug_count++;
                }
            }
        }
    } else {
        // Euler path - optimized for gas dynamics
        constexpr double gamma = 1.4;
        constexpr double gamma_m1 = gamma - 1.0;
        constexpr double rho_floor = 1e-10;

        for (int i = ng; i < nx + ng; i++) {
            for (int j = ng; j < ny + ng; j++) {
                for (int k = ng; k < nz + ng; k++) {
                    const double rho = std::max((*state_)(0, i, j, k), rho_floor);
                    const double rho_u = (*state_)(1, i, j, k);
                    const double rho_v = (*state_)(2, i, j, k);
                    const double rho_w = (*state_)(3, i, j, k);
                    const double E = (*state_)(4, i, j, k);

                    const double inv_rho = 1.0 / rho;
                    const double u = rho_u * inv_rho;
                    const double v = rho_v * inv_rho;
                    const double w = rho_w * inv_rho;

                    const double ke = 0.5 * (u*u + v*v + w*w);
                    const double p = gamma_m1 * (E - rho * ke);

                    double max_speed = 1.0;
                    if (p > 0) {
                        const double cs = std::sqrt(gamma * p * inv_rho);
                        const double vel_mag = std::sqrt(u*u + v*v + w*w);
                        max_speed = vel_mag + cs;
                    }

                    max_wave_speed_global = std::max(max_wave_speed_global, max_speed);
                    const double dt_cell = min_dx / (max_speed + 1e-10);
                    min_dt = std::min(min_dt, dt_cell);
                }
            }
        }
    }

    static int debug_call_count = 0;
    if (parallel::MPIUtils::is_root() && debug_call_count < 100) {
        std::cout << "compute_dt call #" << debug_call_count << ": min_dt=" << min_dt
                  << ", max_wave_speed=" << max_wave_speed_global << std::endl;
        if (min_dt <= 0 || std::isnan(min_dt) || std::isinf(min_dt)) {
            std::cout << "  ERROR: Invalid min_dt detected!" << std::endl;
        }
        debug_call_count++;
    }

    // Global reduction to get minimum across all ranks
    const double global_min_dt = global_reduction_->reduce_cfl_timestep(min_dt);

    return config_.cfl * global_min_dt;
}

void MPIFVMSolver3D::compute_statistics() {
    int nx = local_grid_->nx_local();
    int ny = local_grid_->ny_local();
    int nz = local_grid_->nz_local();
    int ng = local_grid_->nghost();
    int nvars = config_.num_vars;

    // Compute local statistics
    double local_min_rho = 1e10;
    double local_max_rho = -1e10;
    double local_min_p = 1e10;
    double local_max_p = -1e10;
    double local_min_speed = 1e10;
    double local_max_speed = -1e10;
    double local_kinetic_energy = 0.0;
    double local_magnetic_energy = 0.0;
    double local_internal_energy = 0.0;

    const auto& geom = local_grid_->geometry();
    const double dV = geom.dx * geom.dy * geom.dz;

    Eigen::VectorXd U(nvars);
    Eigen::VectorXd V(nvars);

    for (int i = ng; i < nx + ng; i++) {
        for (int j = ng; j < ny + ng; j++) {
            for (int k = ng; k < nz + ng; k++) {
                // Load conservative variables
                for (int v = 0; v < nvars; v++) {
                    U(v) = (*state_)(v, i, j, k);
                }

                // Convert to primitive
                V = physics_->conservative_to_primitive(U);
                const double rho = V(0);
                const double u = V(1);
                const double v = V(2);
                const double w = V(3);
                const double p = V(4);

                // Update min/max
                local_min_rho = std::min(local_min_rho, rho);
                local_max_rho = std::max(local_max_rho, rho);
                local_min_p = std::min(local_min_p, p);
                local_max_p = std::max(local_max_p, p);

                const double speed = std::sqrt(u*u + v*v + w*w);
                local_min_speed = std::min(local_min_speed, speed);
                local_max_speed = std::max(local_max_speed, speed);

                // Compute energies for MHD
                if (nvars == 9) {
                    // Kinetic energy
                    local_kinetic_energy += 0.5 * rho * (u*u + v*v + w*w) * dV;

                    // Magnetic energy
                    const double Bx = V(5);
                    const double By = V(6);
                    const double Bz = V(7);
                    local_magnetic_energy += 0.5 * (Bx*Bx + By*By + Bz*Bz) * dV;

                    // Internal energy
                    const double gamma = 5.0/3.0;
                    local_internal_energy += p / (gamma - 1.0) * dV;
                }
            }
        }
    }

    // Global reductions
    stats_.min_rho = global_reduction_->global_min(local_min_rho);
    stats_.max_rho = global_reduction_->global_max(local_max_rho);
    stats_.min_p = global_reduction_->global_min(local_min_p);
    stats_.max_p = global_reduction_->global_max(local_max_p);
    stats_.min_speed = global_reduction_->global_min(local_min_speed);
    stats_.max_speed = global_reduction_->global_max(local_max_speed);
    stats_.kinetic_energy = global_reduction_->global_sum(local_kinetic_energy);
    stats_.magnetic_energy = global_reduction_->global_sum(local_magnetic_energy);
    stats_.internal_energy = global_reduction_->global_sum(local_internal_energy);
    stats_.total_energy = stats_.kinetic_energy + stats_.magnetic_energy + stats_.internal_energy;

    // Compute div(B) with global reduction
    double local_max_div_B = compute_max_div_B();
    stats_.max_div_B = global_reduction_->global_max(local_max_div_B);
}

void MPIFVMSolver3D::print_progress() {
    if (config_.num_vars == 9) {
        // MHD output with energies and div(B)
        std::cout << std::scientific << std::setprecision(3)
                  << "Step " << std::setw(6) << step_count_
                  << " | t=" << std::setw(9) << t_current_
                  << " | dt=" << std::setw(9) << dt_
                  << " | KE=" << std::setw(9) << stats_.kinetic_energy
                  << " | BE=" << std::setw(9) << stats_.magnetic_energy
                  << " | E_tot=" << std::setw(9) << stats_.total_energy
                  << " | div(B)=" << std::setw(9) << stats_.max_div_B
                  << " | p:[" << stats_.min_p << "," << stats_.max_p << "]\n";
    } else {
        // Euler output
        std::cout << std::scientific << std::setprecision(6)
                  << "Step " << std::setw(6) << step_count_
                  << " | Time " << std::setw(12) << t_current_
                  << " | dt " << std::setw(12) << dt_
                  << " | rho:[" << stats_.min_rho << "," << stats_.max_rho << "]"
                  << " | p:[" << stats_.min_p << "," << stats_.max_p << "]\n";
    }
}

double MPIFVMSolver3D::compute_max_div_B() const {
    // Only compute for MHD simulations
    if (config_.num_vars != 9) {
        return 0.0;
    }

    int nx = local_grid_->nx_local();
    int ny = local_grid_->ny_local();
    int nz = local_grid_->nz_local();
    int ng = local_grid_->nghost();

    const auto& geom = local_grid_->geometry();
    const double inv_dx = 1.0 / geom.dx;
    const double inv_dy = 1.0 / geom.dy;
    const double inv_dz = 1.0 / geom.dz;

    double max_div_B = 0.0;

    // Compute div(B) = dBx/dx + dBy/dy + dBz/dz using centered differences
    for (int i = ng + 1; i < nx + ng - 1; i++) {
        for (int j = ng + 1; j < ny + ng - 1; j++) {
            for (int k = ng + 1; k < nz + ng - 1; k++) {
                const double dBx_dx = ((*state_)(5, i+1, j, k) - (*state_)(5, i-1, j, k)) * 0.5 * inv_dx;
                const double dBy_dy = ((*state_)(6, i, j+1, k) - (*state_)(6, i, j-1, k)) * 0.5 * inv_dy;
                const double dBz_dz = ((*state_)(7, i, j, k+1) - (*state_)(7, i, j, k-1)) * 0.5 * inv_dz;

                const double div_B = std::abs(dBx_dx + dBy_dy + dBz_dz);
                max_div_B = std::max(max_div_B, div_B);
            }
        }
    }

    return max_div_B;
}

void MPIFVMSolver3D::apply_positivity_limiter() {
    int nx = local_grid_->nx_local();
    int ny = local_grid_->ny_local();
    int nz = local_grid_->nz_local();
    int ng = local_grid_->nghost();
    int nvars = config_.num_vars;

    // More physical floor values to reduce energy injection
    // Higher floors mean fewer fixes and less artificial energy added
    const double rho_floor = 0.01;   // 1% of typical density scale
    const double p_floor = 0.001;    // Small but physical pressure floor

    int local_num_fixes = 0;
    double local_energy_added = 0.0;

    Eigen::VectorXd U(nvars);
    Eigen::VectorXd V(nvars);

    for (int i = ng; i < nx + ng; i++) {
        for (int j = ng; j < ny + ng; j++) {
            for (int k = ng; k < nz + ng; k++) {
                // Load conservative variables
                for (int v = 0; v < nvars; v++) {
                    U(v) = (*state_)(v, i, j, k);
                }

                // Store energy before limiting
                double E_before = U(4);

                // Convert to primitive to check positivity
                V = physics_->conservative_to_primitive(U);

                bool needs_fix = false;

                // Check and fix density with gentler blending approach
                // Instead of hard floor, blend towards floor to minimize mass injection
                if (V(0) < rho_floor) {
                    // Blend 50% old value with 50% floor
                    V(0) = 0.5 * (V(0) + rho_floor);
                    // Safety check - ensure we're above floor
                    if (V(0) < rho_floor) {
                        V(0) = rho_floor;
                    }
                    needs_fix = true;
                }

                // Check and fix pressure with very gentle approach
                // Use 90% old value to minimize energy injection
                if (nvars >= 5 && V(4) < p_floor) {
                    double p_old = V(4);
                    // Blend 90% old with 10% floor - minimizes energy injection
                    V(4) = 0.9 * p_old + 0.1 * p_floor;
                    // Safety check - allow lower floor for extreme cases
                    if (V(4) < p_floor * 0.1) {
                        V(4) = p_floor * 0.1;
                    }
                    needs_fix = true;
                }

                // If fixes were applied, convert back to conservative
                if (needs_fix) {
                    U = physics_->primitive_to_conservative(V);

                    // Write back to state
                    for (int v = 0; v < nvars; v++) {
                        (*state_)(v, i, j, k) = U(v);
                    }

                    // Track energy injection
                    double E_after = U(4);
                    local_energy_added += (E_after - E_before);

                    local_num_fixes++;
                }
            }
        }
    }

    // Global reductions for diagnostics
    int global_num_fixes = static_cast<int>(global_reduction_->global_sum(static_cast<double>(local_num_fixes)));
    double global_energy_added = global_reduction_->global_sum(local_energy_added);

    // Report if fixes were applied (rank 0 only, and only if verbose > 1)
    // Reduced verbosity to avoid output spam - most users don't need this detail
    if (global_num_fixes > 0 && config_.verbose > 1 && parallel::MPIUtils::is_root()) {
        const auto& geom = local_grid_->geometry();
        const double dV = geom.dx * geom.dy * geom.dz;
        std::cout << "  [Positivity] Fixed " << global_num_fixes
                  << " cells, ΔE = " << std::scientific << std::setprecision(2)
                  << (global_energy_added * dV) << "\n";
    }
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
