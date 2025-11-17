#include "core/fvm_solver3d.hpp"
#include "io/hdf5_checkpoint.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

namespace fvm3d::core {

FVMSolver3D::FVMSolver3D(const FVMSolverConfig& config)
    : config_(config),
      grid_(GridGeometry3D(config.xmin, config.ymin, config.zmin,
                          config.Lx, config.Ly, config.Lz,
                          config.nx, config.ny, config.nz),
            config.nghost),
      state_(config.num_vars, grid_.nx_total(), grid_.ny_total(), grid_.nz_total()),
      rhs_(config.num_vars, grid_.nx_total(), grid_.ny_total(), grid_.nz_total()),
      u_temp_(config.num_vars, grid_.nx_total(), grid_.ny_total(), grid_.nz_total()),
      flux_x_(config.num_vars, grid_.nx_total(), grid_.ny_total(), grid_.nz_total()),
      flux_y_(config.num_vars, grid_.nx_total(), grid_.ny_total(), grid_.nz_total()),
      flux_z_(config.num_vars, grid_.nx_total(), grid_.ny_total(), grid_.nz_total()),
      t_current_(0.0),
      step_count_(0),
      dt_(0.0),
      max_wave_speed_x_(0.0),
      max_wave_speed_y_(0.0),
      max_wave_speed_z_(0.0)
{
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

    if (config.verbose > 0) {
        std::cout << "FVM3D Solver initialized:\n"
                  << "  Grid: " << config.nx << " x " << config.ny << " x " << config.nz << "\n"
                  << "  Domain: [" << config.xmin << "," << (config.xmin + config.Lx) << "] x "
                  << "[" << config.ymin << "," << (config.ymin + config.Ly) << "] x "
                  << "[" << config.zmin << "," << (config.zmin + config.Lz) << "]\n"
                  << "  Flux Calculator: " << flux_calculator_->name() << "\n"
                  << "  Reconstruction: " << reconstruction_->name() << "\n"
                  << "  Time Integrator: " << time_integrator_->name() << "\n"
                  << "  Boundary Condition: " << boundary_condition_->name() << "\n"
                  << "  CFL: " << config.cfl << "\n";
    }
}

void FVMSolver3D::initialize(const InitFunction& init_func) {
    int nx_total = grid_.nx_total();
    int ny_total = grid_.ny_total();
    int nz_total = grid_.nz_total();
    int nvars = config_.num_vars;

    // Set initial conditions at all grid points
    for (int i = 0; i < nx_total; i++) {
        for (int j = 0; j < ny_total; j++) {
            for (int k = 0; k < nz_total; k++) {
                double x = grid_.cell_center_x(i);
                double y = grid_.cell_center_y(j);
                double z = grid_.cell_center_z(k);

                Eigen::VectorXd U(nvars);
                init_func(x, y, z, U);

                for (int v = 0; v < nvars; v++) {
                    state_(v, i, j, k) = U(v);
                }
            }
        }
    }

    // Apply boundary conditions
    apply_boundary_conditions();

    if (config_.verbose > 0) {
        std::cout << "Initial condition applied.\n";
    }
}

void FVMSolver3D::run() {
    if (config_.verbose > 0) {
        std::cout << "\nStarting main time loop...\n"
                  << "Final time: " << config_.t_final << "\n"
                  << "Max steps: " << config_.num_steps << "\n\n";
    }

    while (step_count_ < config_.num_steps && t_current_ < config_.t_final) {
        step();

        // Print progress
        if (config_.verbose > 0 && step_count_ % config_.output_interval == 0) {
            print_progress();
        }

        // Check for early termination
        if (std::isnan(t_current_) || std::isinf(t_current_)) {
            std::cerr << "Error: NaN or Inf detected in simulation. Aborting.\n";
            break;
        }
    }

    if (config_.verbose > 0) {
        std::cout << "\nSimulation completed.\n"
                  << "Total steps: " << step_count_ << "\n"
                  << "Final time: " << std::scientific << std::setprecision(6) << t_current_ << "\n";
    }
}

void FVMSolver3D::step() {
    // Compute adaptive time step
    dt_ = compute_dt();

    // Ensure we don't overstep
    if (t_current_ + dt_ > config_.t_final) {
        dt_ = config_.t_final - t_current_;
    }

    // Time integration using RHS callback
    auto rhs_func = [this](const StateField3D& U, StateField3D& dUdt) {
        this->compute_rhs();
        dUdt.assign(this->rhs_);
    };

    time_integrator_->step(state_, dt_, rhs_func);

    // Apply boundary conditions
    apply_boundary_conditions();

    // Update time
    t_current_ += dt_;
    step_count_++;

    // Compute statistics
    compute_statistics();
}

void FVMSolver3D::compute_rhs() {
    // Zero out RHS
    rhs_.fill(0.0);

    // Reset cached wave speeds
    max_wave_speed_x_ = 0.0;
    max_wave_speed_y_ = 0.0;
    max_wave_speed_z_ = 0.0;

    // Compute flux divergence: -dF/dx - dG/dy - dH/dz
    // Also cache max wave speeds for CFL calculation
    for (int direction = 0; direction < 3; direction++) {
        StateField3D& flux_out = (direction == 0) ? flux_x_ :
                                 (direction == 1) ? flux_y_ : flux_z_;
        double& max_speed_cache = (direction == 0) ? max_wave_speed_x_ :
                                  (direction == 1) ? max_wave_speed_y_ : max_wave_speed_z_;

        compute_fluxes(direction, flux_out, max_speed_cache);
        add_flux_divergence(flux_out, direction);
    }

    // Add source terms: S(U, x, y, z)
    // (e.g., for resistive MHD: Ohmic dissipation, GLM divergence cleaning)
    const auto& geom = grid_.geometry();
    const int i_begin = grid_.i_begin();
    const int i_end = grid_.i_end();
    const int j_begin = grid_.j_begin();
    const int j_end = grid_.j_end();
    const int k_begin = grid_.k_begin();
    const int k_end = grid_.k_end();
    const int nvars = state_.nvars();

    #pragma omp parallel for collapse(2) if((i_end - i_begin) * (j_end - j_begin) > 100)
    for (int i = i_begin; i < i_end; i++) {
        for (int j = j_begin; j < j_end; j++) {
            for (int k = k_begin; k < k_end; k++) {
                // Get cell position
                double x = grid_.cell_center_x(i);
                double y = grid_.cell_center_y(j);
                double z = grid_.cell_center_z(k);

                // Get state at this cell
                Eigen::VectorXd U(nvars);
                for (int v = 0; v < nvars; v++) {
                    U(v) = state_(v, i, j, k);
                }

                // Compute source term (default returns zero for physics without sources)
                Eigen::VectorXd S = physics_->compute_source(
                    U, x, y, z, geom.dx, geom.dy, geom.dz
                );

                // Add source to RHS
                for (int v = 0; v < nvars; v++) {
                    rhs_(v, i, j, k) += S(v);
                }
            }
        }
    }
}

void FVMSolver3D::compute_fluxes(int direction, StateField3D& flux_out, double& max_wave_speed_out) {
    flux_out.fill(0.0);

    const int i_begin = grid_.i_begin();
    const int i_end = grid_.i_end();
    const int j_begin = grid_.j_begin();
    const int j_end = grid_.j_end();
    const int k_begin = grid_.k_begin();
    const int k_end = grid_.k_end();
    const int nvars = state_.nvars();

    // Initialize max wave speed for this direction
    max_wave_speed_out = 0.0;

    // OpenMP parallelization with thread-private Eigen vectors
    // Each thread gets its own vectors to avoid race conditions
    // Use reduction to compute max wave speed efficiently
    if (direction == 0) {
        // X-direction fluxes: F_{i+1/2,j,k}
        #pragma omp parallel reduction(max:max_wave_speed_out)
        {
            Eigen::VectorXd U_L(nvars), U_R(nvars), F(nvars);

            #pragma omp for collapse(2)
            for (int i = i_begin - 1; i < i_end; i++) {
                for (int j = j_begin; j < j_end; j++) {
                    for (int k = k_begin; k < k_end; k++) {
                        reconstruct_1d(state_, direction, i, j, k, U_L, U_R);
                        F = flux_calculator_->compute_flux(U_L, U_R, *physics_, direction);

                        // Compute and cache max wave speed for CFL
                        double local_speed = flux_calculator_->compute_max_wave_speed(U_L, U_R, *physics_, direction);
                        max_wave_speed_out = std::max(max_wave_speed_out, local_speed);

                        // SIMD vectorization for flux storage
                        #pragma omp simd
                        for (int v = 0; v < nvars; v++) {
                            flux_out(v, i, j, k) = F(v);
                        }
                    }
                }
            }
        }
    } else if (direction == 1) {
        // Y-direction fluxes: G_{i,j+1/2,k}
        #pragma omp parallel reduction(max:max_wave_speed_out)
        {
            Eigen::VectorXd U_L(nvars), U_R(nvars), F(nvars);

            #pragma omp for collapse(2)
            for (int i = i_begin; i < i_end; i++) {
                for (int j = j_begin - 1; j < j_end; j++) {
                    for (int k = k_begin; k < k_end; k++) {
                        reconstruct_1d(state_, direction, i, j, k, U_L, U_R);
                        F = flux_calculator_->compute_flux(U_L, U_R, *physics_, direction);

                        // Compute and cache max wave speed for CFL
                        double local_speed = flux_calculator_->compute_max_wave_speed(U_L, U_R, *physics_, direction);
                        max_wave_speed_out = std::max(max_wave_speed_out, local_speed);

                        // SIMD vectorization for flux storage
                        #pragma omp simd
                        for (int v = 0; v < nvars; v++) {
                            flux_out(v, i, j, k) = F(v);
                        }
                    }
                }
            }
        }
    } else {
        // Z-direction fluxes: H_{i,j,k+1/2}
        #pragma omp parallel reduction(max:max_wave_speed_out)
        {
            Eigen::VectorXd U_L(nvars), U_R(nvars), F(nvars);

            #pragma omp for collapse(2)
            for (int i = i_begin; i < i_end; i++) {
                for (int j = j_begin; j < j_end; j++) {
                    for (int k = k_begin - 1; k < k_end; k++) {
                        reconstruct_1d(state_, direction, i, j, k, U_L, U_R);
                        F = flux_calculator_->compute_flux(U_L, U_R, *physics_, direction);

                        // Compute and cache max wave speed for CFL
                        double local_speed = flux_calculator_->compute_max_wave_speed(U_L, U_R, *physics_, direction);
                        max_wave_speed_out = std::max(max_wave_speed_out, local_speed);

                        // SIMD vectorization for flux storage
                        #pragma omp simd
                        for (int v = 0; v < nvars; v++) {
                            flux_out(v, i, j, k) = F(v);
                        }
                    }
                }
            }
        }
    }
}

void FVMSolver3D::add_flux_divergence(const StateField3D& flux, int direction) {
    const auto& geom = grid_.geometry();
    const double inv_dx = 1.0 / geom.dx;
    const double inv_dy = 1.0 / geom.dy;
    const double inv_dz = 1.0 / geom.dz;

    const int i_begin = grid_.i_begin();
    const int i_end = grid_.i_end();
    const int j_begin = grid_.j_begin();
    const int j_end = grid_.j_end();
    const int k_begin = grid_.k_begin();
    const int k_end = grid_.k_end();
    const int nvars = rhs_.nvars();

    // Hybrid OpenMP + SIMD vectorization:
    // - OpenMP parallel for on outer loops (thread-level parallelism)
    // - SIMD on innermost loop (instruction-level parallelism)

    if (direction == 0) {
        // dF/dx: flux at right - flux at left
        #pragma omp parallel for collapse(2) if(nvars * (i_end - i_begin) * (j_end - j_begin) > 1000)
        for (int v = 0; v < nvars; v++) {
            for (int i = i_begin; i < i_end; i++) {
                for (int j = j_begin; j < j_end; j++) {
                    // Vectorize innermost loop (k direction - contiguous in memory)
                    #pragma omp simd
                    for (int k = k_begin; k < k_end; k++) {
                        rhs_(v, i, j, k) -= (flux(v, i, j, k) - flux(v, i - 1, j, k)) * inv_dx;
                    }
                }
            }
        }
    } else if (direction == 1) {
        // dG/dy
        #pragma omp parallel for collapse(2) if(nvars * (i_end - i_begin) * (j_end - j_begin) > 1000)
        for (int v = 0; v < nvars; v++) {
            for (int i = i_begin; i < i_end; i++) {
                for (int j = j_begin; j < j_end; j++) {
                    #pragma omp simd
                    for (int k = k_begin; k < k_end; k++) {
                        rhs_(v, i, j, k) -= (flux(v, i, j, k) - flux(v, i, j - 1, k)) * inv_dy;
                    }
                }
            }
        }
    } else {
        // dH/dz
        #pragma omp parallel for collapse(2) if(nvars * (i_end - i_begin) * (j_end - j_begin) > 1000)
        for (int v = 0; v < nvars; v++) {
            for (int i = i_begin; i < i_end; i++) {
                for (int j = j_begin; j < j_end; j++) {
                    #pragma omp simd
                    for (int k = k_begin; k < k_end; k++) {
                        rhs_(v, i, j, k) -= (flux(v, i, j, k) - flux(v, i, j, k - 1)) * inv_dz;
                    }
                }
            }
        }
    }
}

void FVMSolver3D::apply_boundary_conditions() {
    boundary_condition_->apply(state_, grid_);
}

double FVMSolver3D::compute_dt() {
    // Use cached wave speeds from flux computation (computed during compute_rhs())
    // This avoids redundant wave speed calculation for all cells
    const double max_speed = std::max({max_wave_speed_x_, max_wave_speed_y_, max_wave_speed_z_});

    // Safety check to avoid division by zero
    if (max_speed < 1e-10) {
        return 1e-3;  // Return small default time step
    }

    const auto& geom = grid_.geometry();
    const double min_dx = std::min({geom.dx, geom.dy, geom.dz});

    return config_.cfl * min_dx / max_speed;
}

void FVMSolver3D::compute_statistics() {
    // Initialize statistics
    double min_rho = 1e10;
    double max_rho = -1e10;
    double min_p = 1e10;
    double max_p = -1e10;
    double min_speed = 1e10;
    double max_speed = -1e10;

    const int i_begin = grid_.i_begin();
    const int i_end = grid_.i_end();
    const int j_begin = grid_.j_begin();
    const int j_end = grid_.j_end();
    const int k_begin = grid_.k_begin();
    const int k_end = grid_.k_end();
    const int nvars = config_.num_vars;

    // OpenMP parallel reduction for min/max statistics
    // Each thread computes local min/max, then reduce to global
    #pragma omp parallel reduction(min:min_rho,min_p,min_speed) reduction(max:max_rho,max_p,max_speed)
    {
        // Thread-private Eigen vectors
        Eigen::VectorXd U(nvars);
        Eigen::VectorXd V(nvars);

        #pragma omp for collapse(2)
        for (int i = i_begin; i < i_end; i++) {
            for (int j = j_begin; j < j_end; j++) {
                for (int k = k_begin; k < k_end; k++) {
                    // Load conservative variables
                    for (int v = 0; v < nvars; v++) {
                        U(v) = state_(v, i, j, k);
                    }

                    // Convert to primitive variables
                    V = physics_->conservative_to_primitive(U);
                    const double rho = V(0);
                    const double u = V(1);
                    const double v = V(2);
                    const double w = V(3);
                    const double p = V(4);

                    // Update local min/max (will be reduced across threads)
                    min_rho = std::min(min_rho, rho);
                    max_rho = std::max(max_rho, rho);
                    min_p = std::min(min_p, p);
                    max_p = std::max(max_p, p);

                    const double speed = std::sqrt(u*u + v*v + w*w);
                    min_speed = std::min(min_speed, speed);
                    max_speed = std::max(max_speed, speed);
                }
            }
        }
    }

    // Store results in stats struct
    stats_.min_rho = min_rho;
    stats_.max_rho = max_rho;
    stats_.min_p = min_p;
    stats_.max_p = max_p;
    stats_.min_speed = min_speed;
    stats_.max_speed = max_speed;
}

void FVMSolver3D::print_progress() {
    std::cout << std::scientific << std::setprecision(6)
              << "Step " << std::setw(6) << step_count_
              << " | Time " << std::setw(12) << t_current_
              << " | dt " << std::setw(12) << dt_
              << " | rho:[" << stats_.min_rho << "," << stats_.max_rho << "]"
              << " | p:[" << stats_.min_p << "," << stats_.max_p << "]\n";
}

// Note: reconstruct_1d() is now inherited from FVMSolverBase

void FVMSolver3D::save_checkpoint(
    const std::string& filename,
    const std::string& description
) {
    io::HDF5Checkpoint checkpoint(filename);
    checkpoint.save(state_, grid_, t_current_, step_count_, description);
}

bool FVMSolver3D::load_checkpoint(const std::string& filename) {
    io::HDF5Checkpoint checkpoint(filename);
    return checkpoint.load(state_, t_current_, step_count_);
}

} // namespace fvm3d::core
