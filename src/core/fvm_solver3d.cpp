#include "core/fvm_solver3d.hpp"
#include "physics/euler3d.hpp"
#include "physics/resistive_mhd3d_advanced.hpp"
#include "spatial/flux_calculation/flux_calculator_factory.hpp"
#include "temporal/time_integrator_factory.hpp"
#include "spatial/reconstruction/reconstruction_factory.hpp"
#include "spatial/reconstruction/reconstruction_base.hpp"
#include "boundary/periodic_bc.hpp"
#include "boundary/reflective_bc.hpp"
#include "boundary/transmissive_bc.hpp"
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
      dt_(0.0)
{
    // Initialize physics based on physics_type
    if (config.physics_type == "euler") {
        physics_ = std::make_shared<physics::EulerEquations3D>();
    } else if (config.physics_type == "mhd_advanced" || config.physics_type == "mhd") {
        physics::AdvancedResistiveMHD3D::ResistivityModel resistivity;
        resistivity.eta0 = 1e-3;
        resistivity.eta1 = 0.01667;
        resistivity.localization_scale = 1.0;
        physics::AdvancedResistiveMHD3D::GLMParameters glm(0.2, 0.2);
        physics_ = std::make_shared<physics::AdvancedResistiveMHD3D>(resistivity, glm);
    } else {
        throw std::invalid_argument("Unknown physics type: " + config.physics_type);
    }

    // Initialize flux calculator
    flux_calculator_ = spatial::FluxCalculatorFactory::create(config.flux_calculator, config.physics_type, config.num_vars);

    // Initialize reconstruction scheme
    reconstruction_ = spatial::ReconstructionFactory::create(
        config.reconstruction,
        config.num_vars,
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

    // Set initial conditions at all grid points
    for (int i = 0; i < nx_total; i++) {
        for (int j = 0; j < ny_total; j++) {
            for (int k = 0; k < nz_total; k++) {
                double x = grid_.cell_center_x(i);
                double y = grid_.cell_center_y(j);
                double z = grid_.cell_center_z(k);

                Eigen::VectorXd U(5);
                init_func(x, y, z, U);

                for (int v = 0; v < 5; v++) {
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

    // Compute flux divergence: -dF/dx - dG/dy - dH/dz
    for (int direction = 0; direction < 3; direction++) {
        StateField3D& flux_out = (direction == 0) ? flux_x_ :
                                 (direction == 1) ? flux_y_ : flux_z_;
        compute_fluxes(direction, flux_out);
        add_flux_divergence(flux_out, direction);
    }
}

void FVMSolver3D::compute_fluxes(int direction, StateField3D& flux_out) {
    flux_out.fill(0.0);

    const int i_begin = grid_.i_begin();
    const int i_end = grid_.i_end();
    const int j_begin = grid_.j_begin();
    const int j_end = grid_.j_end();
    const int k_begin = grid_.k_begin();
    const int k_end = grid_.k_end();
    const int nvars = state_.nvars();

    // OpenMP parallelization with thread-private Eigen vectors
    // Each thread gets its own vectors to avoid race conditions
    if (direction == 0) {
        // X-direction fluxes: F_{i+1/2,j,k}
        #pragma omp parallel
        {
            Eigen::VectorXd U_L(nvars), U_R(nvars), F(nvars);

            #pragma omp for collapse(2)
            for (int i = i_begin - 1; i < i_end; i++) {
                for (int j = j_begin; j < j_end; j++) {
                    for (int k = k_begin; k < k_end; k++) {
                        reconstruct_1d(state_, direction, i, j, k, U_L, U_R);
                        F = flux_calculator_->compute_flux(U_L, U_R, *physics_, direction);

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
        #pragma omp parallel
        {
            Eigen::VectorXd U_L(nvars), U_R(nvars), F(nvars);

            #pragma omp for collapse(2)
            for (int i = i_begin; i < i_end; i++) {
                for (int j = j_begin - 1; j < j_end; j++) {
                    for (int k = k_begin; k < k_end; k++) {
                        reconstruct_1d(state_, direction, i, j, k, U_L, U_R);
                        F = flux_calculator_->compute_flux(U_L, U_R, *physics_, direction);

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
        #pragma omp parallel
        {
            Eigen::VectorXd U_L(nvars), U_R(nvars), F(nvars);

            #pragma omp for collapse(2)
            for (int i = i_begin; i < i_end; i++) {
                for (int j = j_begin; j < j_end; j++) {
                    for (int k = k_begin - 1; k < k_end; k++) {
                        reconstruct_1d(state_, direction, i, j, k, U_L, U_R);
                        F = flux_calculator_->compute_flux(U_L, U_R, *physics_, direction);

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
    const int i_begin = grid_.i_begin();
    const int i_end = grid_.i_end();
    const int j_begin = grid_.j_begin();
    const int j_end = grid_.j_end();
    const int k_begin = grid_.k_begin();
    const int k_end = grid_.k_end();
    const int nvars = state_.nvars();

    double max_speed = 1e-10;

    // OpenMP parallel reduction for maximum wave speed
    // Each thread computes local max, then reduce to global max
    #pragma omp parallel reduction(max:max_speed)
    {
        // Thread-local Eigen vector to avoid conflicts
        Eigen::VectorXd U(nvars);

        #pragma omp for collapse(2)
        for (int i = i_begin; i < i_end; i++) {
            for (int j = j_begin; j < j_end; j++) {
                for (int k = k_begin; k < k_end; k++) {
                    // Load state vector - accessing contiguous memory in k direction
                    for (int v = 0; v < nvars; v++) {
                        U(v) = state_(v, i, j, k);
                    }

                    // Check wave speed in each direction
                    // Unroll direction loop manually for better optimization
                    const double speed_x = physics_->max_wave_speed(U, 0);
                    const double speed_y = physics_->max_wave_speed(U, 1);
                    const double speed_z = physics_->max_wave_speed(U, 2);

                    // Local max to reduce number of comparisons
                    const double local_max = std::max({speed_x, speed_y, speed_z});
                    max_speed = std::max(max_speed, local_max);
                }
            }
        }
    }

    const auto& geom = grid_.geometry();
    const double min_dx = std::min({geom.dx, geom.dy, geom.dz});

    return config_.cfl * min_dx / max_speed;
}

void FVMSolver3D::compute_statistics() {
    stats_.min_rho = 1e10;
    stats_.max_rho = -1e10;
    stats_.min_p = 1e10;
    stats_.max_p = -1e10;
    stats_.min_speed = 1e10;
    stats_.max_speed = -1e10;

    int i_begin = grid_.i_begin();
    int i_end = grid_.i_end();
    int j_begin = grid_.j_begin();
    int j_end = grid_.j_end();
    int k_begin = grid_.k_begin();
    int k_end = grid_.k_end();

    for (int i = i_begin; i < i_end; i++) {
        for (int j = j_begin; j < j_end; j++) {
            for (int k = k_begin; k < k_end; k++) {
                Eigen::VectorXd U(config_.num_vars);
                for (int v = 0; v < config_.num_vars; v++) {
                    U(v) = state_(v, i, j, k);
                }

                Eigen::VectorXd V = physics_->conservative_to_primitive(U);
                double rho = V(0);
                double u = V(1);
                double v = V(2);
                double w = V(3);
                double p = V(4);

                stats_.min_rho = std::min(stats_.min_rho, rho);
                stats_.max_rho = std::max(stats_.max_rho, rho);
                stats_.min_p = std::min(stats_.min_p, p);
                stats_.max_p = std::max(stats_.max_p, p);

                double speed = std::sqrt(u*u + v*v + w*w);
                stats_.min_speed = std::min(stats_.min_speed, speed);
                stats_.max_speed = std::max(stats_.max_speed, speed);
            }
        }
    }
}

void FVMSolver3D::print_progress() {
    std::cout << std::scientific << std::setprecision(6)
              << "Step " << std::setw(6) << step_count_
              << " | Time " << std::setw(12) << t_current_
              << " | dt " << std::setw(12) << dt_
              << " | rho:[" << stats_.min_rho << "," << stats_.max_rho << "]"
              << " | p:[" << stats_.min_p << "," << stats_.max_p << "]\n";
}

void FVMSolver3D::reconstruct_1d(
    const StateField3D& state,
    int direction,
    int i, int j, int k,
    Eigen::VectorXd& U_L,
    Eigen::VectorXd& U_R
) {
    // Use new reconstruction API that operates directly on Field3D
    // Reconstructs left and right states at interface (i+1/2, j, k) for direction=0
    reconstruction_->reconstruct(state, i, j, k, direction, U_L, U_R);
}

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
