#include "core/fvm_solver3d.hpp"
#include "physics/resistive_mhd3d_advanced.hpp"
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

    // ================================================================
    // STRANG SPLITTING: D^(dt/2) ∘ L^(dt) ∘ D^(dt/2)
    // ================================================================
    // Reference: OpenMHD glm_ss2.f90, Dedner et al. (2002)
    //
    // This implements operator splitting for GLM divergence cleaning:
    // - D: Damping operator (exponential decay of ψ)
    // - L: Hyperbolic operator (flux evolution via RK)
    //
    // The splitting is second-order accurate and unconditionally stable.
    // ================================================================

    // ========== GLM DAMPING: First half-step ==========
    // Apply ψ *= exp(-0.5 × dt × ch/cr) to all cells
    auto mhd = std::dynamic_pointer_cast<physics::AdvancedResistiveMHD3D>(physics_);
    if (mhd) {
        const int i_begin = grid_.i_begin();
        const int i_end = grid_.i_end();
        const int j_begin = grid_.j_begin();
        const int j_end = grid_.j_end();
        const int k_begin = grid_.k_begin();
        const int k_end = grid_.k_end();
        const int nvars = state_.nvars();

        #pragma omp parallel
        {
            Eigen::VectorXd U(nvars);

            #pragma omp for collapse(3)
            for (int i = i_begin; i < i_end; i++) {
                for (int j = j_begin; j < j_end; j++) {
                    for (int k = k_begin; k < k_end; k++) {
                        // Extract state vector at (i,j,k)
                        for (int v = 0; v < nvars; v++) {
                            U(v) = state_(v, i, j, k);
                        }

                        // Apply GLM damping
                        mhd->apply_glm_damping(U, dt_, 0.5);

                        // Write back
                        for (int v = 0; v < nvars; v++) {
                            state_(v, i, j, k) = U(v);
                        }
                    }
                }
            }
        }
    }

    // ========== MHD EVOLUTION: Full time step (RK2/RK3) ==========
    // Time integration using RHS callback
    auto rhs_func = [this](const StateField3D& U, StateField3D& dUdt) {
        this->compute_rhs();
        dUdt.assign(this->rhs_);
    };

    time_integrator_->step(state_, dt_, rhs_func);

    // ========== GLM DAMPING: Second half-step ==========
    // Apply ψ *= exp(-0.5 × dt × ch/cr) to all cells
    if (mhd) {
        const int i_begin = grid_.i_begin();
        const int i_end = grid_.i_end();
        const int j_begin = grid_.j_begin();
        const int j_end = grid_.j_end();
        const int k_begin = grid_.k_begin();
        const int k_end = grid_.k_end();
        const int nvars = state_.nvars();

        #pragma omp parallel
        {
            Eigen::VectorXd U(nvars);

            #pragma omp for collapse(3)
            for (int i = i_begin; i < i_end; i++) {
                for (int j = j_begin; j < j_end; j++) {
                    for (int k = k_begin; k < k_end; k++) {
                        // Extract state vector at (i,j,k)
                        for (int v = 0; v < nvars; v++) {
                            U(v) = state_(v, i, j, k);
                        }

                        // Apply GLM damping
                        mhd->apply_glm_damping(U, dt_, 0.5);

                        // Write back
                        for (int v = 0; v < nvars; v++) {
                            state_(v, i, j, k) = U(v);
                        }
                    }
                }
            }
        }
    }

    // ========== POSITIVITY LIMITER ==========
    // Ensure ρ > 0 and p > 0 after time evolution
    apply_positivity_limiter();

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
    // Initialize statistics
    double min_rho = 1e10;
    double max_rho = -1e10;
    double min_p = 1e10;
    double max_p = -1e10;
    double min_speed = 1e10;
    double max_speed = -1e10;
    double kinetic_energy = 0.0;
    double magnetic_energy = 0.0;
    double internal_energy = 0.0;

    const int i_begin = grid_.i_begin();
    const int i_end = grid_.i_end();
    const int j_begin = grid_.j_begin();
    const int j_end = grid_.j_end();
    const int k_begin = grid_.k_begin();
    const int k_end = grid_.k_end();
    const int nvars = config_.num_vars;

    const auto& geom = grid_.geometry();
    const double dV = geom.dx * geom.dy * geom.dz;  // Cell volume

    // OpenMP parallel reduction for statistics
    #pragma omp parallel reduction(min:min_rho,min_p,min_speed) \
                         reduction(max:max_rho,max_p,max_speed) \
                         reduction(+:kinetic_energy,magnetic_energy,internal_energy)
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

                    // Update local min/max
                    min_rho = std::min(min_rho, rho);
                    max_rho = std::max(max_rho, rho);
                    min_p = std::min(min_p, p);
                    max_p = std::max(max_p, p);

                    const double speed = std::sqrt(u*u + v*v + w*w);
                    min_speed = std::min(min_speed, speed);
                    max_speed = std::max(max_speed, speed);

                    // Compute energies for MHD
                    if (nvars == 9) {  // MHD with GLM
                        // Kinetic energy: 0.5 * rho * v^2
                        kinetic_energy += 0.5 * rho * (u*u + v*v + w*w) * dV;

                        // Magnetic energy: 0.5 * B^2
                        const double Bx = V(5);
                        const double By = V(6);
                        const double Bz = V(7);
                        magnetic_energy += 0.5 * (Bx*Bx + By*By + Bz*Bz) * dV;

                        // Internal energy: p / (gamma - 1)
                        const double gamma = 5.0/3.0;  // Adiabatic index
                        internal_energy += p / (gamma - 1.0) * dV;
                    }
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
    stats_.kinetic_energy = kinetic_energy;
    stats_.magnetic_energy = magnetic_energy;
    stats_.internal_energy = internal_energy;
    stats_.total_energy = kinetic_energy + magnetic_energy + internal_energy;
    stats_.max_div_B = compute_max_div_B();
}

double FVMSolver3D::compute_max_div_B() const {
    // Only compute for MHD simulations
    if (config_.num_vars != 9) {
        return 0.0;
    }

    const int i_begin = grid_.i_begin();
    const int i_end = grid_.i_end();
    const int j_begin = grid_.j_begin();
    const int j_end = grid_.j_end();
    const int k_begin = grid_.k_begin();
    const int k_end = grid_.k_end();

    const auto& geom = grid_.geometry();
    const double inv_dx = 1.0 / geom.dx;
    const double inv_dy = 1.0 / geom.dy;
    const double inv_dz = 1.0 / geom.dz;

    double max_div_B = 0.0;

    // Compute div(B) = dBx/dx + dBy/dy + dBz/dz using centered differences
    #pragma omp parallel reduction(max:max_div_B)
    {
        #pragma omp for collapse(2)
        for (int i = i_begin; i < i_end; i++) {
            for (int j = j_begin; j < j_end; j++) {
                for (int k = k_begin; k < k_end; k++) {
                    // Central differences for interior points
                    if (i > i_begin && i < i_end-1 &&
                        j > j_begin && j < j_end-1 &&
                        k > k_begin && k < k_end-1) {

                        const double dBx_dx = (state_(5, i+1, j, k) - state_(5, i-1, j, k)) * 0.5 * inv_dx;
                        const double dBy_dy = (state_(6, i, j+1, k) - state_(6, i, j-1, k)) * 0.5 * inv_dy;
                        const double dBz_dz = (state_(7, i, j, k+1) - state_(7, i, j, k-1)) * 0.5 * inv_dz;

                        const double div_B = std::abs(dBx_dx + dBy_dy + dBz_dz);
                        max_div_B = std::max(max_div_B, div_B);
                    }
                }
            }
        }
    }

    return max_div_B;
}

void FVMSolver3D::apply_positivity_limiter() {
    const int i_begin = grid_.i_begin();
    const int i_end = grid_.i_end();
    const int j_begin = grid_.j_begin();
    const int j_end = grid_.j_end();
    const int k_begin = grid_.k_begin();
    const int k_end = grid_.k_end();
    const int nvars = config_.num_vars;

    // Density and pressure floors - use small but reasonable values
    // Set based on expected minimum values in Harris sheet
    const double rho_floor = 0.01;   // 1% of typical density
    const double p_floor = 0.001;    // Small but physical

    int num_fixes = 0;
    double total_energy_added = 0.0;

    #pragma omp parallel reduction(+:num_fixes,total_energy_added)
    {
        Eigen::VectorXd U(nvars);
        Eigen::VectorXd V(nvars);

        #pragma omp for collapse(3)
        for (int i = i_begin; i < i_end; i++) {
            for (int j = j_begin; j < j_end; j++) {
                for (int k = k_begin; k < k_end; k++) {
                    // Load conservative variables
                    for (int v = 0; v < nvars; v++) {
                        U(v) = state_(v, i, j, k);
                    }

                    // Store original energy for tracking
                    double E_before = U(4);  // Total energy (conservative)

                    // Convert to primitive to check positivity
                    V = physics_->conservative_to_primitive(U);

                    bool needs_fix = false;

                    // Check density - use gentler approach
                    if (V(0) < rho_floor) {
                        // Gently push towards floor, not hard set
                        V(0) = 0.5 * (V(0) + rho_floor);
                        if (V(0) < rho_floor) V(0) = rho_floor;
                        needs_fix = true;
                    }

                    // Check pressure - critical for energy conservation
                    if (nvars >= 5 && V(4) < p_floor) {
                        // Very gentle pressure fix to minimize energy injection
                        double p_old = V(4);
                        V(4) = 0.9 * p_old + 0.1 * p_floor;
                        if (V(4) < p_floor * 0.1) V(4) = p_floor * 0.1;
                        needs_fix = true;
                    }

                    // If fixes were applied, convert back to conservative
                    if (needs_fix) {
                        U = physics_->primitive_to_conservative(V);

                        // Track energy added
                        double E_after = U(4);
                        total_energy_added += (E_after - E_before);

                        // Write back to state
                        for (int v = 0; v < nvars; v++) {
                            state_(v, i, j, k) = U(v);
                        }

                        num_fixes++;
                    }
                }
            }
        }
    }

    // Report if fixes were applied
    if (num_fixes > 0 && config_.verbose > 1) {  // Only verbose > 1
        const auto& geom = grid_.geometry();
        const double dV = geom.dx * geom.dy * geom.dz;
        std::cout << "  [Positivity] Fixed " << num_fixes << " cells, ΔE = "
                  << std::scientific << std::setprecision(2)
                  << (total_energy_added * dV) << "\n";
    }
}

void FVMSolver3D::print_progress() {
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
