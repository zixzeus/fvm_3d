#include "core/mpi_fvm_solver3d.hpp"
#include "spatial/riemann_solver_factory.hpp"
#include "spatial/reconstruction.hpp"
#include "temporal/time_integrator_factory.hpp"
#include "boundary/periodic_bc.hpp"
#include "boundary/reflective_bc.hpp"
#include "boundary/transmissive_bc.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <mpi.h>

namespace fvm3d::core {

MPIFVMSolver3D::MPIFVMSolver3D(const MPIFVMSolverConfig& config)
    : config_(config),
      t_current_(0.0),
      step_count_(0),
      dt_(0.0)
{
    // Get MPI rank and size
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    // Create domain decomposer
    decomposer_ = std::make_unique<parallel::MPIDomainDecomposer>(
        config.nx, config.ny, config.nz,
        config.px, config.py, config.pz,  // 0 = auto-determine
        MPI_COMM_WORLD
    );

    // Get local domain size
    int local_nx = decomposer_->local_nx();
    int local_ny = decomposer_->local_ny();
    int local_nz = decomposer_->local_nz();

    // Get local domain offset in global grid
    int offset_x = decomposer_->local_offset_x();
    int offset_y = decomposer_->local_offset_y();
    int offset_z = decomposer_->local_offset_z();

    // Compute local domain physical extent
    double dx_global = config.Lx / config.nx;
    double dy_global = config.Ly / config.ny;
    double dz_global = config.Lz / config.nz;

    double local_xmin = config.xmin + offset_x * dx_global;
    double local_ymin = config.ymin + offset_y * dy_global;
    double local_zmin = config.zmin + offset_z * dz_global;

    double local_Lx = local_nx * dx_global;
    double local_Ly = local_ny * dy_global;
    double local_Lz = local_nz * dz_global;

    // Create local grid
    GridGeometry3D local_geom(
        local_xmin, local_ymin, local_zmin,
        local_Lx, local_Ly, local_Lz,
        local_nx, local_ny, local_nz
    );
    local_grid_ = Grid3D(local_geom, config.nghost);

    // Allocate local state arrays
    int num_vars = config.num_vars;
    int nx_total = local_grid_.nx_total();
    int ny_total = local_grid_.ny_total();
    int nz_total = local_grid_.nz_total();

    state_ = StateField3D(num_vars, nx_total, ny_total, nz_total);
    rhs_ = StateField3D(num_vars, nx_total, ny_total, nz_total);
    u_temp_ = StateField3D(num_vars, nx_total, ny_total, nz_total);
    flux_x_ = StateField3D(num_vars, nx_total, ny_total, nz_total);
    flux_y_ = StateField3D(num_vars, nx_total, ny_total, nz_total);
    flux_z_ = StateField3D(num_vars, nx_total, ny_total, nz_total);

    // Create MPI halo exchange
    halo_exchange_ = std::make_unique<parallel::MPIHaloExchange>(
        *decomposer_,
        local_nx, local_ny, local_nz,
        config.nghost,
        num_vars
    );

    // Create global reduction handler
    global_reduction_ = std::make_unique<parallel::MPIGlobalReduction>(MPI_COMM_WORLD);

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

    // Print initialization info (rank 0 only)
    if (rank_ == 0 && config.verbose > 0) {
        std::cout << "=== MPI FVM3D Solver Initialized ===\n"
                  << "Global Grid: " << config.nx << " x " << config.ny << " x " << config.nz << "\n"
                  << "Global Domain: [" << config.xmin << "," << (config.xmin + config.Lx) << "] x "
                  << "[" << config.ymin << "," << (config.ymin + config.Ly) << "] x "
                  << "[" << config.zmin << "," << (config.zmin + config.Lz) << "]\n"
                  << "MPI Processes: " << size_ << " ("
                  << decomposer_->px() << " x "
                  << decomposer_->py() << " x "
                  << decomposer_->pz() << ")\n"
                  << "Physics: " << config.physics_type << " (" << num_vars << " variables)\n"
                  << "Riemann Solver: " << riemann_solver_->name() << "\n"
                  << "Reconstruction: " << reconstruction_->name() << "\n"
                  << "Time Integrator: " << time_integrator_->name() << "\n"
                  << "Boundary Condition: " << boundary_condition_->name() << "\n"
                  << "CFL: " << config.cfl << "\n"
                  << "====================================\n";
    }

    if (config.verbose > 1) {
        std::cout << "[Rank " << rank_ << "] Local grid: "
                  << local_nx << " x " << local_ny << " x " << local_nz
                  << " (offset: " << offset_x << ", " << offset_y << ", " << offset_z << ")\n";
    }
}

void MPIFVMSolver3D::initialize(const InitFunction& init_func) {
    int nx_total = local_grid_.nx_total();
    int ny_total = local_grid_.ny_total();
    int nz_total = local_grid_.nz_total();
    int num_vars = config_.num_vars;

    // Initialize local subdomain with global coordinates
    for (int i = 0; i < nx_total; i++) {
        for (int j = 0; j < ny_total; j++) {
            for (int k = 0; k < nz_total; k++) {
                double x = local_grid_.cell_center_x(i);
                double y = local_grid_.cell_center_y(j);
                double z = local_grid_.cell_center_z(k);

                Eigen::VectorXd U(num_vars);
                init_func(x, y, z, U);

                for (int v = 0; v < num_vars; v++) {
                    state_(v, i, j, k) = U(v);
                }
            }
        }
    }

    // Apply boundary conditions and halo exchange
    apply_boundary_conditions();
    halo_exchange_->exchange(state_);

    if (rank_ == 0 && config_.verbose > 0) {
        std::cout << "Initial conditions set.\n";
    }
}

void MPIFVMSolver3D::run() {
    if (rank_ == 0 && config_.verbose > 0) {
        std::cout << "\n=== Starting Simulation ===\n";
    }

    while (t_current_ < config_.t_final && step_count_ < config_.num_steps) {
        step();

        if (config_.output_interval > 0 && step_count_ % config_.output_interval == 0) {
            compute_statistics();
            print_progress();
        }

        if (config_.checkpoint_interval > 0 && step_count_ % config_.checkpoint_interval == 0) {
            // TODO: Implement parallel checkpoint
            if (rank_ == 0 && config_.verbose > 0) {
                std::cout << "[Checkpoint at step " << step_count_ << "]\n";
            }
        }
    }

    if (rank_ == 0 && config_.verbose > 0) {
        std::cout << "\n=== Simulation Complete ===\n"
                  << "Final time: " << t_current_ << "\n"
                  << "Total steps: " << step_count_ << "\n";
    }
}

void MPIFVMSolver3D::step() {
    // Compute adaptive time step (global reduction)
    dt_ = compute_dt();

    // Ensure we don't overshoot t_final
    if (t_current_ + dt_ > config_.t_final) {
        dt_ = config_.t_final - t_current_;
    }

    // Apply boundary conditions before RHS computation
    apply_boundary_conditions();

    // Perform time integration step
    // The RHS function is a lambda that computes RHS with halo exchange
    auto rhs_function = [this](const StateField3D& U, StateField3D& dUdt) {
        // Update state for RHS computation
        state_.assign(U);

        // Exchange ghost cells
        halo_exchange_->exchange(state_);

        // Compute RHS
        compute_rhs();

        // Copy RHS to output
        dUdt.assign(rhs_);
    };

    // Perform integration step
    time_integrator_->step(state_, rhs_, dt_, rhs_function);

    // Update time and step count
    t_current_ += dt_;
    step_count_++;
}

void MPIFVMSolver3D::compute_rhs() {
    // Zero out RHS
    rhs_.fill(0.0);

    // Compute fluxes in each direction
    compute_fluxes(0, flux_x_);  // X-direction
    compute_fluxes(1, flux_y_);  // Y-direction
    compute_fluxes(2, flux_z_);  // Z-direction

    // Add flux divergence: dU/dt = -(dF/dx + dG/dy + dH/dz)
    add_flux_divergence(flux_x_, 0);
    add_flux_divergence(flux_y_, 1);
    add_flux_divergence(flux_z_, 2);

    // TODO: Add source terms for MHD (resistivity, GLM)
}

void MPIFVMSolver3D::compute_fluxes(int direction, StateField3D& flux_out) {
    int nghost = local_grid_.nghost();
    int nx = local_grid_.nx();
    int ny = local_grid_.ny();
    int nz = local_grid_.nz();
    int num_vars = config_.num_vars;

    // Interior cell range
    int i_start = nghost;
    int j_start = nghost;
    int k_start = nghost;
    int i_end = nghost + nx;
    int j_end = nghost + ny;
    int k_end = nghost + nz;

    // Compute fluxes at all interfaces in the specified direction
    for (int i = i_start; i <= i_end; i++) {
        for (int j = j_start; j <= j_end; j++) {
            for (int k = k_start; k <= k_end; k++) {
                // Reconstruct left and right states at interface
                Eigen::VectorXd U_L(num_vars), U_R(num_vars);
                reconstruct_1d(state_, direction, i, j, k, U_L, U_R);

                // Get interface position
                double x = (direction == 0) ? local_grid_.cell_center_x(i) : local_grid_.cell_center_x(i);
                double y = (direction == 1) ? local_grid_.cell_center_y(j) : local_grid_.cell_center_y(j);
                double z = (direction == 2) ? local_grid_.cell_center_z(k) : local_grid_.cell_center_z(k);

                // Solve Riemann problem
                Eigen::VectorXd F = compute_interface_flux(U_L, U_R, direction, x, y, z);

                // Store flux
                for (int v = 0; v < num_vars; v++) {
                    flux_out(v, i, j, k) = F(v);
                }
            }
        }
    }
}

void MPIFVMSolver3D::add_flux_divergence(const StateField3D& flux, int direction) {
    int nghost = local_grid_.nghost();
    int nx = local_grid_.nx();
    int ny = local_grid_.ny();
    int nz = local_grid_.nz();
    int num_vars = config_.num_vars;

    double dx = local_grid_.geometry().dx;
    double dy = local_grid_.geometry().dy;
    double dz = local_grid_.geometry().dz;
    double inv_d = (direction == 0) ? (1.0 / dx) : ((direction == 1) ? (1.0 / dy) : (1.0 / dz));

    // Update RHS in interior cells
    for (int i = nghost; i < nghost + nx; i++) {
        for (int j = nghost; j < nghost + ny; j++) {
            for (int k = nghost; k < nghost + nz; k++) {
                // Compute flux difference
                int ip = (direction == 0) ? (i + 1) : i;
                int jp = (direction == 1) ? (j + 1) : j;
                int kp = (direction == 2) ? (k + 1) : k;

                for (int v = 0; v < num_vars; v++) {
                    double flux_diff = flux(v, ip, jp, kp) - flux(v, i, j, k);
                    rhs_(v, i, j, k) -= flux_diff * inv_d;
                }
            }
        }
    }
}

void MPIFVMSolver3D::apply_boundary_conditions() {
    // Apply physical boundary conditions only at domain boundaries
    // MPI boundaries are handled by halo exchange

    // Check if this rank is at domain boundary
    bool at_xmin = decomposer_->at_boundary_xmin();
    bool at_xmax = decomposer_->at_boundary_xmax();
    bool at_ymin = decomposer_->at_boundary_ymin();
    bool at_ymax = decomposer_->at_boundary_ymax();
    bool at_zmin = decomposer_->at_boundary_zmin();
    bool at_zmax = decomposer_->at_boundary_zmax();

    // Apply boundary conditions using the boundary condition object
    // The BC object handles interior boundaries; we just call it
    boundary_condition_->apply(state_, local_grid_);
}

double MPIFVMSolver3D::compute_dt() {
    // Compute local maximum wave speed
    double max_speed = 0.0;
    int nghost = local_grid_.nghost();
    int nx = local_grid_.nx();
    int ny = local_grid_.ny();
    int nz = local_grid_.nz();

    // Scan interior cells
    for (int i = nghost; i < nghost + nx; i++) {
        for (int j = nghost; j < nghost + ny; j++) {
            for (int k = nghost; k < nghost + nz; k++) {
                Eigen::VectorXd U(config_.num_vars);
                for (int v = 0; v < config_.num_vars; v++) {
                    U(v) = state_(v, i, j, k);
                }

                // Compute max wave speed in each direction
                // For now, use simple estimate (will need physics-specific method)
                double speed_x = 0.0, speed_y = 0.0, speed_z = 0.0;

                // TODO: Use physics-specific wave speed calculation
                // For Euler: speed = |u| + a
                // For MHD: speed = |u| + c_f (fast magnetosonic speed)

                double speed = std::max({speed_x, speed_y, speed_z});
                max_speed = std::max(max_speed, speed);
            }
        }
    }

    // Global reduction to find maximum wave speed across all ranks
    double global_max_speed = global_reduction_->max(max_speed);

    // Compute CFL-limited time step
    double dx = local_grid_.geometry().dx;
    double dy = local_grid_.geometry().dy;
    double dz = local_grid_.geometry().dz;
    double min_d = std::min({dx, dy, dz});

    double dt_cfl = (global_max_speed > 1e-12) ? (config_.cfl * min_d / global_max_speed) : 1e-3;

    return dt_cfl;
}

void MPIFVMSolver3D::compute_statistics() {
    // Compute local statistics
    double local_min_rho = 1e10, local_max_rho = -1e10;
    int nghost = local_grid_.nghost();
    int nx = local_grid_.nx();
    int ny = local_grid_.ny();
    int nz = local_grid_.nz();

    for (int i = nghost; i < nghost + nx; i++) {
        for (int j = nghost; j < nghost + ny; j++) {
            for (int k = nghost; k < nghost + nz; k++) {
                double rho = state_(0, i, j, k);
                local_min_rho = std::min(local_min_rho, rho);
                local_max_rho = std::max(local_max_rho, rho);
            }
        }
    }

    // Global reduction
    stats_.min_rho = global_reduction_->min(local_min_rho);
    stats_.max_rho = global_reduction_->max(local_max_rho);

    // TODO: Add more statistics (pressure, velocity, magnetic field)
}

void MPIFVMSolver3D::print_progress() {
    if (rank_ == 0) {
        std::cout << "Step " << std::setw(6) << step_count_
                  << "  t=" << std::setw(10) << std::scientific << std::setprecision(4) << t_current_
                  << "  dt=" << std::setw(10) << dt_
                  << "  rho=[" << stats_.min_rho << ", " << stats_.max_rho << "]\n";
    }
}

void MPIFVMSolver3D::reconstruct_1d(
    const StateField3D& state,
    int direction,
    int i, int j, int k,
    Eigen::VectorXd& U_L,
    Eigen::VectorXd& U_R)
{
    int num_vars = config_.num_vars;

    // Get 5-point stencil for reconstruction
    int im2 = (direction == 0) ? (i - 2) : i;
    int jm2 = (direction == 1) ? (j - 2) : j;
    int km2 = (direction == 2) ? (k - 2) : k;

    int im1 = (direction == 0) ? (i - 1) : i;
    int jm1 = (direction == 1) ? (j - 1) : j;
    int km1 = (direction == 2) ? (k - 1) : k;

    int ip1 = (direction == 0) ? (i + 1) : i;
    int jp1 = (direction == 1) ? (j + 1) : j;
    int kp1 = (direction == 2) ? (k + 1) : k;

    Eigen::VectorXd Um2(num_vars), Um1(num_vars), U0(num_vars), Up1(num_vars);

    for (int v = 0; v < num_vars; v++) {
        Um2(v) = state(v, im2, jm2, km2);
        Um1(v) = state(v, im1, jm1, km1);
        U0(v) = state(v, i, j, k);
        Up1(v) = state(v, ip1, jp1, kp1);
    }

    // Use reconstruction scheme
    reconstruction_->reconstruct(Um2, Um1, U0, Up1, U_L, U_R);
}

Eigen::VectorXd MPIFVMSolver3D::compute_interface_flux(
    const Eigen::VectorXd& U_L,
    const Eigen::VectorXd& U_R,
    int direction,
    double x, double y, double z)
{
    // Use Riemann solver to compute flux
    return riemann_solver_->solve(U_L, U_R, direction);
}

void MPIFVMSolver3D::save_checkpoint(const std::string& filename, const std::string& description) {
    // TODO: Implement parallel HDF5 checkpoint
    if (rank_ == 0) {
        std::cout << "Parallel checkpoint not yet implemented.\n";
    }
}

bool MPIFVMSolver3D::load_checkpoint(const std::string& filename) {
    // TODO: Implement parallel HDF5 checkpoint
    if (rank_ == 0) {
        std::cout << "Parallel checkpoint not yet implemented.\n";
    }
    return false;
}

} // namespace fvm3d::core
