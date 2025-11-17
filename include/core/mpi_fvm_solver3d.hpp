#pragma once

#include "grid3d.hpp"
#include "field3d.hpp"
#include "physics/physics_base.hpp"
#include "spatial/flux_calculation/flux_calculator_base.hpp"
#include "temporal/time_integrator.hpp"
#include "spatial/reconstruction/reconstruction_base.hpp"
#include "boundary/boundary_condition.hpp"
#include "parallel/mpi_utils.hpp"
#include "parallel/mpi_domain_decomposer.hpp"
#include "parallel/mpi_halo_exchange.hpp"
#include "parallel/mpi_global_reduction.hpp"
#include <memory>
#include <string>
#include <functional>

namespace fvm3d::core {

/**
 * Configuration structure for MPI-parallel FVM solver.
 */
struct MPIFVMSolverConfig {
    // Global grid parameters
    double xmin, ymin, zmin;    // Global domain minimum coordinates
    double Lx, Ly, Lz;          // Global domain dimensions
    int nx, ny, nz;             // Global number of cells
    int nghost;                 // Number of ghost cells (typically 2)

    // MPI decomposition (optional - auto-determined if not specified)
    int px, py, pz;             // Process grid dimensions (0 = auto)

    // Physics type
    std::string physics_type;   // "euler", "mhd", "mhd_advanced"
    int num_vars;               // Number of variables (5 for Euler, 8/9 for MHD)

    // MHD-specific parameters (only used when physics_type = "mhd" or "mhd_advanced")
    double mhd_resistivity = 0.0;            // Basic resistivity for simple MHD
    double mhd_eta0 = 1e-3;                  // Background resistivity for advanced MHD
    double mhd_eta1 = 1.67e-2;               // Enhanced resistivity for advanced MHD
    double mhd_localization_scale = 1.0;     // Resistivity localization scale
    double mhd_glm_ch = 0.2;                 // GLM divergence wave speed
    double mhd_glm_cr = 0.2;                 // GLM parabolic dissipation ratio

    // Numerical schemes
    std::string flux_calculator;     // "laxfriedrichs", "hll", "hllc", "hlld"
    std::string reconstruction;     // "constant", "muscl"
    std::string reconstruction_limiter;  // "minmod", "van_leer", "superbee"
    std::string time_integrator;    // "euler", "rk2", "rk3"
    std::string boundary_condition; // "periodic", "reflective", "transmissive"

    // Boundary condition flags
    bool bc_x, bc_y, bc_z;

    // Time stepping
    double cfl;        // CFL number (typically 0.3-0.4)
    double t_final;    // Final simulation time
    int num_steps;     // Maximum number of steps
    int output_interval; // Output every N steps
    int checkpoint_interval; // Checkpoint every N steps (0 = disabled)

    // Optional
    int verbose;       // Verbosity level (0=silent, 1=normal, 2=debug)
};

/**
 * MPI-parallel FVM Solver for 3D compressible flow and MHD.
 *
 * Features:
 * - Domain decomposition across MPI ranks
 * - Non-blocking halo exchange for ghost cells
 * - Global reductions for CFL, statistics, convergence
 * - Parallel checkpoint/restart (planned)
 * - Support for Euler and MHD equations
 *
 * Architecture:
 * - Each rank owns a local subdomain with ghost cells
 * - Parallel RHS computation with halo exchange
 * - Global CFL constraint via MPI_Allreduce
 * - Synchronized time stepping across all ranks
 */
class MPIFVMSolver3D {
public:
    /**
     * Constructor: Initialize MPI-parallel solver.
     * Must be called after MPI_Init() and within an MPIContext.
     */
    MPIFVMSolver3D(const MPIFVMSolverConfig& config);

    /**
     * Destructor.
     */
    ~MPIFVMSolver3D() = default;

    /**
     * Initialize with initial condition function.
     * @param init_func: function that sets U at each global grid point (x, y, z)
     *
     * Note: init_func is called for each local cell using global coordinates.
     */
    using InitFunction = std::function<void(double, double, double, Eigen::VectorXd&)>;
    void initialize(const InitFunction& init_func);

    /**
     * Run the main simulation loop.
     * All ranks execute synchronously.
     */
    void run();

    /**
     * Perform a single time step.
     * Includes halo exchange and global CFL reduction.
     */
    void step();

    /**
     * Get current simulation time (same on all ranks).
     */
    double time() const { return t_current_; }

    /**
     * Get current step number (same on all ranks).
     */
    int step_count() const { return step_count_; }

    /**
     * Get local state field (this rank's subdomain only).
     */
    StateField3D& state() { return *state_; }
    const StateField3D& state() const { return *state_; }

    /**
     * Get local grid information.
     */
    const Grid3D& grid() const { return *local_grid_; }

    /**
     * Get domain decomposer (for global/local index mapping).
     */
    const parallel::MPIDomainDecomposer& decomposer() const { return *decomposer_; }

    /**
     * Get MPI rank.
     */
    int rank() const { return rank_; }

    /**
     * Get MPI size.
     */
    int size() const { return size_; }

    /**
     * Save parallel checkpoint to HDF5 file (planned).
     * Each rank writes its subdomain.
     */
    void save_checkpoint(
        const std::string& filename,
        const std::string& description = ""
    );

    /**
     * Load parallel checkpoint from HDF5 file (planned).
     */
    bool load_checkpoint(const std::string& filename);

private:
    // Configuration
    MPIFVMSolverConfig config_;
    int rank_, size_;

    // MPI components
    std::unique_ptr<parallel::MPIDomainDecomposer> decomposer_;
    std::unique_ptr<parallel::MPIHaloExchange> halo_exchange_;
    std::unique_ptr<parallel::MPIGlobalReduction> global_reduction_;

    // Local grid and state (using unique_ptr for deferred initialization)
    std::unique_ptr<Grid3D> local_grid_;               // This rank's subdomain
    std::unique_ptr<StateField3D> state_;              // Local conservative variables
    std::unique_ptr<StateField3D> rhs_;                // Local RHS
    std::unique_ptr<StateField3D> u_temp_;             // Temporary storage
    std::unique_ptr<StateField3D> flux_x_, flux_y_, flux_z_;  // Local fluxes

    // Physics and numerical schemes (polymorphic)
    std::shared_ptr<physics::PhysicsBase> physics_;  // Physics equation system
    std::unique_ptr<spatial::FluxCalculator> flux_calculator_;
    std::unique_ptr<spatial::ReconstructionMethod> reconstruction_;
    std::unique_ptr<temporal::TimeIntegrator> time_integrator_;
    std::unique_ptr<boundary::BoundaryCondition> boundary_condition_;

    // Time stepping
    double t_current_;
    int step_count_;
    double dt_;

    // Monitoring
    struct Stats {
        double min_rho, max_rho;
        double min_p, max_p;
        double min_speed, max_speed;
        double min_B, max_B;  // For MHD
    } stats_;

    /**
     * Compute the right-hand side (flux divergence + source terms).
     * Includes halo exchange before flux computation.
     */
    void compute_rhs();

    /**
     * Compute fluxes in a given direction on local subdomain.
     * @param direction: 0=X, 1=Y, 2=Z
     */
    void compute_fluxes(int direction, StateField3D& flux_out);

    /**
     * Compute flux divergence on local subdomain.
     */
    void add_flux_divergence(const StateField3D& flux, int direction);

    /**
     * Apply boundary conditions to ghost cells.
     * For domain boundaries and MPI boundaries (after halo exchange).
     */
    void apply_boundary_conditions();

    /**
     * Compute adaptive time step using CFL condition.
     * Performs global reduction to get minimum dt across all ranks.
     */
    double compute_dt();

    /**
     * Compute local statistics and perform global reduction.
     */
    void compute_statistics();

    /**
     * Print progress information (rank 0 only).
     */
    void print_progress();

    /**
     * Reconstruct interface values for a 1D slice (local subdomain).
     */
    void reconstruct_1d(
        const StateField3D& state,
        int direction,
        int i, int j, int k,
        Eigen::VectorXd& U_L,
        Eigen::VectorXd& U_R
    );

    /**
     * Compute flux for a single interface using physics-dependent method.
     * Handles Euler and MHD cases.
     */
    Eigen::VectorXd compute_interface_flux(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction,
        double x, double y, double z
    );

    /**
     * Helper: Initialize physics object based on config.
     */
    void initialize_physics();

    /**
     * Helper: Map global to local indices.
     */
    void global_to_local(int gi, int gj, int gk, int& li, int& lj, int& lk) const;

    /**
     * Helper: Map local to global indices.
     */
    void local_to_global(int li, int lj, int lk, int& gi, int& gj, int& gk) const;
};

} // namespace fvm3d::core
