#pragma once

#include "grid3d.hpp"
#include "field3d.hpp"
#include "fvm_solver_base.hpp"
#include <memory>
#include <string>
#include <functional>

namespace fvm3d::core {

/**
 * Configuration structure for FVM solver.
 */
struct FVMSolverConfig {
    // Grid parameters
    double xmin, ymin, zmin;    // Domain minimum coordinates
    double Lx, Ly, Lz;          // Domain dimensions
    int nx, ny, nz;             // Number of cells
    int nghost;                 // Number of ghost cells

    // Physics type
    std::string physics_type;       // "euler", "mhd", "mhd_advanced"
    int num_vars;                   // Number of variables (5 for Euler, 8/9 for MHD)

    // MHD-specific parameters (only used when physics_type = "mhd" or "mhd_advanced")
    double mhd_resistivity = 0.0;            // Basic resistivity for simple MHD
    double mhd_eta0 = 1e-3;                  // Background resistivity for advanced MHD
    double mhd_eta1 = 1.67e-2;               // Enhanced resistivity for advanced MHD
    double mhd_localization_scale = 1.0;     // Resistivity localization scale
    double mhd_glm_ch = 0.2;                 // GLM divergence wave speed
    double mhd_glm_cr = 0.2;                 // GLM parabolic dissipation ratio

    // Numerical schemes
    std::string flux_calculator;    // "laxfriedrichs", "hll", "hllc", "hlld"
    std::string reconstruction;     // "constant", "muscl"
    std::string reconstruction_limiter;  // "minmod", "van_leer", "superbee"
    std::string time_integrator;    // "euler", "rk2", "rk3"
    std::string boundary_condition; // "periodic", "reflective", "transmissive"

    // Boundary condition flags
    bool bc_x, bc_y, bc_z;

    // Time stepping
    double cfl;        // CFL number (typically 0.4-0.5)
    double t_final;    // Final simulation time
    int num_steps;     // Maximum number of steps
    int output_interval; // Output every N steps

    // Optional
    int verbose;       // Verbosity level (0=silent, 1=normal, 2=debug)
};

/**
 * Main FVM Solver for 3D compressible flow.
 *
 * Orchestrates:
 * - Spatial discretization (reconstruction + flux calculator)
 * - Time integration (explicit Runge-Kutta)
 * - Boundary condition application
 * - Main simulation loop
 * - Output and monitoring
 */
class FVMSolver3D : public FVMSolverBase {
public:
    /**
     * Constructor: Initialize solver with configuration.
     */
    FVMSolver3D(const FVMSolverConfig& config);

    /**
     * Destructor.
     */
    ~FVMSolver3D() = default;

    /**
     * Initialize with initial condition function.
     * @param init_func: function that sets U at each grid point
     */
    using InitFunction = std::function<void(double, double, double, Eigen::VectorXd&)>;
    void initialize(const InitFunction& init_func);

    /**
     * Run the main simulation loop.
     */
    void run();

    /**
     * Perform a single time step.
     */
    void step();

    /**
     * Get current simulation time.
     */
    double time() const { return t_current_; }

    /**
     * Get current step number.
     */
    int step_count() const { return step_count_; }

    /**
     * Get current state.
     */
    StateField3D& state() { return state_; }
    const StateField3D& state() const { return state_; }

    /**
     * Get grid information.
     */
    const Grid3D& grid() const { return grid_; }

    /**
     * Save checkpoint to HDF5 file.
     * @param filename: path to output HDF5 file
     * @param description: optional description string
     */
    void save_checkpoint(
        const std::string& filename,
        const std::string& description = ""
    );

    /**
     * Load checkpoint from HDF5 file.
     * @param filename: path to input HDF5 file
     * @return true if load successful, false otherwise
     */
    bool load_checkpoint(const std::string& filename);

private:
    // Configuration and grid
    FVMSolverConfig config_;
    Grid3D grid_;

    // State and scratch arrays
    StateField3D state_;           // Conservative variables: U = [rho, rho_u, rho_v, rho_w, E]
    StateField3D rhs_;             // RHS of dU/dt = RHS(U)
    StateField3D u_temp_;          // Temporary storage for interface reconstruction
    StateField3D flux_x_, flux_y_, flux_z_;  // Fluxes in each direction

    // Note: physics_, flux_calculator_, reconstruction_, time_integrator_,
    // and boundary_condition_ are now inherited from FVMSolverBase

    // Time stepping
    double t_current_;
    int step_count_;
    double dt_;

    // Monitoring
    struct Stats {
        double min_rho, max_rho;
        double min_p, max_p;
        double min_speed, max_speed;
    } stats_;

    /**
     * Compute the right-hand side (flux divergence + source terms).
     * dU/dt = -dF/dx - dG/dy - dH/dz + S
     */
    void compute_rhs();

    /**
     * Compute fluxes in a given direction.
     * @param direction: 0=X, 1=Y, 2=Z
     */
    void compute_fluxes(int direction, StateField3D& flux_out);

    /**
     * Compute flux divergence: -dF/dx (and similar for Y, Z).
     */
    void add_flux_divergence(const StateField3D& flux, int direction);

    /**
     * Apply boundary conditions to all ghost cells.
     */
    void apply_boundary_conditions();

    /**
     * Compute adaptive time step using CFL condition.
     */
    double compute_dt();

    /**
     * Compute statistics (min/max values).
     */
    void compute_statistics();

    /**
     * Print progress information.
     */
    void print_progress();

    // Note: reconstruct_1d() is now inherited from FVMSolverBase
};

} // namespace fvm3d::core
