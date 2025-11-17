#pragma once

#include "grid3d.hpp"
#include "field3d.hpp"
#include "physics/physics_base.hpp"
#include "spatial/flux_calculation/flux_calculator_base.hpp"
#include "temporal/time_integrator.hpp"
#include "spatial/reconstruction/reconstruction_base.hpp"
#include "boundary/boundary_condition.hpp"
#include <memory>
#include <string>

namespace fvm3d::core {

/**
 * Base class for FVM solvers (serial and parallel).
 *
 * Provides common functionality for:
 * - Physics initialization
 * - Numerical scheme setup (flux calculator, reconstruction, time integrator)
 * - Boundary condition setup
 * - Common utility methods
 *
 * This class reduces code duplication between FVMSolver3D and MPIFVMSolver3D.
 */
class FVMSolverBase {
protected:
    // Physics and numerical schemes (shared by all solvers)
    std::shared_ptr<physics::PhysicsBase> physics_;
    std::unique_ptr<spatial::FluxCalculator> flux_calculator_;
    std::unique_ptr<spatial::ReconstructionMethod> reconstruction_;
    std::unique_ptr<temporal::TimeIntegrator> time_integrator_;
    std::unique_ptr<boundary::BoundaryCondition> boundary_condition_;

    /**
     * Initialize physics object based on configuration parameters.
     *
     * @param physics_type Type of physics: "euler", "mhd", "mhd_advanced"
     * @param mhd_resistivity Basic resistivity coefficient for simple MHD
     * @param mhd_eta0 Background resistivity for advanced MHD
     * @param mhd_eta1 Enhanced resistivity for advanced MHD
     * @param mhd_localization_scale Resistivity localization scale
     * @param mhd_glm_ch GLM divergence wave speed
     * @param mhd_glm_cr GLM parabolic dissipation ratio
     */
    void initialize_physics(
        const std::string& physics_type,
        double mhd_resistivity,
        double mhd_eta0,
        double mhd_eta1,
        double mhd_localization_scale,
        double mhd_glm_ch,
        double mhd_glm_cr
    );

    /**
     * Initialize flux calculator based on configuration.
     *
     * @param flux_calculator_name Name of flux calculator (e.g., "hll", "hllc", "hlld")
     */
    void initialize_flux_calculator(const std::string& flux_calculator_name);

    /**
     * Initialize reconstruction scheme.
     *
     * @param reconstruction_name Name of reconstruction method (e.g., "constant", "muscl")
     * @param num_vars Number of variables
     * @param limiter_name Name of limiter (e.g., "minmod", "van_leer", "superbee")
     */
    void initialize_reconstruction(
        const std::string& reconstruction_name,
        int num_vars,
        const std::string& limiter_name
    );

    /**
     * Initialize time integrator.
     *
     * @param integrator_name Name of time integrator (e.g., "euler", "rk2", "rk3")
     */
    void initialize_time_integrator(const std::string& integrator_name);

    /**
     * Initialize boundary conditions.
     *
     * @param bc_type Type of boundary condition (e.g., "periodic", "reflective", "transmissive")
     * @param bc_x Apply in X direction
     * @param bc_y Apply in Y direction
     * @param bc_z Apply in Z direction
     */
    void initialize_boundary_conditions(
        const std::string& bc_type,
        bool bc_x, bool bc_y, bool bc_z
    );

    /**
     * Reconstruct interface values for a 1D slice.
     * This method is identical in both serial and parallel solvers.
     *
     * @param state State field
     * @param direction Direction (0=X, 1=Y, 2=Z)
     * @param i Grid index i
     * @param j Grid index j
     * @param k Grid index k
     * @param U_L Left state (output)
     * @param U_R Right state (output)
     */
    void reconstruct_1d(
        const StateField3D& state,
        int direction,
        int i, int j, int k,
        Eigen::VectorXd& U_L,
        Eigen::VectorXd& U_R
    ) const;

public:
    /**
     * Virtual destructor for proper cleanup in derived classes.
     */
    virtual ~FVMSolverBase() = default;

    /**
     * Get physics object.
     */
    const std::shared_ptr<physics::PhysicsBase>& physics() const {
        return physics_;
    }
};

} // namespace fvm3d::core
