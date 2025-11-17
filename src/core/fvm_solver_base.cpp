#include "core/fvm_solver_base.hpp"
#include "physics/physics_factory.hpp"
#include "spatial/flux_calculation/flux_calculator_factory.hpp"
#include "temporal/time_integrator_factory.hpp"
#include "spatial/reconstruction/reconstruction_factory.hpp"
#include "boundary/periodic_bc.hpp"
#include "boundary/reflective_bc.hpp"
#include "boundary/transmissive_bc.hpp"

namespace fvm3d::core {

void FVMSolverBase::initialize_physics(
    const std::string& physics_type,
    double mhd_resistivity,
    double mhd_eta0,
    double mhd_eta1,
    double mhd_localization_scale,
    double mhd_glm_ch,
    double mhd_glm_cr
) {
    // Initialize physics using PhysicsFactory
    if (physics_type == "mhd_advanced") {
        // Advanced MHD with configurable resistivity model and GLM parameters
        physics::AdvancedResistiveMHD3D::ResistivityModel resistivity;
        resistivity.eta0 = mhd_eta0;
        resistivity.eta1 = mhd_eta1;
        resistivity.localization_scale = mhd_localization_scale;

        physics::AdvancedResistiveMHD3D::GLMParameters glm(
            mhd_glm_ch,
            mhd_glm_cr
        );

        physics_ = physics::PhysicsFactory::create_advanced_mhd(resistivity, glm);
    } else {
        // Euler or basic MHD
        physics_ = physics::PhysicsFactory::create(
            physics_type,
            mhd_resistivity
        );
    }
}

void FVMSolverBase::initialize_flux_calculator(const std::string& flux_calculator_name) {
    flux_calculator_ = spatial::FluxCalculatorFactory::create(flux_calculator_name, physics_);
}

void FVMSolverBase::initialize_reconstruction(
    const std::string& reconstruction_name,
    int num_vars,
    const std::string& limiter_name
) {
    reconstruction_ = spatial::ReconstructionFactory::create(
        reconstruction_name,
        num_vars,
        limiter_name
    );
}

void FVMSolverBase::initialize_time_integrator(const std::string& integrator_name) {
    time_integrator_ = temporal::TimeIntegratorFactory::create(integrator_name);
}

void FVMSolverBase::initialize_boundary_conditions(
    const std::string& bc_type,
    bool bc_x, bool bc_y, bool bc_z
) {
    if (bc_type == "periodic") {
        boundary_condition_ = std::make_unique<boundary::PeriodicBC>(
            physics_, bc_x, bc_y, bc_z
        );
    } else if (bc_type == "reflective") {
        boundary_condition_ = std::make_unique<boundary::ReflectiveBC>(
            physics_, bc_x, bc_y, bc_z
        );
    } else if (bc_type == "transmissive") {
        boundary_condition_ = std::make_unique<boundary::TransmissiveBC>(
            physics_, bc_x, bc_y, bc_z
        );
    } else {
        throw std::invalid_argument("Unknown boundary condition: " + bc_type);
    }
}

void FVMSolverBase::reconstruct_1d(
    const StateField3D& state,
    int direction,
    int i, int j, int k,
    Eigen::VectorXd& U_L,
    Eigen::VectorXd& U_R
) const {
    // Use reconstruction API that operates directly on Field3D
    // Reconstructs left and right states at interface (i+1/2, j, k) for direction=0
    reconstruction_->reconstruct(state, i, j, k, direction, U_L, U_R);
}

} // namespace fvm3d::core
