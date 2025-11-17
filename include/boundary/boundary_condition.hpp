#pragma once

#include "core/field3d.hpp"
#include "core/grid3d.hpp"
#include "physics/physics_base.hpp"
#include <string>
#include <memory>

namespace fvm3d::boundary {

using StateField3D = fvm3d::core::StateField3D;
using Grid3D = fvm3d::core::Grid3D;

/**
 * Abstract base class for boundary conditions.
 * Fills ghost cells based on interior values.
 */
class BoundaryCondition {
public:
    /**
     * Constructor with physics object.
     * @param physics: Physics equation system for getting constants
     */
    explicit BoundaryCondition(const std::shared_ptr<physics::PhysicsBase>& physics)
        : physics_(physics) {}

    virtual ~BoundaryCondition() = default;

    /**
     * Apply boundary condition to all ghost layers.
     */
    virtual void apply(
        StateField3D& state,
        const Grid3D& grid
    ) = 0;

    /**
     * Get name of boundary condition.
     */
    virtual std::string name() const = 0;

protected:
    std::shared_ptr<physics::PhysicsBase> physics_;  ///< Physics object for constants

    /**
     * Convert conservative to primitive variables using physics object.
     * Uses physics object's conversion method for accuracy.
     */
    void conservative_to_primitive(
        const StateField3D& state, int i, int j, int k,
        double& rho, double& u, double& v, double& w, double& p
    ) const {
        // Extract conservative state at (i,j,k)
        int nvars = physics_->num_variables();
        Eigen::VectorXd U(nvars);
        for (int var = 0; var < nvars; ++var) {
            U(var) = state(var, i, j, k);
        }

        // Convert using physics object
        Eigen::VectorXd V = physics_->conservative_to_primitive(U);

        // Extract primitive variables (common to Euler and MHD)
        rho = V(0);
        u = V(1);
        v = V(2);
        w = V(3);
        p = V(4);
    }
};

} // namespace fvm3d::boundary

#include <algorithm>
#include <cmath>
