#pragma once

#include "core/field3d.hpp"
#include "core/grid3d.hpp"
#include <string>

namespace fvm3d::boundary {

using StateField3D = fvm3d::core::StateField3D;
using Grid3D = fvm3d::core::Grid3D;

/**
 * Abstract base class for boundary conditions.
 * Fills ghost cells based on interior values.
 */
class BoundaryCondition {
public:
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
    static constexpr double GAMMA = 1.4;
    static constexpr double RHO_FLOOR = 1e-10;
    static constexpr double P_FLOOR = 1e-11;

    void conservative_to_primitive(
        const StateField3D& state, int i, int j, int k,
        double& rho, double& u, double& v, double& w, double& p
    ) const {
        rho = std::max(state(0, i, j, k), RHO_FLOOR);
        u = state(1, i, j, k) / rho;
        v = state(2, i, j, k) / rho;
        w = state(3, i, j, k) / rho;

        double kinetic_energy = 0.5 * (state(1, i, j, k)*state(1, i, j, k) +
                                       state(2, i, j, k)*state(2, i, j, k) +
                                       state(3, i, j, k)*state(3, i, j, k)) / rho;
        double internal_energy = state(4, i, j, k) / rho - kinetic_energy;
        p = std::max((GAMMA - 1.0) * rho * internal_energy, P_FLOOR);
    }
};

} // namespace fvm3d::boundary

#include <algorithm>
#include <cmath>
