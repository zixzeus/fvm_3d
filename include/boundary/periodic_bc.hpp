#pragma once

#include "boundary_condition.hpp"

namespace fvm3d::boundary {

/**
 * Periodic boundary condition.
 * Ghost cells filled from opposite side of domain.
 */
class PeriodicBC : public BoundaryCondition {
public:
    /**
     * Constructor specifies which directions have periodic boundaries.
     * @param physics: Physics object for variable conversion
     * @param periodic_x/y/z: Enable periodic BC in each direction
     */
    PeriodicBC(
        const std::shared_ptr<physics::PhysicsBase>& physics,
        bool periodic_x = true,
        bool periodic_y = false,
        bool periodic_z = false
    )
        : BoundaryCondition(physics),
          periodic_x_(periodic_x),
          periodic_y_(periodic_y),
          periodic_z_(periodic_z) {}

    void apply(StateField3D& state, const Grid3D& grid) override;
    std::string name() const override { return "Periodic"; }

private:
    bool periodic_x_, periodic_y_, periodic_z_;

    void apply_periodic_x(StateField3D& state, const Grid3D& grid);
    void apply_periodic_y(StateField3D& state, const Grid3D& grid);
    void apply_periodic_z(StateField3D& state, const Grid3D& grid);
};

} // namespace fvm3d::boundary
