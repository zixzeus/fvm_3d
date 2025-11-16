#pragma once

#include "boundary_condition.hpp"

namespace fvm3d::boundary {

/**
 * Transmissive (zero-gradient outflow) boundary condition.
 * Ghost cells are filled with values from nearest interior cells.
 * Allows material to leave the domain freely with minimal reflections.
 *
 * Simple but effective for open boundaries in outflow regions.
 * Also known as "free-slip" or "zero-gradient" boundary condition.
 */
class TransmissiveBC : public BoundaryCondition {
public:
    /**
     * Constructor specifies which directions have transmissive boundaries.
     * @param physics: Physics object for variable conversion
     * @param transmissive_x/y/z: Enable transmissive BC in each direction
     */
    TransmissiveBC(
        const std::shared_ptr<physics::PhysicsBase>& physics,
        bool transmissive_x = false,
        bool transmissive_y = false,
        bool transmissive_z = false
    )
        : BoundaryCondition(physics),
          transmissive_x_(transmissive_x),
          transmissive_y_(transmissive_y),
          transmissive_z_(transmissive_z) {}

    void apply(StateField3D& state, const Grid3D& grid) override;
    std::string name() const override { return "Transmissive"; }

private:
    bool transmissive_x_, transmissive_y_, transmissive_z_;

    void apply_transmissive_x(StateField3D& state, const Grid3D& grid);
    void apply_transmissive_y(StateField3D& state, const Grid3D& grid);
    void apply_transmissive_z(StateField3D& state, const Grid3D& grid);
};

} // namespace fvm3d::boundary
