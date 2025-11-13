#pragma once

#include "boundary_condition.hpp"

namespace fvm3d::boundary {

/**
 * Reflective (slip wall) boundary condition.
 * Ghost cells are filled by mirroring interior with reversed normal velocity.
 *
 * Used for solid walls where normal velocity = 0 but tangential slip is allowed.
 * Perfect for inviscid flow against rigid walls.
 */
class ReflectiveBC : public BoundaryCondition {
public:
    /**
     * Constructor specifies which directions have reflective boundaries.
     */
    ReflectiveBC(bool reflect_x = true, bool reflect_y = false, bool reflect_z = false)
        : reflect_x_(reflect_x), reflect_y_(reflect_y), reflect_z_(reflect_z) {}

    void apply(StateField3D& state, const Grid3D& grid) override;
    std::string name() const override { return "Reflective"; }

private:
    bool reflect_x_, reflect_y_, reflect_z_;

    void apply_reflective_x(StateField3D& state, const Grid3D& grid);
    void apply_reflective_y(StateField3D& state, const Grid3D& grid);
    void apply_reflective_z(StateField3D& state, const Grid3D& grid);
};

} // namespace fvm3d::boundary
