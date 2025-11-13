#pragma once

#include "boundary_condition.hpp"

namespace fvm3d::boundary {

/**
 * Periodic boundary condition.
 * Ghost cells filled from opposite side of domain.
 */
class PeriodicBC : public BoundaryCondition {
public:
    PeriodicBC(bool periodic_x = true, bool periodic_y = false, bool periodic_z = false)
        : periodic_x_(periodic_x), periodic_y_(periodic_y), periodic_z_(periodic_z) {}

    void apply(StateField3D& state, const Grid3D& grid) override;
    std::string name() const override { return "Periodic"; }

private:
    bool periodic_x_, periodic_y_, periodic_z_;

    void apply_periodic_x(StateField3D& state, const Grid3D& grid);
    void apply_periodic_y(StateField3D& state, const Grid3D& grid);
    void apply_periodic_z(StateField3D& state, const Grid3D& grid);
};

} // namespace fvm3d::boundary
