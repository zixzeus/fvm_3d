#include "core/grid3d.hpp"

namespace fvm3d::core {

Grid3D::Grid3D(const GridGeometry3D& geom, int nghost)
    : geom_(geom), nghost_(nghost) {
    if (nghost < 0) {
        throw std::invalid_argument("Number of ghost cells must be non-negative");
    }
}

double Grid3D::cell_center_x(int i) const {
    // Cell centers are at xmin + (0.5 + i - nghost) * dx for i in [0, nx_total)
    int i_interior = i - nghost_;
    return geom_.xmin + (0.5 + i_interior) * geom_.dx;
}

double Grid3D::cell_center_y(int j) const {
    // Cell centers are at ymin + (0.5 + j - nghost) * dy for j in [0, ny_total)
    int j_interior = j - nghost_;
    return geom_.ymin + (0.5 + j_interior) * geom_.dy;
}

double Grid3D::cell_center_z(int k) const {
    // Cell centers are at zmin + (0.5 + k - nghost) * dz for k in [0, nz_total)
    int k_interior = k - nghost_;
    return geom_.zmin + (0.5 + k_interior) * geom_.dz;
}

bool Grid3D::is_ghost(int i, int j, int k) const {
    return i < nghost_ || i >= geom_.nx + nghost_ ||
           j < nghost_ || j >= geom_.ny + nghost_ ||
           k < nghost_ || k >= geom_.nz + nghost_;
}

} // namespace fvm3d::core
