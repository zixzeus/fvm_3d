#pragma once

#include <stdexcept>

namespace fvm3d::core {

/**
 * Geometric description of a 3D Cartesian grid.
 * Domain spans [xmin, xmin+Lx] x [ymin, ymin+Ly] x [zmin, zmin+Lz]
 */
struct GridGeometry3D {
    double xmin, ymin, zmin;  // Minimum coordinates (domain origin)
    double Lx, Ly, Lz;        // Domain dimensions (extent in each direction)
    int nx, ny, nz;           // Number of cells in each direction

    double dx, dy, dz;        // Cell spacing (computed from dimensions and cell counts)

    /**
     * Constructor with explicit minimum coordinates.
     */
    GridGeometry3D(double xmin_val, double ymin_val, double zmin_val,
                   double lx, double ly, double lz,
                   int nx_cells, int ny_cells, int nz_cells)
        : xmin(xmin_val), ymin(ymin_val), zmin(zmin_val),
          Lx(lx), Ly(ly), Lz(lz), nx(nx_cells), ny(ny_cells), nz(nz_cells) {
        if (nx <= 0 || ny <= 0 || nz <= 0) {
            throw std::invalid_argument("Grid dimensions must be positive");
        }
        if (Lx <= 0 || Ly <= 0 || Lz <= 0) {
            throw std::invalid_argument("Domain size must be positive");
        }
        dx = Lx / nx;
        dy = Ly / ny;
        dz = Lz / nz;
    }

    /**
     * Convenience constructor: assumes domain starts at origin (0, 0, 0).
     */
    GridGeometry3D(double lx, double ly, double lz,
                   int nx_cells, int ny_cells, int nz_cells)
        : xmin(0.0), ymin(0.0), zmin(0.0),
          Lx(lx), Ly(ly), Lz(lz), nx(nx_cells), ny(ny_cells), nz(nz_cells) {
        if (nx <= 0 || ny <= 0 || nz <= 0) {
            throw std::invalid_argument("Grid dimensions must be positive");
        }
        if (Lx <= 0 || Ly <= 0 || Lz <= 0) {
            throw std::invalid_argument("Domain size must be positive");
        }
        dx = Lx / nx;
        dy = Ly / ny;
        dz = Lz / nz;
    }

    /**
     * Get maximum coordinates.
     */
    double xmax() const { return xmin + Lx; }
    double ymax() const { return ymin + Ly; }
    double zmax() const { return zmin + Lz; }
};

/**
 * Cartesian grid in 3D with ghost cells.
 */
class Grid3D {
public:
    /**
     * Constructor: Create a grid with specified geometry and ghost cell width.
     */
    Grid3D(const GridGeometry3D& geom, int nghost = 2);

    // Accessors
    const GridGeometry3D& geometry() const { return geom_; }
    int nx_local() const { return geom_.nx; }
    int ny_local() const { return geom_.ny; }
    int nz_local() const { return geom_.nz; }
    int nx_total() const { return geom_.nx + 2 * nghost_; }
    int ny_total() const { return geom_.ny + 2 * nghost_; }
    int nz_total() const { return geom_.nz + 2 * nghost_; }
    int nghost() const { return nghost_; }

    // Index range accessors (for interior cells, excluding ghost)
    int i_begin() const { return nghost_; }
    int i_end() const { return geom_.nx + nghost_; }
    int j_begin() const { return nghost_; }
    int j_end() const { return geom_.ny + nghost_; }
    int k_begin() const { return nghost_; }
    int k_end() const { return geom_.nz + nghost_; }

    /**
     * Get cell center coordinates.
     */
    double cell_center_x(int i) const;
    double cell_center_y(int j) const;
    double cell_center_z(int k) const;

    /**
     * Check if a cell is in the ghost region.
     */
    bool is_ghost(int i, int j, int k) const;

private:
    GridGeometry3D geom_;
    int nghost_;
};

} // namespace fvm3d::core
