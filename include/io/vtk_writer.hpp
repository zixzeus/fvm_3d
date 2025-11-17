#pragma once

#include "core/grid3d.hpp"
#include "core/field3d.hpp"
#include <string>
#include <vector>
#include <fstream>

namespace fvm3d::io {

/**
 * VTK structured grid writer for 3D FVM data.
 *
 * Supports:
 * - XML VTK format (.vti files)
 * - Binary or ASCII output
 * - Multiple scalar and vector fields
 * - Cell-centered data (typical for FVM)
 *
 * Usage:
 *   VTKWriter writer("output.vti");
 *   writer.set_grid(grid);
 *   writer.add_scalar_field("density", rho_data);
 *   writer.add_vector_field("velocity", vel_data);
 *   writer.write();
 *
 * Output can be visualized with ParaView, VisIt, etc.
 */
class VTKWriter {
public:
    /**
     * Constructor.
     * @param filename: Output VTK file path (.vti extension)
     * @param binary: Use binary format (default: true for smaller files)
     */
    explicit VTKWriter(const std::string& filename, bool binary = true);

    /**
     * Set the grid geometry.
     * Must be called before adding fields.
     */
    void set_grid(const core::Grid3D& grid);

    /**
     * Add a scalar field (e.g., density, pressure, temperature).
     * @param name: Field name for visualization
     * @param data: 3D array data (nx x ny x nz)
     * @param nx, ny, nz: Grid dimensions (interior cells only, no ghosts)
     */
    void add_scalar_field(
        const std::string& name,
        const double* data,
        int nx, int ny, int nz
    );

    /**
     * Add a vector field (e.g., velocity, magnetic field).
     * @param name: Field name for visualization
     * @param data_x, data_y, data_z: Components of the vector field
     * @param nx, ny, nz: Grid dimensions (interior cells only)
     */
    void add_vector_field(
        const std::string& name,
        const double* data_x,
        const double* data_y,
        const double* data_z,
        int nx, int ny, int nz
    );

    /**
     * Add scalar field from StateField3D.
     * @param name: Field name
     * @param field: State field
     * @param var_index: Variable index in the field
     * @param grid: Grid (to determine interior region)
     */
    void add_scalar_from_state(
        const std::string& name,
        const core::StateField3D& field,
        int var_index,
        const core::Grid3D& grid
    );

    /**
     * Add vector field from StateField3D.
     * @param name: Field name
     * @param field: State field
     * @param var_x, var_y, var_z: Variable indices for vector components
     * @param grid: Grid (to determine interior region)
     */
    void add_vector_from_state(
        const std::string& name,
        const core::StateField3D& field,
        int var_x, int var_y, int var_z,
        const core::Grid3D& grid
    );

    /**
     * Write VTK file to disk.
     * All fields must be added before calling this.
     */
    void write();

    /**
     * Clear all fields (to reuse writer for different timestep).
     */
    void clear_fields();

private:
    struct ScalarField {
        std::string name;
        std::vector<double> data;
    };

    struct VectorField {
        std::string name;
        std::vector<double> data;  // Interleaved: [x0,y0,z0, x1,y1,z1, ...]
    };

    std::string filename_;
    bool binary_;

    // Grid information
    int nx_, ny_, nz_;
    double xmin_, ymin_, zmin_;
    double dx_, dy_, dz_;
    bool grid_set_;

    // Field data
    std::vector<ScalarField> scalar_fields_;
    std::vector<VectorField> vector_fields_;

    /**
     * Write XML header.
     */
    void write_header(std::ofstream& file);

    /**
     * Write grid points.
     */
    void write_points(std::ofstream& file);

    /**
     * Write scalar field data.
     */
    void write_scalar_data(std::ofstream& file, const ScalarField& field);

    /**
     * Write vector field data.
     */
    void write_vector_data(std::ofstream& file, const VectorField& field);

    /**
     * Encode binary data to base64 (for XML VTK binary format).
     */
    std::string encode_base64(const void* data, size_t size);
};

/**
 * Utility function: Export FVM solver state to VTK file.
 *
 * @param filename: Output file path
 * @param state: Conservative state field
 * @param grid: Grid geometry
 * @param num_vars: Number of variables (5 for Euler, 8/9 for MHD)
 * @param physics_type: "euler" or "mhd"
 * @param time: Current simulation time (added as field data)
 * @param binary: Use binary format
 */
void export_state_to_vtk(
    const std::string& filename,
    const core::StateField3D& state,
    const core::Grid3D& grid,
    int num_vars,
    const std::string& physics_type,
    double time = 0.0,
    bool binary = true
);

} // namespace fvm3d::io
