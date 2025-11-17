#include "io/vtk_writer.hpp"
#include "physics/euler3d.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstring>

namespace fvm3d::io {

VTKWriter::VTKWriter(const std::string& filename, bool binary)
    : filename_(filename),
      binary_(binary),
      nx_(0), ny_(0), nz_(0),
      xmin_(0.0), ymin_(0.0), zmin_(0.0),
      dx_(0.0), dy_(0.0), dz_(0.0),
      grid_set_(false)
{
}

void VTKWriter::set_grid(const core::Grid3D& grid) {
    const auto& geom = grid.geometry();

    nx_ = geom.nx;
    ny_ = geom.ny;
    nz_ = geom.nz;

    xmin_ = geom.xmin;
    ymin_ = geom.ymin;
    zmin_ = geom.zmin;

    dx_ = geom.dx;
    dy_ = geom.dy;
    dz_ = geom.dz;

    grid_set_ = true;
}

void VTKWriter::add_scalar_field(
    const std::string& name,
    const double* data,
    int nx, int ny, int nz)
{
    if (!grid_set_) {
        throw std::runtime_error("Grid must be set before adding fields");
    }

    if (nx != nx_ || ny != ny_ || nz != nz_) {
        throw std::runtime_error("Field dimensions must match grid dimensions");
    }

    ScalarField field;
    field.name = name;
    field.data.resize(nx * ny * nz);

    // Copy data (row-major order: k varies fastest, then j, then i)
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int idx = k + nz * (j + ny * i);  // VTK order
                int src_idx = k + nz * (j + ny * i);
                field.data[idx] = data[src_idx];
            }
        }
    }

    scalar_fields_.push_back(field);
}

void VTKWriter::add_vector_field(
    const std::string& name,
    const double* data_x,
    const double* data_y,
    const double* data_z,
    int nx, int ny, int nz)
{
    if (!grid_set_) {
        throw std::runtime_error("Grid must be set before adding fields");
    }

    if (nx != nx_ || ny != ny_ || nz != nz_) {
        throw std::runtime_error("Field dimensions must match grid dimensions");
    }

    VectorField field;
    field.name = name;
    field.data.resize(3 * nx * ny * nz);

    // Interleave vector components: [x0,y0,z0, x1,y1,z1, ...]
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int idx = k + nz * (j + ny * i);
                int vec_idx = 3 * idx;

                field.data[vec_idx + 0] = data_x[idx];
                field.data[vec_idx + 1] = data_y[idx];
                field.data[vec_idx + 2] = data_z[idx];
            }
        }
    }

    vector_fields_.push_back(field);
}

void VTKWriter::add_scalar_from_state(
    const std::string& name,
    const core::StateField3D& field,
    int var_index,
    const core::Grid3D& grid)
{
    int nx = grid.nx_local();
    int ny = grid.ny_local();
    int nz = grid.nz_local();
    int ng = grid.nghost();

    std::vector<double> data(nx * ny * nz);

    // Extract interior cells only
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int idx = k + nz * (j + ny * i);
                data[idx] = field(var_index, i + ng, j + ng, k + ng);
            }
        }
    }

    add_scalar_field(name, data.data(), nx, ny, nz);
}

void VTKWriter::add_vector_from_state(
    const std::string& name,
    const core::StateField3D& field,
    int var_x, int var_y, int var_z,
    const core::Grid3D& grid)
{
    int nx = grid.nx_local();
    int ny = grid.ny_local();
    int nz = grid.nz_local();
    int ng = grid.nghost();

    std::vector<double> data_x(nx * ny * nz);
    std::vector<double> data_y(nx * ny * nz);
    std::vector<double> data_z(nx * ny * nz);

    // Extract interior cells only
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int idx = k + nz * (j + ny * i);
                data_x[idx] = field(var_x, i + ng, j + ng, k + ng);
                data_y[idx] = field(var_y, i + ng, j + ng, k + ng);
                data_z[idx] = field(var_z, i + ng, j + ng, k + ng);
            }
        }
    }

    add_vector_field(name, data_x.data(), data_y.data(), data_z.data(), nx, ny, nz);
}

void VTKWriter::write() {
    if (!grid_set_) {
        throw std::runtime_error("Grid must be set before writing");
    }

    std::ofstream file(filename_);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename_);
    }

    // Write XML header
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";

    // ImageData element with grid extents and spacing
    file << "  <ImageData WholeExtent=\""
         << "0 " << nx_ << " "
         << "0 " << ny_ << " "
         << "0 " << nz_ << "\" "
         << "Origin=\"" << xmin_ << " " << ymin_ << " " << zmin_ << "\" "
         << "Spacing=\"" << dx_ << " " << dy_ << " " << dz_ << "\">\n";

    // Piece element (single piece for serial output)
    file << "    <Piece Extent=\""
         << "0 " << nx_ << " "
         << "0 " << ny_ << " "
         << "0 " << nz_ << "\">\n";

    // CellData (for cell-centered FVM data)
    file << "      <CellData>\n";

    // Write scalar fields
    for (const auto& field : scalar_fields_) {
        write_scalar_data(file, field);
    }

    // Write vector fields
    for (const auto& field : vector_fields_) {
        write_vector_data(file, field);
    }

    file << "      </CellData>\n";
    file << "    </Piece>\n";
    file << "  </ImageData>\n";
    file << "</VTKFile>\n";

    file.close();

    std::cout << "VTK file written: " << filename_ << "\n";
    std::cout << "  Grid: " << nx_ << " x " << ny_ << " x " << nz_ << "\n";
    std::cout << "  Scalar fields: " << scalar_fields_.size() << "\n";
    std::cout << "  Vector fields: " << vector_fields_.size() << "\n";
}

void VTKWriter::clear_fields() {
    scalar_fields_.clear();
    vector_fields_.clear();
}

void VTKWriter::write_scalar_data(std::ofstream& file, const ScalarField& field) {
    file << "        <DataArray type=\"Float64\" Name=\"" << field.name
         << "\" format=\"ascii\">\n";

    file << "          ";
    for (size_t i = 0; i < field.data.size(); i++) {
        file << std::scientific << std::setprecision(8) << field.data[i];
        if (i < field.data.size() - 1) {
            file << " ";
        }
        // Line break every 6 values for readability
        if ((i + 1) % 6 == 0 && i < field.data.size() - 1) {
            file << "\n          ";
        }
    }
    file << "\n";

    file << "        </DataArray>\n";
}

void VTKWriter::write_vector_data(std::ofstream& file, const VectorField& field) {
    file << "        <DataArray type=\"Float64\" Name=\"" << field.name
         << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    file << "          ";
    for (size_t i = 0; i < field.data.size(); i++) {
        file << std::scientific << std::setprecision(8) << field.data[i];
        if (i < field.data.size() - 1) {
            file << " ";
        }
        // Line break every 6 values (2 complete vectors) for readability
        if ((i + 1) % 6 == 0 && i < field.data.size() - 1) {
            file << "\n          ";
        }
    }
    file << "\n";

    file << "        </DataArray>\n";
}

// ========== Utility Function ==========

void export_state_to_vtk(
    const std::string& filename,
    const core::StateField3D& state,
    const core::Grid3D& grid,
    int num_vars,
    const std::string& physics_type,
    double time,
    bool binary)
{
    VTKWriter writer(filename, binary);
    writer.set_grid(grid);

    int nx = grid.nx_local();
    int ny = grid.ny_local();
    int nz = grid.nz_local();
    int ng = grid.nghost();

    if (physics_type == "euler" && num_vars == 5) {
        // Extract conservative variables
        writer.add_scalar_from_state("density", state, 0, grid);
        writer.add_scalar_from_state("momentum_x", state, 1, grid);
        writer.add_scalar_from_state("momentum_y", state, 2, grid);
        writer.add_scalar_from_state("momentum_z", state, 3, grid);
        writer.add_scalar_from_state("energy", state, 4, grid);

        // Compute and add primitive variables
        std::vector<double> pressure(nx * ny * nz);
        std::vector<double> vel_x(nx * ny * nz);
        std::vector<double> vel_y(nx * ny * nz);
        std::vector<double> vel_z(nx * ny * nz);
        std::vector<double> temperature(nx * ny * nz);

        physics::EulerEquations3D euler;
        const double gamma = 1.4;
        const double R = 1.0;  // Gas constant

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int idx = k + nz * (j + ny * i);

                    double rho = state(0, i + ng, j + ng, k + ng);
                    double rho_u = state(1, i + ng, j + ng, k + ng);
                    double rho_v = state(2, i + ng, j + ng, k + ng);
                    double rho_w = state(3, i + ng, j + ng, k + ng);
                    double E = state(4, i + ng, j + ng, k + ng);

                    if (rho > 1e-10) {
                        vel_x[idx] = rho_u / rho;
                        vel_y[idx] = rho_v / rho;
                        vel_z[idx] = rho_w / rho;

                        double ke = 0.5 * (rho_u * rho_u + rho_v * rho_v + rho_w * rho_w) / rho;
                        pressure[idx] = (gamma - 1.0) * (E - ke);
                        temperature[idx] = pressure[idx] / (rho * R);
                    } else {
                        vel_x[idx] = 0.0;
                        vel_y[idx] = 0.0;
                        vel_z[idx] = 0.0;
                        pressure[idx] = 0.0;
                        temperature[idx] = 0.0;
                    }
                }
            }
        }

        writer.add_scalar_field("pressure", pressure.data(), nx, ny, nz);
        writer.add_scalar_field("temperature", temperature.data(), nx, ny, nz);
        writer.add_vector_field("velocity", vel_x.data(), vel_y.data(), vel_z.data(), nx, ny, nz);

    } else if (physics_type == "mhd" && (num_vars == 8 || num_vars == 9)) {
        // Conservative variables
        writer.add_scalar_from_state("density", state, 0, grid);
        writer.add_scalar_from_state("momentum_x", state, 1, grid);
        writer.add_scalar_from_state("momentum_y", state, 2, grid);
        writer.add_scalar_from_state("momentum_z", state, 3, grid);
        writer.add_scalar_from_state("energy", state, 4, grid);
        writer.add_scalar_from_state("Bx", state, 5, grid);
        writer.add_scalar_from_state("By", state, 6, grid);
        writer.add_scalar_from_state("Bz", state, 7, grid);

        if (num_vars == 9) {
            writer.add_scalar_from_state("psi", state, 8, grid);  // GLM divergence cleaning
        }

        // Compute primitive variables
        std::vector<double> pressure(nx * ny * nz);
        std::vector<double> vel_x(nx * ny * nz);
        std::vector<double> vel_y(nx * ny * nz);
        std::vector<double> vel_z(nx * ny * nz);
        std::vector<double> B_magnitude(nx * ny * nz);
        std::vector<double> current_density_x(nx * ny * nz);
        std::vector<double> current_density_y(nx * ny * nz);
        std::vector<double> current_density_z(nx * ny * nz);

        const double gamma = 5.0/3.0;  // Typical for MHD

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int idx = k + nz * (j + ny * i);

                    double rho = state(0, i + ng, j + ng, k + ng);
                    double rho_u = state(1, i + ng, j + ng, k + ng);
                    double rho_v = state(2, i + ng, j + ng, k + ng);
                    double rho_w = state(3, i + ng, j + ng, k + ng);
                    double E = state(4, i + ng, j + ng, k + ng);
                    double Bx = state(5, i + ng, j + ng, k + ng);
                    double By = state(6, i + ng, j + ng, k + ng);
                    double Bz = state(7, i + ng, j + ng, k + ng);

                    if (rho > 1e-10) {
                        vel_x[idx] = rho_u / rho;
                        vel_y[idx] = rho_v / rho;
                        vel_z[idx] = rho_w / rho;

                        double ke = 0.5 * (rho_u * rho_u + rho_v * rho_v + rho_w * rho_w) / rho;
                        double B2 = Bx * Bx + By * By + Bz * Bz;
                        pressure[idx] = (gamma - 1.0) * (E - ke - 0.5 * B2);
                        B_magnitude[idx] = std::sqrt(B2);
                    } else {
                        vel_x[idx] = 0.0;
                        vel_y[idx] = 0.0;
                        vel_z[idx] = 0.0;
                        pressure[idx] = 0.0;
                        B_magnitude[idx] = std::sqrt(Bx*Bx + By*By + Bz*Bz);
                    }

                    // Current density J = curl(B) (simple finite difference)
                    // For now, set to zero (proper calculation requires neighbor data)
                    current_density_x[idx] = 0.0;
                    current_density_y[idx] = 0.0;
                    current_density_z[idx] = 0.0;
                }
            }
        }

        writer.add_scalar_field("pressure", pressure.data(), nx, ny, nz);
        writer.add_scalar_field("B_magnitude", B_magnitude.data(), nx, ny, nz);
        writer.add_vector_field("velocity", vel_x.data(), vel_y.data(), vel_z.data(), nx, ny, nz);

        // Extract magnetic field components
        std::vector<double> Bx_array(nx * ny * nz);
        std::vector<double> By_array(nx * ny * nz);
        std::vector<double> Bz_array(nx * ny * nz);

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    int idx = k + nz * (j + ny * i);
                    Bx_array[idx] = state(5, i + ng, j + ng, k + ng);
                    By_array[idx] = state(6, i + ng, j + ng, k + ng);
                    Bz_array[idx] = state(7, i + ng, j + ng, k + ng);
                }
            }
        }

        writer.add_vector_field("magnetic_field", Bx_array.data(), By_array.data(), Bz_array.data(), nx, ny, nz);
    }

    writer.write();
}

} // namespace fvm3d::io
