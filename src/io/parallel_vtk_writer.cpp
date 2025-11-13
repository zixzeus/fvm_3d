#include "io/parallel_vtk_writer.hpp"
#include "physics/euler3d.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>

namespace fvm3d::io {

ParallelVTKWriter::ParallelVTKWriter(
    const std::string& base_filename,
    const parallel::MPIDomainDecomposer& decomposer,
    MPI_Comm comm)
    : base_filename_(base_filename),
      decomposer_(decomposer),
      comm_(comm),
      local_writer_("", false),  // Will be set later
      grid_set_(false)
{
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);

    // Get global grid information
    auto global_grid = decomposer_.global_grid();
    global_nx_ = global_grid[0];
    global_ny_ = global_grid[1];
    global_nz_ = global_grid[2];
}

void ParallelVTKWriter::set_local_grid(const core::Grid3D& local_grid) {
    const auto& geom = local_grid.geometry();

    // Set up local VTK writer with rank-specific filename
    std::string piece_filename = get_piece_filename(rank_);
    local_writer_ = VTKWriter(piece_filename, false);  // ASCII for now
    local_writer_.set_grid(local_grid);

    // Store global grid spacing (assume uniform)
    global_dx_ = geom.dx;
    global_dy_ = geom.dy;
    global_dz_ = geom.dz;

    // Compute global origin from first cell of rank 0
    auto idx_range = decomposer_.global_index_range();
    global_xmin_ = geom.xmin - idx_range.i_min * global_dx_;
    global_ymin_ = geom.ymin - idx_range.j_min * global_dy_;
    global_zmin_ = geom.zmin - idx_range.k_min * global_dz_;

    grid_set_ = true;
}

void ParallelVTKWriter::add_scalar_field(
    const std::string& name,
    const double* data,
    int nx, int ny, int nz)
{
    local_writer_.add_scalar_field(name, data, nx, ny, nz);

    // Track field name for .pvti metadata (avoid duplicates)
    if (std::find(scalar_field_names_.begin(), scalar_field_names_.end(), name)
        == scalar_field_names_.end()) {
        scalar_field_names_.push_back(name);
    }
}

void ParallelVTKWriter::add_vector_field(
    const std::string& name,
    const double* data_x,
    const double* data_y,
    const double* data_z,
    int nx, int ny, int nz)
{
    local_writer_.add_vector_field(name, data_x, data_y, data_z, nx, ny, nz);

    // Track field name
    if (std::find(vector_field_names_.begin(), vector_field_names_.end(), name)
        == vector_field_names_.end()) {
        vector_field_names_.push_back(name);
    }
}

void ParallelVTKWriter::add_scalar_from_state(
    const std::string& name,
    const core::StateField3D& field,
    int var_index,
    const core::Grid3D& grid)
{
    local_writer_.add_scalar_from_state(name, field, var_index, grid);

    if (std::find(scalar_field_names_.begin(), scalar_field_names_.end(), name)
        == scalar_field_names_.end()) {
        scalar_field_names_.push_back(name);
    }
}

void ParallelVTKWriter::add_vector_from_state(
    const std::string& name,
    const core::StateField3D& field,
    int var_x, int var_y, int var_z,
    const core::Grid3D& grid)
{
    local_writer_.add_vector_from_state(name, field, var_x, var_y, var_z, grid);

    if (std::find(vector_field_names_.begin(), vector_field_names_.end(), name)
        == vector_field_names_.end()) {
        vector_field_names_.push_back(name);
    }
}

void ParallelVTKWriter::write_parallel() {
    if (!grid_set_) {
        throw std::runtime_error("Grid must be set before writing");
    }

    // Each rank writes its local subdomain
    local_writer_.write();

    // Synchronize before rank 0 writes master file
    MPI_Barrier(comm_);

    // Rank 0 writes the master .pvti file
    if (rank_ == 0) {
        write_master_file();
    }

    MPI_Barrier(comm_);

    if (rank_ == 0) {
        std::cout << "Parallel VTK files written:\n";
        std::cout << "  Master file: " << base_filename_ << ".pvti\n";
        std::cout << "  " << size_ << " subdomain files (.vti)\n";
    }
}

void ParallelVTKWriter::clear_fields() {
    local_writer_.clear_fields();
    scalar_field_names_.clear();
    vector_field_names_.clear();
}

std::string ParallelVTKWriter::get_piece_filename(int rank) const {
    std::ostringstream oss;
    oss << base_filename_ << "_rank" << std::setw(4) << std::setfill('0') << rank << ".vti";
    return oss.str();
}

ParallelVTKWriter::Extent ParallelVTKWriter::get_rank_extent(int rank) const {
    // Query decomposer for rank's global index range
    // For now, compute it manually (would need to extend decomposer API)

    auto proc_grid = decomposer_.process_grid();
    int px = proc_grid[0];
    int py = proc_grid[1];
    int pz = proc_grid[2];

    // Compute process coordinates from rank
    int proc_x = rank / (py * pz);
    int proc_y = (rank / pz) % py;
    int proc_z = rank % pz;

    // Compute local sizes (simple division)
    auto compute_local_size = [](int global_size, int num_procs, int my_rank) {
        int base_size = global_size / num_procs;
        int remainder = global_size % num_procs;
        return base_size + (my_rank < remainder ? 1 : 0);
    };

    int local_nx = compute_local_size(global_nx_, px, proc_x);
    int local_ny = compute_local_size(global_ny_, py, proc_y);
    int local_nz = compute_local_size(global_nz_, pz, proc_z);

    // Compute offsets
    int offset_x = 0;
    for (int p = 0; p < proc_x; p++) {
        offset_x += compute_local_size(global_nx_, px, p);
    }

    int offset_y = 0;
    for (int p = 0; p < proc_y; p++) {
        offset_y += compute_local_size(global_ny_, py, p);
    }

    int offset_z = 0;
    for (int p = 0; p < proc_z; p++) {
        offset_z += compute_local_size(global_nz_, pz, p);
    }

    Extent extent;
    extent.i_min = offset_x;
    extent.i_max = offset_x + local_nx;
    extent.j_min = offset_y;
    extent.j_max = offset_y + local_ny;
    extent.k_min = offset_z;
    extent.k_max = offset_z + local_nz;

    return extent;
}

void ParallelVTKWriter::write_master_file() {
    std::string master_filename = base_filename_ + ".pvti";
    std::ofstream file(master_filename);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open master file: " + master_filename);
    }

    // Write XML header
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"PImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";

    // PImageData element with global grid info
    file << "  <PImageData WholeExtent=\""
         << "0 " << global_nx_ << " "
         << "0 " << global_ny_ << " "
         << "0 " << global_nz_ << "\" "
         << "GhostLevel=\"0\" "
         << "Origin=\"" << global_xmin_ << " " << global_ymin_ << " " << global_zmin_ << "\" "
         << "Spacing=\"" << global_dx_ << " " << global_dy_ << " " << global_dz_ << "\">\n";

    // PCellData (metadata for cell data)
    file << "    <PCellData>\n";

    // List all scalar fields
    for (const auto& name : scalar_field_names_) {
        file << "      <PDataArray type=\"Float64\" Name=\"" << name << "\"/>\n";
    }

    // List all vector fields
    for (const auto& name : vector_field_names_) {
        file << "      <PDataArray type=\"Float64\" Name=\"" << name
             << "\" NumberOfComponents=\"3\"/>\n";
    }

    file << "    </PCellData>\n";

    // List all piece files (subdomains)
    for (int r = 0; r < size_; r++) {
        Extent extent = get_rank_extent(r);
        std::string piece_file = get_piece_filename(r);

        // Extract just the filename (not full path)
        size_t last_slash = piece_file.find_last_of("/\\");
        std::string piece_basename = (last_slash != std::string::npos) ?
            piece_file.substr(last_slash + 1) : piece_file;

        file << "    <Piece Extent=\""
             << extent.i_min << " " << extent.i_max << " "
             << extent.j_min << " " << extent.j_max << " "
             << extent.k_min << " " << extent.k_max << "\" "
             << "Source=\"" << piece_basename << "\"/>\n";
    }

    file << "  </PImageData>\n";
    file << "</VTKFile>\n";

    file.close();
}

// ========== Utility Function ==========

void export_parallel_state_to_vtk(
    const std::string& base_filename,
    const core::StateField3D& state,
    const core::Grid3D& local_grid,
    const parallel::MPIDomainDecomposer& decomposer,
    int num_vars,
    const std::string& physics_type,
    double time,
    MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    ParallelVTKWriter writer(base_filename, decomposer, comm);
    writer.set_local_grid(local_grid);

    int nx = local_grid.nx_local();
    int ny = local_grid.ny_local();
    int nz = local_grid.nz_local();
    int ng = local_grid.nghost();

    if (physics_type == "euler" && num_vars == 5) {
        // Conservative variables
        writer.add_scalar_from_state("density", state, 0, local_grid);
        writer.add_scalar_from_state("momentum_x", state, 1, local_grid);
        writer.add_scalar_from_state("momentum_y", state, 2, local_grid);
        writer.add_scalar_from_state("momentum_z", state, 3, local_grid);
        writer.add_scalar_from_state("energy", state, 4, local_grid);

        // Compute primitive variables
        std::vector<double> pressure(nx * ny * nz);
        std::vector<double> vel_x(nx * ny * nz);
        std::vector<double> vel_y(nx * ny * nz);
        std::vector<double> vel_z(nx * ny * nz);

        const double gamma = 1.4;

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
                    } else {
                        vel_x[idx] = 0.0;
                        vel_y[idx] = 0.0;
                        vel_z[idx] = 0.0;
                        pressure[idx] = 0.0;
                    }
                }
            }
        }

        writer.add_scalar_field("pressure", pressure.data(), nx, ny, nz);
        writer.add_vector_field("velocity", vel_x.data(), vel_y.data(), vel_z.data(), nx, ny, nz);

    } else if (physics_type == "mhd" && (num_vars == 8 || num_vars == 9)) {
        // Conservative variables
        writer.add_scalar_from_state("density", state, 0, local_grid);
        writer.add_scalar_from_state("energy", state, 4, local_grid);

        // Magnetic field
        writer.add_vector_from_state("magnetic_field", state, 5, 6, 7, local_grid);

        if (num_vars == 9) {
            writer.add_scalar_from_state("psi", state, 8, local_grid);
        }

        // Compute primitive variables
        std::vector<double> pressure(nx * ny * nz);
        std::vector<double> vel_x(nx * ny * nz);
        std::vector<double> vel_y(nx * ny * nz);
        std::vector<double> vel_z(nx * ny * nz);
        std::vector<double> B_magnitude(nx * ny * nz);

        const double gamma = 5.0/3.0;

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
                }
            }
        }

        writer.add_scalar_field("pressure", pressure.data(), nx, ny, nz);
        writer.add_scalar_field("B_magnitude", B_magnitude.data(), nx, ny, nz);
        writer.add_vector_field("velocity", vel_x.data(), vel_y.data(), vel_z.data(), nx, ny, nz);
    }

    writer.write_parallel();
}

} // namespace fvm3d::io
