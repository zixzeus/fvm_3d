#include "parallel/mpi_global_reduction.hpp"
#include "parallel/mpi_utils.hpp"
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace fvm3d::parallel {

MPIGlobalReduction::MPIGlobalReduction() {
    // Verify MPI is initialized
    int is_initialized = 0;
    MPI_Initialized(&is_initialized);
    if (!is_initialized) {
        throw std::runtime_error("MPI must be initialized before MPIGlobalReduction");
    }
}

MPIGlobalReduction::~MPIGlobalReduction() {
    // Nothing to clean up - MPI finalization handled by MPIUtils
}

double MPIGlobalReduction::global_min(double local_value) const {
    double result = 0.0;
    int error = MPI_Allreduce(
        &local_value,
        &result,
        1,
        MPI_DOUBLE,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    MPIUtils::check_mpi_error(error, "MPI_Allreduce (MIN)");
    return result;
}

double MPIGlobalReduction::global_max(double local_value) const {
    double result = 0.0;
    int error = MPI_Allreduce(
        &local_value,
        &result,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    MPIUtils::check_mpi_error(error, "MPI_Allreduce (MAX)");
    return result;
}

double MPIGlobalReduction::global_sum(double local_value) const {
    double result = 0.0;
    int error = MPI_Allreduce(
        &local_value,
        &result,
        1,
        MPI_DOUBLE,
        MPI_SUM,
        MPI_COMM_WORLD
    );
    MPIUtils::check_mpi_error(error, "MPI_Allreduce (SUM)");
    return result;
}

double MPIGlobalReduction::reduce_cfl_timestep(double local_cfl) const {
    // CFL condition requires dt = min over all processes
    // Safety: ensure all processes have positive dt
    if (local_cfl <= 0.0) {
        throw std::runtime_error("CFL time step must be positive: " + std::to_string(local_cfl));
    }

    // Global minimum determines the time step
    return global_min(local_cfl);
}

MPIGlobalReduction::Statistics MPIGlobalReduction::gather_statistics(
    const Statistics& local_stats
) const {
    Statistics global_stats;

    // Gather minimum density
    {
        struct MinVal {
            double value;
            int rank;
        } local_min, global_min;

        local_min.value = local_stats.rho_min;
        local_min.rank = MPIUtils::rank();

        int error = MPI_Allreduce(
            &local_min, &global_min, 1,
            MPI_DOUBLE_INT, MPI_MINLOC,
            MPI_COMM_WORLD
        );
        MPIUtils::check_mpi_error(error, "MPI_Allreduce (rho_min)");
        global_stats.rho_min = global_min.value;
    }

    // Gather maximum density
    {
        struct MaxVal {
            double value;
            int rank;
        } local_max, global_max;

        local_max.value = local_stats.rho_max;
        local_max.rank = MPIUtils::rank();

        int error = MPI_Allreduce(
            &local_max, &global_max, 1,
            MPI_DOUBLE_INT, MPI_MAXLOC,
            MPI_COMM_WORLD
        );
        MPIUtils::check_mpi_error(error, "MPI_Allreduce (rho_max)");
        global_stats.rho_max = global_max.value;
    }

    // Gather min/max pressure
    global_stats.p_min = global_min(local_stats.p_min);
    global_stats.p_max = global_max(local_stats.p_max);

    // Gather maximum velocity and magnetic field
    global_stats.u_max = global_max(local_stats.u_max);
    global_stats.B_max = global_max(local_stats.B_max);

    return global_stats;
}

void MPIGlobalReduction::global_barrier() const {
    int error = MPI_Barrier(MPI_COMM_WORLD);
    MPIUtils::check_mpi_error(error, "MPI_Barrier");
}

bool MPIGlobalReduction::gather_to_root(
    const std::vector<double>& local_vec,
    std::vector<double>& global_vec
) const {
    int rank = MPIUtils::rank();
    int size = MPIUtils::size();

    // Get sizes from all processes
    int local_size = local_vec.size();
    std::vector<int> all_sizes(size);
    int error = MPI_Gather(
        &local_size, 1, MPI_INT,
        all_sizes.data(), 1, MPI_INT,
        0, MPI_COMM_WORLD
    );
    MPIUtils::check_mpi_error(error, "MPI_Gather (sizes)");

    // Compute displacements on root
    std::vector<int> displacements(size);
    if (rank == 0) {
        displacements[0] = 0;
        for (int i = 1; i < size; i++) {
            displacements[i] = displacements[i-1] + all_sizes[i-1];
        }
        int total_size = displacements[size-1] + all_sizes[size-1];
        global_vec.resize(total_size);
    }

    // Gather data to root
    error = MPI_Gatherv(
        const_cast<double*>(local_vec.data()), local_size, MPI_DOUBLE,
        global_vec.data(), all_sizes.data(), displacements.data(), MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );
    MPIUtils::check_mpi_error(error, "MPI_Gatherv");

    return (rank == 0);
}

void MPIGlobalReduction::broadcast_from_root(double& value) const {
    int error = MPI_Bcast(&value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPIUtils::check_mpi_error(error, "MPI_Bcast (double)");
}

void MPIGlobalReduction::broadcast_from_root(std::vector<double>& vec) const {
    // Broadcast size first
    int size = vec.size();
    int error = MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPIUtils::check_mpi_error(error, "MPI_Bcast (size)");

    // Resize on non-root processes
    if (!MPIUtils::is_root()) {
        vec.resize(size);
    }

    // Broadcast data
    error = MPI_Bcast(vec.data(), size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPIUtils::check_mpi_error(error, "MPI_Bcast (data)");
}

bool MPIGlobalReduction::is_globally_converged(
    double local_residual,
    double tolerance
) const {
    // Check if local process has converged
    bool local_converged = (local_residual <= tolerance);

    // Use MPI_Allreduce with MPI_LAND (logical AND)
    int local_int = local_converged ? 1 : 0;
    int global_int = 0;

    int error = MPI_Allreduce(
        &local_int, &global_int, 1,
        MPI_INT, MPI_LAND,
        MPI_COMM_WORLD
    );
    MPIUtils::check_mpi_error(error, "MPI_Allreduce (convergence)");

    return (global_int == 1);
}

void MPIGlobalReduction::print_statistics(const Statistics& stats) const {
    if (!MPIUtils::is_root()) {
        return;
    }

    std::cout << "\n=== Global Statistics ===" << std::endl;
    std::cout << "Density (ρ):       [" << stats.rho_min << ", " << stats.rho_max << "]" << std::endl;
    std::cout << "Pressure (p):      [" << stats.p_min << ", " << stats.p_max << "]" << std::endl;
    std::cout << "Max velocity:      " << stats.u_max << std::endl;
    std::cout << "Max B-field:       " << stats.B_max << std::endl;
}

void MPIGlobalReduction::mpi_minloc(void* invec, void* inoutvec, int* len, MPI_Datatype* dtype) {
    // Custom reduction: minimum with location tracking
    // Not used in current implementation but available for future optimization
}

void MPIGlobalReduction::mpi_maxloc(void* invec, void* inoutvec, int* len, MPI_Datatype* dtype) {
    // Custom reduction: maximum with location tracking
    // Not used in current implementation but available for future optimization
}

} // namespace fvm3d::parallel
