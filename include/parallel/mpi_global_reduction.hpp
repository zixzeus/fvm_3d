#pragma once

#include <mpi.h>
#include <vector>
#include <memory>

namespace fvm3d::parallel {

/**
 * Global MPI reduction operations for distributed FVM solver.
 *
 * Implements collective MPI operations for:
 * - CFL time step reduction (minimum across all processes)
 * - Global statistics gathering (min/max density, pressure)
 * - Global synchronization barriers
 * - Aggregated communication (data gathering to root)
 */
class MPIGlobalReduction {
public:
    /**
     * Statistics structure for physics variables.
     */
    struct Statistics {
        double rho_min, rho_max;
        double p_min, p_max;
        double u_max;  // Maximum velocity magnitude
        double B_max;  // Maximum magnetic field magnitude
    };

    /**
     * Default constructor.
     * Assumes MPI is initialized.
     */
    explicit MPIGlobalReduction();

    /**
     * Destructor.
     */
    ~MPIGlobalReduction();

    /**
     * Perform global minimum reduction on double value.
     * @param local_value: Local process value
     * @return: Minimum across all processes
     */
    double global_min(double local_value) const;

    /**
     * Perform global maximum reduction on double value.
     * @param local_value: Local process value
     * @return: Maximum across all processes
     */
    double global_max(double local_value) const;

    /**
     * Perform global sum reduction on double value.
     * @param local_value: Local process value
     * @return: Sum across all processes
     */
    double global_sum(double local_value) const;

    /**
     * Compute CFL time step restriction from local data.
     * @param local_cfl: Local CFL value (dt suggested by this process)
     * @return: Global minimum dt (most restrictive)
     */
    double reduce_cfl_timestep(double local_cfl) const;

    /**
     * Gather local statistics and compute global statistics.
     * @param local_stats: Statistics from this process
     * @return: Global statistics (min/max across all processes)
     */
    Statistics gather_statistics(const Statistics& local_stats) const;

    /**
     * Perform global synchronization barrier.
     */
    void global_barrier() const;

    /**
     * Gather all local vectors to root process.
     * @param local_vec: Local vector on this process
     * @param global_vec: [OUT] Global concatenated vector (only valid on root)
     * @return: true if this is root process and has valid data
     */
    bool gather_to_root(
        const std::vector<double>& local_vec,
        std::vector<double>& global_vec
    ) const;

    /**
     * Broadcast single double from root to all processes.
     * @param value: [IN/OUT] Value on root, receives broadcast on others
     */
    void broadcast_from_root(double& value) const;

    /**
     * Broadcast vector from root to all processes.
     * @param vec: [IN/OUT] Vector on root, receives broadcast on others
     */
    void broadcast_from_root(std::vector<double>& vec) const;

    /**
     * Check global convergence criterion across all processes.
     * @param local_residual: Residual norm on this process
     * @param tolerance: Convergence tolerance
     * @return: true if all processes satisfy residual <= tolerance
     */
    bool is_globally_converged(double local_residual, double tolerance) const;

    /**
     * Print global statistics to console (root only).
     */
    void print_statistics(const Statistics& stats) const;

private:
    /**
     * Helper for min/max reduction on complex structures.
     */
    struct MinMaxValue {
        double value;
        int rank;
    };

    /**
     * Custom MPI operation for min with rank tracking.
     */
    static void mpi_minloc(void* invec, void* inoutvec, int* len, MPI_Datatype* dtype);

    /**
     * Custom MPI operation for max with rank tracking.
     */
    static void mpi_maxloc(void* invec, void* inoutvec, int* len, MPI_Datatype* dtype);
};

} // namespace fvm3d::parallel
