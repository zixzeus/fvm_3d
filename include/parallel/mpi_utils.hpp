#pragma once

#include <mpi.h>
#include <vector>
#include <array>
#include <string>

namespace fvm3d::parallel {

/**
 * MPI utilities and helper functions.
 *
 * Provides wrappers around MPI functions for:
 * - Process rank and size queries
 * - Cartesian topology creation
 * - Data type registration
 * - Error checking
 */
class MPIUtils {
public:
    /**
     * Initialize MPI (call once at program start).
     * @param argc, argv: Command line arguments
     */
    static void init(int& argc, char**& argv);

    /**
     * Finalize MPI (call once at program end).
     */
    static void finalize();

    /**
     * Check if MPI is initialized.
     */
    static bool is_initialized();

    /**
     * Get global communicator (MPI_COMM_WORLD).
     */
    static MPI_Comm world_comm();

    /**
     * Get process rank.
     */
    static int rank();

    /**
     * Get total number of processes.
     */
    static int size();

    /**
     * Check if this is the root process.
     */
    static bool is_root();

    /**
     * Print message from root process only.
     */
    static void print_root(const std::string& message);

    /**
     * Print message from all processes with rank prefix.
     */
    static void print_all(const std::string& message);

    /**
     * Global barrier synchronization.
     */
    static void barrier();

    /**
     * MPI error checking helper.
     */
    static void check_mpi_error(int error_code, const std::string& context);
};

/**
 * RAII wrapper for MPI initialization/finalization.
 */
class MPIGuard {
public:
    /**
     * Constructor: Initialize MPI.
     */
    MPIGuard(int& argc, char**& argv);

    /**
     * Destructor: Finalize MPI.
     */
    ~MPIGuard();

    /**
     * Deleted copy operations.
     */
    MPIGuard(const MPIGuard&) = delete;
    MPIGuard& operator=(const MPIGuard&) = delete;
};

} // namespace fvm3d::parallel
