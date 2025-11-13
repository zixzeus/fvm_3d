#include "parallel/mpi_utils.hpp"
#include <iostream>
#include <stdexcept>

namespace fvm3d::parallel {

void MPIUtils::init(int& argc, char**& argv) {
    int initialized;
    MPI_Initialized(&initialized);

    if (!initialized) {
        int error = MPI_Init(&argc, &argv);
        check_mpi_error(error, "MPI_Init");
    }
}

void MPIUtils::finalize() {
    int initialized;
    MPI_Initialized(&initialized);

    if (initialized) {
        int finalized;
        MPI_Finalized(&finalized);

        if (!finalized) {
            int error = MPI_Finalize();
            check_mpi_error(error, "MPI_Finalize");
        }
    }
}

bool MPIUtils::is_initialized() {
    int initialized;
    MPI_Initialized(&initialized);
    return initialized != 0;
}

MPI_Comm MPIUtils::world_comm() {
    return MPI_COMM_WORLD;
}

int MPIUtils::rank() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int MPIUtils::size() {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

bool MPIUtils::is_root() {
    return rank() == 0;
}

void MPIUtils::print_root(const std::string& message) {
    if (is_root()) {
        std::cout << message << std::endl;
    }
}

void MPIUtils::print_all(const std::string& message) {
    std::cout << "[Rank " << rank() << "] " << message << std::endl;
}

void MPIUtils::barrier() {
    MPI_Barrier(MPI_COMM_WORLD);
}

void MPIUtils::check_mpi_error(int error_code, const std::string& context) {
    if (error_code != MPI_SUCCESS) {
        char error_string[MPI_MAX_ERROR_STRING];
        int error_len;
        MPI_Error_string(error_code, error_string, &error_len);

        std::string msg = "MPI Error in " + context + ": " +
                         std::string(error_string, error_len);
        throw std::runtime_error(msg);
    }
}

// MPIGuard implementation
MPIGuard::MPIGuard(int& argc, char**& argv) {
    MPIUtils::init(argc, argv);
}

MPIGuard::~MPIGuard() {
    MPIUtils::finalize();
}

} // namespace fvm3d::parallel
