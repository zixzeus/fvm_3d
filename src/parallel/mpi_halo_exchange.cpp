#include "parallel/mpi_halo_exchange.hpp"
#include "parallel/mpi_utils.hpp"
#include <algorithm>
#include <stdexcept>

namespace fvm3d::parallel {

MPIHaloExchange::MPIHaloExchange(
    const std::shared_ptr<MPIDomainDecomposer>& decomposer,
    int nghost
)
    : decomposer_(decomposer), nghost_(nghost), exchange_pending_(false) {

    if (!decomposer_) {
        throw std::invalid_argument("Domain decomposer required");
    }

    // Pre-allocate buffers for each direction
    auto local_cells = decomposer_->local_total_cells();
    int ny = local_cells[1];
    int nz = local_cells[2];

    // X-direction buffers (send/recv nghost x ny x nz * nvars)
    size_t sz_x = nghost_ * ny * nz * 8;  // 8 variables for MHD
    send_buffer_xm_.resize(sz_x);
    recv_buffer_xm_.resize(sz_x);
    send_buffer_xp_.resize(sz_x);
    recv_buffer_xp_.resize(sz_x);

    // Y-direction buffers
    int nx = local_cells[0];
    size_t sz_y = nx * nghost_ * nz * 8;
    send_buffer_ym_.resize(sz_y);
    recv_buffer_ym_.resize(sz_y);
    send_buffer_yp_.resize(sz_y);
    recv_buffer_yp_.resize(sz_y);

    // Z-direction buffers
    size_t sz_z = nx * ny * nghost_ * 8;
    send_buffer_zm_.resize(sz_z);
    recv_buffer_zm_.resize(sz_z);
    send_buffer_zp_.resize(sz_z);
    recv_buffer_zp_.resize(sz_z);
}

MPIHaloExchange::~MPIHaloExchange() {
    // Wait for pending communications
    if (exchange_pending_) {
        wait_exchange();
    }
}

void MPIHaloExchange::exchange(StateField3D& state) {
    start_exchange(state);
    wait_exchange();
}

void MPIHaloExchange::start_exchange(StateField3D& state) {
    if (exchange_pending_) {
        throw std::runtime_error("Exchange already in progress");
    }

    requests_.clear();

    // Post non-blocking exchanges for all directions
    post_exchange_direction(state, 0);  // X
    post_exchange_direction(state, 1);  // Y
    post_exchange_direction(state, 2);  // Z

    exchange_pending_ = true;
}

void MPIHaloExchange::wait_exchange() {
    if (!exchange_pending_) {
        return;
    }

    // Wait for all requests
    if (!requests_.empty()) {
        std::vector<MPI_Status> statuses(requests_.size());
        int error = MPI_Waitall(
            requests_.size(),
            requests_.data(),
            statuses.data()
        );
        MPIUtils::check_mpi_error(error, "MPI_Waitall");
    }

    exchange_pending_ = false;
}

void MPIHaloExchange::post_exchange_direction(
    StateField3D& state,
    int direction
) {
    // Pack and post minus-side (lower boundary)
    if (decomposer_->has_neighbor(direction, -1)) {
        pack_send_buffer(state, direction, -1);

        int neighbor = decomposer_->neighbor_rank(direction, -1);
        int tag = direction * 2;

        // Non-blocking send/recv pair
        MPI_Request req_send, req_recv;

        // Send interior data to neighbor
        int error = MPI_Isend(
            send_buffer_xm_.data(),
            send_buffer_xm_.size(),
            MPI_DOUBLE,
            neighbor,
            tag,
            decomposer_->cartesian_comm(),
            &req_send
        );
        MPIUtils::check_mpi_error(error, "MPI_Isend (minus)");
        requests_.push_back(req_send);

        // Receive ghost data from neighbor
        error = MPI_Irecv(
            recv_buffer_xm_.data(),
            recv_buffer_xm_.size(),
            MPI_DOUBLE,
            neighbor,
            tag + 1,
            decomposer_->cartesian_comm(),
            &req_recv
        );
        MPIUtils::check_mpi_error(error, "MPI_Irecv (minus)");
        requests_.push_back(req_recv);
    }

    // Pack and post plus-side (upper boundary)
    if (decomposer_->has_neighbor(direction, +1)) {
        pack_send_buffer(state, direction, +1);

        int neighbor = decomposer_->neighbor_rank(direction, +1);
        int tag = direction * 2 + 1;

        MPI_Request req_send, req_recv;

        // Send interior data to neighbor
        int error = MPI_Isend(
            send_buffer_xp_.data(),
            send_buffer_xp_.size(),
            MPI_DOUBLE,
            neighbor,
            tag,
            decomposer_->cartesian_comm(),
            &req_send
        );
        MPIUtils::check_mpi_error(error, "MPI_Isend (plus)");
        requests_.push_back(req_send);

        // Receive ghost data from neighbor
        error = MPI_Irecv(
            recv_buffer_xp_.data(),
            recv_buffer_xp_.size(),
            MPI_DOUBLE,
            neighbor,
            tag + 1,
            decomposer_->cartesian_comm(),
            &req_recv
        );
        MPIUtils::check_mpi_error(error, "MPI_Irecv (plus)");
        requests_.push_back(req_recv);
    }
}

void MPIHaloExchange::pack_send_buffer(
    const StateField3D& state,
    int direction,
    int side
) {
    // Extract interior cells for sending
    const int nx = state.nx();
    const int ny = state.ny();
    const int nz = state.nz();
    const int nvars = state.nvars();  // Get actual number of variables

    // Vectorized approach: maintain interleaved buffer layout but optimize reads
    // Buffer layout: [v0,v1,...,v7] for each (i,j,k) point

    if (direction == 0) {
        // X-direction: send layers from x=nghost or x=nx-2*nghost
        std::vector<double>& buffer = (side < 0) ? send_buffer_xm_ : send_buffer_xp_;
        const int i_start = (side < 0) ? nghost_ : (nx - 2*nghost_);

        size_t idx = 0;
        for (int i = i_start; i < i_start + nghost_; i++) {
            for (int j = 0; j < ny; j++) {
                // Vectorize k loop with manual gather from SoA layout
                for (int k = 0; k < nz; k++) {
                    // Pack all variables for this (i,j,k) point
                    // Use pointer arithmetic for better compiler optimization
                    double* buf_ptr = &buffer[idx];
                    #pragma omp simd
                    for (int v = 0; v < nvars; v++) {
                        buf_ptr[v] = state(v, i, j, k);
                    }
                    idx += nvars;
                }
            }
        }
    } else if (direction == 1) {
        // Y-direction
        std::vector<double>& buffer = (side < 0) ? send_buffer_ym_ : send_buffer_yp_;
        const int j_start = (side < 0) ? nghost_ : (ny - 2*nghost_);

        size_t idx = 0;
        for (int i = 0; i < nx; i++) {
            for (int j = j_start; j < j_start + nghost_; j++) {
                for (int k = 0; k < nz; k++) {
                    double* buf_ptr = &buffer[idx];
                    #pragma omp simd
                    for (int v = 0; v < nvars; v++) {
                        buf_ptr[v] = state(v, i, j, k);
                    }
                    idx += nvars;
                }
            }
        }
    } else {
        // Z-direction
        std::vector<double>& buffer = (side < 0) ? send_buffer_zm_ : send_buffer_zp_;
        const int k_start = (side < 0) ? nghost_ : (nz - 2*nghost_);

        size_t idx = 0;
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = k_start; k < k_start + nghost_; k++) {
                    double* buf_ptr = &buffer[idx];
                    #pragma omp simd
                    for (int v = 0; v < nvars; v++) {
                        buf_ptr[v] = state(v, i, j, k);
                    }
                    idx += nvars;
                }
            }
        }
    }
}

void MPIHaloExchange::unpack_recv_buffer(
    StateField3D& state,
    int direction,
    int side
) {
    // Write received ghost cells
    const int nx = state.nx();
    const int ny = state.ny();
    const int nz = state.nz();
    const int nvars = state.nvars();  // Get actual number of variables

    // Vectorized approach: maintain interleaved buffer layout but optimize writes
    // Buffer layout: [v0,v1,...,v7] for each (i,j,k) point

    if (direction == 0) {
        // X-direction: write to x=0 or x=nx-nghost
        std::vector<double>& buffer = (side < 0) ? recv_buffer_xm_ : recv_buffer_xp_;
        const int i_start = (side < 0) ? 0 : (nx - nghost_);

        size_t idx = 0;
        for (int i = i_start; i < i_start + nghost_; i++) {
            for (int j = 0; j < ny; j++) {
                // Vectorize k loop with manual scatter to SoA layout
                for (int k = 0; k < nz; k++) {
                    // Unpack all variables for this (i,j,k) point
                    const double* buf_ptr = &buffer[idx];
                    #pragma omp simd
                    for (int v = 0; v < nvars; v++) {
                        state(v, i, j, k) = buf_ptr[v];
                    }
                    idx += nvars;
                }
            }
        }
    } else if (direction == 1) {
        // Y-direction
        std::vector<double>& buffer = (side < 0) ? recv_buffer_ym_ : recv_buffer_yp_;
        const int j_start = (side < 0) ? 0 : (ny - nghost_);

        size_t idx = 0;
        for (int i = 0; i < nx; i++) {
            for (int j = j_start; j < j_start + nghost_; j++) {
                for (int k = 0; k < nz; k++) {
                    const double* buf_ptr = &buffer[idx];
                    #pragma omp simd
                    for (int v = 0; v < nvars; v++) {
                        state(v, i, j, k) = buf_ptr[v];
                    }
                    idx += nvars;
                }
            }
        }
    } else {
        // Z-direction
        std::vector<double>& buffer = (side < 0) ? recv_buffer_zm_ : recv_buffer_zp_;
        const int k_start = (side < 0) ? 0 : (nz - nghost_);

        size_t idx = 0;
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = k_start; k < k_start + nghost_; k++) {
                    const double* buf_ptr = &buffer[idx];
                    #pragma omp simd
                    for (int v = 0; v < nvars; v++) {
                        state(v, i, j, k) = buf_ptr[v];
                    }
                    idx += nvars;
                }
            }
        }
    }
}

size_t MPIHaloExchange::buffer_size(int direction) const {
    auto local = decomposer_->local_total_cells();
    int nvars = 8;  // MHD

    if (direction == 0) {
        return nghost_ * local[1] * local[2] * nvars;
    } else if (direction == 1) {
        return local[0] * nghost_ * local[2] * nvars;
    } else {
        return local[0] * local[1] * nghost_ * nvars;
    }
}

} // namespace fvm3d::parallel
