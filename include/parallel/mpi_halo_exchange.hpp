#pragma once

#include "core/field3d.hpp"
#include "parallel/mpi_domain_decomposer.hpp"
#include <mpi.h>
#include <vector>
#include <memory>

namespace fvm3d::parallel {

/**
 * Halo (ghost cell) exchange for distributed 3D data.
 *
 * Manages non-blocking MPI communication of boundary data between
 * neighboring subdomains. Enables stencil-based computations on
 * distributed arrays while maintaining data consistency.
 *
 * Data layout: [variables, x, y, z] with ghost cells on boundaries
 */
class MPIHaloExchange {
public:
    using StateField3D = core::StateField3D;

    /**
     * Constructor.
     * @param decomposer: Domain decomposition information
     * @param nghost: Width of ghost cell layer (typically 2)
     */
    explicit MPIHaloExchange(
        const std::shared_ptr<MPIDomainDecomposer>& decomposer,
        int nghost = 2
    );

    /**
     * Destructor: Free allocated buffers.
     */
    ~MPIHaloExchange();

    /**
     * Perform halo exchange (blocking).
     * Exchanges all ghost cell values with neighboring processes.
     */
    void exchange(StateField3D& state);

    /**
     * Start non-blocking halo exchange.
     * Must be paired with wait_exchange().
     */
    void start_exchange(StateField3D& state);

    /**
     * Wait for pending non-blocking exchange to complete.
     */
    void wait_exchange();

    /**
     * Check if exchange is in progress.
     */
    bool is_exchange_pending() const { return exchange_pending_; }

    /**
     * Get number of ghost cell layers.
     */
    int nghost() const { return nghost_; }

private:
    std::shared_ptr<MPIDomainDecomposer> decomposer_;
    int nghost_;
    bool exchange_pending_;

    // State pointer for async exchange (valid only during pending exchange)
    StateField3D* exchange_state_;

    // MPI request handles for non-blocking communication
    std::vector<MPI_Request> requests_;

    // Send/receive buffers
    std::vector<double> send_buffer_xm_, recv_buffer_xm_;  // X-minus
    std::vector<double> send_buffer_xp_, recv_buffer_xp_;  // X-plus
    std::vector<double> send_buffer_ym_, recv_buffer_ym_;  // Y-minus
    std::vector<double> send_buffer_yp_, recv_buffer_yp_;  // Y-plus
    std::vector<double> send_buffer_zm_, recv_buffer_zm_;  // Z-minus
    std::vector<double> send_buffer_zp_, recv_buffer_zp_;  // Z-plus

    /**
     * Pack data from interior region into send buffer.
     */
    void pack_send_buffer(
        const StateField3D& state,
        int direction,  // 0=X, 1=Y, 2=Z
        int side        // -1=minus, +1=plus
    );

    /**
     * Unpack received data into ghost cell region.
     */
    void unpack_recv_buffer(
        StateField3D& state,
        int direction,  // 0=X, 1=Y, 2=Z
        int side        // -1=minus, +1=plus
    );

    /**
     * Post non-blocking send/receive pairs for a direction.
     */
    void post_exchange_direction(
        StateField3D& state,
        int direction   // 0=X, 1=Y, 2=Z
    );

    /**
     * Calculate buffer size needed for a halo layer.
     */
    size_t buffer_size(int direction) const;
};

} // namespace fvm3d::parallel
