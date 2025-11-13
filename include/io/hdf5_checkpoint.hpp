#pragma once

#include "core/field3d.hpp"
#include "core/grid3d.hpp"
#include <string>
#include <memory>

namespace fvm3d::io {

using StateField3D = fvm3d::core::StateField3D;
using Grid3D = fvm3d::core::Grid3D;

/**
 * HDF5-based checkpoint/restart functionality.
 *
 * Provides save and load capabilities for FVM simulations:
 * - Saves complete solution state (all variables)
 * - Saves grid geometry and metadata
 * - Saves simulation time and step count
 * - Allows restart from checkpoint
 *
 * File format:
 * ├─ /grid
 * │  ├─ geometry (dataset with domain bounds and spacing)
 * │  └─ nghost (attribute)
 * ├─ /state
 * │  ├─ rho (dataset: nx × ny × nz)
 * │  ├─ rho_u (dataset)
 * │  ├─ rho_v (dataset)
 * │  ├─ rho_w (dataset)
 * │  └─ E (dataset)
 * └─ /metadata
 *    ├─ time (attribute)
 *    ├─ step_count (attribute)
 *    └─ description (attribute)
 */
class HDF5Checkpoint {
public:
    /**
     * Constructor with filename.
     * @param filename: path to HDF5 file (will be created or overwritten on save)
     */
    HDF5Checkpoint(const std::string& filename);

    ~HDF5Checkpoint() = default;

    /**
     * Save checkpoint to HDF5 file.
     * @param state: solution state field
     * @param grid: grid geometry information
     * @param time: current simulation time
     * @param step_count: current step number
     * @param description: optional description string
     */
    void save(
        const StateField3D& state,
        const Grid3D& grid,
        double time,
        int step_count,
        const std::string& description = ""
    );

    /**
     * Load checkpoint from HDF5 file.
     * @param state: solution state field (modified in-place)
     * @param time: current simulation time (output)
     * @param step_count: current step number (output)
     * @return true if load successful, false otherwise
     */
    bool load(
        StateField3D& state,
        double& time,
        int& step_count
    );

    /**
     * Get the filename.
     */
    const std::string& filename() const { return filename_; }

    /**
     * Check if file exists and is readable.
     */
    static bool file_exists(const std::string& filename);

    /**
     * Get metadata from checkpoint file without loading full state.
     * @param filename: HDF5 file to read
     * @param time: current simulation time (output)
     * @param step_count: current step number (output)
     * @param description: description string (output)
     * @return true if metadata read successfully, false otherwise
     */
    static bool read_metadata(
        const std::string& filename,
        double& time,
        int& step_count,
        std::string& description
    );

private:
    std::string filename_;

    /**
     * Save grid information to HDF5 group.
     */
    void save_grid(void* file_id, const Grid3D& grid);

    /**
     * Save state variables to HDF5 group.
     */
    void save_state(void* file_id, const StateField3D& state);

    /**
     * Save metadata to HDF5 group.
     */
    void save_metadata(
        void* file_id,
        double time,
        int step_count,
        const std::string& description
    );

    /**
     * Load grid information from HDF5 group.
     * @param file_id: HDF5 file ID
     * @return true if successful
     */
    bool load_grid(void* file_id);

    /**
     * Load state variables from HDF5 group.
     * @param file_id: HDF5 file ID
     * @param state: state field to populate
     * @return true if successful
     */
    bool load_state(void* file_id, StateField3D& state);

    /**
     * Load metadata from HDF5 group.
     * @param file_id: HDF5 file ID
     * @param time: simulation time (output)
     * @param step_count: step count (output)
     * @param description: description (output)
     * @return true if successful
     */
    bool load_metadata(
        void* file_id,
        double& time,
        int& step_count,
        std::string& description
    );
};

} // namespace fvm3d::io
