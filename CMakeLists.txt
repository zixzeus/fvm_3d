cmake_minimum_required(VERSION 3.10)
project(FVM3D)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# If not specified, use Release mode
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Find Eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Find HDF5 with MPI support
set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 REQUIRED COMPONENTS C)
include_directories(${HDF5_INCLUDE_DIRS})

# Find MPI
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_DIRS})

# Include directories
include_directories(include)

# Source files
set(CORE_SOURCES
    src/core/grid3d.cpp
    src/core/fvm_solver_base.cpp
    src/core/fvm_solver3d.cpp
    src/core/mpi_fvm_solver3d.cpp
)

set(PHYSICS_SOURCES
    src/physics/euler3d.cpp
    src/physics/resistive_mhd3d.cpp
    src/physics/resistive_mhd3d_advanced.cpp
)

set(SPATIAL_SOURCES
    src/spatial/flux_calculation/lax_friedrichs.cpp
    src/spatial/flux_calculation/riemann_flux.cpp
    src/spatial/riemann_solvers/riemann_hll.cpp
    src/spatial/riemann_solvers/riemann_hllc.cpp
    src/spatial/riemann_solvers/riemann_hlld.cpp
    src/spatial/riemann_solvers/riemann_solver_factory.cpp
    src/spatial/flux_calculation/flux_calculator_factory.cpp
    src/spatial/reconstruction/reconstruction_base.cpp
    src/spatial/reconstruction/constant_reconstruction.cpp
    src/spatial/reconstruction/muscl_reconstruction.cpp
    src/spatial/reconstruction/reconstruction_factory.cpp
)

set(TEMPORAL_SOURCES
    src/temporal/forward_euler.cpp
    src/temporal/rk_integrators.cpp
    src/temporal/time_integrator_factory.cpp
)

set(BOUNDARY_SOURCES
    src/boundary/periodic_bc.cpp
    src/boundary/reflective_bc.cpp
    src/boundary/transmissive_bc.cpp
)

set(IO_SOURCES
    src/io/hdf5_checkpoint.cpp
    src/io/mpi_hdf5_checkpoint.cpp
    src/io/vtk_writer.cpp
    src/io/parallel_vtk_writer.cpp
    src/io/vtk_physics_helper.cpp
)

set(PARALLEL_SOURCES
    src/parallel/mpi_utils.cpp
    src/parallel/mpi_domain_decomposer.cpp
    src/parallel/mpi_halo_exchange.cpp
    src/parallel/mpi_global_reduction.cpp
)

# Create static library
add_library(fvm3d_lib
    ${CORE_SOURCES}
    ${PHYSICS_SOURCES}
    ${SPATIAL_SOURCES}
    ${TEMPORAL_SOURCES}
    ${BOUNDARY_SOURCES}
    ${IO_SOURCES}
    ${PARALLEL_SOURCES}
)

target_link_libraries(fvm3d_lib Eigen3::Eigen ${HDF5_C_LIBRARIES} MPI::MPI_CXX)

# Example executables
add_executable(blast_wave_example examples/blast_wave_3d.cpp)
target_link_libraries(blast_wave_example fvm3d_lib Eigen3::Eigen ${HDF5_C_LIBRARIES} MPI::MPI_CXX)

add_executable(blast_wave_simulation examples/blast_wave_simulation.cpp)
target_link_libraries(blast_wave_simulation fvm3d_lib Eigen3::Eigen ${HDF5_C_LIBRARIES} MPI::MPI_CXX)

add_executable(checkpoint_restart_example examples/checkpoint_restart_example.cpp)
target_link_libraries(checkpoint_restart_example fvm3d_lib Eigen3::Eigen ${HDF5_C_LIBRARIES} MPI::MPI_CXX)

add_executable(test_hdf5_checkpoint examples/test_hdf5_checkpoint.cpp)
target_link_libraries(test_hdf5_checkpoint fvm3d_lib Eigen3::Eigen ${HDF5_C_LIBRARIES} MPI::MPI_CXX)

add_executable(harris_sheet_3d examples/harris_sheet_3d.cpp)
target_link_libraries(harris_sheet_3d fvm3d_lib Eigen3::Eigen ${HDF5_C_LIBRARIES} MPI::MPI_CXX)

# New MPI-parallel examples
add_executable(mpi_blast_wave_3d examples/mpi_blast_wave_3d.cpp)
target_link_libraries(mpi_blast_wave_3d fvm3d_lib Eigen3::Eigen ${HDF5_C_LIBRARIES} MPI::MPI_CXX)

add_executable(mpi_magnetic_reconnection examples/mpi_magnetic_reconnection.cpp)
target_link_libraries(mpi_magnetic_reconnection fvm3d_lib Eigen3::Eigen ${HDF5_C_LIBRARIES} MPI::MPI_CXX)

add_executable(vtk_visualization_example examples/vtk_visualization_example.cpp)
target_link_libraries(vtk_visualization_example fvm3d_lib Eigen3::Eigen ${HDF5_C_LIBRARIES} MPI::MPI_CXX)

add_executable(mpi_vtk_visualization_example examples/mpi_vtk_visualization_example.cpp)
target_link_libraries(mpi_vtk_visualization_example fvm3d_lib Eigen3::Eigen ${HDF5_C_LIBRARIES} MPI::MPI_CXX)

# Print build information
message(STATUS "FVM3D Build Configuration:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  MPI Support: Enabled")
