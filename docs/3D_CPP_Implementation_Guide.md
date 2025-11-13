# 3D C++ MPI框架实现指南

**文档类型**: 开发者指南
**目标读者**: C++开发工程师
**难度级别**: 中等

---

## 目录

1. [环境配置](#环境配置)
2. [项目初始化](#项目初始化)
3. [核心模块实现](#核心模块实现)
4. [MPI并行实现](#mpi并行实现)
5. [常见问题与调试](#常见问题与调试)
6. [性能优化检查清单](#性能优化检查清单)

---

## 环境配置

### 依赖安装（Ubuntu 20.04/22.04）

```bash
# 更新包管理器
sudo apt-get update

# 1. GCC/G++ 11+
sudo apt-get install gcc-11 g++-11
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 100
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 100

# 2. CMake 3.20+
sudo apt-get install cmake

# 3. MPI (OpenMPI)
sudo apt-get install libopenmpi-dev openmpi-bin

# 4. Eigen3
sudo apt-get install libeigen3-dev

# 5. HDF5（带MPI支持）
sudo apt-get install libhdf5-openmpi-dev

# 6. Git
sudo apt-get install git

# 7. 可选：Catch2（单元测试）
sudo apt-get install catch2

# 8. 可选：ccache（加快编译）
sudo apt-get install ccache
```

### 验证安装

```bash
# 检查编译器
g++ --version  # >= 11.0

# 检查MPI
mpicc --version
mpirun --version

# 检查CMake
cmake --version  # >= 3.20

# 检查Eigen3
pkg-config --cflags --libs eigen3

# 检查HDF5
h5cc -show
```

---

## 项目初始化

### 克隆与初始设置

```bash
# 创建项目目录
mkdir -p fvm3d_project
cd fvm3d_project

# 初始化Git
git init
git config user.name "Your Name"
git config user.email "your.email@example.com"

# 创建目录结构
mkdir -p include/fvm3d/{core,parallel,physics,spatial,temporal,boundary,pipeline,solver,utils}
mkdir -p src/{core,parallel,physics,spatial,temporal,boundary,pipeline,solver,utils}
mkdir -p tests
mkdir -p examples
mkdir -p benchmarks
mkdir -p docs
mkdir -p scripts
mkdir -p build

# 复制现有Python框架文件供参考
cp -r /path/to/fvm_framework fvm_framework_python_ref
```

### 根目录CMakeLists.txt（完整版）

```cmake
# CMakeLists.txt
cmake_minimum_required(VERSION 3.20)
project(FVM3D
    VERSION 1.0.0
    DESCRIPTION "3D Finite Volume Method Framework with MPI"
    LANGUAGES CXX)

# ============================================================================
# 编译选项
# ============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 默认编译类型
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# 优化标志
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG -ffast-math")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -DDEBUG")

# 启用ccache加速
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

# ============================================================================
# 依赖查找
# ============================================================================
find_package(MPI REQUIRED COMPONENTS CXX)
find_package(Eigen3 REQUIRED CONFIG)
find_package(HDF5 REQUIRED COMPONENTS CXX HL)
find_package(OpenMP)  # 可选

# Catch2测试框架
include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.4.0
)
FetchContent_MakeAvailable(Catch2)

# ============================================================================
# 头文件包含
# ============================================================================
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${MPI_CXX_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

# ============================================================================
# 库定义
# ============================================================================
add_library(fvm3d SHARED)

# 添加源文件
file(GLOB_RECURSE FVM3D_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/**/*.cpp")
target_sources(fvm3d PRIVATE ${FVM3D_SOURCES})

# 链接依赖
target_link_libraries(fvm3d
    PUBLIC
        MPI::MPI_CXX
        Eigen3::Eigen
        ${HDF5_LIBRARIES}
        ${HDF5_HL_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(fvm3d PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(fvm3d PUBLIC FVM3D_OPENMP)
endif()

# ============================================================================
# 编译信息输出
# ============================================================================
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "MPI: ${MPI_CXX_LIBRARIES}")
message(STATUS "Eigen3: ${EIGEN3_INCLUDE_DIRS}")
message(STATUS "HDF5: ${HDF5_LIBRARIES}")

# ============================================================================
# 测试
# ============================================================================
enable_testing()
add_subdirectory(tests)

# ============================================================================
# 示例程序
# ============================================================================
add_subdirectory(examples)

# ============================================================================
# 安装
# ============================================================================
install(TARGETS fvm3d
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
        DESTINATION include)
```

---

## 核心模块实现

### 1. Grid3D实现

#### 头文件（include/fvm3d/core/grid3d.hpp）

```cpp
#ifndef FVM3D_CORE_GRID3D_HPP
#define FVM3D_CORE_GRID3D_HPP

namespace fvm3d {

/**
 * @class GridGeometry3D
 * @brief 3D网格几何信息
 */
struct GridGeometry3D {
    double Lx, Ly, Lz;      ///< 域尺寸
    int nx, ny, nz;         ///< 单元数
    double dx, dy, dz;      ///< 单元尺寸

    GridGeometry3D() = default;

    GridGeometry3D(double lx, double ly, double lz,
                   int nx_, int ny_, int nz_)
        : Lx(lx), Ly(ly), Lz(lz), nx(nx_), ny(ny_), nz(nz_)
        , dx(lx/nx), dy(ly/ny), dz(lz/nz) {
        if (dx <= 0 || dy <= 0 || dz <= 0) {
            throw std::invalid_argument("Grid spacing must be positive");
        }
    }

    // 单元中心坐标
    double cell_center_x(int i) const { return (i + 0.5) * dx; }
    double cell_center_y(int j) const { return (j + 0.5) * dy; }
    double cell_center_z(int k) const { return (k + 0.5) * dz; }

    // 单元体积
    double cell_volume() const { return dx * dy * dz; }
};

/**
 * @class Grid3D
 * @brief 3D均匀笛卡尔网格
 */
class Grid3D {
public:
    /// 构造函数
    /// @param geom 网格几何
    /// @param nghost 幽灵层数
    Grid3D(const GridGeometry3D& geom, int nghost = 2)
        : geom_(geom), nghost_(nghost)
        , nx_total_(geom.nx + 2*nghost)
        , ny_total_(geom.ny + 2*nghost)
        , nz_total_(geom.nz + 2*nghost) {
        if (nghost < 0) {
            throw std::invalid_argument("nghost must be non-negative");
        }
    }

    // 访问器
    int nx() const { return geom_.nx; }
    int ny() const { return geom_.ny; }
    int nz() const { return geom_.nz; }
    int nghost() const { return nghost_; }

    int nx_total() const { return nx_total_; }
    int ny_total() const { return ny_total_; }
    int nz_total() const { return nz_total_; }

    // 几何信息
    double Lx() const { return geom_.Lx; }
    double Ly() const { return geom_.Ly; }
    double Lz() const { return geom_.Lz; }
    double dx() const { return geom_.dx; }
    double dy() const { return geom_.dy; }
    double dz() const { return geom_.dz; }
    double volume() const { return geom_.cell_volume(); }

    // 单元中心坐标
    double x(int i) const { return geom_.cell_center_x(i - nghost_); }
    double y(int j) const { return geom_.cell_center_y(j - nghost_); }
    double z(int k) const { return geom_.cell_center_z(k - nghost_); }

    // 索引范围
    int i_min() const { return nghost_; }
    int i_max() const { return nx_total_ - nghost_; }
    int j_min() const { return nghost_; }
    int j_max() const { return ny_total_ - nghost_; }
    int k_min() const { return nghost_; }
    int k_max() const { return nz_total_ - nghost_; }

    // 检查索引是否在幽灵层内
    bool is_ghost_x(int i) const {
        return i < nghost_ || i >= nx_total_ - nghost_;
    }
    bool is_ghost_y(int j) const {
        return j < nghost_ || j >= ny_total_ - nghost_;
    }
    bool is_ghost_z(int k) const {
        return k < nghost_ || k >= nz_total_ - nghost_;
    }

private:
    GridGeometry3D geom_;
    int nghost_;
    int nx_total_, ny_total_, nz_total_;
};

}  // namespace fvm3d

#endif  // FVM3D_CORE_GRID3D_HPP
```

---

#### 实现文件（src/core/grid3d.cpp）

```cpp
// src/core/grid3d.cpp
#include <fvm3d/core/grid3d.hpp>
#include <iostream>

namespace fvm3d {

// Grid3D的所有方法都内联，无需额外实现
// 但可以添加序列化等高级功能

std::ostream& operator<<(std::ostream& os, const GridGeometry3D& geom) {
    os << "GridGeometry3D: "
       << "Lx=" << geom.Lx << ", Ly=" << geom.Ly << ", Lz=" << geom.Lz << " "
       << "nx=" << geom.nx << ", ny=" << geom.ny << ", nz=" << geom.nz << " "
       << "dx=" << geom.dx << ", dy=" << geom.dy << ", dz=" << geom.dz;
    return os;
}

}  // namespace fvm3d
```

---

### 2. Field3D实现

#### 头文件（include/fvm3d/core/field3d.hpp）

```cpp
#ifndef FVM3D_CORE_FIELD3D_HPP
#define FVM3D_CORE_FIELD3D_HPP

#include <memory>
#include <vector>
#include <algorithm>
#include <cstring>

namespace fvm3d {

/**
 * @class Field3D
 * @brief 3D多维数组容器（SoA布局）
 * @tparam T 数据类型
 */
template<typename T>
class Field3D {
public:
    /**
     * 构造函数
     * @param nvars 变量数
     * @param nx X方向单元数
     * @param ny Y方向单元数
     * @param nz Z方向单元数
     */
    Field3D(int nvars, int nx, int ny, int nz)
        : nvars_(nvars), nx_(nx), ny_(ny), nz_(nz)
        , size_(static_cast<size_t>(nvars) * nx * ny * nz)
        , data_(new T[size_]) {
        std::fill(data_.get(), data_.get() + size_, T(0));
    }

    // 移动构造函数
    Field3D(Field3D&& other) noexcept
        : nvars_(other.nvars_), nx_(other.nx_), ny_(other.ny_), nz_(other.nz_)
        , size_(other.size_), data_(std::move(other.data_)) {
        other.nvars_ = other.nx_ = other.ny_ = other.nz_ = 0;
        other.size_ = 0;
    }

    // 移动赋值
    Field3D& operator=(Field3D&& other) noexcept {
        if (this != &other) {
            nvars_ = other.nvars_;
            nx_ = other.nx_;
            ny_ = other.ny_;
            nz_ = other.nz_;
            size_ = other.size_;
            data_ = std::move(other.data_);

            other.nvars_ = other.nx_ = other.ny_ = other.nz_ = 0;
            other.size_ = 0;
        }
        return *this;
    }

    // 禁用复制（使用move代替）
    Field3D(const Field3D&) = delete;
    Field3D& operator=(const Field3D&) = delete;

    /**
     * 访问操作符：data(var, i, j, k)
     * @param var 变量索引（0到nvars-1）
     * @param i X方向索引
     * @param j Y方向索引
     * @param k Z方向索引
     * @return 引用
     */
    T& operator()(int var, int i, int j, int k) {
        return data_[index(var, i, j, k)];
    }

    const T& operator()(int var, int i, int j, int k) const {
        return data_[index(var, i, j, k)];
    }

    /**
     * 原始指针访问（用于MPI或C库）
     */
    T* data() { return data_.get(); }
    const T* data() const { return data_.get(); }

    /**
     * 尺寸信息
     */
    int nvars() const { return nvars_; }
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }
    size_t size() const { return size_; }
    size_t bytes() const { return size_ * sizeof(T); }

    /**
     * 填充数据
     */
    void fill(T value) {
        std::fill(data_.get(), data_.get() + size_, value);
    }

    /**
     * 获取某个变量的2D切片（固定k值）
     */
    std::vector<T> get_slice_xy(int var, int k) const {
        std::vector<T> slice(nx_ * ny_);
        for (int i = 0; i < nx_; ++i) {
            for (int j = 0; j < ny_; ++j) {
                slice[i*ny_ + j] = operator()(var, i, j, k);
            }
        }
        return slice;
    }

private:
    int nvars_, nx_, ny_, nz_;
    size_t size_;
    std::unique_ptr<T[]> data_;

    /**
     * 内存布局：行优先（C风格）
     * index = var*(nx*ny*nz) + i*(ny*nz) + j*nz + k
     */
    inline size_t index(int var, int i, int j, int k) const {
        return static_cast<size_t>(var) * (nx_ * ny_ * nz_)
             + static_cast<size_t>(i) * (ny_ * nz_)
             + static_cast<size_t>(j) * nz_
             + static_cast<size_t>(k);
    }
};

// 常见别名
using StateField3D = Field3D<double>;
using IntField3D = Field3D<int>;

}  // namespace fvm3d

#endif  // FVM3D_CORE_FIELD3D_HPP
```

---

### 3. 物理方程实现（Euler3D）

#### 头文件（include/fvm3d/physics/euler3d.hpp）

```cpp
#ifndef FVM3D_PHYSICS_EULER3D_HPP
#define FVM3D_PHYSICS_EULER3D_HPP

#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <string>

namespace fvm3d {

/**
 * @class EulerEquations3D
 * @brief 三维可压缩Euler方程
 *
 * 守恒变量: U = [ρ, ρu, ρv, ρw, E]
 * 原始变量: V = [ρ, u, v, w, p]
 */
class EulerEquations3D {
public:
    /// 构造函数
    /// @param gamma 比热比 (默认 1.4)
    explicit EulerEquations3D(double gamma = 1.4)
        : gamma_(gamma) {
        if (gamma <= 1.0) {
            throw std::invalid_argument("gamma must be > 1.0");
        }
    }

    // 物理信息
    static constexpr int num_vars() { return 5; }

    std::vector<std::string> var_names() const {
        return {"rho", "rho_u", "rho_v", "rho_w", "E"};
    }

    // 属性
    double gamma() const { return gamma_; }

    /**
     * 守恒变量 → 原始变量
     */
    void conservative_to_primitive(
        const Eigen::VectorXd& U,
        Eigen::VectorXd& V) const {

        double rho = U(0);
        double rho_u = U(1);
        double rho_v = U(2);
        double rho_w = U(3);
        double E = U(4);

        // 避免rho过小导致数值不稳定
        constexpr double rho_floor = 1e-10;
        rho = std::max(rho, rho_floor);

        // 原始变量
        double u = rho_u / rho;
        double v = rho_v / rho;
        double w = rho_w / rho;

        // 动能（稳定计算方式）
        double ke = 0.5 * (rho_u*rho_u + rho_v*rho_v + rho_w*rho_w) / rho;

        // 内能和压强
        double e_internal = E - ke;
        double p = (gamma_ - 1.0) * e_internal;
        p = std::max(p, 1e-11);  // 压强下限

        V.resize(5);
        V << rho, u, v, w, p;
    }

    /**
     * 原始变量 → 守恒变量
     */
    void primitive_to_conservative(
        const Eigen::VectorXd& V,
        Eigen::VectorXd& U) const {

        double rho = V(0);
        double u = V(1);
        double v = V(2);
        double w = V(3);
        double p = V(4);

        U.resize(5);
        U(0) = rho;
        U(1) = rho * u;
        U(2) = rho * v;
        U(3) = rho * w;

        // 总能量
        double ke = 0.5 * rho * (u*u + v*v + w*w);
        double e_internal = p / (gamma_ - 1.0);
        U(4) = e_internal + ke;
    }

    /**
     * X方向通量 F(U)
     */
    void compute_flux_x(
        const Eigen::VectorXd& U,
        Eigen::VectorXd& F) const {

        Eigen::VectorXd V;
        conservative_to_primitive(U, V);

        double rho = V(0);
        double u = V(1);
        double v = V(2);
        double w = V(3);
        double p = V(4);
        double E = U(4);

        F.resize(5);
        F(0) = rho * u;                    // 质量
        F(1) = rho*u*u + p;                // x动量
        F(2) = rho*u*v;                    // y动量
        F(3) = rho*u*w;                    // z动量
        F(4) = u * (E + p);                // 能量
    }

    /**
     * Y方向通量 G(U)
     */
    void compute_flux_y(
        const Eigen::VectorXd& U,
        Eigen::VectorXd& G) const {

        Eigen::VectorXd V;
        conservative_to_primitive(U, V);

        double rho = V(0);
        double u = V(1);
        double v = V(2);
        double w = V(3);
        double p = V(4);
        double E = U(4);

        G.resize(5);
        G(0) = rho * v;                    // 质量
        G(1) = rho*v*u;                    // x动量
        G(2) = rho*v*v + p;                // y动量
        G(3) = rho*v*w;                    // z动量
        G(4) = v * (E + p);                // 能量
    }

    /**
     * Z方向通量 H(U)
     */
    void compute_flux_z(
        const Eigen::VectorXd& U,
        Eigen::VectorXd& H) const {

        Eigen::VectorXd V;
        conservative_to_primitive(U, V);

        double rho = V(0);
        double u = V(1);
        double v = V(2);
        double w = V(3);
        double p = V(4);
        double E = U(4);

        H.resize(5);
        H(0) = rho * w;                    // 质量
        H(1) = rho*w*u;                    // x动量
        H(2) = rho*w*v;                    // y动量
        H(3) = rho*w*w + p;                // z动量
        H(4) = w * (E + p);                // 能量
    }

    /**
     * 最大波速（用于CFL）
     * @param U 守恒变量
     * @param direction 方向 (0=x, 1=y, 2=z)
     * @return 该方向的最大波速
     */
    double max_wave_speed(
        const Eigen::VectorXd& U,
        int direction) const {

        Eigen::VectorXd V;
        conservative_to_primitive(U, V);

        double rho = V(0);
        double u = V(1);
        double v = V(2);
        double w = V(3);
        double p = V(4);

        // 声速
        double c = std::sqrt(gamma_ * p / rho);

        // 速度分量
        double vel_component = 0.0;
        if (direction == 0) {
            vel_component = std::abs(u);
        } else if (direction == 1) {
            vel_component = std::abs(v);
        } else {  // direction == 2
            vel_component = std::abs(w);
        }

        return vel_component + c;
    }

private:
    double gamma_;
};

}  // namespace fvm3d

#endif  // FVM3D_PHYSICS_EULER3D_HPP
```

---

## MPI并行实现

### MPIGrid3D

#### 头文件（include/fvm3d/parallel/mpi_grid3d.hpp）

```cpp
#ifndef FVM3D_PARALLEL_MPI_GRID3D_HPP
#define FVM3D_PARALLEL_MPI_GRID3D_HPP

#include <mpi.h>
#include <fvm3d/core/grid3d.hpp>
#include <stdexcept>

namespace fvm3d {

/**
 * @class MPIGrid3D
 * @brief MPI分布式3D网格（1D分解）
 *
 * 沿X方向进行区域分解
 * 每个MPI进程负责 [x_start, x_end) × [0, Ly) × [0, Lz)
 */
class MPIGrid3D {
public:
    /**
     * 构造函数
     * @param comm MPI通信器
     * @param global_geom 全局网格几何
     * @param nghost 幽灵层数
     */
    MPIGrid3D(MPI_Comm comm, const GridGeometry3D& global_geom, int nghost = 2);

    ~MPIGrid3D();

    // 全局信息
    int global_nx() const { return global_geom_.nx; }
    int global_ny() const { return global_geom_.ny; }
    int global_nz() const { return global_geom_.nz; }

    // 本地信息
    int local_nx() const { return local_nx_; }
    int local_ny() const { return local_ny_; }
    int local_nz() const { return local_nz_; }

    // 总数（包含幽灵层）
    int nx_total() const { return local_nx_ + 2*nghost_; }
    int ny_total() const { return local_ny_ + 2*nghost_; }
    int nz_total() const { return local_nz_ + 2*nghost_; }

    int nghost() const { return nghost_; }

    // 全局索引偏移
    int global_i_start() const { return i_start_; }
    int global_j_start() const { return 0; }
    int global_k_start() const { return 0; }

    // MPI信息
    MPI_Comm comm() const { return comm_; }
    int rank() const { return rank_; }
    int size() const { return size_; }

    // 邻居
    int neighbor_xm() const { return neighbor_xm_; }
    int neighbor_xp() const { return neighbor_xp_; }

    // 几何信息
    double Lx() const { return global_geom_.Lx; }
    double Ly() const { return global_geom_.Ly; }
    double Lz() const { return global_geom_.Lz; }
    double dx() const { return global_geom_.dx; }
    double dy() const { return global_geom_.dy; }
    double dz() const { return global_geom_.dz; }

    // 本地网格对象
    Grid3D local_grid() const {
        GridGeometry3D local_geom(
            global_geom_.Lx * local_nx_ / global_geom_.nx,
            global_geom_.Ly,
            global_geom_.Lz,
            local_nx_,
            local_ny_,
            local_nz_
        );
        return Grid3D(local_geom, nghost_);
    }

private:
    MPI_Comm comm_;
    int rank_, size_;

    GridGeometry3D global_geom_;
    int nghost_;

    // 本地域大小
    int local_nx_, local_ny_, local_nz_;

    // 全局索引起点
    int i_start_;

    // 邻居rank
    int neighbor_xm_, neighbor_xp_;
};

}  // namespace fvm3d

#endif  // FVM3D_PARALLEL_MPI_GRID3D_HPP
```

#### 实现文件（src/parallel/mpi_grid3d.cpp）

```cpp
// src/parallel/mpi_grid3d.cpp
#include <fvm3d/parallel/mpi_grid3d.hpp>
#include <iostream>

namespace fvm3d {

MPIGrid3D::MPIGrid3D(MPI_Comm comm, const GridGeometry3D& global_geom, int nghost)
    : comm_(comm), nghost_(nghost), global_geom_(global_geom) {

    // 获取MPI信息
    MPI_Comm_rank(comm, &rank_);
    MPI_Comm_size(comm, &size_);

    // 简单1D分解：沿X方向均匀分割
    local_ny_ = global_geom.ny;
    local_nz_ = global_geom.nz;

    if (size_ > global_geom.nx) {
        throw std::runtime_error(
            "Number of MPI processes exceeds number of cells in X direction"
        );
    }

    // 计算本进程的本地X范围
    int cells_per_rank = global_geom.nx / size_;
    int remainder = global_geom.nx % size_;

    if (rank_ < remainder) {
        local_nx_ = cells_per_rank + 1;
        i_start_ = rank_ * (cells_per_rank + 1);
    } else {
        local_nx_ = cells_per_rank;
        i_start_ = remainder * (cells_per_rank + 1) + (rank_ - remainder) * cells_per_rank;
    }

    // 邻居信息
    neighbor_xm_ = (rank_ > 0) ? rank_ - 1 : MPI_PROC_NULL;
    neighbor_xp_ = (rank_ < size_ - 1) ? rank_ + 1 : MPI_PROC_NULL;

    if (rank_ == 0) {
        std::cout << "[MPI] Domain decomposition (1D):\n";
        std::cout << "  Global grid: " << global_geom.nx << " × "
                  << global_geom.ny << " × " << global_geom.nz << "\n";
        std::cout << "  MPI processes: " << size_ << "\n";
        std::cout << "  Local grid per rank: " << local_nx_ << " × "
                  << local_ny_ << " × " << local_nz_ << "\n";
    }
}

MPIGrid3D::~MPIGrid3D() {
    // MPI_Comm 不需要主动释放（除非是自己创建的）
}

}  // namespace fvm3d
```

---

### Halo交换实现

#### 头文件（include/fvm3d/parallel/halo_exchange.hpp）

```cpp
#ifndef FVM3D_PARALLEL_HALO_EXCHANGE_HPP
#define FVM3D_PARALLEL_HALO_EXCHANGE_HPP

#include <mpi.h>
#include <fvm3d/core/field3d.hpp>
#include <fvm3d/parallel/mpi_grid3d.hpp>
#include <vector>

namespace fvm3d {

/**
 * @class MPIHaloExchange
 * @brief MPI幽灵单元交换
 *
 * 处理MPI进程间的幽灵层通信
 */
class MPIHaloExchange {
public:
    explicit MPIHaloExchange(const MPIGrid3D& grid);

    /**
     * 交换幽灵层数据
     * @param field 需要交换的场
     */
    void exchange(StateField3D& field);

private:
    const MPIGrid3D& grid_;
    std::vector<MPI_Request> requests_;

    // 缓冲区（避免重复分配）
    std::vector<double> send_buffer_xm_, recv_buffer_xm_;
    std::vector<double> send_buffer_xp_, recv_buffer_xp_;

    void allocate_buffers(const StateField3D& field);
    void pack_send_buffers(const StateField3D& field);
    void unpack_recv_buffers(StateField3D& field);
};

}  // namespace fvm3d

#endif  // FVM3D_PARALLEL_HALO_EXCHANGE_HPP
```

#### 实现文件（src/parallel/halo_exchange.cpp）

```cpp
// src/parallel/halo_exchange.cpp
#include <fvm3d/parallel/halo_exchange.hpp>
#include <algorithm>

namespace fvm3d {

MPIHaloExchange::MPIHaloExchange(const MPIGrid3D& grid)
    : grid_(grid) {
}

void MPIHaloExchange::allocate_buffers(const StateField3D& field) {
    // 缓冲大小：所有变量 × ny × nz × nghost层
    int buffer_size = field.nvars() * grid_.local_ny() * grid_.local_nz() * grid_.nghost();

    send_buffer_xm_.resize(buffer_size);
    recv_buffer_xm_.resize(buffer_size);
    send_buffer_xp_.resize(buffer_size);
    recv_buffer_xp_.resize(buffer_size);
}

void MPIHaloExchange::pack_send_buffers(const StateField3D& field) {
    int ng = grid_.nghost();
    int nx = grid_.local_nx();
    int ny = grid_.local_ny();
    int nz = grid_.local_nz();
    int nvars = field.nvars();

    // 打包左边界（i = ng to ng+1）发送给左邻居
    int idx = 0;
    for (int var = 0; var < nvars; ++var) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                for (int ig = 0; ig < ng; ++ig) {
                    send_buffer_xm_[idx++] = field(var, ng+ig, ng+j, ng+k);
                }
            }
        }
    }

    // 打包右边界（i = ng+nx-ng to ng+nx）发送给右邻居
    idx = 0;
    for (int var = 0; var < nvars; ++var) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                for (int ig = 0; ig < ng; ++ig) {
                    send_buffer_xp_[idx++] = field(var, ng+nx-ng+ig, ng+j, ng+k);
                }
            }
        }
    }
}

void MPIHaloExchange::exchange(StateField3D& field) {
    allocate_buffers(field);
    pack_send_buffers(field);

    int ng = grid_.nghost();
    int buffer_size = send_buffer_xm_.size();

    requests_.clear();

    // 非阻塞通信
    // 发送/接收左边界
    if (grid_.neighbor_xm() != MPI_PROC_NULL) {
        MPI_Request req;
        MPI_Isend(send_buffer_xm_.data(), buffer_size, MPI_DOUBLE,
                  grid_.neighbor_xm(), 0, grid_.comm(), &req);
        requests_.push_back(req);

        MPI_Irecv(recv_buffer_xm_.data(), buffer_size, MPI_DOUBLE,
                  grid_.neighbor_xm(), 0, grid_.comm(), &req);
        requests_.push_back(req);
    }

    // 发送/接收右边界
    if (grid_.neighbor_xp() != MPI_PROC_NULL) {
        MPI_Request req;
        MPI_Isend(send_buffer_xp_.data(), buffer_size, MPI_DOUBLE,
                  grid_.neighbor_xp(), 1, grid_.comm(), &req);
        requests_.push_back(req);

        MPI_Irecv(recv_buffer_xp_.data(), buffer_size, MPI_DOUBLE,
                  grid_.neighbor_xp(), 1, grid_.comm(), &req);
        requests_.push_back(req);
    }

    // 等待所有通信完成
    if (!requests_.empty()) {
        MPI_Waitall(requests_.size(), requests_.data(), MPI_STATUSES_IGNORE);
    }

    // 解包接收的数据
    unpack_recv_buffers(field);
}

void MPIHaloExchange::unpack_recv_buffers(StateField3D& field) {
    int ng = grid_.nghost();
    int ny = grid_.local_ny();
    int nz = grid_.local_nz();
    int nvars = field.nvars();

    // 解包左ghost
    if (grid_.neighbor_xm() != MPI_PROC_NULL) {
        int idx = 0;
        for (int var = 0; var < nvars; ++var) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    for (int ig = 0; ig < ng; ++ig) {
                        field(var, ig, ng+j, ng+k) = recv_buffer_xm_[idx++];
                    }
                }
            }
        }
    }

    // 解包右ghost
    if (grid_.neighbor_xp() != MPI_PROC_NULL) {
        int idx = 0;
        int nx_total = grid_.nx_total();
        for (int var = 0; var < nvars; ++var) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    for (int ig = 0; ig < ng; ++ig) {
                        field(var, nx_total-ng+ig, ng+j, ng+k) = recv_buffer_xp_[idx++];
                    }
                }
            }
        }
    }
}

}  // namespace fvm3d
```

---

## 常见问题与调试

### 问题1: MPI通信死锁

**症状**: 程序挂起，无错误信息

**原因**:
- 发送/接收顺序不匹配
- 缓冲区大小不一致

**调试步骤**:
```cpp
// 添加调试输出
MPI_Barrier(MPI_COMM_WORLD);  // 同步所有进程

if (rank == 0) {
    std::cout << "Rank " << rank << ": 进入Halo交换\n";
}

halo_exchange.exchange(field);

MPI_Barrier(MPI_COMM_WORLD);
if (rank == 0) {
    std::cout << "Rank " << rank << ": Halo交换完成\n";
}
```

### 问题2: 数值不稳定（NaN/Inf）

**症状**: 几步后出现NaN或Inf

**原因**:
- 密度或压强过小（除数过小）
- CFL条件违反
- 数值问题在某个进程的边界处放大

**修复**:
```cpp
// 在每一步后进行稳健性检查
void check_solution_validity(const StateField3D& field, int rank) {
    for (int var = 0; var < field.nvars(); ++var) {
        for (int i = 0; i < field.nx(); ++i) {
            for (int j = 0; j < field.ny(); ++j) {
                for (int k = 0; k < field.nz(); ++k) {
                    double val = field(var, i, j, k);
                    if (!std::isfinite(val)) {
                        std::cerr << "NaN detected at rank " << rank
                                  << " var=" << var << " (" << i << "," << j << "," << k << ")\n";
                        std::cerr << "Value: " << val << "\n";
                    }
                }
            }
        }
    }
}
```

---

### 问题3: 性能瓶颈（MPI通信 vs 计算）

**诊断**:
```bash
# 使用MPI profiler
mpirun -x MPICH_INSTRUMENTATION=1 -np 4 ./fvm3d 2>&1 | grep "MPI time"

# 或使用Tau/Scalasca
tau_exec mpirun -np 4 ./fvm3d
```

**优化**:
- 增加计算密度（使用更高阶格式）
- 减少MPI调用频率（每N步通信一次）
- 使用MPI+OpenMP混合并行

---

## 性能优化检查清单

### 编译优化

- [ ] `-O3` 优化级别
- [ ] `-march=native` 架构优化
- [ ] `-ffast-math` 浮点优化（如果允许）
- [ ] `-flto` 链接时优化（LTO）
- [ ] 启用 ccache 加速增量编译

```bash
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-O3 -march=native -flto"
```

### 内存优化

- [ ] 使用SoA布局（已实现）
- [ ] 一次分配所有临时缓冲（内存池）
- [ ] 避免频繁的malloc/free
- [ ] NUMA感知内存绑定（多socket系统）

```cpp
// 内存绑定示例（Linux）
#ifdef __linux__
#include <numaif.h>
void bind_memory_to_numa(void* ptr, size_t size, int numa_node) {
    // 将内存绑定到特定NUMA节点
    numa_tonode_memory(ptr, size, numa_node);
}
#endif
```

### MPI优化

- [ ] 使用非阻塞通信（Isend/Irecv）
- [ ] 最小化通信次数
- [ ] 合并多个MPI调用
- [ ] 使用MPI-IO并行输出

### 计算优化

- [ ] 向量化循环（编译器自动或显式SIMD）
- [ ] 减少浮点运算除法次数
- [ ] 使用快速数学库函数

```cpp
// 快速平方根近似
inline double fast_rsqrt(double x) {
    // Newton-Raphson 方法
    double y = x;
    y = 1.5 * y - 0.5 * x * y * y;
    y = 1.5 * y - 0.5 * x * y * y;
    return y;
}
```

---

## 总结

这份实现指南涵盖了：

✅ **完整的开发环境配置**
✅ **核心数据结构的完整实现**
✅ **物理方程的标准模板**
✅ **MPI并行的关键代码**
✅ **调试和优化建议**

**下一步**:
1. 完整实现所有Riemann求解器
2. 实现时间积分器
3. 集成Pipeline框架
4. 开发测试用例

---

**文档版本**: 1.0
**最后更新**: 2025年1月12日
