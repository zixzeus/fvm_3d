# 3D C++ MPI并行框架设计文档

**项目名称**: FVM3D - 3D Finite Volume Method Framework with MPI
**基于**: 现有Python 2D框架
**目标**: 高性能3D MPI并行计算框架
**工期**: 20天
**预算**: 4000元

---

## 目录

1. [项目概述](#项目概述)
2. [技术选型与架构决策](#技术选型与架构决策)
3. [整体架构设计](#整体架构设计)
4. [MPI并行策略](#mpi并行策略)
5. [模块设计](#模块设计)
6. [内存布局与性能优化](#内存布局与性能优化)
7. [算法迁移清单](#算法迁移清单)
8. [20天实施计划](#20天实施计划)
9. [关键代码框架](#关键代码框架)
10. [测试策略](#测试策略)

---

## 项目概述

### 核心目标

**从**: Python 2D串行框架
**到**: C++ 3D MPI并行框架

**预期性能提升**:
- **10-50x**: Python → C++（编译优化）
- **N倍**: MPI并行（N = CPU核数）
- **总计**: 100-500x+（在多节点集群上）

### 关键特性

✅ **3D几何**: 完整3维空间离散
✅ **MPI并行**: 区域分解并行计算
✅ **高性能**: C++17 + 现代编译器优化
✅ **算法完整**: 迁移所有Python算法
✅ **可扩展**: 模块化设计，易于添加新算法

---

## 技术选型与架构决策

### 编程语言与标准

**C++17**
- 原因: 现代特性（结构化绑定、if constexpr、并行算法）
- 兼容性: 所有主流编译器支持

### 核心依赖库

| 库 | 版本 | 用途 | 必需性 |
|-------|---------|------|----------|
| **MPI** | 3.0+ | 并行通信 | **必须** |
| **Eigen3** | 3.4+ | 线性代数 | **必须** |
| **HDF5** | 1.12+ | 数据I/O | **推荐** |
| **Catch2** | 3.x | 单元测试 | **推荐** |
| **CMake** | 3.15+ | 构建系统 | **必须** |
| **OpenMP** | 4.5+ | 节点内并行 | 可选 |

**依赖说明**:

**MPI实现选择**:
- **OpenMPI** 4.x: 开源，功能全，推荐开发
- **MPICH** 4.x: 稳定，HPC常用
- **Intel MPI**: 商业，性能最优（如有许可）

**Eigen3**:
- 轻量级纯头文件库
- 优秀的编译期优化
- 向量化支持（SSE/AVX）

**HDF5**:
- 大规模科学数据标准格式
- MPI-IO并行写入
- Python/VisIt/ParaView兼容

---

## 整体架构设计

### 分层架构

```
┌─────────────────────────────────────────────────────────┐
│  应用层 (Applications)                                  │
│  drivers/, tests/, examples/                            │
├─────────────────────────────────────────────────────────┤
│  求解器层 (Solvers)                                     │
│  FVMSolver3D, DGSolver3D                                │
├─────────────────────────────────────────────────────────┤
│  并行层 (Parallel Layer) - MPI                          │
│  MPIDomainDecomposer, MPIHaloExchange, MPIReducer       │
├─────────────────────────────────────────────────────────┤
│  流水线层 (Pipeline)                                    │
│  ComputePipeline3D, StageExecutor                       │
├─────────────────────────────────────────────────────────┤
│  计算模块层 (Computation Modules)                       │
│  - Boundary (边界条件)                                  │
│  - Reconstruction (空间重构: WENO, DG)                  │
│  - Riemann Solvers (HLL, HLLC, HLLD)                    │
│  - TimeIntegrators (Euler, RK2/3/4)                     │
│  - Source Terms (源项)                                   │
├─────────────────────────────────────────────────────────┤
│  物理层 (Physics)                                       │
│  EulerEquations3D, MHDEquations3D, ResistiveMHD3D       │
├─────────────────────────────────────────────────────────┤
│  数据层 (Data Structures)                               │
│  Grid3D, Field3D, StateVector, MPIGrid3D                │
├─────────────────────────────────────────────────────────┤
│  基础设施层 (Infrastructure)                            │
│  Logger, Timer, Config, I/O (HDF5)                      │
└─────────────────────────────────────────────────────────┘
```

### 目录结构

```
fvm3d/
├── CMakeLists.txt                 # 根CMake配置
├── README.md
├── LICENSE
│
├── include/                        # 头文件（公共接口）
│   └── fvm3d/
│       ├── core/
│       │   ├── grid3d.hpp         # 3D网格
│       │   ├── field3d.hpp        # 场数据容器
│       │   ├── state_vector.hpp   # 状态向量
│       │   └── types.hpp          # 类型定义
│       ├── parallel/               # MPI并行
│       │   ├── mpi_grid3d.hpp     # MPI分布式网格
│       │   ├── halo_exchange.hpp  # 幽灵单元通信
│       │   ├── domain_decomposer.hpp
│       │   └── mpi_io.hpp         # 并行I/O
│       ├── physics/
│       │   ├── physics_base.hpp
│       │   ├── euler3d.hpp
│       │   ├── mhd3d.hpp
│       │   └── resistive_mhd3d.hpp
│       ├── spatial/
│       │   ├── reconstruction/
│       │   │   ├── reconstruction_base.hpp
│       │   │   ├── constant.hpp
│       │   │   ├── muscl.hpp
│       │   │   ├── weno.hpp
│       │   │   └── dg_modal.hpp   # DG P0-P2
│       │   └── riemann/
│       │       ├── riemann_base.hpp
│       │       ├── lax_friedrichs.hpp
│       │       ├── tvdlf.hpp
│       │       ├── hll.hpp
│       │       ├── hllc.hpp
│       │       └── hlld.hpp
│       ├── temporal/
│       │   ├── time_integrator_base.hpp
│       │   ├── forward_euler.hpp
│       │   ├── rk2.hpp
│       │   ├── rk3.hpp
│       │   └── rk4.hpp
│       ├── boundary/
│       │   ├── boundary_base.hpp
│       │   ├── periodic.hpp
│       │   ├── reflective.hpp
│       │   └── transmissive.hpp
│       ├── pipeline/
│       │   ├── compute_pipeline3d.hpp
│       │   └── stage_executor.hpp
│       ├── solver/
│       │   ├── fvm_solver3d.hpp
│       │   └── dg_solver3d.hpp
│       └── utils/
│           ├── logger.hpp
│           ├── timer.hpp
│           ├── config.hpp
│           └── hdf5_writer.hpp
│
├── src/                            # 实现文件
│   ├── core/
│   ├── parallel/
│   ├── physics/
│   ├── spatial/
│   ├── temporal/
│   ├── boundary/
│   ├── pipeline/
│   ├── solver/
│   └── utils/
│
├── tests/                          # 单元测试
│   ├── test_grid3d.cpp
│   ├── test_mpi_halo_exchange.cpp
│   ├── test_riemann_solvers.cpp
│   ├── test_euler3d.cpp
│   └── ...
│
├── examples/                       # 示例程序
│   ├── blast_wave_3d.cpp
│   ├── magnetic_reconnection.cpp
│   ├── cme_eruption.cpp
│   ├── kh_instability.cpp
│   └── rt_instability.cpp
│
├── benchmarks/                     # 性能测试
│   ├── weak_scaling.cpp
│   └── strong_scaling.cpp
│
├── python/                         # Python绑定（可选）
│   └── pybind11_wrapper.cpp
│
└── docs/                           # 文档
    ├── API.md
    ├── MPI_GUIDE.md
    └── PERFORMANCE.md
```

---

## MPI并行策略

### 区域分解（Domain Decomposition）

#### 1D分解（初期，简单）

```
3D域按 X 方向切分:

全局域: [0, Lx] × [0, Ly] × [0, Lz]
N个MPI进程

Rank 0: [0, Lx/N] × [0, Ly] × [0, Lz]
Rank 1: [Lx/N, 2Lx/N] × [0, Ly] × [0, Lz]
...
Rank N-1: [(N-1)Lx/N, Lx] × [0, Ly] × [0, Lz]

优点: 实现简单，通信少（仅2个邻居）
缺点: 扩展性受限（最多N = Nx个进程）
```

#### 3D分解（高级，强扩展）

```
3D域按 X, Y, Z 三个方向切分:

进程拓扑: Px × Py × Pz

例如 64个进程 = 4×4×4

Rank (i,j,k):
  X: [i*Lx/Px, (i+1)*Lx/Px]
  Y: [j*Ly/Py, (j+1)*Ly/Py]
  Z: [k*Lz/Pz, (k+1)*Lz/Pz]

邻居: 最多26个（3D情况，含角点）

优点: 最佳扩展性
缺点: 通信复杂（需优化）
```

**推荐策略**:
- **阶段1（开发）**: 1D分解
- **阶段2（优化）**: 3D分解

---

### 幽灵单元通信（Halo Exchange）

#### 数据结构

```cpp
// 每个MPI进程的本地域
class MPIGrid3D {
    // 全局索引范围
    int global_nx, global_ny, global_nz;

    // 本进程的本地索引范围（不含ghost）
    int local_nx, local_ny, local_nz;

    // 本地数据（含ghost层）
    int nx_total = local_nx + 2*nghost;
    int ny_total = local_ny + 2*nghost;
    int nz_total = local_nz + 2*nghost;

    // 邻居进程rank
    int neighbor_xm, neighbor_xp;  // X方向邻居
    int neighbor_ym, neighbor_yp;  // Y方向邻居
    int neighbor_zm, neighbor_zp;  // Z方向邻居

    // 本地数据
    Field3D<double> state;  // [nvars, nx_total, ny_total, nz_total]
};
```

#### 通信模式

**非阻塞通信（推荐）**:

```cpp
void MPIHaloExchange::exchange(Field3D<double>& field) {
    std::vector<MPI_Request> requests;

    // X方向
    if (has_neighbor_xm) {
        // 发送左边界 → neighbor_xm
        MPI_Isend(left_send_buffer, ..., neighbor_xm, tag, comm, &req);
        requests.push_back(req);

        // 接收左ghost ← neighbor_xm
        MPI_Irecv(left_recv_buffer, ..., neighbor_xm, tag, comm, &req);
        requests.push_back(req);
    }

    if (has_neighbor_xp) {
        // 发送右边界 → neighbor_xp
        // 接收右ghost ← neighbor_xp
        ...
    }

    // Y, Z方向类似

    // 等待所有通信完成
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}
```

**优化**:
- **Pack/Unpack**: 连续内存复制
- **重叠计算**: 边界通信时，内部单元可以先计算

---

### 负载均衡

#### 均匀网格

对于**均匀笛卡尔网格**:
```
每个进程分配相同数量的单元
→ 天然负载均衡
```

#### 自适应网格（未来扩展）

- 动态负载平衡（Zoltan库）
- 按计算量重新分配

---

### MPI I/O（并行输出）

**HDF5 + MPI-IO**:

```cpp
void MPIWriter::write_snapshot(const Field3D<double>& field,
                                const std::string& filename) {
    // 创建MPI-IO文件访问属性
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);

    // 打开文件（所有进程）
    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC,
                           H5P_DEFAULT, plist);

    // 每个进程写入自己的数据块
    // 使用hyperslab选择
    hsize_t offset[3] = {global_i_start, global_j_start, global_k_start};
    hsize_t count[3] = {local_nx, local_ny, local_nz};

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    // 集体写入
    H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
             xfer_plist, field.data());

    H5Fclose(file);
}
```

**优点**:
- 单个文件（无需后处理合并）
- 高性能并行写入
- 标准格式（Python/VisIt可读）

---

## 模块设计

### 核心数据结构

#### Grid3D（基础3D网格）

```cpp
// include/fvm3d/core/grid3d.hpp
namespace fvm3d {

struct GridGeometry3D {
    double Lx, Ly, Lz;          // 域尺寸
    int nx, ny, nz;             // 单元数
    double dx, dy, dz;          // 单元尺寸

    // 构造函数
    GridGeometry3D(double Lx_, double Ly_, double Lz_,
                   int nx_, int ny_, int nz_)
        : Lx(Lx_), Ly(Ly_), Lz(Lz_)
        , nx(nx_), ny(ny_), nz(nz_)
        , dx(Lx/nx), dy(Ly/ny), dz(Lz/nz) {}

    // 单元中心坐标
    double cell_center_x(int i) const { return (i + 0.5) * dx; }
    double cell_center_y(int j) const { return (j + 0.5) * dy; }
    double cell_center_z(int k) const { return (k + 0.5) * dz; }
};

class Grid3D {
public:
    Grid3D(const GridGeometry3D& geom, int nghost = 2)
        : geom_(geom), nghost_(nghost)
        , nx_total_(geom.nx + 2*nghost)
        , ny_total_(geom.ny + 2*nghost)
        , nz_total_(geom.nz + 2*nghost) {}

    // 访问器
    int nx() const { return geom_.nx; }
    int ny() const { return geom_.ny; }
    int nz() const { return geom_.nz; }
    int nghost() const { return nghost_; }

    // 单元中心坐标
    double x(int i) const { return geom_.cell_center_x(i); }
    double y(int j) const { return geom_.cell_center_y(j); }
    double z(int k) const { return geom_.cell_center_z(k); }

private:
    GridGeometry3D geom_;
    int nghost_;
    int nx_total_, ny_total_, nz_total_;
};

} // namespace fvm3d
```

---

#### Field3D（3D场数据容器）

```cpp
// include/fvm3d/core/field3d.hpp
namespace fvm3d {

// 3D多维数组（基于Eigen或原始指针）
template<typename T>
class Field3D {
public:
    Field3D(int nvars, int nx, int ny, int nz)
        : nvars_(nvars), nx_(nx), ny_(ny), nz_(nz)
        , size_(nvars * nx * ny * nz)
        , data_(new T[size_]) {
        std::fill(data_.get(), data_.get() + size_, T(0));
    }

    // 访问：data(var, i, j, k)
    T& operator()(int var, int i, int j, int k) {
        return data_[index(var, i, j, k)];
    }

    const T& operator()(int var, int i, int j, int k) const {
        return data_[index(var, i, j, k)];
    }

    // 原始指针（用于MPI通信）
    T* data() { return data_.get(); }
    const T* data() const { return data_.get(); }

    // 尺寸
    int nvars() const { return nvars_; }
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }
    size_t size() const { return size_; }

private:
    int nvars_, nx_, ny_, nz_;
    size_t size_;
    std::unique_ptr<T[]> data_;

    // 内存布局：行优先（C风格）
    // index = var*(nx*ny*nz) + i*(ny*nz) + j*nz + k
    inline size_t index(int var, int i, int j, int k) const {
        return var*(nx_*ny_*nz_) + i*(ny_*nz_) + j*nz_ + k;
    }
};

// 别名
using StateField3D = Field3D<double>;

} // namespace fvm3d
```

**内存布局选择**:
- **行优先（C风格）**: `[var][i][j][k]` 连续
- 优点: 与C++数组习惯一致，MPI通信友好

---

#### MPIGrid3D（MPI分布式网格）

```cpp
// include/fvm3d/parallel/mpi_grid3d.hpp
namespace fvm3d {

class MPIGrid3D {
public:
    MPIGrid3D(MPI_Comm comm, const GridGeometry3D& global_geom,
              int nghost = 2);

    // 全局信息
    int global_nx() const { return global_geom_.nx; }
    int global_ny() const { return global_geom_.ny; }
    int global_nz() const { return global_geom_.nz; }

    // 本地信息
    int local_nx() const { return local_nx_; }
    int local_ny() const { return local_ny_; }
    int local_nz() const { return local_nz_; }

    // 本地数据尺寸（含ghost）
    int nx_total() const { return local_nx_ + 2*nghost_; }
    int ny_total() const { return local_ny_ + 2*nghost_; }
    int nz_total() const { return local_nz_ + 2*nghost_; }

    // 全局索引转换
    int global_i_start() const { return i_start_; }
    int global_j_start() const { return j_start_; }
    int global_k_start() const { return k_start_; }

    // MPI信息
    MPI_Comm comm() const { return comm_; }
    int rank() const { return rank_; }
    int size() const { return size_; }

    // 邻居信息
    int neighbor_xm() const { return neighbor_xm_; }
    int neighbor_xp() const { return neighbor_xp_; }
    // ... ym, yp, zm, zp

private:
    MPI_Comm comm_;
    int rank_, size_;

    GridGeometry3D global_geom_;
    int nghost_;

    // 本地域大小
    int local_nx_, local_ny_, local_nz_;

    // 全局索引起点
    int i_start_, j_start_, k_start_;

    // 邻居rank（-1表示无邻居或边界）
    int neighbor_xm_, neighbor_xp_;
    int neighbor_ym_, neighbor_yp_;
    int neighbor_zm_, neighbor_zp_;
};

} // namespace fvm3d
```

---

### 物理方程

#### PhysicsBase（抽象基类）

```cpp
// include/fvm3d/physics/physics_base.hpp
namespace fvm3d {

class PhysicsBase {
public:
    virtual ~PhysicsBase() = default;

    // 必须实现的接口
    virtual int num_vars() const = 0;
    virtual std::vector<std::string> var_names() const = 0;

    // 守恒变量 ↔ 原始变量
    virtual void conservative_to_primitive(
        const Eigen::VectorXd& U, Eigen::VectorXd& V) const = 0;

    virtual void primitive_to_conservative(
        const Eigen::VectorXd& V, Eigen::VectorXd& U) const = 0;

    // 通量计算
    virtual void compute_flux_x(
        const Eigen::VectorXd& U, Eigen::VectorXd& F) const = 0;

    virtual void compute_flux_y(
        const Eigen::VectorXd& U, Eigen::VectorXd& G) const = 0;

    virtual void compute_flux_z(
        const Eigen::VectorXd& U, Eigen::VectorXd& H) const = 0;

    // 最大波速（CFL）
    virtual double max_wave_speed(const Eigen::VectorXd& U, int direction) const = 0;

    // 源项（默认无）
    virtual void compute_source(
        const Eigen::VectorXd& U,
        double x, double y, double z,
        Eigen::VectorXd& S) const {
        S.setZero();
    }
};

} // namespace fvm3d
```

---

#### EulerEquations3D

```cpp
// include/fvm3d/physics/euler3d.hpp
namespace fvm3d {

class EulerEquations3D : public PhysicsBase {
public:
    EulerEquations3D(double gamma = 1.4) : gamma_(gamma) {}

    int num_vars() const override { return 5; }

    std::vector<std::string> var_names() const override {
        return {"rho", "rho_u", "rho_v", "rho_w", "E"};
    }

    void conservative_to_primitive(
        const Eigen::VectorXd& U, Eigen::VectorXd& V) const override {
        double rho = U(0);
        double u = U(1) / rho;
        double v = U(2) / rho;
        double w = U(3) / rho;
        double E = U(4);

        double ke = 0.5 * rho * (u*u + v*v + w*w);
        double p = (gamma_ - 1.0) * (E - ke);

        V.resize(5);
        V << rho, u, v, w, p;
    }

    void compute_flux_x(const Eigen::VectorXd& U, Eigen::VectorXd& F) const override {
        Eigen::VectorXd V;
        conservative_to_primitive(U, V);

        double rho = V(0), u = V(1), v = V(2), w = V(3), p = V(4);
        double E = U(4);

        F.resize(5);
        F << rho*u,
             rho*u*u + p,
             rho*u*v,
             rho*u*w,
             u*(E + p);
    }

    // compute_flux_y, compute_flux_z 类似

    double max_wave_speed(const Eigen::VectorXd& U, int direction) const override {
        Eigen::VectorXd V;
        conservative_to_primitive(U, V);

        double rho = V(0), u = V(1), v = V(2), w = V(3), p = V(4);
        double c = std::sqrt(gamma_ * p / rho);  // 声速

        if (direction == 0) return std::abs(u) + c;
        if (direction == 1) return std::abs(v) + c;
        return std::abs(w) + c;
    }

private:
    double gamma_;
};

} // namespace fvm3d
```

---

### Riemann求解器

#### RiemannSolverBase

```cpp
// include/fvm3d/spatial/riemann/riemann_base.hpp
namespace fvm3d {

class RiemannSolverBase {
public:
    RiemannSolverBase(std::shared_ptr<PhysicsBase> physics)
        : physics_(physics) {}

    virtual ~RiemannSolverBase() = default;

    // 求解Riemann问题
    virtual void solve(
        const Eigen::VectorXd& U_L,
        const Eigen::VectorXd& U_R,
        int direction,
        Eigen::VectorXd& F) const = 0;

protected:
    std::shared_ptr<PhysicsBase> physics_;
};

} // namespace fvm3d
```

---

#### HLLSolver

```cpp
// include/fvm3d/spatial/riemann/hll.hpp
namespace fvm3d {

class HLLSolver : public RiemannSolverBase {
public:
    using RiemannSolverBase::RiemannSolverBase;

    void solve(const Eigen::VectorXd& U_L,
               const Eigen::VectorXd& U_R,
               int direction,
               Eigen::VectorXd& F) const override {
        // 计算波速
        double S_L = compute_wave_speed_left(U_L, U_R, direction);
        double S_R = compute_wave_speed_right(U_L, U_R, direction);

        // 通量
        Eigen::VectorXd F_L, F_R;
        physics_->compute_flux(U_L, direction, F_L);
        physics_->compute_flux(U_R, direction, F_R);

        // HLL通量
        if (S_L >= 0.0) {
            F = F_L;
        } else if (S_R <= 0.0) {
            F = F_R;
        } else {
            F = (S_R*F_L - S_L*F_R + S_L*S_R*(U_R - U_L)) / (S_R - S_L);
        }
    }

private:
    double compute_wave_speed_left(...) const {
        // 估计左波速
        double lambda_L = physics_->max_wave_speed(U_L, direction);
        double lambda_R = physics_->max_wave_speed(U_R, direction);
        return std::min(lambda_L, lambda_R);
    }

    double compute_wave_speed_right(...) const {
        // 估计右波速
        ...
    }
};

} // namespace fvm3d
```

**其他求解器**（HLLC, HLLD等）类似实现。

---

### 时间积分器

#### TimeIntegratorBase

```cpp
// include/fvm3d/temporal/time_integrator_base.hpp
namespace fvm3d {

class TimeIntegratorBase {
public:
    TimeIntegratorBase(std::shared_ptr<PhysicsBase> physics)
        : physics_(physics) {}

    virtual ~TimeIntegratorBase() = default;

    // 推进一个时间步
    virtual void step(StateField3D& state, double dt, double t) = 0;

protected:
    std::shared_ptr<PhysicsBase> physics_;

    // 计算RHS（通用）
    void compute_rhs(const StateField3D& state, StateField3D& rhs);
};

} // namespace fvm3d
```

---

#### RungeKutta3

```cpp
// include/fvm3d/temporal/rk3.hpp
namespace fvm3d {

class RungeKutta3 : public TimeIntegratorBase {
public:
    using TimeIntegratorBase::TimeIntegratorBase;

    void step(StateField3D& state, double dt, double t) override {
        // 保存初始状态
        StateField3D U0 = state;
        StateField3D rhs(state.nvars(), state.nx(), state.ny(), state.nz());

        // Stage 1
        compute_rhs(state, rhs);
        state = U0 + dt * rhs;

        // Stage 2
        compute_rhs(state, rhs);
        state = 0.75*U0 + 0.25*state + 0.25*dt*rhs;

        // Stage 3
        compute_rhs(state, rhs);
        state = (1.0/3.0)*U0 + (2.0/3.0)*state + (2.0/3.0)*dt*rhs;
    }
};

} // namespace fvm3d
```

---

## 内存布局与性能优化

### SoA vs AoS

**Structure of Arrays (SoA)** - 推荐:
```cpp
// 内存连续，缓存友好
Field3D<double> state(nvars, nx, ny, nz);
state(0, i, j, k) = rho;   // 所有rho连续
state(1, i, j, k) = rho_u; // 所有rho_u连续
```

**Array of Structures (AoS)** - 不推荐:
```cpp
struct Cell { double rho, rho_u, rho_v, rho_w, E; };
Cell cells[nx][ny][nz];  // 单元内变量连续，跨单元不连续
```

**为什么SoA更快？**
1. SIMD向量化友好
2. 缓存行利用率高
3. MPI通信连续

---

### 循环优化

#### 缓存友好的循环顺序

**对于行优先存储 `[var][i][j][k]`**:

```cpp
// ✓ 推荐（k最内层，内存连续）
for (int var = 0; var < nvars; ++var)
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
                // 访问 state(var, i, j, k)

// ✗ 不推荐（i最内层，跨步访问）
for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            for (int var = 0; var < nvars; ++var)
                // 跨步访问，缓存miss多
```

---

#### SIMD向量化

**编译器自动向量化**:
```cpp
// 简单循环，编译器可自动向量化
for (int k = 0; k < nz; ++k) {
    state(0, i, j, k) = state(0, i, j, k) * factor;  // ✓ 可向量化
}
```

**显式SIMD（可选，高级）**:
```cpp
#include <immintrin.h>  // AVX intrinsics

// 手动AVX-256（8个double）
__m256d factor_vec = _mm256_set1_pd(factor);
for (int k = 0; k < nz; k += 4) {
    __m256d data = _mm256_load_pd(&state(0, i, j, k));
    data = _mm256_mul_pd(data, factor_vec);
    _mm256_store_pd(&state(0, i, j, k), data);
}
```

---

### 内存池（Memory Pool）

**避免频繁分配/释放**:
```cpp
class MemoryPool {
    std::vector<StateField3D> temp_fields_;

public:
    StateField3D& get_temp(int nvars, int nx, int ny, int nz) {
        // 重用已分配的内存
        if (temp_fields_.empty()) {
            temp_fields_.emplace_back(nvars, nx, ny, nz);
        }
        return temp_fields_.back();
    }
};
```

---

### OpenMP节点内并行（可选）

**MPI + OpenMP混合**:
```cpp
// 外层MPI（跨节点）
// 内层OpenMP（节点内多核）

#pragma omp parallel for collapse(3)
for (int i = nghost; i < nx+nghost; ++i)
    for (int j = nghost; j < ny+nghost; ++j)
        for (int k = nghost; k < nz+nghost; ++k) {
            // 计算单元(i,j,k)的RHS
            compute_cell_rhs(state, i, j, k);
        }
```

**编译**: `mpic++ -fopenmp ...`

**运行**: `export OMP_NUM_THREADS=4; mpirun -np 16 ./fvm3d`
- 16个MPI进程
- 每个进程4个OpenMP线程
- 总计64个线程

---

## 算法迁移清单

### 空间离散算法

| 算法 | Python文件 | C++目标文件 | 难度 | 优先级 |
|------|-----------|------------|------|--------|
| **Lax-Friedrichs** | `spatial/flux_calculation/lax_friedrichs_flux.py` | `spatial/riemann/lax_friedrichs.cpp` | 低 | P0 |
| **TVDLF** | （可能在drivers中）| `spatial/riemann/tvdlf.cpp` | 低 | P1 |
| **HLL** | `spatial/riemann_solvers.py:HLLSolver` | `spatial/riemann/hll.cpp` | 中 | P0 |
| **HLLC** | `spatial/riemann_solvers.py:HLLCSolver` | `spatial/riemann/hllc.cpp` | 中 | P0 |
| **HLLD** | `spatial/riemann_solvers.py:HLLDSolver` | `spatial/riemann/hlld.cpp` | 高 | P1 |
| **DG P0** | `core/tensor_dg_2d.py` | `spatial/reconstruction/dg_p0.cpp` | 中 | P2 |
| **DG P1** | `core/tensor_dg_2d.py` | `spatial/reconstruction/dg_p1.cpp` | 高 | P2 |
| **DG P2** | `core/tensor_dg_2d.py` | `spatial/reconstruction/dg_p2.cpp` | 高 | P3 |

**优先级说明**:
- P0: 第一周完成（基础）
- P1: 第二周完成（重要）
- P2: 第三周完成（高级）
- P3: 时间允许则实现

---

### 时间离散算法

| 算法 | Python文件 | C++目标文件 | 难度 | 优先级 |
|------|-----------|------------|------|--------|
| **Forward Euler** | `temporal/time_integrators.py:ForwardEuler` | `temporal/forward_euler.cpp` | 低 | P0 |
| **RK2** | `temporal/time_integrators.py:RungeKutta2` | `temporal/rk2.cpp` | 低 | P0 |
| **RK3** | `temporal/time_integrators.py:RungeKutta3` | `temporal/rk3.cpp` | 中 | P0 |
| **RK4** | `temporal/time_integrators.py:RungeKutta4` | `temporal/rk4.cpp` | 中 | P1 |

---

### 物理方程

| 方程 | Python文件 | C++目标文件 | 难度 | 优先级 |
|------|-----------|------------|------|--------|
| **Euler 3D** | `physics/euler_equations.py` | `physics/euler3d.cpp` | 中 | P0 |
| **MHD 3D** | `physics/mhd_equations.py` | `physics/mhd3d.cpp` | 高 | P1 |
| **Resistive MHD 3D** | `physics/resistive_mhd_equations_flux.py` | `physics/resistive_mhd3d.cpp` | 高 | P2 |

---

## 20天实施计划

### 阶段划分

| 阶段 | 天数 | 任务 | 交付物 |
|------|------|------|--------|
| **阶段1: 基础架构** | 5天 | 搭建框架骨架 | 编译通过的空框架 |
| **阶段2: 核心算法** | 7天 | 实现基础算法 | 串行版本可运行 |
| **阶段3: MPI并行** | 5天 | 实现MPI并行 | 并行版本通过测试 |
| **阶段4: 测试优化** | 3天 | 测试用例和性能优化 | 完整可发布版本 |

---

### 详细日程

#### 第1周（Day 1-7）: 基础框架 + 核心算法

**Day 1-2: 环境搭建与基础架构**
- [ ] 创建CMake项目结构
- [ ] 配置MPI、Eigen3、HDF5依赖
- [ ] 实现 `Grid3D`, `Field3D` 基础数据结构
- [ ] 实现 `Logger`, `Timer` 工具类
- [ ] 编写第一个测试（test_grid3d）

**Day 3-4: 物理方程层**
- [ ] 实现 `PhysicsBase` 抽象基类
- [ ] 实现 `EulerEquations3D`
- [ ] 单元测试: 守恒变量↔原始变量转换
- [ ] 单元测试: 通量计算

**Day 5-7: Riemann求解器**
- [ ] 实现 `RiemannSolverBase`
- [ ] 实现 `LaxFriedrichsSolver`
- [ ] 实现 `HLLSolver`
- [ ] 实现 `HLLCSolver`（Euler专用）
- [ ] 单元测试: Sod激波管问题

**交付物**: 串行版本3D Euler求解器（无MPI）

---

#### 第2周（Day 8-14）: 完整串行实现

**Day 8-9: 时间积分器**
- [ ] 实现 `TimeIntegratorBase`
- [ ] 实现 `ForwardEuler`
- [ ] 实现 `RungeKutta2`
- [ ] 实现 `RungeKutta3`
- [ ] 单元测试: 时间精度收敛性

**Day 10-11: 边界条件**
- [ ] 实现 `BoundaryBase`
- [ ] 实现 `PeriodicBC`, `ReflectiveBC`, `TransmissiveBC`
- [ ] 单元测试: 各种边界条件

**Day 12-14: 求解器整合**
- [ ] 实现 `FVMSolver3D`
- [ ] 整合：网格 + 物理 + Riemann + 时间积分 + 边界
- [ ] 第一个完整示例: 3D Blast Wave
- [ ] 串行性能测试

**交付物**: 完整串行3D FVM求解器

---

#### 第3周（Day 15-19）: MPI并行化

**Day 15-16: MPI基础设施**
- [ ] 实现 `MPIGrid3D`（1D区域分解）
- [ ] 实现 `MPIDomainDecomposer`
- [ ] 单元测试: 域分解正确性

**Day 17-18: Halo通信**
- [ ] 实现 `MPIHaloExchange`（非阻塞通信）
- [ ] 实现边界条件的MPI扩展（跨进程周期BC）
- [ ] 单元测试: 幽灵单元一致性

**Day 19: MPI I/O**
- [ ] 实现 `MPIHDFWriter`（并行HDF5输出）
- [ ] 测试: 多进程并行写入
- [ ] 验证: Python读取HDF5文件

**交付物**: MPI并行版本（1D分解）

---

#### 第4周（Day 20）: 测试与收尾

**Day 20: 最终测试与文档**
- [ ] 运行所有测试用例
  - Blast wave 3D
  - KH instability 3D
  - RT instability 3D
- [ ] 弱扩展性测试（固定每进程负载）
- [ ] 强扩展性测试（固定总负载）
- [ ] 编写 README.md
- [ ] 编写 MPI_GUIDE.md
- [ ] 代码清理和注释

**交付物**: 完整可发布的v1.0版本

---

### 里程碑检查点

| 检查点 | 日期 | 标准 |
|--------|------|------|
| **Checkpoint 1** | Day 7 | 串行Euler求解器通过Sod测试 |
| **Checkpoint 2** | Day 14 | 3D Blast Wave串行模拟成功 |
| **Checkpoint 3** | Day 19 | MPI并行通过2进程测试 |
| **最终验收** | Day 20 | 所有测试通过 + 文档完整 |

---

## 关键代码框架

### CMakeLists.txt（根目录）

```cmake
cmake_minimum_required(VERSION 3.15)
project(FVM3D VERSION 1.0.0 LANGUAGES CXX)

# C++17标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 优化选项
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")

# 查找依赖
find_package(MPI REQUIRED)
find_package(Eigen3 3.4 REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX HL)

# 可选：OpenMP
find_package(OpenMP)

# 包含目录
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${MPI_CXX_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

# 库目录
add_subdirectory(src)

# 测试
enable_testing()
add_subdirectory(tests)

# 示例
add_subdirectory(examples)

# 安装
install(DIRECTORY include/ DESTINATION include)
```

---

### src/CMakeLists.txt

```cmake
# 收集所有源文件
file(GLOB_RECURSE FVM3D_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/core/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/parallel/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/physics/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/spatial/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/temporal/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pipeline/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/solver/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/*.cpp
)

# 创建库
add_library(fvm3d SHARED ${FVM3D_SOURCES})

# 链接依赖
target_link_libraries(fvm3d
    PUBLIC
        MPI::MPI_CXX
        Eigen3::Eigen
        ${HDF5_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(fvm3d PUBLIC OpenMP::OpenMP_CXX)
endif()

# 安装
install(TARGETS fvm3d
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
```

---

### examples/blast_wave_3d.cpp

```cpp
#include <fvm3d/solver/fvm_solver3d.hpp>
#include <fvm3d/physics/euler3d.hpp>
#include <fvm3d/spatial/riemann/hllc.hpp>
#include <fvm3d/temporal/rk3.hpp>
#include <fvm3d/boundary/periodic.hpp>
#include <fvm3d/utils/hdf5_writer.hpp>
#include <mpi.h>

int main(int argc, char** argv) {
    // 初始化MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 全局域
    fvm3d::GridGeometry3D global_geom(1.0, 1.0, 1.0, 100, 100, 100);

    // MPI分布式网格
    fvm3d::MPIGrid3D mpi_grid(MPI_COMM_WORLD, global_geom, 2);

    // 物理方程
    auto physics = std::make_shared<fvm3d::EulerEquations3D>(1.4);

    // Riemann求解器
    auto riemann = std::make_shared<fvm3d::HLLCSolver>(physics);

    // 时间积分器
    auto integrator = std::make_shared<fvm3d::RungeKutta3>(physics);

    // 边界条件
    auto bc = std::make_shared<fvm3d::PeriodicBC>();

    // 求解器
    fvm3d::FVMSolver3D solver(mpi_grid, physics, riemann, integrator, bc);

    // 初始条件：Blast Wave
    solver.set_initial_condition([](double x, double y, double z, auto& U) {
        double r = std::sqrt(x*x + y*y + z*z);
        if (r < 0.1) {
            U << 1.0, 0.0, 0.0, 0.0, 10.0;  // 高压区
        } else {
            U << 1.0, 0.0, 0.0, 0.0, 0.1;   // 低压区
        }
    });

    // 运行模拟
    double t_end = 0.2;
    double cfl = 0.4;

    solver.run(t_end, cfl);

    // 输出结果
    fvm3d::HDF5Writer writer(MPI_COMM_WORLD);
    writer.write("blast_wave_output.h5", solver.get_state(), mpi_grid);

    if (rank == 0) {
        std::cout << "Simulation completed successfully!\n";
    }

    MPI_Finalize();
    return 0;
}
```

**编译**:
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8
```

**运行**:
```bash
mpirun -np 4 ./examples/blast_wave_3d
```

---

## 测试策略

### 单元测试（Catch2）

```cpp
// tests/test_euler3d.cpp
#include <catch2/catch_test_macros.hpp>
#include <fvm3d/physics/euler3d.hpp>

TEST_CASE("EulerEquations3D: Conservative to Primitive", "[euler3d]") {
    fvm3d::EulerEquations3D euler(1.4);

    Eigen::VectorXd U(5);
    U << 1.0, 0.5, 0.0, 0.0, 2.5;  // rho=1, u=0.5, v=0, w=0, E=2.5

    Eigen::VectorXd V;
    euler.conservative_to_primitive(U, V);

    REQUIRE(V(0) == Approx(1.0));      // rho
    REQUIRE(V(1) == Approx(0.5));      // u
    REQUIRE(V(2) == Approx(0.0));      // v
    REQUIRE(V(3) == Approx(0.0));      // w
    REQUIRE(V(4) > 0.0);               // p > 0
}

TEST_CASE("EulerEquations3D: Flux calculation", "[euler3d]") {
    // 测试通量计算
    ...
}
```

---

### 集成测试（MPI）

```cpp
// tests/test_mpi_halo_exchange.cpp
#include <catch2/catch_test_macros.hpp>
#include <fvm3d/parallel/mpi_grid3d.hpp>
#include <fvm3d/parallel/halo_exchange.hpp>
#include <mpi.h>

TEST_CASE("MPI Halo Exchange", "[mpi]") {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    REQUIRE(size >= 2);  // 至少2个进程

    // 创建MPI网格
    fvm3d::GridGeometry3D geom(1.0, 1.0, 1.0, 10, 10, 10);
    fvm3d::MPIGrid3D mpi_grid(MPI_COMM_WORLD, geom, 1);

    // 创建数据
    fvm3d::StateField3D field(1,
                              mpi_grid.nx_total(),
                              mpi_grid.ny_total(),
                              mpi_grid.nz_total());

    // 填充本地数据（rank作为值）
    for (int i = 0; i < mpi_grid.nx_total(); ++i)
        for (int j = 0; j < mpi_grid.ny_total(); ++j)
            for (int k = 0; k < mpi_grid.nz_total(); ++k)
                field(0, i, j, k) = static_cast<double>(rank);

    // Halo通信
    fvm3d::MPIHaloExchange exchanger(mpi_grid);
    exchanger.exchange(field);

    // 验证：ghost单元应该是邻居的rank值
    if (rank == 0 && size > 1) {
        // Rank 0的右ghost应该是Rank 1的值
        double right_ghost = field(0, mpi_grid.nx_total()-1, 0, 0);
        REQUIRE(right_ghost == 1.0);
    }
}
```

---

### 物理测试用例

#### 1. Blast Wave 3D

**初始条件**:
```
高压球: r < 0.1, p = 10.0
低压背景: r >= 0.1, p = 0.1
密度均匀: rho = 1.0
```

**预期结果**:
- 球形激波向外传播
- 对称性保持（球形）

---

#### 2. Kelvin-Helmholtz Instability 3D

**初始条件**:
```
上层: ρ=2, u=0.5, v=0
下层: ρ=1, u=-0.5, v=0
界面: y=0.5，加小扰动
```

**预期结果**:
- 界面涡旋形成
- KH不稳定性特征结构

---

#### 3. Rayleigh-Taylor Instability 3D

**初始条件**:
```
重流体在上: ρ=2, p=2.5
轻流体在下: ρ=1, p=1.5
重力: g=-0.1（向下）
界面扰动: Δy = 0.01*sin(2πx)
```

**预期结果**:
- 重流体下沉，轻流体上升
- 蘑菇云结构形成

---

### 性能测试

#### 弱扩展性（Weak Scaling）

**固定每进程负载**:
```
1进程: 100³ 单元
2进程: 200×100×100 单元
4进程: 400×100×100 单元
...
```

**理想**: 运行时间不变（效率100%）

---

#### 强扩展性（Strong Scaling）

**固定总负载**:
```
总网格: 400³ 单元

1进程: 400³
2进程: 200³ × 2
4进程: 100³ × 4
8进程: 50³ × 8
```

**目标**: 加速比 > 70%（8进程时）

---

## 风险与应对

### 技术风险

| 风险 | 概率 | 影响 | 应对策略 |
|------|------|------|----------|
| MPI通信bug导致结果错误 | 中 | 高 | 充分单元测试，逐步验证 |
| 性能不达预期 | 低 | 中 | Profile分析，针对优化 |
| 3D算法理解偏差 | 低 | 高 | 先实现2D验证，再扩展3D |
| HDF5并行I/O问题 | 中 | 低 | 提供后备的串行I/O |

---

### 进度风险

| 风险 | 概率 | 影响 | 应对策略 |
|------|------|------|----------|
| 某个算法卡住超时 | 中 | 中 | 降低优先级，先完成核心 |
| 测试用例不通过 | 低 | 高 | 预留最后3天调试时间 |
| 依赖库编译问题 | 低 | 中 | 使用conda/spack统一环境 |

---

## 交付清单

### 代码交付

- [ ] 完整的C++源代码（GitHub仓库）
- [ ] CMake构建系统
- [ ] 单元测试（覆盖率>80%）
- [ ] 至少5个示例程序
- [ ] MPI运行脚本

### 文档交付

- [ ] README.md（快速入门）
- [ ] API.md（API文档）
- [ ] MPI_GUIDE.md（MPI使用指南）
- [ ] PERFORMANCE.md（性能测试报告）
- [ ] 算法验证报告（测试用例结果）

### 性能数据

- [ ] 串行vs Python性能对比（至少10x）
- [ ] MPI弱扩展性曲线
- [ ] MPI强扩展性曲线
- [ ] 内存使用分析

---

## 未来扩展方向

### v2.0（3个月后）

- [ ] 3D区域分解（Px × Py × Pz）
- [ ] AMR自适应网格
- [ ] GPU加速（CUDA/HIP）
- [ ] 更多物理：辐射、热传导

### v3.0（6个月后）

- [ ] Python绑定（pybind11）
- [ ] 在线可视化（VTK）
- [ ] 动态负载均衡（Zoltan）

---

## 总结

本设计文档为**3D C++ MPI并行FVM框架**提供了完整的技术方案：

✅ **明确的架构设计**：分层模块化，易于扩展
✅ **详细的实施计划**：20天分4阶段，里程碑清晰
✅ **完整的代码框架**：关键类和接口定义
✅ **严格的测试策略**：单元测试 + 物理验证 + 性能测试
✅ **MPI并行核心**：区域分解 + Halo通信 + 并行I/O

**预期成果**:
- 相比Python框架：**100-500x**加速
- MPI扩展性：在64核上效率>**70%**
- 代码质量：测试覆盖率>**80%**

---

**文档版本**: 1.0
**创建日期**: 2025年1月12日
**作者**: FVM3D Development Team
