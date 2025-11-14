# 时间积分器向量化优化

## 改进概述

将时间积分器从4层嵌套循环改为单层向量化操作，显著提高性能和可读性。

## 主要变化

### 1. Field3D 新增向量化操作

在 `include/core/field3d.hpp` 中添加了以下BLAS风格的操作：

- **axpy(a, x, y)**: `this = a*x + y`
- **linear_combination(a, x, b, y)**: `this = a*x + b*y`
- **linear_combination_3(a, x, b, y, c, z)**: `this = a*x + b*y + c*z`
- **add_scaled(a, x)**: `this = this + a*x` (原地AXPY)
- **scale(a)**: `this = a * this`

所有操作都使用 `#pragma omp simd` 指令以启用SIMD向量化。

### 2. 时间积分器重写

#### ForwardEuler (之前: 4层循环)
```cpp
// 之前: 20行代码，4层循环
for (int v = 0; v < nvars; v++) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                U_current(v, i, j, k) += dt * (*temp_rhs_)(v, i, j, k);
            }
        }
    }
}

// 现在: 1行代码，向量化操作
U_current.add_scaled(dt, *temp_rhs_);
```

#### RK2 (之前: 45行代码 → 现在: 27行代码)
```cpp
// Stage 1: U^(1) = U^n + dt * RHS(U^n)
rhs(U_current, *temp_rhs_);
temp_stage1_->axpy(dt, *temp_rhs_, U_current);

// Stage 2: U^(n+1) = 0.5 * U^n + 0.5 * U^(1) + 0.5 * dt * RHS(U^(1))
rhs(*temp_stage1_, *temp_rhs_);
U_current.linear_combination_3(
    0.5, U_current,
    0.5, *temp_stage1_,
    0.5 * dt, *temp_rhs_
);
```

#### RK3 (之前: 100行代码 → 现在: 58行代码)
完全消除了所有嵌套循环，使用3个向量化操作。

## 性能优势

### 1. 内存访问优化
- **之前**: 4层循环导致差的缓存局部性
- **现在**: 单次顺序遍历整个数组，最优缓存利用

### 2. SIMD向量化
- **之前**: 编译器难以自动向量化深层嵌套循环
- **现在**: `#pragma omp simd` 明确指示编译器向量化，利用AVX/AVX512指令

### 3. 代码简洁性
- **之前**: RK3积分器100行代码
- **现在**: RK3积分器58行代码 (减少42%)
- 更清晰的数学表达，更易维护

### 4. 预期性能提升
对于典型的3D网格 (nx=ny=nz=100, nvars=5):
- **数据量**: 5 × 100³ = 5M doubles = 40 MB
- **理论加速比**:
  - 消除循环开销: 1.5-2×
  - SIMD向量化 (AVX2): 4×
  - 更好的缓存利用: 1.2-1.5×
  - **总体预期**: 5-10× 加速

## 编译器优化建议

为了获得最佳性能，使用以下编译器标志：

```bash
-O3                    # 最高优化等级
-march=native          # 针对当前CPU优化 (启用AVX2/AVX512)
-ffast-math            # 快速浮点运算
-fopenmp-simd          # OpenMP SIMD支持
```

在 `CMakeLists.txt` 中:
```cmake
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math")
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(fvm3d_lib PUBLIC OpenMP::OpenMP_CXX)
endif()
```

## 验证

所有测试用例编译通过，无错误：
- ✅ Forward Euler
- ✅ RK2 (SSP-RK2)
- ✅ RK3 (SSP-RK3)
- ✅ 所有示例程序

## 未来改进

1. **多线程并行**: 添加 `#pragma omp parallel for` 用于多核CPU
2. **GPU加速**: 使用CUDA/HIP将操作移至GPU
3. **BLAS库集成**: 可选使用Intel MKL或OpenBLAS的daxpy/dgemv
4. **自适应向量长度**: 根据数据大小自动选择最优块大小
