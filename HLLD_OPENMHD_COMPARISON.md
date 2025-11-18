# HLLD求解器与OpenMHD 3D重联实现详细对比分析

**日期**: 2025-11-18
**对比版本**:
- **fvm_3d**: 当前实现 (commit e2cdd23)
- **OpenMHD**: 3D_reconnection_gpu 版本

---

## 执行摘要 📋

### ✅ 正确实现的部分
1. **GLM Strang Splitting**: 完全正确，与OpenMHD一致
2. **电阻率模型**: 位置依赖sech²函数，参数匹配
3. **Harris Sheet磁场拓扑**: 正确的tanh(y)配置
4. **代码架构**: C++现代化设计优于Fortran
5. **并行化策略**: MPI + OpenMP + SIMD多级并行

### ❌ 严重问题
1. **HLLD求解器不完整**: 只实现4波简化版本（应为5波）
2. **能量更新公式错误**: 缺少Poynting flux项
3. **切向磁场/速度演化错误**: 公式缺少Bn²修正项
4. **缺少双星态**: 无U**L和U**R计算
5. **GLM源项系数错误**: 使用了cr/ch而非ch/cr

---

## 一、OpenMHD 3D重联GPU版本架构

### 1.1 文件结构
```
3D_reconnection_gpu/
├── main.cuf              # 主程序（CUDA Fortran）
├── model.f90             # Harris Sheet初始化
├── flux_solver.cuf       # HLLD Riemann求解器
├── flux_resistive.cuf    # 电阻MHD源项
├── glm_ss2.cuf           # GLM Strang splitting
├── rk.cuf                # TVD RK2时间积分
└── mpibc.cuf             # MPI边界通信
```

### 1.2 物理模型参数（OpenMHD标准配置）

```fortran
! Harris Sheet
beta = 0.2                    ! 等离子体beta
b1 = 0.03                     ! 扰动振幅
L_sheet = 1.0                 ! 电流片厚度

! 电阻率
eta0 = 1.0/1000.0            ! 背景: Rm₀ = 1000
eta01 = eta0 * (1000.0/60.0 - 1.0)  ! 增强: Rm₁ = 60
! η(r) = η₀ + (η₁-η₀)·sech²(r)

! GLM参数
ch = 0.2                      ! 波速
cr = 0.2                      ! 阻尼比
```

---

## 二、HLLD Riemann求解器详细对比

### 2.1 波结构对比

#### OpenMHD实现（5波HLLD）
```
      SL      SL*     SM      SR*     SR
  UL ──── U*L ──── U**L ─ U**R ──── U*R ──── UR
       ↑      ↑       ↑       ↑      ↑
     快波   Alfvén  接触   Alfvén   快波
```

**关键特性**:
- 5个独立状态：UL, U*L, U** (中心), U*R, UR
- 7个通量计算分支
- 双星态U**处理Alfvén波相互作用

#### fvm_3d实现（4波简化）
```
      SL       SM       SR
  UL ──── U*L ──── U*R ──── UR
       ↑       ↑       ↑
     快波    接触     快波
```

**问题**:
- ❌ 缺少Alfvén波（SL*, SR*）
- ❌ 无双星态U**
- ❌ 只有4个通量分支

---

### 2.2 波速计算对比

#### OpenMHD（正确）
```fortran
! 快磁声波速度（平方）
f1 = gamma * p_L
f2 = 4.0 * Bn**2
B2 = Bx**2 + By**2 + Bz**2
aL = ((f1 + B2) + sqrt(max((f1+B2)**2 - f1*f2, 0.d0))) / (2*rho_L)
aR = ((f1 + B2) + sqrt(max((f1+B2)**2 - f1*f2, 0.d0))) / (2*rho_R)

! 取最大值作为统一波速
f1 = sqrt(max(aL, aR))
SL = min(vn_L, vn_R) - f1
SR = max(vn_L, vn_R) + f1

! Alfvén波速
aVL = abs(Bn_hll) / sqrt(rho_L*)
aVR = abs(Bn_hll) / sqrt(rho_R*)
SL* = SM - aVL
SR* = SM + aVR
```

#### fvm_3d实现
```cpp
// riemann_hlld.cpp:49-54
double cf_L = mhd_->fast_speed(rho_L, p_L, B_x, By_L, Bz_L);
double cf_R = mhd_->fast_speed(rho_R, p_R, B_x, By_R, Bz_R);

double S_L = std::min(u_L - cf_L, u_R - cf_R);  // ⚠️ 不对称
double S_R = std::max(u_L + cf_L, u_R + cf_R);
```

**问题**:
1. ✅ 快磁声波速度计算正确
2. ❌ 波速选择不对称（应该用max(cf_L, cf_R)）
3. ❌ **完全缺少Alfvén波速计算**

---

### 2.3 中间态计算对比

#### 2.3.1 左中间态 U*L

##### OpenMHD（正确）
```fortran
! 密度
rho_L* = rho_L * (SL - vn_L) / (SL - SM)

! 分母因子（包含磁压修正）
f1 = 1.0 / (rho_L*(SL - vn_L)*(SL - SM) - Bn²)

! 切向磁场
Bty_L* = Bty_L * f1 * (rho_L*(SL - vn_L)² - Bn²)
Btz_L* = Btz_L * f1 * (rho_L*(SL - vn_L)² - Bn²)

! 切向速度（磁张力修正）
vty_L* = vty_L - f1 * Bn * Bty_L * (SM - vn_L)
vtz_L* = vtz_L - f1 * Bn * Btz_L * (SM - vn_L)

! 能量（包含Poynting flux）
e_L* = ((SL - vn_L)*e_L - pt_L*vn_L + pt*·SM
       + Bn*(v·B - SM*Bn - vty_L**Bty_L* - vtz_L**Btz_L*)) / (SL - SM)

其中:
  pt_L = p_L + 0.5*B²_L        # 总压力
  pt* = p* + 0.5*(Bn² + Bty_L*² + Btz_L*²)
  v·B = vn_L*Bn + vty_L*Bty_L + vtz_L*Btz_L
```

##### fvm_3d实现
```cpp
// riemann_hlld.cpp:187-254
double rho_Lstar = rho_L * (S_L - u_L) / (S_L - S_M);  // ✅ 正确

// ❌ 错误的切向磁场公式（缺少Bn²）
By_Lstar = By_L * (S_L - u_L) / (S_L - S_M);
Bz_Lstar = Bz_L * (S_L - u_L) / (S_L - S_M);

// ❌ 错误的切向速度公式
double delta_p_fact = px_jump / (rho_L * (S_L - u_L) * (S_L - S_M));
double v_Lstar = v_L - B_x * By_L / (rho_L * (S_L - u_L)) * delta_p_fact;

// ❌ 完全错误的能量更新
double p_L_scaled = p_L + 0.5 * (By_L * By_L + Bz_L * Bz_L);
double p_m_scaled = p_m + 0.5 * (By_Lstar * By_Lstar + Bz_Lstar * Bz_Lstar);
U_star(4) = U_L(4) + (p_m_scaled - p_L_scaled) * (S_L - u_L) / (S_L - S_M);
```

**错误分析**:

| 量 | OpenMHD公式 | fvm_3d公式 | 正确性 |
|---|------------|-----------|--------|
| ρ* | ✅ 正确 | ✅ 正确 | ✅ |
| By*, Bz* | f1·(ρ(S-v)² - Bn²) | (S-v)/(S-SM) | ❌ **缺少Bn²** |
| vy*, vz* | v - f1·Bn·B·(SM-v) | 错误公式 | ❌ **公式错误** |
| E* | 完整公式(7项) | 简化公式(2项) | ❌ **严重错误** |

---

#### 2.3.2 双星态 U** (fvm_3d完全缺失)

##### OpenMHD实现
```fortran
! Roe平均（跨越Alfvén波）
sqrt_rho_L = sqrt(rho_L*)
sqrt_rho_R = sqrt(rho_R*)
f1 = 1.0 / (sqrt_rho_L + sqrt_rho_R)
f2 = sign(1.0, Bn)

! 平均切向速度
vty_** = f1 * (sqrt_rho_L*vty_L* + sqrt_rho_R*vty_R*
              + (Bty_R* - Bty_L*)*f2)
vtz_** = f1 * (sqrt_rho_L*vtz_L* + sqrt_rho_R*vtz_R*
              + (Btz_R* - Btz_L*)*f2)

! 平均切向磁场
Bty_** = f1 * (sqrt_rho_L*Bty_R* + sqrt_rho_R*Bty_L*
              + sqrt_rho_L*sqrt_rho_R*(vty_R* - vty_L*)*f2)
Btz_** = f1 * (sqrt_rho_L*Btz_R* + sqrt_rho_R*Btz_L*
              + sqrt_rho_L*sqrt_rho_R*(vtz_R* - vtz_L*)*f2)

! 能量修正
e_L** = e_L* - sqrt_rho_L*(vty_L**Bty_L* + vtz_L**Btz_L*
                           - vty_***Bty_** - vtz_***Btz_**)*f2
```

##### fvm_3d实现
```cpp
// ❌ 完全缺失！没有双星态计算
```

**影响**:
- 无法正确处理Alfvén波相互作用
- 磁场旋转间断处理错误
- **磁重联时Alfvén波传播失真**

---

### 2.4 通量选择逻辑对比

#### OpenMHD（7分支）
```fortran
if (SL >= 0) then
   F = FL                          ! 分支1
else if (SR <= 0) then
   F = FR                          ! 分支2
else
   ! 计算HLL平均态
   U_hll = ...
   SM = U_hll(mn) / U_hll(ro)

   ! 计算U*L, U*R
   ...

   if (SL* >= 0) then
      F = FL + SL*(U*L - UL)       ! 分支3: U*L区域
   else if (SR* <= 0) then
      F = FR + SR*(U*R - UR)       ! 分支4: U*R区域
   else
      ! 计算U**
      ...

      if (SM >= 0) then
         F = FL + SL*U*L + SL**U**L  ! 分支5: U**L区域
      else
         F = FR + SR*U*R + SR**U**R  ! 分支6: U**R区域
      end if
   end if
end if
```

#### fvm_3d实现（4分支）
```cpp
// riemann_hlld.cpp:78-94
if (0.0 <= S_L) {
    F_hlld = F_L;                  // 分支1
} else if (0.0 >= S_R) {
    F_hlld = F_R;                  // 分支2
} else if (0.0 <= S_M) {
    F_hlld = F_L + S_L * (U_Lstar - UL_rot);  // 分支3
} else {
    F_hlld = F_R + S_R * (U_Rstar - UR_rot);  // 分支4
}
```

**❌ 缺少3个分支**: 无法正确处理Alfvén波区域！

---

## 三、电阻率模型对比 ✅

### 3.1 公式对比

#### OpenMHD
```fortran
r = sqrt(x**2 + y**2)
r = min(r, 25.0)  ! 防止溢出
eta = eta0 + eta01 * (cosh(r))**(-2)

其中:
  eta0 = 1.0/1000.0
  eta01 = eta0 * (1000.0/60.0 - 1.0) ≈ 0.01567
```

#### fvm_3d
```cpp
// resistive_mhd3d_advanced.hpp:68-78
double r = std::sqrt(x*x + y*y);
r = std::min(r, 25.0);
double sech_sq = 1.0 / std::cosh(r / localization_scale);
sech_sq *= sech_sq;
return eta0 + (eta1 - eta0) * sech_sq;

其中:
  eta0 = 1e-3
  eta1 = 1.667e-2
```

**✅ 完全一致**（数学等价）

---

### 3.2 电阻源项实现

#### OpenMHD (flux_resistive.cuf)
```fortran
! 电流密度（界面处）
JyS = -(U(i+1,j,k,bz) - U(i,j,k,bz)) / dx + ...
JzS = (U(i+1,j,k,by) - U(i,j,k,by)) / dx - ...

! Ohmic加热
F(i,j,k,en) += 0.5*eta_dx * (JyS*Bz - JzS*By)

! 磁场扩散
F(i,j,k,by) -= eta_dx * JzS
F(i,j,k,bz) += eta_dx * JyS
```

#### fvm_3d (resistive_mhd3d_advanced.cpp:278-324)
```cpp
// 电流密度J = ∇×B（单元中心）
Eigen::Vector3d J = compute_current_density(dx, dy, dz, ...);
double J_sq = J.squaredNorm();

// Ohmic加热
S(4) = eta * J_sq;

// 磁场扩散
Eigen::Vector3d lapl_B = compute_laplacian_B(dx, dy, dz, ...);
S(5) = eta * lapl_B(0);
S(6) = eta * lapl_B(1);
S(7) = eta * lapl_B(2);
```

**差异分析**:
- OpenMHD: 在通量处计算（界面值）
- fvm_3d: 在源项处计算（体积值）
- ✅ 两种方法都物理正确，数值差异很小

---

## 四、Harris Sheet初始条件对比

### 4.1 磁场配置

#### OpenMHD
```fortran
! 主磁场
Bx = tanh(y(j))
By = 0
Bz = 0

! Gaussian扰动
r2 = x(i)**2 + y(j)**2
dBx = -b1 * y(j) * exp(-r2/4.0)
dBy = +b1 * x(i) * exp(-r2/4.0)

Bx = Bx + dBx
By = By + dBy
```

#### fvm_3d
```cpp
// resistive_mhd3d_advanced.cpp:361-422
double Bx = harris.B0 * std::tanh(y_norm);

// m=1模扰动
double By = harris.perturbation_amplitude *
            std::sin(M_PI * x) *
            std::exp(-y_norm * y_norm);
```

**差异**:
- OpenMHD: Gaussian扰动 `exp(-r²/4)`
- fvm_3d: 正弦调制 `sin(πx)·exp(-y²)`

两种都是标准做法，fvm_3d的**更适合周期边界**。

---

### 4.2 密度分布

#### OpenMHD
```fortran
rho = 1.0 + (cosh(y))**(-2) / beta

其中 beta = 0.2
```

这确保了压力平衡：
```
p = 0.5 * beta * rho
∇p = ∇(B²/2μ₀) + ρ∇Φ  ! 力平衡
```

#### fvm_3d
```cpp
// 当前: 使用常密度（简化）
double rho = harris.n0;

// ⚠️ 注释中提到应该用:
// double sech_y_sq = 1.0 / std::cosh(y_norm);
// sech_y_sq *= sech_y_sq;
// double rho = harris.n0 * (1.0 + (1.0/harris.beta - 1.0) * sech_y_sq);
```

**❌ 问题**: 当前实现使用**常密度**，不满足MHD平衡！

**影响**:
- 初始条件不是严格平衡态
- 会产生声波振荡（但很快衰减）
- 对磁重联主要物理影响**较小**

---

### 4.3 压力分布

#### OpenMHD
```fortran
p = 0.5 * beta * rho
```

#### fvm_3d
```cpp
// resistive_mhd3d_advanced.cpp:406-410
double p_mag = 0.5 * Bx * Bx / MU0;
double p = harris.p0 - p_mag;
p = std::max(p, P_FLOOR);
```

**差异**:
- OpenMHD: `p ∝ ρ`（温度均匀）
- fvm_3d: `p = p₀ - B²/(2μ₀)`（总压力平衡）

fvm_3d的方法**更物理正确**（满足力平衡）。

---

## 五、GLM散度清理对比

### 5.1 Strang Splitting实现

#### OpenMHD (glm_ss2.f90)
```fortran
! D^(dt/2) ∘ L^(dt) ∘ D^(dt/2)

subroutine glm_ss2(U, ch, dt, ix, jx)
   f1 = exp(-0.5 * dt * ch / cr)
   U(:,:,ps) = U(:,:,ps) * f1
end subroutine
```

主循环中:
```fortran
call glm_ss2(U, ch, 0.5*dt, ix, jx)  ! D^(dt/2)
call tvdrk2(U, dt, ...)               ! L^(dt)
call glm_ss2(U, ch, 0.5*dt, ix, jx)  ! D^(dt/2)
```

#### fvm_3d (fvm_solver3d.cpp:140-235)
```cpp
// ========== GLM DAMPING: First half-step ==========
mhd->apply_glm_damping(U, dt_, 0.5);  // ψ *= exp(-0.5·dt·ch/cr)

// ========== MHD EVOLUTION: Full time step ==========
time_integrator_->step(state_, dt_, rhs_func);

// ========== GLM DAMPING: Second half-step ==========
mhd->apply_glm_damping(U, dt_, 0.5);
```

**✅ 完全一致**！包括:
- 分裂顺序: D-L-D
- 阻尼因子: `exp(-0.5·dt·ch/cr)`
- 与RK积分器的接口

---

### 5.2 GLM源项

#### OpenMHD
```fortran
! GLM演化方程（在通量中）
F(ps) = ch**2 * Bn        ! ψ通量
! 注意: 阻尼项通过Strang splitting处理，不在源项中
```

#### fvm_3d
```cpp
// resistive_mhd3d_advanced.cpp:326-333
double psi_source = -ch * div_B - (cr / ch) * psi;  // ❌ 错误！
```

**❌ 问题**:
- 源项中的阻尼系数应该是 `-(ch/cr)·ψ`
- 当前实现用了 `-(cr/ch)·ψ`，系数反了！

但由于Strang splitting中单独处理了阻尼（且公式正确），这个源项实际上**可能没有使用**。需要检查是否被调用。

---

### 5.3 GLM通量

#### OpenMHD
```fortran
! X方向通量
F(bx) = ps                  ! F(Bx) = ψ
F(ps) = ch**2 * bx          ! F(ψ) = ch²·Bx
```

#### fvm_3d
```cpp
// resistive_mhd3d_advanced.cpp:88-91
F(5) = psi;                 // F(Bx) = ψ
F(8) = glm_params_.ch * glm_params_.ch * Bx;  // F(ψ) = ch²·Bx
```

**✅ 完全正确**！

---

## 六、数值方法对比

### 6.1 时间积分

| 特性 | OpenMHD | fvm_3d |
|------|---------|--------|
| **方法** | TVD RK2 | TVD RK2/RK3（可配置） |
| **实现** | 固定 | 灵活 |
| **CFL** | 0.4 | 0.3（可配置） |

✅ fvm_3d更灵活

---

### 6.2 空间重构

#### OpenMHD
```fortran
! MUSCL重构（原始变量）
! minmod限制器
phi = max(0.0, min(1.0, r))  ! minmod
dq = phi * (q(i+1) - q(i))

! 重构
q_L = q(i) + 0.5 * dq
q_R = q(i+1) - 0.5 * dq
```

#### fvm_3d
```cpp
// muscl_reconstruction.cpp:84-124
// 混合变量重构
if (var_type == PRIMITIVE) {
    // 原始变量: u, v, w, p
    val = V(var);
} else {
    // 守恒变量: ρ, Bx, By, Bz, ψ
    val = U(var);
}

// Positivity limiter
if (V_L(0) < rho_floor || V_R(0) < rho_floor) {
    V_L(0) = V_center(0);
    V_R(0) = V_plus1(0);
}
```

**✅ fvm_3d优势**:
- 混合变量重构减少数值振荡
- Positivity-preserving limiter（参考Waagan 2009）
- 比OpenMHD更先进

---

### 6.3 性能优化

| 特性 | OpenMHD | fvm_3d |
|------|---------|--------|
| **CPU并行** | MPI + OpenMP | MPI + OpenMP |
| **GPU加速** | CUDA Fortran | 无 |
| **SIMD** | 编译器自动 | 显式指令 |
| **内存布局** | SoA | SoA |

fvm_3d虽然没有GPU支持，但**SIMD优化更明确**。

---

## 七、问题总结与修复优先级

### 🔴 **P0: 紧急修复（阻碍科学可信度）**

#### 1. HLLD能量公式错误
**文件**: `src/spatial/riemann_solvers/riemann_hlld.cpp:245-247, 313-315`

**当前**（错误）:
```cpp
double p_L_scaled = p_L + 0.5 * (By_L * By_L + Bz_L * Bz_L);
U_star(4) = U_L(4) + (p_m_scaled - p_L_scaled) * (S_L - u_L) / (S_L - S_M);
```

**应改为**（OpenMHD公式）:
```cpp
double vBL = u_L*B_x + v_L*By_L + w_L*Bz_L;
double pt_L = p_L + 0.5 * (B_x*B_x + By_L*By_L + Bz_L*Bz_L);
double pt_star = p_m + 0.5 * (B_x*B_x + By_Lstar*By_Lstar + Bz_Lstar*Bz_Lstar);

U_star(4) = ((S_L - u_L)*U_L(4) - pt_L*u_L + pt_star*S_M
            + B_x*(vBL - S_M*B_x - v_Lstar*By_Lstar - w_Lstar*Bz_Lstar))
            / (S_L - S_M);
```

**影响**: 能量不守恒，磁重联时非物理能量耗散

---

#### 2. 切向磁场演化公式错误
**文件**: `src/spatial/riemann_solvers/riemann_hlld.cpp:222-223, 290-291`

**当前**（错误）:
```cpp
By_Lstar = By_L * (S_L - u_L) / (S_L - S_M);
```

**应改为**:
```cpp
double f1 = 1.0 / (rho_L * (S_L - u_L) * (S_L - S_M) - B_x * B_x);
By_Lstar = By_L * f1 * (rho_L * (S_L - u_L) * (S_L - u_L) - B_x * B_x);
```

**影响**: 强磁场时数值不稳定

---

### 🟠 **P1: 高优先级（完整性）**

#### 3. 实现Alfvén波和双星态
**文件**: `src/spatial/riemann_solvers/riemann_hlld.cpp`

需要添加:
```cpp
// Alfvén波速
double aVL = std::abs(B_x) / std::sqrt(rho_Lstar);
double aVR = std::abs(B_x) / std::sqrt(rho_Rstar);
double S_Lstar = S_M - aVL;
double S_Rstar = S_M + aVR;

// 双星态计算
if (S_Lstar < 0.0 && S_Rstar > 0.0) {
    Eigen::VectorXd U_LstarStar = compute_double_star_L(...);
    Eigen::VectorXd U_RstarStar = compute_double_star_R(...);
    // ...
}
```

**影响**: 无法正确模拟Alfvén波，磁场旋转间断失真

---

### 🟡 **P2: 中优先级（精度改进）**

#### 4. GLM源项系数修正
**文件**: `src/physics/resistive_mhd3d_advanced.cpp:331`

**当前**: `-(cr/ch)*psi`
**应改**: `-(ch/cr)*psi`

但需要**先检查**这个源项是否被使用，因为Strang splitting已经处理了阻尼。

---

#### 5. Harris Sheet密度平衡
**文件**: `src/physics/resistive_mhd3d_advanced.cpp:379-392`

从常密度改为:
```cpp
double sech_y_sq = 1.0 / std::cosh(y_norm);
sech_y_sq *= sech_y_sq;
double rho = harris.n0 * (1.0 + (1.0/harris.beta - 1.0) * sech_y_sq);
```

**影响**: 初始平衡更精确，减少声波振荡

---

## 八、基准测试建议

### 8.1 必须通过的测试

1. **Brio-Wu激波管**（1D MHD Riemann问题）
   - 目的: 验证HLLD求解器正确性
   - 预期: 当前实现**可能失败**（无Alfvén波）

2. **Orszag-Tang涡旋**（2D MHD湍流）
   - 目的: 验证能量守恒和数值稳定性
   - 预期: 当前实现**可能失败**（能量公式错误）

3. **Harris重联**（3D磁重联）
   - 目的: 与OpenMHD对比主要物理结果
   - 诊断量:
     - 重联率 dΨ/dt
     - 磁能/动能演化
     - div(B) < 1e-9
     - 岛链形成时间

---

### 8.2 与OpenMHD对比流程

```bash
# 1. OpenMHD运行
cd OpenMHD/3D_reconnection_gpu
# 编辑参数匹配
./a.out

# 2. fvm_3d运行（相同参数）
cd fvm_3d/build
./harris_sheet_3d 64 32 32

# 3. 提取诊断量
# - t vs. max|By|（重联指标）
# - t vs. KE, ME（能量演化）
# - t vs. max|div B|（散度约束）

# 4. 对比分析
python compare_results.py openmhd_data.txt fvm_3d_data.txt
```

---

## 九、总体评价

### ✅ **优秀之处（超越OpenMHD）**

1. **代码架构**: 现代C++设计，模块化清晰
2. **GLM实现**: Strang splitting完全正确
3. **重构方法**: 混合变量 + positivity limiter更先进
4. **并行优化**: MPI + OpenMP + SIMD三级并行
5. **可维护性**: 比Fortran CUDA代码更易维护

### ❌ **严重缺陷（亟需修复）**

1. **HLLD求解器不完整**: 只实现了4波，缺少Alfvén波
2. **能量公式错误**: 缺少Poynting flux，可能导致能量不守恒
3. **切向分量演化错误**: 强磁场时不稳定
4. **缺少双星态**: 无法正确处理Alfvén波相互作用

### 📊 **综合评分**

| 方面 | 分数 | 评语 |
|------|------|------|
| **GLM清理** | 10/10 | 完美实现 |
| **电阻率模型** | 10/10 | 与OpenMHD一致 |
| **Riemann求解器** | **3/10** | 不完整且有错误 |
| **重构方法** | 9/10 | 优于OpenMHD |
| **代码质量** | 9/10 | 现代化设计 |
| **并行性能** | 8/10 | 缺GPU但CPU优化好 |
| **Harris Sheet** | 7/10 | 缺密度平衡 |

**总体**: **7.0/10**（受HLLD求解器严重拖累）

---

## 十、修复路线图

### Week 1: P0紧急修复
- [ ] 修正HLLD能量公式（添加Poynting flux）
- [ ] 修正切向磁场/速度演化（添加Bn²项）
- [ ] 单元测试: Sod激波管（验证能量守恒）

### Week 2: P1完整性实现
- [ ] 实现Alfvén波速计算
- [ ] 实现双星态U**, U**
- [ ] 实现7分支通量选择
- [ ] 单元测试: Brio-Wu激波管

### Week 3: 验证与基准测试
- [ ] Orszag-Tang涡旋
- [ ] Harris重联（与OpenMHD对比）
- [ ] 性能基准测试
- [ ] 文档更新

### Week 4: P2精度改进（可选）
- [ ] Harris密度平衡
- [ ] GLM源项系数检查
- [ ] 波速估算对称化

---

## 十一、结论

**当前状态**: fvm_3d在工程架构和代码质量上**优于OpenMHD**，但核心物理求解器**存在严重缺陷**。

**关键问题**: HLLD求解器是**简化的4波版本**（更像HLLC），而非完整的5波HLLD。能量更新公式**完全错误**。

**建议**:
1. **立即修复P0问题**，否则无法用于科学研究
2. **添加标准基准测试**，建立持续验证机制
3. **与OpenMHD交叉验证**，确保物理正确性
4. **考虑临时方案**: 在修复前可退回到HLL求解器

**最终判断**: 项目架构优秀，但**当前不建议用于发表科学论文**，必须先修复HLLD求解器。

---

**参考文献**:
1. Miyoshi & Kusano (2005), JCP 208:315-344
2. Dedner et al. (2002), JCP 175:645-673
3. Waagan (2009), JCP 228:8609-8626
4. OpenMHD: https://github.com/zenitani/OpenMHD
