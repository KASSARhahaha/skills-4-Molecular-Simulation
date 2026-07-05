# Patterns: Reusable Algorithms & Techniques

从 "Understanding Molecular Simulation" (Frenkel & Smit, 3rd ed.) 提炼的可复用算法模式与决策规则。

## 1. Sampling Construction Patterns (采样构造)

### 1.1 Four-step MC construction recipe (Ch 6)
任意系综的 MC 构造四步法:
1. 写出目标分布 $\pi(\xi)$ (取决于系综)
2. 应用细致平衡: $\pi(o)\, \alpha(o \to n)\, \text{acc}(o \to n) = \pi(n)\, \alpha(n \to o)\, \text{acc}(n \to o)$
3. 选择试探矩阵 $\alpha$
4. 导出接受规则 $\text{acc}$

**用途**: NPT, $\mu$VT, semigrand, Gibbs, expanded ensemble, parallel tempering。

### 1.2 Bias sampling + weight correction (Ch 12)
- 用偏置 $\alpha(o \to n) = f[U(n)]$ 产生试探
- 接受规则修正偏置: $\text{acc} = \min\!\big[1, \dfrac{f(U_o)}{f(U_n)} e^{-\beta \Delta U}\big]$
- 例: CBMC 用 Rosenbluth 权重 $W$

### 1.3 Super-detailed balance (Ch 12)
对连续偏置: 对每组固定试探方向 $\Omega^*$ 强制细致平衡。证明更简单; 结果与试探数 $k$ 无关。

## 2. Reference State Patterns (参考态模式)

### 2.1 Choose analytic reference state
| System 体系 | Reference 参考 | Why 原因 |
|---|---|---|
| Liquid 液体 | Ideal gas (B_2 known) | Density TI from low $\rho$ |
| Solid 固体 | Einstein crystal ($F$ analytic) | Springs anchored to lattice |
| Hard spheres 硬球 | Carnahan-Starling $F$ | Closed form |
| Polymer chain 高分子 | Non-self-avoiding ($Q_{\text{int}}$ analytic) | Boltzmann-sample angles |

### 2.2 Thermodynamic integration template
```
1. Pick reference state with F_ref analytic
2. Define U(λ) interpolating U_ref → U_real
3. Sample ⟨∂U/∂λ⟩_λ at λ = 0, 0.1, ..., 1.0
4. Numerical integration: ΔF = ∫₀¹ ⟨∂U/∂λ⟩_λ dλ
5. Consistency check: Gibbs-Bogoliubov ∂²F/∂λ² ≤ 0
```

公式: $\Delta F = \int_0^1 \langle \partial U / \partial \lambda \rangle_\lambda\, d\lambda$

### 2.3 Avoid first-order phase transition
- 不跨越一阶相变 (hysteresis)
- 绕过方法: (a) 外场, (b) single-occupancy cell, (c) lattice switching, (d) 在 $(T, \rho)$ 平面绕过临界点

## 3. Force/Energy Patterns (力与能量)

### 3.1 Cutoff + tail correction
- 短程 (LJ-like): 截断于 $r_c$，加尾部修正
  $$U_{\text{tail}} = \dfrac{N \rho}{2} \int_{r_c}^\infty 4 \pi r^2\, u(r)\, dr$$
- Shifted-force: $u_{\text{sf}}(r) = u(r) - u(r_c) - (r - r_c)\, u'(r_c)$，使 $u$ 和 $u'$ 在 $r_c$ 处都连续

### 3.2 Long-range (Coulomb, dipole)
- Ewald 分解: $\dfrac{1}{r} = \dfrac{\text{erfc}(\alpha r)}{r} + \dfrac{\text{erf}(\alpha r)}{r}$
- PME/SPME: FFT 加速 Fourier 部分, $O(N \log N)$
- Reaction field: 截断 + 外部 $\varepsilon_{\text{RF}}$ 连续介质

### 3.3 Neighbor lists (App I)
- Verlet list + skin: 重建触发 $\max_i |\mathbf{r}_i(t) - \mathbf{r}_i(t_0)| > \delta/2$
- Cell list: 格大小 $\geq r_c$，只查相邻 26 格 (3D)
- 组合: cell list 用于重建 Verlet list

## 4. Integration Patterns (积分模式)

### 4.1 Verlet family (Ch 4)
- Position Verlet: $O(\Delta t^4)$ position, $O(\Delta t^2)$ velocity
- Velocity Verlet: 等价轨迹，速度定义在整数步
- Leap-frog: 速度在半整数步
- 全部 symplectic (辛)、time-reversible (可逆)、长期漂移小

Velocity Verlet 公式:
$$\mathbf{r}(t + \Delta t) = \mathbf{r}(t) + \mathbf{v}(t)\, \Delta t + \dfrac{\mathbf{F}(t)}{2 m}\, \Delta t^2$$
$$\mathbf{v}(t + \Delta t) = \mathbf{v}(t) + \dfrac{\mathbf{F}(t) + \mathbf{F}(t + \Delta t)}{2 m}\, \Delta t$$

### 4.2 Multiple time step (r-RESPA, Ch 14)
```
Forces split: F_fast (bonded), F_mid (short non-bonded), F_slow (long non-bonded)
Outer step Δt = n · δt
Trotter: exp(iL_slow Δt/2) · exp(iL_mid Δt/2) · [exp(iL_fast δt)]^n · exp(iL_mid Δt/2) · exp(iL_slow Δt/2)
```

### 4.3 Constraints (Ch 14)
- SHAKE: 不动点迭代，逐键修正直至 $\sigma_k < \text{tol}$
- RATTLE/LINCS: velocity-Verlet 形式，更稳
- Fixman 修正: $\Delta F_{\text{corr}} = k_B T \ln |H|^{1/2}$

## 5. Thermostat/Barostat Patterns (Ch 7)

### 5.1 Choice based on what you trade
| Want 需求 | Use 用 | Don't use 不用 |
|---|---|---|
| Static properties only | Andersen, Bussi | — |
| Transport ($D, \eta, \lambda_T$) | NVE, Nosé-Hoover, Bussi | Andersen, Berendsen, Langevin |
| Equilibration only | Berendsen | (never for sampling) |
| Polymer / Brownian | Langevin | (loses hydrodynamics) |

### 5.2 Extended Lagrangian framework
- 在 $L$ 中加扩展变量 ($s$ for $T$, $V$ for $P$, $\mathbf{h}$ for stress)
- $H' = K_{\text{extra}} + U_{\text{extra}} + (\text{target constraint})$
- 质量 $Q$ / 惯性 调节耦合强度

Nosé-Hoover Hamiltonian: $\mathcal{H}_{\text{Nosé}} = \sum_i \dfrac{\mathbf{p}_i^2}{2 m_i s^2} + U(\mathbf{r}^N) + \dfrac{p_s^2}{2 Q} + g k_B T \ln s$

## 6. Measurement Patterns (Ch 5)

### 6.1 Fluctuation → response
- $C_V = \dfrac{\langle E^2 \rangle - \langle E \rangle^2}{k_B T^2}$ (NVT)
- $\kappa_T = \dfrac{\langle V^2 \rangle - \langle V \rangle^2}{\langle V \rangle k_B T}$ (NPT)
- $\varepsilon = 1 + \dfrac{4 \pi}{3 V k_B T} \langle M^2 \rangle$ (NVE, dielectric)

### 6.2 Time correlation → transport
| Coefficient 系数 | Formula 公式 |
|---|---|
| $D_s$ (self-diffusion) | $\int_0^\infty \langle \mathbf{v}(0) \cdot \mathbf{v}(t) \rangle / 3\, dt$ |
| $\eta$ (shear viscosity) | $\dfrac{V}{k_B T} \int_0^\infty \langle \sigma_{xy}(0)\, \sigma_{xy}(t) \rangle\, dt$ |
| $\lambda_T$ (thermal conductivity) | $\dfrac{V}{k_B T^2} \int_0^\infty \langle j_z(0)\, j_z(t) \rangle\, dt$ |
| $\sigma_e$ (electric conductivity) | $\dfrac{V}{k_B T} \int_0^\infty \langle j_x^e(0)\, j_x^e(t) \rangle\, dt$ |

### 6.3 Error estimation
- Flyvbjerg-Petersen block average: $\sigma$ from plateau
- $O(n)$ block summing for long correlations
- Statistical error $\sim 1/\sqrt{N}$ (independent samples)

## 7. Rare Event Patterns (Ch 15)

### 7.1 Bennett-Chandler recipe
```
1. Identify reaction coordinate q
2. Sample P(q) via umbrella sampling → k^TST = ⟨|q̇|⟩ e^{-βF(q*)} / ∫_{-∞}^{q*} e^{-βF(q)} dq
3. Release MD from q = q* → κ = lim_{t→∞} ⟨f_B(t)⟩_release
4. k_{A→B} = k^TST × κ
```

### 7.2 Blue-Moon alternative
- 在 MD 中约束 $q = q^*$ (现有约束代码即可)
- Fixman 修正: $P(q) \propto |H|^{-1/2} \langle \delta(\sigma) \rangle_c$
- 平均力 $\langle \partial H / \partial q \rangle_{\sigma=0}$ → 自由能梯度

### 7.3 Reaction coordinate validation
- Committor test: 在 $q = q^*$ 处释放 MD，统计到达 B 的比例
- Ideal: 直方图峰值在 0.5
- Bad: 双峰 0/1 → 反应坐标错误

### 7.4 TPS / TIS / FFS
- TPS: 路径系综 $C(t) = Z_{AB}(t) / Z_A$
- TIS: 分界面 $\lambda_i$，逐面计算条件概率 $P(\lambda_{i+1} | \lambda_i)$
- FFS: $k_{A \to 1} = 1 / \tau_{\text{MFPT}}$，$\tau_{\text{MFPT}}$ 从 Brute-force MFPT 反演

## 8. Acceleration Patterns (Ch 13)

### 8.1 Parallel tempering setup
- 选 $T_1 < T_2 < \dots < T_n$
- 间距: $\Delta T / T \sim 1/\sqrt{N}$
- 交换接受: $\text{acc} = \min\!\big[1, e^{(\beta_i - \beta_j)(U_j - U_i)}\big]$

### 8.2 Hamiltonian PT
- 每副本不同势能 (例: 高分子柔软度)
- 接受: $\text{acc} = \min\!\big[1, e^{-\beta(U_i^{(j)} - U_i^{(i)}) + \beta(U_j^{(i)} - U_j^{(j)})}\big]$

### 8.3 Expanded ensemble
- 单系统在 $\beta_i$ 间跳变
- 权重 $w_i = e^{\beta_i F - \gamma_i}$; 迭代 $\gamma$ 使采样均等
- 副产品: $F(T_i) - F(T_j)$

## 9. Mesoscopic Patterns (Ch 16)

### 9.1 Preserve momentum for hydrodynamics
- Langevin: $\times$ (breaks momentum)
- DPD: $\checkmark$ pairwise dissipative + random; FDT $\sigma^2 = 2 k_B T \gamma^D$
- Lowe-Andersen: $\checkmark$ local pairwise

DPD 力分解:
$$\mathbf{F}_i = \sum_{j \neq i} \big(\mathbf{f}_{ij}^C + \mathbf{f}_{ij}^D + \mathbf{f}_{ij}^R\big)$$
- $\mathbf{f}_{ij}^D = -\gamma\, \omega^D(r_{ij})\, (\mathbf{v}_{ij} \cdot \hat{\mathbf{r}}_{ij})\, \hat{\mathbf{r}}_{ij}$
- $\mathbf{f}_{ij}^R = \sigma\, \omega^R(r_{ij})\, \xi_{ij}\, \hat{\mathbf{r}}_{ij}$
- Español-Warren: $\omega^D = [\omega^R]^2$, $\sigma^2 = 2 k_B T \gamma$

### 9.2 Mesoscopic solvent selection
| Need 需求 | Method 方法 |
|---|---|
| Static only | Implicit solvent + effective potential |
| Hydrodynamic interactions | DPD or MPC |
| Heat transport | SDPD |
| Dilute gas | DSMC |

## 10. Linear Response Patterns (App F)

### 10.1 Kubo response function
$$\chi_{AA}(t) = -\beta\, \langle A(0)\, \dot{A}(t) \rangle\, \Theta(t)$$

### 10.2 Dissipation-fluctuation (energy dissipation)
$$\dot{E} = -\pi \omega\, |f_\omega|^2\, \text{Im}[\chi_{AA}(\omega)]$$

### 10.3 Static fluctuation-response
$$\chi_{AA}(0) = \beta\, \sigma_A^2 = \beta\, \big(\langle A^2 \rangle - \langle A \rangle^2\big)$$

### 10.4 Kramers-Kronig (causality)
$$\text{Re}\, \chi(\omega) = \dfrac{1}{\pi}\, \mathcal{P}\!\int_{-\infty}^\infty \dfrac{\text{Im}\, \chi(\omega')}{\omega' - \omega}\, d\omega'$$

## 11. Non-Equilibrium Patterns (App D, E)

### 11.1 Entropy production (second law)
$$\sigma = \mathbf{J}_q \cdot \nabla(1/T) - \sum_i \mathbf{J}_i \cdot \nabla(\mu_i / T) \geq 0$$

### 11.2 Onsager reciprocal + Green-Kubo
- $J_i = \sum_j L_{ij}\, X_j$, $L_{ij} = L_{ji}$
- $L_{ij} = \dfrac{1}{k_B} \int_0^\infty \langle J_j(0)\, J_i(t) \rangle\, dt$

### 11.3 Jarzynski / Crooks
- Jarzynski: $\langle e^{-\beta W} \rangle = e^{-\beta \Delta F}$
- Crooks: $\dfrac{P_F(W)}{P_R(-W)} = e^{\beta(W - \Delta F)}$

## 12. Validation Patterns (验证)

### 12.1 Always check
- Total energy conservation in NVE (drift = bug indicator)
- $\langle P \rangle$ matches imposed $P$ in NPT-MC (virial self-consistency)
- $g(r)$ from force method vs binning (equilibrium diagnostic)
- Committor distribution at proposed $q^*$ (coordinate quality)

### 12.2 Common anti-patterns to avoid
- Cutoff Coulomb (use Ewald/PME)
- Berendsen for sampling (only equilibration)
- Force-field without Fixman correction
- Naive Verlet + constraints (use SHAKE)
- Random insertion of long chains (use CBMC)
- Linear response in strong perturbation
- Box-Muller in hot loop (use Ziggurat)
- Hand-written FFT (use FFTW)

## 13. RNG and Special Function Patterns (App J)

### 13.1 Random number generators
- 推荐: Mersenne Twister (MT19937), PCG, Xoshiro
- 反对: `rand()`, linear congruential

### 13.2 Gaussian sampling
- Box-Muller: $g_1 = \sqrt{-2 \ln u_1}\, \cos(2 \pi u_2)$, $g_2 = \sqrt{-2 \ln u_1}\, \sin(2 \pi u_2)$
- Ziggurat (Marsaglia): 比 Box-Muller 快 3-5×

### 13.3 Special functions
- $\text{erfc}(x)$: Abramowitz-Stegun 近似 ($\sim 10^{-7}$) 或 Chebyshev 展开
- B-spline (PME): Cox-de Boor 递归, $W_{2p}$ 用于 SPME 电荷分配
- $\text{erfc}(\sqrt{\alpha}\, r) / r$ 近似决定 Ewald 实空间性能

### 13.4 Optimization
- Steepest descent: 简单但收敛慢
- Conjugate gradient (CG): 二次收敛
- L-BFGS: 内存受限, 大系统首选
- Simulated annealing: 全局优化, 跳出局部极小
