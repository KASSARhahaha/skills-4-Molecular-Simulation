# Cheatsheet: Quick Decision Tables

分子模拟方法选择的快速参考决策表 (Frenkel & Smit, 3rd ed.)。

## 1. Integrator Selection (积分器选择)

| Goal 目标 | Pick 选择 | Why 原因 |
|---|---|---|
| General MD 通用 MD | Velocity Verlet | Symplectic, simple, low drift 辛、简单、能量漂移小 |
| Constrained bonds 含约束键 | Verlet + SHAKE/LINCS | Allows $2 \times$ larger $\Delta t$ 允许 2 倍时间步 |
| Multi-timescale 多时间步 | r-RESPA | Bonded fast, non-bonded slow 键合快、非键慢 |
| High accuracy 高精度 (rare) | Runge-Kutta | Not for long MD (drift) 不适合长时间 MD |

**Never**: Forward Euler (drift), high-order predictor-corrector (Lyapunov defeats accuracy)。

## 2. Ensemble Selection (系综选择)

| Need 需求 | Ensemble 系综 | Method 方法 |
|---|---|---|
| Equation of state 状态方程 | NPT | NPT-MC or NPT-MD |
| Chemical potential 化学势 | $\mu$VT | $\mu$VT-MC (small molecules) |
| Gas-liquid coexistence 气液共存 | Gibbs | Gibbs ensemble MC |
| Adsorption in pores 孔内吸附 | $\mu$VT | $\mu$VT-MC (PBC + framework) |
| Phase diagram tracking 相图追踪 | NPT + Gibbs-Duhem | Trace coexistence curve |
| Transport properties 输运性质 | NVE | No thermostat disturbance 无恒温器干扰 |

## 3. Thermostat Decision (恒温器决策)

```
Need transport (D, η, λ_T)?
├─ YES → NVE (best) or Nosé-Hoover / Bussi (momentum-preserving)
│        NEVER Andersen, Langevin, Berendsen
└─ NO (static properties only)
    ├─ Equilibration → Berendsen (fast, weak coupling)
    │                  BUT never sample with Berendsen
    └─ Sampling → Bussi (recommended) or Nosé-Hoover
```

## 4. Free Energy Method Selection (自由能方法)

| System 体系 | Reference state 参考 | Method 方法 |
|---|---|---|
| Liquid (any $T, \rho$) | Ideal gas | Density TI (low $\rho$ → target) |
| Solid 固体 | Einstein crystal | Hamiltonian TI |
| Hard sphere 硬球 | Carnahan-Starling | Coupling-softened TI |
| Crystal polymorphs (fcc/hcp) | Either | Lattice switching (Bruce) |
| Polymer 高分子 | Non-self-avoid chain | Rosenbluth sampling |
| Small $\Delta\lambda$ steps | — | FEP (Zwanzig) |
| Many $\lambda$ steps + reuse | — | Histogram reweighting (MBAR) |

## 5. Long-Range Force Method (长程力求解)

| System size $N$ | Method 方法 |
|---|---|
| $N < 10^3$ | Plain Ewald (gold standard) |
| $10^3 \leq N \leq 10^5$ | PME / SPME (FFT accelerated) |
| $N > 10^5$ | FMM ($O(N)$) |
| Isotropic + low precision | Wolf / Reaction field |
| Surface systems | 2D Ewald (special form) |

Ewald 分解: $\dfrac{1}{r} = \dfrac{\text{erfc}(\alpha r)}{r} + \dfrac{\text{erf}(\alpha r)}{r}$

## 6. Reaction Coordinate Quality (反应坐标质量)

Test: 在候选 $q = q^*$ 处取样，释放 MD，统计到达 B 的比例 (committor $P_B$)。

| Committor histogram | Verdict 判定 |
|---|---|
| Single peak at 0.5 (单峰 0.5) | Excellent $q^*$ |
| Broad distribution (宽分布) | Marginal, may still work |
| Bimodal at 0 and 1 (双峰 0/1) | Wrong coordinate, find new $q$ |

## 7. Acceptance Rule Patterns (接受规则)

| Move type 移动类型 | Acceptance 接受概率 |
|---|---|
| NVT displacement | $\min\big(1, e^{-\beta \Delta U}\big)$ |
| NPT volume change | $\min\big(1, e^{-\beta[\Delta U + P \Delta V - N k_B T \ln(V'/V)]}\big)$ |
| $\mu$VT insertion | $\min\big(1, \dfrac{f V}{N+1} e^{-\beta \Delta U}\big)$ |
| $\mu$VT deletion | $\min\big(1, \dfrac{N}{f V} e^{-\beta \Delta U}\big)$ |
| Semigrand species swap | $\min\big(1, e^{-\beta(\Delta\mu + \Delta U)}\big)$ |
| Parallel tempering swap | $\min\big(1, e^{(\beta_i - \beta_j)(U_j - U_i)}\big)$ |
| CBMC regrow | $\min\big(1, W(n)/W(o)\big)$ |

## 8. Common Pitfalls → Fix (常见陷阱与修复)

| Symptom 症状 | Likely cause 原因 | Fix 修复 |
|---|---|---|
| Energy drift in NVE | Non-symplectic integrator or bug | Use Verlet; check force code |
| Wrong $T$ from kinetic energy | Using $v$ instead of $p \cdot v$ | Use $\langle p \cdot v\rangle$ for $T$ |
| Cutoff Coulomb artifacts | Truncating $1/r$ | Use Ewald/PME |
| Berendsen wrong fluctuations | Not a valid ensemble | Use Bussi/Nosé-Hoover |
| Chain insertion rejected | Random in dense fluid | Use CBMC |
| Wrong constraint statistics | Missing Fixman correction | Add $+k_B T \ln |H|^{1/2}$ |
| $g(r)$ noisy at small $r$ | Binning artifacts | Use Borgis force method |
| $S(q)$ negative/oscillatory | Fourier from truncated $g(r)$ | Compute $S(q)$ from $\rho(q)$ |
| Low PT acceptance | $T$ spacing too wide | $\Delta T / T \sim 1/\sqrt{N}$ |

## 9. Key Formulas to Memorize (核心公式)

- **Boltzmann 分布**: $\rho \propto e^{-\beta U}, \quad \beta = 1/(k_B T)$
- **Detailed balance**: $\pi(o)\, \alpha(o \to n)\, \text{acc}(o \to n) = \pi(n)\, \alpha(n \to o)\, \text{acc}(n \to o)$
- **Metropolis**: $\text{acc} = \min\big(1, e^{-\beta \Delta U}\big)$
- **Einstein 扩散系数**: $D = \lim_{t \to \infty} \dfrac{\langle |\Delta \mathbf{r}(t)|^2 \rangle}{2 d t}$
- **Green-Kubo 扩散系数**: $D = \int_0^\infty \langle v_x(0)\, v_x(t) \rangle\, dt$
- **热力学积分 (TI)**: $\Delta F = \int_0^1 \langle \partial H / \partial \lambda \rangle_\lambda\, d\lambda$
- **FEP (Widom/Zwanzig)**: $\Delta F = -k_B T \ln \langle e^{-\beta \Delta U} \rangle$
- **Gibbs-Bogoliubov 不等式**: $\dfrac{\partial^2 F}{\partial \lambda^2} \leq 0$ (linear interpolation)
- **Clausius-Clapeyron**: $\dfrac{dP}{dT} = \dfrac{\Delta h}{T \Delta v}$
- **Fixman 修正**: $\Delta F_{\text{corr}} = k_B T \ln |H|^{1/2}$
- **Jarzynski 等式**: $\langle e^{-\beta W} \rangle = e^{-\beta \Delta F}$
- **响应函数 (Kubo)**: $\chi_{AA}(t) = -\beta\, \langle A(0)\, \dot{A}(t) \rangle\, \Theta(t)$
- **剪切黏度 (Green-Kubo)**: $\eta = \dfrac{1}{V k_B T} \int_0^\infty \langle \sigma_{xy}(0)\, \sigma_{xy}(t) \rangle\, dt$
- **TST 速率**: $k^{\text{TST}} = \langle |\dot{q}| \rangle\, \dfrac{e^{-\beta F(q^*)}}{\int_{-\infty}^{q^*} e^{-\beta F(q)}\, dq}$
- **Bennett-Chandler**: $k_{A \to B} = \kappa\, k^{\text{TST}}$

## 10. Conversion Constants (常用常数)

- $k_B T$ at 300 K $\approx 4.116 \times 10^{-21}$ J $\approx 0.593$ kcal/mol $\approx 2.479$ kJ/mol
- $1$ fs $= 10^{-15}$ s; LJ $\tau^* \approx 1$ ps
- Bjerrum length $\lambda_B$ in water at 300 K $\approx 7.14$ Å
- $1$ atm $\approx 1.013 \times 10^5$ Pa; $1$ bar $= 10^5$ Pa
- LJ reduced units: $T^* = k_B T / \epsilon$, $\rho^* = \rho \sigma^3$, $P^* = P \sigma^3 / \epsilon$
- Ewald 参数: $\alpha \approx 5/L$ (box size $L$), accuracy $10^{-5}$

## 11. Phase Coexistence Relations (相共存关系)

- **Clausius-Clapeyron**: $\dfrac{dP}{dT} = \dfrac{\Delta h}{T \Delta v}$
- **Gibbs 相律**: $f = c - p + 2$ (组分 $c$, 相 $p$, 自由度 $f$)
- **Gibbs-Duhem**: $S\, dT - V\, dP + \sum_i N_i\, d\mu_i = 0$
- **Maxwell 关系**: $\dfrac{\partial^2 F}{\partial T \partial V} = -\dfrac{\partial S}{\partial V} = \dfrac{\partial P}{\partial T}$

## 12. Statistical Mechanics Essentials (统计力学要点)

- **配分函数 (NVT)**: $Q_{NVT} = \dfrac{1}{N! h^{3N}} \int d\mathbf{p}^N d\mathbf{r}^N\, e^{-\beta H}$
- **配分函数 (NVE)**: $\Omega(E) = \int d\mathbf{p}^N d\mathbf{r}^N\, \delta(H - E)$
- **Helmholtz 自由能**: $F = -k_B T \ln Q_{NVT}$
- **熵**: $S = -(\partial F / \partial T)_{V,N}$
- **化学势**: $\mu = (\partial F / \partial N)_{T,V}$
- **压力**: $P = -(\partial F / \partial V)_{T,N}$
