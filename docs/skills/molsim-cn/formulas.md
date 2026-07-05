# 核心公式速查 (Core Formulas)

> 手工策展的分子模拟核心方程。每条包含：公式 + 物理含义 + 边界条件 + 单位/约定。
> 来源：Frenkel & Smit《Understanding Molecular Simulation》第三版。

---

## 1. 统计力学基础 (Ch. 2)

### 1.1 正则配分函数 (Canonical partition function)

$$
Z(N, V, T) = \frac{1}{N!\, h^{3N}} \int \mathrm{d}\mathbf{r}^N\, \mathrm{d}\mathbf{p}^N\;
\exp\!\big[-\beta H(\mathbf{r}^N, \mathbf{p}^N)\big]
$$

- **含义**：$NVT$ 系综的所有微观状态按 Boltzmann 权重的总和。
- **经典极限**：动能 $\mathbf{p}$ 部分解析积分 → 只剩构型积分 $Q_N$。
- **单位**：$\beta = 1/(k_B T)$；$h$ 是 Planck 常数（仅在用量子修正时显式出现）。

### 1.2 构型积分 (Configurational integral)

$$
Q_N(V, T) = \int \mathrm{d}\mathbf{r}^N\, \exp\!\big[-\beta U(\mathbf{r}^N)\big]
$$

- **用途**：所有 NVT MC 的核心求积对象。
- **灾难性**：$N \sim 100$ 时直接求积不可行（$10^{3N}$ 维），必须用重要抽样。

### 1.3 系综平均 (Ensemble average)

$$
\langle A \rangle_{NVT} = \frac{\displaystyle\int \mathrm{d}\mathbf{r}^N\,
A(\mathbf{r}^N)\, e^{-\beta U(\mathbf{r}^N)}}{Q_N}
$$

- **重要抽样实现**：$\langle A \rangle \approx \frac{1}{L}\sum_{i=1}^{L} A(\mathbf{r}_i^N)$，其中 $\mathbf{r}_i^N$ 按 $\pi(\mathbf{r}^N) \propto e^{-\beta U}$ 采样。

### 1.4 能量涨落 (Energy fluctuation)

$$
\langle (\delta E)^2 \rangle = k_B T^2\, C_V
\quad\Rightarrow\quad
C_V = \frac{\langle E^2\rangle - \langle E\rangle^2}{k_B T^2}
$$

- **用途**：从 MC/MD 轨迹直接算热容，无需数值微分。
- **边界**：仅 $NVT$ 严格成立；其他系综需用对应的响应函数公式。

### 1.5 压力涨落 → 体积弹性模量

$$
\kappa_T = \frac{1}{k_B T}\, \frac{\langle V^2\rangle - \langle V\rangle^2}{\langle V\rangle}
\quad (NPT)
$$

---

## 2. Monte Carlo 核心 (Ch. 3, 6)

### 2.1 Metropolis 接受准则

$$
\mathrm{acc}(o \to n) = \min\!\left[1,\; \exp\!\big(-\beta [U(n) - U(o)]\big)\right]
$$

- **前提**：试探概率对称（$\alpha(o\to n) = \alpha(n\to o)$）。
- **非对称试探**：用 Metropolis–Hastings：
  $\mathrm{acc}(o\to n) = \min\!\left[1,\; \dfrac{\alpha(n\to o)\, \pi(n)}{\alpha(o\to n)\, \pi(o)}\right]$。

### 2.2 细致平衡 (Detailed balance)

$$
\pi(o)\, K(o\to n) = \pi(n)\, K(n\to o)
$$

- $K = \alpha \cdot \mathrm{acc}$：转移概率 = 试探 × 接受。
- **必要性**：保证 $\pi$ 是马尔可夫链的不动分布；非细致平衡的采样器会偏。

### 2.3 NPT 接受准则

$$
\mathrm{acc}(o \to n) = \min\!\left[1,\;
\exp\!\big(-\beta[\Delta U + P\Delta V] - \beta^{-1} N \ln(V_n/V_o)\big)\right]
$$

- **附加项**：体积改变的熵贡献 $N\ln(V_n/V_o)$。
- **常见 bug**：漏掉 $\ln V$ 项导致 NPT 漂移。

### 2.4 Gibbs 系综 (两箱无化学势)

$$
\mathrm{acc}(\text{particle swap}) = \min\!\left[1,\;
\frac{N_2\, V_1}{N_1\, V_2}\,
\exp\!\big(-\beta[\Delta U_1 + \Delta U_2]\big)\right]
$$

- **用途**：直接计算气液相平衡，无需迭代化学势。

### 2.5 化学势插入 (Widom test particle)

$$
\mu_{\mathrm{ex}} = -k_B T \ln \langle \exp(-\beta \Delta U) \rangle_{NVT}
$$

- **限制**：高密度下被积量 $\langle e^{-\beta\Delta U}\rangle$ 几乎全由稀有构型贡献 → 方差爆炸。

---

## 3. 分子动力学核心 (Ch. 4, 7)

### 3.1 Velocity Verlet 积分器

$$
\begin{aligned}
\mathbf{r}(t + \Delta t) &= \mathbf{r}(t) + \mathbf{v}(t)\,\Delta t
+ \frac{1}{2}\,\mathbf{a}(t)\,\Delta t^2 \\[4pt]
\mathbf{v}(t + \Delta t) &= \mathbf{v}(t) + \frac{1}{2}\big[\mathbf{a}(t)
+ \mathbf{a}(t + \Delta t)\big]\,\Delta t
\end{aligned}
$$

- **优点**：时间可逆、辛结构近似保持、仅需一次额外力计算。
- **典型步长**：$\Delta t \sim 0.005\tau$（约 5–10 fs for Ar）；$\tau = \sigma\sqrt{m/\epsilon}$。

### 3.2 Lennard-Jones 势

$$
U_{\mathrm{LJ}}(r) = 4\epsilon\!\left[
\left(\frac{\sigma}{r}\right)^{12}
- \left(\frac{\sigma}{r}\right)^{6}
\right]
$$

- **能量零点**：$r = \sigma$ 时 $U = 0$；最低点 $r_{\min} = 2^{1/6}\sigma$，$U_{\min} = -\epsilon$。
- **力**：$\mathbf{F} = -\nabla U = \dfrac{24\epsilon}{r}\!\left[2\!\left(\dfrac{\sigma}{r}\right)^{12}
- \left(\dfrac{\sigma}{r}\right)^6\right]\hat{\mathbf{r}}$。
- **截断**：$r_c = 2.5\sigma$ 是工业默认；临界点附近必须加尾部校正。

### 3.3 LJ 尾部校正

$$
u_{\mathrm{tail}} = \frac{8\pi}{3}\, \rho\, \epsilon\, \sigma^3\!\left[
\tfrac{1}{3}\!\left(\frac{\sigma}{r_c}\right)^{9}
- \left(\frac{\sigma}{r_c}\right)^{3}
\right]
$$

- **假设**：$r > r_c$ 处 $\rho(r) \approx \rho_{\text{bulk}}$（均匀流体）。
- **失效**：高各向异性体系（界面、孔道、晶体）。

### 3.4 Nosé–Hoover 恒温器

$$
\dot{\mathbf{p}}_i = \mathbf{F}_i - \xi\, \mathbf{p}_i,
\qquad
\dot{\xi} = \frac{1}{Q}\!\left[\sum_i \frac{p_i^2}{m_i} - N_f k_B T\right]
$$

- $Q$：恒温器"质量"；$N_f$：自由度数（通常 $3N - 3$）。
- **陷阱**：$Q$ 太小 → 温度震荡；$Q$ 太大 → 实际不变温。
- **非遍历体系**：谐振子下 Nosé–Hoover 不能遍历相空间，需 Nosé–Hoover chains。

### 3.5 Berendsen 恒温器（非平衡专用）

$$
\frac{\mathrm{d}T}{\mathrm{d}t} = \frac{T_0 - T}{\tau_T}
\quad\Rightarrow\quad
\lambda = \sqrt{1 + \frac{\Delta t}{\tau_T}\!\left(\frac{T_0}{T} - 1\right)}
$$

- **用途**：快速平衡（前 10–100 ps）；**不用于产出**（不采样任何正确系综）。

---

## 4. 自由能计算 (Ch. 8–10)

### 4.1 自由能微扰 (FEP)

$$
\Delta F = -k_B T \ln \langle \exp[-\beta (U_1 - U_0)] \rangle_0
$$

- **要求**：$\lambda = 0$ 和 $\lambda = 1$ 的构型分布重叠；否则需中间 $\lambda$ 点。
- **重叠规则**：若 $\langle\cdot\rangle_0$ 和 $\langle\cdot\rangle_1$ 的能量分布峰几乎不相交，必须分段。

### 4.2 热力学积分 (TI)

$$
\Delta F = \int_{\lambda=0}^{1} \mathrm{d}\lambda\;
\left\langle \frac{\partial U(\lambda)}{\partial \lambda} \right\rangle_\lambda
$$

- **数值实现**：在 $\lambda_k \in \{0, 0.05, \ldots, 1\}$ 各跑独立模拟，梯形/Simpson 积分。
- **优于 FEP**：对每个 $\lambda$ 局部平均，不要求全局重叠。

### 4.3 Bennett 接受率 (BAR)

$$
\sum_i \frac{1}{1 + \exp[\beta(U_1 - U_0) + C]}
= \sum_j \frac{1}{1 + \exp[-\beta(U_1 - U_0) - C]}
$$

迭代解 $C$，然后 $\Delta F = C + k_B T \ln(N_0/N_1)$。
- **最小方差**：在给定采样下 BAR 是 FEP 的最优无偏估计。
- **MBAR**：多状态扩展，处理 $\lambda$ 网格时推荐。

### 4.4 Jarzynski 等式（非平衡）

$$
e^{-\beta \Delta F} = \big\langle e^{-\beta W} \big\rangle_{\mathrm{neq}}
$$

- **陷阱**：$W$ 的尾部分布决定收敛；有限样本会高估 $\Delta F$（Crooks 修正在 [Eqs. 4.5–4.6]）。

### 4.5 RDF 与结构因子

$$
g(r) = \frac{V}{N^2}\!\left\langle \sum_{i\neq j} \delta(\mathbf{r}_i - \mathbf{r}_j - \mathbf{r}) \right\rangle,
\qquad
S(k) = 1 + \rho \int \mathrm{d}\mathbf{r}\; [g(r) - 1]\, e^{i\mathbf{k}\cdot\mathbf{r}}
$$

---

## 5. 长程相互作用 (Ch. 11)

### 5.1 Ewald 求和

$$
U_{\mathrm{Coul}} = \frac{1}{2}\sum_{i<j} q_i q_j
\left[
\underbrace{\sum_{\mathbf{n}} \frac{\mathrm{erfc}(\alpha\,|\mathbf{r}_{ij}+\mathbf{n}L|)}{|\mathbf{r}_{ij}+\mathbf{n}L|}}_{\text{real-space}}
+ \underbrace{\frac{1}{\pi V}\sum_{\mathbf{k}\neq 0} \frac{4\pi^2}{k^2}\, e^{-k^2/(4\alpha^2)}\, |\rho(\mathbf{k})|^2}_{\text{reciprocal-space}}
\right] - U_{\mathrm{self}}
$$

- $\alpha$：精度参数（典型 $\alpha L \approx 5$–7）。
- **复杂度**：$O(N^{3/2})$；用 PME / SPME 可降到 $O(N\log N)$。

### 5.2 反作用场 (Reaction field, 截断替代)

$$
U_{\mathrm{rf}}(r < r_c) = -\frac{1}{2}\, \mu_i\cdot
\frac{2(\epsilon_{\mathrm{rf}}-1)}{2\epsilon_{\mathrm{rf}}+1}\,
\frac{\mathbf{r}_{ij}}{r_{ij}^3}
$$

- **便宜**：$O(N)$；中性极性液体（TIP3P 水）足够精确。
- **失败**：高电荷体系（离子液体、带电表面）必须用 Ewald。

---

## 6. 误差估计与采样质量

### 6.1 阻塞平均 (Block averaging)

$$
\sigma_{\bar{A}}^2 = \frac{1}{N_b(N_b - 1)} \sum_{b=1}^{N_b} (\bar{A}_b - \bar{A})^2
\cdot \tau_{\mathrm{int}}\;/\;\Delta t_{\mathrm{block}}
$$

- **要点**：相邻帧高度相关 → 必须按统计独立时长分块，否则严重低估误差。

### 6.2 集成自相关时间

$$
\tau_{\mathrm{int}} = 1 + 2\sum_{t=1}^{\infty} \rho(t),
\qquad
\sigma^2_{\bar{A}} \approx \frac{2\tau_{\mathrm{int}}}{T_{\mathrm{sim}}}\, \mathrm{Var}(A)
$$

- **规则**：有效独立样本数 $\approx T_{\mathrm{sim}} / (2\tau_{\mathrm{int}})$。

---

## 7. 约化单位约定 (Reduced units)

| 量 | 约化形式 | Ar 60 K 对应 SI |
|---|---|---|
| 长度 | $r^* = r/\sigma$ | $\sigma = 3.40$ Å |
| 能量 | $U^* = U/\epsilon$ | $\epsilon/k_B = 119.8$ K |
| 时间 | $t^* = t/\tau$，$\tau = \sigma\sqrt{m/\epsilon}$ | $\tau \approx 2.16$ ps |
| 温度 | $T^* = k_B T/\epsilon$ | $T^* = 0.5$ ↔ 60 K |
| 密度 | $\rho^* = \rho\sigma^3$ | $\rho^* = 0.8442$ ↔ 1350 kg/m³ |
| 压力 | $P^* = P\sigma^3/\epsilon$ | $P^* = 1$ ↔ 41 MPa |

**为什么用约化单位**：同一模拟结果可对应多种物质（Ar 60 K ≡ Kr 84 K ≡ Xe 112 K 的对应状态点）。

---

## 8. 算法常数与经验值

| 参数 | 推荐值 | 说明 |
|---|---|---|
| LJ 截断 $r_c$ | $2.5\sigma$–$3.0\sigma$ | 临界点附近需更大 |
| MD 步长 $\Delta t$ | $0.005\tau$（Verlet） | 刚性键可用 $0.01\tau$ |
| Ewald $\alpha$ | $5/L$（实空间截断到 $r_c = L/2$） | 平衡实/倒空间开销 |
| MC 位移步长 | 调到接受率 30–50% | 自适应在前 1000 步用，之后冻结 |
| 平衡所需步数 | $\geq 10\tau_{\mathrm{int}}$ | 经验值；测量是唯一可靠办法 |
