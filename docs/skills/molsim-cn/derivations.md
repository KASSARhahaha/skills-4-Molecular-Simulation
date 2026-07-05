# 核心算法推导 (Core Derivations)

> 从第一性原理推导分子模拟的核心算法。每条推导：物理/数学起点 → 关键代数步骤 → 最终公式 → 物理含义与边界。
> 与 `formulas.md` 配套：formulas 是速查，derivations 是教学。
> 来源：Frenkel & Smit《Understanding Molecular Simulation》第三版 + 标准统计力学。

---

## 1. Metropolis 接受准则

**问题**：如何采样 $N$ 粒子构型空间 $\mathbf{r}^N$ 上的 Boltzmann 分布
$\pi(\mathbf{r}^N) \propto e^{-\beta U(\mathbf{r}^N)}$？

### 起点：马尔可夫链设计

构造一条以 $\pi$ 为不动分布的马尔可夫链 $K(o \to n)$。**全局平衡** (global balance) 是充要：

$$
\sum_n \pi(o) K(o \to n) = \sum_n \pi(n) K(n \to o).
$$

全局平衡太自由，难以设计；强化为 **细致平衡** (detailed balance)：

$$
\pi(o) K(o \to n) = \pi(n) K(n \to o). \tag{1}
$$

细致平衡蕴含全局平衡（对 $n$ 求和即得），但反过来不必。

### 分解 $K$：试探 + 接受

把转移拆成两步：以概率 $\alpha(o \to n)$ **试探** 一个新构型，再以概率 $\mathrm{acc}(o \to n)$ **接受**：

$$
K(o \to n) = \alpha(o \to n) \cdot \mathrm{acc}(o \to n).
$$

代入 (1)：

$$
\frac{\mathrm{acc}(o \to n)}{\mathrm{acc}(n \to o)} = \frac{\pi(n)\, \alpha(n \to o)}{\pi(o)\, \alpha(o \to n)} = \frac{e^{-\beta U(n)}\, \alpha(n \to o)}{e^{-\beta U(o)}\, \alpha(o \to n)}. \tag{2}
$$

### 对称试探 → Metropolis 形式

选试探分布对称 $\alpha(o \to n) = \alpha(n \to o)$（如均匀随机位移），(2) 简化为

$$
\frac{\mathrm{acc}(o \to n)}{\mathrm{acc}(n \to o)} = e^{-\beta[U(n) - U(o)]} \equiv e^{-\beta \Delta U}. \tag{3}
$$

为最大化接受率（不让链卡住），取满足 (3) 的最大值：

$$
\boxed{\;\mathrm{acc}(o \to n) = \min\!\left[1,\, e^{-\beta \Delta U}\right]\;}
$$

- $\Delta U \leq 0$（下游）：$e^{-\beta\Delta U} \geq 1$ → 接受率 = 1。
- $\Delta U > 0$（上游）：接受率 = $e^{-\beta\Delta U}$，按概率抽签。

### 非对称试探 → Metropolis–Hastings

若 $\alpha$ 不对称（如偏向低能方向的智能试探），不能用 (3) 简化。**Metropolis–Hastings** 直接用 (2)：

$$
\mathrm{acc}(o \to n) = \min\!\left[1,\; \frac{\alpha(n \to o)\, \pi(n)}{\alpha(o \to n)\, \pi(o)}\right].
$$

### 物理含义

- $\pi \propto e^{-\beta U}$ 是构型在 $NVT$ 下的自然权重。
- MC 不需要算 $Q_N = \int e^{-\beta U}$（无量纲化的"配分函数"），只需相对权重 $e^{-\beta\Delta U}$。
- 对 $\mathbf{r}^N \in \mathbb{R}^{3N}$ 的求积（$N \sim 100$ 时直接求积需 $10^{3N}$ 点）化为按权重的离散平均，**这是 MC 的核心魔力**。

### 边界与陷阱

- **遍历性 (irreducibility)**：细致平衡不够，还需要任何两态互达；否则链被困在亚稳态。
- **自相关**：相邻构型高度相关；有效独立样本数 $\sim T_{\mathrm{sim}} / (2\tau_{\mathrm{int}})$。
- **非各向同性试探**：在临界点附近，位移步长需随方向自适应（如构型偏置 MC）。

---

## 2. Velocity Verlet 积分器

**问题**：数值求解牛顿方程 $m_i \ddot{\mathbf{r}}_i = \mathbf{F}_i = -\nabla_i U$，要求长时间能量稳定、计算便宜。

### 起点：Taylor 展开到三阶

$\mathbf{r}(t)$ 在 $t$ 附近 Taylor 展开：

$$
\mathbf{r}(t \pm \Delta t) = \mathbf{r}(t) \pm \dot{\mathbf{r}}(t)\, \Delta t + \tfrac{1}{2}\, \ddot{\mathbf{r}}(t)\, \Delta t^2 \pm \tfrac{1}{6}\, \mathbf{r}^{(3)}(t)\, \Delta t^3 + \mathcal{O}(\Delta t^4).
$$

两式相加（消除一阶项）：

$$
\mathbf{r}(t + \Delta t) + \mathbf{r}(t - \Delta t) = 2\mathbf{r}(t) + \ddot{\mathbf{r}}(t)\, \Delta t^2 + \mathcal{O}(\Delta t^4). \tag{4}
$$

记 $\mathbf{a} = \ddot{\mathbf{r}} = \mathbf{F}/m$，得到 **Verlet 蛙跳公式**：

$$
\boxed{\;\mathbf{r}(t + \Delta t) = 2\mathbf{r}(t) - \mathbf{r}(t - \Delta t) + \mathbf{a}(t)\, \Delta t^2 + \mathcal{O}(\Delta t^4)\;}
$$

误差 $\mathcal{O}(\Delta t^4)$，但速度没显式出现 — 数值上不方便（速度需用中央差分重构）。

### Velocity Verlet：单步形式

把速度显式纳入。等价地写为

$$
\begin{aligned}
\mathbf{r}(t + \Delta t) &= \mathbf{r}(t) + \mathbf{v}(t)\, \Delta t + \tfrac{1}{2}\, \mathbf{a}(t)\, \Delta t^2, \\[4pt]
\mathbf{v}(t + \tfrac{1}{2}\Delta t) &= \mathbf{v}(t) + \tfrac{1}{2}\, \mathbf{a}(t)\, \Delta t, \\[4pt]
\mathbf{a}(t + \Delta t) &= \mathbf{F}(\mathbf{r}(t + \Delta t))/m, \\[4pt]
\mathbf{v}(t + \Delta t) &= \mathbf{v}(t + \tfrac{1}{2}\Delta t) + \tfrac{1}{2}\, \mathbf{a}(t + \Delta t)\, \Delta t.
\end{aligned}
$$

合并成单步：

$$
\boxed{\;\mathbf{r}(t + \Delta t) = \mathbf{r}(t) + \mathbf{v}(t)\, \Delta t + \tfrac{1}{2}\mathbf{a}(t)\, \Delta t^2,\quad
\mathbf{v}(t + \Delta t) = \mathbf{v}(t) + \tfrac{1}{2}\big[\mathbf{a}(t) + \mathbf{a}(t + \Delta t)\big]\, \Delta t\;}
$$

### 关键性质

1. **时间可逆**：交换 $\Delta t \leftrightarrow -\Delta t$ 公式不变（因为 Taylor 展开只含 $\Delta t^2$）。
2. **辛结构 (symplectic)**：保持相空间体积 $\mathrm{d}\mathbf{r}\,\mathrm{d}\mathbf{p}$ 守恒。
3. **能量稳定**：长时间能量漂移为有界振荡（不是发散），误差在 $10^{-5}$ 量级。
4. **每步一次力计算**：$\mathbf{a}(t + \Delta t)$ 需算 $\mathbf{F}$（昂贵），但 $\mathbf{r}$ 和 $\mathbf{v}$ 更新便宜。

### 为何不用 Runge-Kutta？

RK4 等高阶方法**不是辛的** — 长时间下能量单调漂移，最终系统被"加热"。辛积分器（Verlet 族）的"阴影哈密顿量"原理：数值轨迹是某个扰动哈密顿量的精确轨迹 → 误差有界。

### 边界

- **步长上限**：$\Delta t \lesssim \tau_{\mathrm{fastest}}/10$，否则能量漂移。LJ 流体 $\tau = \sigma\sqrt{m/\epsilon}$ → $\Delta t \sim 0.005\tau$。
- **刚性键**：含 H 原子的键振动周期 $\sim 10$ fs → 必须 $\Delta t \leq 0.5$ fs；用约束 (SHAKE/RATTLE) 把高频键冻结后可放回 2 fs。
- **变步长会破坏辛性** → 不要在 MD 中用自适应步长。

---

## 3. Nosé–Hoover 恒温器

**问题**：MD 默认走 $NVE$；要做 $NVT$ 必须把系统耦合到热浴。

### 起点：扩展拉氏量 (Nosé 1984)

引入一个额外自由度 $s$（"热浴坐标"）和有效质量 $Q$：

$$
\mathcal{L}_{\mathrm{Nosé}} = \sum_{i=1}^{N} \frac{m_i s^2 \dot{\mathbf{r}}_i^2}{2} - U(\mathbf{r}^N) + \frac{Q \dot{s}^2}{2} - N_f k_B T_0 \ln s,
$$

其中 $N_f$ 是自由度数。物理意义：
- $s$ 像一个"时间缩放器"：真实时间 $\mathrm{d}t_{\mathrm{real}} = s\, \mathrm{d}t_{\mathrm{virt}}$；
- $Q\dot{s}^2/2$ 是热浴的"动能"；
- $-N_f k_B T_0 \ln s$ 保证 $s$ 平均有界（否则 $\ln s \to \pm\infty$）。

### 推导运动方程

对 $\mathbf{r}_i$ 变分：

$$
\frac{\mathrm{d}}{\mathrm{d}t}\!\left(m_i s^2 \dot{\mathbf{r}}_i\right) = -\nabla_i U.
$$

对 $s$ 变分：

$$
Q\ddot{s} = \sum_i m_i s \dot{\mathbf{r}}_i^2 - \frac{N_f k_B T_0}{s}.
$$

### Hoover 改写 (1985)

引入"虚拟速度" $\mathbf{v}_i = s\dot{\mathbf{r}}_i$ 和摩擦系数 $\xi = \dot{s}/s$，重写为

$$
\boxed{\;
\begin{aligned}
\dot{\mathbf{r}}_i &= \mathbf{v}_i, \\
\dot{\mathbf{v}}_i &= \frac{\mathbf{F}_i}{m_i} - \xi\, \mathbf{v}_i, \\
\dot{\xi} &= \frac{1}{Q}\!\left[\sum_i m_i v_i^2 - N_f k_B T_0\right]
= \frac{N_f k_B}{Q}(T - T_0).
\end{aligned}
\;}
$$

物理直觉：
- $T > T_0$ → $\dot{\xi} > 0$ → $\xi$ 增大 → 速度衰减 → 系统降温。
- $T < T_0$ → 反之 → 系统升温。
- $Q$ 是反馈增益：太小 → 震荡；太大 → 反应慢。

### 系综正确性

Nosé–Hoover 的不变分布是

$$
\rho(\mathbf{r}^N, \mathbf{p}^N, \xi) \propto
\exp\!\big[-\beta H_{\mathrm{Nosé}}\big],\qquad
H_{\mathrm{Nosé}} = H(\mathbf{r}^N, \mathbf{p}^N) + \tfrac{1}{2}Q\xi^2 + N_f k_B T_0 \ln s.
$$

按 $\xi$ 积分后得到 $NVT$ Boltzmann 分布。

### 边界与陷阱

- **非遍历 (non-ergodic)**：谐振子等可积系统下，$\xi$ 不能访问所有 $(\mathbf{r}, \mathbf{p}, \xi)$ 空间 → 实际温度偏离目标。**修复**：Nosé–Hoover chains (Martyna et al. 1992)，把 $\xi$ 自身耦合到第二个热浴。
- **远平衡**：温度梯度大时 Nosé–Hoover 给出错误的热导率（因为它干扰动量流的自然演化）。**修复**：用反向非平衡 MD (RNEMD) 或周期性扰动。
- **Berendsen 不是真恒温器**：$\dot{T}/(T_0 - T) = 1/\tau$ 是一阶松弛，不采样任何已知系综 — 仅用于平衡期。

---

## 4. Lennard-Jones 尾部校正

**问题**：势能在 $r_c = 2.5\sigma$ 处截断，丢失 $r > r_c$ 的吸引贡献。如何加回？

### 起点：均匀流体近似

假设 $r > r_c$ 处径向分布 $g(r) \approx 1$（即粒子分布如理想气体）。$N-1$ 个粒子在距离 $r$、壳厚 $\mathrm{d}r$ 内的数目为 $4\pi r^2 \rho\, \mathrm{d}r$。单粒子的"漏掉"能量：

$$
u_{\mathrm{tail}} = \int_{r_c}^{\infty} 4\pi r^2 \rho\, U_{\mathrm{LJ}}(r)\, \mathrm{d}r.
$$

代入 $U_{\mathrm{LJ}}(r) = 4\epsilon[(\sigma/r)^{12} - (\sigma/r)^6]$：

$$
u_{\mathrm{tail}} = 16\pi\epsilon\rho \int_{r_c}^{\infty}\!\big(\sigma^{12} r^{-10} - \sigma^6 r^{-4}\big)\, \mathrm{d}r.
$$

逐项积分 $\int_{r_c}^\infty r^{-n}\mathrm{d}r = r_c^{1-n}/(n-1)$：

$$
\begin{aligned}
\int_{r_c}^\infty r^{-10}\,\mathrm{d}r &= \frac{r_c^{-9}}{9}, \\
\int_{r_c}^\infty r^{-4}\,\mathrm{d}r &= \frac{r_c^{-3}}{3}.
\end{aligned}
$$

合并：

$$
\boxed{\;
u_{\mathrm{tail}} = \frac{16\pi\epsilon\rho}{3}\!\left[\frac{1}{3}\frac{\sigma^{12}}{r_c^9} - \frac{\sigma^6}{r_c^3}\right]
= \frac{8\pi}{3}\rho\epsilon\sigma^3\!\left[\frac{1}{3}\!\left(\frac{\sigma}{r_c}\right)^9 - \left(\frac{\sigma}{r_c}\right)^3\right]
\;}
$$

### 压力校正

类似推导对 $\mathbf{F} \cdot \mathbf{r}$：

$$
P_{\mathrm{tail}} = \frac{16\pi}{3}\rho^2 \epsilon\sigma^3\!\left[\frac{2}{3}\!\left(\frac{\sigma}{r_c}\right)^9 - \left(\frac{\sigma}{r_c}\right)^3\right].
$$

注意压力校正比能量校正**大三倍**（系数 $\frac{2}{3}$ vs $\frac{1}{3}$），且符号相反主导。$r_c = 2.5\sigma$、$\rho^* = 0.8$ 时 $P_{\mathrm{tail}}^* \approx -1.0$ — 漏掉会使压力严重高估。

### 边界

- **气液界面**：$g(r)$ 各向异性 → 校正过强。
- **晶体**：晶格位置非均匀 → 校正可能错符号。
- **临界点**：长程涨落被截断 → 需 finite-size scaling。

---

## 5. Ewald 求和

**问题**：周期边界下库仑势 $1/r$ 的条件收敛（$\sum 1/r^3$ 衰减慢），直接求和收敛极慢。

### 起点：屏蔽 + 反屏蔽分解

把 $1/r$ 拆成"短程 + 长程"：

$$
\frac{1}{r} = \underbrace{\frac{\mathrm{erfc}(\alpha r)}{r}}_{\text{短程，实空间}}
+ \underbrace{\frac{\mathrm{erf}(\alpha r)}{r}}_{\text{长程，倒空间}},
$$

其中 $\mathrm{erfc}(x) = 1 - \mathrm{erf}(x)$，$\alpha$ 是高斯屏蔽参数。

- $\mathrm{erfc}(\alpha r)/r$ 在 $r \gtrsim 5/\alpha$ 处指数衰减 → 实空间求和截断到 $r_c$ 即可。
- $\mathrm{erf}(\alpha r)/r$ 在 $r \to 0$ 平滑 → 倒空间收敛快。

### 实空间部分

$$
U_{\mathrm{real}} = \frac{1}{2} \sum_{i<j} q_i q_j \sum_{\mathbf{n}} \frac{\mathrm{erfc}(\alpha |\mathbf{r}_{ij} + \mathbf{n}L|)}{|\mathbf{r}_{ij} + \mathbf{n}L|}.
$$

$\mathbf{n}L$ 是周期镜像位移。截断到 $|\mathbf{r}_{ij} + \mathbf{n}L| < r_c \approx 5/\alpha$。

### 倒空间部分

利用 Poisson 求和：把 $\mathrm{erf}(\alpha r)/r$ 的 Fourier 变换代入。倒空间点 $\mathbf{k} = 2\pi\mathbf{m}/L$（$\mathbf{m} \in \mathbb{Z}^3 \setminus \{0\}$），结构因子

$$
S(\mathbf{k}) = \sum_i q_i e^{i\mathbf{k}\cdot\mathbf{r}_i}.
$$

倒空间贡献：

$$
U_{\mathrm{recip}} = \frac{1}{2\pi V} \sum_{\mathbf{k} \neq 0} \frac{4\pi^2}{k^2}\, e^{-k^2/(4\alpha^2)}\, |S(\mathbf{k})|^2.
$$

### 自相互作用

每个电荷 $q_i$ 与自身的高斯屏蔽云有相互作用，必须扣除：

$$
U_{\mathrm{self}} = -\frac{\alpha}{\sqrt{\pi}} \sum_i q_i^2.
$$

### 总势能

$$
\boxed{\;
U_{\mathrm{Ewald}} = U_{\mathrm{real}} + U_{\mathrm{recip}} + U_{\mathrm{self}}
+ U_{\mathrm{bg}}\;\text{(中性背景项)}
\;}
$$

### 复杂度与 $\alpha$ 选取

实空间求和 $O(N \cdot N_{\text{邻居}}) \sim O(N)$；倒空间求和 $O(N \cdot N_k)$，$N_k \sim k_{\max}^3 \sim \alpha^3 L^3$。平衡两者得 $O(N^{3/2})$。

**PME / SPME** (Hockney-Eastwood / Essmann) 把倒空间求和换成 FFT → $O(N \log N)$。

### 边界

- **非中性体系**：需"中性背景"项 $U_{\mathrm{bg}}$，否则能量发散。
- **金属 vs 真空边界条件**：tin-foil (Surrounding conductor, $\epsilon = \infty$) 比 vacuum ($\epsilon = 1$) 多 $2\pi/3V |\sum q_i \mathbf{r}_i|^2$ 项。
- **$\alpha$ 调参**：$\alpha L \approx 5$–7 通常是好选择；对极性体系需小心。

---

## 6. 自由能微扰 (FEP) 公式

**问题**：从已知势 $U_0$ 计算未知势 $U_1$ 的自由能差 $\Delta F = F_1 - F_0$。

### 起点：配分函数比

$$
\Delta F = -k_B T \ln \frac{Z_1}{Z_0}
= -k_B T \ln \frac{\int e^{-\beta U_1(\mathbf{r}^N)}\, \mathrm{d}\mathbf{r}^N}{\int e^{-\beta U_0(\mathbf{r}^N)}\, \mathrm{d}\mathbf{r}^N}.
$$

### 关键技巧：乘 $1 = e^{-\beta U_0} e^{+\beta U_0}$

$$
\int e^{-\beta U_1}\, \mathrm{d}\mathbf{r}^N
= \int e^{-\beta (U_1 - U_0)}\, e^{-\beta U_0}\, \mathrm{d}\mathbf{r}^N
= Z_0 \cdot \langle e^{-\beta (U_1 - U_0)} \rangle_0.
$$

代入 $\Delta F$：

$$
\boxed{\;
\Delta F = -k_B T \ln \big\langle \exp[-\beta (U_1 - U_0)] \big\rangle_0
\;}
$$

### 重要抽样解释

按 $0$ 系综抽样 $\mathbf{r}^N \sim e^{-\beta U_0}$，对每个构型算 $\Delta U = U_1 - U_0$，最后对 $e^{-\beta\Delta U}$ 求平均。

### 边界（重叠灾难）

如果 $\pi_0$ 和 $\pi_1$ 重叠不足，被积量 $e^{-\beta\Delta U}$ 几乎全由 $\pi_0$ 的尾部（罕见构型）贡献 → 方差爆炸、系统偏差。

**解决**：插入中间 $\lambda$ 点 $U(\lambda) = (1-\lambda)U_0 + \lambda U_1$，分段 FEP 或用 TI / BAR。

---

## 7. 热力学积分 (TI)

**起点**：自由能的全微分

$$
\mathrm{d}F = -S\, \mathrm{d}T - P\, \mathrm{d}V + \mu\, \mathrm{d}N
+ \left\langle \frac{\partial U}{\partial \lambda} \right\rangle_\lambda \mathrm{d}\lambda,
$$

其中 $\lambda$ 是耦合参数。在等温等容（$NVT$）下：

$$
\frac{\partial F}{\partial \lambda}\bigg|_{N,V,T} = \left\langle \frac{\partial U(\lambda)}{\partial \lambda} \right\rangle_\lambda.
$$

从 $\lambda = 0$ 到 $1$ 积分：

$$
\boxed{\;
\Delta F = \int_0^1 \left\langle \frac{\partial U(\lambda)}{\partial \lambda} \right\rangle_\lambda\, \mathrm{d}\lambda
\;}
$$

### 数值实现

- 选 10–20 个 $\lambda$ 点（均匀或软核间距）；
- 每点独立 $NVT$ 模拟，记录 $\partial U/\partial\lambda$ 的平均；
- 梯形/Simpson 积分。

### 软核势 (soft-core)

线性混合 $U(\lambda) = (1-\lambda)U_0 + \lambda U_1$ 在 $\lambda \to 0, 1$ 端点有奇点（粒子重叠 $r \to 0$ 时 $U \to \infty$）。

软核替换：

$$
U(\lambda) = (1-\lambda) U_0 + \lambda U_1 + 4\epsilon\!\left[\frac{1}{(\alpha + (r/\sigma)^6)^2} - \frac{1}{\alpha + (r/\sigma)^6}\right],
$$

避免端点奇点，使 $\partial U/\partial\lambda$ 有限。

### FEP vs TI 谁更优

- FEP：单系综采样足够，但重叠不足时失败。
- TI：每 $\lambda$ 点局部采样，需多点但每点可靠。
- BAR/MBAR：在固定样本下统计最优（最小方差），但需多 $\lambda$ 同时模拟。

---

## 8. Bennett 接受率 (BAR)

**问题**：在 $\lambda_0, \lambda_1$ 各采了 $N_0, N_1$ 个样本，最小方差估计 $\Delta F$？

### 起点：极大似然

Bennett (1976) 把 FEP 的"前向"和"反向"估计结合，引入偏移常数 $C$：

$$
\Delta F = C + k_B T \ln\frac{N_0}{N_1},
$$

其中 $C$ 由自洽方程确定：

$$
\boxed{\;
\sum_{i \in \mathcal{D}_0} \frac{1}{1 + \exp[\beta(\Delta U_i - C)]}
= \sum_{j \in \mathcal{D}_1} \frac{1}{1 + \exp[-\beta(\Delta U_j - C)]}
\;}
$$

- $\Delta U_i = U_1(\mathbf{r}_i) - U_0(\mathbf{r}_i)$，$i$ 来自 0 系综（前向）；
- $\Delta U_j$ 同定义，$j$ 来自 1 系综（反向）。

### 推导思路

 Bennett 把 $F_1 - F_0$ 的混合估计写成

$$
e^{-\beta\Delta F} = \frac{\langle f(\Delta U - C) \rangle_0}{\langle f(C - \Delta U) \rangle_1}\, e^{-\beta C},
$$

其中 $f(x) = 1/(1 + e^{\beta x})$ 是 Fermi 函数。

最小化相对估计的方差，对 $C$ 求导得 BAR 方程。

### 边界

- **MBAR** (Multistate BAR)：扩展到 $K > 2$ 个状态，统一处理 $\lambda$ 网格；
- **采样要求**：两系综的能量分布必须有重叠（否则等式两边各为 0 或 1，无信息）；
- **vs FEP**：BAR 在固定样本下是渐近最小方差，FEP 单方向采样通常较差。

---

## 9. Crooks 等式 与 Jarzynski 等式

**问题**：从非平衡过程（如快速拉伸分子）能提取平衡自由能吗？

### 起点：微观可逆性

哈密顿 $H(\mathbf{r}, \mathbf{p}; \lambda)$ 依赖外参 $\lambda(t)$。前向过程 $\lambda: 0 \to 1$，逆向 $\lambda: 1 \to 0$。微观可逆性 → 前向轨迹 $\Gamma$ 与反向轨迹 $\tilde{\Gamma}$ 的概率比：

$$
\frac{\mathcal{P}_F[\Gamma]}{\mathcal{P}_R[\tilde{\Gamma}]} = \exp\!\big[\beta(W[\Gamma] - \Delta F)\big],
$$

其中 $W[\Gamma] = \int_0^1 \partial H/\partial\lambda\, \mathrm{d}\lambda$ 是沿轨迹做的功。

### Crooks 等式 (1999)

对耗散功 $W_{\mathrm{diss}} = W - \Delta F$ 的分布：

$$
\boxed{\;
\frac{P_F(W)}{P_R(-W)} = \exp\!\big[\beta(W - \Delta F)\big]
\;}
$$

两分布的交点 $W^*$ 满足 $P_F(W^*) = P_R(-W^*)$ → $W^* = \Delta F$。

### Jarzynski 等式 (1997)

把 Crooks 等式对 $P_F(W)$ 求期望：

$$
\big\langle e^{-\beta W} \big\rangle_F
= \int P_F(W) e^{-\beta W}\, \mathrm{d}W
= \int P_R(-W) e^{-\beta\Delta F}\, \mathrm{d}W
= e^{-\beta \Delta F}.
$$

得 **Jarzynski 等式**：

$$
\boxed{\;
\big\langle e^{-\beta W} \big\rangle = e^{-\beta \Delta F}
\quad\Longleftrightarrow\quad
\Delta F = -k_B T \ln \big\langle e^{-\beta W} \big\rangle
\;}
$$

### 实际陷阱：稀有样本偏差

被积量 $e^{-\beta W}$ 由稀有小 $W$ 轨迹主导。**第二定律 $\langle W\rangle \geq \Delta F$** 严格成立，但样本估计 $\widehat{\Delta F} = -k_B T \ln \langle e^{-\beta W} \rangle$ 偏高（Jensen 不等式）。

**修复**：
- 大量重复（典型 100–1000 次拉伸）；
- 用 BAR/Crooks 交叉点法（更稳健）；
- 用二阶累积量展开 (Hummer 2001)。

---

## 10. Green-Kubo 关系

**问题**：从平衡 MD 轨迹的涨落计算输运系数（黏度 $\eta$、热导 $\kappa$、扩散 $D$）。

### 起点：线性响应理论

对可观测量 $A$ 加微扰 $f(t)$，扰动后期望值（Kubo 公式，附录 F）：

$$
\delta \langle A(t) \rangle = \beta \int_{-\infty}^{t}\!\mathrm{d}t'\, f(t')\, C_{AJ}(t - t'),
$$

其中 $C_{AJ}(\tau) = \langle A(0)\, J(\tau) \rangle_{\mathrm{eq}}$ 是平衡相关函数，$J = \dot{A}$。

### 黏度

剪切应力 $P_{xy}$ 与流动 $v_x(y)$ 的耦合：取 $A = \sum_i m_i v_{i,x} y_i$（横流动量），$J = \dot{A} = V P_{xy}$。

线性响应 + 周期性剪切外场 $f(t) = F_x \delta(t)$ → 流动响应 → 黏度：

$$
\boxed{\;
\eta = \frac{\beta}{V}\!\int_0^\infty \langle P_{xy}(0)\, P_{xy}(t) \rangle\, \mathrm{d}t
\;}
$$

### 热导

类似地，能流 $J_q$ 与温度梯度的耦合：

$$
\kappa = \frac{\beta}{V}\!\int_0^\infty \langle J_q(0)\, J_q(t) \rangle\, \mathrm{d}t.
$$

### 自扩散

特殊情形：速度自相关积分

$$
D = \int_0^\infty \langle v_x(0)\, v_x(t) \rangle\, \mathrm{d}t = \frac{1}{3}\!\int_0^\infty \langle \mathbf{v}(0) \cdot \mathbf{v}(t) \rangle\, \mathrm{d}t.
$$

等价的 Einstein 关系：

$$
D = \lim_{t \to \infty} \frac{1}{6t} \big\langle |\mathbf{r}(t) - \mathbf{r}(0)|^2 \big\rangle.
$$

### 实际陷阱

- **长时尾**：$C(t) \sim t^{-d/2}$（$d$ 维），积分慢收敛；
- **有限尺寸**：$D(L) = D(\infty) - a/L$ 是经典校正；
- **统计噪声**：相关函数长时被噪声淹没 → 用多个时间原点平均 + 阻塞平均。

---

## 推导索引

| 推导 | 起点 | 关键技巧 | 公式 |
|---|---|---|---|
| Metropolis | 细致平衡 | 试探 + 接受分解 | $\mathrm{acc} = \min(1, e^{-\beta\Delta U})$ |
| Velocity Verlet | Taylor 到 $\Delta t^3$ | 前后两式相加 | $r(t+\Delta t) = 2r(t) - r(t-\Delta t) + a\Delta t^2$ |
| Nosé–Hoover | 扩展拉氏量 | 引入 $s, \xi$ | $\dot{\xi} = (N_f k_B/Q)(T - T_0)$ |
| LJ 尾部校正 | 均匀流体 $g(r)=1$ | $4\pi r^2 \rho\, U(r)\, dr$ 积分 | $u_{\mathrm{tail}} = \frac{8\pi}{3}\rho\epsilon\sigma^3[\frac{1}{3}(\sigma/r_c)^9 - (\sigma/r_c)^3]$ |
| Ewald 求和 | erf + erfc 分裂 | Fourier 变换 + 自能扣除 | $U = U_{\mathrm{real}} + U_{\mathrm{recip}} + U_{\mathrm{self}}$ |
| FEP | $Z_1/Z_0$ | 乘 $e^{-\beta U_0} e^{+\beta U_0}$ | $\Delta F = -k_BT \ln\langle e^{-\beta\Delta U}\rangle_0$ |
| TI | $\mathrm{d}F$ 全微分 | $\partial F/\partial\lambda$ | $\Delta F = \int_0^1 \langle \partial U/\partial\lambda\rangle \mathrm{d}\lambda$ |
| BAR | 极大似然 + Fermi 函数 | 自洽方程 | $C + k_BT \ln(N_0/N_1) = \Delta F$ |
| Crooks/Jarzynski | 微观可逆性 | 轨迹概率比 | $\langle e^{-\beta W}\rangle = e^{-\beta\Delta F}$ |
| Green-Kubo | 线性响应 | 平衡相关函数积分 | $\eta = (\beta/V)\!\int\langle P_{xy}(0)P_{xy}(t)\rangle\mathrm{d}t$ |

---

## 推荐学习路径

1. **入门**：先读 Metropolis（§1）和 Velocity Verlet（§2），这两个是 MC/MD 的全部基础。
2. **进阶**：Nosé–Hoover（§3）+ LJ 尾部校正（§4）→ 能写出生产级 NVT MD。
3. **自由能**：FEP（§6）→ TI（§7）→ BAR（§8）→ Crooks/Jarzynski（§9）。这条线是现代自由能计算的完整谱系。
4. **高级**：Ewald（§5）+ Green-Kubo（§10）→ 处理带电体系和输运系数。
