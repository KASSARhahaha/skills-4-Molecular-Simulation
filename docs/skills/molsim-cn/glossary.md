# Glossary: Understanding Molecular Simulation (CN)

Alphabetical reference of key terms (Chinese ↔ English ↔ math symbol).

## A
- **Andersen 恒温器** (Andersen thermostat): 随机重抽速度 → Maxwell-Boltzmann; 破坏动量守恒
- **Andersen 恒压器** (Andersen barostat): V 当作动力学变量, 扩展拉格朗日

## B
- **Bennett-Chandler 方法**: k = k^TST × κ; 静态 × 动态分解
- **Berendsen 恒温器**: 弱耦合, 不对应明确系综, 仅用于平衡
- **Bjerrum 长度** λ_B = e²/(4πε₀ε_r k_BT): 静电能 = 热能的距离
- **Blue-Moon 系综**: 用约束 MD 把体系锁在 q=q*; 硬约束需 Fixman 校正
- **Boltzmann 分布**: ρ ∝ exp(−βU), β = 1/(k_BT)
- **Borgis 力方法**: g(r) 从力场直接积分, 无分箱噪声
- **Bussi 恒温器** (随机速度标度): 现代推荐, 比 Nosé-Hoover 稳健
- **Brownian 动力学**: Langevin 过阻尼极限, m·v̇ 忽略

## C
- **CBMC (构型偏置 MC)**: 用 Rosenbluth 因子校正偏置试探, 链分子采样
- **Clausius-Clapeyron**: dP/dT = Δh/(T Δv), 共存曲线斜率
- **committor q_B(x)**: 从 x 出发首次到 B 的概率; 过渡态 = q_B = 0.5
- **Coulomb 势** u = q²/(4πεr): 必须 Ewald/PME, 不能截断
- **Crooks 涨落定理**: P_F(W)/P_R(−W) = exp(β(W − ΔF))

## D
- **DPD (耗散粒子动力学)**: 软势 + 成对耗散 + 噪声; FDT: σ² = 2k_BT γ^D
- **DSMC (直接模拟 MC)**: Bird, 稀薄气体动力学的 MC 方法
- **detailed balance** (细致平衡): π(o)α(o→n)acc(o→n) = π(n)α(n→o)acc(n→o)
- **delta-correlated noise** ⟨R(0)R(t)⟩ = 2γk_BT δ(t)

## E
- **Einstein 关系**: D = lim_{t→∞} ⟨|Δr(t)|²⟩/(2dt)
- **Einstein 晶体参考**: U(λ) = U(r₀) + (1−λ)[U(r) − U(r₀)] + λΣ αᵢ|rᵢ−r₀,ᵢ|²
- **Ewald 求和**: 1/r = erfc(αr)/r + erf(αr)/r; 短程 + Fourier 分解
- **erfc** (complementary error function): 互补误差函数, Ewald 实空间部分

## F
- **FDT** (涨落耗散定理): ⟨R(0)R(t)⟩ = 2γk_BT δ(t)
- **FEP** (自由能微扰, Zwanzig): ΔF = −k_BT ln ⟨exp(−βΔU)⟩
- **Fixman 校正**: 硬约束统计权重差; +k_BT ln |H|^{1/2}
- **Flyvbjerg-Petersen 块平均**: 误差估计的标准方法, 从平台区读 σ
- **FMM** (快速多极方法): O(N) 长程力, 常数大, 适合 N > 10⁵

## G
- **Gibbs 系综** (Panagiotopoulos): 两盒 + 三类移动; 中密度气-液平衡
- **Gibbs-Bogoliubov 不等式**: 线性插值下 ∂²F/∂λ² ≤ 0
- **Gibbs-Duhem 积分** (Kofke): dP/dT = Δh/(T Δv), 沿共存曲线积分
- **Green-Kubo**: 输运系数 = 平衡时间关联积分, 如 η = (V/k_BT) ∫⟨σ_xy(0)σ_xy(t)⟩dt

## H
- **Hamiltonian TI** (Kirkwood): U(λ) = (1−λ)U_I + λU_II
- **Hamilton 量** H = K + U; Hamilton 方程 q̇ = ∂H/∂p, ṗ = −∂H/∂q
- **hard sphere** (硬球): u(r) = ∞ (r < σ), 0 otherwise; F 解析 (Carnahan-Starling)

## I
- **isotropic ensemble** (等温等压 NPT): 见 NPT-MC

## J
- **Jarzynski 等式**: ⟨exp(−βW)⟩ = exp(−βΔF); 非平衡功 → 自由能

## K
- **k_BT** (Boltzmann 常数 × 温度): 热能尺度; ~ 4.1×10⁻²¹ J at 300 K
- **Kirkwood-Buff 积分** G_ij = 4π∫₀^∞ [g_ij(r) − 1]r²dr
- **Kubo 公式**: χ(t) = −β⟨A(0)B(t)⟩_0 Θ(t); 响应函数

## L
- **Lagrange 乘子** λ: 约束力大小; SHAKE 迭代求解
- **Langevin 动力学**: m·v̇ = −γv − ∇U + R(t); 破坏动量守恒
- **Legendre 变换**: L ↔ H, H(q,p) = Σ p_i q̇_i − L
- **LINCS**: 比 SHAKE 稳定的约束算法
- **Lowe-Andersen 恒温器**: 重置粒子对相对径向速度; 局部动量守恒

## M
- **Markov chain**: 状态链 + 转移概率; 平衡分布 = 不变分布
- **Martyna-Tobias-Klein (MTK)**: 同时恒 T,P 的辛时间可逆组合
- **Maxwell-Boltzmann 分布**: f(v) ∝ v² exp(−mv²/(2k_BT))
- **Metropolis 准则**: acc = min[1, exp(−βΔU)]
- **MPC** (多粒子碰撞动力学, Malevanets-Kapral): 弹道 + 网格碰撞
- **minimum image convention** (最小镜像): PBC 下两粒子距离 = 最近镜像距离

## N
- **Nosé-Hoover 恒温器**: 引入 ξ 变量; 最常用确定性恒温器
- **N_f** (自由度数): d(N−1) − 1 (扣掉总动量 + 能量守恒)

## O
- **Onsager 倒易**: L_ij = L_ji; 来自微观可逆性
- **O(n) 关联算法**: 多级块求和, 长时关联内存友好

## P
- **PBC** (周期性边界条件): r → r mod L; 模拟无限晶格
- **Parrinello-Rahman 恒压器**: h 矩阵当变量, 允许盒形变化
- **PME / SPME** (粒子-网格 Ewald): FFT 加速 Fourier 部分; O(N log N)
- **parallel tempering** (平行回火, PT): 多温度并行 + 构型交换

## Q
- **q_B** (committor): 见 committor

## R
- **RATTLE**: SHAKE 在速度 Verlet 形式下的实现
- **reaction field** (反应场 RFC): 截断外用 ε_RF 连续介质近似
- **replica exchange MD** (REMD): PT 在 MD 中的实现
- **r-RESPA** (多时间步): Trotter 分解按力频率分层
- **Rosenbluth 因子** W = Π_i w_i / k^ℓ; 偏置采样校正

## S
- **SDPD**: 平滑 DPD (Español-Revenga); SPH + 不可逆热力学
- **SHAKE** (Ryckaert-Ciccotti): 顺序迭代约束, 修正键长
- **SRD** (随机旋转动力学 = MPC): 见 MPC
- **semigrand**: 恒 Δμ, N, V, T; 物种身份交换
- **shadow Hamiltonian** (影子 H'): Verlet 守恒的伪 H' = H + O(Δt^{2n})
- **symplectic** (辛): 保持相空间 2-form; Hamilton 流的离散化

## T
- **tail correction** (尾部修正): ∫_{r_c}^∞ 4πr²ρu(r)dr
- **thermal quantity** (热学量): F, G, S, μ; 相空间体积, 非 Boltzmann 平均
- **thermodynamic integration (TI)**: ΔF = ∫₀¹ ⟨∂H/∂λ⟩_λ dλ
- **Tin-foil boundary** (tin-foil): ε' = ∞; 无表面极化, 最常用
- **TPS** (过渡路径采样): 路径马尔可夫链, 不需预设反应坐标
- **Trotter 分解**: exp(iLΔt) ≈ exp(iL_p Δt/2)·exp(iL_r Δt)·exp(iL_p Δt/2)
- **TST** (过渡态理论): k^TST = ⟨δ(q−q*)q̇Θ(q̇)⟩/⟨f_A⟩; 无再交叉

## V
- **Verlet 算法**: r(t+Δt) = 2r(t) − r(t−Δt) + F(t)/m Δt²; 辛, 时间可逆
- **virial pressure** (维里压力): P = ρk_BT + (1/(dV))⟨Σ r_ij·F_ij⟩

## W
- **Wolf method**: 截断 + 阻尼的 Coulomb 近似
- **Widom 试探粒子**: μ_ex = −k_BT ln⟨e^{−βΔU}⟩
- **Wiener-Khinchin**: S(ω) = Fourier[⟨A(0)A(t)⟩]

## X-Y-Z
- **Yukawa potential**: u(r) ∝ exp(−κr)/r; 屏蔽 Coulomb
