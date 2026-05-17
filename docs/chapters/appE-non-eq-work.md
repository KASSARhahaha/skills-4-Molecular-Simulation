# 非平衡功与细致平衡

第 8.7 节中描述的自由能差与（非平衡）功之间的关系，对于所有满足马尔可夫性且每一步都满足细致平衡的方案或运动方程都成立。

为了理解这意味着什么，考虑一个通过改变控制参数 $\lambda$ 来做功的过程。例如，$\lambda$ 可以是系统的体积，也可以是系统哈密顿量中的一个参数。我们将 $\lambda$ 的初始值/最终值记为 $\lambda_0 / \lambda_K$（我们使用下标 $K$ 以保持后续记号的一致性）。

将 $\lambda$ 从 $\lambda_0$ 改变到 $\lambda_K$ 的方案可以分解为两类基本步骤：在第一类步骤中，系统的所有相空间坐标（$\Gamma$）保持固定，$\lambda$ 改变了 $\Delta \lambda_i$，其中 $i$ 标记步骤编号；如果共有 $K$ 个这样的步骤，则 $i = \{1, 2, \cdots, K\}$。注意，仅改变 $\lambda$ 的步骤是确定性的。然而，方案中还有第二类步骤，即 $\lambda$ 保持不变，但系统通过其自然动力学进行演化，与恒温器交换能量。恒定 $\lambda$ 下演化的例子包括一系列一个或多个蒙特卡罗移动，或恒定温度 MD 模拟中的一个或多个时间步长。^1 在 Crooks ^[\ref{389}] 的术语中，我们将与"热库"（即恒温器）交换的能量记为 $Q$：它可以被解释为热库传递给系统的热量。由于细致平衡，系统在恒定 $\lambda$（比如 $\lambda_i$）下从相空间坐标 $\Gamma$ 演化到 $\Gamma'$ 的概率与逆向移动概率之比由下式给出：

$$
\frac{P(\Gamma \to \Gamma'; \lambda_i)}{P(\Gamma' \to \Gamma; \lambda_i)} = e^{-\beta \Delta E(\lambda_i)},
$$

其中 $\Delta E(\lambda_i) \equiv E(\Gamma'; \lambda_i) - E(\Gamma; \lambda_i)$。

当我们在恒定 $\Gamma$ 下改变 $\lambda$ 时，我们对系统做功。我们将 $\Gamma$ 固定时 $\lambda$ 从 $\lambda_{i-1}$ 变到 $\lambda_i$ 对应的功记为 $w_i$。当 $\lambda$ 从 $\lambda_0$ 增加到 $\lambda_K$ 时，对系统所做的总功 $W$ 等于

$$
W = \sum_{i=1}^{K} w_i.
$$

如果系统的时间演化满足马尔可夫性，我们可以写出系统从 $\lambda_0$ 处的 $\Gamma_0$ 演化到 $\lambda_K$ 处的 $\Gamma_K$ 的概率为

$$
\prod_{i=0}^{K-1} \left[ P(\Gamma_i \to \Gamma_{i+1}; \lambda_i) \times 1 \right].
$$

我们在每一步中包含了一个因子一，以表示在恒定 $\Gamma$ 下改变 $\lambda$ 是确定性的。在接下来的讨论中，我们省略这个平凡的因子。我们也可以写出系统沿相同路径从 $\lambda_K$ 处的 $\Gamma_K$ 逆向演化到 $\lambda_0$ 处的 $\Gamma_0$ 的概率：

$$
\prod_{i=K-1}^{0} P(\Gamma_{i+1} \to \Gamma_i; \lambda_i).
$$

由于式 (E.1.1)，我们有

$$
\prod_{i=0}^{K-1} P(\Gamma_i \to \Gamma_{i+1}; \lambda_i) = \prod_{i=K-1}^{0} \left[ P(\Gamma_{i+1} \to \Gamma_i; \lambda_i) \, e^{-\beta \Delta E(\lambda_i)} \right].
$$

注意 $\prod_{i=K-1}^{0} e^{-\beta \Delta E(\lambda_i)}$ 等于 $e^{-\beta Q(\{\Gamma\})}$，其中 $Q$ 是状态序列 $\{\Gamma\} \equiv \Gamma_0, \Gamma_1, \cdots, \Gamma_K$ 的总能量——即从热库传递给系统的能量。注意 $Q$ 不等于 $E(\Gamma_K; \lambda_K) - E(\Gamma_0; \lambda_0)$，因为系统的能量也因做功而改变：

$$
E(\Gamma_K; \lambda_K) - E(\Gamma_0; \lambda_0) = W(\{\Gamma\}) + Q(\{\Gamma\}),
$$

这可以看作热力学第一定律的微观版本。特别需要注意的是，$W$ 和 $Q$ 都依赖于路径，但它们的和与路径无关。如果我们从玻尔兹曼分布中采样初始条件，则可以将 $e^{-\beta W}$ 的平均值表示为

$$
\langle e^{-\beta W} \rangle = \sum_{\Gamma_0, \cdots, \Gamma_K} P_B(\Gamma_0, \lambda_0) \prod_{i=0}^{K-1} P(\Gamma_i \to \Gamma_{i+1}; \lambda_i) \, e^{-\beta W(\{\Gamma\})},
$$

其中 $P_B(\Gamma_0, \lambda_0) = \exp\left(-\beta [E(\Gamma_0; \lambda_0) - F(\lambda_0)]\right)$，而Helmholtz自由能 $F$ 如通常一样由下式给出：

$$
\beta F(\lambda) = -\ln \sum_{\Gamma} e^{-\beta E(\Gamma; \lambda)}.
$$

利用式 (E.1.2)，我们可以写出

$$
\frac{\prod_{i=0}^{K-1} P(\Gamma_i \to \Gamma_{i+1}; \lambda_i)}{\prod_{i=K-1}^{0} P(\Gamma_{i+1} \to \Gamma_i; \lambda_i)} = e^{-\beta Q(\{\Gamma\})}.
$$

或者，利用玻尔兹曼权重的定义：

$$
\frac{P_B(\Gamma_0, \lambda_0) \prod_{i=0}^{K-1} P(\Gamma_i \to \Gamma_{i+1}; \lambda_i)}{P_B(\Gamma_K, \lambda_K) \prod_{i=K-1}^{0} P(\Gamma_{i+1} \to \Gamma_i; \lambda_i)} = e^{-\beta Q(\{\Gamma\})} \, e^{+\beta[\Delta E - \Delta F]},
$$

其中 $\Delta E \equiv E(\Gamma_K; \lambda_K) - E(\Gamma_0; \lambda_0)$，$\Delta F \equiv F(\lambda_K) - F(\lambda_0)$。由此我们可以写出：

$$
\begin{aligned}
P_B(\Gamma_0, \lambda_0) \prod_{i=0}^{K-1} P(\Gamma_i \to \Gamma_{i+1}; \lambda_i) \, e^{-\beta W(\{\Gamma\})}
&= P_B(\Gamma_K, \lambda_K) \prod_{i=K-1}^{0} P(\Gamma_{i+1} \to \Gamma_i; \lambda_i) \, e^{-\beta[Q(\{\Gamma\}) + W(\{\Gamma\}) - \Delta E]} \, e^{-\beta \Delta F} \\
&= P_B(\Gamma_K, \lambda_K) \prod_{i=K-1}^{0} P(\Gamma_{i+1} \to \Gamma_i; \lambda_i) \, e^{-\beta \Delta F},
\end{aligned}
$$

其中我们利用了式 (E.1.3)。由此可得，式 (E.1.4) 可以写为

$$
\langle e^{-\beta W} \rangle = e^{-\beta \Delta F} \sum_{\Gamma_0, \cdots, \Gamma_K} \prod_{i=K-1}^{0} P(\Gamma_{i+1} \to \Gamma_i; \lambda_i) = e^{-\beta \Delta F},
$$

其中最后一个等号源于所有转移概率都是归一化的这一事实。这就是 Jarzynski 等式：对非平衡功的指数取平均，可以得到自由能差 $\Delta F$。

^1 该论证可以推广到系统也与热库交换体积或粒子的情况，但此处我们考虑最简单的情况。