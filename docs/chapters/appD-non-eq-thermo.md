# 非平衡热力学

## 熵产生

当一个系统从非平衡态弛豫时，会产生熵。由昂萨格 [[55,56]](references.md#ref-55)建立的非平衡热力学（不可逆热力学）确定了熵产生的不同贡献之间的关系。因此，我们需要一个熵产生的表达式作为出发点。熵产生的“经典”推导可以在 De Groot 和 Mazur 所著的《Non-equilibrium Thermodynamics》[[57]](references.md#ref-57)一书中找到。然而，文献[[57]](references.md#ref-57)中熵产生表达式的推导虽然完整，但有些令人生畏。在这里，我们选择了一条不同的路线：它不够完整，但更加简洁，也更容易理解其中的物理本质。

让我们首先简要考察将一个封闭系统的两个子系统（I 和 II）接触时所产生的熵产生。这两个子系统可以交换能量和粒子。我们忽略介质可能发生流动这一事实——这是为保持问题简洁所付出的代价。实际上，粘性流是一个重要的非平衡现象，将其纳入熵产生表达式的方法在文献[[57]](references.md#ref-57)中有详细描述。

子系统 I 的熵$S$由于一个无穷小量交换而发生的变化为：

$$
\mathrm{d}S^I = \frac{1}{T^I}\mathrm{d}U^I - \sum_i \frac{\mu_i^I}{T^I}\mathrm{d}N_i^I .

\tag{D.1.1}
$$

由于整个系统是封闭的，我们有$\mathrm{d}U^I = -\mathrm{d}U^{II}$，$\mathrm{d}N_i^{II} = -\mathrm{d}N_i^I$。因此，熵的总变化为

$$
\mathrm{d}S_{\mathrm{total}} = \mathrm{d}S^I + \mathrm{d}S^{II} = \left(\frac{1}{T^I} - \frac{1}{T^{II}}\right)\mathrm{d}U^I - \sum_i \left(\frac{\mu_i^I}{T^I} - \frac{\mu_i^{II}}{T^{II}}\right)\mathrm{d}N_i^I .

\tag{D.1.2}
$$

该方程表明，如果系统中存在能量流且伴有$1/T$的梯度，或者存在粒子流且伴有$\mu_i/T$的梯度，则会产生熵。现在考虑一个厚度为$\mathrm{d}x$、截面积为$A$的薄板，热量在温度$T(x)$处流入系统，在温度$T(x+\mathrm{d}x)$处流出。我们用$\dot{Q}$表示单位时间内传输的热量。热流（heat flux）则为$j_q = \dot{Q}/A$。

类似地，流过系统的第$i$种粒子的数量为$\dot{N}_i$，第$i$种物种的粒子流（particle flux）为$j_i = \dot{N}_i/A$。注意，我们假设系统整体处于静止状态，因此所有粒子流的总和必须为零[^1]。

于是我们可以将熵产生率$\dot{S}$写为

$$
\dot{S} = \left[\frac{\partial 1/T(x)}{\partial x}j_q - \sum_i \frac{\partial \mu_i(x)/T(x)}{\partial x}j_i\right]\mathrm{d}x \, A ,

\tag{D.1.3}
$$

或者，以$\sigma$（单位体积的熵产生）表示：

$$
\sigma = \frac{\partial 1/T(x)}{\partial x}j_q - \sum_i \frac{\partial \mu_i(x)/T(x)}{\partial x}j_i .

\tag{D.1.4}
$$

我们可以将此稍作推广，取消输运是一维的假设，并将化学势梯度与作用于物种$i$的其他外力$F_i$（由外势梯度引起）结合起来。于是得到

$$
\sigma = \mathbf{J}_q \cdot \nabla \frac{1}{T} - \frac{1}{T}\sum_i \mathbf{J}_i \cdot \nabla \frac{\mu_i}{T} - \mathbf{F}_i .

\tag{D.1.5}
$$

出于下文将要讨论的原因，将$\mu_i/T$的梯度分离为依赖于温度梯度的部分和不依赖于温度梯度的部分是方便的。如果化学势依赖于温度和压力，我们可以利用

$$
\left.\frac{\partial \beta \mu_i}{\partial \beta}\right|_{P,\{N_j\}} = h_i ,

\tag{D.1.6}
$$

其中$\beta \equiv 1/k_B T$，以及

$$
\sigma = \left(\mathbf{J}_q - \sum_i \mathbf{J}_i h_i\right) \cdot \nabla \frac{1}{T} - \frac{1}{T}\sum_i \mathbf{J}_i \cdot \left(\nabla \mu_i\right)_T - \mathbf{F}_i ,

\tag{D.1.7}
$$

这使我们能够定义“非扩散性”的不可逆热流$\mathbf{J}'_q$：

$$
\mathbf{J}'_q \equiv \mathbf{J}_q - \sum_i \mathbf{J}_i h_i .

\tag{D.1.8}
$$

在下文中，我们将用$\mathbf{J}_h$表示焓流（enthalpy flux）：

$$
\mathbf{J}_h \equiv \sum_i \mathbf{J}_i h_i .

\tag{D.1.9}
$$

### 焓流

为什么要减去焓流？首先：正如 de Groot 和 Mazur 所指出的，混合物中的热流不是唯一确定的。焓流也是如此。焓流的值取决于我们如何选择粒子能量的零点。例如，我们可以将与粒子$i$的静止质量相关的$m_i c^2$包括在内。这并不像看起来那么荒谬，因为例如，如果我们正在泵送$\mathrm{UF}_6$或类似的核燃料，能量流确实会将静止质量考虑在内。

关键在于，“内能”的选择对于“在恒定能量下转移一个粒子”的含义有着巨大的影响。如果我们将粒子$i$从储库 1 移动到储库 2，而不移动能量，那么从储库 1 中移除粒子$i$将导致熵的大幅增加（因为曾经存在于粒子$i$静止质量中的能量现在被其他粒子的库所吸收），反过来，系统 2 的熵将会大幅减少以补偿与引入粒子$i$相关的能量增加。可以方便地将这个过程看作：我们允许粒子$i$保留其焓，然后添加一个从 2 到 1 的焓流来精确补偿这一点。在恒定能量下转移粒子的过程如图 D.1 所示。

![图 D.1](../images/fig_D_1.png)

*图 D.1　在不改变任一体系能量的前提下，将粒子从温度为 $T_1$ 的体系 1 转移到温度为 $T_2$ 的体系 2，这要求每个粒子把全部能量（甚至动能）都留在体系 1 中，并以消耗体系 2 的能量为代价获得其新的能量。*

具体而言，我们可以考虑这样的情况：除了粒子转移之外，还存在一个从 1 到 2 的热流$\sum_i \mathbf{J}_i h_i$。这对应于我们不强制粒子在离开 1 时提取能量或在到达 2 时添加能量的情况。在这种特殊情况下，扩散焓输运对净“不可逆”热流$\mathbf{J}'_q$的贡献为零。重要的是，如果我们允许粒子携带其相关的焓进行传输，那么参考态的整个问题就消失了。当然，在现实中，仍然会存在热流，但这是与热运动和分子间相互作用相关的热。

总结如下：熵产生率由下式给出：

$$
\sigma = \mathbf{J}'_q \cdot \nabla \frac{1}{T} - \frac{1}{T}\sum_i \mathbf{J}_i \cdot \left[\left(\nabla \mu_i\right)_T - \mathbf{F}_i\right] ,

\tag{D.1.10}
$$

其中不可逆热流不依赖于粒子流$\mathbf{J}_i$所携带的焓。

## 涨落

在前面的章节中，我们论述了系统的平衡态对应于具有最多微观实现数目的状态。此外，我们论证了可以将这一概率图景与封闭系统的熵达到最大值的实验观测联系起来，通过将熵$S$与$k_B \ln \Omega$等同：

$$
S = k_B \ln \Omega .

\tag{D.2.1}
$$

现在我们假设熵是$n$个线性无关广延变量$\{A_1, A_2, \cdots, A_n\}$的唯一函数。在平衡态下，孤立系统的熵必须达到最大值。如果$A_i$表示一个非守恒量（例如，结晶度），则根据第二定律，在平衡态时，

$$
\left.\frac{\partial S_t}{\partial A_i}\right|_{\mathrm{eq}} = 0 ,

\tag{D.2.2}
$$

其中$S_t$表示整个系统的熵。然而，如果$A_i$表示一个守恒量，则它在封闭系统中不能改变。在这种情况下，我们可以考虑当$\mathrm{d}A_i$从子系统 1 转移到子系统 2 时熵的变化。于是

$$
\left.\frac{\partial S_t}{\partial A_i^{(2)}}\right|_{\mathrm{eq}} = \left.\frac{\partial S^{(2)}}{\partial A_i^{(2)}}\right|_{\mathrm{eq}} - \left.\frac{\partial S^{(1)}}{\partial A_i^{(1)}}\right|_{\mathrm{eq}} = 0 ,

\tag{D.2.3}
$$

因为$\mathrm{d}A_i^{(1)} = -\mathrm{d}A_i^{(2)}$。更一般地，如果我们考虑$m$个子系统，若$A_i$是守恒量，则有$m-1$个独立变量$A_i^{(n)}$；若$A_i$不是守恒量，则有$m$个独立变量。对于守恒量，式 (D.2.3) 简单地表达了两个子系统中$T$、$P/T$或$\mu_i/T$相等的条件。因此，$S$在$A_i$的线性阶不发生变化。然而，在$A_i$的二次阶，$S$确实会变化[^2]。因此，到二次阶，我们可以写出

$$
S = S_0 + \frac{1}{2}\sum_{i,j} \left.\frac{\partial^2 S_t}{\partial A_i \partial A_j}\right|_{\mathrm{eq}} \alpha_i \alpha_j ,

\tag{D.2.4}
$$

其中$\alpha_i \equiv A_i - A_i^0$。为了与 De Groot 和 Mazur [[57]](references.md#ref-57)的记号保持一致，我们写

$$
g_{ij} \equiv -\left.\frac{\partial^2 S_t}{\partial A_i \partial A_j}\right|_{\mathrm{eq}} .

\tag{D.2.5}
$$

如果$\alpha_i$是线性无关的，则$g_{ij}$是一个对称正定矩阵。系统处于偏离（但接近）最概然状态、由变量$\{\alpha_1, \alpha_2, \cdots, \alpha_k\}$表征的状态的概率为

$$
P(\{\alpha_1, \alpha_2, \cdots, \alpha_k\}) \propto \exp\!\left(-\frac{1}{2k_B}\sum_{i,j} g_{ij}\alpha_i \alpha_j\right) .

\tag{D.2.6}
$$

由于第二定律表明系统将从概率较低的状态演化到概率较高的状态，我们可以定义使系统返回其最概然状态的“驱动力”（driving force）。与每个变量$\alpha_i$相关联的驱动力$X_i$由下式给出：

$$
X_i = -\left.\frac{\partial}{\partial \alpha_i}\frac{1}{2}\sum_{j} g_{ij}\alpha_i \alpha_j \right|_{\mathrm{eq}} = -\sum_j g_{ij}\alpha_j .

\tag{D.2.7}
$$

为了定义驱动力$X_i$，我们利用了局部热力学平衡假设，即局部地，熵$S$与基本广延热力学量$U$、$V$和$M_i$（组分$i$的质量）之间的关系与准静态过程中相同。我们先验地不知道系统返回平衡的速度有多快。需要一组额外的本构方程（constitutive equations）来描述驱动力$X_j$与“流”$\mathbf{J}_i \equiv \dot{\alpha}_i$之间的关系。我们假设，在最低阶，这些关系具有如下形式：

$$
J_i = \sum_j L_{ij} X_j .

\tag{D.2.8}
$$

在这一阶段，我们对输运系数（transport coefficients）$L_{ij}$一无所知。Onsager [[55]](references.md#ref-55)假设，描述变量$\alpha_i$向其平衡值（零）衰减速度的定律对于任意小的$\alpha_i$都成立，因此也描述了平衡附近自发涨落向其平均值衰减的速度。这一“昂萨格回归假设”（Onsager Regression Hypothesis）提供了宏观输运系数与平衡态系统微观动力学之间的联系。回归假设可以被视为爱因斯坦关于扩散输运可以处理为布朗运动的宏观表现这一假设的推广。

由式 (D.2.4) 可以容易地导出熵产生$\dot{S}$的表达式：

$$
\dot{S} = -\sum_{i,j} g_{ij}\dot{\alpha}_i \alpha_j = \sum_i X_i \cdot \dot{\alpha}_i = \sum_{i,j} L_{ij} X_i X_j .

\tag{D.2.9}
$$

## 昂萨格倒易关系

为了与平衡态中涨落的衰减建立联系，我们现在证明平衡涨落$\alpha_i$仅与其共轭力$X_j$相关：

$$
\langle \alpha_i X_j \rangle = -k_B \delta_{ij} ,

\tag{D.3.1}
$$

其中$k_B$是玻尔兹曼常数，$\delta_{ij}$是克罗内克 $\delta$。式 (D.3.1) 可由以下事实推导：

$$
\langle \alpha_i X_j \rangle = k_B \int \mathrm{d}\{\alpha\}\, \alpha_i \frac{\partial P(\{\alpha\})}{\partial \alpha_j} = -k_B \int \mathrm{d}\{\alpha\}\, \frac{\partial \alpha_i}{\partial \alpha_j} P(\{\alpha\}) = -k_B \delta_{ij} ,
\tag{D.3.2}
$$

其中$\{\alpha\}$表示集合$\{\alpha_1, \alpha_2, \cdots, \alpha_k\}$。利用$J_i = \sum_j L_{ij} X_j$，可以得到：

$$
\langle \alpha_j(t) J_i(t) \rangle = -k_B L_{ij} .

\tag{D.3.3}
$$

式 (D.3.3) 使我们能够推导昂萨格倒易关系（Onsager reciprocal relations）。但在推导之前，请注意，从微观上看，式 (D.3.3) 有些奇怪，因为时刻$t$的涨落只会导致$t > 0$时的流。事实上，由于$\alpha_j$和$J_j$具有不同的时间反演对称性，等时积$\langle \alpha_j(t) J_i(t) \rangle$实际上为零。非零的量是$\langle \alpha_j(t) J_i(t+\epsilon) \rangle$。我们稍后将回到这一点。目前，我们继续使用式 (D.3.3) 并将其重写为：

$$
\int_0^\infty \mathrm{d}t\, \langle \alpha_j(0) \dot{J}_i(t) \rangle = -k_B L_{ij} ,

\tag{D.3.4}
$$

其中我们利用了$\langle \alpha_j(0) \dot{J}_i(0) \rangle$为零这一事实。接下来，利用时间平移不变性：

$$
\int_0^\infty \mathrm{d}t\, \langle \dot{\alpha}_j(0) J_i(t) \rangle = +k_B L_{ij} ,

\tag{D.3.5}
$$

并利用$J_j = \dot{\alpha}_j$这一事实：

$$
\int_0^\infty \mathrm{d}t\, \langle J_j(0) J_i(t) \rangle = k_B L_{ij} .

\tag{D.3.6}
$$

重要的是（经典地）：

$$
\langle J_j(0) J_i(t) \rangle = \langle J_i(t) J_j(0) \rangle = \langle J_i(0) J_j(-t) \rangle .
$$

在此及下文中，我们将讨论限制在所有流具有相同时间反演对称性的情况，此时$\langle J_j(0) J_i(t) \rangle = \langle J_i(0) J_j(t) \rangle$，由此可得

$$
L_{ij} = L_{ji} .

\tag{D.3.7}
$$

这就是昂萨格倒易关系：它表明输运系数矩阵$L_{ij}$是对称的。这一关系是非平衡热力学的核心结果之一，在计算如互扩散系数等耦合输运性质时具有深远的意义。

---

[^1]: 粒子流最优雅的选择是质量流，在这种情况下，所有粒子流之和为零的条件意味着质心是静止的。然而，在实际应用中，使用哪种流差别不大。我们将使用数密度。需要记住的是，化学势的定义（例如每粒子、每单位质量或每摩尔）必须与流的选取保持一致。
[^2]: 我们假设$S$是$A_i$的解析函数。这似乎是合理的，但并不总是成立。