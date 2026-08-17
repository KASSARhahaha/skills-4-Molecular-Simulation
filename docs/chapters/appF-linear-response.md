# 线性响应：示例

## 耗散

许多实验技术通过测量外加场（如可见光、红外辐射、微波辐射）的吸收来探测多体系统的动力学。线性响应理论使我们能够建立吸收谱与时间关联函数傅里叶变换之间的简单关系。为了理解这一点，让我们再次考虑一个与动力学变量$A(\mathbf{p}^N, \mathbf{q}^N)$耦合的外场。系统含时哈密顿量$H$为

$$
H(t) = H_0 - f(t) A(\mathbf{p}^N, \mathbf{q}^N).
$$

注意唯一显含时间的量是$f(t)$。由于哈密顿量依赖于时间，系统的总能量$E$也随时间变化：

$$
E(t) = \langle H(t) \rangle.
$$

让我们计算系统能量变化的平均速率。这是系统吸收（或发射）的能量，单位时间内为：

$$
\frac{\partial E}{\partial t} = \left\langle \frac{\mathrm{d}H}{\mathrm{d}t} \right\rangle
= \left\langle \sum_i \left( \dot{q}_i \frac{\partial H}{\partial q_i} + \dot{p}_i \frac{\partial H}{\partial p_i} \right) + \frac{\partial H}{\partial t} \right\rangle.
$$

但由哈密顿方程，我们有

$$
\dot{q}_i = \frac{\partial \mathcal{H}}{\partial p_i} \quad \text{和} \quad \dot{p}_i = -\frac{\partial \mathcal{H}}{\partial q_i}.
\tag{F.1.1}
$$

因此，式 (F.1.1) 简化为

$$
\frac{\partial E}{\partial t} = \left\langle \frac{\partial H}{\partial t} \right\rangle = -\left\langle \dot{f}(t) A(\mathbf{p}^N, \mathbf{q}^N) \right\rangle = -\dot{f}(t) \langle A(t) \rangle.
\tag{F.1.2}
$$

但请注意，$\langle A(t) \rangle$本身就是对所加场$f$的响应：

$$
\langle A(t) \rangle = \int_{-\infty}^{\infty} \mathrm{d}t' \, \chi_{AA}(t - t') f(t').
$$

现在考虑$f(t)$是频率为$\omega$的周期场（如单色光）的情况。在这种情况下，我们可以将$f(t)$写为

$$
f(t) = \mathrm{Re} \left[ f_{\omega} e^{i\omega t} \right]
$$

和

$$
\dot{f}(t) = \frac{i\omega}{2} \left[ f_{\omega} e^{i\omega t} - f_{\omega}^* e^{-i\omega t} \right].
$$

平均能量耗散速率为

$$
\frac{\partial E}{\partial t} = -\dot{f}(t) \langle A(t) \rangle = -\dot{f}(t) \int_{-\infty}^{\infty} \mathrm{d}t' \, \chi_{AA}(t - t') f(t').
\tag{F.1.3}
$$

对于周期场，我们有

$$
\begin{aligned}
\int_{-\infty}^{\infty} \mathrm{d}t' \, \chi_{AA}(t - t') f(t') &= \frac{f_{\omega} e^{i\omega t}}{2} \int_{-\infty}^{\infty} \mathrm{d}t' \, \chi_{AA}(t - t') e^{i\omega(t' - t)}\\
&\quad + \frac{f_{\omega}^* e^{-i\omega t}}{2} \int_{-\infty}^{\infty} \mathrm{d}t' \, \chi_{AA}(t - t') e^{-i\omega(t' - t)},
\end{aligned}
$$

其中

$$
\chi_{AA}(\omega) \equiv \frac{1}{2\pi} \int_{0}^{\infty} \mathrm{d}t \, \chi_{AA}(t) e^{-i\omega t}.
$$

为了计算能量变化速率$\dot{E}$，我们必须将$\langle \partial H / \partial t \rangle$在一个周期$T$（$= 2\pi/\omega$）内取平均：

$$
\begin{aligned}
\dot{E} &= -\frac{\pi}{2T} \int_0^T \mathrm{d}t \left[ i\omega \left( f_{\omega} e^{i\omega t} - f_{\omega}^* e^{-i\omega t} \right) \times \left( f_{\omega} e^{i\omega t} \chi_{AA}(\omega) + f_{\omega}^* e^{-i\omega t} \chi_{AA}(-\omega) \right) \right]\\
&= -\pi \omega |f_{\omega}|^2 \frac{\chi_{AA}(\omega) - \chi_{AA}(-\omega)}{2i}\\
&= -\pi \omega |f_{\omega}|^2 \, \mathrm{Im}[\chi_{AA}(\omega)].
\end{aligned}
\tag{F.1.4}
$$

我们利用$\chi_{AA}(t)$与$A$的自关联函数（式 2.5.17）之间的关系：

$$
\chi_{AA}(\omega) = \frac{1}{2\pi} \int_0^{\infty} \mathrm{d}t \, e^{-i\omega t} \left[ -\beta \langle A(0) \dot{A}(t) \rangle \right].
\tag{F.1.5}
$$

$\chi_{AA}(\omega)$的虚部为

$$
\mathrm{Im}[\chi_{AA}(\omega)] = \frac{\beta}{2\pi} \int_0^{\infty} \mathrm{d}t \, \sin(\omega t) \langle A(0) \dot{A}(t) \rangle = -\frac{\beta}{4\pi} \int_{-\infty}^{\infty} \mathrm{d}t \, \omega \cos(\omega t) \langle A(0) A(t) \rangle.
\tag{F.1.6}
$$

最终，我们得到

$$
\dot{E} = |f_{\omega}|^2 \frac{\beta \omega^2}{4} \int_{-\infty}^{\infty} \mathrm{d}t \, \cos(\omega t) \langle A(0) A(t) \rangle.
\tag{F.1.7}
$$

因此，通过了解与外场耦合的量的自关联函数，我们可以计算吸收谱的形状。这个关系是在假设经典动力学的情况下推导的，因此仅在$\hbar \omega \ll k_B T$时才有效。然而，也有可能推导出适用于任意频率的量子力学版本的线性响应理论（参见例如文献[[53]](references.md#ref-53)）。

举一个具体的例子，让我们计算稀薄极性分子气体的吸收谱形状。在这种情况下，相关的关联函数是偶极（$\boldsymbol{\mu}$）自关联函数：

$$
\langle M_x(0) M_x(t) \rangle = \frac{N}{3} \langle \boldsymbol{\mu}(0) \cdot \boldsymbol{\mu}(t) \rangle.
$$

对于几乎自由旋转的分子（“几乎”，否则就没有耗散），$\boldsymbol{\mu}(0) \cdot \boldsymbol{\mu}(t)$依赖于时间，因为每个分子都在旋转。对于一个旋转频率为$\omega$的分子，我们有

$$
\boldsymbol{\mu}(0) \cdot \boldsymbol{\mu}(t) = \mu^2 \cos(\omega t),
$$

对于具有旋转速度热分布$P(\omega)$的分子集合，我们有

$$
\langle \boldsymbol{\mu}(0) \cdot \boldsymbol{\mu}(t) \rangle = \mu^2 \int \mathrm{d}\omega \, P(\omega) \cos(\omega t).
$$

辐射吸收速率则由下式给出：

$$
\dot{E} = \frac{\pi \beta \omega^2 N \mu^2}{12} P(\omega) |f_{\omega}|^2.
\tag{F.1.8}
$$

关于光谱性质与时间关联函数之间关系的更多细节，读者可参阅 Madden 在文献[[44]](references.md#ref-44)中的文章。

## 电导率

在第 2 章线性响应理论的推导中，我们假设系统在扰动开启时被制备为平衡态，然后让系统弛豫到扰动关闭的新平衡态。然而，这并不总是可行的。例如，考虑电导率。在这种情况下，扰动是一个电场，它将引起系统中电流的流动。因此，我们用电场开启时制备的系统的状态不是平衡态，而是稳定的非平衡态。恒定剪切下的系统也是如此。在这些情况下，似乎不能用最简单形式的线性响应理论框架来推导输运系数，如电导率$\sigma_e$或黏度$\eta$。幸运的是，情况并没有那么糟糕。以电导率为例，确实，如果我们将导电系统置于外场中，我们将产生一个非平衡稳态。但是，我们可以通过开启一个弱的均匀矢势$\mathbf{A}$来扰动系统。开启矢势后系统的哈密顿量为

$$
\mathcal{H}' = \sum_{i=1}^{N} \frac{1}{2m_i} \left( \mathbf{p}_i - \frac{e_i}{c} \mathbf{A} \right)^2 + \mathcal{U}_{\mathrm{pot}}.
\tag{F.2.1}
$$

由该哈密顿量描述的系统满足与未扰动系统相同的运动方程（$\mathbf{A}$是规范场），系统在$t = 0$时处于平衡态。然后我们突然关闭矢势。从电动力学可知，时变矢势产生电场：

$$
\mathbf{E} = -\frac{1}{c} \dot{\mathbf{A}}.
\tag{F.2.2}
$$

在目前的情况下，电场将是$t = 0$处的一个无穷小$\delta$尖峰：

$$
\mathbf{E}(t) = \frac{1}{c} \mathbf{A} \, \delta(t).
\tag{F.2.3}
$$

我们可以用标准方法计算由此产生的电流。注意到我们可以将式 (F.2.1) 中的$H'$写为

$$
\begin{aligned}
\mathcal{H}' &= \mathcal{H}_0 - \sum_{i=1}^{N} \frac{e_i}{c m_i} \mathbf{p}_i \cdot \mathbf{A} + \mathcal{O}(A^2)\\
&= \mathcal{H}_0 - \frac{\mathbf{A}}{c} \int \mathrm{d}\mathbf{r} \sum_{i=1}^{N} \frac{e_i}{m_i} \mathbf{p}_i \, \delta(\mathbf{r}_i - \mathbf{r})\\
&= \mathcal{H}_0 - \frac{\mathbf{A}}{c} \int \mathrm{d}\mathbf{r} \, \mathbf{j}(\mathbf{r}),
\end{aligned}
\tag{F.2.4}
$$

其中$\mathbf{j}(\mathbf{r})$表示点$\mathbf{r}$处的电流密度。扰动在时间$t$引起的平均电流密度为

$$
\langle \mathbf{j}(t) \rangle = \frac{\mathbf{A}}{c V k_B T} \int \mathrm{d}\mathbf{r} \, \mathrm{d}\mathbf{r}' \, \langle \mathbf{j}(\mathbf{r}, 0) \cdot \mathbf{j}(\mathbf{r}', t) \rangle.
\tag{F.2.5}
$$

对所加$\delta$函数电场尖峰的电流响应的唯象表达式为（参见式 2.5.14）

$$
\langle \mathbf{j}(t) \rangle = \int_{-\infty}^{t} \mathrm{d}t' \, \boldsymbol{\sigma}(t - t') \mathbf{E}(t') = \boldsymbol{\sigma}(t) \frac{\mathbf{A}}{c}.
\tag{F.2.6}
$$

由此立即可得

$$
\boldsymbol{\sigma}(t) = \frac{1}{V k_B T} \int \mathrm{d}\mathbf{r} \, \mathrm{d}\mathbf{r}' \, \langle \mathbf{j}(\mathbf{r}, 0) \cdot \mathbf{j}(\mathbf{r}', t) \rangle.
\tag{F.2.7}
$$

直流电导率则为

$$
\boldsymbol{\sigma}(\omega = 0) = \frac{1}{V k_B T} \int_0^{\infty} \mathrm{d}t \int \mathrm{d}\mathbf{r} \, \mathrm{d}\mathbf{r}' \, \langle \mathbf{j}(\mathbf{r}, 0) \cdot \mathbf{j}(\mathbf{r}', t) \rangle.
\tag{F.2.8}
$$

## 黏度

黏度相应的线性响应表达式似乎更加微妙，因为剪切通常不被解释为作用于所有分子的外场。尽管如此，我们可以通过与电导率情况类比，使用正则变换——即对应于均匀剪切的时间导数。为此，我们考虑一个由$N$个粒子组成的系统，其坐标为$\mathbf{r}^N$，哈密顿量为

$$
\mathcal{H}_0 = \sum_{i=1}^{N} \frac{p_i^2}{2m_i} + \mathcal{U}(\mathbf{r}^N).
\tag{F.3.1}
$$

现在考虑另一个系统，由一组坐标$\mathbf{r}'^N$描述，与$\mathbf{r}^N$通过线性变换相关：

$$
\mathbf{r}'_i = \mathbf{h} \mathbf{r}_i.
\tag{F.3.2}
$$

新系统的哈密顿量可以写为

$$
\mathcal{H}_1 = \sum_{i=1}^{N} \frac{1}{2m_i} \mathbf{p}'_i \cdot \mathbf{G}^{-1} \cdot \mathbf{p}'_i + \mathcal{U}(\mathbf{r}'^N),
\tag{F.3.3}
$$

其中$\mathbf{G}$是度量张量，定义为

$$
\mathbf{G} \equiv \mathbf{h}^T \cdot \mathbf{h}.
\tag{F.3.4}
$$

我们假设$\mathbf{h}$与单位矩阵$\mathbf{I}$相差无穷小：

$$
\mathbf{h} = \mathbf{I} + \boldsymbol{\epsilon}.
\tag{F.3.5}
$$

在我们感兴趣的均匀剪切效应的情况下，例如可以选择$\epsilon_{xy} = \epsilon$，而$\boldsymbol{\epsilon}$的所有其他元素为 0。现在考虑我们用哈密顿量$H_1$使系统达到平衡，在$t = 0$时关闭无穷小变形$\boldsymbol{\epsilon}$。这意味着在$t = 0$时，系统经历剪切速率的$\delta$函数尖峰：

$$
\frac{\partial v_x}{\partial y} = -\epsilon \, \delta(t).
\tag{F.3.6}
$$

我们可以计算剪应力$\sigma_{xy}(t)$对从$H_1$到$H_0$突然变化的时间依赖响应：

$$
\langle \sigma_{xy}(t) \rangle = -\epsilon \frac{1}{V k_B T} \langle \sigma_{xy}(0) \sigma_{xy}(t) \rangle.
\tag{F.3.7}
$$

将式 (F.3.6) 和 (F.3.7) 与式 (2.5.14) 结合，我们立即看到由稳定剪切产生的稳态应力$\sigma_{xy}$为

$$
\sigma_{xy} = \frac{\partial v_x}{\partial y} \times \frac{1}{V k_B T} \int_0^{\infty} \mathrm{d}t \, \langle \sigma_{xy}(0) \sigma_{xy}(t) \rangle,
\tag{F.3.8}
$$

剪切黏度$\eta$的表达式为

$$
\eta = \frac{1}{V k_B T} \int_0^{\infty} \mathrm{d}t \, \langle \sigma_{xy}(0) \sigma_{xy}(t) \rangle.
\tag{F.3.9}
$$

## 弹性常数

液体在剪切力作用下会流动。固体则不会。相反，固体的任何小变形都会引起弹性响应（应力）来抵消它。这种弹性应力与所加变形（应变）成正比。应力与应变（将在下面更精确定义）之间的比例常数称为弹性常数。下面我们讨论如何通过计算机模拟来测量这些常数。为简单起见，我们将讨论限制在各向同性（静水）压力下的晶体。

当考虑应变对固体自由能的影响时，必须引入所谓的拉格朗日应变张量（参见例如文献[[182]](references.md#ref-182)）。[^1] 原因是在局部尺度上，自由能的所有变化都是由构成固体的粒子之间距离的变化引起的。而测量这种变化的量正是拉格朗日应变。我们从弹性变形引起的新旧坐标关系开始：

$$
\mathbf{r}' = (\mathbf{1} + \boldsymbol{\epsilon}) \mathbf{r},
\tag{F.4.1}
$$

其中

$$
\epsilon_{\alpha\beta} \equiv \frac{\partial u_{\alpha}}{\partial x_{\beta}}
\tag{F.4.2}
$$

是（常规）应变张量。它度量位移场$\mathbf{u}$随原始坐标$\mathbf{r}$的变化。由于应变，固体中两点$i$和$j$之间的距离$r_{ij}$发生了变化。新的平方距离与旧距离的关系为

$$
r_{ij}'^2 = \mathbf{r}_{ij} (\mathbf{1} + \boldsymbol{\epsilon}^T) (\mathbf{1} + \boldsymbol{\epsilon}) \mathbf{r}_{ij} = \mathbf{r}_{ij} (\mathbf{1} + \boldsymbol{\epsilon}^T + \boldsymbol{\epsilon} + \boldsymbol{\epsilon}^T \boldsymbol{\epsilon}) \mathbf{r}_{ij} \equiv \mathbf{r}_{ij} (\mathbf{1} + 2\boldsymbol{\eta}) \mathbf{r}_{ij}.
$$

这定义了拉格朗日应变张量$\boldsymbol{\eta}$。系统的新体积$V'$与原始体积$V_0$的关系为

$$
V' = V_0 \det(\mathbf{1} + \boldsymbol{\epsilon})
\tag{F.4.3}
$$

或

$$
V' = V_0 \sqrt{\det(\mathbf{1} + 2\boldsymbol{\eta})}.
\tag{F.4.4}
$$

现在我们将单位（未变形）体积$(V)$的亥姆霍兹自由能$(F)$按拉格朗日应变参数$\boldsymbol{\eta}$的幂次展开：

$$
\begin{aligned}
F(\boldsymbol{\eta})/V &= V^{-1} \left[ F(0) + \frac{\partial F}{\partial \eta_{\alpha\beta}} \eta_{\alpha\beta} + \frac{1}{2} \frac{\partial^2 F}{\partial \eta_{\alpha\beta} \partial \eta_{\gamma\delta}} \eta_{\alpha\beta} \eta_{\gamma\delta} + \cdots \right]\\
&= V^{-1} F(0) + C^{(1)}_{\alpha\beta} \eta_{\alpha\beta} + \frac{1}{2} C^{(2)}_{\alpha\beta\gamma\delta} \eta_{\alpha\beta} \eta_{\gamma\delta} + \cdots.
\end{aligned}
\tag{F.4.5}
$$

该方程定义了（二阶）弹性常数$C^{(2)}_{\alpha\beta\gamma\delta}$。为了数值计算弹性常数，我们需要$F$对$\boldsymbol{\eta}$依赖关系的微观表达式。为了推导这种关系，我们必须详细考虑系统变形对配分函数的影响。让我们首先考虑变形系统。该系统的配分函数（忽略常数，如$h^{-3N}$）为

$$
Q(\boldsymbol{\eta}) = \int \mathrm{d}\mathbf{p}^N \mathrm{d}\mathbf{r}^N \exp\left[ -\beta \mathcal{H}\left( \mathbf{p}^N, \mathbf{r}^N \right) \right].
\tag{F.4.6}
$$

该配分函数通过坐标积分的边界条件依赖于变形。这在计算对应变的导数时不太方便。因此，我们首先用原始未变形系统的坐标和动量来表示变形系统的配分函数。我们可以用应变张量$\mathbf{h} \equiv (\mathbf{1} + \boldsymbol{\epsilon})$和原始坐标$(\mathbf{r}_{0,i})$及速度$(\dot{\mathbf{r}}_{0,i})$来表示该系统中的坐标$(\mathbf{r}_i)$和速度$(\dot{\mathbf{r}}_i)$：

$$
\begin{align}
\mathbf{r}_i &= \mathbf{h} \mathbf{r}_{0,i} \nonumber \\
\dot{\mathbf{r}}_i &= \mathbf{h} \dot{\mathbf{r}}_{0,i}
\tag{F.4.7}
\end{align}
$$

动能$K = \sum \frac{1}{2} m_i \dot{\mathbf{r}}_i^2$可以写为

$$
\mathcal{K} = \sum \frac{1}{2} m_i \dot{\mathbf{r}}_i^2 = \sum \frac{1}{2} m_i \dot{\mathbf{r}}_{0,i} (\mathbf{h}^T \mathbf{h}) \dot{\mathbf{r}}_{0,i} \equiv \sum \frac{1}{2} m_i \dot{\mathbf{r}}_{0,i} \cdot \mathbf{G} \cdot \dot{\mathbf{r}}_{0,i},
\tag{F.4.8}
$$

其中$\mathbf{h}^T = (\mathbf{1} + \boldsymbol{\epsilon}^T)$是$\mathbf{h}$的转置，$\mathbf{G} = \mathbf{h}^T \mathbf{h}$是度量张量。由$\mathbf{h}$的定义可知$\mathbf{G} = (\mathbf{1} + 2\boldsymbol{\eta})$。我们现在可以写出与坐标$\mathbf{r}_{0,i}$共轭的广义动量$\mathbf{p}_{0,i}$（见附录 A）：

$$
p_{0,i}^{\alpha} = \left( \frac{\partial \mathcal{K}}{\partial \dot{r}_{0,i}^{\alpha}} \right) = m_i G^{\alpha\beta} \dot{r}_{0,i}^{\beta}
\tag{F.4.9}
$$

因此

$$
\mathcal{K} = \sum \frac{1}{2} m_i \dot{\mathbf{r}}_{0,i} \cdot \mathbf{G} \cdot \dot{\mathbf{r}}_{0,i} = \sum \frac{1}{2m_i} \mathbf{p}_{0,i} \cdot \mathbf{G}^{-1} \cdot \mathbf{p}_{0,i} = \sum \frac{1}{2m_i} \mathbf{p}_{0,i} \cdot (\mathbf{1} + 2\boldsymbol{\eta})^{-1} \cdot \mathbf{p}_{0,i}.
\tag{F.4.10}
$$

由于

$$
\begin{align}
\mathbf{p}_i &= m_i \dot{\mathbf{r}}_i = m_i \mathbf{h} \dot{\mathbf{r}}_{0,i} = (\mathbf{h}^T)^{-1} \mathbf{p}_{0,i} \nonumber \\
\mathbf{r}_i &= \mathbf{h} \mathbf{r}_{0,i},
\tag{F.4.11}
\end{align}
$$

$\{\mathbf{p}^N, \mathbf{r}^N\}$与$\{\mathbf{p}_0^N, \mathbf{r}_0^N\}$之间变换的雅可比行列式等于 1。因此，我们可以写

$$
\begin{aligned}
Q(\boldsymbol{\eta}) &= \int \mathrm{d}\mathbf{p}^N \mathrm{d}\mathbf{r}^N \exp\left[ -\beta \mathcal{H}\left( \mathbf{p}^N, \mathbf{r}^N \right) \right]\\
&= \int \mathrm{d}\mathbf{p}_0^N \mathrm{d}\mathbf{r}_0^N \exp\left[ -\beta \left\{ \sum \frac{1}{2m_i} \mathbf{p}_{0,i} \cdot (\mathbf{1} + 2\boldsymbol{\eta})^{-1} \cdot \mathbf{p}_{0,i} + \mathcal{U}\left( \mathbf{r}_0^N; \boldsymbol{\eta} \right) \right\} \right].
\end{aligned}
\tag{F.4.12}
$$

现在$Q(\boldsymbol{\eta})$对$\boldsymbol{\eta}$的依赖仅包含在哈密顿量中。我们现在可以显式地完成对$\boldsymbol{\eta}$的微分。利用

$$
\begin{aligned}
\left( \frac{\partial U}{\partial \eta_{\alpha\beta}} \right) &= \sum_{i<j} \left( \frac{\partial U}{\partial r_{ij}^2} \right) \left( \frac{\partial r_{ij}^2}{\partial \eta_{\alpha\beta}} \right) = \sum_{i<j} \left( \frac{\partial U}{\partial r_{ij}} \right) \frac{r_{0,ij}^{\alpha} r_{0,ij}^{\beta}}{r_{ij}}\\
&= \left( \mathbf{h}^{-1} \sum_{i<j} \frac{\partial U}{\partial r_{ij}} \frac{\mathbf{r}_{ij} \mathbf{r}_{ij}}{r_{ij}} (\mathbf{h}^T)^{-1} \right)_{\alpha\beta}
\end{aligned}
\tag{F.4.13}
$$

和

$$
\begin{aligned}
\sum \frac{1}{2m_i} \mathbf{p}_{0,i} \cdot \left( \frac{\partial \mathbf{G}^{-1}}{\partial \eta_{\alpha\beta}} \right) \cdot \mathbf{p}_{0,i} &= -\sum \frac{1}{m_i} \left( \mathbf{p}_{0,i} \cdot \mathbf{G}^{-1} \right)_{\alpha} \left( \mathbf{G}^{-1} \cdot \mathbf{p}_{0,i} \right)_{\beta}\\
&= -\left( \mathbf{h}^{-1} \sum \frac{1}{m_i} \mathbf{p} \mathbf{p} (\mathbf{h}^T)^{-1} \right)_{\alpha\beta},
\end{aligned}
\tag{F.4.14}
$$

我们得到

$$
\left( \frac{\partial F}{\partial \eta_{\alpha\beta}} \right) = -\left( \mathbf{h}^{-1} \left[ \sum \frac{1}{m_i} \mathbf{p} \mathbf{p} + \sum_{j<i} \mathbf{r}_{ij} \mathbf{f}_{ij} \right] (\mathbf{h}^T)^{-1} \right)_{\alpha\beta}.
\tag{F.4.15}
$$

由此立即可得

$$
\begin{aligned}
C^{(1)}_{\alpha\beta} &\equiv \left( \frac{\partial F}{\partial \eta_{\alpha\beta}} \right)\\
&= \frac{V'}{V} \left[ (\mathbf{1} + \boldsymbol{\epsilon})^{-1} \boldsymbol{\sigma} (\mathbf{1} + \boldsymbol{\epsilon}^T)^{-1} \right]_{\alpha\beta}\\
&= \sqrt{\det(\mathbf{1} + 2\boldsymbol{\eta})} \left[ (\mathbf{1} + \boldsymbol{\epsilon})^{-1} \boldsymbol{\sigma} (\mathbf{1} + \boldsymbol{\epsilon}^T)^{-1} \right]_{\alpha\beta},
\end{aligned}
\tag{F.4.16}
$$

其中

$$
\sigma_{\gamma\delta} = -\frac{1}{V'} \sum_i \left[ \frac{p_i^{\gamma} p_i^{\delta}}{m_i} + \sum_{j<i} r_{ij}^{\gamma} f_{ij}^{\delta} \right]
\tag{F.4.17}
$$

表示变形系统中的微观应力。注意$\boldsymbol{\sigma}$可以在模拟中测量，而$\boldsymbol{\epsilon}$由所加应变固定。对于未变形的系统，$C^{(1)}$简单等于$-P$，其中$P$是静水压力。由式 (F.4.16) 还可得应力$\sigma_{ij}$与线性应变$\epsilon_{rs}$之间的比例常数为

$$
\displaystyle
\left( \frac{\partial \sigma_{\alpha\beta}{\partial \epsilon_{\gamma\delta}} \right) = (\sigma_{\alpha\delta} \delta_{\beta\gamma} + \sigma_{\beta\delta} \delta_{\alpha\gamma} - \sigma_{\alpha\beta} \delta_{\gamma\delta}) + C^{(2)}_{\alpha\beta\gamma\delta}.
}
\tag{F.4.18}
$$

为了确定二阶弹性常数$C^{(2)}_{\alpha\beta\gamma\delta}$，我们必须确定$C^{(1)}_{\alpha\beta}$对$\eta_{\gamma\delta}$的初始线性依赖。这种测量弹性常数的技术简单且相当精确（参见例如文献[[720]](references.md#ref-720)）。然而，需要多次计算来测量不同的弹性常数。晶体对称性越低，需要的计算次数越多。这可以通过直接考虑$C^{(2)}_{\alpha\beta\gamma\delta}$的微观表达式来避免。这样的表达式由 Squire 等人[[721]](references.md#ref-721)推导：

$$
\begin{aligned}
C^{(2)}_{\alpha\beta\gamma\delta} = {}& -\frac{1}{V k_B T} \left\langle \Delta\sigma_{\alpha\beta} \, \Delta\sigma_{\gamma\delta} \right\rangle + 2\rho k_B T (\delta_{\alpha\gamma} \delta_{\beta\delta} + \delta_{\alpha\delta} \delta_{\beta\gamma})\\
& + 4 \sum_{i<j,\,k<l} \left\langle \frac{\partial^2 U}{\partial r_{ij}^2 \partial r_{kl}^2} r_{ij\alpha} r_{ij\beta} r_{kl\gamma} r_{kl\delta} \right\rangle.
\end{aligned}
\tag{F.4.19}
$$

使用式 (F.4.19)，只需要一次模拟即可测量所有弹性常数。不幸的是，评估这个涨落表达式的统计误差通常大于计算式 (F.4.16) 时的误差。在恒应力 MD 模拟中统计问题更加严重，其中弹性柔度（而非模量）由盒子形状的涨落确定[[722]](references.md#ref-722)。式 (F.4.19) 仅在分子间势处处连续时才能使用。然而，Farago 和 Kantor [[723]](references.md#ref-723)已经开发了一种适用于硬核系统的涨落表达式。关于弹性常数数值评估的更多细节可以在文献[[436,721,723,724]](references.md#ref-436)中找到，而有限温度下弹性性质的一般框架在文献[[182]](references.md#ref-182)中讨论。

---

[^1]: 面对希腊字母的有限性，我们使用符号$\boldsymbol{\eta}$表示拉格朗日应变张量。这个符号很容易与标量黏度$\eta$混淆。