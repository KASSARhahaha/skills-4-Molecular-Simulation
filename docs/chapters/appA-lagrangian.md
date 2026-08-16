# 拉格朗日与哈密顿运动方程

了解牛顿运动方程足以理解分子动力学方法的基础。然而，许多更高级的模拟技术利用了经典力学的拉格朗日表述或哈密顿表述。在这里，我们简要勾勒这些不同方法之间的关系（另见[[54]](references.md#ref-54)）。关于经典力学更详细和更严格的描述，读者可参考 Goldstein 的著作[[54]](references.md#ref-54)。

## 作用量

经典力学的拉格朗日表述基于变分原理。经典系统在时间间隔$\{t_b, t_e\}$内、从初始位置$\mathbf{x}_b$到最终位置$\mathbf{x}_e$之间所遵循的实际轨迹，是使作用量$S$取极值（通常为极小值）的轨迹。经典作用量$S$对于任意轨迹的定义为系统动能$K$与势能$U$之差沿该轨迹的时间积分：

$$
S = \int_{t_b}^{t_e} \mathrm{d}t\,[K - U].
$$

在考虑由这个极值原理导出的一般拉格朗日运动方程之前，让我们先考虑几个简单的例子。

第一种情况是在没有外势的情况下运动的单个粒子，即$U = 0$。由于粒子必须在时间间隔$t_e - t_b$内从$\mathbf{x}_b$运动到$\mathbf{x}_e$，我们已经知道其平均速度：$\mathbf{v}_{\mathrm{av}}$。如果粒子始终以这个平均速度运动，它将遵循一条直线轨迹，我们用$\bar{\mathbf{x}}(t)$表示。设粒子的真实轨迹为$\mathbf{x}(t) = \bar{\mathbf{x}}(t) + \boldsymbol{\eta}(t)$，其中$\boldsymbol{\eta}(t)$是尚待确定的偏差。则粒子的速度是平均速度$\mathbf{v}_{\mathrm{av}}$和偏差$\dot{\boldsymbol{\eta}}(t)$之和：

$$
\mathbf{v}(t) = \mathbf{v}_{\mathrm{av}} + \dot{\boldsymbol{\eta}}(t).
$$

根据构造，

$$
\int \mathrm{d}t\,\dot{\boldsymbol{\eta}}(t) = 0.
$$

在本例中，势能始终为零，因此作用量$S$由动能的时间积分决定：

$$
S = \frac{1}{2}m \int \mathrm{d}t\,[\mathbf{v}_{\mathrm{av}} + \dot{\boldsymbol{\eta}}(t)]^2 = S_{\mathrm{av}} + \frac{1}{2}m \int \mathrm{d}t\,\dot{\boldsymbol{\eta}}^2(t).
$$

由于最后一项不可能小于零，作用量在$\dot{\boldsymbol{\eta}}(t) = 0$（对所有$t$）时取极小值。换言之，我们恢复了众所周知的结果：在没有外力的情况下，粒子以恒定速度运动。这就是牛顿第一定律。

接下来，考虑在一维势$U(x)$中运动的粒子。此时作用量为

$$
S = \int_{t_b}^{t_e} \mathrm{d}t \left[ \frac{1}{2}m \left(\frac{\mathrm{d}x(t)}{\mathrm{d}t}\right)^2 - \mathcal{U}(x) \right].
\tag{A.1.1}
$$

任意路径$x(t)$可以写成经典粒子将遵循的实际路径$\bar{x}(t)$加上偏离该路径的$\eta(t)$之和：

$$
x(t) = \bar{x}(t) + \eta(t).
$$

如前所述，我们施加粒子的初始和最终位置，因此$\eta(t_b) = \eta(t_e) = 0$。对于接近实际路径的路径$x(t)$，我们可以将作用量展开为（小）量$\eta(t)$的幂级数。实际上，由于$\eta(t)$本身就是$t$的函数，这样的展开称为泛函展开。如果此泛函展开中的前导（线性）项消失，则作用量取极值。现在让我们考虑作用量在真实路径的作用量附近的线性$\eta(t)$阶泛函展开：

$$
\begin{aligned}
S &= \int_{t_b}^{t_e} \mathrm{d}t\,\frac{1}{2}m \left(\frac{\mathrm{d}\bar{x}(t)}{\mathrm{d}t} + \frac{\mathrm{d}\eta(t)}{\mathrm{d}t}\right)^2 - U[\bar{x}(t) + \eta(t)] \\
&\approx \int_{t_b}^{t_e} \mathrm{d}t\,\frac{1}{2}m \left[\left(\frac{\mathrm{d}\bar{x}(t)}{\mathrm{d}t}\right)^2 + 2\frac{\mathrm{d}\bar{x}(t)}{\mathrm{d}t}\frac{\mathrm{d}\eta(t)}{\mathrm{d}t}\right] - \left[U(\bar{x}(t)) + \frac{\partial U(\bar{x})}{\partial x}\eta(t)\right] \\
&= \bar{S} + \int_{t_b}^{t_e} \mathrm{d}t \left[\frac{md\bar{x}(t)}{\mathrm{d}t}\frac{\mathrm{d}\eta(t)}{\mathrm{d}t} - \frac{\partial U(\bar{x})}{\partial x}\eta(t)\right] \\
&= \bar{S} + \left.m\frac{\mathrm{d}\bar{x}(t)}{\mathrm{d}t}\eta(t)\right|_{t_b}^{t_e} - \int_{t_b}^{t_e} \mathrm{d}t \left[\frac{md^2\bar{x}(t)}{\mathrm{d}t^2} + \frac{\partial U(\bar{x})}{\partial x}\right]\eta(t),
\end{aligned}
$$

其中最后一步通过分部积分得到。由于根据定义，$\eta(t)$在边界处为零，等式右边的第二项消失。如果上述方程最后一行中的被积函数对任意$\eta(t)$都为零，则作用量取极值。这个条件当且仅当以下关系成立时才能满足：

$$
m\frac{d^2\bar{x}(t)}{\mathrm{d}t^2} = -\frac{\partial U(\bar{x})}{\partial x},
\tag{A.1.2}
$$

这正是牛顿第二定律。换言之，牛顿运动方程可以从粒子遵循使作用量取极值的路径这一陈述中推导出来。

## 拉格朗日量

如果这种经典力学定律的替代表达方式不能使我们做比简单重新推导$F = ma$更多的事情，那么引入它就没有什么意义。事实上，经典力学的拉格朗日表述被证明是非常强大的。例如，拉格朗日方法使得在非笛卡尔坐标系中推导运动方程变得容易。假设我们希望使用某些广义坐标$q$代替笛卡尔坐标$x$。例如，考虑均匀重力场中长度为$l$的摆。摆与竖直方向（即与重力场方向）所成的角度可以用来指定其取向。由于摆所遵循的路径显然与我们碰巧用来指定其状态的坐标无关，作用量$S$应该是相同的：

$$
S = \int \mathrm{d}t\,\mathcal{L}(\mathbf{x},\dot{\mathbf{x}}) = \int \mathrm{d}t\,\mathcal{L}(q,\dot{q}),
\tag{A.2.1}
$$

其中量$L$称为拉格朗日量。拉格朗日量定义为动能减去势能[^1]：

$$
\mathcal{L} \equiv \mathcal{K}(\dot{q}) - \mathcal{U}(q).
\tag{A.2.2}
$$

我们再次引入实际路径$\bar{q}(t)$和偏离它的$\eta(t)$：

$$
\begin{aligned}
q(t) &= \bar{q}(t) + \eta(t) \\
\dot{q}(t) &= \dot{\bar{q}}(t) + \dot{\eta}(t).
\end{aligned}
$$

我们可以将拉格朗日量$L$写为：

$$
L(q,\dot{q}) = L(\bar{q},\dot{\bar{q}}) + \frac{\partial L(\bar{q},\dot{\bar{q}})}{\partial \dot{q}}\dot{\eta}(t) + \frac{\partial L(\bar{q},\dot{\bar{q}})}{\partial q}\eta(t).
$$

如前一节一样，我们使用$S$的$\eta(t)$幂次的泛函展开来推导经典路径的表达式。为此，我们将拉格朗日量代入作用量的表达式 (A.2.1) 中。接下来，我们将粒子可能的路径写为实际路径和修正$\eta(t)$之和。如前所述，我们使用分部积分并利用$\eta(t)$在积分边界处为零的事实。由此可得，作用量取极值的条件是：

$$
\int \mathrm{d}t \left[-\frac{d}{\mathrm{d}t}\left(\frac{\partial \mathcal{L}(\bar{q},\dot{\bar{q}})}{\partial \dot{q}}\right) + \frac{\partial \mathcal{L}(\bar{q},\dot{\bar{q}})}{\partial q}\right]\eta(t) = 0,
\tag{A.2.3}
$$

该式对任意$\eta(t)$成立当且仅当：

$$
-\frac{d}{\mathrm{d}t}\left(\frac{\partial \mathcal{L}(\bar{q},\dot{\bar{q}})}{\partial \dot{q}}\right) + \frac{\partial \mathcal{L}(\bar{q},\dot{\bar{q}})}{\partial q} = 0.
\tag{A.2.4}
$$

这就是拉格朗日运动方程。为了将此运动方程写成更熟悉的形式，我们引入与广义坐标$q$关联的广义动量$p$：

$$
p \equiv \frac{\partial \mathcal{L}(q,\dot{q})}{\partial \dot{q}}.
\tag{A.2.5}
$$

将此表达式代入式 (A.2.4) 得到：

$$
\dot{p} = \frac{\partial \mathcal{L}(q,\dot{q})}{\partial q}.
\tag{A.2.6}
$$

由于上述表述对任何坐标系都成立，它当然对笛卡尔坐标也成立。在这些坐标中，拉格朗日量为：

$$
L(x,\dot{x}) = \frac{1}{2}m\dot{x}^2 - U(x).
$$

与$x$关联的动量为：

$$
p_x = \frac{\partial L(x,\dot{x})}{\partial \dot{x}} = m\dot{x}
$$

运动方程为：

$$
m\ddot{x} = -\frac{\partial U(x)}{\partial x},
$$

这确实是我们从牛顿运动方程得到的结果。

???+ example "例证 24（重力场中的摆）"

    考虑长度为$l$、质量为$m$的简单摆（见图 A.1）。均匀重力场作用在摆上，势能是摆与竖直方向所成角度$\theta$的简单函数：

    $$
    U(\theta) = mgl[1 - \cos(\theta)].
    $$

    我们希望用广义坐标$\theta$表达运动方程。拉格朗日量$L$为：

    $$
    L = K - U = \frac{1}{2}m\left[\dot{x}^2(t) + \dot{y}^2(t)\right] - U(\theta) = \frac{ml^2}{2}\dot{\theta}^2 - U(\theta).
    $$

    广义动量定义为：

    $$
    p_\theta = \frac{\partial L}{\partial \dot{q}} = ml^2\dot{\theta}
    $$

    运动方程由式 (A.2.6) 得到：

    $$
    \dot{p}_\theta = -\frac{\partial U(\theta)}{\partial \theta}
    $$

    或

    $$
    \ddot{\theta} = -\frac{1}{ml^2}\frac{\partial U(\theta)}{\partial \theta}.
    $$

![图 A.1](../images/fig_A_1.png)

*图 A.1　长度为$l$、质量为$m$的简单摆。*

## 哈密顿量

使用拉格朗日量，我们得到以$q$和$\dot{q}$表示的运动方程。通常，以$q$及其共轭动量$p$表示运动方程更为方便。为此，我们可以进行勒让德变换[^2]：

$$
\mathcal{H}(q,p) \equiv p\dot{q} - \mathcal{L}(q,\dot{q},t).
\tag{A.3.1}
$$

该方程定义了系统的哈密顿量$H$。由于$H$是$q$、$p$的函数，通常也是$t$的函数，显然我们可以将$H$的微小变化写为：

$$
\mathrm{d}\mathcal{H}(q,p) = \frac{\partial \mathcal{H}}{\partial p}\mathrm{d}p + \frac{\partial \mathcal{H}}{\partial q}\mathrm{d}q + \frac{\partial \mathcal{H}}{\partial t}\mathrm{d}t.
\tag{A.3.2}
$$

但是，利用$H$的定义，我们也可以写：

$$
\begin{aligned}
\mathrm{d}H(q,p) &= d(p\dot{q}) - \mathrm{d}L(q,\dot{q}) \\
&= p\,\mathrm{d}\dot{q} + \dot{q}\,\mathrm{d}p - \left[\frac{\partial L}{\partial q}\mathrm{d}q + \frac{\partial L}{\partial \dot{q}}\mathrm{d}\dot{q} + \frac{\partial L}{\partial t}\mathrm{d}t\right] \\
&= p\,\mathrm{d}\dot{q} + \dot{q}\,\mathrm{d}p - \dot{p}\,\mathrm{d}q - p\,\mathrm{d}\dot{q} - \frac{\partial L}{\partial t}\mathrm{d}t \\
&= \dot{q}\,\mathrm{d}p - \dot{p}\,\mathrm{d}q - \frac{\partial L}{\partial t}\mathrm{d}t,
\end{aligned}
$$

其中我们分别使用了$p$和$\dot{p}$的定义，即式 (A.2.5) 和 (A.2.6)。由此直接得到：

$$
\begin{align}
\frac{\partial \mathcal{H}}{\partial p} &= \dot{q} \tag{A.3.3}\\
\frac{\partial \mathcal{H}}{\partial q} &= -\dot{p}.
\tag{A.3.4}
\end{align}
$$

这就是以$q$、$p$表示的所期望的运动方程。对于我们本书中考虑的大多数系统，拉格朗日量不显含时间。在这些情况下，哈密顿量是守恒的。这可以从运动方程直接得出：

$$
\frac{\mathrm{d}H(q,p)}{\mathrm{d}t} = \frac{\partial H}{\partial p}\dot{p} + \frac{\partial H}{\partial q}\dot{q} = -\frac{\partial H}{\partial p}\frac{\partial H}{\partial q} + \frac{\partial H}{\partial q}\frac{\partial H}{\partial p} = 0.
$$

这个守恒定律表达了在封闭系统中总能量守恒的事实。在笛卡尔坐标中，哈密顿量可以写为：

$$
H(x,p_x) = \dot{x}p_x - L(x,\dot{x}) = m\dot{x}^2 - \frac{1}{2}m\dot{x}^2 + U(x) = \frac{1}{2m}p_x^2 + U(x),
$$

哈密顿运动方程简化为牛顿方程：

$$
\begin{aligned}
\dot{x} &= \frac{\partial H}{\partial p_x} = \frac{p_x}{m} \\
\dot{p}_x &= -\frac{\partial H}{\partial x} = -\frac{\partial U(x)}{\partial x}.
\end{aligned}
$$

哈密顿运动方程是两个一阶微分方程——一个是关于$p$的，另一个是关于$q$的。相比之下，拉格朗日形式产生一个二阶微分方程。然而，两种形式产生相同的结果。两者之间的选择取决于数学便利性的考虑。

???+ example "例 28（重力场中的摆：第二部分）"

    我们再次考虑均匀重力场中的简单摆，如例 24 中所引入的：

    $$
    U(\theta) = mgl[1 - \cos(\theta)],
    $$

    其中$\theta$是摆与竖直方向的夹角，$g$是重力加速度。

    在例 24 中，我们从拉格朗日量导出了关于$\theta$的二阶微分方程形式的运动方程。现在我们将使用哈密顿表述。

    拉格朗日量为：

    $$
    L(\theta,\dot{\theta}) = U_K - U_P = \frac{ml^2}{2}\dot{\theta}^2 - U(\theta).
    $$

    拉格朗日量依赖于变量$\theta$和$\dot{\theta}$，而在哈密顿语言中我们希望以$\theta$及其共轭动量$p_\theta$表达运动方程。该共轭动量由式 (A.2.5) 定义：

    $$
    p_\theta \equiv \frac{\partial L(\theta,\dot{\theta})}{\partial \dot{\theta}} = ml^2\dot{\theta}.
    $$

    哈密顿量由勒让德变换 (A.3.1) 得到：

    $$
    H = p_\theta\dot{\theta} - L(\theta,\dot{\theta}) = \frac{p_\theta^2}{2ml^2} + U(\theta) = \frac{1}{2}ml^2\dot{\theta}^2 + U(\theta),
    $$

    当然，这等于摆的总能量。

    运动方程由式 (A.3.3) 和 (A.3.4) 得到：

    $$
    \begin{aligned}
    \dot{\theta} &= \frac{\partial H}{\partial p_\theta} = \frac{P_\theta}{ml^2} \\
    \dot{p}_\theta &= -\frac{\partial H}{\partial \theta} = -\frac{\mathrm{d}U(\theta)}{\mathrm{d}\theta},
    \end{aligned}
    $$

    这就是以两个一阶微分方程表示的所期望的运动方程。

## 哈密顿动力学与统计力学

在经典力学的哈密顿表述和拉格朗日表述之间的选择取决于便利性的考虑。拉格朗日形式更方便的一个例子是在推导带有约束的系统的运动方程时（见第 14.1 节）。另一方面，当建立与统计力学的联系时（见第 2 章），应使用哈密顿表达式。

### 正则变换

在哈密顿表述中，广义坐标和动量是独立变量。因此可以同时引入两个变量的变换。例如，坐标$q$、$p$到$Q$、$P$的变换记为：

$$
\begin{align}
Q &= Q(q,p) \\
P &= P(q,p)
\tag{A.4.1}
\end{align}
$$

逆变换，$Q$、$P$到$q$、$p$，记为：

$$
\begin{align}
q &= q(Q,P) \\
p &= p(Q,P).
\tag{A.4.2}
\end{align}
$$

显然，相空间坐标的任何函数的值不受坐标变换的影响。对于哈密顿量，这意味着：

$$
\mathcal{H}(q,p) \equiv \mathcal{H}[Q(p,q),P(q,p)] \equiv \mathcal{H}'(Q,P).
\tag{A.4.3}
$$

通常，新坐标中的运动方程不具有正则形式，除非坐标变换是正则的[^3]。如果坐标变换是正则的，新相空间坐标$Q$、$P$的运动方程为：

$$
\begin{align}
\dot{Q} &= \left(\frac{\partial \mathcal{H}'(Q,P)}{\partial P}\right) \tag{A.4.4}\\
\dot{P} &= -\left(\frac{\partial \mathcal{H}'(Q,P)}{\partial Q}\right).
\tag{A.4.5}
\end{align}
$$

从式 (A.4.1) 和坐标$q$、$p$的哈密顿运动方程可得：

$$
\dot{Q} = \left(\frac{\partial Q(q,p)}{\partial q}\right)\dot{q} + \left(\frac{\partial Q(q,p)}{\partial p}\right)\dot{p} = \left(\frac{\partial Q(q,p)}{\partial q}\right)\left(\frac{\partial H(q,p)}{\partial p}\right) - \left(\frac{\partial Q(q,p)}{\partial p}\right)\left(\frac{\partial H(q,p)}{\partial q}\right).
$$

利用式 (A.4.3)，我们可以写：

$$
\left(\frac{\partial H'(Q,P)}{\partial P}\right) = \left(\frac{\partial H(q,p)}{\partial p}\right)\left(\frac{\partial p(P,Q)}{\partial P}\right) + \left(\frac{\partial H(q,p)}{\partial q}\right)\left(\frac{\partial q(P,Q)}{\partial P}\right).
$$

该方程只有在以下条件满足时才能与$\dot{Q}$的表达式 (A.4.4) 相等：

$$
\begin{align}
\left(\frac{\partial Q(q,p)}{\partial q}\right) &= \left(\frac{\partial p(Q,P)}{\partial P}\right) \\
\left(\frac{\partial Q(q,p)}{\partial p}\right) &= -\left(\frac{\partial q(Q,P)}{\partial P}\right).
\tag{A.4.6}
\end{align}
$$

类似地，我们可以从$\dot{P}$出发，推导出另外两个条件：

$$
\begin{align}
\left(\frac{\partial P(q,p)}{\partial q}\right) &= -\left(\frac{\partial p(Q,P)}{\partial Q}\right) \\
\left(\frac{\partial P(q,p)}{\partial p}\right) &= \left(\frac{\partial q(Q,P)}{\partial Q}\right).
\tag{A.4.7}
\end{align}
$$

这两个方程定义了正则变换的条件。

### 辛条件

我们可以通过使用矩阵记号将上述正则变换的条件表达为单个方程。设$\boldsymbol{\xi}$为包含$N$个粒子在$d$维中的广义坐标$q_i$和动量$p_i$的$2\mathrm{d}N$维向量（见第 2.5.1 节）。哈密顿运动方程 (A.3.3) 和 (A.3.4) 可以写为：

$$
\dot{\boldsymbol{\xi}} = \boldsymbol{\omega}\frac{\partial \mathcal{H}}{\partial \boldsymbol{\xi}},
\tag{A.4.8}
$$

其中$\boldsymbol{\omega}$是反对称矩阵，定义为：

$$
\boldsymbol{\omega} = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}.
$$

类似地，我们可以定义$\boldsymbol{\xi}$为包含广义坐标$Q_i$和$P_i$的$2N$维向量。使用矩阵记号，从$Q$、$P$到$q$、$p$的变换 (A.4.1) 写为：

$$
\boldsymbol{\xi} = \boldsymbol{\xi}(\boldsymbol{\zeta}).
$$

对于$\boldsymbol{\xi}$的时间导数，我们可以写：

$$
\dot{\boldsymbol{\xi}} = \mathbf{M}\dot{\boldsymbol{\zeta}},
$$

其中$\mathbf{M}$是变换的雅可比矩阵。该矩阵的元素为：

$$
M_{ij} = \frac{\partial \xi_i}{\partial \zeta_j}.
\tag{A.4.9}
$$

利用式 (A.4.8)，我们可以写$\boldsymbol{\xi}$的时间导数为：

$$
\dot{\boldsymbol{\xi}} = \mathbf{M}\boldsymbol{\omega}\frac{\partial \mathcal{H}}{\partial \boldsymbol{\zeta}}.
\tag{A.4.10}
$$

类似地，我们可以定义逆变换 (A.4.2)：

$$
\boldsymbol{\zeta} = \boldsymbol{\zeta}(\boldsymbol{\xi}).
$$

由于$H(\mathbf{p},\mathbf{q}) = H(\mathbf{P},\mathbf{Q})$，我们可以写：

$$
\frac{\partial \mathcal{H}(\boldsymbol{\zeta})}{\partial \zeta_i} = \sum_j \frac{\partial \mathcal{H}(\boldsymbol{\xi})}{\partial \xi_j}\frac{\partial \xi_j}{\partial \zeta_i}.
\tag{A.4.11}
$$

如果我们定义式 (A.4.9) 中定义的$\mathbf{M}$的转置矩阵[^4]：

$$
\tilde{M}_{ij} = \frac{\partial \xi_j}{\partial \zeta_i}.
$$

这使我们可以将式 (A.4.11) 用矩阵记号重写为：

$$
\frac{\partial \mathcal{H}(\boldsymbol{\zeta})}{\partial \boldsymbol{\zeta}} = \tilde{\mathbf{M}}\frac{\partial \mathcal{H}(\boldsymbol{\xi})}{\partial \boldsymbol{\xi}}.
\tag{A.4.12}
$$

结合式 (A.4.10) 和 (A.4.12)，我们有：

$$
\dot{\boldsymbol{\xi}} = \mathbf{M}\boldsymbol{\omega}\tilde{\mathbf{M}}\frac{\partial H}{\partial \boldsymbol{\xi}}.
$$

这个运动方程的表达式对于任何从$\boldsymbol{\zeta}$集合变换（与时间无关地）得到的变量集$\boldsymbol{\xi}$都是有效的。这种变换在新坐标中具有正则形式的运动方程时是正则的：

$$
\dot{\boldsymbol{\xi}} = \boldsymbol{\omega}\frac{\partial H}{\partial \boldsymbol{\xi}}.
$$

这只有在$\mathbf{M}$满足以下条件时才成立：

$$
\mathbf{M}\boldsymbol{\omega}\tilde{\mathbf{M}} = \boldsymbol{\omega}.
\tag{A.4.13}
$$

这个条件通常被称为辛条件。满足此条件的矩阵$\mathbf{M}$称为辛矩阵[^5]。

### 统计力学

使用正则变换的辛记号，我们考虑其对统计力学的意义。在微正则系综中，三维原子系统的经典配分函数$\Omega$定义为：

$$
\Omega_{N,V,E} = \frac{1}{h^{3N}N!}\int \mathrm{d}p^N \mathrm{d}q^N\,\delta(\mathcal{H}(\mathbf{p},\mathbf{q}) - E),
\tag{A.4.14}
$$

其中$h$是普朗克常数，$\delta$函数将积分限制在由$H(\mathbf{p},\mathbf{q}) = E$定义的相空间超曲面上。我们可以用其他相空间坐标重新表达这个积分，但此时必须考虑两种坐标集中的体积元不一定相同。与$\boldsymbol{\zeta}$关联的体积元为：

$$
\mathrm{d}\boldsymbol{\zeta} = \mathrm{d}q_1... \mathrm{d}q_N\,\mathrm{d}p_1... \mathrm{d}p_N
$$

与$\boldsymbol{\xi}$关联的为：

$$
\mathrm{d}\boldsymbol{\xi} = \mathrm{d}Q_1... \mathrm{d}Q_N\,\mathrm{d}P_1... \mathrm{d}P_N.
$$

这两个体积元通过变换矩阵的雅可比矩阵相关联：

$$
\mathrm{d}\boldsymbol{\zeta} = |\mathrm{Det}(\mathbf{M})|\,\mathrm{d}\boldsymbol{\xi}.
\tag{A.4.15}
$$

该方程表明，通常坐标变换将导致配分函数中出现雅可比行列式：

$$
\Omega_{N,V,E} = \frac{1}{h^{3N}N!}\int \mathrm{d}P^N \mathrm{d}Q^N\,|\mathrm{Det}(\mathbf{M})|\,\delta\left[\mathcal{H}'(\mathbf{P},\mathbf{Q}) - E\right].
\tag{A.4.16}
$$

在计算非原始笛卡尔坐标系中的系综平均时，变换的雅可比矩阵 $\mathbf{M}$可能不等于 1，应当加以考虑。在下文中，我们用符号$\omega$表示雅可比行列式 $|\mathrm{Det}(\mathbf{M})|$。

对于正则变换，即满足条件 (A.4.13) 的变换，雅可比行列式的绝对值为 1。为了推导这个结果，我们对辛条件 (A.4.13) 两边取行列式：

$$
\begin{aligned}
\mathrm{Det}(\mathbf{M}\boldsymbol{\omega}\tilde{\mathbf{M}}) &= \mathrm{Det}(\boldsymbol{\omega}) \\
\mathrm{Det}^2(\mathbf{M})\mathrm{Det}(\boldsymbol{\omega}) &= \mathrm{Det}(\boldsymbol{\omega}).
\end{aligned}
$$

这个方程只有在$\mathbf{M}$的行列式为$\pm 1$时才成立，这意味着对于正则变换，与该变换关联的雅可比行列式的绝对值必须为 1。

经典系统在相空间中的自然时间演化可以被视为一种坐标变换：

$$
\boldsymbol{\zeta}(t_0) \to \boldsymbol{\zeta}(t).
$$

哈密顿系统的一个重要性质是自然时间演化对应于辛坐标变换。我们可以将$\boldsymbol{\zeta}(t_0)$到$\boldsymbol{\zeta}(t)$的变换视为具有时间步$\delta t$的无穷小变换序列。假设我们定义时间间隔$\delta t$内坐标的演化为从$\boldsymbol{\zeta}$到$\boldsymbol{\xi}$的坐标变换：

$$
\boldsymbol{\xi} = \boldsymbol{\xi}(\boldsymbol{\zeta}) = \boldsymbol{\zeta}(t + \delta t) = \boldsymbol{\zeta}(t) + \dot{\boldsymbol{\zeta}}(t)\delta t.
$$

此变换的雅可比行列式为：

$$
\mathbf{M} \equiv \frac{\partial \boldsymbol{\xi}}{\partial \boldsymbol{\zeta}} = \mathbf{1} + \delta t\frac{\partial}{\partial \boldsymbol{\zeta}}\left[\boldsymbol{\omega}\frac{\partial H}{\partial \boldsymbol{\zeta}}\right] = \mathbf{1} + \delta t\,\boldsymbol{\omega}\frac{\partial^2 H}{\partial \boldsymbol{\zeta}\partial \boldsymbol{\zeta}},
$$

其中：

$$
\left[\frac{\partial^2 H}{\partial \boldsymbol{\zeta}\partial \boldsymbol{\zeta}}\right]_{ij} = \frac{\partial^2 H}{\partial \zeta_i \partial \zeta_j}.
$$

考虑到$\boldsymbol{\omega}$是反对称矩阵，我们可以写出矩阵$\mathbf{M}$的转置：

$$
\tilde{\mathbf{M}} = \mathbf{1} - \frac{\partial^2 H}{\partial \boldsymbol{\zeta}\partial \boldsymbol{\zeta}}\boldsymbol{\omega}.
$$

将雅可比行列式的这个表达式代入辛条件 (A.4.13) 得到（在$\delta t$的一阶近似下）：

$$
\begin{aligned}
\mathbf{M}\boldsymbol{\omega}\tilde{\mathbf{M}} &= \left[\mathbf{1} + \delta t\,\boldsymbol{\omega}\frac{\partial^2 H}{\partial \boldsymbol{\zeta}\partial \boldsymbol{\zeta}}\right]\boldsymbol{\omega}\left[\mathbf{1} - \delta t\frac{\partial^2 H}{\partial \boldsymbol{\zeta}\partial \boldsymbol{\zeta}}\boldsymbol{\omega}\right] \\
&\approx \boldsymbol{\omega} + \delta t\,\boldsymbol{\omega}\frac{\partial^2 H}{\partial \boldsymbol{\zeta}\partial \boldsymbol{\zeta}}\boldsymbol{\omega} - \boldsymbol{\omega}\,\delta t\frac{\partial^2 H}{\partial \boldsymbol{\zeta}\partial \boldsymbol{\zeta}}\boldsymbol{\omega} \\
&= \boldsymbol{\omega}.
\end{aligned}
$$

因此辛条件在无穷小时间间隔内$\boldsymbol{\zeta}$的演化中成立。由于我们可以将有限时间间隔内$\boldsymbol{\zeta}$的演化视为无穷小步长正则变换的序列，总的时间演化也满足辛条件。

可以将哈密顿量视为作用于相空间所有点的正则变换的生成元。由于正则变换的雅可比行列式等于 1，相空间中体积元的大小在哈密顿系统的自然时间演化过程中不变。此外，相空间中任意点周围的密度$f(\mathbf{q}(t),\mathbf{p}(t))$在时间演化过程中也保持不变。要理解这一点，考虑相空间中由曲面$S$包围的体积$V$。在时间演化过程中，曲面移动，曲面内的所有点也随之移动。然而，点不能穿过曲面。原因很简单：如果相空间中的两条轨迹相交，将意味着两条轨迹从同一相空间点出发。但这是不可能的，因为这将意味着从该点出发的轨迹不由其初始条件唯一确定。因此，任意体积内的相空间点数不随时间变化。由于体积本身也是恒定的，这意味着相空间密度（即单位体积的点数）是恒定的。换言之：哈密顿系统的相空间密度表现得像不可压缩流体：

$$
\frac{\mathrm{d}f}{\mathrm{d}t} = 0.
\tag{A.4.17}
$$

虽然哈密顿运动方程的精确解将满足不可压缩性条件，但离散的数值格式——通常——会违反它。如前所述，我们可以将任何数值 MD 算法（例如 Verlet、速度 Verlet 等）视为从$(\mathbf{q}(t),\mathbf{p}(t))$到$(\mathbf{q}(t + \Delta t),\mathbf{q}(t + \Delta t))$的变换。然后我们可以计算此变换的雅可比行列式，并检查它是否等于 1（见第 4.3 节和第 4.3.4 节）。对于所有求解牛顿运动方程的“好的”算法，从$(\mathbf{q}(t),\mathbf{p}(t))$到$(\mathbf{q}(t + \Delta t),\mathbf{q}(t + \Delta t))$的变换的雅可比行列式等于 1——这种算法被称为“保面积”的。应当注意，辛条件所蕴含的不仅仅是保面积性质。不幸的是，这些其他后果没有如此简单的直观解释。当我们说一个算法应该是辛的时候，我们的意思不仅仅是它应该是保面积的——它应该真正满足辛条件。幸运的是，在许多情况下，利用任何一组经典哈密顿运动方程都满足辛条件这一事实，算法的辛性质很容易证明。可以写成由简单哈密顿量生成的精确时间演化序列的算法，因此必然是辛的。一个例子是 Verlet 算法。正如第 4.3.4 节中所讨论的，该算法可以被视为使用哈密顿量的动能部分或势能部分的一系列精确传播。两种传播都满足辛条件。因此，Verlet 算法整体上是辛的。关于辛动力学的通俗讨论，参见文献[[713]](references.md#ref-713)。关于分子动力学模拟中辛积分器的讨论可在文献[[714]](references.md#ref-714)中找到。

---

[^1]: 正确定义更为严格；详见[[54]](references.md#ref-54)。
[^2]: 在热力学中，勒让德变换用于推导各种热力学势。例如，能量$E$是熵$S$和体积$V$的自然函数：$E = E(S,V)$，即在这些变量中，$E$是热力学势。在大多数实际应用中，以温度$T$而非熵$S$作为独立变量更为方便。由于温度是与熵共轭的变量（$\partial E/\partial S = T$），我们可以进行勒让德变换来消除$S$的依赖性：
$A \equiv E - TS$，
得到：
$\mathrm{d}A = \mathrm{d}E - d(TS) = -S\mathrm{d}T - p\mathrm{d}V$。
由于历史原因，将拉格朗日量与哈密顿量联系起来的勒让德变换具有相反的符号。
[^3]: 由于我们假设时间在这些方程中不显式出现，我们定义的是所谓的受限正则变换。
[^4]: 给定矩阵$\mathbf{A}$的转置矩阵可以通过交换行和列得到，即$\tilde{a}_{ij} = a_{ji}$。
[^5]: 要看出此条件与式 (A.4.6) 和 (A.4.7) 等价，我们必须从右边将此方程乘以$\tilde{\mathbf{M}}$的逆矩阵：
$\mathbf{M}\boldsymbol{\omega} = \boldsymbol{\omega}\tilde{\mathbf{M}}^{-1}$。