# 非哈密顿动力学

将经典统计力学方法推广到非哈密顿系统的系统性方法由 Tuckerman 等人[[267,715,716]](references.md#ref-267)提出。在本附录中，我们简要概述分析扩展拉格朗日系统的一般方法。然而，我们将跳过大部分推导。对于使用微分几何数学方法的更完整、更严格的推导，读者可参考原始文献。

一般而言，求解非哈密顿运动方程所产生的动力学不是保面积的。正如我们在附录 A.4.3 中所看到的，求解运动方程可以被视为从时间$t_0$处的相空间坐标到时间$t$处的相空间坐标的坐标变换。如果系统是哈密顿的，系统的时间演化将改变相空间中无穷小体积元的形状，但不会改变其体积$\mathrm{d}\boldsymbol{\Gamma}$。相比之下，对于非哈密顿系统，我们必须考虑与$\mathrm{d}\boldsymbol{\Gamma}(t_0) \rightarrow \mathrm{d}\boldsymbol{\Gamma}(t)$演化相关的变换的雅可比行列式：

$$
\mathrm{d}\boldsymbol{\Gamma}_t = J(\boldsymbol{\Gamma}_t ; \boldsymbol{\Gamma}_0)\,\mathrm{d}\boldsymbol{\Gamma}_0,
$$

其中下标$0$表示$t = t_0$时的相空间体积，$J$是变换的雅可比矩阵$\mathbf{M}$的行列式。为方便起见，我们选择$t_0 = 0$。

哈密顿系统在相空间中的运动类似于不可压缩液体的运动：这种“液体”的体积随时间不发生变化。相比之下，非哈密顿系统是可压缩的。在考虑将刘维尔定理推广到非哈密顿系统时，必须考虑这种可压缩性。

可压缩性可以从雅可比行列式的时间依赖性推导出来：

$$
\frac{\mathrm{d}J(\boldsymbol{\Gamma}_t ; \boldsymbol{\Gamma}_0)}{\mathrm{d}t} = \kappa(\boldsymbol{\Gamma}_t, t)\,J(\boldsymbol{\Gamma}_t ; \boldsymbol{\Gamma}_0),
\tag{B.0.1}
$$

其中$\kappa(\boldsymbol{\Gamma}_t, t)$，即动力学系统的相空间压缩率，定义为：

$$
\kappa(\boldsymbol{\Gamma}_t, t) \equiv \nabla_{\boldsymbol{\Gamma}} \cdot \dot{\boldsymbol{\Gamma}}.
\tag{B.0.2}
$$

式 (B.0.1) 的解为

$$
J(\boldsymbol{\Gamma}_t ; \boldsymbol{\Gamma}_0) = \exp\!\left[\int_0^t \kappa(\boldsymbol{\Gamma}_s, s)\,\mathrm{d}s\right].
$$

如果我们将$w(\boldsymbol{\Gamma}_t, t)$定义为$\kappa(\boldsymbol{\Gamma}_t, t)$的原函数，则可以写为

$$
J(\boldsymbol{\Gamma}_t ; \boldsymbol{\Gamma}_0) = \exp[w(\boldsymbol{\Gamma}_t, t) - w(\boldsymbol{\Gamma}_0, 0)] \equiv \frac{\sqrt{g(\boldsymbol{\Gamma}_0, 0)}}{\sqrt{g(\boldsymbol{\Gamma}_t, t)}},
$$

其中最后一行定义了量$g$。回顾

$$
\mathrm{d}\boldsymbol{\Gamma}_t = J(\boldsymbol{\Gamma}_t ; \boldsymbol{\Gamma}_0)\,\mathrm{d}\boldsymbol{\Gamma}_0 = \frac{\sqrt{g(\boldsymbol{\Gamma}_0, 0)}}{\sqrt{g(\boldsymbol{\Gamma}_t, t)}}\,\mathrm{d}\boldsymbol{\Gamma}_0,
$$

因此，

$$
\int \sqrt{g(\boldsymbol{\Gamma}_t, t)}\,\mathrm{d}\boldsymbol{\Gamma}_t = \int \sqrt{g(\boldsymbol{\Gamma}_0, 0)}\,\mathrm{d}\boldsymbol{\Gamma}_0,
$$

这定义了相空间中的一个不变测度。这一结果可以用来推导非哈密顿系统刘维尔方程的新形式。这里的关键点是，相空间分布$f(\boldsymbol{\Gamma})$（我们感兴趣的函数，给出相空间中的概率密度）应当与相空间度量$g$（确保非哈密顿系统的相空间体积在时间演化下不变）区分开来：

$$
\frac{\partial \left(f\sqrt{g}\right)}{\partial t} + \nabla \cdot \left(f\sqrt{g}\,\dot{\boldsymbol{\Gamma}}\right) = 0.
\tag{B.0.3}
$$

与系综平均对应的表达式为

$$
\langle A \rangle = \frac{\int \mathrm{d}\boldsymbol{\Gamma}\,\sqrt{g(\boldsymbol{\Gamma})}\,A(\boldsymbol{\Gamma})\,f(\boldsymbol{\Gamma})}{\int \mathrm{d}\boldsymbol{\Gamma}\,\sqrt{g(\boldsymbol{\Gamma})}\,f(\boldsymbol{\Gamma})}.
\tag{B.0.4}
$$

假设存在$n_c$个守恒定律$\Phi_k(\boldsymbol{\Gamma}') = C_k$（$k = 1, ..., n_c$），则非哈密顿系统的配分函数为

$$
\Xi(C_1, ..., C_{n_c}) = \int \mathrm{d}\boldsymbol{\Gamma}'\,\sqrt{g(\boldsymbol{\Gamma}')}\;\prod_{k=1}^{n_c}\delta\!\left[\Phi_k(\boldsymbol{\Gamma}') - C_k\right].
\tag{B.0.5}
$$

在许多应用中，可以通过对引入的、用于表示恒温器或恒压器效应的非物理变量进行积分，从上述“微正则”配分函数获得正确的（$NVT$或$NPT$）配分函数。为了正确地进行这一步，必须识别所有的守恒定律。此外，从分析中消除所有线性依赖于其他变量的坐标和“被驱动”的变量是有用的。当一个变量不影响（且不通过守恒定律耦合）系统中感兴趣的物理变量的时间演化时，即使其自身的时间演化可能依赖于这些物理变量，该变量就被称为“被驱动”的。