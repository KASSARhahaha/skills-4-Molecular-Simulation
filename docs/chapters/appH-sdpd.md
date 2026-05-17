# 光滑耗散粒子动力学

## Navier-Stokes方程与Fourier定律

我们简要总结支撑Español和Revenga 的光滑耗散粒子动力学（SDPD）方法的最重要关系式。SDPD的描述适合这本基于粒子的模拟书籍，因为它提供了流体输运的基于粒子的图像，该图像与非平衡热力学  兼容，正确地考虑了熵产生。此外，与DPD一样，它正确地考虑了热涨落。

如正文所述，该方法是对20世纪70年代光滑粒子流体动力学方法  的扩展，旨在通过一组粒子轨迹来表示连续的流体动力学流场，这些粒子遵循运动方程，在足够大的尺度上模拟粘性流的Navier-Stokes方程，并结合热传导方程。

流体"粒子"$i$的能量通过形如$E_i = E(m, S_i, V_i)$的状态方程与基本热力学参数相关联，该方程仍待确定。每个粒子具有温度$T_i$和压力$P_i$，由下式给出：

$$
T_i = \left( \frac{\partial E}{\partial S} \right)_V \quad \text{和} \quad P_i = -\left( \frac{\partial E}{\partial V} \right)_S .
$$

在不存在热涨落的情况下，宏观输运由Navier-Stokes方程描述：

$$
m\rho \frac{d\mathbf{v}}{dt} = -\nabla P + \eta \nabla^2 \mathbf{v} + \left(\zeta + \frac{\eta}{3}\right) \nabla (\nabla \cdot \mathbf{v}) ,
\quad (H.1.1)
$$

并辅以由于粘性流动和热输运导致的熵产生的（简化）表达式：

$$
T\rho \frac{ds}{dt} = \phi + \kappa \nabla^2 T
\quad (H.1.2)
$$

其中

$$
\phi = 2\eta \nabla \mathbf{v} : \nabla \mathbf{v} + \zeta (\nabla \cdot \mathbf{v})^2 ,
$$

最后是质量守恒定律：

$$
\frac{d\rho}{dt} = -\rho \nabla \cdot \mathbf{v} .
\quad (H.1.3)
$$

在上述方程中，$\rho$表示数密度：为与本书其他部分的记号保持一致，我们使用的记号与文献  略有不同。$\mathbf{v}$表示流动速度，$P$表示静水压力，$\eta$表示剪切粘度，$\zeta$表示体积粘度，$s$表示每个粒子的熵，$\kappa$表示热导率。

## 离散化SDPD方程

式(H.1.1)至(H.1.3)的离散化形式为：

$$
m\dot{\mathbf{v}}_i = -\frac{(\nabla P)_i}{\rho_i} + \frac{\eta (\nabla^2 \mathbf{v})_i}{\rho_i} + \frac{(\zeta + \eta/3)(\nabla \nabla \cdot \mathbf{v})_i}{\rho_i} ,
\quad (H.2.1)
$$

$$
T_i \dot{S}_i = \frac{\phi + \kappa(\nabla^2 T)_i}{\rho_i} ,
\quad (H.2.2)
$$

以及

$$
\dot{\rho}_i = -\rho_i (\nabla \cdot \mathbf{v})_i ,
\quad (H.2.3)
$$

其中我们现在使用密度、压力和速度梯度的局部表达式。然而，为了将连续形式的Navier-Stokes方程联系起来，我们现在必须给出在粒子$i$附近计算梯度项的规定。

为了在Navier-Stokes方程和基于粒子的描述之间建立这种联系，Español和Revenga假设了粒子$i$的体积$V_i$与该粒子周围的局部密度$\rho_i$之间的如下关系：

$$
V_i \equiv 1/\rho_i
$$

其中

$$
\rho_i \equiv \int \mathrm{d}\mathbf{r} \, \rho(\mathbf{r}) W(\mathbf{r} - \mathbf{r}_i) ,
$$

其中$\rho(\mathbf{r}) \equiv \sum_j \delta(\mathbf{r} - \mathbf{r}_j)$，而$W(\mathbf{r} - \mathbf{r}_i)$是归一化的权重函数：

$$
\int \mathrm{d}\mathbf{r} \, W(\mathbf{r} - \mathbf{r}_i) = 1 .
$$

注意粒子$i$周围的密度包含$i$自身：因此$\rho_i > 0$。函数形式$W(\mathbf{r})$有多种选择；选择具有有限支撑的$W$形式较为方便，即对于$r > h$有$W(r) = 0$，其中$h$仍待确定。文献  使用了：

$$
W(r) = \frac{105}{16\pi h^3} \left(1 - \frac{r}{h}\right)^3 \left(1 + 3\frac{r}{h}\right) \quad \text{当} \quad r \leq h
$$

$$
= 0 \quad \text{其他情况} .
$$

在以下的推导中，进入方程的不是$W(r)$本身，而是通过以下关系定义的量$F(r)$：

$$
\nabla W(r) = -\mathbf{r} F(r) .
$$

对于文献  的选择，

$$
F(r) = \frac{315}{4\pi h^5} \left(1 - r/h\right)^2 .
$$

对于$F(\mathbf{r}_i - \mathbf{r}_j)$我们将使用简写符号$F_{ij}$，类似地$T_{ij} \equiv T_i - T_j$，$\mathbf{r}_{ij} \equiv (\mathbf{r}_i - \mathbf{r}_j)$，$\mathbf{v}_{ij} \equiv (\mathbf{v}_i - \mathbf{v}_j)$以及$\hat{\mathbf{e}}_{ij} \equiv \mathbf{r}_{ij}/|\mathbf{r}_{ij}|$。

由此可得（经过一些略微繁琐的代数运算——见文献 ）：

$$
m\dot{\mathbf{v}}_i = \sum_j \left[
\left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) F_{ij} \mathbf{r}_{ij}
- \frac{5\eta}{\rho_i \rho_j} F_{ij} \mathbf{v}_{ij}
- \frac{\zeta + \eta/3}{\rho_i \rho_j} F_{ij} (\mathbf{r}_{ij} \cdot \mathbf{v}_{ij}) \hat{\mathbf{e}}_{ij}
\right] ,
\quad (H.2.4)
$$

式(H.2.4)与

$$
\dot{\mathbf{r}}_i = \mathbf{v}_i ,
$$

一起提供了流体"粒子"位置随时间演化的显式表达式。重要的是，该方程满足牛顿第三定律。

我们仍需要给出熵产生的离散化形式：

$$
T_i \dot{S}_i = (\phi)_i - 2\kappa \sum_j \frac{F_{ij}}{\rho_i \rho_j} T_{ij} ,
\quad (H.2.5)
$$

其中

$$
(\phi)_i = \sum_j \frac{F_{ij}}{\rho_i \rho_j} \left[
\left( \frac{5\eta}{2} - \zeta + \frac{\zeta + \eta/3}{3} \right) \mathbf{v}_{ij}^2
+ 5 \left(\frac{\eta}{2} + \frac{\zeta}{3}\right) (\hat{\mathbf{e}}_{ij} \cdot \mathbf{v}_{ij})^2
\right] .
$$

为了求解这些方程，我们需要一个状态方程来在每个位置将压力与密度联系起来：这就是我们的关系式$E_i = E(m, S_i, V_i)$，我们可以选择它来重现所考虑流体的性质。

## 添加噪声

唯一剩下的事情是添加噪声项，以考虑流体中的热涨落。这涉及以正确的幅度向粒子的熵和速度添加随机噪声，同时满足能量守恒和动量守恒。有关细节，我们请读者参阅文献 。