# Kirkwood-Buff 理论

## 混合物的结构因子

在第 5.1.7.1 节中，我们讨论了单组分系统的结构因子$S(\mathbf{q})$（式 (5.1.40)）与粒子密度傅里叶变换的均方值（式 (5.1.38)）之间的关系：

$$
\rho(\mathbf{q}) = \sum_{i=1}^{N} \mathrm{e}^{\mathrm{i}\mathbf{q}\cdot\mathbf{r}_i} = \int_V \mathrm{d}\mathbf{r}\,\rho(\mathbf{r})\,\mathrm{e}^{\mathrm{i}\mathbf{q}\cdot\mathbf{r}}.
$$

对于$n$组分系统，我们可以将这一关系推广，得到偏结构因子$S_{ab}(\mathbf{q})$的表达式，该因子度量了组分$a$和$b$的密度傅里叶变换涨落之间的交叉关联：

$$
\begin{align}
S_{ab}(\mathbf{q}) &= \frac{1}{\sqrt{\langle N_a \rangle \langle N_b \rangle}}\left\langle \delta\rho_a(\mathbf{q})\,\delta\rho_b(-\mathbf{q}) \right\rangle \nonumber \\
&= \frac{1}{\sqrt{\langle N_a \rangle \langle N_b \rangle}} \frac{1}{V^2}\int_V \mathrm{d}\mathbf{r}\int_V \mathrm{d}\mathbf{r}'\,\delta\rho_a(\mathbf{r})\,\delta\rho_b(\mathbf{r}')\,\mathrm{e}^{\mathrm{i}\mathbf{q}\cdot(\mathbf{r}-\mathbf{r}')},
\tag{C.1.1}
\end{align}
$$

其中$\delta\rho_i(\mathbf{q})$表示$\rho(\mathbf{q})$围绕其平均值的涨落。如果我们在式 (C.1.1b) 的后半部分取$q \to 0$的极限，可得：

$$
\lim_{q\to 0} S_{ab}(q) = \sqrt{\frac{\langle \delta N_a\,\delta N_b \rangle}{\langle N_a \rangle \langle N_b \rangle}},
\tag{C.1.2}
$$

其中$\delta N_a$表示系统中物种$a$的粒子总数的涨落。沿用 Kirkwood 和 Buff 的方法[[717]](references.md#ref-717)，我们将证明式 (C.1.2) 具有直接的热力学解释，这为确定溶液中化学势随组成变化的依赖关系提供了一条强有力的途径。

对于多组分系统，我们可以将巨正则配分函数$\Xi$（式 (6.3.1)）推广为

$$
\Xi(\{\mu\}, V, T) \equiv \sum_{N_1,N_2,\cdots,N_n=0}^{\infty} \prod_{a=1}^{n} \exp(\beta\mu_a N_a)\,\mathrm{e}^{-\beta F(\{N\},V,T)},
\tag{C.1.3}
$$

其中$\{\mu\}$表示$\mu_1, \mu_2, \cdots, \mu_n$，$\{N\}$代表$N_1, N_2, \cdots, N_n$。巨势$\Omega = \Omega(\{\mu\}, V, T)$由下式给出：

$$
\Omega = -k_B T \ln \Xi(\{\mu\}, V, T).
$$

由式 (C.1.3) 可得：

$$
\left(\frac{\partial \Omega}{\partial \mu_a}\right) = -\langle N_a \rangle
\tag{C.1.4}
$$

以及

$$
\frac{\partial^2 \Omega}{\partial \mu_a\, \partial \mu_b} = -\frac{\partial \langle N_a \rangle}{\partial \mu_b} = -\beta\left\langle \delta N_a\,\delta N_b \right\rangle_{\{\mu\},V,T}.
\tag{C.1.5}
$$

将式 (C.1.5) 与式 (C.1.2) 进行比较，可以看出结构因子$S_{ab}(\mathbf{q})$在$q \to 0$极限下的行为与巨势的热力学导数之间存在密切关系。

Kirkwood 和 Buff 首先提出了这些关系[[717]](references.md#ref-717)。然而，他们也将结果表示为对径向分布函数$g_{ab}(r)$的积分形式，许多模拟研究使用$g(r)$方法来计算化学势随组成的变化。基于$g(r)$的关系在原理上是正确的，但正如式 (5.1.42) 下方所解释的，在（甚至不非常小的）小系统模拟中使用它们是非常危险的。因此，最好坚持使用式 (C.1.2) [[718]](references.md#ref-718)。[^1]

我们还需要解释为什么式 (C.1.5) 中的关系很重要。关键原因是，在粒子插入方法（见第 8.5.1 节）失效的条件下，它们使我们能够计算多组分混合物中各物种化学势随组成的变化。

我们首先注意到，在恒定$T$和$V$下：

$$
\mathrm{d}\Omega = -\sum_{r=1}^{n} \langle N_r \rangle\,\mathrm{d}\mu_r
$$

因此：

$$
\left(\frac{\partial \Omega}{\partial N_a}\right)_{T,V,N'} = -\sum_{r=1}^{n} \langle N_r \rangle \left(\frac{\partial \mu_r}{\partial N_a}\right)_{T,V,N'},
$$

其中$N'$和$\mu'$表示在该特定微分操作中保持不变的所有$\{N\}$和$\{\mu\}$的集合（即我们对不同集合使用相同的符号）。由此可得：

$$
\begin{align}
\frac{\partial^2 \Omega}{\partial N_a\,\partial N_b}\bigg|_{T,V,N'} &= \sum_{r,s=1}^{n} \frac{\partial^2 \Omega}{\partial \mu_r\,\partial \mu_s}\bigg|_{T,V,\mu'} \left(\frac{\partial \mu_r}{\partial N_a}\right)_{T,V,N'}\left(\frac{\partial \mu_s}{\partial N_b}\right)_{T,V,N'} \notag \\
&\quad - \sum_{r=1}^{n} \langle N_r \rangle \frac{\partial^2 \mu_r}{\partial N_a\,\partial N_b}\bigg|_{T,V,N'}.
\tag{C.1.6}
\end{align}
$$

由吉布斯-杜亥姆方程（$\mathrm{d}\Omega|_{V,T} = -\sum_r \langle N_r \rangle\,\mathrm{d}\mu_r$）可得：

$$
\frac{\partial^2 \Omega}{\partial N_a\,\partial N_b}\bigg|_{T,V,N'} = -\frac{\partial \mu_a}{\partial N_b}\bigg|_{T,V,N'} - \sum_{r=1}^{n} \langle N_r \rangle \frac{\partial^2 \mu_r}{\partial N_a\,\partial N_b}\bigg|_{T,V,N'}.
\tag{C.1.7}
$$

将式 (C.1.6) 和式 (C.1.7) 结合，我们得到：

$$
\frac{\partial \mu_a}{\partial N_b}\bigg|_{T,V,N'} = -\sum_{r,s=1}^{n} \frac{\partial^2 \Omega}{\partial \mu_r\,\partial \mu_s}\bigg|_{T,V,\mu'} \left(\frac{\partial \mu_r}{\partial N_a}\right)_{T,V,N'}\left(\frac{\partial \mu_s}{\partial N_b}\right)_{T,V,N'}.
\tag{C.1.8}
$$

注意式 (C.1.8) 具有矩阵方程的形式：

$$
A_{ab} = \sum_{r,s} A_{ar} B_{rs} A_{sb},
$$

其中

$$
A_{ab} = A_{ba} = \left(\frac{\partial \mu_a}{\partial N_b}\right)_{T,V,N'},
$$

以及

$$
B_{rs} = B_{sr} = \beta\left\langle \delta N_r\,\delta N_s \right\rangle_{\{\mu\},V,T}.
$$

因此，$\mathbf{A} = \mathbf{B}^{-1}$，或者用紧凑的记号表示为：

$$
\beta\left(\frac{\partial \mu_r}{\partial N_b}\right)_{T,V,N'} \left\langle \delta N_b\,\delta N_s \right\rangle_{T,V,\mu'} = \delta_{rs},
\tag{C.1.9}
$$

其中$\delta_{rs}$是克罗内克 $\delta$函数。换言之，一旦我们计算出了粒子数涨落的交叉关联矩阵，我们就知道了所有物种化学势随组成的变化规律。需要注意的是，$\mathbf{A}$矩阵对应于$T$、$V$和$N'$固定的系综，而$\mathbf{B}$则表达了恒定$\mu$、$V$、$T$系综中的涨落。

## 模拟中的 Kirkwood-Buff 方法

为了计算混合物的相行为，我们需要知道化学势随组成的变化。对于稠密液体，我们不能使用粒子插入方法（第 8.5.1 节）或巨正则模拟（第 6.3 节）来计算$\mu_i$，因为成功插入的概率变得非常小。正是在这样的条件下，式 (C.1.2) 变得有用：如果不同物种的粒子数是固定的，$S_{ab}(\mathbf{q}=0)$恒等于零，但我们可以计算$q \to 0$时的$S_{ab}(\mathbf{q}=0)$极限，该极限定义良好，且除可能的有限尺寸效应外，等于所需的值。

---

[^1]: 文献[[719]](references.md#ref-719)描述了在各物种分子总数固定的系统中利用 Kirkwood-Buff 关系的另一种方法。