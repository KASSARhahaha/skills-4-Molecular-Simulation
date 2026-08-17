# 算法

本附录描述了正文中使用的几个通用算法。

## 高斯分布

**算法 36　高斯分布**

<table class="algbox" markdown="1">
<tbody markdown="1">
<tr markdown="1"><td class="algcode" markdown="span"><code>若&nbsp;</code>$R_1$<code>&nbsp;与&nbsp;</code>$R_2$<code>&nbsp;是区间&nbsp;</code>$\{0,1\}$<code>&nbsp;上均匀分布的两个随机数，则由</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span">$X_g = X_{\mathrm{avg}} + \sigma\sqrt{-\ln(R_1)}\,\cos(2\pi R_2)$</td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>给出的&nbsp;</code>$X_g$<code>&nbsp;服从均值为&nbsp;</code>$X_{\mathrm{avg}}$<code>、方差为&nbsp;</code>$\sigma^2$<code>&nbsp;的高斯分布。</code></td><td class="algcom" markdown="span"></td></tr>
</tbody>
</table>

**注释** 上述算法只是其中一个示例。它简单，但不一定是最快的。

## 试探取向的选择

**算法 37　试探取向的选择**

在构型偏倚 Monte Carlo 方法中，我们经常需要从一组 $k$ 个试探方向中选择下一个键方向。下面，我们假设各个试探方向的（玻尔兹曼）权重 $w(n)$ 已知。

<table class="algbox" markdown="1">
<tbody markdown="1">
<tr markdown="1"><td class="algcode" markdown="span"><code>function&nbsp;select(w,sumw)</code></td><td class="algcom" markdown="span">以概率 $p(n) = w(n)/\sum_j w(j)$ 选择试探取向 $n$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>ws=R*sumw</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>cumw=w(1)</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>n=1</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>while&nbsp;cumw&nbsp;&lt;&nbsp;ws&nbsp;do</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>n=n+1</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>cumw=cumw+w(n)</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>enddo</code></td><td class="algcom" markdown="span">函数返回 $n$，</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>end&nbsp;function</code></td><td class="algcom" markdown="span">即所选试探位置的索引</td></tr>
</tbody>
</table>

**具体说明**（一般说明见第 1 章「算法」）：

1. 对于较大的 $k$ 值，二分法[[38]](references.md#ref-38) 可能更高效。

## 在球面上生成随机向量

**算法 38　单位球面上的随机向量**

<table class="algbox" markdown="1">
<tbody markdown="1">
<tr markdown="1"><td class="algcode" markdown="span"><code>function&nbsp;ranor</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"></td><td class="algcom" markdown="span">生成一个三维随机单位向量，分量为 bx, by, bz</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>ransq=2.</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>do&nbsp;while&nbsp;(ransq.ge.1.0)</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>ran1=1.-2.*R</code></td><td class="algcom" markdown="span">继续直到向量在单位球内</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>ran2=1.-2.*R</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>ran3=1.-2.*R</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>ransq=ran1*ran1+ran2*ran2+ran3*ran3</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>enddo</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>or&nbsp;=&nbsp;1.0/sqrt(ransq)</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>bx=ran1*or</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>by=ran2*or</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>bz=ran3*or</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>end&nbsp;function</code></td><td class="algcom" markdown="span"></td></tr>
</tbody>
</table>

**具体说明**（一般说明见第 1 章「算法」）：

1. 上述算法只是其中一个示例。它简单，但不一定是最快的。

## 生成键长

**算法 39　使用谐振弹簧生成键长（三维）**

<table class="algbox" markdown="1">
<tbody markdown="1">
<tr markdown="1"><td class="algcode" markdown="span"><code>function&nbsp;bondl</code></td><td class="algcom" markdown="span">返回键长 $\ell$；假设为谐振弹簧，弹簧常数 $k_v$。$\ell_0$：$T=0$ 时的键长</td></tr>
<tr markdown="1"><td class="algcode" markdown="span">$\alpha = k_v/(k_B T)$</td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span">$\ell_M = (\ell_0/2) * (1 + \sqrt{1 + 8/(\alpha \ell_0^2)})$</td><td class="algcom" markdown="span">$T$ 时的极大值位置</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>ready=.false.</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>while&nbsp;ready&nbsp;==&nbsp;.false.&nbsp;do</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code>$\ell$<code>&nbsp;=&nbsp;gauss(</code>$\alpha$<code>,</code>$\ell_M$<code>)</code></td><td class="algcom" markdown="span">以高斯分布生成 $\ell$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>aux&nbsp;</code>$= 2 * [-(\ell/\ell_M - 1) + \ln(\ell/\ell_M)]$</td><td class="algcom" markdown="span">辅助量</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>if&nbsp;R&nbsp;</code>$\leq$<code>&nbsp;exp([aux])&nbsp;then</code></td><td class="algcom" markdown="span">是否接受？</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</code><code>ready=.true.</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>endif</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>enddo</code></td><td class="algcom" markdown="span">拒绝步骤</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>end&nbsp;function</code></td><td class="algcom" markdown="span"></td></tr>
</tbody>
</table>

**具体说明**（一般说明见第 1 章「算法」）：

1. 键长具有以下分布：
   $$
   p(l) \propto \exp[-\beta \cdot 0.5 k_v (l - l_0)^2] \mathrm{d}l \propto l^2 \exp[-\beta \cdot 0.5 k_v (l - l_0)^2] \mathrm{d}l
   $$
1. 我们利用了 $x - 1 \geq \ln x$ 这一事实。
1. $\text{gauss}(\alpha, l_M)$ 是一维正态分布，见第 J.1 节。

三维线性谐振分子的键长分布接近高斯分布，可以从 Box-Muller 算法[[66]](references.md#ref-66) 出发生成。考虑一个可以沿其轴向振动的线性分子，振动力常数记为 $\kappa$。对于固定取向，平衡键长记为 $l_0$，逆温度记为 $\beta$。我们忽略旋转与振动之间的耦合。

在这些条件下，分子的长度分布为

$$
P(l) \sim l^2 \exp[-0.5\beta\kappa(l - l_0)^2] = \exp[-0.5\beta\kappa(l - l_0)^2 + 2\ln l] \tag{J.4.1}
$$

下文中，归一化并不重要。我们无法直接对这个分布进行采样。然而，我们可以使用拒绝法。首先，我们确定 $\ln P(l)$ 的极大值位置 $l_M$。于是我们得到：

$$
\beta\kappa(l_M - l_0) - 2/l_M = 0 \tag{J.4.2}
$$

将 $\beta\kappa$ 记为 $\alpha$，我们得到

$$
l_M^2 - l_M l_0 - 2/\alpha = 0 \tag{J.4.3}
$$

即

$$
l_M = 0.5(l_0 + \sqrt{l_0^2 + 8/\alpha}) \tag{J.4.4}
$$

我们可以用在 $l_0$ 附近的高斯分布来近似公式 (J.4.1)：

$$
P'(l) \sim \exp[-0.5\alpha(l - l_M)^2] \tag{J.4.5}
$$

$P'(l)$ 与 $P(l)$ 的关系为

$$
P(l)/P'(l) = \exp\left[2\ln(l/l_M) - \frac{2(l - l_M)}{l_M}\right] \leq 1 \quad \text{。}
$$

因此我们可以从 $P'(l)$ 中抽取 $l$ 的值，然后当

$$
R > \exp\left[2\ln(l/l_M) - \frac{2(l - l_M)}{l_M}\right]
$$

时拒绝该 $l$ 值。

## 生成键角

**算法 40　生成键角**

<table class="algbox" markdown="1">
<tbody markdown="1">
<tr markdown="1"><td class="algcode" markdown="span"><code>function&nbsp;bonda(xn,i)</code></td><td class="algcom" markdown="span">以键弯曲势给定的玻尔兹曼概率生成键取向向量 $\hat{\mathbf b}$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>ready=.false.</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>while&nbsp;ready&nbsp;==&nbsp;.false.&nbsp;do</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>b&nbsp;=&nbsp;ranor</code></td><td class="algcom" markdown="span">单位球上的单位向量</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>dx1x2=xn(i-1)-xn(i-2)</code></td><td class="algcom" markdown="span">向量 $\mathbf r_{21} = \mathbf r_{i-1} - \mathbf r_{i-2}$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>u12=dx1x2/|dx1x2|</code></td><td class="algcom" markdown="span">归一化向量</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>theta=acos(b*dx1x2)</code></td><td class="algcom" markdown="span">弯曲角 $\theta = \arccos(\hat{\mathbf u}_{12} \cdot \hat{\mathbf b})$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>bu&nbsp;=&nbsp;ubb(theta)</code></td><td class="algcom" markdown="span">键弯曲能量</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>if&nbsp;R&nbsp;&lt;&nbsp;exp(-beta*bu)&nbsp;then</code></td><td class="algcom" markdown="span">拒绝检验</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</code><code>ready=.true.</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>endif</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>enddo</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>end&nbsp;function</code></td><td class="algcom" markdown="span"></td></tr>
</tbody>
</table>

**具体说明**（一般说明见第 1 章「算法」）：

1. 本算法使用朴素的拒绝方案来生成取向 $\hat{\mathbf b}$ 的玻尔兹曼分布。函数 `ranor` 在单位球上生成随机向量（算法 38）。函数 `ubb`（未具体说明）给出给定角度的键弯曲能量。
1. 译注：原书此框末尾连写两个 `endif`，第二个用来收 `while` 循环，与全书体例（`while $...$ do` 以 `enddo` 收，参见算法 41）不合，应为 `enddo`；本书已改正。

## 生成键角和扭转角

**算法 41　生成键角和扭转角**

<table class="algbox" markdown="1">
<tbody markdown="1">
<tr markdown="1"><td class="algcode" markdown="span"><code>function&nbsp;tors_bonda(xn,i)</code></td><td class="algcom" markdown="span">生成一个具有由扭转势和键弯曲势决定的取向玻尔兹曼分布的单位向量</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>ready=.false.</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>while&nbsp;ready&nbsp;==&nbsp;.false.&nbsp;do</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>b&nbsp;=&nbsp;ranor</code></td><td class="algcom" markdown="span">生成随机单位向量 $\hat{\mathbf b}$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>dx1x2=xn(i-1)-xn(i-2)</code></td><td class="algcom" markdown="span">向量 $\mathbf r_{21} = \mathbf r_{i-1} - \mathbf r_{i-2}$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>dx1x2=dx1x2/|dx1x2|</code></td><td class="algcom" markdown="span">归一化 $\mathbf r_{12}$：$\hat{\mathbf u}_{12} \equiv \mathbf r_{12}/|\mathbf r_{12}|$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>dx2x3=xn(i-2)-xn(i-3)</code></td><td class="algcom" markdown="span">向量 $\mathbf r_{23} = \mathbf r_{i-2} - \mathbf r_{i-3}$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>dx2x3=dx2x3/|dx2x3|</code></td><td class="algcom" markdown="span">归一化 $\mathbf r_{23}$：$\hat{\mathbf u}_{23} \equiv \mathbf r_{23}/|\mathbf r_{23}|$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>theta=acos(b&nbsp;*&nbsp;dx1x2)</code></td><td class="algcom" markdown="span">弯曲角 $\theta = \arccos(\hat{\mathbf u}_{12} \cdot \hat{\mathbf b})$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>ubb=ubb(theta)</code></td><td class="algcom" markdown="span">键弯曲能量</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>xx1=b&nbsp;</code>$\times$<code>&nbsp;dx1x2</code></td><td class="algcom" markdown="span">叉积：$\mathbf{xx1} = \hat{\mathbf b} \times \hat{\mathbf u}_{12}$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>xx2=dx1x2&nbsp;</code>$\times$<code>&nbsp;dx2x3</code></td><td class="algcom" markdown="span">叉积：$\mathbf{xx2} = \hat{\mathbf u}_{12} \times \hat{\mathbf u}_{23}$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>[...&nbsp;归一化&nbsp;xx1&nbsp;和&nbsp;xx2&nbsp;...]</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>phi=acos(xx1&nbsp;*&nbsp;xx2)</code></td><td class="algcom" markdown="span">扭转角 $\phi = \arccos(\hat{\mathbf{xx1}} \cdot \hat{\mathbf{xx2}})$</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>utors=utors(phi)</code></td><td class="algcom" markdown="span">确定扭转能量</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>usum=ubb+utors</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>if&nbsp;R&nbsp;&lt;&nbsp;exp(-beta*usum)&nbsp;then</code></td><td class="algcom" markdown="span">拒绝检验</td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</code><code>ready=.true.</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>&nbsp;&nbsp;&nbsp;&nbsp;</code><code>endif</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>enddo</code></td><td class="algcom" markdown="span"></td></tr>
<tr markdown="1"><td class="algcode" markdown="span"><code>end&nbsp;function</code></td><td class="algcom" markdown="span"></td></tr>
</tbody>
</table>

**具体说明**（一般说明见第 1 章「算法」）：

1. 本算法使用朴素的拒绝方案来生成取向 $\hat{b}$ 的玻尔兹曼分布。
1. 在文献中，扭转角有不同的定义。
1. 函数 `ranor` 在单位球上均匀生成向量（算法 38），函数 `ubb`（未具体说明）给出给定角度 $\theta$ 的键弯曲能量，函数 `utors`（也未具体说明）给出二面角 $\phi$ 的扭转能量。