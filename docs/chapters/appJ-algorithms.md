# 算法

本附录描述了正文中使用的几个通用算法。

## 高斯分布

**算法 36　高斯分布**

```
如果 R1 和 R2 是在区间 {0, 1} 上均匀分布的两个随机数，
则由下式给出的 Xg
  Xg = Xavg + sigma * sqrt(-ln(R1)) * cos(2*pi*R2)
服从均值为 Xavg、方差为 sigma^2 的高斯分布。
```

**注释** 上述算法只是其中一个示例。它简单，但不一定是最快的。

## 试探取向的选择

**算法 37　试探取向的选择**

```
在构型偏倚Monte Carlo方法中，我们经常需要从一组 k 个
试探方向中选择下一个键方向。下面，我们假设各个
试探方向的（玻尔兹曼）权重 w(n) 已知。

function select(w,sumw)

  ws=R*sumw
  cumw=w(1)
  n=1
  while cumw < ws do
    n=n+1
    cumw=cumw+w(n)
  enddo
end function
```

**具体注释**

1. 对于较大的 $k$ 值，二分法[[38]](references.md#ref-38) 可能更高效。

## 在球面上生成随机向量

**算法 38　单位球面上的随机向量**

```
function ranor

  ransq=2.
  do while (ransq.ge.1.0)
    ran1=1.-2.*R      % 继续直到向量在单位球内
    ran2=1.-2.*R
    ran3=1.-2.*R
    ransq=ran1*ran1+ran2*ran2+ran3*ran3
  enddo
  or = 1.0/sqrt(ransq)
  bx=ran1*or
  by=ran2*or
  bz=ran3*or
end function
```

**具体注释**

1. 上述算法只是其中一个示例。它简单，但不一定是最快的。

## 生成键长

**算法 39　使用谐振弹簧生成键长（三维）**

```
function bondl

  alpha = kv/(kB*T)
  lM = (l0/2)*(1+sqrt(1+8/(alpha*l0^2)))

  ready=.false.
  while ready == .false. do
    l = gauss(alpha,lM)         % 以高斯分布生成 l
    aux = 2*(-(l/lM-1)+ln(l/lM))  % 辅助量
    if R <= exp(aux) then       % 是否接受？
      ready=.true.
    endif
  enddo                       % 拒绝步骤
end function
```

**具体注释**

1. 键长具有以下分布：
   $$
   p(l) \propto \exp[-\beta \cdot 0.5 k_v (l - l_0)^2] dl \propto l^2 \exp[-\beta \cdot 0.5 k_v (l - l_0)^2] dl
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

```
function bonda(xn,i)
  ready=.false.
  while ready == .false. do
    b = ranor                  % 以键弯曲势给定的玻尔兹曼概率

    dx1x2=xn(i-1)-xn(i-2)     % 单位球上的单位向量
    u12=dx1x2/|dx1x2|          % 向量 r21 = r_{i-1} - r_{i-2}
    theta=acos(b*dx1x2)        % 归一化向量
    bu = ubb(theta)            % 弯曲角 theta = arccos(u12_hat . b_hat)

    if R < exp(-beta*bu) then  % 拒绝检验
      ready=.true.
    endif
  enddo
end function
```

**具体注释**

1. 本算法使用朴素的拒绝方案来生成取向 $\hat{b}$ 的玻尔兹曼分布。函数 `ranor` 在单位球上生成随机向量（算法 38）。函数 `ubb`（未具体说明）给出给定角度的键弯曲能量。

## 生成键角和扭转角

**算法 41　生成键角和扭转角**

```
function tors_bonda(xn,i)

  ready=.false.
  while ready == .false. do
    b = ranor                  % 生成随机单位向量 b_hat
    dx1x2=xn(i-1)-xn(i-2)     % 向量 r21 = r_{i-1} - r_{i-2}
    dx1x2=dx1x2/|dx1x2|       % 归一化 r12: u12_hat = r12/|r12|
    dx2x3=xn(i-2)-xn(i-3)     % 向量 r23 = r_{i-2} - r_{i-3}
    dx2x3=dx2x3/|dx2x3|       % 归一化 r23: u23_hat = r23/|r23|
    theta=acos(b*dx1x2)        % 弯曲角 theta = arccos(u12_hat . b_hat)
    ubb=ubb(theta)             % 键弯曲能量
    xx1=b x dx1x2              % 叉积: xx1 = b_hat x u12_hat
    xx2=dx1x2 x dx2x3          % 叉积: xx2 = u12_hat x u23_hat
    [... 归一化 xx1 和 xx2 ...]
    phi=acos(xx1*xx2)          % 扭转角 phi = arccos(xx1_hat . xx2_hat)
    utors=utors(phi)           % 确定扭转能量
    usum=ubb+utors             % 拒绝检验

    if R < exp(-beta*usum) then
      ready=.true.
    endif
  enddo
end function
```

**具体注释**

1. 本算法使用朴素的拒绝方案来生成取向 $\hat{b}$ 的玻尔兹曼分布。
1. 在文献中，扭转角有不同的定义。
1. 函数 `ranor` 在单位球上均匀生成向量（算法 38），函数 `ubb`（未具体说明）给出给定角度 $\theta$ 的键弯曲能量，函数 `utors`（也未具体说明）给出二面角 $\phi$ 的扭转能量。