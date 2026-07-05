# Glossary — Aspen Adsorption 在气体吸附过程模拟方面的应用

**Aspen Adsorption** — AspenTech 公司开发的气液相吸附脱附过程模拟软件，原名 Aspen Adsim；内置稳态与动态模拟、Cycle Organizer 循环控制器、User Procedure 自定义接口 (Ch 1)。

**Aspen Plus / Hysys / EDR** — 同公司流程模拟软件，可与 Aspen Adsorption 物性同源联用，做能量/有效能/换热几何分析 (Ch 4)。

**BET 方程** — 多层吸附等温方程，常用于比表面积测定；Aspen Adsorption 内置可选 (Ch 1)。

**CSS（Cyclic Steady State，循环稳态）** — PSA/TSA 循环多次运行后床层温度/浓度分布达到周期性重复的状态；循环优化的目标终点 (Ch 3)。

**CMT（传质系数）** — 单位 s⁻¹，线性推动力模型中的常数；LDF 方程 ∂qi/∂t = CMT·(q*i−qi) (Ch 1)。

**Cycle Organizer** — Aspen Adsorption 循环控制模块，按时序参数自动切换流向/流速/阀门，实现 PSA/TSA 周期操作 (Ch 1, Ch 3)。

**DR-PSA（Dual Reflux PSA，双回流变压吸附）** — 在吸附基础上加塔顶轻组分回流与塔底重组分回流，同时富集轻、重组分 (Ch 3)。

**EDR（Exchanger Design Rating）** — Aspen 换热器与床层几何评估软件；用于新型吸附床结构对比 (Ch 4)。

**Ergun 方程** — 默认动量守恒/压降方程；含黏性项与惯性项，需 ψ、Rp、εi、μ、vg、ρg (Ch 1)。

**ESA（Electric Swing Adsorption，变电吸附）** — 通过电加热带/焦耳热直接加热床层实现脱附；升温快但能耗高于 TSA (Ch 3)。

**Freundlich 方程** — 经验型吸附等温方程，中压段适用；高压时不具物理意义 (Ch 1)。

**Henry 方程** — 极低压线性吸附等温近似 (Ch 1)。

**IP1–IP6** — Aspen Adsorption 中吸附质 i 的 1–6 个等温线参数；Langmuir 用 IP1/IP2，Langmuir-Freundlich 用 IP1–IP6 含温度修正 (Ch 1)。

**LDF（Linear Driving Force，线性推动力）** — ∂qi/∂t = CMT·(q*i−qi)；Aspen Adsorption 默认且最常用的传质速率模型，简洁且有物理基础 (Ch 1)。

**Langmuir 方程** — 单层均匀表面吸附等温方程；低压至中压常用；只用到 IP1、IP2 (Ch 1)。

**Langmuir-Freundlich 方程** — 含温度修正项的扩展 Langmuir；适用于非均匀表面、多孔分子筛；用 IP1–IP6 (Ch 1)。

**MOF（Metal-Organic Framework，金属有机框架）** — 高比表面积多孔吸附材料 (Ch 1)。

**PSA（Pressure Swing Adsorption，变压吸附）** — 高压吸附、低压脱附的循环分离过程；4 基本步骤为升压-吸附-降压-吹扫置换 (Ch 3)。

**TSA（Temperature Swing Adsorption，变温吸附）** — 低温吸附、高温脱附；适合溶剂回收，可利用低等级废热 (Ch 3)。

**User Procedure** — Aspen Adsorption 自定义模型接口，可挂载自定义等温方程或传递模型 (Ch 1)。

**VTSA（Vacuum + TSA，真空变温吸附）** — 抽真空与适度加热联合脱附；减少 TSA 升温滞后 (Ch 3)。

**VPSA（Vacuum Pressure Swing Adsorption，真空变压吸附）** — 在 PSA 降压/吹扫步骤抽真空的变型；真空度与压力均衡方式是关键变量 (Ch 3)。

**穿透曲线（Breakthrough Curve）** — 出口吸附质浓度-时间曲线；越陡吸附速率越大；Aspen Adsorption 内置作图工具 (Ch 2)。

**穿透时间（Breakthrough Time）** — 出口浓度达到阈值（典型 5%–10% 入口浓度）的时刻，即床层有效工作时间 (Ch 2)。

**动态吸附能力** — 动态流动条件下的实际吸附量，通常略小于静态吸附能力；直接套静态等温线会引起系统性偏差 (Ch 2)。

**压力均衡（Pressure Equalization, PE）** — PSA 多床循环中用降压床的余压给另一床升压，节能关键步骤 (Ch 3)。

**床层分层（Layered Bed）** — 一床内填充多种吸附剂（如活性炭+分子筛），针对不同组分；双层吸附剂 PSA 是常见配置 (Ch 3)。
