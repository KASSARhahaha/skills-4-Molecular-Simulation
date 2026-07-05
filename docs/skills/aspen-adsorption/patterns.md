# Patterns — Aspen Adsorption 在气体吸附过程模拟方面的应用

## Pattern: 等温线选型决策
**When to use**: 拟合实验吸附数据选 Aspen Adsorption 等温线方程时。
**How**:
1. 单层均匀表面、低中压 → Langmuir（IP1/IP2）。
2. 多孔分子筛、非均匀表面、需温度修正 → Langmuir-Freundlich/Toth（IP1–IP6）。
3. 中压经验拟合、无物理外推需求 → Freundlich。
4. 多层吸附、比表面积测定场景 → BET。
5. 极低压线性区 → Henry。
6. 内置全部不匹配 → User Procedure 自定义。
**Trade-offs**: 越复杂的方程拟合越好但参数更难独立验证；优先用最简且物理意义明确的 Langmuir/LDF 组合。

## Pattern: LDF 优先传质模型
**When to use**: 99% 的 Aspen Adsorption 动态模拟。
**How**: 选 LDF（线性推动力），∂qi/∂t = CMT_i·(q*i−qi)；从实验或经验式得到 CMT_i。
**Trade-offs**: LDF 是形式最简且有合理物理基础的模型；其他模型（如固体扩散）精度相当但参数更难取，仅特殊体系使用。

## Pattern: 穿透曲线 → 模型验证闭环
**When to use**: 引入新材料、新体系、新床型时。
**How**:
1. 用静态实验数据回归等温线 IP。
2. 在 Aspen Adsorption 建立固定床 + 等温线 + LDF + Ergun。
3. 运行动态求解，画出口浓度-时间曲线。
4. 与实验穿透曲线对齐。
5. 一致 → 进入 PSA/TSA 设计；不一致 → 反推机理偏差（孔道、热效应、动态 vs 静态吸附能力、CMT 来源）。
**Trade-offs**: 单点验证不够；至少 2 个压力/温度交叉验证。

## Pattern: 操作参数单变量敏感性扫描
**When to use**: 床层设计或操作条件未定时。
**How**: 在 Aspen Adsorption 中对 5 个变量逐项扫描——
- 操作压力 P
- 吸附剂装量 m
- 传质系数 CMT
- 床层温度 T
- 进料速率 F
记录穿透时间与陡度变化，找主导变量。
**Trade-offs**: 单变量扫描不揭示耦合；找完主导变量后可做 2 变量网格扫描。

## Pattern: 动态等温线纠偏
**When to use**: 模拟穿透时间系统性偏长于实验。
**How**:
1. 怀疑静态等温线高估动态吸附能力。
2. 用动态吸附实验重新拟合等温参数（Zhang 等[21] 改进方案）。
3. 重新跑穿透曲线验证。
**Trade-offs**: 动态实验成本高于静态；只在偏差超出容忍域时启动。

## Pattern: PSA 循环时序设计
**When to use**: 设计普通 PSA、VPSA、DR-PSA。
**How**:
1. 从 4 步循环母版开始：升压 → 吸附 → 降压 → 吹扫/置换。
2. 加入压力均衡（PE）步节能。
3. 床数选择：1 床机理验证 → 4 床工程化 → 6–12 床规模化。
4. VPSA 在降压/吹扫步骤加抽真空；扫真空度。
5. DR-PSA 在吸附基础上加塔顶+塔底回流。
**Trade-offs**: 多床与多步可提升纯度与回收率，但时序复杂度非线性上升，CSS 收敛时间变长。

## Pattern: TSA / VTSA 选型
**When to use**: 溶剂回收、CO₂ 燃烧尾气捕集、有利废热场景。
**How**:
1. 默认 TSA：低温吸附、高温脱附（热惰性气体或蒸汽吹扫）。
2. TSA 升温慢 → 加抽真空形成 VTSA。
3. 无废热源、需快速升温 → ESA（电加热）。
**Trade-offs**: TSA 利废热节能但循环慢；ESA 快但能耗高；VTSA 在两者间平衡。

## Pattern: 多软件联用分工
**When to use**: 单软件出现盲区——Aspen Adsorption 不擅长能量/有效能/换热几何。
**How**:
- 物性 → Aspen Hysys 预测/评估
- 动力学与循环 → Aspen Adsorption
- 全流程能耗/有效能 → Aspen Plus
- 床层/换热几何 → Aspen EDR
- Property Method 必须保持一致（如都选 Peng-Robinson）
**Trade-offs**: 联用复杂度上升，需保证联接面物性一致；典型失真来源是 Property Method 不匹配。

## Pattern: 偏差诊断三件套
**When to use**: 模拟结果与实验系统性偏差，且参数已校准。
**How**: 按顺序排查——
1. 等温参数 IP 是否来自实验回归（不是默认值）。
2. CMT 来源是否独立（不是硬凑）。
3. 床层均匀性、单一孔径、热效应、径向扩散——这四项是 Aspen Adsorption 默认假设的失真源。
**Trade-offs**: 每补一项复杂度都上升；先验证再加复杂模型。
