# Glossary — Simulation of Adsorption Processes

**Adsorbate** — species retained on the solid surface (Ch 1, Ch 2).
**Adsorbent** — solid phase that retains adsorbate (e.g., zeolite 13X, activated carbon) (Ch 1, Ch 5).
**Adsorption** — gas/liquid molecules binding to a solid surface (Ch 1).
**Adsorption isotherm** — q* = f(C) at constant T; equilibrium relation (Ch 2).
**Bed voidage (ε)** — inter-particle void volume / total bed volume (Ch 3).
**BET isotherm** — multilayer extension of Langmuir (Ch 2).
**Bi-Langmuir (Dual-site Langmuir)** — sum of two Langmuir terms for heterogeneous sites (Ch 2).
**Blowdown** — PSA step that depressurizes the bed to release adsorbate (Ch 5, Ch 9).
**Breakthrough curve** — C_product(t) showing the front arrival at the column exit (Ch 8).
**Bulk density** — mass of particles / total bed volume (Ch 3).
**Convection-only balance** — column PDE with axial dispersion term dropped (Ch 4, Ch 7).
**CSS (Cyclic Steady State)** — bed profile periodicity across cycles (Ch 9).
**Cycle Organizer** — Aspen's state machine for PSA step sequencing (Ch 6, Ch 9).
**Cyclic steady state** — see CSS.
**Countercurrent (cc) blowdown / purge** — regeneration steps flowing opposite to feed direction; standard in Skarstrom cycle (Ch 5).
**Density hierarchy** — bulk < particle < skeletal (Ch 3).
**Discretization** — converting the column PDE to ODEs over spatial nodes (Ch 7).
**Dual-site Langmuir** — see Bi-Langmuir (Ch 2).
**Equilibrium loading (q*)** — solid loading at equilibrium with local fluid (Ch 2, Ch 7).
**Equilibrium step (PE)** — see Pressure equalization.
**Ergun equation** — packed-bed ΔP correlation spanning laminar + turbulent (Ch 4).
**Extract / Heavy Product (HP)** — more-strongly-adsorbed species released during desorption step (Ch 5, Ch 8).
**Extended Langmuir** — multicomponent Langmuir with Σ_j b_j c_j denominator (Ch 2).
**Film model** — gas-film vs solid-film LDF driving force choice (Ch 7).
**Freundlich isotherm** — empirical power-law q = h c^(1/n) (Ch 2).
**gas_bed** — Aspen's PDE-based adsorption column unit (Ch 6, Ch 8).
**gas_interaction** — Aspen's simplified gas adsorber unit (Ch 6, Ch 10).
**Henry's law (linear) isotherm** — q = Hc; low-C limit (Ch 2).
**Kinetic model** — pore-diffusion / LDF / Fickian choice in Aspen (Ch 7).
**LDF (Linear Driving Force)** — ∂q/∂t = k(q* − q); lumped mass-transfer model (Ch 7).
**LiLSX (Li-exchanged low-Si X zeolite)** — modern high-capacity adsorbent for O2 PSA (Ch 5).
**Linear isotherm** — see Henry's law.
**Loading (q)** — adsorbed-phase concentration (mol/m³ adsorbent) (Ch 2).
**Mass transfer coefficient (k_i)** — LDF rate constant, 1/s (Ch 7).
**Mass transfer zone (MTZ)** — axial region where loading changes from feed to clean (Ch 8).
**Momentum balance** — pressure-velocity coupling; Ergun is the default (Ch 4).
**Multi-bed PSA** — staggered bed trains for continuous product (Ch 5).
**Naphtha upgrading** — adsorptive removal of trace species from naphtha (Ch 10).
**Nodes** — spatial discretization points in the column (Ch 7).
**Particle density** — mass / particle envelope volume (incl. pores) (Ch 3).
**Peng-Robinson** — default EOS for gas-phase property method (Ch 6).
**P/F ratio (purge-to-feed)** — purge gas volume / feed volume; key purity/recovery lever (Ch 9).
**Pressure equalization (PE)** — PSA step sharing pressure between beds to save compression energy (Ch 5, Ch 9).
**Raffinate / Light Product (LP)** — less-strongly-adsorbed species collected at column top during adsorption step (Ch 5, Ch 8).
**PSA (Pressure Swing Adsorption)** — regeneration via pressure change (Ch 1).
**PTSA** — hybrid pressure/temperature swing (Ch 1).
**Purge** — low-pressure clean-up step using product or exhaust gas (Ch 9).
**q*** — see Equilibrium loading.
**Recovery** — fraction of feed component recovered in product (Ch 9).
**Regeneration** — desorption step that releases adsorbate (Ch 1).
**Skeletal (solid) density** — mass / solid-only volume (Ch 3).
**Sips (Langmuir–Freundlich)** — hybrid isotherm with monolayer ceiling (Ch 2).
**Skarstrom cycle** — Charles Skarstrom's 1960s 4-step PSA: pressurize / adsorb / countercurrent blowdown / countercurrent purge; optional 5th pressure-equalization step (Ch 5).
**Solid film** — LDF driving force written in terms of q* − q (Ch 7).
**Superficial velocity (u)** — volumetric flow / empty-column cross-section (Ch 4).
**Toth isotherm** — empirical gas-phase isotherm with monolayer ceiling (Ch 2).
**TSA (Temperature Swing Adsorption)** — regeneration via temperature change (Ch 1).
**UDS1 (Upwind Differencing Scheme 1)** — default first-order spatial discretization (Ch 7).
**TSA thermal coupling** — heating/cooling path slope on q-P diagram depends on C_ps, C_pg, m_s, m_g (Ch 1).
**Valve schedule** — 2-bed PSA valve table using codes 0 (closed), 1 (open), 2 (flow-controlled), 3 (pressure-controlled) (Ch 5).
**Vacuum PSA (VPSA)** — PSA with vacuum-assisted regeneration (Ch 1).
**Void fraction** — see Bed voidage.
**Working capacity (Δq)** — q at feed P_high minus q at regeneration P_low; the useful swing per cycle, expressed per mass of adsorbent (Ch 1, Ch 5, Ch 9).
**Zeolite 5A** — original Skarstrom-cycle adsorbent for O2 PSA; also molecular-sieve for n-paraffins (Ch 5).
**Zeolite 13X** — equilibrium-selective adsorbent for N2, CO2 (Ch 5, Ch 8).

## Bilingual EN ↔ CN Cross-Reference (from DOCX source)
- Adsorption → 吸附 (xī fù)
- Adsorbent → 吸附剂
- Adsorbate → 吸附质
- Desorption / Regeneration → 解吸 / 再生
- Pressure swing adsorption → 变压吸附
- Temperature swing adsorption → 变温吸附
- Vacuum pressure swing adsorption → 真空变压吸附
- Isotherm → 等温线
- Working capacity → 工作容量
- Packed bed → 填料床
- Void fraction → 空隙率
- Pressure equalization → 均压
- Raffinate / light product → 顶部产品 / 轻组分产品
- Extract / heavy product → 提取物 / 重组分产品
- Skeletal density → 骨骼密度
- Particle density → 颗粒密度
- Linear driving force → 线性推动力
- Mass transfer coefficient → 传质系数
- Cycle organizer → 循环组织器
- Breakthrough curve → 穿透曲线
