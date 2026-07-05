# Glossary — Adsorption Simulation Notes

**Adsorbate** — gas molecule retained on a solid surface (Note 1).
**Adsorbent** — porous solid that retains adsorbate (Note 1).
**Apparent (particle) density ρ_p** — mass over particle envelope volume, includes internal pores (Note 3).
**Axial dispersion coefficient D_L** — effective back-mixing along the bed axis (Note 3).
**Bed voidage (total) φ** — all voids over total bed volume (Note 3).
**Blowdown** — depressurization step in PSA that releases the extract (Note 1).
**Breakthrough time t_bt** — feed time after which raffinate impurity exceeds a threshold (Note 6).
**Bulk density ρ_b** — particles poured in a beaker, divided by beaker volume (Note 3).
**Burke–Plummer equation** — turbulent-limit pressure-drop relation (Re > ~2000) (Note 3).
**Carman–Kozeny equation** — laminar-limit pressure-drop relation (Note 3).
**CSS_Info** — Aspen metadata struct required by gCSS: step time, step-time node count (Note 6, 7).
**Cyclic steady state (CSS)** — periodic solution where bed loading at cycle start = cycle end (Note 6).
**Cycle Organizer** — Aspen block that schedules valve state changes per step (Note 4).
**Ergun equation** — combined laminar+turbulent pressure-drop correlation for packed beds (Note 3).
**Extract (Heavy Product, HP)** — gas released during regeneration, enriched in strongly-adsorbed species (Note 1).
**Flowsheet Constraints** — Aspen domain scripting language; RealVariable declarations and assignments (Note 7).
**Freundlich isotherm** — empirical power-law `q = K·p^(1/n)`; no finite saturation (Note 2).
**gas_bed block** — Aspen 1D plug-flow column with axial dispersion, Ergun ΔP, configurable kinetics (Note 4).
**gas_interaction block** — Aspen tank/mixer; used as buffer in single-bed PSA variants (Note 4).
**gCSS block** — gas Cyclic Steady State block; discretizes time+space, finds CSS directly (Note 6).
**Henry's law (linear isotherm)** — `q = H·c`; low-loading limit of every isotherm (Note 2).
**Interstitial (external) voidage ε_e** — volume outside particles over total bed volume (Note 3).
**Intra-particle (internal) voidage ε_i** — pore volume inside particles; two conventions exist (Note 3).
**Isotherm** — equilibrium loading q* as function of partial pressure at fixed T (Note 1, 2).
**Langmuir isotherm** — `q = q_max · b·p/(1+b·p)`; monolayer on homogeneous sites (Note 2).
**LDF (Linear Driving Force)** — `dq_i/dt = k_i·(q_i* − q_i)`; default lumped kinetic model (Note 5).
**Loading q** — adsorbate mass (or mol) per adsorbent mass (or per bed volume) (Note 1).
**Mass-transfer coefficient k_i** — fitted constant in LDF; species- and adsorbent-specific (Note 5).
**Peng–Robinson EOS** — cubic EOS for moderately nonideal gas; default for H₂/CO/CO₂ streams (Note 3).
**Pressure equalization** — PSA step connecting high-P and low-P beds to reuse compression work (Note 1, 6).
**PSA (Pressure Swing Adsorption)** — cycles pressure to drive adsorption/desorption (Note 1).
**Purge** — low-pressure step using LP to sweep desorbed species out of bed (Note 1).
**Raffinate (Light Product, LP)** — gas leaving during adsorption, depleted in strongly-adsorbed species (Note 1).
**Redlich–Kwong–Soave (RK-Soave) EOS** — cubic EOS for hydrocarbons up to ~10 bar (Note 3).
**Selectivity** — ratio of loadings or isotherm slopes for two species (Note 1).
**Shape factor ψ** — particle sphericity, 1 = sphere (Note 3).
**Sips isotherm** — Langmuir–Freundlich hybrid; Freundlich at low p, Langmuir saturation at high p (Note 2).
**Skeletal (solid) density ρ_s** — mass over volume of pore-free solid (Note 3).
**Skarstrom cycle** — canonical 4-step two-bed PSA: pressurize → adsorb → blowdown → purge (Note 1).
**Superficial velocity u** — gas volumetric flow over total bed cross-section (Note 3).
**Toth isotherm** — empirical, gas-phase, finite saturation; common for hydrocarbon/zeolite systems (Note 2).
**TSA (Temperature Swing Adsorption)** — cycles temperature to drive adsorption/desorption (Note 1, 7).
**Working capacity Δq** — difference between max and min loading reached in one cycle (Note 1).
