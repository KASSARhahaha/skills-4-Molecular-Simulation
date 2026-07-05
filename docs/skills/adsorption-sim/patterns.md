# Patterns — Simulation of Adsorption Processes

## Pattern: Build-Then-Cycle Aspen Workflow
**When to use**: every new Aspen Adsorption simulation, single-bed or multi-bed.
**How**:
1. Components + property method
2. Flowsheet topology (units + streams)
3. Per-unit specs (geometry, adsorbent, isotherm)
4. Cycle organizer (step sequence + valve states)
5. Numerical options (discretization, nodes)
6. Run breakthrough first, then assemble cycle
**Trade-offs**: skipping steps forces rework — fixing property method late invalidates isotherm parameters.

## Pattern: Single-Bed Breakthrough Validation
**When to use**: before assembling a multi-step PSA cycle.
**How**:
1. Configure gas_bed with feed at process P, T
2. Clean bed initial condition (q = 0)
3. Run for 2–3 estimated breakthrough times
4. Compare breakthrough time × flow against bed working capacity (mass balance)
5. If mismatch > 20 %, re-fit isotherm or k_i
**Trade-offs**: costs one extra run; saves 10+ cycles of debugging a non-converging PSA.

## Pattern: LDF Coefficient Estimation from Uptake Experiment
**When to use**: when literature k_i is unavailable for your adsorbent/gas pair.
**How**:
1. Batch uptake experiment, record q(t)
2. Identify half-uptake time t_0.5
3. k ≈ 0.4 / t_0.5 (Glueckauf approximation for spherical particles)
4. Use as starting value; refine against breakthrough curve
**Trade-offs**: lumped LDF cannot separate film vs pore resistances; use Fickian model if that distinction matters.

## Pattern: Skarstrom PSA Cycle Design
**When to use**: simplest two-bed PSA for bulk gas separation; default for O2 PSA, gas drying.
**How**:
1. Two beds in anti-phase; each cycle = pressurize → adsorb → countercurrent blowdown → countercurrent purge
2. Feed at high pressure (e.g., 5 bar), purge at near-atmospheric using raffinate slipstream
3. Choose step durations so mass-transfer zone stays inside bed during adsorption step
4. Set P/F ratio ≈ 0.1–0.3 as starting point; tune from there
5. Optional 5th pressure-equalization step shares gas from depressurizing bed to repressurizing bed before blowdown completes
**Trade-offs**: simple 4-step has lower recovery; adding PE saves 20–40 % compression work but complicates CSS convergence.
**Directionality matters**: blowdown and purge are **countercurrent** to feed direction — this strips the strongly-adsorbed species from the feed end where it accumulated, preventing contamination of the next adsorption step.

## Pattern: Pressure Equalization in 2-Bed PSA
**When to use**: high-recovery H2 / CO2 PSA where compression energy matters.
**How**:
1. Add a co-current depressurization step that pushes gas from high-P bed to low-P bed
2. Bed A: feed → PE-A (provide) → blowdown → purge → PE-B (receive) → pressurize
3. Bed B does the mirror image
4. Saves 20–40 % of feed compression work
**Trade-offs**: more cycle steps → harder CSS convergence; longer commissioning time.

## Pattern: TSA Indirect Heating Cycle
**When to use**: trace contaminant removal (drying, sulfur) where PSA cannot regenerate cleanly.
**How**:
1. Bed with internal coil or external jacket
2. Feed at ambient T; switching to hot purge (150–250 °C, often steam) regenerates
3. Cooling step before re-adsorption (avoid hot bed in feed step)
4. Cycle times minutes to hours (much slower than PSA)
5. Track the q-P cycle path: heating/cooling lines slope depends on (C_ps·m_s + C_pg·m_g) — not pure isotherms
**Trade-offs**: lower energy efficiency than PSA but cleaner regeneration for strongly-bound species.

## Pattern: Adsorbent Upgrade Path (Classic → Modern)
**When to use**: improving an existing PSA design's working capacity.
**How**:
1. Baseline: classic 5A or 13X zeolite
2. Upgrade to LiLSX (Li-exchanged low-Si X) for O2 PSA — typically ~2× working capacity
3. Re-fit isotherm parameters; verify breakthrough before re-converging CSS
**Trade-offs**: LiLSX has higher adsorbent cost; payback only justified for larger plants or higher-purity targets.

## Pattern: Layered Bed Design
**When to use**: feed contains both fast-adsorbing contaminants and the target component.
**How**:
1. Front layer: cheap, high-capacity adsorbent for the contaminant
2. Back layer: selective adsorbent for the target
3. Configure as Layer 1 + Layer 2 in the same gas_bed
4. Tune layer thicknesses so contaminant never reaches Layer 2
**Trade-offs**: complicates bed packing; saves cost by using less of the expensive adsorbent.

## Pattern: CSS Convergence Tuning
**When to use**: PSA cycle that drifts or fails to converge.
**How**:
1. Run 5–10 cycles, plot q(z, t) at end of each cycle
2. If profile drifts monotonically → bed is overloaded; shorten feed step or lengthen purge
3. If profile oscillates → step timing mismatch; balance feed loading with purge unloading
4. Reduce P/F ratio if product purity collapses after a few cycles
5. Tighten time-step at valve transitions if integrator struggles
**Trade-offs**: longer cycle runs cost compute; under-converged CSS gives wrong purity/recovery.

## Pattern: Isotherm Selection by Data Shape
**When to use**: fitting q* = f(C) to experimental data.
**How**:
- Linear at low C only → use Henry (linear)
- Monolayer ceiling visible → Langmuir
- Two-step ceiling → Bi-Langmuir (two site populations)
- Smooth power-law at low C, ceiling at high C → Sips or Toth
- Multilayer (gas/solid surface) → BET
**Trade-offs**: more parameters fit better but lose physical meaning; prefer the simplest model that residuals accept.

## Pattern: Multi-Bed Parallel Synchronization
**When to use**: continuous product flow required (most industrial PSA).
**How**:
1. N beds (typically 2–12) running the same cycle, phase-offset by cycle_time / N
2. Cycle organizer drives each bed independently
3. Headers + valves sequence feed/product/purge across beds
4. Pressure-equalization steps pair beds at complementary cycle phases
**Trade-offs**: more beds smooth product flow but increase capital + valve complexity.
