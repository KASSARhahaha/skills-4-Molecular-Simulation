# Patterns — Adsorption Simulation Notes

## Pattern: Aspen PSA Build Order
**When to use**: every new Aspen Adsorption model, PSA or TSA, dynamic or gCSS.
**How**:
1. Components + property method → generate `*.aprbkp`.
2. Lay out feeds, beds, valves, tanks, products on the flowsheet.
3. Specify feed (T, P, y, F) and bed geometry (L, D, ε, ρ_b).
4. Pick isotherm family and enter IP parameters per component.
5. Pick kinetic model (default LDF solid film) and k_i per component.
6. Specify steady Cv values on flow-setting valves.
7. Add Cycle Organizer and define step table.
8. Initialize with small non-zero loading; run.
**Trade-offs**: skipping step 1 (property file) is the #1 cause of inactive Run buttons; skipping step 7 gives breakthrough, not a cycle.

## Pattern: Skarstrom 4-Step Cycle with Optional Equalization
**When to use**: designing a basic two-bed PSA for a new gas separation.
**How**: define 4 steps (adsorb, blowdown, purge, repressurize) per Note 4 Table; add a 5th pressure-equalization step between blowdown and repressurization on each bed.
**Trade-offs**: equalization adds ~10–30% compressor-work saving but extends cycle time and adds a control valve.

## Pattern: One-Bed Dynamic Then Two-Bed gCSS
**When to use**: H₂-PSA and any new PSA design study.
**How**: run a one-bed breakthrough first to read t_bt; size the adsorption step as 0.5–0.8·t_bt; transfer schedule to two-bed gCSS for design-grade CSS purity/recovery.
**Trade-offs**: dynamic gives intuition and step sizing; gCSS gives design-grade CSS but converges only with a sensible initial guess.

## Pattern: Reduce Vendor Isotherm Algebraically Before Constraining
**When to use**: whenever a vendor supplies a multi-step isotherm "recipe" (P_sat, P_rel, exponentials).
**How**: identify that `Y_i·P_total = p_i`; merge constant prefactors; use `e^(a+b)=e^a·e^b` to peel constants; aim for a single-line form matching a built-in isotherm family.
**Trade-offs**: ~30 min of algebra vs hours of debugging a Flowsheet Constraint implementation that Aspen solves slower and less stably.

## Pattern: LDF Solid Film as Default Kinetic Model
**When to use**: every gas_bed block unless you have breakthrough data showing LDF fails.
**How**: select "solid form" in the Kinetic Model tab; enter one k_i per component; keep k_i constant.
**Trade-offs**: simplest stable choice; hides the real film/pore/surface resistances but only matters when fitting breakthrough data.

## Pattern: Throttle Product/Purge Valves with Cv, Never On/Off
**When to use**: every PSA cycle.
**How**: set VP and Vpurge as AS=2 (Cv valve) with a small Cv (order 1e-5 to 1e-6 kmol/(s·bar)); let VF and VW be larger Cv.
**Trade-offs**: on/off valves (AS=0/1) cause pressure spikes and back-mixing of impurities into the product.

## Pattern: Two Conventions for Intra-particle Voidage
**When to use**: any time you read a voidage from a paper or vendor datasheet.
**How**: chromatography convention is pore volume over total column volume; chemical engineering convention is pore volume over particle volume. State which one you are using before any balance.
**Trade-offs**: silently mixing conventions introduces a factor-of-several holdup error.

## Pattern: Pressure-Drop Model by Reynolds Number
**When to use**: choosing between Carman–Kozeny, Burke–Plummer, and Ergun.
**How**: Re_p < ~10 → Carman–Kozeny; Re_p > ~2000 → Burke–Plummer; otherwise Ergun (combines both).
**Trade-offs**: Ergun is always safe and is the Aspen default; do not overthink the choice.

## Pattern: EOS Selection by Pressure and Composition
**When to use**: specifying the gas property method in Aspen.
**How**: P < 2–3 bar → ideal gas; up to ~10 bar hydrocarbon → RK-Soave; refinery/H₂/CO/CO₂ streams or higher P → Peng–Robinson.
**Trade-offs**: stiffer EOS costs solver time but avoids 10–20% volumetric-flow errors at high P.

## Pattern: Fixing TSA Step-Change Temperature Spikes
**When to use**: when TSA simulations show periodic T spikes at step transitions.
**How**: enlarge tank volumes (default is tiny); disable gas compressibility in the tanks (keep it in the beds).
**Trade-offs**: slight loss of physical fidelity in the tank; large gain in numerical stability.

## Pattern: Initialize gCSS from a Dynamic Run
**When to use**: gCSS refuses to converge from a cold start.
**How**: run the same flowsheet in dynamic mode for a few cycles; export the loading profile; import it as the gCSS initial state.
**Trade-offs**: extra setup time; converges cases that otherwise fail.

## Pattern: Define Repeated Values as Flowsheet Constraints
**When to use**: when N blocks share the same Cv, k_i, or T.
**How**: declare `X as RealVariable;` once, assign `X = value;`, reference `X` in each block's spec.
**Trade-offs**: one change point instead of N; small overhead per constraint.
