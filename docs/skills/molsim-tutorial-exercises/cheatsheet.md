# Cheatsheet — Molsim Practical Decisions

Decision-focused quick reference. Every line should help you *choose* under constraints.

## Ensemble Picker

| Want to sample... | Use ensemble | Constant | Typical method |
|---|---|---|---|
| Equilibrium structure at fixed T | NVT | N, V, T | MC (Ch 2), Langevin/Nosé-Hoover MD (Ch 3, 5) |
| Phase behavior, compressibility | NPT | N, p, T | MC volume moves (Ch 4), MTTK MD |
| Adsorption isotherm | μVT | μ, V, T | GC MC; CBMC/CFCMC for chains (Ch 4, 7) |
| Microcanonical dynamics, transport | NVE | N, V, E | Velocity Verlet MD (Ch 3) |

## Cutoff & Tail-Correction Decision

| Potential | Decays as | Tail correction? | Default |
|---|---|---|---|
| LJ `1/r⁶` | faster than `r⁻³` | Yes, analytic | `r_c = 2.5 σ`, add `u_tail`, `p_tail` |
| Coulomb `1/r` | `r⁻¹` | **Diverges** — use Ewald/PME | `α ≈ 3.5/L`, real + reciprocal |
| Shifted LJ | — | No (continuous at `r_c`) | Use for MD |

## Integrator Picker

| Goal | Integrator | Thermostat | Notes |
|---|---|---|---|
| NVE dynamics, transport | Velocity Verlet | none | Check energy drift `< 1e-4` |
| NVT, easy | Langevin (BAOAB) | friction γ | Destroy momentum; great for equilibrium |
| NVT, deterministic | Nosé-Hoover | friction ζ | Conserve `H_NH`; tune `Q` |
| Enhanced sampling | Metadynamics on top of any | same | Bias `V_bias → −F` |

## Step-Size & Frequency Defaults

| Move / parameter | Default | Diagnostic |
|---|---|---|
| MD timestep | `Δt = 0.001–0.005` (LJ units) | Energy drift per 10⁴ steps `< 1e-4` |
| MC max translation | tuned to 50% acceptance | disable tuning during measurement |
| MC volume move frequency | `1 / N` per step | acceptance ~0.2–0.3 |
| Replica exchange attempt | every 10–100 MD steps | swap acceptance 20–30% |
| Langevin friction γ | `~ 1 / τ_phys` (e.g. `γ ≈ 1` LJ) | too large → overdamped |
| Nosé-Hoover `Q` | `Q = g τ² Δt`, `τ ~ 100 Δt` | `H_NH` conserved |

## Common Bugs → First Check

| Symptom | Likely cause | Fix |
|---|---|---|
| Energy drifts linearly in NVE | `Δt` too large | halve `Δt`, recheck |
| `<v²>` ≠ `k_B T/m` | noise variance wrong | include `1−θ²` in Langevin σ |
| Acceptance ≈ 0 | step too large or overlapping init | reduce `Δ_max`, re-equilibrate |
| Tail-corrected pressure negative | density too low or `g(r)` not unity at `r_c` | extend `r_c`, verify `g(r)` |
| WHAM gives jagged F(s) | windows do not overlap | add windows, lower `k` |
| Replica swap never accepted | adjacent T gap too large | insert intermediate replicas |
| K-Means output unstable | features not standardized | `StandardScaler` |

## Sampling Method Picker (Barrier-Crossing)

| Barrier height | Method | Chapter |
|---|---|---|
| `≤ 2 k_B T` | plain MD/Langevin | Ch 3 |
| `3 – 8 k_B T`, known CV | Umbrella + WHAM | Ch 5 |
| `5 – 15 k_B T`, multi-CV | Metadynamics | Ch 9 |
| Complex, many CVs | Path-metadynamics | Ch 9 |
| Need true dynamics | TPS (OpenPathSampling) | Ch 6 |

## CV / Feature Choice Rules

- Use coordinates invariant under molecular symmetry (dihedrals, distances, R_g, RMSD).
- Avoid raw Cartesian for clustering — global rotation/translation dominates.
- For chains/molecules: backbone φ, ψ angles + radius of gyration are a good starter set.
- For reactions: distances that distinguish reactant, TS, product.

## Adsorption (μVT) Rules

- Always compute `W^IG` once per molecule and `T` before running CBMC.
- Convert absolute → excess adsorption before comparing to experiment.
- Fugacity (not pressure) couples to the reservoir: `βμ = βμ⁰_IG + ln(βf)`.
- At high loading, switch plain μVT → CFCMC; insertion acceptance in plain μVT can be `< 1e-6`.

## Visualization Quick Picks

| File type | Use |
|---|---|
| LJ/MD trajectory (PDB/XYZ) | NGLView in notebook (Ch 8), OVITO for large |
| MOF framework (CIF) | Mercury, OVITO |
| Protein structure (PDB) | PyMOL for renders, VS Code Protein Viewer for quick peek |
| Free-energy surface | matplotlib `pcolormesh` + contour labels |

## Cost-vs-Accuracy Tells

- "My MC result has no error bars" → unreliable.
- "My MD run has drift > 1e-3" → halve `Δt`.
- "Adjacent umbrella histograms don't overlap" → WHAM will fail.
- "100 transitions seen in 10⁶ steps" → barrier crossing resolved; <10 → undersampled.
- "Cluster populations jump on rerun" → wrong K or insufficient sampling.
