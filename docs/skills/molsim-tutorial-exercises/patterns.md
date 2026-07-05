# Patterns — Molsim Tutorial Exercises 2026

Reusable patterns from the practicals. Each pattern lists when to use it, how, and trade-offs.

## Pattern: Metropolis Accept/Reject Loop

**When**: any equilibrium MC simulation (NVT, NPT, μVT) with symmetric trial moves.

**How**:
1. Propose a trial move (translation, rotation, insertion) — proposal distribution must be symmetric.
2. Compute `ΔU = U_new − U_old`.
3. Accept with probability `min(1, exp(−βΔU))` by drawing a uniform `r ∈ [0,1]`.
4. Always count the move for the acceptance-rate statistic, whether accepted or not.

**Trade-offs**: easy to implement; acceptance collapses at high density / large steps → adaptive step tuning needed.

## Pattern: Minimum-Image Distance (Orthorhombic Box)

**When**: any periodic system with pairwise interactions.

**How** (per axis):
```
dx = x[j] - x[i]
dx -= L * round(dx / L)   # rint/round returns nearest integer
```
Keep coordinates "unwrapped" to avoid bookkeeping of how often a particle crossed a boundary.

**Trade-offs**: requires smallest perpendicular box width `> 2 r_c` or pairs are double-counted. For triclinic boxes, work in fractional coordinates `s = h⁻¹ r`, apply MIC in `s`-space, transform back.

## Pattern: Tail Correction for Truncated LJ

**When**: simulating LJ with cutoff `r_c` in any ensemble; compensates for the missing `r > r_c` interactions.

**How**: assuming `g(r) = 1` beyond `r_c`,
- `u_tail/N = (8/3) π ρ ε σ³ [ (2/3)(σ/r_c)³ − (σ/r_c)⁹ ]`
- `p_tail = (16/3) π ρ² ε σ³ [ (2/3)(σ/r_c)³ − (σ/r_c)⁹ ]`

**Trade-offs**: only valid for potentials decaying faster than `r⁻³`. For Coulomb use Ewald/PME instead.

## Pattern: Adaptive Step Tuning

**When**: any MC code where you want stable ~50% acceptance across densities.

**How**:
1. Every `K` steps compute `f = accepted / K`.
2. Scale `Δ_max *= f / 0.5` (multiplicative).
3. Reset counter; disable during measurement runs to preserve detailed balance.

**Trade-offs**: target 0.5 is empirical for translations; rotations and volume moves want different targets (often 0.2–0.3).

## Pattern: Block-Averaged Error Bars

**When**: any long simulation; you need a standard error on an observable.

**How**:
1. Drop the first `~10%` of the trajectory for equilibration.
2. Split the rest into `M ≈ 10–20` blocks.
3. Compute the per-block mean `x_i`.
4. Standard error = `std(x_i) / √(M − 1)`.

**Trade-offs**: blocks must be longer than the autocorrelation time or the error is underestimated. If SE decreases as you cut fewer blocks, the blocks are too short.

## Pattern: Conserved-Quantity Diagnostic

**When**: validating an integrator or thermostat.

**How**:
- NVE: total energy `E = K + U`; drift should be `< 1e-4 · |E|` over the run.
- Nosé-Hoover NVT: extended Hamiltonian `H_NH = K + U + ½ p_ζ²/Q + g k_B T ln ζ`.
- NPT: includes `+ P V` and box-kinetic terms.

**Trade-offs**: a perfectly conserved `H` does not guarantee correct *sampling*; always also check the velocity distribution.

## Pattern: Replica-Exchange (Parallel Tempering) Move

**When**: free-energy landscape with high barriers; one temperature alone cannot cross.

**How**:
1. Run `N` replicas at temperatures `T_1 < T_2 < ... < T_N`.
2. Every `K` steps pick adjacent pair `(i, i+1)` at random.
3. Accept swap with `min(1, exp((β_i − β_{i+1})(U_i − U_{i+1})))`.

**Trade-offs**: requires enough replicas that adjacent energy histograms overlap (~20–30% swap acceptance). Cheaper on shared memory; on clusters needs MPI.

## Pattern: Widom Insertion for μ_ex

**When**: low-to-moderate density NVT; need the chemical potential.

**How**:
1. Every `K` steps, place a test particle at a random position.
2. Compute `ΔU` of the (virtual) insertion.
3. Accumulate `⟨exp(−β ΔU)⟩`; then `μ_ex = −k_B T ln ⟨exp(−β ΔU)⟩`.

**Trade-offs**: variance explodes at high density → switch to CBMC or CFCMC (Ch 7).

## Pattern: Langevin Thermostat (Half-Step)

**When**: cheap, robust NVT sampling where momentum conservation is not required.

**How**:
- `θ = exp(−γ Δt)`
- `σ = √(k_B T (1 − θ²))`
- velocity update: `v ← θ v + (Δt/2) F/m + σ ξ`, with `ξ ~ N(0, 1)`.

**Trade-offs**: overdamped at large `γ`; destroys hydrodynamics. Use Nosé-Hoover when dynamics matter.

## Pattern: Ewald Split of Coulomb

**When**: periodic charged system.

**How**:
- Real-space sum (cutoff ~`5/α`): `Σ_i<j erfc(α r_ij)/r_ij`.
- Reciprocal-space Fourier sum (k-cutoff ~`5 α / V^(1/3)`): `Σ_k (4π/V) exp(−k²/(4α²)) |Σ_j q_j exp(i k · r_j)|² / k²`.
- Add self and surface (tin-foil) corrections.

**Trade-offs**: `O(N^(3/2))` cost; for large systems PME/PPPM at `O(N log N)` is preferred.

## Pattern: Umbrella Sampling + WHAM

**When**: free-energy profile along a known CV with barriers `5–15 k_B T`.

**How**:
1. Pick `M` windows with bias centers `s_i` and spring constant `k`.
2. In window `i` run biased sampling with `U_bias = ½ k (s − s_i)²`.
3. Collect histograms `P_i(s)`.
4. Run WHAM (or MBAR) to combine them into unbiased `F(s)`.

**Trade-offs**: requires overlap between adjacent histograms; for unknown CVs use metadynamics (Ch 9).

## Pattern: K-Means Conformational Clustering

**When**: discovery of metastable states from a trajectory.

**How**:
1. Compute per-frame features (`R_g`, RMSD, dihedrals).
2. Standardize features.
3. Run sklearn KMeans for `K = 1..10`; pick `K` from the elbow + silhouette.
4. Inspect the highest-populated clusters' representative structures.

**Trade-offs**: Euclidean on dihedrals ignores periodicity; use `min(|d|, 2π − |d|)` or a circular distance. For time-aware clustering use GMM/TICA.

## Pattern: Adaptive Path-Metadynamics

**When**: complex transition requiring many CVs (>3); plain metadynamics infeasible.

**How**:
1. Define N CVs `z_i(q)` and states A, B.
2. Initialize M path nodes `s_j` between A and B.
3. Run metadynamics biasing only the 1-D progress `σ(q)` along the path.
4. Periodically reparametrize nodes to be equidistant and shift them toward `⟨z⟩_σ`.

**Trade-offs**: convergence of the path can be slow; check stability by re-running with different initial guesses.

## Pattern: Two-Way MD Shooting for TPS

**When**: unbiased sampling of rare transition paths.

**How**:
1. From a frame in the current path, perturb momenta slightly.
2. Integrate forward and backward in time.
3. Accept the new path if it still connects A → B.

**Trade-offs**: low acceptance if the perturbation is too large or the frame is far from the transition state; warm-start from a committor analysis.

## Pattern: WHAM / MBAR Histogram Combination

**When**: combining biased simulations (umbrella, replica exchange) into an unbiased FES.

**How**: iteratively solve self-consistent equations mixing biased histograms with their bias factors. MBAR is the modern, lower-variance successor.

**Trade-offs**: WHAM needs overlap; MBAR is more robust but heavier to set up.

## Pattern: EDA Before ML

**When**: applying any clustering/regression to trajectory features.

**How**: ranges, NaNs, outliers, correlations (`pandas.plotting.scatter_matrix`), histogram per feature.

**Trade-offs**: 5 minutes of EDA often saves hours of misapplied ML.
