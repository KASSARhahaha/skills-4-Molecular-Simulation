# Glossary — Molsim Tutorial Exercises 2026

**Acceptance ratio / acceptance probability** — probability with which an MC trial move is accepted; for Metropolis it is `min(1, exp(−βΔU))` (Ch 2, 4).

**Adaptive step tuning** — periodically rescaling `Δ_max` so the MC acceptance rate stays near 0.5 (Ch 2).

**Adsorption isotherm** — equilibrium relation `q = q(c₁, ..., c_NC)` between adsorbed amount and bulk pressure/concentration at constant `T` (Ch 7).

**Block averaging** — splitting a trajectory into blocks, computing per-block means, and using their standard deviation to estimate the standard error of an observable (Ch 2, 3).

**Boltzmann distribution** — `P(s) ∝ exp(−β U(s))`; the equilibrium distribution over configurations in the canonical ensemble (Ch 1, 2).

**CBMC (Configurational-Bias Monte Carlo)** — segment-by-segment chain growth with biased trial orientations and Rosenbluth-weight correction (Ch 7).

**CFCMC (Continuous Fractional Component MC)** — extended-ensemble scheme introducing a continuous λ-coupled fractional molecule to make insertion/deletion tractable at high density (Ch 7).

**Cheat patterns / Verlet neighbor list** — precomputed lists of nearby particle pairs to avoid O(N²) force loops (Ch 3, 5).

**Chemical potential (excess)** — `μ_ex = −k_B T ln ⟨exp(−β ΔU)⟩` from Widom insertion; what is added on top of the ideal-gas contribution (Ch 4).

**Committor** — probability `p_B(s)` that a trajectory initiated from configuration `s` with random momenta reaches state B before A; the transition state is the `p_B = 0.5` surface (Ch 6).

**Coulomb potential** — `U = q_i q_j / (4 π ε₀ r_ij)`; diverges in extended systems and needs Ewald/PME (Ch 5).

**Detailed balance** — `acc(o→n) N(o) = acc(n→o) N(n)`; the master equation condition that guarantees the equilibrium distribution is sampled (Ch 2, 4, 7).

**Ewald summation** — splitting `1/r = erfc(αr)/r + erf(αr)/r` so the real-space and reciprocal-space sums both converge (Ch 5).

**Fugacity** — effective pressure used in μVT; related to chemical potential by `βμ = βμ⁰_IG + ln(βf)` (Ch 4, 7).

**Heat capacity C_V** — fluctuation formula `C_V = (⟨U²⟩ − ⟨U⟩²) / (N k_B T²)` (Ch 2).

**Ideal-gas limit** — `β p = ρ`; verified in LJ at high `T`, low `ρ` (Ch 2).

**Langevin integrator** — `m r̈ = −γ m ṙ − ∇U + η(t)` with `⟨η(t)η(t')⟩ = 2γ k_B T δ(t−t')` (Ch 3).

**Lennard-Jones (LJ) potential** — `U(r) = 4ε[(σ/r)¹² − (σ/r)⁶]`; standard reduced units `σ` (length) and `ε` (energy) (Ch 2, 3, 5).

**Maxwell-Boltzmann distribution** — `P(v_x) ∝ exp(−m v_x² / (2 k_B T))`; check of any thermostat (Ch 3, 5).

**Metropolis algorithm** — propose + accept/reject scheme producing a Boltzmann-distributed Markov chain (Ch 2).

**Minimum image convention (MIC)** — for each pair, take the shortest distance across periodic replicas; requires box width `> 2 r_c` (Ch 1, 2).

**Müller-Brown potential** — 2-D analytic surface with multiple minima and barriers; standard test for enhanced sampling (Ch 3, 4, 5).

**Nosé-Hoover thermostat** — deterministic friction variable `ζ` coupled to kinetic energy; extended Hamiltonian `H_NH` is conserved (Ch 5).

**NVT / NVE / NPT / μVT ensembles** — fixed (N, V, T) / (N, V, E) / (N, p, T) / (μ, V, T) respectively (Ch 1–5, 7).

**Parallel tempering (replica exchange)** — N replicas at increasing `T`; swap adjacent replicas with `acc = min(1, exp((β_i − β_j)(U_i − U_j)))` (Ch 4).

**Path-CV (path collective variable)** — parametrized curve `s(σ)` from state A to B in CV-space; `σ` is the progress along the path (Ch 9).

**Periodic boundary conditions (PBC)** — particles leaving one side re-enter from the opposite side; mimics an infinite system (Ch 1, 2).

**Rosenbluth weight** — `W = Π_i w(i)` for a grown chain in CBMC; cancels the growth bias in the acceptance rule (Ch 7).

**Stirling's approximation** — `ln x! ≈ x ln x − x + ½ ln(2πx)`; valid for `x ≳ 10` (Ch 1).

**Tail correction** — analytic integral of `g(r) ≈ 1` beyond `r_c`; gives `u_tail/N`, `p_tail`, `μ_tail` for LJ (Ch 2, 5).

**Thermodynamic integration (TI)** — integrate `⟨∂U/∂λ⟩` over a coupling parameter to obtain a free-energy difference (Ch 5, 7).

**Transition path sampling (TPS)** — importance sampling of dynamical trajectories connecting two stable states via shooting/shifting moves (Ch 6).

**Umbrella sampling** — harmonic bias on a CV per window; combined via WHAM into a free-energy profile (Ch 5).

**Velocity Verlet** — symplectic, time-reversible integrator; standard for NVE MD (Ch 3).

**Well-tempered metadynamics** — variant where Gaussian height decreases with time, avoiding over-filling (Ch 9).

**Widom insertion** — virtual test-particle insertion giving `μ_ex` from the average Boltzmann factor (Ch 4).

**`g(r)` radial distribution function** — pair-correlation; `g(r) → 1` in the dilute limit and at large `r` (Ch 2).
