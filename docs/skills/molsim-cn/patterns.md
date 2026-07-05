# Patterns: Reusable Algorithms & Techniques

Curated algorithmic patterns and decision rules from "Understanding Molecular Simulation."

## 1. Sampling Construction Patterns

### 1.1 Four-step MC construction recipe (Ch 6)
For any ensemble:
1. Write target distribution π(ξ) (depends on ensemble)
2. Apply detailed balance: π(o)·α(o→n)·acc(o→n) = π(n)·α(n→o)·acc(n→o)
3. Choose trial matrix α
4. Derive acceptance rule acc

**Usage**: NPT, μVT, semigrand, Gibbs, expanded ensemble, parallel tempering.

### 1.2 Bias sampling + weight correction (Ch 12)
- Generate trial move with bias α(o→n) = f[U(n)]
- Acceptance corrects bias: acc = min[1, (f(U_o)/f(U_n))·exp(−βΔU)]
- Example: CBMC uses Rosenbluth W as weight

### 1.3 Super-detailed balance (Ch 12)
For continuous bias: impose detailed balance per fixed set of trial directions Ω*.
Simpler proof; result independent of trial count k.

## 2. Reference State Patterns (Free Energy)

### 2.1 Choose analytic reference state
| System | Reference | Why |
|---|---|---|
| Liquid | Ideal gas (B_2 known) | Density TI from low density |
| Solid | Einstein crystal (F analytic) | Springs anchored to lattice |
| Hard spheres | Carnahan-Starling F | Known closed form |
| Polymer chain | Non-self-avoiding (Q_int analytic) | Boltzmann-sample angles |

### 2.2 Thermodynamic integration template
```
1. Pick reference state with F_ref analytic
2. Define U(λ) interpolating U_ref → U_real
3. Sample ⟨∂U/∂λ⟩_λ at λ = 0, 0.1, ..., 1.0
4. Numerical integration: ΔF = ∫₀¹ ⟨∂U/∂λ⟩_λ dλ
5. Consistency check: Gibbs-Bogoliubov ∂²F/∂λ² ≤ 0
```

### 2.3 Avoid first-order phase transition
- Don't cross first-order transition (hysteresis)
- Bypass: (a) external field, (b) single-occupancy cell, (c) lattice switching, (d) bypass critical point in (T, ρ) plane

## 3. Force/Energy Patterns

### 3.1 Cutoff + tail correction
- Short-range (LJ-like): truncate at r_c, add tail ∫_{r_c}^∞ 4πr²ρu(r)dr
- Use shifted-force: u'(r) = u(r) − u(r_c) − (r − r_c)(du/dr)|_{r_c} (better energy conservation)

### 3.2 Long-range (Coulomb, dipole)
- Ewald: 1/r = erfc(αr)/r + erf(αr)/r
- PME/SPME: FFT-accelerated Fourier part, O(N log N)
- Reaction field: cutoff + ε_RF continuum outside

### 3.3 Neighbor lists (App I)
- Verlet list + skin: rebuild when max displacement > skin/2
- Cell list: divide into cells of size ≥ r_c, only check neighbors
- Combination: cell list to build Verlet list

## 4. Integration Patterns

### 4.1 Verlet family (Ch 4)
- Position Verlet: O(Δt⁴) position, O(Δt²) velocity
- Velocity Verlet: equivalent trajectory, velocity defined at integer steps
- Leap-frog: velocity at half-integer steps
- All symplectic, time-reversible, low long-term drift

### 4.2 Multiple time step (r-RESPA, Ch 14)
```
Forces split: F_fast (bonded), F_mid (short non-bonded), F_slow (long non-bonded)
Outer step Δt = n · δt
Trotter: exp(iL_slow Δt/2) · exp(iL_mid Δt/2) · [exp(iL_fast δt)]^n · exp(iL_mid Δt/2) · exp(iL_slow Δt/2)
```

### 4.3 Constraints (Ch 14)
- SHAKE: iterate per-bond correction until σ_k < tolerance
- RATTLE/LINCS: velocity-Verlet form, more stable
- Fixman correction required: +k_BT ln |H|^{1/2}

## 5. Thermostat/Barostat Patterns (Ch 7)

### 5.1 Choice based on what you trade
| Want | Use | Don't use |
|---|---|---|
| Static properties only | Andersen, Bussi | — |
| Transport (D, η, λ_T) | NVE, Nosé-Hoover, Bussi | Andersen, Berendsen, Langevin |
| Equilibration only | Berendsen | (never for sampling) |
| Polymer / Brownian | Langevin | (loses hydrodynamics) |

### 5.2 Extended Lagrangian framework
- Add variable (s for T, V for P, h for stress) to L
- New H' = K_extra + U_extra + (target constraint)
- Mass Q / inertia tunes coupling strength

## 6. Measurement Patterns (Ch 5)

### 6.1 Fluctuation → response
- C_V = (⟨E²⟩ − ⟨E⟩²)/(k_BT²) in NVT
- κ_T = (⟨V²⟩ − ⟨V⟩²)/(⟨V⟩k_BT) in NPT
- ε = 1 + (4π/(3Vk_BT))⟨M²⟩

### 6.2 Time correlation → transport
| Coefficient | Formula |
|---|---|
| D_s | ∫₀^∞ ⟨v(0)·v(t)⟩/(3) dt |
| η (shear) | (V/k_BT) ∫ ⟨σ_xy(0)σ_xy(t)⟩ dt |
| λ_T | (V/k_BT²) ∫ ⟨j_z(0)j_z(t)⟩ dt |
| σ_e | (V/k_BT) ∫ ⟨j_x^e(0)j_x^e(t)⟩ dt |

### 6.3 Error estimation
- Flyvbjerg-Petersen block average: σ from plateau
- O(n) block summing for long correlations
- Statistical error ~ 1/√N (independent samples)

## 7. Rare Event Patterns (Ch 15)

### 7.1 Bennett-Chandler recipe
```
1. Identify reaction coordinate q
2. Sample P(q) via umbrella sampling → k^TST = P(q*) D / (k_BT P(q_A))
3. Release MD from q = q* → κ = lim_{t→∞} ⟨f_B(t)⟩_release
4. k_{A→B} = k^TST × κ
```

### 7.2 Blue-Moon alternative
- Constrain q = q* in MD (existing codes: just enable constraint)
- Fixman correction: P(q) ∝ |H|^{−1/2} ⟨δ(σ)⟩_c
- Release MD from constrained configurations

### 7.3 Reaction coordinate validation
- Committor test: at q = q*, release MD, count fraction reaching B
- Ideal: histogram peaks at 0.5
- Bad: bimodal 0/1 → wrong coordinate

## 8. Acceleration Patterns (Ch 13)

### 8.1 Parallel tempering setup
- Choose T_1 < T_2 < ... < T_n
- Spacing: ΔT/T ~ 1/√N (energy fluctuation matching)
- Accept exchange: acc = min[1, exp((β_i − β_j)(U_j − U_i))]

### 8.2 Hamiltonian PT
- Each replica has different potential (e.g., polymer softness)
- Accept: acc = min[1, exp(−β(U_i^(j) − U_i^(i)) + β(U_j^(i) − U_j^(j)))]

### 8.3 Expanded ensemble
- Single system hops between β_i values
- Weights w_i = exp(β_i F − γ_i); iterate γ to equalize sampling
- Side product: F(T_i) − F(T_j)

## 9. Mesoscopic Patterns (Ch 16)

### 9.1 Preserve momentum for hydrodynamics
- Langevin: ✗ (breaks momentum)
- DPD: ✓ pairwise dissipative + random; FDT σ² = 2k_BT γ^D
- Lowe-Andersen: ✓ local pairwise

### 9.2 Mesoscopic solvent selection
| Need | Method |
|---|---|
| Static only | Implicit solvent + effective potential |
| Hydrodynamic interactions | DPD or MPC |
| Heat transport | SDPD |
| Dilute gas | DSMC |

## 10. Validation Patterns

### 10.1 Always check
- Total energy conservation in NVE (drift = bug indicator)
- ⟨P⟩ matches imposed P in NPT-MC (virial self-consistency)
- g(r) from force method vs binning (equilibrium diagnostic)
- Committor distribution at proposed q* (coordinate quality)

### 10.2 Common anti-patterns to avoid
- Cutoff Coulomb (use Ewald/PME)
- Berendsen for sampling (only equilibration)
- Force-field without Fixman correction
- Naive Verlet + constraints (use SHAKE)
- Random insertion of long chains (use CBMC)
