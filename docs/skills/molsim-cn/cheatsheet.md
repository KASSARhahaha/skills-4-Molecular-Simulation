# Cheatsheet: Quick Decision Tables

## 1. Integrator Selection

| Goal | Pick | Why |
|---|---|---|
| General MD | Velocity Verlet | Symplectic, simple, low drift |
| Constrained bonds | Verlet + SHAKE/LINCS | Allows 2× larger Δt |
| Multi-timescale | r-RESPA | Bonded fast, non-bonded slow |
| High accuracy needed | (rare) Runge-Kutta | Not for long MD (drift) |

**Never**: Forward Euler (drift), high-order predictor-corrector (Lyapunov defeats accuracy).

## 2. Ensemble Selection

| Need | Ensemble | Method |
|---|---|---|
| Equation of state | NPT | NPT-MC or NPT-MD |
| Chemical potential | μVT | μVT-MC (small molecules) |
| Gas-liquid coexistence | Gibbs | Gibbs ensemble MC |
| Adsorption in pores | μVT | μVT-MC (PBC + framework) |
| Phase diagram tracking | NPT + Gibbs-Duhem | Trace coexistence curve |
| Transport properties | NVE | No thermostat disturbance |

## 3. Thermostat Decision

```
Do you need transport (D, η, λ_T)?
├─ YES → Use NVE (best) or Nosé-Hoover / Bussi (momentum-preserving)
│        NEVER use Andersen, Langevin, Berendsen
└─ NO (static properties only)
    ├─ Equilibration phase → Berendsen (fast, weak coupling)
    │                       BUT never sample with Berendsen
    └─ Sampling → Bussi (recommended) or Nosé-Hoover
```

## 4. Free Energy Method Selection

| System | Reference state | Method |
|---|---|---|
| Liquid (any T, ρ) | Ideal gas | Density TI (low ρ → target) |
| Solid | Einstein crystal | Hamiltonian TI |
| Hard sphere | Carnahan-Starling | Coupling-softened TI |
| Crystal polymorphs (fcc/hcp) | Either | Lattice switching (Bruce) |
| Polymer | Non-self-avoid chain | Rosenbluth sampling |
| Small Δλ steps | — | FEP (Zwanzig) |
| Many λ steps + reuse | — | Histogram reweighting (MBAR) |

## 5. Long-Range Force Method

| System size N | Method |
|---|---|
| N < 10³ | Plain Ewald (gold standard) |
| 10³ ≤ N ≤ 10⁵ | PME / SPME (FFT accelerated) |
| N > 10⁵ | FMM (O(N)) |
| Isotropic + low precision | Wolf / Reaction field |
| Surface systems | 2D Ewald (special form) |

## 6. Reaction Coordinate Quality

Test: at proposed q = q*, sample configurations, release MD, count fraction reaching B.

| Committor histogram | Verdict |
|---|---|
| Single peak at 0.5 | Excellent q* |
| Broad distribution | Marginal, may still work |
| Bimodal at 0 and 1 | Wrong coordinate, find new q |

## 7. Acceptance Rule Patterns

| Move type | Acceptance |
|---|---|
| NVT displacement | min(1, e^{−βΔU}) |
| NPT volume change | min(1, e^{−β[ΔU + PΔV − Nk_BT ln(V'/V)]}) |
| μVT insertion | min(1, (fV/(N+1))·e^{−βΔU}) |
| μVT deletion | min(1, (N/(fV))·e^{−βΔU}) |
| Semigrand species swap | min(1, e^{−β(Δμ + ΔU)}) |
| Parallel tempering swap | min(1, e^{(β_i − β_j)(U_j − U_i)}) |
| CBMC regrow | min(1, W(n)/W(o)) |

## 8. Common Pitfalls → Fix

| Symptom | Likely cause | Fix |
|---|---|---|
| Energy drift in NVE | Non-symplectic integrator or bug | Use Verlet; check force code |
| Wrong T from kinetic energy | Using v instead of p·v | Use ⟨p·v⟩ for T |
| Cutoff Coulomb artifacts | Truncating 1/r | Use Ewald/PME |
| Berendsen gives wrong fluctuations | Not a valid ensemble | Use Bussi/Nosé-Hoover for sampling |
| Chain insertion always rejected | Random insertion in dense fluid | Use CBMC |
| Wrong constraint statistics | Missing Fixman correction | Add +k_BT ln |H|^{1/2} |
| g(r) noisy at small r | Binning artifacts | Use Borgis force method |
| S(q) negative or oscillatory | Fourier from truncated g(r) | Compute S(q) directly from ρ(q) |
| Low PT acceptance | Temperature spacing too wide | ΔT/T ~ 1/√N |

## 9. Key Formulas to Memorize

- **Boltzmann**: ρ ∝ e^{−βU}, β = 1/(k_BT)
- **Detailed balance**: π(o)α(o→n)acc(o→n) = π(n)α(n→o)acc(n→o)
- **Metropolis**: acc = min(1, e^{−βΔU})
- **Einstein D**: D = lim_{t→∞} ⟨|Δr(t)|²⟩/(2dt)
- **Green-Kubo D**: D = ∫₀^∞ ⟨v_x(0)v_x(t)⟩dt
- **TI**: ΔF = ∫₀¹ ⟨∂H/∂λ⟩_λ dλ
- **FEP (Widom/Zwanzig)**: ΔF = −k_BT ln⟨e^{−βΔU}⟩
- **Gibbs-Bogoliubov**: ∂²F/∂λ² ≤ 0 (linear interpolation)
- **Clausius-Clapeyron**: dP/dT = Δh/(T Δv)
- **Fixman**: ΔF_correction = k_BT ln |H|^{1/2}

## 10. Conversion Constants

- k_BT at 300 K ≈ 4.116 × 10⁻²¹ J ≈ 0.593 kcal/mol ≈ 2.479 kJ/mol
- 1 fs = 10⁻¹⁵ s; LJ τ* ≈ 1 ps
- Bjerrum length λ_B in water at 300 K ≈ 7.14 Å
- 1 atm ≈ 1.013 × 10⁵ Pa; 1 bar = 10⁵ Pa
