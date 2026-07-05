# Cheatsheet — Simulation of Adsorption Processes

## Cycle Configuration Decision
| If regeneration driver is… | Use | Don't use because… |
|---|---|---|
| Pressure swing | **PSA** / VPSA | TSA too slow for bulk gas |
| Temperature swing | **TSA** | PSA cannot displace strongly-bound species |
| Both | **PTSA** | (extra complexity only when one driver is insufficient) |
| Vacuum available | **VPSA** | higher CAPEX but better recovery for O2 / CO2 |

## Adsorbent Selection
| Target species | First choice | Why |
|---|---|---|
| N2 from air | Zeolite 13X | equilibrium selectivity N2 ≫ O2 |
| O2 from air (modern) | **LiLSX** | ~2× working capacity vs classic 5A/13X |
| O2 from air (classic) | Zeolite 5A or 13X | original Skarstrom-cycle adsorbents |
| H2O (drying) | Silica gel, activated alumina | polar surface binds water |
| VOC / organics | Activated carbon | hydrophobic, large surface area |
| H2S / sulfur | Zeolite, promoted carbon | strong chemisorption |
| n-paraffins | Zeolite 5A | molecular sieve |

## Isotherm Selection Decision
```
Data shape →
  Linear only at low C      → Henry (linear)
  Monolayer ceiling         → Langmuir            (DEFAULT)
  Two distinct ceilings     → Bi-Langmuir
  Power-law low, ceiling high → Sips or Toth
  Multilayer (surface area) → BET
```
**Rule**: prefer the simplest model whose residuals are random. Adding parameters always lowers SSE but may not improve physics.

## LDF vs Fickian Kinetics
| Situation | Pick | Driving force |
|---|---|---|
| Single dominant resistance (pore OR film) | **LDF (lumped)** | q* − q (solid film) or C_bulk − C_surface (gas film) |
| Pore diffusion explicitly resolved | Fickian | ∂q/∂t = D_eff ∇²q inside particle |
| Default for column simulation | LDF | — |
| Default k_i starting estimate | k ≈ 0.4 / t_0.5 | — |

## Numerical Defaults (Aspen Adsorption)
| Setting | Default | Escalate when… |
|---|---|---|
| Discretization | UDS1 | breakthrough front smeared > 10 % of bed length |
| Nodes | 20–60 | long bed (> 2 m) or sharp front |
| Time integrator | Implicit DAE | valve transitions cause spikes → smaller dt |
| Material balance | Convection only | axial dispersion matters for short beds |
| Momentum | Ergun | always for packed beds |
| Mass balance term | include (1−ε)/ε · ∂q/∂t | never drop this |

## PSA Cycle Tuning Levers (ordered by impact)
1. **P/F ratio** — biggest purity/recovery trade-off. Start 0.1–0.3; raise for purity, lower for recovery.
2. **Feed step duration** — must keep MTZ inside bed; too long → breakthrough.
3. **Pressure ratio** (P_feed / P_purge) — higher → sharper regeneration.
4. **Pressure equalization steps** — improve recovery 20–40 % in 2-bed systems.
5. **Bed length / diameter** — increases working capacity but raises ΔP.

## Skarstrom Cycle Directionality
| Step | Direction | Why |
|---|---|---|
| Pressurization | co-current with feed | raises P using feed gas |
| Adsorption | co-current with feed | produces raffinate (LP) |
| Blowdown | **countercurrent** to feed | strips strongly-adsorbed species from feed end |
| Purge | **countercurrent** to feed | slips raffinate gas to push adsorbate out feed end |
| Pressure equalization | bed-to-bed | shares pressure; co-current for receiver, countercurrent for donor |

## PSA Valve Code Convention (Aspen / Table 1.5)
| Code | Meaning |
|---|---|
| 0 | Closed |
| 1 | Open continuously |
| 2 | Flow-controlled |
| 3 | Pressure-controlled |

## CSS Convergence Tells
| Symptom | Likely cause | Fix |
|---|---|---|
| q profile drifts each cycle | bed overloaded | shorten feed step, lengthen purge |
| Purity collapses after 5 cycles | insufficient purge | raise P/F ratio |
| T spikes at step transitions | numerical, not physical | tighten dt at transitions |
| Oscillating q profile | step timing mismatch | balance feed vs purge loading |
| No convergence after 30 cycles | wrong isotherm sign or unit | re-validate against breakthrough |

## Workshop Build Order (Aspen)
```
1. Components  →  2. Property method  →  3. Flowsheet  →
4. gas_bed specs (geometry, adsorbent, isotherm)  →
5. Kinetic (LDF, k_i)  →  6. Discretization (UDS1, nodes)  →
7. Single-feed breakthrough (validate)  →
8. Cycle organizer (multi-step PSA)  →  9. CSS convergence
```
**Never** reverse this order — late component/property changes invalidate downstream specs.

## Common Pitfalls (Aspen)
- Mixing units for k_i: Aspen wants 1/s; using 1/min silently scales by 60×.
- Forgetting (1−ε)/ε factor in solid term: breaks mass balance.
- Using constant ΔP instead of Ergun: distorts velocity during blowdown.
- Skipping the breakthrough run before assembling a cycle: hides model errors.
- Ideal-gas property method at P > 10 bar: distorts gas density, especially H2.

## Quick Numbers
- ε_bed for spheres: 0.36–0.42
- Particle shape factor ψ: 1.0 (sphere), ~0.65 (Raschig rings), ~0.75 (crushed
- Zeolite 13X N2/O2 selectivity: ~5–6 at 25 °C
- LiLSX working capacity (O2 PSA): ~2× of classic 5A/13X
- LDF k_i for N2 on 13X: ~0.08 s⁻¹ (typical, 2 mm beads)
- Typical PSA cycle time: 30 s – 5 min
- Typical TSA cycle time: 30 min – several hours
- Pressure ratio for PSA: 4–10
- P/F ratio starting value: 0.1–0.3
- Skarstrom O2 PSA scale: 1 L/min (medical) to tens of tons/day (industrial)
