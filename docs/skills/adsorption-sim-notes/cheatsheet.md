# Cheatsheet — Problem-Solver's Companion

## Decision: PSA vs TSA
| Symptom / Requirement | Choose |
|---|---|
| Target moderately adsorbed; fast cycles needed | **PSA** |
| Target strongly adsorbed; P-swing alone cannot desorb enough | **TSA** |
| Heat integration available (steam, waste heat) | **TSA** |
| Cycles/min required (small units) | **PSA** |
| Refinery H₂, O₂ from air, N₂ from air | **PSA** |
| Air drying, solvent removal, CO₂ with amine-modified solids | **TSA** |

## Decision: Which Isotherm
| Data shape / situation | Pick |
|---|---|
| Homogeneous surface, monolayer, default start | **Langmuir** |
| Two distinct site populations | **Dual-site Langmuir** |
| Heterogeneous, finite saturation, gas/zeolite | **Toth** or **Sips** |
| Heterogeneous, bioseparation, narrow p range | **Freundlich** (never extrapolate) |
| Surface area measurement (N₂, 77 K) | **BET** |
| Trace/dilute region of any isotherm | **Linear (Henry)** |

## Decision: Which Aspen Block
| Goal | Block |
|---|---|
| Learn breakthrough, size step from t_bt | gas_bed + dynamic |
| Design-grade CSS for two-bed PSA | **gCSS** |
| Buffer tank in single-bed variant | gas_interaction |
| Schedule valve state changes | Cycle Organizer |
| Hold gCSS step-time metadata | CSS_Info |

## Decision: Valve Spec Code (Aspen AS)
| Need | Code |
|---|---|
| Always closed | 0 |
| Always open, flow by mass balance | 1 |
| Flow ∝ ΔP (Cv valve) — **default for product/purge** | 2 |
| Fixed flow setter | 3 |

## Decision: Kinetic Model
| Situation | Pick |
|---|---|
| Default, no breakthrough data | **LDF solid film**, constant k_i |
| Film-controlled system | LDF gas film |
| LDF fails to fit breakthrough | Micropore + macropore split |
| Have k_i(T,q) data | Variable k_i |

## Decision: EOS by Pressure
| Max P | EOS |
|---|---|
| < 2–3 bar | Ideal gas |
| ≤ 10 bar, hydrocarbons | RK-Soave |
| Refinery, H₂/CO/CO₂, higher P | **Peng–Robinson** |

## Quick Diagnosis Map
| Symptom | First thing to check |
|---|---|
| Run button inactive | Property file `*.aprbkp` not generated |
| Solver fails at t=0 | Initial loading exactly zero; use ~1e-6 mol/kg |
| PSA purity below target | Under-purging → lengthen purge step |
| T spikes at step boundaries (TSA) | Tank volumes too small + compressibility on |
| Custom isotherm returns NaN | Vendor expression has div-by-~0; reduce algebraically |
| ΔP/L way off vs expectation | ε or d_p in wrong units; verify SI vs cgs vs Aspen |
| gCSS diverges from cold start | Initialize from a dynamic run |
| Simulated purity differs from paper | Valve type (Cv vs on/off) — small Cv raises purity |
| Bed saturates far slower than experiment | k_i too small (τ_i too large) |
| Breakthrough too fast vs experiment | Don't increase q_max — reduce k_i or check isotherm |

## Thresholds & Defaults
| Quantity | Rule of thumb |
|---|---|
| Adsorption step duration | 0.5–0.8 × t_bt |
| Pressure-equalization duration | 10–20% of adsorption step |
| Cv for VP/Vpurge valves | order 1e-5 to 1e-6 kmol/(s·bar) |
| Pressure ratio (PSA) | feed P / purge P ≥ 4–10 |
| Initial loading | ~1e-6 mol/kg (avoid exact zero) |
| ΔP/L sanity ceiling (low-drop design) | ~0.1 bar/m |
| Re boundary (Ergun terms) | laminar < ~10; turbulent > ~2000 |
| Pressure-equalization energy saving | ~10–30% of compressor work |

## Build Order Checklist (every PSA/TSA model)
1. Components + property method → `*.aprbkp`
2. Flowsheet: feeds, beds, valves, tanks, products
3. Feed T, P, y, F
4. Bed geometry (L, D, ε, ρ_b)
5. Isotherm family + IP parameters per component
6. Kinetic model + k_i per component
7. Steady Cv on product/purge valves
8. Cycle Organizer step table
9. Initialize with non-zero loading
10. Run → check CSS purity/recovery

## Cycle Step Tables (Skarstrom + Pressure Equalization)
**4-step**: adsorb → blowdown → purge → repressurize.
**6-step (with EQ)**: adsorb → EQ → blowdown → purge → EQ → repressurize.
