# CCS Patterns and Anti-patterns

Recurring practitioner heuristics, structured as Patterns (do this) and Anti-patterns (don't do this).

## Cross-cutting decision patterns

### Pattern: Match capture technology to CO2 partial pressure
- 0.04% (air) -> DAC; 10-15% (coal flue) -> MEA or VSA; 5-8% (gas flue) -> advanced solvents; 30-60% (IGCC) -> Selexol/Rectisol/membranes.
- Reason: thermodynamic minimum work scales with ln(1/y_CO2). Higher partial pressure -> cheaper separation.
- Apply at: first-pass technology selection.

### Pattern: Use scale parity to test plausibility
- USA point-source CO2 ~3.4 Gm^3/yr vs USA oilfield water injection ~3.0 Gm^3/yr. Comparable.
- Conclusion: CCS injection at Mt/yr scale is engineering, not science fiction.
- Apply at: portfolio-level feasibility arguments.

### Pattern: Sandstone-shale dichotomy for site design
- Store in sandstone (k ~10^-12 to 10^-14 m^2); seal with shale (k ~10^-18 to 10^-21 m^2).
- Mixing both into one formation ruins both functions.
- Apply at: reservoir screening.

### Pattern: Read Henry coefficients in context
- Pure water kappa ~1600 at 25 C, 1 atm -> useless for CO2.
- 30% MEA kappa ~50 -> 30x better.
- Compare solvents at the same partial pressure and temperature, not at saturation.
- Apply at: solvent selection.

## Quantitative pattern library

### Pattern: Minimum work of separation scales with ln(1/y)
- W_min = -R*T*[ ln(y/(1-eta)) + ((1-eta*y)/(1-eta))*ln((1-y)/(1-eta*y)) ] / M_CO2.
- y = 0.12 (coal flue), eta = 0.9: W_min ~ 158 kJ/kg CO2.
- y = 0.0004 (air), eta = 0.5: W_min ~ 500 kJ/mol CO2.
- Use to bound capture energy before detailed simulation.

### Pattern: McCabe-Thiele triangle count
- Step off triangles between equilibrium line (y = kappa*x) and operating line.
- Solvent with lower kappa -> fewer stages -> cheaper column.
- L/G ratio shift changes slope; over-recompense increases steam cost.
- Use for first-pass absorber sizing.

### Pattern: Capillary entry pressure sets column height
- h_max = 2*gamma*cos(theta) / [(rho_brine - rho_CO2) * g * r_throat].
- Shale throat ~10 nm -> entry pressure ~6 MPa -> holds km-scale columns.
- Sandstone throat ~10 um -> m-scale columns -> not a seal.
- Use to screen caprocks.

### Pattern: Plume vs pressure footprint asymmetry
- Pressure pulse propagates ~100x farther than CO2 plume.
- Pressure monitoring watches seal integrity; 4D seismic watches plume extent.
- Use to design MVA strategy.

### Pattern: Capillary heterogeneity controls plume architecture
- Cm-scale shale laminae within sandstone act as local capillary barriers.
- CO2 piles up beneath them, increasing residual trapping far above lab-scale prediction.
- Always grid heterogeneity explicitly; don't average it away.

## Capture technology anti-patterns

### Anti-pattern: Citing saturation capacity for cyclic operation
- Saturation capacity (sigma_max) at 1 bar pure CO2 is irrelevant for 12% flue gas.
- Working capacity = sigma(adsorb) - sigma(regenerate). Use this.
- Solvent analog: capacity 1.0 mol/mol is reachable only at stoichiometric regeneration; real systems operate at 0.4-0.7.

### Anti-pattern: Forgetting solvents are not neat CO2-water systems
- Real flue gas has SOx, NOx, O2, particulates.
- MEA in such streams loses 1-2 kg/t CO2 -> major OPEX and nitrosamine by-products.
- Pretreat SOx and NOx BEFORE the amine unit.

### Anti-pattern: Chasing selectivity at low partial pressure
- For post-combustion flue gas at p_F/p_P = 5, alpha > 50 is wasted (Merkel).
- The pressure ratio, not the material, sets the floor on purity.
- Fix pressure ratio first; pick material second.

### Anti-pattern: Equating EOR with CCS
- CO2-EOR is optimized to recycle CO2 and maximize oil, not store CO2.
- CCS accounting must track net: (CO2 injected) - (CO2 produced with oil).
- Some EOR projects are net emitters, not sinks.

## Storage anti-patterns

### Anti-pattern: Treating caprock as impermeable
- Caprock has finite k (~10^-20 m^2). "Seal" means capillary entry pressure exceeds buoyancy pressure at planned column height.
- Always test via MICP (Mercury Injection Capillary Pressure).
- One high-k lens in caprock short-circuits the seal -> require stochastic realizations.

### Anti-pattern: Single-phase injection simulation
- Ignoring brine displacement underestimates pressure buildup 2-3x.
- Pressure limit (fracture gradient) often binds before plume footprint.

### Anti-pattern: No hysteresis in relative permeability
- Symmetric k_r curves overestimate CO2 mobility after imbibition.
- Underpredicts residual trapping by 30-50%.
- Always use Land trapping model or equivalent for hysteretic k_r.

### Anti-pattern: Forgetting legacy wells
- Each pre-existing well through the caprock is a potential leakage path.
- Inventory abandoned wells in the area of review; model them explicitly.

## Removal and geoengineering anti-patterns

### Anti-pattern: Citing GPP as sequestration
- Gross primary productivity is mostly returned by respiration within a year.
- Net ecosystem productivity (NEP) is what accumulates. NEP is typically 10-30% of GPP.

### Anti-pattern: Counting avoided emissions as negative emissions
- Avoided deforestation prevents a positive flux; it is not a negative flux.
- Removal requires actual CO2 take-up from the atmosphere.

### Anti-pattern: Treating SRM as substitute for mitigation
- SRM masks warming but CO2 keeps accumulating.
- Ocean acidification continues; termination shock is a real risk.
- Treat SRM as a complement, not a substitute.

## Cross-chapter modeling discipline

### Pattern: Use depth-density diagram to predict CO2 phase
- Above 800 m: gas (200-700 kg/m^3).
- Below 800 m: supercritical (600-750 kg/m^3).
- Determines injectivity, plume behavior, storage efficiency.

### Pattern: Use Robeson upper bound for polymer membranes
- New material must lie on or above the bound to be a real improvement.
- Otherwise it is an incremental tradeoff, not a breakthrough.

### Pattern: Use Otway blind prediction as simulator benchmark
- Any new code must reproduce Otway plume geometry within 2x.
- Else: model is suspect.

### Pattern: Sort removal options by permanence
- Mineral millennia > biochar centuries > forests decades > SRM not-at-all.
- Permanence hierarchy helps triage when budgets are limited.
