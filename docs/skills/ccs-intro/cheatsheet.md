# CCS Cheatsheet

One-line facts and equations to retrieve fast. Grouped by topic.

## Quick Numbers

| Quantity | Value |
|---|---|
| CO2 critical point | 31 C, 7.4 MPa |
| Supercritical threshold depth | ~800 m |
| CO2 density at 2 km, 20 MPa, 60 C | 725 kg/m^3 |
| Brine density (typical) | 1200 kg/m^3 |
| CO2 viscosity (reservoir) | 50-70 uPa*s |
| Brine viscosity (reservoir) | 800-1000 uPa*s |
| Air CO2 | 380-420 ppm (atmospheric) |
| Coal flue gas CO2 | 10-15% |
| Natural gas flue CO2 | 5-8% |
| IGCC syngas CO2 | 40-60% |
| MEA heat of absorption | ~85 kJ/mol CO2 |
| MEA capacity | 0.5 mol CO2/mol amine |
| MDEA capacity | 1.0 mol/mol (slower) |
| Sandstone porosity | 0.15-0.35 |
| Shale porosity | 0.06-0.15 |
| Sandstone permeability | 10^-12 to 10^-14 m^2 |
| Shale permeability | 10^-18 to 10^-21 m^2 |
| Smectite surface area | ~800 m^2/g |
| Quartz surface area | ~0.02 m^2/g |
| Sleipner injection rate | 1 Mt/yr |
| USA CO2 point-source total | ~3.4 Gm^3/yr |
| USA oilfield water injection | ~3.0 Gm^3/yr |
| USA coal fleet derate (90% capture, MEA) | ~30% |
| Committed emissions | 496 Gt CO2 (Davis 2010) |

## Capture Energy Bounds

- Minimum work of separation (90% capture, 12% CO2): **158 kJ/kg CO2**.
- Compression work to 150 bar: **218 kJ/kg CO2**.
- Total parasitic energy: **E_par = 0.75*eta_final*Q + W_comp**.
- Energy penalty for 90% capture, PC plant: **~30% of output**.
- DAC thermodynamic minimum: **~500 kJ/mol CO2** (10x post-combustion).

## Solvent Selection

| Application | Solvent | Reason |
|---|---|---|
| Coal flue gas (1 atm, 12% CO2) | MEA, KS-1, piperazine | Low partial pressure demands chemical solvent |
| Natural gas sweetening (high P) | Selexol | Physical solvent, no water |
| IGCC (40-60% CO2) | Selexol, Rectisol | High partial pressure favors physical |
| Syngas cleanup (low T) | Rectisol (-40 C) | Refrigerated MeOH; very clean |
| Benchmark for new solvents | MEA (30 wt%) | Industrial standard since 1930 |

## Membrane Quick Facts

- Solution-diffusion: **j = (P'/L)*(pR - pP)**, P' = D*H.
- 1 Barrer = 3.348x10^-19 kmol*m/(m^2*s*Pa).
- 1 GPU = 10^-6 cm^3(STP)/(cm^2*s*cmHg).
- Merkel optimum alpha for post-combustion: **20-40** (p_F/p_P = 5).
- Robeson upper bound: alpha vs P' tradeoff; new polymer must beat the line.

## Adsorption Quick Facts

- Langmuir: **theta = bp/(1+bp)**; Henry: sigma = sigma_max*b*p.
- Selectivity (Henry): **S = H_CO2/H_N2**.
- Van't Hoff: **d(ln H)/dT = delta_h_ads/(R*T^2)**.
- Zeolite 13X: delta_h_ads = 35-50 kJ/mol, S = 30-100 (CO2/N2).
- Working capacity = sigma(adsorb) - sigma(regenerate). Use this, not saturation.
- Cycle choice: **TSA for purity, VSA for throughput**.

## Storage Quick Facts

- CO2 below 800 m is supercritical; density 600-750 kg/m^3.
- Trapping cascade: structural (0-10 yr) -> residual (~20%) -> solubility -> mineral (>1000 yr).
- Capillary column height: h_max = 2*gamma*cos(theta) / [(rho_brine - rho_CO2)*g*r_throat].
- Pressure pulse radius ~100x plume radius.
- Capacity (IPCC 2005): oil/gas 675-900 Gt; coal seams 3-200 Gt; saline 1000-10^6 Gt.

## Removal Capacity Stack

| Approach | Capacity (Gt C) | Permanence |
|---|---|---|
| Avoided deforestation | 50-100 | decades |
| Afforestation | 100-200 | decades |
| Biochar | 50-100 | centuries |
| Enhanced weathering | 1000+ | millennia |
| DAC + storage | unlimited | millennia |
| Ocean fertilization | 1-10 (uncertain) | decades |
| SRM (mask only) | n/a | months |

## Key Equations

- Two-phase Darcy: **u_i = -(k*k_r,i(S_i)/mu_i)*(grad P_i - rho_i*g)**.
- Ideal selectivity: **alpha* = P'_CO2/P'_N2**.
- Stage cut: **theta = j_P/j_F**.
- Kozeny-Carman: **k ~ phi^3/(1-phi)^2**.
- Redfield: **C:N:P:Fe = 106:16:1:0.001**.
- Stefan-Boltzmann: **F = sigma*T^4** (blackbody); Earth equivalent T = 255 K.

## Project Quick Reference

| Site | Country | Mass | Depth | Notes |
|---|---|---|---|---|
| Sleipner | Norway | 1 Mt/yr | 800-1000 m | Utsira Fm; since 1996 |
| Snohvit | Norway | 0.7 Mt/yr | sub-sea | since 2008 |
| In Salah | Algeria | 1.2 Mt/yr | 1800 m | suspended 2011 |
| Weyburn | Canada | EOR | - | anthropogenic CO2 |
| Gorgon | Australia | 3.4-4.1 Mt/yr | saline | since 2015 |
| Otway | Australia | 65,000 t | 2 km | CO2CRC pilot |
| Frio I | USA | 1,600 t | 1.5 km | pilot 2004 |
