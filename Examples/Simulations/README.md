# Simulation Examples

`01_*` folders are generic sweep templates. They all use a 600 l/mm blazed Au/Si grating as a worked example; swap the geometry and stack block to adapt them to any grating.

- `01_Fixed_Angle_Energy_Sweep` — energy sweep at fixed α
- `01_Fixed_cff_Energy_Sweep` — energy sweep at fixed C_ff; α is recomputed from the grating equation at each energy step
- `01_Fixed_Energy_Angle_Sweep` — angle sweep at fixed energy

`02_*` folders simulate specific physical  measured gratings, benchmarked against `../Reference_Data/`.

- `02_Blazed Grating 600L` — HORIBA 600 l/mm blazed, Au + optional C contamination; multiple result sets for different polarisations, angles, and contamination states
- `02_Blazed Multilayer 2400L` — 2400 l/mm blazed Cr/C multilayer, Bragg sweep along a lookup table, 3–5 keV range
- `02_Laminar Grating 400L` — 400 l/mm laminar (trapezoidal) Pt coating, results for TE/TM and contaminated/clean at two incidence angles

Each `Results/` folder contains `simulation_results.csv` and PNG plots where available.
