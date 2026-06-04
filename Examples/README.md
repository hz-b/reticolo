# Examples

`Gratings/` — scripts that build a stack and render a cross-section meshgrid. No sweep is run; use these to inspect geometry before simulating.

`Simulations/` — complete RCWA runners with pre-computed `Results/`.

`Optical_Constants/` — CXRO `.txt` refractive-index files (columns: Energy eV, Delta, Beta). `ELISA/` holds instrument-specific CSV files.

`Reference_Data/` — experimental measurements and third-party results (REFLEC, DiffractMod) for benchmarking.

All runners resolve paths relative to their own location and add `helpers/` automatically.
