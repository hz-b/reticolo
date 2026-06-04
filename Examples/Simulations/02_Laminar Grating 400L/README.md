# Laminar Grating 400 l/mm

Physical 400 l/mm laminar (trapezoidal) grating measured at HZB/ELISA. Compare against `../../Reference_Data/HORIBA_ELISA_Experimental/`.

Geometry: depth 14.9 nm, sidewall 15°, groove/period 0.67. Stack: 29 nm Pt on Si. Uncomment the `add_layer` call for C to add a 1 nm contamination layer.

The runner does an energy sweep (50–1000 eV) at fixed C_ff = 2.25, TM, order −1. To switch to a fixed angle instead, comment out `sweep.Cff` and uncomment `sweep.alpha_deg`.

Result sets in `Results/`:
- `TE_4deg` — TE, α = 4°
- `TM_4deg` — TM, α = 4°
- `TM_4deg_C_top_layer` — TM, α = 4°, contaminated
- `TM_2.25cff` — TM, C_ff = 2.25
- `TM_2.25cff_C_top_layer` — TM, C_ff = 2.25, contaminated
