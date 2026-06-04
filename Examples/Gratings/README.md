# Grating Definitions

Each script builds a stack and saves a cross-section meshgrid PNG to `Results/`. No RCWA sweep is run. Change `photon_eV` near the bottom of any script to render at a different energy.

`Single_Layer_Blazed_Grating` — 600 l/mm, blaze 0.73° / anti-blaze 5.60°, 31 nm Au on Si.

`Single_Layer_Laminar_Grating` — 400 l/mm trapezoidal, depth 14.9 nm, sidewall 15°, groove/period 0.67, Pt + 1 nm C contamination on Si.

`Multilayer_Blazed_Grating` — 2400 l/mm, blaze 1.37° / anti-blaze 3.25°, 60× Cr(1.9 nm)/C(2.9 nm) on Si. Use `x_res_nm = 0.1` or finer to resolve individual multilayer periods.
