# Helpers

Functions for building grating stacks, running RCWA sweeps, and plotting results. Add this folder to the path; runner scripts do this automatically.

---

**`build_grating(type, period_nm, depth_nm, ...)`**  
Creates a grating geometry descriptor. Two profile types:
- `'blazed'` — triangular groove; extra args: `blaze_deg, antiblaze_deg, x_res_nm`
- `'trapezoidal'` — flat-bottomed groove; extra args: `width_ratio, angle_deg, x_res_nm` where `width_ratio` is groove floor / period and `angle_deg` is sidewall angle from horizontal

**`build_stack(grating)`**  
Initialises an empty layer stack from a grating descriptor.

**`add_layer(stack, material_file, thickness_nm)`**  
Appends a conformal coating layer above the existing stack. `material_file` is a CXRO `.txt` file (columns: Energy eV, Delta, Beta; non-numeric header lines are skipped). Layers follow the grating surface shape.

**`run_rcwa(stack, substrate_file, sweep, options)`**  
Main solver. Builds the meshgrid, calls RETICOLO at each sweep point, and writes `simulation_results.csv` to `options.output_dir`.

Sweep modes controlled by `sweep.type`:
- `'energy'` with `sweep.alpha_deg` — energy sweep at fixed grazing angle
- `'energy'` with `sweep.Cff` — energy sweep at fixed C_ff; α recomputed from the grating equation at each step
- `'alpha'` with `sweep.energy_eV` — angle sweep at fixed energy
- `'bragg'` — follows a (energy, α) lookup table; use `load_bragg_table` to build the sweep struct

Key options: `FourierOrders` (default 11), `pol` (+1 TE / −1 TM), `GR_Order`, `z_res_nm`, `reticolo_path`, `output_dir`, `oc_path`.

**`load_bragg_table(filepath, energy_range)`**  
Loads a CSV/TSV with `Energy` and `alpha` columns and returns a sweep struct for `run_rcwa`. If `energy_range` is supplied, interpolates onto that grid; errors if any requested energies fall outside the table.

**`plot_meshgrid(stack, substrate_file, photon_eV, z_res_nm, save_png)`**  
Renders the Im(n) cross-section of the grating stack at a given photon energy and saves a PNG. Shows exactly the refractive-index map the solver receives.

**`plot_results(stack, substrate_file, sweep, options, results)`**  
Convenience wrapper: saves a meshgrid PNG at the mid-sweep energy and an efficiency-vs-sweep-variable plot. Called at the end of most runner scripts.
