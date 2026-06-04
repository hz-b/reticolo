# Blazed Multilayer Grating 2400 l/mm

Physical 2400 l/mm blazed grating with a Cr/C multilayer coating for hard x-ray Bragg reflection. Compare against `../../Reference_Data/`.

Geometry: blaze 1.37°, anti-blaze 3.25°. Stack: 60× Cr(1.9 nm)/C(2.9 nm) on Si.

The runner sweeps energy along the Bragg condition using `bragg_lookup_table.csv`, which maps each energy to the matching α. Range: 3000–5000 eV, TM, order −2. Update the table with `load_bragg_table` from `helpers/` if you change the multilayer period.
