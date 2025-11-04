# RETICOLO-Based Grating Efficiency Examples

This folder contains Octave/MATLAB examples demonstrating how to simulate **diffraction efficiency** of multilayer-coated gratings using the [RETICOLO](https://www.lp2n.institutoptique.fr/equipes/RETICOLO) electromagnetic solver.  
Each example builds the required structure (geometry, materials, and optical parameters), then calls a helper routine that interfaces with RETICOLO to compute diffraction efficiencies and optionally visualize field distributions.

---

## 📁 Folder structure

```
reticolo/
├── Example_MLBG.m               # Multilayer Blazed Grating example
├── Example_SLBG.m               # Single-Layer Blazed Grating example
├── Example_SLAG.m               # Single-Layer Laminar Grating example
└── helpers/
    ├── estimateTheta.m          # Computes approximate Bragg/grazing incidence angle
    └── efficiency_bgrML.m       # Builds structure, runs RETICOLO, returns efficiencies
```

---

## ⚙️ Common dependencies

All examples share the same helper functions and require:

- Material refractive index files named as:
  ```
  n_Si_cxro.txt
  n_Cr_cxro.txt
  n_C_cxro.txt
  n_Au_cxro.txt
  ...
  ```
  Each text file should contain three columns:
  | Photon Energy [eV] | Re(1-n) | Im(1-n) |

The examples automatically add the `helpers/` folder to the MATLAB/Octave path at runtime.

---

## 🧪 Example 1 — `Example_MLBG.m`: Multilayer Blazed Grating

### Overview
This example simulates a **multilayer-coated blazed grating (MLBG)** with a substrate and two materials alternating in the multilayer stack.  
The multilayer period, number of bilayers, blaze angles, and photon energy are user-configurable.

### Workflow
1. **Setup of parameters**
   - Define grating geometry (`grPeriod_lpermm`, `grBA_deg`, `grAntiBA_deg`)
   - Define materials for substrate, high-Z layer, and low-Z layer
   - Define multilayer stack: period (`ML_d_nm`), ratio (`ML_d_HZtod`), and repetitions (`ML_N`)
   - Define simulation parameters: grid resolution, Fourier orders, and polarization

2. **Estimate grazing incidence**
   The helper `estimateTheta.m` computes the expected Bragg condition for multilayer reflection combined with the diffraction equation for the grating.

3. **Run RETICOLO efficiency calculation**
   The helper `efficiency_bgrML.m`:
   - Builds the multilayered, blazed surface geometry
   - Converts it to RETICOLO texture and profile arrays
   - Calls the RETICOLO routines (`res0`, `res1`, `res2`, optionally `res3`)
   - Returns the diffraction efficiency in the top reflected orders

4. **Plot results**
   The diffraction efficiency vs. incidence angle is plotted, showing the peak efficiency for the target diffraction order.

---

## 🧩 Helper functions

### `estimateTheta.m`
Computes the estimated grazing angle by combining Bragg’s law and the grating equation.  
Inputs: materials, multilayer period, order, groove density, photon energy.  
Output: estimated grazing incidence angle in degrees.

### `efficiency_bgrML.m`
Builds the 2D multilayer grating profile, computes its optical response with RETICOLO, and optionally plots the field distribution.  
Inputs: geometry, materials, multilayer parameters, photon energy, polarization, and RETICOLO settings.  
Output: RETICOLO structure `ef` containing efficiencies and angles for each diffracted order.

---

## 🧱 Example 2 — `Example_SLBG.m` (placeholder)
Single-layer **blazed grating** with uniform coating on the grooves.  
To be added.

---

## 🧱 Example 3 — `Example_SLAG.m` (placeholder)
Single-layer **laminar grating** (rectangular profile) with coating.  
To be added.

---

## 🧭 Usage

Run from Octave or MATLAB:
```
>> Example_MLBG
```
This will load the helper functions, compute the diffraction efficiency for a multilayer blazed grating, and produce a figure of `efficiency vs. grazing angle`.

---
