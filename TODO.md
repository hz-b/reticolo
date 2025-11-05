# 📝 TODO — RETICOLO Grating Examples

This document tracks pending tasks and improvements for the Octave/MATLAB examples that interface with **RETICOLO**.

## 🧩 General structure and coding guidelines

All reusable helper routines must be moved into individual `.m` files located in the `helpers/` folder.

### Guidelines
- Each function should be stored in **its own file**.  
  Example:  
  ```text
  helpers/
  ├── estimateTheta.m
  ├── efficiency_bgrML.m
  └── formatElapsedTime.m
  ```
- The filename **must match the function name** (case-sensitive), e.g.:
  - Function: `function thetaEst = estimateTheta(...)`
  - File: `helpers/estimateTheta.m`
- Avoid defining multiple functions in a single `.m` file (except for nested local helpers explicitly needed inside a single algorithm).
- Use clear, descriptive function names that start with lowercase letters (e.g., `computeProfile`, `loadMaterialData`).
- All examples (`Example_MLBG.m`, `Example_SLAG.m`, `Example_SLBG.m`) should load helpers via:
  ```matlab
  addpath(genpath(fullfile(pwd, 'helpers')))
  ```
- When a new helper is created, document its purpose briefly in a comment block at the top of the file.

---

---

## ✅ Example 1 — `Example_MLBG.m` (Multilayer Blazed Grating)

### Current behavior
The current version of `Example_MLBG.m` computes the **diffraction efficiency** of a **multilayer blazed grating (MLBG)** over a **range of photon energies**, while keeping the **incident (grazing) angle** fixed.  
Each simulation step stores efficiency results and angular data in structured outputs (`outputdata`, `outputdata0`, etc.).

### Goal
Modify the workflow so that for **each photon energy**, the **incident angle** is **automatically computed** using the helper function `estimateTheta()` from the `helpers/` folder.  
This ensures that the Bragg condition is satisfied dynamically across the energy range.

### Desired behavior
- The example should loop over a **list of photon energies** (e.g., `photonEnergy_eV = [1000, 1500, 2000, 2500, ...]`).
- For each energy value:
  1. Call `thetaEst = estimateTheta(...)` to compute the corresponding incidence angle.
  2. Run the RETICOLO efficiency calculation (`efficiency_bgrML`) using that `thetaEst`.
- The result will be an **energy-dependent efficiency curve**, where the incidence angle adapts to the multilayer’s Bragg condition.

### Notes
- `estimateTheta()` uses optical constants of the selected high-Z and low-Z materials to estimate the effective Bragg angle for a given energy.
- This approach is expected to produce smoother and more physically meaningful energy-dependent efficiency spectra.
- **Validation:** Computed incidence angles and resulting efficiency profiles should be **checked with Andrey**.

### Tasks
- [ ] For each energy, compute `thetaEst` using `estimateTheta()`.
- [ ] Pass both parameters to `efficiency_bgrML()` inside the loop.
- [ ] Plot efficiency as a function of photon energy.
- [ ] Validate and cross-check results with Andrey’s data.
- [ ] write a function that does the job, we will use it to call it from python. 

---

## ⏳ Example 2 — `Example_SLAG.m` (Single-Layer Laminar Grating)

### Objective
Simulate a **single-layer laminar grating** (SLAG) and compute its **diffraction efficiency** over a range of photon energies or incident angles.

### Validation
The results must be **benchmarked against an independent reference code**, such as **Reflec** or another validated diffraction solver.  
**Ask Andrey** which tool he prefers to use for cross-validation.

### Future development
This script should be **refactored into a callable function**, so that it can later be invoked directly from **Python**, for example via:
```python
subprocess.run(["octave", "--eval", "Example_SLAG"])
```
or through a dedicated Octave-Python bridge such as **oct2py**.

### Tasks
- [ ] Implement SLAG efficiency simulation logic.
- [ ] Validate results with Reflec (or another code chosen by Andrey).
- [ ] Convert the script into a callable function interface for Python integration.

---

## ⏳ Example 3 — `Example_SLBG.m` (Single-Layer Blazed Grating)

### Objective
Simulate a **single-layer blazed grating** (SLBG) and compute the **diffraction efficiency** under varying photon energies and blaze geometries.

### Validation
As with `Example_SLAG.m`, validate all computed diffraction efficiencies against **Reflec** or a similar tool.  
Confirm validation methodology and reference results with **Andrey**.

### Future development
This example should also evolve into a **function-based script** for integration into a future **Python workflow** that automates all three simulations.

### Tasks
- [ ] Implement SLBG simulation setup using RETICOLO.
- [ ] Validate with Reflec or Andrey’s chosen comparison tool.
- [ ] Refactor into a function callable from Python.

---
