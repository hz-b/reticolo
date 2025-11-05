# RETICOLO-Based Grating Efficiency Examples

This folder contains Octave/MATLAB examples demonstrating how to simulate **diffraction efficiency** of multilayer-coated gratings using the [RETICOLO](https://www.lp2n.institutoptique.fr/equipes/RETICOLO) electromagnetic solver.  
Each example builds the required structure (geometry, materials, and optical parameters), then calls a helper routine that interfaces with RETICOLO to compute diffraction efficiencies and optionally visualize field distributions.

---

### 🐧 Linux Installation

#### Ubuntu / Debian
You can install Octave from the official repositories:

```bash
sudo apt update
sudo apt install octave octave-common octave-geometry octave-control
```

Optional: to install the graphical user interface (GUI) and plotting dependencies:

```bash
sudo apt install gnuplot-x11 liboctave-dev
```

Then launch Octave:

```bash
octave --gui
```

#### ❗ Known Issue

When running Octave from the terminal or from Python, you may encounter the following error:

```
/usr/libexec/octave/8.4.0/exec/x86_64-pc-linux-gnu/octave-gui: symbol lookup error: /snap/core20/current/lib/x86_64-linux-gnu/libpthread.so.0: undefined symbol: __libc_pthread_init, version GLIBC_PRIVATE
```

This happens when Octave (often the **Snap** version) links against the wrong `libpthread.so.0` library provided by `/snap/core20/...`, instead of the system library from `/lib/x86_64-linux-gnu`.



##### ✅ Permanent Solution (System-Wide)

To make the fix permanent for all users, modify `/etc/environment`:

```bash
sudo nano /etc/environment
```

and then add the following line

```bash
LD_LIBRARY_PATH="/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu"
```

Then log out and back in (or reboot) to apply the changes.


---

## 📁 Folder structure

```
reticolo/
├── V9                           # Reticolo V9
├── V7-reticolo-blazr            # Reticolo V7 (with blazed grating support?)
├── ZYQ_files                    # examples from ZYQ
├── Example_MLBG.m               # Multilayer Blazed Grating example
├── Example_SLBG.m               # Single-Layer Blazed Grating example
├── Example_SLAG.m               # Single-Layer Laminar Grating example
└── helpers/
    ├── estimateTheta.m          # Computes approximate Bragg/grazing incidence angle
    └── efficiency_bgrML.m       # Builds structure, runs RETICOLO, returns efficiencies
```

---

## Materials refractive index files

The material refractive index files required by the examples must be in the CXRO format. To download file for new material visist the [CXRO - INdex of Refraction](https://henke.lbl.gov/optical_constants/getdb2.html).

- Material refractive index files named as:
  ```
  n_Si_cxro.txt
  n_Cr_cxro.txt
  n_C_cxro.txt
  n_Au_cxro.txt
  ...
  ```
  Each text file should contain three columns:
  | Photon Energy [eV] | Delta | Beta |


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


## 🧱 Example 1 — `Example_MLBG.m` 
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


## 🧱 Example 2 — `Example_SLBG.m`
Single-layer **blazed grating** with uniform coating on the grooves.  


---

## 🧱 Example 3 — `Example_SLAG.m` 
Single-layer **laminar grating** (rectangular profile) with coating.  

---


