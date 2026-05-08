# SLAG Simulation Workflow

This workflow allows you to run SLAG simulations in Octave and visualize the results using Python.

## Files

### Octave/MATLAB Scripts

1. **Example_SLAG_save_data.m** - Main simulation script that:
   - Runs the SLAG simulation with configurable parameters
   - Saves results to a CSV file with embedded metadata
   - Outputs summary statistics

2. **Example_SLAG.m** - Original simulation script (reference)

### Python Scripts

1. **plot_simulation_data.py** - Visualization script that:
   - Reads CSV files with embedded metadata
   - Displays metadata in a formatted way
   - Creates efficiency vs energy plots
   - Saves plots with timestamps

## Usage

### Step 1: Run the Octave Simulation

```bash
octave Example_SLAG_save_data.m
```

This will:
- Run the simulation with the parameters defined in the script
- Create a CSV file named like: `SLAG_simulation_400lmm_alpha4.0deg_order1_Pt_Au.csv`
- Display simulation summary

### Step 2: Modify Parameters (Optional)

Edit `Example_SLAG_save_data.m` to change simulation parameters:

```matlab
metadata.grPeriod_lpermm = 400;        % grating period (lines/mm)
metadata.grDepth_nm = 14.9;            % grating depth (nm)
metadata.grazing_angle_deg = 4;        % grazing incidence angle (degrees)
metadata.photonEnergy_min = 100;       % minimum photon energy (eV)
metadata.photonEnergy_max = 2000;      % maximum photon energy (eV)
metadata.photonEnergy_step = 10;       % energy step size (eV)
```

### Step 3: Plot with Python

```bash
python3 plot_simulation_data.py
```

Or specify a specific file:

```bash
python3 plot_simulation_data.py SLAG_simulation_400lmm_alpha4.0deg_order1_Pt_Au.csv
```

This will:
- Parse the CSV file and extract metadata
- Display all simulation parameters
- Create a plot with efficiency vs energy
- Save the plot as `SLAG_plot_YYYYMMDD_HHMMSS.png`

## Output Files

### CSV Data File

The simulation creates a CSV file with:
- Metadata as comments (lines starting with #)
- Data columns: `photon_energy_eV;diffraction_efficiency`
- European number format support

Example structure:
```
# SLAG Simulation Data
# grating_period_lpermm: 400
# grazing_angle_deg: 4
# photon_energy_min_eV: 100
# photon_energy_max_eV: 2000
#
# photon_energy_eV;diffraction_efficiency
51.031;0.03528338
51.535;0.03446764
...
```

### Plot Files

The Python script creates PNG plots with:
- Main plot: Efficiency vs Photon Energy
- Distribution plot: Efficiency histogram
- Statistics: Max, Min, Mean efficiency
- Metadata in title

## Requirements

### Octave/MATLAB
- GNU Octave 4.0+ or MATLAB
- RETICOLO package (in V9/reticolo_allege_v9/)

### Python
- Python 3.6+
- matplotlib
- numpy

Install Python dependencies:
```bash
pip install matplotlib numpy
```

## Customization

### Changing Simulation Parameters

Edit the following sections in `Example_SLAG_save_data.m`:

**Grating Structure:**
```matlab
metadata.grPeriod_lpermm = 400;
metadata.grDepth_nm = 14.9;
metadata.grWidthtoD = 0.67;
```

**Incidence Conditions:**
```matlab
grazing_angle_deg = 4;
metadata.photonEnergy_eV = 100:10:2000;
```

**Materials:**
```matlab
metadata.material_sub = 'Pt';
metadata.material_layer = 'Au';
```

### Customizing Plots

Edit `plot_simulation_data.py` to:
- Change plot styles
- Add additional analysis
- Modify output formats (PDF, SVG)
- Add comparison with experimental data

## Workflow Diagram

```
Example_SLAG_save_data.m
         ↓
    [Octave]
         ↓
SLAG_simulation_*.csv (with metadata)
         ↓
plot_simulation_data.py
         ↓
    [Python]
         ↓
SLAG_plot_*.png (visualization)
```

## Troubleshooting

### Octave Errors
- Ensure RETICOLO path is correct: `V9/reticolo_allege_v9/`
- Check that material index files exist (n_Pt_cxro.txt, n_Au_cxro.txt, etc.)

### Python Errors
- Install dependencies: `pip install matplotlib numpy`
- Check that CSV file exists and is readable
- Ensure file has proper metadata format

## Advanced Usage

### Batch Processing

To process multiple simulation files:

```bash
for file in SLAG_simulation_*.csv; do
    python3 plot_simulation_data.py "$file"
done
```

### Comparing with Experimental Data

1. Run simulation: `octave Example_SLAG_save_data.m`
2. Load experimental data in Python
3. Modify `plot_simulation_data.py` to overlay experimental data
4. Calculate comparison metrics (MAE, RMSE)

## Support

For issues or questions, check:
- README.md (main project documentation)
- TODO.md (planned improvements)
- Original Example_SLAG.m (reference implementation)
