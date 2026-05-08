from ast import Load

import pandas as pd
import matplotlib.pyplot as plt

# List your CSV files here
files = [
    "rcwa_trap_dc0.67_400lmm_d15nm_Si_SiO21nm_Cr6nm_Pt29nm_C1nm_TM_order-1_alpha86.00deg_50-1000eV.csv"
]

plt.figure()


labels = ["new constants", "old constants"]  # Adjust labels as needed
for f in files:
    # Load CSV (adjust if you have headers)
    data = pd.read_csv(f)

    # Assume first column = x, second = y
    x = data["PhotonEnergy_eV"]  # or data.iloc[:, 0] if no headers
    y = data["DiffractionEfficiency"]  # or data.iloc[:, 1] if no headers
    

    plt.plot(x, y, label=labels[files.index(f)],linewidth=.5)  # Use label if you want a legend






# RetiPy= pd.read_csv("ELISA_400lmm_top_C.csv")
# plt.plot(RetiPy["energy_ev"], RetiPy["efficiency_order_-1"], label='RETICOLO-PYTHON',linewidth=.5)





# Load experimental data
experiment = pd.read_csv(
    "Re__ELISA,_400l_mm_laminar_grating_from_HORIBA/lG400-HZB-ELISA_ascan-energy_alpha-4deg_1-order.csv",
    sep=";",
    decimal=",",
    skiprows=3,        
    header=None
)


experiment = experiment.dropna(axis=1, how='all')

experiment.columns = ["PhotonEnergy_eV", "DiffractionEfficiency"]

plt.plot(experiment["PhotonEnergy_eV"], experiment["DiffractionEfficiency"], label="experiment",linewidth=.5)








# # Load reflec file
# reflec = pd.read_csv(
#     "Reflec_simulations/Simulations_REFLEC_REFLEC(SPECS).txt",
#     sep="\t",
#     skiprows=3,     # skip header rows
#     header=None
# )

# # Drop completely empty columns (caused by uneven headers)
# reflec = reflec.dropna(axis=1, how='all')

# # Extract columns (0-based indexing)
# x = reflec.iloc[:, 0]   # Energy
# y = reflec.iloc[:, 3]   # 4th column

# # Plot REFLEC simulation
# plt.plot(x, y, label='REFLEC',linewidth=.5)







plt.xlabel("Photon Energy (eV)")
plt.ylabel("Diffraction Efficiency")
plt.title(" Diffraction Efficiency vs Photon Energy")
plt.legend()
plt.grid()
plt.savefig("csv_plots.png", dpi=300, bbox_inches="tight")
plt.show()
