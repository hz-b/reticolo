import pandas as pd
import matplotlib.pyplot as plt

# List your CSV files here
files = [
    "rcwa_blazed_b0.73deg_ab5.60deg_600lmm_d19nm_Si_Au31nm_C1nm_TM_order-1_Cff2.25_50-2000eV.csv",
    "rcwa_blazed_b0.73deg_ab5.60deg_600lmm_d19nm_Si_Au31nm_TM_order-1_Cff2.25_50-2000eV.csv"
]

plt.figure()
labels = ["carbon coating", "bare gold"]
for f in files:
    # Load CSV (adjust if you have headers)
    data = pd.read_csv(f)

    # Assume first column = x, second = y
    x = data["PhotonEnergy_eV"]  # or data.iloc[:, 0] if no headers
    y = data["DiffractionEfficiency"]  # or data.iloc[:, 1] if no headers
    

    plt.plot(x, y, label=labels[files.index(f)])

plt.xlabel("Photon Energy (eV)")
plt.ylabel("Diffraction Efficiency")
plt.title("CSV plots")
plt.legend()
plt.grid()
plt.savefig("csv_plots.png", dpi=300, bbox_inches="tight")
plt.show()
