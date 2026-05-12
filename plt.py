import numpy as np
import matplotlib.pyplot as plt



eff = [45.54,  47.95,  50.22,  52.37,  54.31,  55.24,  57.72,
       59.33,57.89,61.49,62.98,64.00,58.39,64.91,
       66.29,67.04,67.45,65.18,66.13,68.04,68.35,
       68.21,67.80,67.14,65.51,60.35,63.05,61.42,
       58.53,55.57,51.61,46.83,41.05,34.85,29.47,
       21.64,16.09,3.87,1.86,3.61,1.33]

ev = np.linspace(2000, 6000, len(eff))

for e, eV in zip(eff, ev):
    print(f'{eV},{e}')
plt.figure(figsize=(8, 6))
plt.plot(ev, eff, marker='o')
plt.xlabel('Diffraction Efficiency (%)')
plt.ylabel('Photon Energy (eV)')
plt.title('Diffraction Efficiency vs Photon Energy')
plt.grid()
plt.show()
plt.savefig('diffraction_efficiency.png')   