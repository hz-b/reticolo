from __future__ import annotations

from pathlib import Path

import numpy as np


def load_cxro_data(material: str, base_dir: str | Path = ".") -> np.ndarray:
    """Load legacy `n_<material>_cxro.txt` data used by Example_SLAG.m."""

    path = Path(base_dir) / f"n_{material}_cxro.txt"
    if not path.exists():
        raise FileNotFoundError(f"CXRO file not found: {path}")

    data = np.loadtxt(path, skiprows=2)
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(f"Unexpected CXRO format in {path}")
    return data[:, :3]


def interpolate_cxro_index(
    material: str,
    photon_energy_ev: float,
    base_dir: str | Path = ".",
) -> complex:
    """Return the complex refractive index used by the MATLAB SLAG example."""

    data = load_cxro_data(material, base_dir=base_dir)
    delta = np.interp(photon_energy_ev, data[:, 0], data[:, 1])
    beta = np.interp(photon_energy_ev, data[:, 0], data[:, 2])
    return 1.0 - delta + 1j * beta
