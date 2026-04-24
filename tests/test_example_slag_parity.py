from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from reticolo_py.slag import SlagConfig, run_example_slag


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.skipif(shutil.which("octave") is None, reason="Octave is required for parity checks")
def test_example_slag_matches_octave_reference(tmp_path: Path) -> None:
    output_path = tmp_path / "octave_reference.csv"
    energies = [140.0, 150.0, 160.0]

    command = [
        "octave",
        "-qf",
        str(REPO_ROOT / "tests" / "octave_example_slag_reference.m"),
        str(output_path),
        *[str(value) for value in energies],
    ]
    subprocess.run(command, cwd=REPO_ROOT, check=True, capture_output=True, text=True)

    octave_rows = []
    with output_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            octave_rows.append((float(row["energy_ev"]), float(row["efficiency"])))

    config = SlagConfig(photon_energy_ev=np.asarray(energies, dtype=float), base_dir=str(REPO_ROOT))
    result = run_example_slag(config=config, save_plot=False)
    python_rows = list(zip(result["energy_ev"], result["efficiency"]))

    assert len(python_rows) == len(octave_rows)
    for (python_energy, python_efficiency), (octave_energy, octave_efficiency) in zip(python_rows, octave_rows):
        assert python_energy == pytest.approx(octave_energy)
        assert python_efficiency == pytest.approx(octave_efficiency, abs=5e-3)
