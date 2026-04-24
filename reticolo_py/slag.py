from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .materials import interpolate_cxro_index
from .rcwa_1d import res0, res1, res2


@dataclass
class SlagConfig:
    grating_period_lpermm: int = 400
    diffraction_order: int = 1
    width_to_period_ratio: float = 0.67
    depth_nm: float = 14.9
    trapezoid_left_deg: float = 15.0
    trapezoid_right_deg: float = 15.0
    substrate_material: str = "Si"
    layer_material: str = "Pt"
    layer_thickness_nm: float = 28.77
    z_resolution_nm: float = 0.1
    x_resolution_nm: float = 0.1
    fourier_orders: int = 25
    photon_energy_ev: np.ndarray = field(default_factory=lambda: np.arange(140.0, 165.0, 5.0))
    grazing_angle_deg: float = 4.0
    base_dir: str = "."
    experimental_csv: str = (
        "Re__ELISA,_400l_mm_laminar_grating_from_HORIBA/"
        "lG400-HZB-ELISA_ascan-energy_alpha-4deg_1-order.csv"
    )


def default_example_slag_config() -> SlagConfig:
    return SlagConfig()


def run_example_slag(
    config: SlagConfig | None = None,
    save_plot: bool = True,
    plot_filename: str = "Example_SLAG_efficiency_4deg.png",
) -> dict[str, np.ndarray]:
    config = config or default_example_slag_config()
    energies = np.asarray(config.photon_energy_ev, dtype=float)
    efficiencies = []
    diffraction_angles = []

    for energy in energies:
        simulation = simulate_single_energy(config, float(energy))
        efficiencies.append(simulation["efficiency"])
        diffraction_angles.append(simulation["diffraction_angle_deg"])

    efficiencies = np.asarray(efficiencies, dtype=float)
    diffraction_angles = np.asarray(diffraction_angles, dtype=float)

    result = {
        "energy_ev": energies,
        "efficiency": efficiencies,
        "diffraction_angle_deg": diffraction_angles,
    }

    experimental = load_experimental_csv(Path(config.base_dir) / config.experimental_csv)
    if save_plot:
        plot_example_slag(
            result["energy_ev"],
            result["efficiency"],
            experimental[:, 0],
            experimental[:, 1],
            plot_filename,
        )

    result["experimental_energy_ev"] = experimental[:, 0]
    result["experimental_efficiency"] = experimental[:, 1]
    return result


def simulate_single_energy(config: SlagConfig, photon_energy_ev: float) -> dict[str, float]:
    wavelength_nm = 1239.8 / photon_energy_ev
    n_sub = interpolate_cxro_index(config.substrate_material, photon_energy_ev, base_dir=config.base_dir)
    n_layer = interpolate_cxro_index(config.layer_material, photon_energy_ev, base_dir=config.base_dir)
    period_nm = 1e6 / config.grating_period_lpermm
    k_parallel = np.sin(np.deg2rad(90.0 - config.grazing_angle_deg))

    textures, profile = build_slag_textures(config, n_sub=n_sub, n_layer=n_layer, n_inc=1.0 + 0.0j, period_nm=period_nm)
    parm = res0(1)
    aa = res1(
        wavelength_nm,
        period_nm,
        textures,
        config.fourier_orders,
        k_parallel,
        parm,
    )
    ef = res2(aa, profile, parm)
    order_index = np.where(ef.inc_top_reflected.order == -config.diffraction_order)[0]
    if len(order_index) != 1:
        raise ValueError(f"Unable to locate diffraction order {config.diffraction_order}")
    idx = int(order_index[0])
    return {
        "efficiency": float(np.real_if_close(ef.inc_top_reflected.efficiency[idx])),
        "diffraction_angle_deg": float(90.0 - ef.inc_top_reflected.theta[idx]),
    }


def build_slag_textures(
    config: SlagConfig,
    *,
    n_sub: complex,
    n_layer: complex,
    n_inc: complex,
    period_nm: float,
) -> tuple[list[object], tuple[np.ndarray, np.ndarray]]:
    thickness_nm = config.depth_nm + config.layer_thickness_nm + 5.0

    if config.trapezoid_left_deg == 0 and config.trapezoid_right_deg == 0:
        edge_z = config.depth_nm
        edge_x = config.width_to_period_ratio * period_nm
        positions = np.array(
            [0.0, period_nm / 2 - edge_x / 2, period_nm / 2 + edge_x / 2, period_nm],
            dtype=float,
        )
        heights = np.array([0.0, 0.0, edge_z, 0.0], dtype=float)
    else:
        pos1 = [
            (period_nm - config.width_to_period_ratio * period_nm) / 2
            - config.depth_nm * np.tan(np.deg2rad(config.trapezoid_left_deg)),
            0.0,
        ]
        pos2 = [(period_nm - config.width_to_period_ratio * period_nm) / 2, config.depth_nm]
        pos3 = [(period_nm + config.width_to_period_ratio * period_nm) / 2, config.depth_nm]
        pos4 = [
            (period_nm + config.width_to_period_ratio * period_nm) / 2
            + config.depth_nm * np.tan(np.deg2rad(config.trapezoid_right_deg)),
            0.0,
        ]
        profile_points = np.array(
            [[0.0, 0.0], pos1, pos2, pos3, pos4, [period_nm, 0.0]],
            dtype=float,
        )
        positions = profile_points[:, 0]
        heights = profile_points[:, 1]

    x = np.linspace(0.0, period_nm, int(round(period_nm / config.x_resolution_nm)) + 1)
    z = np.linspace(thickness_nm, 0.0, int(round(thickness_nm / config.z_resolution_nm)) + 1)
    x_grid, z_grid = np.meshgrid(x, z)
    surface = np.interp(x, positions, heights)
    coating_top = surface + config.layer_thickness_nm

    index_grid = np.full_like(x_grid, n_sub, dtype=complex)
    index_grid[z_grid >= surface[None, :]] = n_inc
    layer_mask = (z_grid >= surface[None, :]) & (z_grid < coating_top[None, :])
    index_grid[layer_mask] = n_layer

    textures: list[object] = [n_inc]
    layer_texture_indices: list[int] = []

    previous_texture: object | None = None
    current_texture_index = -1
    thicknesses: list[float] = [0.0]
    texture_sequence: list[int] = [0]

    for row in range(index_grid.shape[0]):
        row_values = index_grid[row]
        row_changes = np.nonzero(np.diff(row_values) != 0)[0]
        if len(row_changes) == 0:
            texture: object = complex(row_values[0])
        else:
            texture = [x[row_changes + 1], row_values[row_changes]]

        if _textures_equal(previous_texture, texture):
            thicknesses[-1] += config.z_resolution_nm
        else:
            textures.append(texture)
            current_texture_index = len(textures) - 1
            texture_sequence.append(current_texture_index)
            thicknesses.append(config.z_resolution_nm)
            previous_texture = texture
        layer_texture_indices.append(current_texture_index)

    textures.append(n_sub)
    texture_sequence.append(len(textures) - 1)
    thicknesses.append(0.0)
    return textures, (np.asarray(thicknesses, dtype=float), np.asarray(texture_sequence, dtype=int))


def load_experimental_csv(path: str | Path) -> np.ndarray:
    path = Path(path)
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for _ in range(3):
            next(handle)
        for line in handle:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            energy_text, efficiency_text, *_ = line.split(";")
            rows.append(
                [
                    float(energy_text.replace(",", ".")),
                    float(efficiency_text.replace(",", ".")),
                ]
            )
    return np.asarray(rows, dtype=float)


def plot_example_slag(
    energy: np.ndarray,
    efficiency: np.ndarray,
    experimental_energy: np.ndarray,
    experimental_efficiency: np.ndarray,
    output_filename: str,
) -> None:
    figure, axis = plt.subplots(figsize=(10, 7))
    axis.plot(energy, efficiency, "b-o", linewidth=0.5, markersize=2.0, label="Simulation")
    axis.plot(
        experimental_energy,
        experimental_efficiency,
        "r-s",
        linewidth=0.5,
        markersize=2.0,
        label="Experimental Data",
    )
    axis.set_xlabel("Photon Energy (eV)")
    axis.set_ylabel("Diffraction Efficiency")
    axis.set_title("SLAG Simulation vs Experimental Data (400 l/mm Grating, alpha = 4 deg)")
    axis.grid(True, alpha=0.3)
    axis.legend(loc="best")
    figure.tight_layout()
    figure.savefig(output_filename, dpi=150, bbox_inches="tight")
    plt.close(figure)


def _textures_equal(left: object | None, right: object | None) -> bool:
    if left is None or right is None:
        return False
    if isinstance(left, complex) or isinstance(left, float) or isinstance(left, int):
        if isinstance(right, complex) or isinstance(right, float) or isinstance(right, int):
            return complex(left) == complex(right)
        return False
    if not isinstance(left, list) or not isinstance(right, list):
        return False
    return np.array_equal(np.asarray(left[0]), np.asarray(right[0])) and np.array_equal(
        np.asarray(left[1]), np.asarray(right[1])
    )
