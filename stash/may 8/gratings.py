"""Grating profile objects used by RCWA simulations."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .materials import interpolate_cxro_index
from .stacks import BaseStack, MultilayerStack, SingleLayerStack


@dataclass
class BaseGrating(ABC):
    """Base class for grating profiles and material stacks."""

    period_lpermm: int = 400
    width_to_period_ratio: float = 0.67
    depth_nm: float = 14.9
    coating_stack: BaseStack | None = None
    substrate_material: str = "Si"
    layer_material: str = "Pt"
    layer_thickness_nm: float = 28.77
    top_cap_material: str | None = None
    top_cap_thickness_nm: float = 0.0
    z_resolution_nm: float = 0.1
    x_resolution_nm: float = 1.0
    base_dir: str = "."

    @property
    def period_nm(self) -> float:
        """Return the grating period in nanometers."""

        return 1e6 / self.period_lpermm

    @abstractmethod
    def profile_points(self) -> tuple[np.ndarray, np.ndarray]:
        """Return one grating period as x and height arrays."""

    def _tiled_profile_points(self, num_periods: int = 3) -> tuple[np.ndarray, np.ndarray]:
        """Return profile points repeated over multiple periods.

        Args:
            num_periods: Number of periods to tile.

        Returns:
            Tiled x positions and heights.
        """

        if num_periods < 1:
            raise ValueError("num_periods must be at least 1.")

        base_positions, base_heights = self.profile_points()
        tiled_positions = []
        tiled_heights = []
        for period_index in range(num_periods):
            offset = period_index * self.period_nm
            if period_index == 0:
                tiled_positions.extend((base_positions + offset).tolist())
                tiled_heights.extend(base_heights.tolist())
            else:
                tiled_positions.extend((base_positions[1:] + offset).tolist())
                tiled_heights.extend(base_heights[1:].tolist())
        return np.asarray(tiled_positions, dtype=float), np.asarray(tiled_heights, dtype=float)

    def resolved_stack(self) -> BaseStack:
        """Return the common coating stack used by the grating."""

        if self.coating_stack is not None:
            return self.coating_stack
        return SingleLayerStack(
            substrate_material=self.substrate_material,
            layer_material=self.layer_material,
            layer_thickness_nm=self.layer_thickness_nm,
            top_cap_material=self.top_cap_material,
            top_cap_thickness_nm=self.top_cap_thickness_nm,
        )

    def _material_plot_data(
        self,
        num_periods: int = 1,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
        """Return discretized material data for profile visualization.

        Args:
            num_periods: Number of periods to include.

        Returns:
            X grid, z grid, material-code map, and material labels.
        """

        coating_stack = self.resolved_stack()
        x_grid = self._build_x_grid(num_periods=num_periods)
        z_grid = self._build_plot_z_grid(coating_stack)
        surface = self._surface_profile_on_grid(x_grid, num_periods=num_periods)
        material_labels = coating_stack.plot_material_names()
        material_codes = {label: index for index, label in enumerate(material_labels)}
        material_map = self._build_material_code_grid(
            x_grid=x_grid,
            z_grid=z_grid,
            surface=surface,
            coating_stack=coating_stack,
            material_codes=material_codes,
        )

        return x_grid, z_grid, material_map, material_labels

    def plot_profile(self, output_filename: str | Path) -> None:
        """Save a plot showing the profile and the material stack."""

        positions, heights = self._tiled_profile_points(num_periods=3)
        x_grid, z_grid, material_map, material_labels = self._material_plot_data(num_periods=3)

        figure, axes = plt.subplots(
            2,
            1,
            figsize=(10, 7),
            sharex=True,
            gridspec_kw={"height_ratios": [2.0, 1.0]},
        )
        profile_axis, material_axis = axes

        profile_axis.fill_between(
            positions,
            0.0,
            heights,
            color="tab:blue",
            alpha=0.35,
            label="Grating groove",
        )
        profile_axis.plot(positions, heights, color="tab:blue", linewidth=2.0)
        profile_axis.set_xlim(0.0, 3.0 * self.period_nm)
        profile_axis.set_ylim(0.0, max(float(np.max(heights)) * 1.15, 1.0))
        profile_axis.set_ylabel("Depth (nm)")
        profile_axis.set_title(f"{self.__class__.__name__} Profile (Three Periods)")
        profile_axis.grid(True, alpha=0.3)
        profile_axis.legend(loc="best")

        base_colors = [
            "#f3f4f6",
            "#d97706",
            "#16a34a",
            "#2563eb",
            "#9333ea",
            "#dc2626",
        ]
        color_map = plt.matplotlib.colors.ListedColormap(base_colors[: len(material_labels)])
        material_axis.imshow(
            material_map,
            origin="lower",
            aspect="auto",
            extent=[x_grid[0], x_grid[-1], z_grid[0], z_grid[-1]],
            cmap=color_map,
            interpolation="nearest",
            vmin=0,
            vmax=max(len(material_labels) - 1, 1),
        )
        legend_handles = [
            plt.matplotlib.patches.Patch(color=color_map.colors[index], label=label)
            for index, label in enumerate(material_labels)
        ]
        material_axis.set_xlabel("x (nm)")
        material_axis.set_ylabel("z (nm)")
        material_axis.set_title("Material Stack")
        material_axis.legend(handles=legend_handles, loc="upper right")
        material_axis.grid(False)

        figure.tight_layout()
        figure.savefig(output_filename, dpi=150, bbox_inches="tight")
        plt.close(figure)

    def save_structure_debug_data(
        self,
        photon_energy_ev: float,
        output_directory: str | Path,
        *,
        num_periods: int = 1,
        n_inc: complex = 1.0 + 0.0j,
    ) -> None:
        """Write discretized structure debug data to CSV files.

        Args:
            photon_energy_ev: Photon energy used to interpolate optical constants.
            output_directory: Destination directory for the CSV files.
            num_periods: Number of grating periods to export in x.
            n_inc: Refractive index of the incident medium.
        """

        output_path = Path(output_directory)
        output_path.mkdir(parents=True, exist_ok=True)

        coating_stack = self.resolved_stack()
        x_grid = self._build_x_grid(num_periods=num_periods)
        z_grid = self._build_solver_z_grid(coating_stack)
        surface = self._surface_profile_on_grid(x_grid, num_periods=num_periods)
        material_labels = coating_stack.plot_material_names()
        material_codes = {label: index for index, label in enumerate(material_labels)}
        material_map = self._build_material_code_grid(
            x_grid=x_grid,
            z_grid=z_grid,
            surface=surface,
            coating_stack=coating_stack,
            material_codes=material_codes,
        )
        index_grid = self._build_refractive_index_grid(
            x_grid=x_grid,
            z_grid=z_grid,
            surface=surface,
            coating_stack=coating_stack,
            photon_energy_ev=photon_energy_ev,
            n_inc=n_inc,
        )

        np.savetxt(output_path / "x.csv", x_grid, delimiter=",")
        np.savetxt(output_path / "z.csv", z_grid, delimiter=",")
        np.savetxt(output_path / "surface.csv", surface, delimiter=",")
        np.savetxt(output_path / "material_id.csv", material_map, fmt="%d", delimiter=",")
        np.savetxt(output_path / "index_real.csv", np.real(index_grid), delimiter=",")
        np.savetxt(output_path / "index_imag.csv", np.imag(index_grid), delimiter=",")
        (output_path / "material_labels.txt").write_text("\n".join(material_labels), encoding="utf-8")

    def build_textures(
        self,
        photon_energy_ev: float,
        *,
        n_inc: complex = 1.0 + 0.0j,
    ) -> tuple[list[object], tuple[np.ndarray, np.ndarray]]:
        """Build RETICOLO-compatible textures and profile arrays."""

        coating_stack = self.resolved_stack()
        n_sub = interpolate_cxro_index(
            coating_stack.substrate_material,
            photon_energy_ev,
            base_dir=self.base_dir,
        )
        x = self._build_x_grid(num_periods=1)
        z = self._build_solver_z_grid(coating_stack)
        surface = self._surface_profile_on_grid(x, num_periods=1)
        index_grid = self._build_refractive_index_grid(
            x_grid=x,
            z_grid=z,
            surface=surface,
            coating_stack=coating_stack,
            photon_energy_ev=photon_energy_ev,
            n_inc=n_inc,
        )

        textures: list[object] = [n_inc]
        for row in range(index_grid.shape[0]):
            row_values = index_grid[row]
            row_changes = np.nonzero(np.diff(row_values) != 0)[0]
            if len(row_changes) == 0:
                textures.append(complex(row_values[0]))
            else:
                textures.append([x[row_changes + 1], row_values[row_changes]])

        textures.append(n_sub)
        thicknesses = np.concatenate(
            (
                np.asarray([0.0], dtype=float),
                np.full(index_grid.shape[0], self.z_resolution_nm, dtype=float),
                np.asarray([0.0], dtype=float),
            )
        )
        texture_sequence = np.arange(len(textures), dtype=int)
        return textures, (np.asarray(thicknesses, dtype=float), np.asarray(texture_sequence, dtype=int))

    def _build_x_grid(self, *, num_periods: int) -> np.ndarray:
        """Return the x grid used for discretized structure generation."""

        x_span_nm = num_periods * self.period_nm
        return np.linspace(0.0, x_span_nm, int(round(x_span_nm / self.x_resolution_nm)) + 1)

    def _build_solver_z_grid(self, coating_stack: BaseStack) -> np.ndarray:
        """Return the solver z grid, matching the Octave helper orientation."""

        if isinstance(coating_stack, MultilayerStack):
            thickness_nm = self.depth_nm + coating_stack.d_period_nm * coating_stack.n_bilayers
            thickness_nm += coating_stack.top_cap_thickness_nm + 25.0 * self.z_resolution_nm
        else:
            thickness_nm = self.depth_nm + coating_stack.total_thickness_nm + 5.0
        return np.linspace(thickness_nm, 0.0, int(round(thickness_nm / self.z_resolution_nm)) + 1)

    def _build_plot_z_grid(self, coating_stack: BaseStack) -> np.ndarray:
        """Return the plot z grid with increasing depth."""

        solver_z = self._build_solver_z_grid(coating_stack)
        return solver_z[::-1]

    def _surface_profile_on_grid(self, x_grid: np.ndarray, *, num_periods: int) -> np.ndarray:
        """Return the surface profile interpolated on the requested x grid."""

        positions, heights = self._tiled_profile_points(num_periods=num_periods)
        return np.interp(x_grid, positions, heights)

    def _build_material_code_grid(
        self,
        *,
        x_grid: np.ndarray,
        z_grid: np.ndarray,
        surface: np.ndarray,
        coating_stack: BaseStack,
        material_codes: dict[str, int],
    ) -> np.ndarray:
        """Return the discretized material-code grid for plotting."""

        substrate_code = material_codes[coating_stack.substrate_material]
        incident_code = material_codes["Incident Medium"]
        z_mesh = np.repeat(z_grid[:, None], x_grid.size, axis=1)
        material_map = np.full((z_grid.size, x_grid.size), substrate_code, dtype=int)

        if isinstance(coating_stack, MultilayerStack):
            bottom_material, top_material = coating_stack.bilayer_materials_bottom_up
            bottom_thickness_nm, top_thickness_nm = coating_stack.bilayer_thicknesses_bottom_up
            for bilayer_index in range(coating_stack.n_bilayers):
                lower_interface = surface + 1.0 + coating_stack.d_period_nm * bilayer_index
                middle_interface = lower_interface + bottom_thickness_nm
                upper_interface = lower_interface + coating_stack.d_period_nm
                lower_mask = (z_mesh >= lower_interface[None, :]) & (z_mesh < middle_interface[None, :])
                upper_mask = (z_mesh >= middle_interface[None, :]) & (z_mesh < upper_interface[None, :])
                material_map[lower_mask] = material_codes[bottom_material]
                material_map[upper_mask] = material_codes[top_material]

            current_top = surface + 1.0 + coating_stack.d_period_nm * coating_stack.n_bilayers
            if coating_stack.top_cap_material and coating_stack.top_cap_thickness_nm > 0.0:
                top_cap_upper = current_top + coating_stack.top_cap_thickness_nm
                top_cap_mask = (z_mesh >= current_top[None, :]) & (z_mesh < top_cap_upper[None, :])
                material_map[top_cap_mask] = material_codes[coating_stack.top_cap_material]
                current_top = top_cap_upper

            material_map[z_mesh >= current_top[None, :]] = incident_code
            return material_map

        current_lower = surface
        for material_name, layer_thickness_nm in coating_stack.layer_sequence_bottom_up():
            current_upper = current_lower + layer_thickness_nm
            layer_mask = (z_mesh >= current_lower[None, :]) & (z_mesh < current_upper[None, :])
            material_map[layer_mask] = material_codes[material_name]
            current_lower = current_upper
        material_map[z_mesh >= current_lower[None, :]] = incident_code
        return material_map

    def _build_refractive_index_grid(
        self,
        *,
        x_grid: np.ndarray,
        z_grid: np.ndarray,
        surface: np.ndarray,
        coating_stack: BaseStack,
        photon_energy_ev: float,
        n_inc: complex,
    ) -> np.ndarray:
        """Return the refractive-index grid used to derive RETICOLO textures."""

        n_sub = interpolate_cxro_index(
            coating_stack.substrate_material,
            photon_energy_ev,
            base_dir=self.base_dir,
        )
        z_mesh = np.repeat(z_grid[:, None], x_grid.size, axis=1)
        index_grid = np.full((z_grid.size, x_grid.size), n_sub, dtype=complex)

        if isinstance(coating_stack, MultilayerStack):
            n_material_a = interpolate_cxro_index(
                coating_stack.material_a,
                photon_energy_ev,
                base_dir=self.base_dir,
            )
            n_material_b = interpolate_cxro_index(
                coating_stack.material_b,
                photon_energy_ev,
                base_dir=self.base_dir,
            )
            refractive_index_by_material = {
                coating_stack.material_a: n_material_a,
                coating_stack.material_b: n_material_b,
            }
            bottom_material, top_material = coating_stack.bilayer_materials_bottom_up
            bottom_thickness_nm, _ = coating_stack.bilayer_thicknesses_bottom_up
            for bilayer_index in range(coating_stack.n_bilayers):
                lower_interface = surface + 1.0 + coating_stack.d_period_nm * bilayer_index
                middle_interface = lower_interface + bottom_thickness_nm
                upper_interface = lower_interface + coating_stack.d_period_nm
                lower_mask = (z_mesh >= lower_interface[None, :]) & (z_mesh < middle_interface[None, :])
                upper_mask = (z_mesh >= middle_interface[None, :]) & (z_mesh < upper_interface[None, :])
                index_grid[lower_mask] = refractive_index_by_material[bottom_material]
                index_grid[upper_mask] = refractive_index_by_material[top_material]

            current_top = surface + 1.0 + coating_stack.d_period_nm * coating_stack.n_bilayers
            if coating_stack.top_cap_material and coating_stack.top_cap_thickness_nm > 0.0:
                n_top_cap = interpolate_cxro_index(
                    coating_stack.top_cap_material,
                    photon_energy_ev,
                    base_dir=self.base_dir,
                )
                top_cap_upper = current_top + coating_stack.top_cap_thickness_nm
                top_cap_mask = (z_mesh >= current_top[None, :]) & (z_mesh < top_cap_upper[None, :])
                index_grid[top_cap_mask] = n_top_cap
                current_top = top_cap_upper

            index_grid[z_mesh >= current_top[None, :]] = n_inc
            return index_grid

        current_lower = surface
        for material_name, layer_thickness_nm in coating_stack.layer_sequence_bottom_up():
            refractive_index = interpolate_cxro_index(
                material_name,
                photon_energy_ev,
                base_dir=self.base_dir,
            )
            current_upper = current_lower + layer_thickness_nm
            layer_mask = (z_mesh >= current_lower[None, :]) & (z_mesh < current_upper[None, :])
            index_grid[layer_mask] = refractive_index
            current_lower = current_upper
        index_grid[z_mesh >= current_lower[None, :]] = n_inc
        return index_grid


@dataclass
class LaminarGrating(BaseGrating):
    """Laminar or trapezoidal grating profile."""

    left_wall_angle_deg: float = 15.0
    right_wall_angle_deg: float = 15.0

    def profile_points(self) -> tuple[np.ndarray, np.ndarray]:
        """Return the laminar profile for one period."""

        effective_width_ratio = 1.0 - self.width_to_period_ratio
        if effective_width_ratio < 0.0 or effective_width_ratio > 1.0:
            raise ValueError("width_to_period_ratio must be in [0, 1].")

        if self.left_wall_angle_deg == 0.0 and self.right_wall_angle_deg == 0.0:
            edge_x = effective_width_ratio * self.period_nm
            positions = np.array(
                [0.0, self.period_nm / 2 - edge_x / 2, self.period_nm / 2 + edge_x / 2, self.period_nm],
                dtype=float,
            )
            heights = np.array([0.0, 0.0, self.depth_nm, 0.0], dtype=float)
            return positions, heights

        pos1 = [
            (self.period_nm - effective_width_ratio * self.period_nm) / 2
            - self.depth_nm / np.tan(np.deg2rad(self.left_wall_angle_deg)),
            0.0,
        ]
        pos2 = [(self.period_nm - effective_width_ratio * self.period_nm) / 2, self.depth_nm]
        pos3 = [(self.period_nm + effective_width_ratio * self.period_nm) / 2, self.depth_nm]
        pos4 = [
            (self.period_nm + effective_width_ratio * self.period_nm) / 2
            + self.depth_nm / np.tan(np.deg2rad(self.right_wall_angle_deg)),
            0.0,
        ]
        profile_points = np.array(
            [[0.0, 0.0], pos1, pos2, pos3, pos4, [self.period_nm, 0.0]],
            dtype=float,
        )
        return profile_points[:, 0], profile_points[:, 1]


@dataclass
class BlazedGrating(BaseGrating):
    """Blazed sawtooth grating profile."""

    blaze_angle_deg: float = 0.9

    def __post_init__(self) -> None:
        """Set depth from blaze angle and period."""

        self.depth_nm = self.period_nm * np.tan(np.deg2rad(self.blaze_angle_deg))

    def profile_points(self) -> tuple[np.ndarray, np.ndarray]:
        """Return the blazed sawtooth profile for one period."""

        if self.blaze_angle_deg <= 0.0:
            raise ValueError("blaze_angle_deg must be > 0.")

        # Use a one-step drop at period end to emulate the sawtooth reset.
        edge_step_nm = max(self.x_resolution_nm, self.period_nm * 1e-6)
        x_peak = self.period_nm - edge_step_nm

        positions = np.array([0.0, x_peak, self.period_nm], dtype=float)
        heights = np.array([0.0, self.depth_nm, 0.0], dtype=float)
        return positions, heights
