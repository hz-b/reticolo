from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.linalg import expm


ArrayLike = Any


@dataclass
class Res1Parameters:
    trace: int = 0


@dataclass
class Res2Parameters:
    result: int = 1


@dataclass
class Parameters:
    dim: int = 1
    polarization: int = 1
    not_io: int = 0
    res1: Res1Parameters = field(default_factory=Res1Parameters)
    res2: Res2Parameters = field(default_factory=Res2Parameters)


@dataclass
class Texture1D:
    period: float
    breaks: np.ndarray
    refractive_index: np.ndarray
    epsilon_fourier: np.ndarray
    homogeneous_index: complex | None = None

    @property
    def signature(self) -> tuple[Any, ...]:
        if self.homogeneous_index is not None:
            n = self.homogeneous_index
            return ("homogeneous", round(float(np.real(n)), 12), round(float(np.imag(n)), 12))
        return (
            "patterned",
            tuple(np.round(self.breaks, 12)),
            tuple((round(float(np.real(v)), 12), round(float(np.imag(v)), 12)) for v in self.refractive_index),
        )


@dataclass
class Res1Result:
    wavelength: float
    period: float
    orders: np.ndarray
    beta0: float
    polarization: int
    textures: list[Texture1D]


@dataclass
class DiffractionResult:
    order: np.ndarray
    theta: np.ndarray
    efficiency: np.ndarray
    amplitude: np.ndarray


@dataclass
class Res2Result:
    inc_top_reflected: DiffractionResult
    inc_top_transmitted: DiffractionResult


def res0(dim: int) -> Parameters:
    if dim == 0 or abs(dim) != 1:
        raise NotImplementedError("The Python port currently supports 1D TE/TM-style entry only.")
    return Parameters(dim=1, polarization=1 if dim > 0 else -1)


def res1(
    wavelength: float,
    period: float,
    textures: list[ArrayLike],
    nn: int | tuple[int, int] | list[int],
    beta0: float,
    parm: Parameters | None = None,
) -> Res1Result:
    parm = parm or res0(1)
    if parm.polarization != 1:
        raise NotImplementedError("The native Python port currently implements the 1D TE path only.")

    orders = _normalize_orders(nn)
    converted = [_convert_texture(texture, period, orders) for texture in textures]
    return Res1Result(
        wavelength=float(wavelength),
        period=float(period),
        orders=orders,
        beta0=float(beta0),
        polarization=parm.polarization,
        textures=converted,
    )


def res2(aa: Res1Result, profile: tuple[ArrayLike, ArrayLike], parm: Parameters | None = None) -> Res2Result:
    parm = parm or res0(1)
    if parm.polarization != 1:
        raise NotImplementedError("The native Python port currently implements the 1D TE path only.")

    thicknesses = np.asarray(profile[0], dtype=float)
    texture_indices = np.asarray(profile[1], dtype=int)
    if thicknesses.shape != texture_indices.shape:
        raise ValueError("Profile thicknesses and texture indices must have the same length.")
    if len(texture_indices) < 2:
        raise ValueError("The profile must contain at least superstrate and substrate textures.")

    top_texture = aa.textures[int(texture_indices[0])]
    bottom_texture = aa.textures[int(texture_indices[-1])]
    if top_texture.homogeneous_index is None or bottom_texture.homogeneous_index is None:
        raise ValueError("The Python port expects homogeneous superstrate and substrate textures.")

    compressed: list[tuple[float, Texture1D]] = []
    for thickness, texture_index in zip(thicknesses[1:-1], texture_indices[1:-1]):
        if thickness <= 0:
            continue
        texture = aa.textures[int(texture_index)]
        if compressed and compressed[-1][1].signature == texture.signature:
            compressed[-1] = (compressed[-1][0] + float(thickness), texture)
        else:
            compressed.append((float(thickness), texture))

    reflected, transmitted = _solve_te_stack(
        wavelength=aa.wavelength,
        period=aa.period,
        orders=aa.orders,
        beta0=aa.beta0,
        n_top=top_texture.homogeneous_index,
        n_bottom=bottom_texture.homogeneous_index,
        layers=compressed,
    )

    return Res2Result(
        inc_top_reflected=reflected,
        inc_top_transmitted=transmitted,
    )


def _normalize_orders(nn: int | tuple[int, int] | list[int]) -> np.ndarray:
    if np.isscalar(nn):
        n = int(nn)
        return np.arange(-n, n + 1, dtype=int)
    values = np.asarray(nn, dtype=int).ravel()
    if values.size != 2:
        raise ValueError("Only scalar or two-value 1D order definitions are supported.")
    return np.arange(values[0], values[1] + 1, dtype=int)


def _convert_texture(texture: ArrayLike, period: float, orders: np.ndarray) -> Texture1D:
    if not isinstance(texture, (list, tuple)):
        refractive_index = complex(texture)
        return Texture1D(
            period=period,
            breaks=np.array([0.0, period], dtype=float),
            refractive_index=np.array([refractive_index], dtype=complex),
            epsilon_fourier=_piecewise_fourier_coefficients(
                np.array([0.0, period], dtype=float),
                np.array([refractive_index], dtype=complex),
                period=period,
                max_order=2 * int(np.max(np.abs(orders))),
            ),
            homogeneous_index=refractive_index,
        )

    if len(texture) == 1:
        refractive_index = complex(texture[0])
        return _convert_texture(refractive_index, period, orders)

    if len(texture) != 2:
        raise ValueError("1D Python textures must be homogeneous or [x_positions, n_left].")

    x_positions = np.asarray(texture[0], dtype=float).ravel()
    n_left = np.asarray(texture[1], dtype=complex).ravel()
    if x_positions.size != n_left.size:
        raise ValueError("Texture boundary and refractive-index vectors must have the same length.")
    if x_positions.size == 0:
        raise ValueError("Patterned textures require at least one discontinuity.")

    breaks = np.concatenate(([0.0], x_positions, [period]))
    refractive_index = np.concatenate((n_left, [n_left[0]]))
    return Texture1D(
        period=period,
        breaks=breaks,
        refractive_index=refractive_index,
        epsilon_fourier=_piecewise_fourier_coefficients(
            breaks,
            refractive_index,
            period=period,
            max_order=2 * int(np.max(np.abs(orders))),
        ),
    )


def _piecewise_fourier_coefficients(
    breaks: np.ndarray,
    refractive_index: np.ndarray,
    period: float,
    max_order: int,
) -> np.ndarray:
    g = np.arange(-max_order, max_order + 1, dtype=int)
    epsilon = refractive_index**2
    coeffs = np.zeros_like(g, dtype=complex)
    for idx, order in enumerate(g):
        if order == 0:
            coeffs[idx] = np.sum(epsilon * (breaks[1:] - breaks[:-1])) / period
            continue
        phase_right = np.exp(-1j * 2 * np.pi * order * breaks[1:] / period)
        phase_left = np.exp(-1j * 2 * np.pi * order * breaks[:-1] / period)
        coeffs[idx] = np.sum(
            epsilon * (phase_right - phase_left) / (-1j * 2 * np.pi * order)
        )
    return coeffs


def _convolution_matrix(coefficients: np.ndarray, orders: np.ndarray) -> np.ndarray:
    max_order = (len(coefficients) - 1) // 2
    conv = np.zeros((len(orders), len(orders)), dtype=complex)
    for row, m in enumerate(orders):
        for col, n in enumerate(orders):
            conv[row, col] = coefficients[(m - n) + max_order]
    return conv


def _solve_te_stack(
    wavelength: float,
    period: float,
    orders: np.ndarray,
    beta0: float,
    n_top: complex,
    n_bottom: complex,
    layers: list[tuple[float, Texture1D]],
) -> tuple[DiffractionResult, DiffractionResult]:
    k0 = 2 * np.pi / wavelength
    kx = k0 * beta0 + (2 * np.pi * orders / period)
    kx_matrix_sq = np.diag(kx**2)
    basis_size = len(orders)

    transfer = np.eye(2 * basis_size, dtype=complex)
    for thickness, texture in layers:
        epsilon_conv = _convolution_matrix(texture.epsilon_fourier, orders)
        operator = kx_matrix_sq - (k0**2) * epsilon_conv
        state_matrix = np.block(
            [
                [np.zeros((basis_size, basis_size), dtype=complex), np.eye(basis_size, dtype=complex)],
                [operator, np.zeros((basis_size, basis_size), dtype=complex)],
            ]
        )
        transfer = expm(state_matrix * thickness) @ transfer

    kz_top = np.array([_kz_branch((k0 * n_top) ** 2 - value**2) for value in kx], dtype=complex)
    kz_bottom = np.array([_kz_branch((k0 * n_bottom) ** 2 - value**2) for value in kx], dtype=complex)
    derivative_top = 1j * np.diag(kz_top)
    derivative_bottom = 1j * np.diag(kz_bottom)

    incident = np.zeros(basis_size, dtype=complex)
    incident[np.where(orders == 0)[0][0]] = 1.0

    a11 = transfer[:basis_size, :basis_size]
    a12 = transfer[:basis_size, basis_size:]
    a21 = transfer[basis_size:, :basis_size]
    a22 = transfer[basis_size:, basis_size:]

    lhs = (a21 - derivative_bottom @ a11) - (a22 - derivative_bottom @ a12) @ derivative_top
    rhs = -((a21 - derivative_bottom @ a11) + (a22 - derivative_bottom @ a12) @ derivative_top) @ incident
    reflected_amplitude = np.linalg.solve(lhs, rhs)
    transmitted_amplitude = a11 @ (incident + reflected_amplitude) + a12 @ (
        derivative_top @ (incident - reflected_amplitude)
    )

    incident_kz = kz_top[np.where(orders == 0)[0][0]]
    reflected_efficiency = np.real(kz_top / incident_kz) * np.abs(reflected_amplitude) ** 2
    transmitted_efficiency = np.real(kz_bottom / incident_kz) * np.abs(transmitted_amplitude) ** 2

    reflected_theta = _angles_from_kx(kx, k0, n_top)
    transmitted_theta = _angles_from_kx(kx, k0, n_bottom)

    return (
        DiffractionResult(
            order=orders.copy(),
            theta=reflected_theta,
            efficiency=reflected_efficiency,
            amplitude=reflected_amplitude,
        ),
        DiffractionResult(
            order=orders.copy(),
            theta=transmitted_theta,
            efficiency=transmitted_efficiency,
            amplitude=transmitted_amplitude,
        ),
    )


def _angles_from_kx(kx: np.ndarray, k0: float, refractive_index: complex) -> np.ndarray:
    if abs(np.imag(refractive_index)) > 1e-12:
        n_for_angles = np.real(refractive_index)
    else:
        n_for_angles = float(np.real(refractive_index))
    ratio = np.real(kx / (k0 * n_for_angles))
    ratio = np.clip(ratio, -1.0, 1.0)
    return np.degrees(np.arcsin(ratio))


def _kz_branch(value: complex) -> complex:
    kz = np.sqrt(value + 0j)
    if kz.real < 0:
        kz = -kz
    if abs(kz.real) < 1e-12 and kz.imag < 0:
        kz = -kz
    return kz
