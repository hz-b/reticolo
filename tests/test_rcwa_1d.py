from __future__ import annotations

import numpy as np

from reticolo_py.rcwa_1d import res0, res1, res2


def test_homogeneous_reflection_matches_fresnel_te() -> None:
    wavelength_nm = 500.0
    period_nm = 1000.0
    n_top = 1.0 + 0.0j
    n_bottom = 1.5 + 0.0j
    incidence_theta_deg = 15.0
    beta0 = np.sin(np.deg2rad(incidence_theta_deg))

    parm = res0(1)
    aa = res1(wavelength_nm, period_nm, [n_top, n_bottom], 2, beta0, parm)
    ef = res2(aa, (np.array([0.0, 0.0]), np.array([0, 1])), parm)

    order_zero = np.where(ef.inc_top_reflected.order == 0)[0][0]
    reflected_efficiency = ef.inc_top_reflected.efficiency[order_zero]

    cos_top = np.cos(np.deg2rad(incidence_theta_deg))
    sin_bottom = np.real(n_top / n_bottom) * np.sin(np.deg2rad(incidence_theta_deg))
    cos_bottom = np.sqrt(1.0 - sin_bottom**2)
    fresnel = abs((n_top * cos_top - n_bottom * cos_bottom) / (n_top * cos_top + n_bottom * cos_bottom)) ** 2

    assert np.isclose(reflected_efficiency, fresnel, atol=1e-10)
