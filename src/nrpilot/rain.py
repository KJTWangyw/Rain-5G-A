"""Specific rain attenuation coefficients of Recommendation ITU-R P.838-3 (horizontal polarisation)."""
from __future__ import annotations

import numpy as np

# Tables 1 and 3 of ITU-R P.838-3: (a_j, b_j, c_j) and the linear terms.
_KH = ((-5.33980, -0.10008, 1.13098), (-0.35351, 1.26970, 0.45400),
       (-0.23789, 0.86036, 0.15354), (-0.94158, 0.64552, 0.16817))
_KH_M, _KH_C = -0.18961, 0.71147
_AH = ((-0.14318, 1.82442, -0.55187), (0.29591, 0.77564, 0.19822), (0.32177, 0.63773, 0.13164),
       (-5.37610, -0.96230, 1.47828), (16.1721, -3.29980, 3.43990))
_AH_M, _AH_C = 0.67849, -1.95537


def p838_horizontal(f_ghz: float) -> tuple[float, float]:
    """(k_H, alpha_H) at frequency f in GHz, from equations (2) and (3) of ITU-R P.838-3."""
    x = np.log10(f_ghz)
    log_k = sum(a * np.exp(-((x - b) / c) ** 2) for a, b, c in _KH) + _KH_M * x + _KH_C
    alpha = sum(a * np.exp(-((x - b) / c) ** 2) for a, b, c in _AH) + _AH_M * x + _AH_C
    return float(10.0 ** log_k), float(alpha)


def rain_rate(attenuation_db: float, k: float, alpha: float, length_km: float) -> float:
    """Path rain rate (mm/h) producing ``attenuation_db`` over ``length_km`` with A = k L R^alpha."""
    return (attenuation_db / (k * length_km)) ** (1.0 / alpha)
