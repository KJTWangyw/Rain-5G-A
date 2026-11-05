"""Specific rain attenuation coefficients of Recommendation ITU-R P.838-3."""
from __future__ import annotations

import numpy as np

# Tables 1-4 of ITU-R P.838-3: (a_j, b_j, c_j) and the linear terms (m, c).
_KH = ((-5.33980, -0.10008, 1.13098), (-0.35351, 1.26970, 0.45400),
       (-0.23789, 0.86036, 0.15354), (-0.94158, 0.64552, 0.16817))
_KH_M, _KH_C = -0.18961, 0.71147
_KV = ((-3.80595, 0.56934, 0.81061), (-3.44965, -0.22911, 0.51059),
       (-0.39902, 0.73042, 0.11899), (0.50167, 1.07319, 0.27195))
_KV_M, _KV_C = -0.16398, 0.63297
_AH = ((-0.14318, 1.82442, -0.55187), (0.29591, 0.77564, 0.19822), (0.32177, 0.63773, 0.13164),
       (-5.37610, -0.96230, 1.47828), (16.1721, -3.29980, 3.43990))
_AH_M, _AH_C = 0.67849, -1.95537
_AV = ((-0.07771, 2.33840, -0.76284), (0.56727, 0.95545, 0.54039), (-0.20238, 1.14520, 0.26809),
       (-48.2991, 0.791669, 0.116226), (48.5833, 0.791459, 0.116479))
_AV_M, _AV_C = -0.053739, 0.83433


def _regression(f_ghz, k_table, k_lin, a_table, a_lin):
    x = np.log10(np.asarray(f_ghz, dtype=float))
    log_k = sum(a * np.exp(-((x - b) / c) ** 2) for a, b, c in k_table) + k_lin[0] * x + k_lin[1]
    alpha = sum(a * np.exp(-((x - b) / c) ** 2) for a, b, c in a_table) + a_lin[0] * x + a_lin[1]
    return 10.0 ** log_k, alpha


def p838_horizontal(f_ghz):
    """(k_H, alpha_H) at frequency f in GHz, from equations (2) and (3) of ITU-R P.838-3."""
    k, a = _regression(f_ghz, _KH, (_KH_M, _KH_C), _AH, (_AH_M, _AH_C))
    return (float(k), float(a)) if np.ndim(k) == 0 else (k, a)


def p838_vertical(f_ghz):
    """(k_V, alpha_V) at frequency f in GHz."""
    k, a = _regression(f_ghz, _KV, (_KV_M, _KV_C), _AV, (_AV_M, _AV_C))
    return (float(k), float(a)) if np.ndim(k) == 0 else (k, a)


def p838(f_ghz, polarization):
    """(k, alpha) for 'H'/'horizontal' or 'V'/'vertical' polarization (scalars or arrays of equal length)."""
    pol = np.atleast_1d(np.asarray(polarization, dtype=str))
    kh, ah = p838_horizontal(np.atleast_1d(f_ghz))
    kv, av = p838_vertical(np.atleast_1d(f_ghz))
    vertical = np.char.lower(np.char.strip(pol)).astype(str)
    vertical = np.array([p.startswith("v") for p in vertical])
    k, a = np.where(vertical, kv, kh), np.where(vertical, av, ah)
    return (float(k[0]), float(a[0])) if np.ndim(f_ghz) == 0 else (k, a)


def rain_rate(attenuation_db, k, alpha, length_km):
    """Path rain rate (mm/h) producing ``attenuation_db`` over ``length_km`` with A = k L R^alpha."""
    return (np.asarray(attenuation_db) / (k * length_km)) ** (1.0 / alpha)
