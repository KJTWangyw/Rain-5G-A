"""Profile amplitude estimator with a finite clear-sky reference.

Model (per pilot j): r_j = h_j + e_0j, z_tj = a_t exp(i phi_t) h_j + e_tj, with e_0j ~ CN(0, v_0j),
e_tj ~ CN(0, v_tj) and the channel h_j unknown. Profiling h_j out of the joint likelihood gives

    a_hat = argmin_a  sum_j (|z_j|^2 + a^2 |r_j|^2) / V_j(a) - 2 a |sum_j r_j^* z_j / V_j(a)|,
    V_j(a) = v_tj + a^2 v_0j,

and A_hat = -20 log10(a_hat) in dB. With white noise (v_0j = v_0, v_tj = v) the stationary point is the
unique positive root of  G v_0 a^2 + (R v - Z v_0) a - G v = 0  with R = sum|r|^2, Z = sum|z|^2,
G = |sum r^* z|.
"""
from __future__ import annotations

import numpy as np

BOUNDS_DB = (-20.0, 80.0)


def _amp_bounds(bounds_db):
    lo_db, hi_db = bounds_db
    return 10.0 ** (-hi_db / 20.0), 10.0 ** (-lo_db / 20.0)


def profile_white(r, z, v0, vt, bounds_db=BOUNDS_DB):
    """Closed-form profile estimate (dB) for white noise; the last axis indexes pilots."""
    R = np.sum(np.abs(r) ** 2, axis=-1)
    Z = np.sum(np.abs(z) ** 2, axis=-1)
    G = np.abs(np.sum(np.conj(r) * z, axis=-1))
    R, Z, G = np.broadcast_arrays(R, Z, G)
    qb = R * vt - Z * v0
    disc = np.sqrt(qb ** 2 + 4.0 * G ** 2 * v0 * vt)
    a_lo, a_hi = _amp_bounds(bounds_db)
    with np.errstate(divide="ignore", invalid="ignore"):
        root = np.where(qb >= 0, 2.0 * G * vt / (qb + disc), (disc - qb) / (2.0 * G * v0))
    root = np.clip(np.nan_to_num(root, nan=a_lo, posinf=a_hi), a_lo, a_hi)

    def objective(a):
        return (Z - 2.0 * a * G + a * a * R) / (vt + a * a * v0)

    cands = np.stack(np.broadcast_arrays(root, a_lo, a_hi))
    vals = np.stack([objective(c) for c in cands])
    a = np.take_along_axis(cands, np.argmin(vals, axis=0)[None], axis=0)[0]
    return -20.0 * np.log10(a)


def profile_numeric(r, z, v0, vt, bounds_db=BOUNDS_DB, grid_step_db=0.5, iters=40):
    """Numeric profile estimate (dB) for per-pilot noise variances: coarse grid, then golden section."""
    J = r.shape[-1]
    v0 = np.broadcast_to(np.asarray(v0, dtype=float), (J,))
    vt = np.broadcast_to(np.asarray(vt, dtype=float), (J,))
    absr2, absz2, rz = np.abs(r) ** 2, np.abs(z) ** 2, np.conj(r) * z

    def objective(A):
        a = (10.0 ** (-A / 20.0))[..., None]
        V = vt + a * a * v0
        return np.sum((absz2 + a * a * absr2) / V, axis=-1) - 2.0 * a[..., 0] * np.abs(np.sum(rz / V, axis=-1))

    lo, hi = bounds_db
    batch = np.broadcast_shapes(r.shape[:-1], z.shape[:-1])
    best_val, best_A = np.full(batch, np.inf), np.full(batch, lo, dtype=float)
    for A in np.arange(lo, hi + 0.5 * grid_step_db, grid_step_db):
        val = objective(np.full(batch, A))
        better = val < best_val
        best_val, best_A = np.where(better, val, best_val), np.where(better, A, best_A)
    left, right = np.maximum(best_A - grid_step_db, lo), np.minimum(best_A + grid_step_db, hi)
    g = (np.sqrt(5.0) - 1.0) / 2.0
    x1, x2 = right - g * (right - left), left + g * (right - left)
    f1, f2 = objective(x1), objective(x2)
    for _ in range(iters):
        move_left = f1 < f2
        right, left = np.where(move_left, x2, right), np.where(move_left, left, x1)
        x1, x2 = right - g * (right - left), left + g * (right - left)
        f1, f2 = objective(x1), objective(x2)
    return 0.5 * (left + right)
