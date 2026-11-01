"""Parametric bootstrap of the averaged attenuation with a shared reference draw.

Each draw regenerates one reference error, reused by all current observations, and the current noise:
    r*_j = h~_j + e*_0j,   z*_tj = a_hat_t h~_j + e*_tj,
with h~ the reference shrunk by its known noise energy. The estimator profiles the common phase, so the
draws use zero phase. The basic interval is [2 A_bar - Q*_(1-g/2), 2 A_bar - Q*_(g/2)].
"""
from __future__ import annotations

import numpy as np

from .estimator import BOUNDS_DB, profile_white


def crandn(rng, shape, var):
    """Circular complex Gaussian with total variance ``var``."""
    return np.sqrt(np.asarray(var) / 2.0) * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape))


def shrunk_channel(r, v0):
    """Reference scaled so that its energy matches the unbiased estimate of sum |h|^2."""
    R = float(np.sum(np.abs(r) ** 2))
    signal = max(R - r.shape[-1] * float(v0), 0.05 * R)
    return r * np.sqrt(signal / R)


def bootstrap_means(rng, h_plug, a_hat, v0, v, n_avg, draws, chunk=100, bounds_db=BOUNDS_DB):
    """Bootstrap replicates of the mean of the first n estimates, shape (draws, len(n_avg))."""
    T, J = len(a_hat), len(h_plug)
    n_avg = np.asarray(n_avg)
    out = np.empty((draws, len(n_avg)))
    for start in range(0, draws, chunk):
        b = min(chunk, draws - start)
        rb = h_plug[None, None, :] + crandn(rng, (b, 1, J), v0)
        zb = a_hat[None, :, None] * h_plug[None, None, :] + crandn(rng, (b, T, J), v)
        Ab = profile_white(rb, zb, v0, v, bounds_db)
        out[start:start + b] = np.cumsum(Ab, axis=1)[:, n_avg - 1] / n_avg
    return out


def basic_interval(estimate, replicates, level):
    lo, hi = np.quantile(replicates, [0.5 - level / 2.0, 0.5 + level / 2.0])
    return 2.0 * estimate - hi, 2.0 * estimate - lo
