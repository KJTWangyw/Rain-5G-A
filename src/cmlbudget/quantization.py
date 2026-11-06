"""Rounding quantizer of a received level and the mean error left after averaging.

A report y = Q(x) = u + D * round((x - u) / D) has the rounding error eps(x) = Q(x) - x, a sawtooth of period D:

    eps(x) = sum_{k>=1} (D / (pi k)) (-1)^k sin(2 pi k (x - u) / D).

If the level fluctuates as x ~ N(c, s^2) while reports are averaged, the mean error is the sawtooth smoothed by
the Gaussian (each harmonic damped by exp(-2 pi^2 k^2 s^2 / D^2)):

    b(c; s) = E[eps(x)] = sum_k (D / (pi k)) (-1)^k exp(-2 pi^2 k^2 s^2 / D^2) sin(2 pi k (c - u) / D).

b is the error that averaging cannot remove. Its root mean square over a grid offset uniform on [0, D) is

    rms_b(s) = sqrt( sum_k (D / (pi k))^2 exp(-4 pi^2 k^2 s^2 / D^2) / 2 ),

which equals D / sqrt(12) for s = 0 and falls below 0.01 D for s > 0.42 D.
"""
from __future__ import annotations

import numpy as np

TERMS = 60


def quantize(x, step, offset=0.0):
    """Round ``x`` to the grid offset + step * integer."""
    return offset + step * np.round((np.asarray(x, dtype=float) - offset) / step)


def dither_bias(c, sigma, step, offset=0.0, terms=TERMS):
    """Mean rounding error E[Q(x) - x] for x ~ N(c, sigma^2) (arrays broadcast)."""
    c, sigma, offset = np.broadcast_arrays(np.asarray(c, float), np.asarray(sigma, float), np.asarray(offset, float))
    y = (c - offset) / step
    out = np.zeros(c.shape)
    for k in range(1, terms + 1):
        out += (step / (np.pi * k)) * (-1.0) ** k * np.exp(-2.0 * (np.pi * k * sigma / step) ** 2) * np.sin(2 * np.pi * k * y)
    return out


def dither_bias_rms(sigma, step, terms=TERMS):
    """Root mean square of ``dither_bias`` over a uniform grid offset."""
    sigma = np.asarray(sigma, float)
    k = np.arange(1, terms + 1).reshape((-1,) + (1,) * sigma.ndim)
    return np.sqrt(np.sum((step / (np.pi * k)) ** 2 * np.exp(-4.0 * (np.pi * k * sigma / step) ** 2), axis=0) / 2.0)
