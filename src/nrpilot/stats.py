"""Coverage statistics."""
from __future__ import annotations

import math
from statistics import NormalDist

import numpy as np


def z_value(level: float) -> float:
    return NormalDist().inv_cdf(0.5 + level / 2.0)


def wilson(k: int, n: int, conf: float = 0.95) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return float("nan"), float("nan")
    z = z_value(conf)
    p = k / n
    den = 1.0 + z ** 2 / n
    centre = (p + z ** 2 / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z ** 2 / (4 * n * n)) / den
    return centre - half, centre + half


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact McNemar p-value for b and c discordant pairs."""
    m = b + c
    if m == 0:
        return 1.0
    tail = sum(math.comb(m, i) for i in range(0, min(b, c) + 1))
    return min(1.0, 2 * tail / (1 << m))


def chi2_variance_interval(sample_var: float, dof: int, conf: float = 0.95) -> tuple[float, float]:
    """Approximate interval for a variance (Wilson-Hilferty quantiles of chi^2_dof / dof)."""
    z = z_value(conf)
    q = lambda s: (1.0 - 2.0 / (9 * dof) + s * z * math.sqrt(2.0 / (9 * dof))) ** 3
    return sample_var / q(+1), sample_var / q(-1)


def fmt_pct(x: float, digits: int = 1) -> str:
    return f"{100.0 * x:.{digits}f}"


def as_float(x):
    return float(np.asarray(x))
