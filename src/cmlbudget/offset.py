"""Length-independent attenuation offset identified from links of different lengths.

Over many links observing the same rain climate, the wet-period mean attenuation of link i is

    a_i = W + beta * e_i,      e_i = k_i L_i R^alpha_i averaged with the reference rain rate,

so the intercept W of a regression across links is the part of the attenuation that does not scale with the path
(wet antennas and radomes, plus any other wet-period offset). Repeating the regression within bins of the reference
rain rate gives a curve W(R); rain rate is then retrieved by inverting g(R) = W(R) + k L R^alpha.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def fit_offset(att_mean, exp_mean):
    """Least-squares intercept and slope of att_mean = W + beta * exp_mean."""
    X = np.c_[np.ones(len(att_mean)), np.asarray(exp_mean, dtype=float)]
    (w, beta), *_ = np.linalg.lstsq(X, np.asarray(att_mean, dtype=float), rcond=None)
    return float(w), float(beta)


def cluster_bootstrap(df, statistic, draws, seed, cluster="Link"):
    """Statistic recomputed on resamples of whole clusters (both directions of a physical link stay together)."""
    rng = np.random.default_rng(seed)
    groups = [g.index.values for _, g in df.groupby(cluster)]
    out = []
    for _ in range(draws):
        idx = np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
        out.append(statistic(df.loc[idx]))
    return np.array(out)


def bin_means(rain_att, expected, radar, bins, min_count):
    """Per rain-rate bin: count and mean measured attenuation, radar-expected attenuation and radar rain rate."""
    ok = np.isfinite(rain_att) & np.isfinite(expected) & np.isfinite(radar)
    b = np.digitize(radar[ok], bins) - 1
    rows = []
    for k in range(len(bins) - 1):
        s = b == k
        if s.sum() >= min_count:
            rows.append(dict(bin=k, n=int(s.sum()), a=float(np.mean(rain_att[ok][s])),
                             e=float(np.mean(expected[ok][s])), R=float(np.mean(radar[ok][s]))))
    return rows


def offset_curve(binned):
    """W, beta and mean rain rate per bin from a table with columns bin, a, e, R (one row per link and bin)."""
    rows = []
    for k, g in binned.groupby("bin"):
        w, beta = fit_offset(g.a.values, g.e.values)
        rows.append(dict(bin=k, links=len(g), R=float(g.R.mean()), W=w, beta=beta,
                         resid_sd=float(np.std(g.a.values - w - beta * g.e.values))))
    return pd.DataFrame(rows)


def offset_function(R_bins, W_bins):
    """W(R): non-decreasing, piecewise linear in log10 R through the bin estimates, constant outside them."""
    R_bins = np.asarray(R_bins, dtype=float)
    W_mono = np.maximum.accumulate(np.asarray(W_bins, dtype=float))
    x = np.log10(R_bins)

    def W(R):
        return np.interp(np.log10(np.clip(np.asarray(R, dtype=float), R_bins[0], R_bins[-1])), x, W_mono)
    return W


def invert(rain_att, kL, alpha, W, r_max=300.0, n_grid=800):
    """Rain rate solving W(R) + kL R^alpha = rain_att (zero below W(0)); NaN stays NaN."""
    grid = np.r_[0.0, np.logspace(-3, np.log10(r_max), n_grid)]
    g = W(grid) + kL * grid ** alpha
    g[0] = W(0.0)
    g = np.maximum.accumulate(g)
    a = np.asarray(rain_att, dtype=float)
    out = np.interp(np.nan_to_num(a, nan=-np.inf), g, grid, left=0.0, right=r_max)
    out[~np.isfinite(a)] = np.nan
    return out
