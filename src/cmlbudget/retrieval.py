"""Wet/dry classification, event references and rain rate for one sub-link on a 1-min grid.

Every wet event uses one clear-sky reference: the mean attenuation over the dry minutes of the ``window``
minutes before the event starts. The same reference is subtracted from every minute of the event, so its error
is shared by all of them.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def radar_on_minutes(minute_time, radar_time, radar_R, max_age_s):
    """Most recent radar value for every minute (NaN when older than ``max_age_s``)."""
    idx = np.clip(np.searchsorted(radar_time, minute_time, side="right") - 1, 0, len(radar_time) - 1)
    out = np.asarray(radar_R, dtype=float)[idx]
    stale = (minute_time - radar_time[idx] > max_age_s) | (minute_time < radar_time[0])
    out[stale] = np.nan
    return out


def wet_dry(radar_min, threshold, margin):
    """wet: radar rain above ``threshold`` within +-``margin`` minutes; dry: radar available and not above it."""
    s = pd.Series(radar_min)
    win = 2 * margin + 1
    wet = s.rolling(win, center=True, min_periods=1).max().to_numpy() > threshold
    full = s.rolling(win, center=True, min_periods=1).count().to_numpy() == win
    return wet, (~wet) & full


def runs(mask):
    """Start and stop indices of the contiguous True runs of ``mask``."""
    d = np.diff(np.concatenate([[0], np.asarray(mask, dtype=np.int8), [0]]))
    return np.flatnonzero(d == 1), np.flatnonzero(d == -1)


def _cumsum0(x):
    return np.concatenate([[0.0], np.cumsum(x, dtype=float)])


def window_stats(values, ok, lo, hi, var_within=None):
    """Count, mean and (optionally) total 10-s variance of ``values`` over ``ok`` minutes in [lo, hi).

    The total variance combines the variance of the 1-min means with the mean within-minute variance.
    """
    okv = ok & np.isfinite(values)
    if var_within is not None:
        okv &= np.isfinite(var_within)
    x = np.where(okv, values, 0.0)
    n = _cumsum0(okv)[hi] - _cumsum0(okv)[lo]
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = (_cumsum0(x)[hi] - _cumsum0(x)[lo]) / n
        if var_within is None:
            return n, mean
        m2 = (_cumsum0(x * x)[hi] - _cumsum0(x * x)[lo]) / n
        vw = (_cumsum0(np.where(okv, var_within, 0.0))[hi] - _cumsum0(np.where(okv, var_within, 0.0))[lo]) / n
    return n, mean, np.sqrt(np.maximum(m2 - mean ** 2, 0.0) + vw)


def event_reference(att, dry, starts, window, min_dry):
    """Mean dry attenuation over [start - window, start) for each event; NaN with fewer than ``min_dry`` minutes."""
    lo = np.maximum(np.asarray(starts) - window, 0)
    n, ref = window_stats(att, dry, lo, np.asarray(starts))
    ref = np.where(n >= min_dry, ref, np.nan)
    return ref, n


def rain_rate(rain_att, k, alpha, length_km, offset_db=0.0):
    """Path rain rate (mm/h) after removing a length-independent offset; negative values are set to zero."""
    a = np.clip(np.asarray(rain_att, dtype=float) - offset_db, 0.0, None)
    return (a / (k * length_km)) ** (1.0 / alpha)


def equivalent_rate(error_db, k, alpha, length_km):
    """Rain rate that would produce an attenuation of ``error_db`` on the path (mm/h)."""
    return (np.asarray(error_db, dtype=float) / (k * length_km)) ** (1.0 / alpha)


def held(n, starts, stops, values):
    """Array of length n holding values[i] on [starts[i], stops[i]) and NaN elsewhere."""
    out = np.full(n, np.nan)
    for s, e, v in zip(starts, stops, values):
        out[s:e] = v
    return out


def link_events(tsl, rsl, radar_min, threshold, margin, window, min_dry, min_event):
    """Wet/dry masks from radar, valid events and rain-period attenuation of one sub-link.

    Returns a dict with the attenuation ``att`` (TSL - RSL), the masks, the event bounds with their reference and
    reference size, and ``rain_att`` = att - reference inside the events (NaN elsewhere).
    """
    att = np.asarray(tsl, dtype=float) - np.asarray(rsl, dtype=float)
    wet, dry = wet_dry(radar_min, threshold, margin)
    return events_from_masks(att, wet, dry, window, min_dry, min_event)


def rolling_std_mask(att, window, threshold):
    """Wet where the centred rolling standard deviation of the attenuation exceeds ``threshold``.

    Schleiss and Berne (2010) type classification that uses the link signal only; returns (wet, dry, rolling sd).
    """
    sd = pd.Series(att).rolling(window, center=True, min_periods=window // 2).std().to_numpy()
    wet = sd > threshold
    dry = (sd <= threshold) & np.isfinite(att)
    return wet, dry, sd


def events_from_masks(att, wet, dry, window, min_dry, min_event):
    """Valid events (with their clear-sky reference) and rain-period attenuation for given wet/dry masks."""
    att = np.asarray(att, dtype=float)
    starts, stops = runs(wet)
    keep = (stops - starts) >= min_event
    starts, stops = starts[keep], stops[keep]
    ref, n_ref = event_reference(att, dry, starts, window, min_dry)
    valid = np.isfinite(ref)
    starts, stops, ref, n_ref = starts[valid], stops[valid], ref[valid], n_ref[valid]
    return dict(att=att, wet=wet, dry=dry, starts=starts, stops=stops, ref=ref, n_ref=n_ref,
                rain_att=att - held(len(att), starts, stops, ref))
