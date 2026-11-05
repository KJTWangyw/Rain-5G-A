"""Read the OpenMRG data set (Andersson et al., 2022; doi:10.5281/zenodo.7107689) into 1-min arrays.

CML: 10-s transmitted and received signal levels (TSL, RSL, dBm) of 728 sub-links in Gothenburg, June-August 2015.
Radar: SMHI composite, 2 km, 5 min, reflectivity in dBZ. Only the files below are used:

    cml/cml.nc, cml/cml_metadata.csv, radar/radar.nc
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

SAMPLES_PER_MINUTE = 6
MISSING = 1e9


def read_metadata(data_dir):
    meta = pd.read_csv(Path(data_dir) / "cml" / "cml_metadata.csv", encoding="utf-8-sig")
    return meta


def select_sublinks(meta, f_min, f_max):
    sel = meta[(meta.Frequency_GHz >= f_min) & (meta.Frequency_GHz <= f_max)]
    return sel.sort_values("Sublink").reset_index(drop=True)


def _minute_reduce(x, func):
    """Reduce (rows, links) 10-s samples to (rows // 6, links) with ``func`` over finite values (>= 3 required)."""
    xr = x.reshape(-1, SAMPLES_PER_MINUTE, x.shape[1])
    ok = np.isfinite(xr)
    n = ok.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        if func == "mean":
            out = np.where(ok, xr, 0.0).sum(axis=1) / n
        else:  # variance (population) within the minute
            m = np.where(ok, xr, 0.0).sum(axis=1) / n
            out = (np.where(ok, (xr - m[:, None, :]) ** 2, 0.0)).sum(axis=1) / n
    out[n < 3] = np.nan
    return out


def minute_aggregates(data_dir, sublinks, steps_db, offsets_db, block_rows=60_000):
    """1-min means of TSL and RSL, within-minute RSL variance, and 1-min means of RSL re-quantized at each step.

    ``offsets_db[step]`` holds one grid offset per selected sub-link; the 10-s RSL is rounded to
    offset + step * integer before averaging, as a receiver reporting at that resolution would do.
    """
    import netCDF4

    nc = netCDF4.Dataset(Path(data_dir) / "cml" / "cml.nc")
    t = np.asarray(nc["time"][:], dtype=np.int64)
    if not (np.all(np.diff(t) == 10) and t[0] % 60 == 0):
        raise ValueError("expected a regular 10-s time axis starting on a full minute")
    index = {s: i for i, s in enumerate(np.asarray(nc["sublink"][:]))}
    cols = np.array([index[s] for s in sublinks])
    T, nl = len(t), len(cols)
    n_min = -(-T // SAMPLES_PER_MINUTE)
    out = {"time": t[0] + 60 * np.arange(n_min, dtype=np.int64)}
    for name in ("tsl", "rsl", "rsl_var") + tuple(f"rsl_q{s:g}" for s in steps_db):
        out[name] = np.full((n_min, nl), np.nan, dtype=np.float32)
    assert block_rows % SAMPLES_PER_MINUTE == 0
    for s in range(0, T, block_rows):
        e = min(T, s + block_rows)
        pad = (-(e - s)) % SAMPLES_PER_MINUTE
        m0, m1 = s // SAMPLES_PER_MINUTE, (e + pad) // SAMPLES_PER_MINUTE
        for var in ("tsl", "rsl"):
            x = np.asarray(nc[var][s:e, :], dtype=np.float64)[:, cols]
            x[x > MISSING] = np.nan
            if pad:
                x = np.vstack([x, np.full((pad, nl), np.nan)])
            out[var][m0:m1] = _minute_reduce(x, "mean")
            if var == "rsl":
                out["rsl_var"][m0:m1] = _minute_reduce(x, "var")
                for step in steps_db:
                    q = offsets_db[step] + step * np.round((x - offsets_db[step]) / step)
                    out[f"rsl_q{step:g}"][m0:m1] = _minute_reduce(q, "mean")
    nc.close()
    return out


def path_pixels(lat, lon, near, far, n_points):
    """Radar cells under ``n_points`` equally spaced points of a path (repeated cells weight by path length)."""
    frac = np.linspace(0.0, 1.0, n_points)
    pts = [(near[0] + f * (far[0] - near[0]), near[1] + f * (far[1] - near[1])) for f in frac]
    out = []
    for la, lo in pts:
        d = (lat - la) ** 2 + ((lon - lo) * np.cos(np.radians(la))) ** 2
        out.append(np.unravel_index(np.argmin(d), d.shape))
    return np.array(out)


def radar_path_average(data_dir, meta, zr_a, zr_b, min_dbz, n_points, block=2000):
    """5-min radar rain rate (mm/h, Z = a R^b) averaged along every path; returns (times, R[time, link])."""
    import netCDF4

    nc = netCDF4.Dataset(Path(data_dir) / "radar" / "radar.nc")
    times = np.asarray(nc["time"][:], dtype=np.int64)
    lat, lon = np.asarray(nc["lat"][:]), np.asarray(nc["lon"][:])
    pix = [path_pixels(lat, lon, (r.NearLatitude_DecDeg, r.NearLongitude_DecDeg),
                       (r.FarLatitude_DecDeg, r.FarLongitude_DecDeg), n_points) for r in meta.itertuples()]
    R = np.full((len(times), len(meta)), np.nan, dtype=np.float32)
    for s in range(0, len(times), block):
        e = min(len(times), s + block)
        dbz = np.ma.filled(nc["data"][s:e].astype(np.float64), np.nan)
        rr = np.where(np.isfinite(dbz), (10.0 ** (dbz / 10.0) / zr_a) ** (1.0 / zr_b), np.nan)
        rr[dbz < min_dbz] = 0.0
        for j, p in enumerate(pix):
            with np.errstate(invalid="ignore"):
                R[s:e, j] = np.nanmean(rr[:, p[:, 0], p[:, 1]], axis=1)
    nc.close()
    return times, R
