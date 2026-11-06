"""Read the OpenRainER data set (Covi and Roversi; doi:10.5281/zenodo.10593848) into 1-min arrays.

CML: monthly NetCDF files with 1-min RSL and TSL of the Lepida network (Emilia-Romagna, Italy), dimensions
(cml_id, sublink_id, time). The data provider rounded the levels to 1 dB and averaged them before logging; the
data set README recommends rounding to the nearest dBm, which is applied here.
Radar: 15-min radar rainfall accumulations (RADrain, about 1 km grid). Gauges: 15-min accumulations (AWS).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

MISSING = -9999.0


def month_list(first, last):
    """['YYYYMM', ...] from first to last inclusive."""
    p = pd.period_range(pd.Period(first, "M"), pd.Period(last, "M"), freq="M")
    return [f"{x.year}{x.month:02d}" for x in p]


def cml_files(data_dir):
    return {f.name[4:10]: f for f in sorted((Path(data_dir) / "cml").glob("CML_*.nc"))}


def read_cml(data_dir, months, round_db=1.0):
    """Sub-link metadata and 1-min TSL/RSL arrays on a regular grid covering ``months`` (NaN where missing)."""
    import netCDF4

    files = cml_files(data_dir)
    months = [m for m in months if m in files]
    first = netCDF4.Dataset(files[months[0]])
    ids = [str(x) for x in first["cml_id"][:]]
    subs = [str(x) for x in first["sublink_id"][:]]
    freq = np.asarray(first["frequency"][:], dtype=float) / 1000.0
    pol = np.asarray(first["polarization"][:])
    rows = []
    for i, cid in enumerate(ids):
        for s, sid in enumerate(subs):
            rows.append(dict(Sublink=f"{cid}_{sid}", Link=cid, Frequency_GHz=freq[i, s],
                             Polarization=str(pol[i, s]), Length_km=float(first["length"][i]) / 1000.0,
                             NearLatitude_DecDeg=float(first["site_0_lat"][i]), NearLongitude_DecDeg=float(first["site_0_lon"][i]),
                             FarLatitude_DecDeg=float(first["site_1_lat"][i]), FarLongitude_DecDeg=float(first["site_1_lon"][i])))
    first.close()
    meta = pd.DataFrame(rows)
    t0 = int(pd.Timestamp(f"{months[0][:4]}-{months[0][4:]}-01", tz="UTC").timestamp())
    end = pd.Period(months[-1], "M").end_time.tz_localize("UTC").ceil("min")
    n_min = (int(end.timestamp()) - t0) // 60
    out = {"time": t0 + 60 * np.arange(n_min, dtype=np.int64)}
    for v in ("tsl", "rsl"):
        out[v] = np.full((n_min, len(meta)), np.nan, dtype=np.float32)
    position = {cid: i for i, cid in enumerate(ids)}
    for m in months:
        nc = netCDF4.Dataset(files[m])
        these = [str(x) for x in nc["cml_id"][:]]
        if [str(x) for x in nc["sublink_id"][:]] != subs or set(these) != set(ids):
            raise ValueError(f"links or sub-links differ in {files[m].name}")
        rows_here = np.array([position[c] for c in these])           # file row -> master row
        t = np.asarray(nc["time"][:], dtype=np.int64)
        idx = (t - t0) // 60
        ok = (idx >= 0) & (idx < n_min)
        for v in ("tsl", "rsl"):
            x = np.asarray(nc[v][:], dtype=float)                     # (cml, sub, time) in file order
            x[x <= MISSING + 1] = np.nan
            if round_db:
                x = round_db * np.round(x / round_db)
            y = np.empty_like(x)
            y[rows_here] = x                                          # master order
            out[v][idx[ok]] = y.reshape(-1, y.shape[-1]).T[ok].astype(np.float32)
        nc.close()
    return meta, out


def read_aws(data_dir, months):
    """15-min gauge accumulations (mm): times of the interval ends, values (time, station) and station metadata."""
    import gzip
    import netCDF4

    times, vals, meta = [], [], None
    for m in months:
        f = Path(data_dir) / "aws" / f"AWS_{m}.nc"
        if not f.exists():
            gz = f.with_suffix(".nc.gz")
            f.write_bytes(gzip.decompress(gz.read_bytes()))
        nc = netCDF4.Dataset(f)
        ids = [str(x) for x in nc["id"][:]]
        this = pd.DataFrame(dict(id=ids, lat=np.asarray(nc["latitude"][:]), lon=np.asarray(nc["longitude"][:])))
        if meta is None:
            meta = this
        x = np.asarray(nc["rainfall_amount"][:], dtype=float)
        x[x <= MISSING + 1] = np.nan
        frame = pd.DataFrame(x.T, columns=ids).reindex(columns=meta.id)
        times.append(np.asarray(nc["time"][:], dtype=np.int64))
        vals.append(frame.values)
        nc.close()
    return np.concatenate(times), np.vstack(vals), meta


def path_indices(lat_grid, lon_grid, meta, n_points):
    """Nearest (lat, lon) grid indices of ``n_points`` equally spaced points on every path."""
    frac = np.linspace(0.0, 1.0, n_points)
    out = []
    for r in meta.itertuples():
        la = r.NearLatitude_DecDeg + frac * (r.FarLatitude_DecDeg - r.NearLatitude_DecDeg)
        lo = r.NearLongitude_DecDeg + frac * (r.FarLongitude_DecDeg - r.NearLongitude_DecDeg)
        iy = np.abs(lat_grid[None, :] - la[:, None]).argmin(axis=1)
        ix = np.abs(lon_grid[None, :] - lo[:, None]).argmin(axis=1)
        out.append((iy, ix))
    return out


def radar_path_average(data_dir, meta, months, n_points, work_dir, chunk=192):
    """15-min radar accumulations along every path, as rain rate (mm/h): (times, R[time, link]).

    Monthly files are decompressed one at a time into ``work_dir`` and removed after use.
    """
    import gzip
    import shutil
    import tarfile
    import netCDF4

    data_dir, work_dir = Path(data_dir), Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    times, rates = [], []
    with tarfile.open(data_dir / "RADrain.tar") as tar:
        for m in months:
            gz = work_dir / f"RADrain_{m}.nc.gz"
            nc_path = work_dir / f"RADrain_{m}.nc"
            with tar.extractfile(f"RADrain_{m}.nc.gz") as src, open(gz, "wb") as dst:
                shutil.copyfileobj(src, dst)
            with gzip.open(gz, "rb") as src, open(nc_path, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 24)
            gz.unlink()
            nc = netCDF4.Dataset(nc_path)
            lat, lon = np.asarray(nc["lat"][:]), np.asarray(nc["lon"][:])
            pix = path_indices(lat, lon, meta, n_points)
            iy = np.concatenate([p[0] for p in pix])
            ix = np.concatenate([p[1] for p in pix])
            t = np.asarray(nc["time"][:], dtype=np.int64)
            var = nc["rainfall_amount"]
            var.set_auto_mask(False)
            R = np.full((len(t), len(meta)), np.nan, dtype=np.float32)
            for s in range(0, len(t), chunk):
                block = np.asarray(var[s:s + chunk], dtype=np.float32)          # (time, lat, lon)
                vals = block[:, iy, ix]
                vals[vals < 0] = np.nan
                vals = vals.reshape(len(block), len(meta), n_points)
                with np.errstate(invalid="ignore"):
                    R[s:s + chunk] = 4.0 * np.nanmean(vals, axis=2)           # mm per 15 min -> mm/h
            nc.close()
            nc_path.unlink()
            times.append(t)
            rates.append(R)
            print("radar", m, len(t), flush=True)
    return np.concatenate(times), np.vstack(rates)
