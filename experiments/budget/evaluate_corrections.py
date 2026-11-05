#!/usr/bin/env python3
"""Robustness of the offset correction and physical check of the offset (results/budget/).

  evaluation_ext.csv  August hourly comparison per sub-link for
                      mask   radar (as in analyze_openmrg.py) or link (rolling standard deviation of the attenuation,
                             threshold set in the calibration months; no radar in August)
                      corr   none, constant, rain (W(R) from length diversity) or kralt (Pastorek et al. KR-alt
                             model with literature parameters)
                      ref    radar path average, or the nearest city gauge for sub-links within ``gauge_km``
                             (the radar itself is also scored against the gauges)
  drying.csv          attenuation after the radar rain on the path has ended (event reference held), per sub-link
                      and 5-min lag bin: the wet-antenna drying signature, compared across link lengths

Usage:  python experiments/budget/evaluate_corrections.py [--config experiments/budget/config.json]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments" / "budget"))
import analyze_openmrg as A  # noqa: E402
from cmlbudget.offset import invert, offset_function  # noqa: E402
from cmlbudget.retrieval import events_from_masks, rain_rate, rolling_std_mask, wet_dry  # noqa: E402


def kralt(params):
    """Pastorek et al. (2022) KR-alt wet-antenna model WAA(R) = A_max (1 - exp(-d R^zeta)), zero at R = 0."""
    a_max, zeta, d = params["A_max"], params["zeta"], params["d"]

    def W(R):
        R = np.asarray(R, dtype=float)
        return a_max * (1.0 - np.exp(-d * np.clip(R, 0, None) ** zeta))
    return W


def segment_distance_km(lat, lon, near, far):
    """Distance (km) from a point to a path segment, equirectangular approximation."""
    kx, ky = 111.32 * np.cos(np.radians(lat)), 110.57
    p = np.array([lon * kx, lat * ky])
    a = np.array([near[1] * kx, near[0] * ky])
    b = np.array([far[1] * kx, far[0] * ky])
    t = np.clip(np.dot(p - a, b - a) / max(np.dot(b - a, b - a), 1e-12), 0.0, 1.0)
    return float(np.linalg.norm(p - (a + t * (b - a))))


def city_gauges(cfg, minute_time):
    """1-min city-gauge rainfall (mm) on the cache grid (interval start) and the gauge metadata."""
    data = ROOT / cfg["data_dir"]
    g = pd.read_csv(data / "gauges/city/CityGauges-2015JJA.csv")
    meta = pd.read_csv(data / "gauges/city/CityGauges-metadata.csv", encoding="utf-8-sig")
    end = pd.to_datetime(g.Time_UTC, utc=True)
    start = ((end - pd.Timestamp("1970-01-01", tz="UTC")) // pd.Timedelta(seconds=1)).to_numpy() - 60
    idx = (start - minute_time[0]) // 60
    ok = (idx >= 0) & (idx < len(minute_time))
    out = np.full((len(minute_time), len(meta)), np.nan)
    out[idx[ok]] = g[meta.Name.values].values[ok]
    return out, meta


def scores(hour, est, ref, ok, hmin, ref_hourly_sum=False):
    """Hourly scores of an estimate (mm/h per minute) against a reference (mm/h per minute, or mm per minute)."""
    d = pd.DataFrame({"h": hour[ok], "x": est[ok], "y": ref[ok]}).groupby("h")
    n = d.size()
    keep = n >= hmin
    x = d.x.mean()[keep].values
    y = (d.y.sum() if ref_hourly_sum else d.y.mean())[keep].values
    if len(x) < 24:
        return None
    wet = y >= 0.1
    return dict(hours=len(x), acc_est=x.sum(), acc_ref=y.sum(), rmse=np.sqrt(np.mean((x - y) ** 2)),
                r=np.corrcoef(x, y)[0, 1] if x.std() > 0 and y.std() > 0 else np.nan,
                pod=np.mean(x[wet] >= 0.1) if wet.any() else np.nan,
                far=np.mean(y[x >= 0.1] < 0.1) if (x >= 0.1).any() else np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "experiments" / "budget" / "config.json")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "budget")
    args = ap.parse_args()
    cfg = json.loads(args.config.read_text())
    ext = cfg["extension"]
    t0 = time.time()
    meta, arr = A.load(cfg)
    links = pd.read_csv(args.out / "links.csv")
    usable = set(links.j[links.wet_min >= cfg["min_wet_minutes"]])
    curves = pd.read_csv(args.out / "offset_curve.csv")
    cal_curve = curves[(curves.variant == "native") & (curves.period == "cal")].sort_values("bin")
    W_rain = offset_function(cal_curve.R.values, cal_curve.W.values)
    W_const = json.loads((args.out / "offset.json").read_text())["W_cal"]["native"]
    W_kr = kralt(ext["kralt"])
    tmin = arr["time"]
    month = pd.to_datetime(tmin, unit="s").month.values
    cal = np.isin(month, cfg["calibration_months"])
    evm = np.isin(month, cfg["evaluation_months"])
    hour = (tmin // 3600) - tmin[0] // 3600
    gauge, gmeta = city_gauges(cfg, tmin)
    thr, margin = cfg["wet"]["threshold_mmh"], cfg["wet"]["margin_min"]
    window, min_dry, min_event = cfg["reference"]["window_min"], cfg["reference"]["min_dry_min"], cfg["event_min_minutes"]
    lag_edges = np.arange(0, ext["drying_max_min"] + 1, ext["drying_bin_min"])

    rows, dry_rows, mask_rows = [], [], []
    for j, m in meta.iterrows():
        if j not in usable:
            continue
        tsl, rsl = arr["tsl"][:, j].astype(float), arr["rsl"][:, j].astype(float)
        att = tsl - rsl
        radar = arr["radar_min"][:, j]
        wet_r, dry_r = wet_dry(radar, thr, margin)
        masks = {"radar": (wet_r, dry_r)}
        # link-only classification: threshold so that the wet fraction in the calibration months matches radar
        _, _, sd = rolling_std_mask(att, ext["std_window_min"], np.inf)
        cal_ok = cal & np.isfinite(sd) & np.isfinite(radar)
        f_wet = np.mean(wet_r[cal_ok]) if cal_ok.any() else np.nan
        if np.isfinite(f_wet) and 0 < f_wet < 1:
            sd_thr = np.quantile(sd[cal_ok], 1.0 - f_wet)
            wet_l, dry_l, _ = rolling_std_mask(att, ext["std_window_min"], sd_thr)
            masks["link"] = (wet_l, dry_l)
            agree = evm & np.isfinite(sd) & np.isfinite(radar)
            mask_rows.append(dict(j=j, f_wet_cal=f_wet, sd_thr=sd_thr,
                                  hit=np.mean(wet_l[agree & wet_r]) if (agree & wet_r).any() else np.nan,
                                  false=np.mean(wet_l[agree & ~wet_r]) if (agree & ~wet_r).any() else np.nan))
        # nearest city gauge
        dist = [segment_distance_km(g.Latitude_DecDeg, g.Longitude_DecDeg, (m.NearLatitude_DecDeg, m.NearLongitude_DecDeg),
                                    (m.FarLatitude_DecDeg, m.FarLongitude_DecDeg)) for g in gmeta.itertuples()]
        gi = int(np.argmin(dist))
        near_gauge = dist[gi] <= ext["gauge_km"]
        for mname, (wet, dry) in masks.items():
            ev = events_from_masks(att, wet, dry, window, min_dry, min_event)
            a = ev["rain_att"]
            inside = np.isfinite(a)
            for corr in ("none", "constant", "rain", "kralt"):
                rate = np.where(dry, 0.0, np.nan)
                if corr in ("none", "constant"):
                    rate[inside] = rain_rate(a[inside], m.k, m.alpha, m.Length_km, 0.0 if corr == "none" else W_const)
                else:
                    rate[inside] = invert(a[inside], m.kL, m.alpha, W_rain if corr == "rain" else W_kr)
                ok = evm & np.isfinite(rate) & np.isfinite(radar)
                s = scores(hour, rate, radar, ok, cfg["hour_min_minutes"])
                if s:
                    rows.append(dict(j=j, Link=m.Link, mask=mname, corr=corr, ref="radar", **s))
                if near_gauge and mname == "radar":
                    gk = gauge[:, gi]
                    okg = evm & np.isfinite(rate) & np.isfinite(gk)
                    s = scores(hour, rate, gk, okg, ext["gauge_hour_min_minutes"], ref_hourly_sum=True)
                    if s:
                        rows.append(dict(j=j, Link=m.Link, mask=mname, corr=corr, ref="gauge", gauge=gmeta.Name[gi],
                                         gauge_km=dist[gi], **s))
            if mname == "radar":
                # drying: attenuation after the last rainy minute of each event, while the path stays dry on radar
                rainy = np.nan_to_num(radar, nan=np.inf) > thr
                ref_held = att - a  # reference inside events
                for s0, e0, r0 in zip(ev["starts"], ev["stops"], ev["ref"]):
                    idx = np.flatnonzero(rainy[s0:e0])
                    if not len(idx):
                        continue
                    t_end = s0 + idx[-1]
                    seg = slice(t_end, min(t_end + ext["drying_max_min"] + 1, len(att)))
                    later_rain = np.cumsum(rainy[seg]) - rainy[t_end]   # rain after t_end within the window
                    valid = (later_rain == 0) & np.isfinite(att[seg])
                    lags = np.arange(seg.stop - seg.start)[valid]
                    vals = att[seg][valid] - r0
                    b = np.digitize(lags, lag_edges) - 1
                    for k in np.unique(b):
                        if 0 <= k < len(lag_edges) - 1:
                            dry_rows.append((j, k, vals[b == k].sum(), int((b == k).sum())))
                del ref_held
        # report resolution after the wet-antenna correction (radar mask, references as above)
        for step in cfg["quantization_steps_db"]:
            att_q = tsl - arr[f"rsl_q{step:g}"][:, j].astype(float)
            evq = events_from_masks(att_q, wet_r, dry_r, window, 1, min_event)
            keep = np.isin(evq["starts"], events_from_masks(att, wet_r, dry_r, window, min_dry, min_event)["starts"])
            aq = np.full(len(att_q), np.nan)
            for s0, e0, r0, kp in zip(evq["starts"], evq["stops"], evq["ref"], keep):
                if kp:
                    aq[s0:e0] = att_q[s0:e0] - r0
            inside = np.isfinite(aq)
            for corr in ("none", "kralt"):
                rate = np.where(dry_r, 0.0, np.nan)
                rate[inside] = (rain_rate(aq[inside], m.k, m.alpha, m.Length_km, 0.0) if corr == "none"
                                else invert(aq[inside], m.kL, m.alpha, W_kr))
                ok = evm & np.isfinite(rate) & np.isfinite(radar)
                s = scores(hour, rate, radar, ok, cfg["hour_min_minutes"])
                if s:
                    rows.append(dict(j=j, Link=m.Link, mask=f"radar_q{step:g}", corr=corr, ref="radar", **s))
                if near_gauge:
                    gk = gauge[:, gi]
                    okg = evm & np.isfinite(rate) & np.isfinite(gk)
                    s = scores(hour, rate, gk, okg, ext["gauge_hour_min_minutes"], ref_hourly_sum=True)
                    if s:
                        rows.append(dict(j=j, Link=m.Link, mask=f"radar_q{step:g}", corr=corr, ref="gauge",
                                         gauge=gmeta.Name[gi], gauge_km=dist[gi], **s))
        if near_gauge:
            gk = gauge[:, gi]
            okg = evm & np.isfinite(radar) & np.isfinite(gk)
            s = scores(hour, radar, gk, okg, ext["gauge_hour_min_minutes"], ref_hourly_sum=True)
            if s:
                rows.append(dict(j=j, Link=m.Link, mask="-", corr="radar", ref="gauge", gauge=gmeta.Name[gi],
                                 gauge_km=dist[gi], **s))
        if j % 100 == 0:
            print(f"{j}/{len(meta)} ({time.time() - t0:.0f} s)", flush=True)

    ev_df = pd.DataFrame(rows)
    ev_df.to_csv(args.out / "evaluation_ext.csv", index=False)
    d = pd.DataFrame(dry_rows, columns=["j", "lag_bin", "sum", "n"]).groupby(["j", "lag_bin"]).sum().reset_index()
    d["mean"] = d["sum"] / d["n"]
    d["lag_min"] = lag_edges[d.lag_bin] + ext["drying_bin_min"] / 2.0
    d.drop(columns="sum").to_csv(args.out / "drying.csv", index=False)
    pd.DataFrame(mask_rows).to_csv(args.out / "link_mask.csv", index=False)
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
