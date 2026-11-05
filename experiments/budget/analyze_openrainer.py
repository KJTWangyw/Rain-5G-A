#!/usr/bin/env python3
"""Second region: the same error-budget checks on OpenRainER (Emilia-Romagna, 25-26 GHz, 1-dB levels).

Seasons: May-September 2021 (calibration of the constant and binned offsets) and May-September 2022 (evaluation).
Outputs (results/budget/openrainer/): timing.json (radar and gauge time-stamp conventions found from the CML data),
links.csv, offset_curve.csv, offset.json, drying.csv, dry_lag.csv, evaluation.csv.

Usage:  python experiments/budget/analyze_openrainer.py [--config experiments/budget/config.json]
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
from cmlbudget.offset import (bin_means, cluster_bootstrap, fit_offset, invert, offset_curve,  # noqa: E402
                              offset_function)
from cmlbudget.openrainer import month_list  # noqa: E402
from cmlbudget.retrieval import held, link_events, rain_rate, wet_dry, window_stats  # noqa: E402
from evaluate_corrections import kralt, scores, segment_distance_km  # noqa: E402
from nrpilot.rain import p838  # noqa: E402

STEP = 900  # 15-min radar and gauge accumulations


def interval_index(stamps, minute_time, convention):
    """Index of the 15-min accumulation that contains each minute ('end': stamp closes the interval)."""
    if convention == "end":
        idx = np.searchsorted(stamps, minute_time, side="right")
        ok = (idx < len(stamps))
        ok[ok] &= (stamps[idx[ok]] - minute_time[ok]) <= STEP
    else:
        idx = np.searchsorted(stamps, minute_time, side="right") - 1
        ok = idx >= 0
        ok[ok] &= (minute_time[ok] - stamps[idx[ok]]) < STEP
    return np.clip(idx, 0, len(stamps) - 1), ok


def convention_score(stamps, values, minute_time, att, convention, max_links=80):
    """Median correlation between 15-min mean attenuation and the 15-min series under a time-stamp convention."""
    idx, ok = interval_index(stamps, minute_time, convention)
    cors = []
    for j in range(min(att.shape[1], max_links)):
        a = att[:, j]
        good = ok & np.isfinite(a)
        if good.sum() < 10000:
            continue
        a_b = pd.Series(a[good]).groupby(idx[good]).mean()
        v = values[a_b.index.values, j] if values.ndim == 2 else values[a_b.index.values]
        m = np.isfinite(v) & np.isfinite(a_b.values)
        if m.sum() > 500 and np.nanstd(v[m]) > 0:
            cors.append(np.corrcoef(a_b.values[m], v[m])[0, 1])
    return float(np.median(cors)) if cors else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "experiments" / "budget" / "config.json")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "budget" / "openrainer")
    args = ap.parse_args()
    cfg = json.loads(args.config.read_text())
    oc = cfg["openrainer"]
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    cache = ROOT / oc["cache"]
    z = np.load(cache)
    arr = {k: z[k] for k in ("time", "tsl", "rsl", "radar_time", "radar_R", "gauge_time", "gauge")}
    meta = pd.read_csv(cache.with_name("openrainer_sublinks.csv"))
    gmeta = pd.read_csv(cache.with_name("openrainer_gauges.csv"))
    k, alpha = p838(meta.Frequency_GHz.values, meta.Polarization.values)
    meta = meta.assign(k=k, alpha=alpha, kL=k * meta.Length_km.values)
    tmin = arr["time"]
    month = pd.to_datetime(tmin, unit="s").strftime("%Y%m").values
    cal_m = month_list(*oc["calibration_months"])
    ev_m = month_list(*oc["evaluation_months"])
    cal, evm = np.isin(month, cal_m), np.isin(month, ev_m)
    att_all = arr["tsl"].astype(float) - arr["rsl"].astype(float)

    # ---- time-stamp conventions, judged by the correlation with the 1-min CML attenuation
    rt, rR = arr["radar_time"], arr["radar_R"].astype(float)
    order = np.argsort(rt)
    rt, rR = rt[order], rR[order]
    conv_r = {c: convention_score(rt, rR, tmin, att_all, c) for c in ("end", "start")}
    radar_conv = max(conv_r, key=conv_r.get)
    ridx, rok = interval_index(rt, tmin, radar_conv)
    radar_min = rR[ridx]
    radar_min[~rok] = np.nan
    gt, g = arr["gauge_time"], arr["gauge"].astype(float)
    order = np.argsort(gt)
    gt, g = gt[order], g[order]
    # nearest gauge (path-segment distance) for each sub-link
    near = []
    for m in meta.itertuples():
        d = [segment_distance_km(r.lat, r.lon, (m.NearLatitude_DecDeg, m.NearLongitude_DecDeg),
                                 (m.FarLatitude_DecDeg, m.FarLongitude_DecDeg)) for r in gmeta.itertuples()]
        near.append((int(np.argmin(d)), float(np.min(d))))
    gauge_series = np.stack([g[:, i] * 4.0 for i, _ in near], axis=1)       # mm/h per 15 min, per sub-link
    conv_g = {c: convention_score(gt, gauge_series, tmin, att_all, c) for c in ("end", "start")}
    gauge_conv = max(conv_g, key=conv_g.get)
    gidx, gok = interval_index(gt, tmin, gauge_conv)
    (args.out / "timing.json").write_text(json.dumps(dict(radar=conv_r, radar_convention=radar_conv,
                                                          gauge=conv_g, gauge_convention=gauge_conv), indent=1))
    print("conventions", radar_conv, conv_r, gauge_conv, conv_g, flush=True)

    thr, margin = cfg["wet"]["threshold_mmh"], cfg["wet"]["margin_min"]
    window, min_dry, min_event = cfg["reference"]["window_min"], cfg["reference"]["min_dry_min"], cfg["event_min_minutes"]
    bins = np.asarray(cfg["rain_bins_mmh"], dtype=float)
    W_kr = kralt(cfg["extension"]["kralt"])
    hour = (tmin // 3600) - tmin[0] // 3600
    lag_edges = np.arange(0, oc["drying_max_min"] + 1, oc["drying_bin_min"])
    link_rows, bin_rows, dry_rows, lag_rows, rain_att = [], [], [], [], {}
    for j, m in meta.iterrows():
        radar = radar_min[:, j]
        ev = link_events(arr["tsl"][:, j], arr["rsl"][:, j], radar, thr, margin, window, min_dry, min_event)
        a, att = ev["rain_att"], ev["att"]
        rain_att[j] = (a.astype(np.float32), ev["dry"])
        expected = m.kL * np.where(radar > 0, radar, 0.0) ** m.alpha
        use = np.isfinite(a) & np.isfinite(expected)
        c = use & cal
        link_rows.append(dict(j=j, Sublink=m.Sublink, Link=m.Link, length_km=m.Length_km, freq_ghz=m.Frequency_GHz,
                              pol=m.Polarization, k=m.k, alpha=m.alpha, kL=m.kL, events=len(ev["starts"]),
                              wet_min=int(use.sum()), att_mean=np.mean(a[use]) if use.any() else np.nan,
                              exp_mean=np.mean(expected[use]) if use.any() else np.nan, wet_min_cal=int(c.sum()),
                              att_mean_cal=np.mean(a[c]) if c.any() else np.nan,
                              exp_mean_cal=np.mean(expected[c]) if c.any() else np.nan,
                              gauge=gmeta.id[near[j][0]], gauge_km=near[j][1]))
        for pname, pmask in (("full", np.ones(len(a), bool)), ("cal", cal)):
            for r in bin_means(np.where(pmask, a, np.nan), expected, radar, bins, cfg["bin_min_minutes"]):
                bin_rows.append(dict(j=j, Link=m.Link, kL=m.kL, period=pname, **r))
        # drying after the last rainy radar interval of each event
        rainy = np.nan_to_num(radar, nan=np.inf) > thr
        for s0, e0, r0 in zip(ev["starts"], ev["stops"], ev["ref"]):
            idx = np.flatnonzero(rainy[s0:e0])
            if not len(idx):
                continue
            t_end = s0 + idx[-1]
            seg = slice(t_end, min(t_end + oc["drying_max_min"] + 1, len(att)))
            later = np.cumsum(rainy[seg]) - rainy[t_end]
            valid = (later == 0) & np.isfinite(att[seg])
            lags = np.arange(seg.stop - seg.start)[valid]
            vals = att[seg][valid] - r0
            b = np.digitize(lags, lag_edges) - 1
            for kb in np.unique(b):
                if 0 <= kb < len(lag_edges) - 1:
                    dry_rows.append((j, kb, vals[b == kb].sum(), int((b == kb).sum())))
        # dry-reference test (1 h and 6 h after the end of the 6-h reference)
        dry = ev["dry"]
        full = pd.Series(dry & np.isfinite(att)).groupby(hour).sum().values
        hour_start = np.searchsorted(hour, np.arange(hour[-1] + 1))
        h_full = hour_start[np.flatnonzero(full == 60)]
        for lag in (1, 6):
            ref_end = h_full - (lag - 1) * 60
            ok = ref_end - window >= 0
            h0, ref_end = h_full[ok], ref_end[ok]
            n_r, r_mean = window_stats(att, dry, ref_end - window, ref_end)
            _, h_mean = window_stats(att, dry, h0, h0 + 60)
            good = n_r >= min_dry
            err = (h_mean - r_mean)[good]
            lag_rows.append(dict(j=j, lag_h=lag, n=int(good.sum()), sd_h=np.std(err) if good.sum() > 1 else np.nan))
        if j % 100 == 0:
            print(f"pass 1: {j}/{len(meta)} ({time.time() - t0:.0f} s)", flush=True)
    links = pd.DataFrame(link_rows)
    binned = pd.DataFrame(bin_rows)

    usable = links[links.wet_min >= cfg["min_wet_minutes"]].copy()
    stat = lambda d: fit_offset(d.att_mean.values, d.exp_mean.values)  # noqa: E731
    w, beta = stat(usable)
    bs = cfg["bootstrap"]
    boot = cluster_bootstrap(usable.reset_index(drop=True), stat, bs["draws"], bs["seed"])
    rf = lambda d: fit_offset(d.att_mean.values, d.kL.values)  # noqa: E731
    w_rf = rf(usable)[0]
    usable_cal = links[links.wet_min_cal >= cfg["min_wet_minutes"]]
    w_cal = fit_offset(usable_cal.att_mean_cal.values, usable_cal.exp_mean_cal.values)[0]
    curves = []
    for pname, d in binned[binned.j.isin(usable.j)].groupby("period"):
        cc = offset_curve(d).assign(period=pname)
        if pname == "full":
            d = d.reset_index(drop=True)
            bc = cluster_bootstrap(d, lambda x: offset_curve(x).set_index("bin").W.reindex(cc.bin).values,
                                   bs["draws"], bs["seed"] + 9)
            cc["W_lo"], cc["W_hi"] = np.nanpercentile(bc, 2.5, axis=0), np.nanpercentile(bc, 97.5, axis=0)
        curves.append(cc)
    curves = pd.concat(curves, ignore_index=True)
    usable["offset_resid"] = usable.att_mean - (w + beta * usable.exp_mean)
    offset = dict(links=len(usable), phys=int(usable.Link.nunique()), W=w, beta=beta,
                  W_ci=np.percentile(boot[:, 0], [2.5, 97.5]).tolist(), W_radar_free=w_rf, W_cal=w_cal,
                  resid_sd=float(np.std(usable.offset_resid)))

    # ---- evaluation season: hourly scores against radar and nearby gauges
    cal_curve = curves[curves.period == "cal"].sort_values("bin")
    W_rain = offset_function(cal_curve.R.values, cal_curve.W.values)
    rows = []
    gauge_min = np.full(len(tmin), np.nan)
    for j, m in meta.iterrows():
        if j not in set(usable.j):
            continue
        a, dry = rain_att[j]
        a = a.astype(float)
        radar = radar_min[:, j]
        inside = np.isfinite(a)
        gi, gkm = near[j]
        gauge_min[:] = np.nan
        gauge_min[gok] = g[gidx[gok], gi] * 4.0          # mm/h, constant within each 15-min interval
        for corr in ("none", "constant", "rain", "kralt"):
            rate = np.where(dry, 0.0, np.nan)
            if corr in ("none", "constant"):
                rate[inside] = rain_rate(a[inside], m.k, m.alpha, m.Length_km, 0.0 if corr == "none" else w_cal)
            else:
                rate[inside] = invert(a[inside], m.kL, m.alpha, W_rain if corr == "rain" else W_kr)
            ok = evm & np.isfinite(rate) & np.isfinite(radar)
            s = scores(hour, rate, radar, ok, cfg["hour_min_minutes"])
            if s:
                rows.append(dict(j=j, Link=m.Link, corr=corr, ref="radar", **s))
            if gkm <= oc["gauge_km"]:
                okg = evm & np.isfinite(rate) & np.isfinite(gauge_min)
                s = scores(hour, rate, gauge_min, okg, 60)
                if s:
                    rows.append(dict(j=j, Link=m.Link, corr=corr, ref="gauge", gauge_km=gkm, **s))
        if gkm <= oc["gauge_km"]:
            okg = evm & np.isfinite(radar) & np.isfinite(gauge_min)
            s = scores(hour, radar, gauge_min, okg, 60)
            if s:
                rows.append(dict(j=j, Link=m.Link, corr="radar", ref="gauge", gauge_km=gkm, **s))
    evaluation = pd.DataFrame(rows)

    links.to_csv(args.out / "links.csv", index=False)
    curves.to_csv(args.out / "offset_curve.csv", index=False)
    (args.out / "offset.json").write_text(json.dumps(offset, indent=1, default=float))
    d = pd.DataFrame(dry_rows, columns=["j", "lag_bin", "sum", "n"]).groupby(["j", "lag_bin"]).sum().reset_index()
    d["mean"] = d["sum"] / d["n"]
    d["lag_min"] = lag_edges[d.lag_bin] + oc["drying_bin_min"] / 2.0
    d.drop(columns="sum").to_csv(args.out / "drying.csv", index=False)
    pd.DataFrame(lag_rows).to_csv(args.out / "dry_lag.csv", index=False)
    evaluation.to_csv(args.out / "evaluation.csv", index=False)
    print(f"done in {time.time() - t0:.0f} s: {len(usable)} sub-links, W = {w:.3f} dB")


if __name__ == "__main__":
    main()
