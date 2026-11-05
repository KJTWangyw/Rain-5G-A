#!/usr/bin/env python3
"""Error budget of rainfall estimation on the 24-40 GHz OpenMRG sub-links.

Every wet event uses the mean attenuation of the dry minutes in the 6 h before it as its clear-sky reference.
Reports at 1 dB and 2 dB are emulated by rounding the 10-s RSL before averaging (prepare_openmrg.py).

Outputs (results/budget/):
  links.csv          per sub-link: geometry, wet-period mean attenuation and radar-expected attenuation
  bins.csv           per sub-link, rain-rate bin, period and report resolution: mean attenuation, radar-expected
                     attenuation and radar rain rate
  offset_curve.csv   length-independent offset W per rain-rate bin (regression across sub-links)
  offset.json        constant offset with cluster-bootstrap interval, frequency bands, radar-free estimate,
                     calibration-period estimate and sensitivity to the wet margin and reference window
  events.csv         per wet event: re-quantization error of the reference and of the rain period, predicted values
  dry_lag.csv        error of a held reference over a dry hour versus the age of the reference
  evaluation.csv     hourly comparison with radar in the evaluation month, per sub-link, report resolution and
                     offset correction (none, constant, rain-dependent; both calibrated in the calibration months)

Usage:  python experiments/budget/analyze_openmrg.py [--config experiments/budget/config.json]
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
from cmlbudget import openmrg  # noqa: E402
from cmlbudget.offset import (bin_means, cluster_bootstrap, fit_offset, invert, offset_curve,  # noqa: E402
                              offset_function)
from cmlbudget.quantization import dither_bias  # noqa: E402
from cmlbudget.retrieval import event_reference, held, link_events, radar_on_minutes, rain_rate, window_stats  # noqa: E402
from nrpilot.rain import p838  # noqa: E402

LAGS_H = (1, 2, 4, 6)
BANDS = {"26-30": (24.0, 30.5), "32-34": (30.5, 35.0), "37-40": (35.0, 40.5)}


def load(cfg):
    """Cache arrays, sub-link metadata with P.838-3 coefficients, and radar rain rate on the 1-min grid."""
    z = np.load(ROOT / cfg["cache"])
    steps = cfg["quantization_steps_db"]
    keys = ["time", "tsl", "rsl", "rsl_var", "radar_time", "radar_R", "sublink"]
    keys += [f"rsl_q{s:g}" for s in steps] + [f"offset_q{s:g}" for s in steps]
    arr = {k: z[k] for k in keys}
    meta = openmrg.select_sublinks(openmrg.read_metadata(ROOT / cfg["data_dir"]), *cfg["frequency_ghz"])
    if not np.array_equal(meta.Sublink.values, arr["sublink"]):
        raise ValueError("cache and metadata disagree; rerun prepare_openmrg.py")
    k, alpha = p838(meta.Frequency_GHz.values, meta.Polarization.values)
    meta = meta.assign(k=k, alpha=alpha, kL=k * meta.Length_km.values)
    arr["radar_min"] = radar_on_minutes(arr["time"], arr["radar_time"], arr["radar_R"], cfg["radar"]["max_age_s"])
    return meta, arr


def events_of(arr, cfg, j, rsl_key="rsl", margin=None, window=None):
    return link_events(arr["tsl"][:, j], arr[rsl_key][:, j], arr["radar_min"][:, j], cfg["wet"]["threshold_mmh"],
                       cfg["wet"]["margin_min"] if margin is None else margin,
                       cfg["reference"]["window_min"] if window is None else window,
                       cfg["reference"]["min_dry_min"], cfg["event_min_minutes"])


def expected_attenuation(m, radar):
    return m.kL * np.where(radar > 0, radar, 0.0) ** m.alpha


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "experiments" / "budget" / "config.json")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "budget")
    args = ap.parse_args()
    cfg = json.loads(args.config.read_text())
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    meta, arr = load(cfg)
    n_min, n_links = arr["rsl"].shape
    tmin = arr["time"]
    month = pd.to_datetime(tmin, unit="s").month.values
    period = {"full": np.ones(n_min, bool), "cal": np.isin(month, cfg["calibration_months"])}
    evaluation_mask = np.isin(month, cfg["evaluation_months"])
    hour = (tmin // 3600) - tmin[0] // 3600
    steps = cfg["quantization_steps_db"]
    variants = {"native": "rsl", **{f"q{s:g}": f"rsl_q{s:g}" for s in steps}}
    window0, min_dry = cfg["reference"]["window_min"], cfg["reference"]["min_dry_min"]
    bins = np.asarray(cfg["rain_bins_mmh"], dtype=float)

    # ---------------- pass 1 ----------------
    link_rows, bin_rows, event_rows, lag_rows = [], [], [], []
    rain_att = {}
    for j, m in meta.iterrows():
        radar = arr["radar_min"][:, j]
        expected = expected_attenuation(m, radar)
        ev = events_of(arr, cfg, j)
        starts, stops, dry, wet = ev["starts"], ev["stops"], ev["dry"], ev["wet"]
        for name, key in variants.items():
            if name == "native":
                a = ev["rain_att"]
            else:  # the reports carry their own reference, built over the same dry minutes
                att_q = arr["tsl"][:, j].astype(float) - arr[key][:, j].astype(float)
                ref_q, _ = event_reference(att_q, dry, starts, window0, 1)
                a = att_q - held(n_min, starts, stops, ref_q)
            rain_att[(name, j)] = a.astype(np.float32)
            for pname, pmask in period.items():
                if name != "native" and pname == "full":
                    continue
                for r in bin_means(np.where(pmask, a, np.nan), expected, radar, bins, cfg["bin_min_minutes"]):
                    bin_rows.append(dict(j=j, Link=m.Link, kL=m.kL, variant=name, period=pname, **r))
        a = ev["rain_att"]
        use = np.isfinite(a) & np.isfinite(expected)
        cal = use & period["cal"]
        link_rows.append(dict(j=j, Sublink=m.Sublink, Link=m.Link, length_km=m.Length_km, freq_ghz=m.Frequency_GHz,
                              pol=m.Polarization, k=m.k, alpha=m.alpha, kL=m.kL, events=len(starts),
                              wet_min=int(use.sum()), att_mean=np.mean(a[use]) if use.any() else np.nan,
                              exp_mean=np.mean(expected[use]) if use.any() else np.nan, wet_min_cal=int(cal.sum()),
                              att_mean_cal=np.mean(a[cal]) if cal.any() else np.nan,
                              exp_mean_cal=np.mean(expected[cal]) if cal.any() else np.nan,
                              **{f"att_mean_cal_{n}": np.nanmean(np.where(cal, rain_att[(n, j)], np.nan))
                                 for n in variants if n != "native"}))
        # re-quantization errors of the reference and of the rain period, per event
        rsl = arr["rsl"][:, j].astype(float)
        var_w = arr["rsl_var"][:, j].astype(float)
        lo = np.maximum(starts - window0, 0)
        _, c_ref, s_ref = window_stats(rsl, dry, lo, starts, var_w)
        _, c_ev, s_ev = window_stats(rsl, wet, starts, stops, var_w)
        _, att_ev = window_stats(a, wet, starts, stops)
        _, exp_ev = window_stats(expected, wet & np.isfinite(a), starts, stops)
        rec = dict(j=np.full(len(starts), j), start=tmin[starts], minutes=stops - starts, n_ref=ev["n_ref"],
                   c_ref=c_ref, sig_ref=s_ref, c_ev=c_ev, sig_ev=s_ev, att_mean=att_ev, exp_mean=exp_ev)
        for s in steps:
            q = arr[f"rsl_q{s:g}"][:, j].astype(float)
            u = arr[f"offset_q{s:g}"][j]
            eps = q - rsl
            _, rec[f"eps_ref_q{s:g}"] = window_stats(eps, dry, lo, starts)
            _, rec[f"eps_ev_q{s:g}"] = window_stats(eps, wet, starts, stops)
            rec[f"pred_ref_q{s:g}"] = dither_bias(c_ref, s_ref, s, offset=u)
            rec[f"pred_ev_q{s:g}"] = dither_bias(c_ev, s_ev, s, offset=u)
        event_rows.append(pd.DataFrame(rec))
        # dry-reference test: reference held over a following dry hour, as during an event
        att = ev["att"]
        full = pd.Series(dry & np.isfinite(att)).groupby(hour).sum().values
        hour_start = np.searchsorted(hour, np.arange(hour[-1] + 1))
        h_full = hour_start[np.flatnonzero(full == 60)]
        for lag in LAGS_H:
            ref_end = h_full - (lag - 1) * 60
            ok = ref_end - window0 >= 0
            h0, ref_end = h_full[ok], ref_end[ok]
            n_r, r_mean = window_stats(att, dry, ref_end - window0, ref_end)
            _, h_mean = window_stats(att, dry, h0, h0 + 60)
            _, h_sq = window_stats(att ** 2, dry, h0, h0 + 60)
            good = n_r >= min_dry
            err = (h_mean - r_mean)[good]
            ms = (h_sq - 2 * h_mean * r_mean + r_mean ** 2)[good]
            lag_rows.append(dict(j=j, lag_h=lag, n=int(good.sum()), mean_h=np.mean(err) if good.any() else np.nan,
                                 sd_h=np.std(err) if good.sum() > 1 else np.nan,
                                 rms_min=np.sqrt(np.mean(ms)) if good.any() else np.nan))
        if j % 100 == 0:
            print(f"pass 1: {j}/{n_links} ({time.time() - t0:.0f} s)", flush=True)
    links = pd.DataFrame(link_rows)
    binned = pd.DataFrame(bin_rows)
    events = pd.concat(event_rows, ignore_index=True)
    lags = pd.DataFrame(lag_rows)

    # ---------------- length-independent offset ----------------
    bs = cfg["bootstrap"]
    usable = links[links.wet_min >= cfg["min_wet_minutes"]].copy()
    stat = lambda d: fit_offset(d.att_mean.values, d.exp_mean.values)  # noqa: E731
    w_all, beta_all = stat(usable)
    boot = cluster_bootstrap(usable, stat, bs["draws"], bs["seed"])
    rf = lambda d: fit_offset(d.att_mean.values, d.kL.values)  # noqa: E731  (radar-free: regression on kL)
    w_rf, s_rf = rf(usable)
    boot_rf = cluster_bootstrap(usable, rf, bs["draws"], bs["seed"] + 1)
    usable_cal = links[links.wet_min_cal >= cfg["min_wet_minutes"]]
    bands = {}
    for i, (name, (lo_f, hi_f)) in enumerate(BANDS.items()):
        d = usable[(usable.freq_ghz >= lo_f) & (usable.freq_ghz < hi_f)]
        bb = cluster_bootstrap(d, stat, bs["draws"], bs["seed"] + 2 + i)
        bands[name] = dict(n=len(d), W=stat(d)[0], beta=stat(d)[1], W_ci=np.percentile(bb[:, 0], [2.5, 97.5]).tolist())
    usable["offset_resid"] = usable.att_mean - (w_all + beta_all * usable.exp_mean)
    links = links.merge(usable[["j", "offset_resid"]], on="j", how="left")
    const_cal = {"native": fit_offset(usable_cal.att_mean_cal.values, usable_cal.exp_mean_cal.values)[0]}
    for n in variants:
        if n != "native":
            const_cal[n] = fit_offset(usable_cal[f"att_mean_cal_{n}"].values, usable_cal.exp_mean_cal.values)[0]

    curves = []
    for (var, pname), d in binned.groupby(["variant", "period"]):
        c = offset_curve(d).assign(variant=var, period=pname)
        c["W_radar_free"] = [fit_offset(g.a.values, g.kL.values)[0] for _, g in d.groupby("bin")]
        if var == "native" and pname == "full":
            d = d.reset_index(drop=True)
            bc = cluster_bootstrap(d, lambda x: offset_curve(x).set_index("bin").W.reindex(c.bin).values,
                                   bs["draws"], bs["seed"] + 9)
            c["W_lo"], c["W_hi"] = np.nanpercentile(bc, 2.5, axis=0), np.nanpercentile(bc, 97.5, axis=0)
        curves.append(c)
    curves = pd.concat(curves, ignore_index=True)

    sens = []
    combos = sorted({(mg, window0) for mg in cfg["sensitivity"]["margin_min"]}
                    | {(cfg["wet"]["margin_min"], w) for w in cfg["sensitivity"]["window_min"]})
    for margin, window in combos:
        rows = []
        for j, m in meta.iterrows():
            ev = events_of(arr, cfg, j, margin=margin, window=window)
            expected = expected_attenuation(m, arr["radar_min"][:, j])
            use = np.isfinite(ev["rain_att"]) & np.isfinite(expected)
            if use.sum() >= cfg["min_wet_minutes"]:
                rows.append((np.mean(ev["rain_att"][use]), np.mean(expected[use])))
        rows = np.array(rows)
        w_s, beta_s = fit_offset(rows[:, 0], rows[:, 1])
        sens.append(dict(margin_min=margin, window_min=window, links=len(rows), W=w_s, beta=beta_s))
        print(f"sensitivity margin={margin} window={window}: W={w_s:.3f} ({time.time() - t0:.0f} s)", flush=True)

    offset = dict(links=len(usable), W=w_all, beta=beta_all, W_ci=np.percentile(boot[:, 0], [2.5, 97.5]).tolist(),
                  beta_ci=np.percentile(boot[:, 1], [2.5, 97.5]).tolist(), W_radar_free=w_rf,
                  slope_radar_free=s_rf, W_radar_free_ci=np.percentile(boot_rf[:, 0], [2.5, 97.5]).tolist(),
                  bands=bands, resid_sd=float(np.std(usable.offset_resid)), links_cal=len(usable_cal),
                  W_cal=const_cal, sensitivity=sens)

    # ---------------- hourly evaluation in the evaluation month ----------------
    W_fun = {}
    for var in variants:
        c = curves[(curves.variant == var) & (curves.period == "cal")].sort_values("bin")
        W_fun[var] = offset_function(c.R.values, c.W.values)
    eval_rows = []
    hmin = cfg["hour_min_minutes"]
    for j, m in meta.iterrows():
        radar = arr["radar_min"][:, j]
        dry = events_of(arr, cfg, j)["dry"]
        for var in variants:
            a = rain_att[(var, j)].astype(float)
            inside = np.isfinite(a)
            for corr in ("none", "constant", "rain"):
                rate = np.where(dry, 0.0, np.nan)
                if corr == "rain":
                    rate[inside] = invert(a[inside], m.kL, m.alpha, W_fun[var])
                else:
                    w = 0.0 if corr == "none" else const_cal[var]
                    rate[inside] = rain_rate(a[inside], m.k, m.alpha, m.Length_km, w)
                ok = evaluation_mask & np.isfinite(rate) & np.isfinite(radar)
                if ok.sum() < cfg["min_wet_minutes"]:
                    continue
                g = pd.DataFrame({"h": hour[ok], "link": rate[ok], "radar": radar[ok]}).groupby("h")
                hm = g.mean()[g.size() >= hmin]
                x, y = hm.link.values, hm.radar.values
                wet_h = y >= 0.1
                eval_rows.append(dict(j=j, variant=var, correction=corr, hours=len(hm), acc_link=x.sum(),
                                      acc_radar=y.sum(), r=np.corrcoef(x, y)[0, 1] if x.std() > 0 else np.nan,
                                      rmse=np.sqrt(np.mean((x - y) ** 2)),
                                      pod=np.mean(x[wet_h] >= 0.1) if wet_h.any() else np.nan,
                                      far=np.mean(y[x >= 0.1] < 0.1) if (x >= 0.1).any() else np.nan))
        if j % 200 == 0:
            print(f"evaluation: {j}/{n_links} ({time.time() - t0:.0f} s)", flush=True)
    evaluation = pd.DataFrame(eval_rows)

    links.to_csv(args.out / "links.csv", index=False)
    binned.to_csv(args.out / "bins.csv", index=False)
    curves.to_csv(args.out / "offset_curve.csv", index=False)
    events.to_csv(args.out / "events.csv.gz", index=False)
    lags.to_csv(args.out / "dry_lag.csv", index=False)
    evaluation.to_csv(args.out / "evaluation.csv", index=False)
    (args.out / "offset.json").write_text(json.dumps(offset, indent=1, default=float))
    (args.out / "run_info.json").write_text(json.dumps({"runtime_s": round(time.time() - t0), "links": n_links,
                                                         "usable_links": len(usable), "events": len(events)}, indent=1))
    print(f"done in {time.time() - t0:.0f} s: {len(usable)} links, {len(events)} events, W={w_all:.3f} dB")


if __name__ == "__main__":
    main()
