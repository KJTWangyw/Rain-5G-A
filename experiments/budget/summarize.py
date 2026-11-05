#!/usr/bin/env python3
"""Every number quoted in the manuscript, from results/budget/ -> results/budget/summary.json.

Values are stored formatted for the text ({{key}} placeholders in paper/content_zh.md); tables are stored
as lists of rows. Usage:  python experiments/budget/summarize.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from cmlbudget.offset import offset_function  # noqa: E402
from cmlbudget.quantization import dither_bias_rms  # noqa: E402
from nrpilot.uncertainty import K_DB  # noqa: E402

RES = ROOT / "results" / "budget"
CFG = json.loads((ROOT / "experiments" / "budget" / "config.json").read_text())


def f(x, d=2):
    return f"{x:.{d}f}"


def pct(x, d=0):
    return f"{100 * x:.{d}f}"


def ci(lo_hi, d=2):
    return f"{lo_hi[0]:.{d}f}–{lo_hi[1]:.{d}f}"


def length_class(L):
    edges = CFG["length_classes_km"]
    return pd.cut(L, edges, labels=[f"{a:g}–{b:g}" for a, b in zip(edges[:-1], edges[1:])])


def req(db, kL, alpha):
    """Equivalent rain rate of an attenuation of ``db`` on a path with kL, alpha."""
    return (np.asarray(db, dtype=float) / kL) ** (1.0 / alpha)


def cluster_ci(values, clusters, stat, draws=2000, seed=20261108):
    """95% interval of ``stat`` over resamples of whole clusters (physical links)."""
    rng = np.random.default_rng(seed)
    values, clusters = np.asarray(values), np.asarray(clusters)
    groups = [np.flatnonzero(clusters == c) for c in np.unique(clusters)]
    out = [stat(values[np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])])
           for _ in range(draws)]
    return np.percentile(out, [2.5, 97.5])


def sign_test(diff):
    """Two-sided exact sign test p-value (zeros dropped)."""
    from math import comb
    d = np.asarray(diff)
    d = d[d != 0]
    n, k = len(d), int(np.sum(d > 0))
    tail = sum(comb(n, i) for i in range(0, min(k, n - k) + 1))
    return min(1.0, 2.0 * tail / 2 ** n)


def extension_summary(use):
    """Numbers of the drying, correction-model, link-classification and gauge analyses."""
    S = {}
    ext = CFG["extension"]
    tags = ["L05", "L1", "L2", "L4", "L20"]
    # drying after the radar rain on the path has ended
    d = pd.read_csv(RES / "drying.csv").merge(use[["j", "length_km"]], on="j")
    d = d[d.n >= 5]
    med = d.groupby("lag_min")["mean"].median()
    S["dry_first"] = f(med.iloc[0])
    S["dry_30"] = f(float(med[med.index >= 30].iloc[0]))
    S["dry_60"] = f(float(med[med.index >= 60].iloc[0]), 3)
    w = ext["drying_bin_min"] / 2.0
    half = float(med.index[np.argmax(med.values < med.iloc[0] / 2)])
    b01 = float(med.index[np.argmax(med.values < 0.1)])
    S["dry_half_lo"], S["dry_half_hi"] = f"{half - w:g}", f"{half + w:g}"
    S["dry_below01_lo"], S["dry_below01_hi"] = f"{b01 - w:g}", f"{b01 + w:g}"
    d["c3"] = pd.cut(d.length_km, [0, 1, 4, 20], labels=["short", "mid", "long"])
    m3 = d.groupby(["lag_min", "c3"], observed=True)["mean"].median().unstack()
    for c in ("short", "mid", "long"):
        S[f"dry_first_{c}"] = f(m3[c].iloc[0])
    late = m3[m3.index >= 10]
    S["dry_spread_after10"] = f(float((late.max(axis=1) - late.min(axis=1)).max()), 2)
    S["drying_curve"] = m3.round(4).reset_index().to_dict("records")
    S["n_dry_links"] = int(d.j.nunique())
    # KR-alt literature curve
    p = ext["kralt"]
    kr = lambda R: p["A_max"] * (1 - np.exp(-p["d"] * R ** p["zeta"]))  # noqa: E731
    S["kr_Amax"], S["kr_zeta"], S["kr_d"] = f"{p['A_max']:g}", f"{p['zeta']:g}", f"{p['d']:g}"
    S["kr_at_1"], S["kr_at_5"] = f(kr(1.0)), f(kr(5.0))
    full = pd.read_csv(RES / "offset_curve.csv")
    full = full[(full.variant == "native") & (full.period == "full")].sort_values("bin")
    mid = full[(full.R > 0.1) & (full.R < 2)]
    diff = kr(mid.R.values) - mid.W.values
    S["kr_minus_bin_lo"], S["kr_minus_bin_hi"] = f(diff.min()), f(diff.max())
    S["kr_over_bin_top"] = f(kr(full.R.iloc[-1]) / full.W.iloc[-1], 1)
    # evaluation against radar (all usable sub-links) and gauges (sub-links within gauge_km of a city gauge)
    e = pd.read_csv(RES / "evaluation_ext.csv").merge(use[["j", "length_km", "cls"]], on="j")
    e["ratio"] = e.acc_est / e.acc_ref
    rad = e[e.ref == "radar"]
    cls_order = list(rad.cls.cat.categories) if hasattr(rad.cls, "cat") else sorted(rad.cls.unique())
    tab = rad.groupby(["cls", "mask", "corr"], observed=True).agg(rmse=("rmse", "median"), ratio=("ratio", "median"),
                                                                 r=("r", "median"), n=("j", "size")).reset_index()
    S["eval_radar"] = tab.round(4).assign(cls=tab.cls.astype(str)).to_dict("records")
    for (cls, mask, corr), row in tab.set_index(["cls", "mask", "corr"]).iterrows():
        t = tags[cls_order.index(cls)]
        S[f"er_{mask}_{corr}_{t}"] = f(row.rmse)
        S[f"eq_{mask}_{corr}_{t}"] = f(row.ratio)
    piv = rad[rad["mask"] == "radar"].pivot_table(index=["j", "Link", "cls"], columns="corr", values="rmse",
                                                   observed=True).reset_index()
    for i, cls in enumerate(cls_order):
        g = piv[piv.cls == cls].dropna(subset=["none", "kralt"])
        red = 1 - g.kralt / g.none
        lo, hi = cluster_ci(red, g.Link, np.median)
        S[f"red_kr_{tags[i]}"], S[f"red_kr_{tags[i]}_ci"] = pct(np.median(red)), f"{100 * lo:.0f}–{100 * hi:.0f}"
        S[f"impr_kr_{tags[i]}"] = pct(np.mean(g.kralt < g.none))
        red_c = 1 - g.constant / g.none
        S[f"red_const_{tags[i]}"] = pct(np.median(red_c))
    q = rad[(rad["corr"] == "kralt") & rad["mask"].isin(["radar", "radar_q1", "radar_q2"])]
    qp = q.pivot_table(index=["j", "cls"], columns="mask", values="rmse", observed=True).dropna().reset_index()
    for step in ("q1", "q2"):
        inc = qp.groupby("cls", observed=True).apply(lambda g: np.median(g[f"radar_{step}"] / g.radar - 1),
                                                    include_groups=False)
        S[f"kr_{step}_inc_05"], S[f"kr_{step}_inc_4"] = pct(inc.iloc[0], 1), pct(inc.iloc[-1], 1)
    lm = pd.read_csv(RES / "link_mask.csv")
    S["lm_hit"], S["lm_false"] = pct(lm.hit.median()), pct(lm.false.median(), 1)
    S["lm_window"] = ext["std_window_min"]
    gau = e[e.ref == "gauge"].copy()
    gau["c2"] = np.where(gau.length_km <= 1, "le1", "gt1")
    gt = gau.groupby(["c2", "mask", "corr"]).agg(rmse=("rmse", "median"), ratio=("ratio", "median"),
                                                r=("r", "median"), n=("j", "size"), phys=("Link", "nunique"),
                                                km=("gauge_km", "median")).reset_index()
    S["eval_gauge"] = gt.round(4).to_dict("records")
    for _, row in gt.iterrows():
        S[f"eg_{row['mask'].replace('-', 'none')}_{row['corr']}_{row.c2}"] = f(row.rmse)
        S[f"egq_{row['mask'].replace('-', 'none')}_{row['corr']}_{row.c2}"] = f(row.ratio)
        S[f"egr_{row['mask'].replace('-', 'none')}_{row['corr']}_{row.c2}"] = f(row.r)
    for c2 in ("le1", "gt1"):
        g = gau[gau.c2 == c2]
        S[f"g_n_{c2}"] = int(g[(g["mask"] == "radar") & (g["corr"] == "kralt")].j.nunique())
        S[f"g_phys_{c2}"] = int(g.Link.nunique())
        # paired comparison KR-alt-corrected link vs radar, averaged per physical link
        a = g[(g["mask"] == "radar") & (g["corr"] == "kralt")].groupby("Link").rmse.mean()
        b = g[(g["mask"] == "-") & (g["corr"] == "radar")].groupby("Link").rmse.mean()
        both = pd.concat([a, b], axis=1, keys=["kr", "radar"]).dropna()
        S[f"g_kr_better_{c2}"] = int((both.kr < both.radar).sum())
        S[f"g_pairs_{c2}"] = len(both)
        pv = sign_test(both.kr - both.radar) if len(both) else float("nan")
        mant, expo = f"{pv:.1e}".split("e")
        S[f"g_sign_p_{c2}"] = f"{mant}\\times10^{{{int(expo)}}}"
    S["gauge_km"] = f"{ext['gauge_km']:g}"
    S["n_gauges"] = int(gau.gauge.nunique())
    return S


def main():
    links = pd.read_csv(RES / "links.csv")
    offset = json.loads((RES / "offset.json").read_text())
    curves = pd.read_csv(RES / "offset_curve.csv")
    events = pd.read_csv(RES / "events.csv.gz")
    lags = pd.read_csv(RES / "dry_lag.csv")
    evaluation = pd.read_csv(RES / "evaluation.csv")
    pilot = pd.read_csv(RES / "pilot_reference.csv")
    run = json.loads((RES / "run_info.json").read_text())
    S = {}

    use = links[links.wet_min >= CFG["min_wet_minutes"]].copy()
    use["cls"] = length_class(use.length_km)
    S.update(n_sel=len(links), n_use=len(use), n_phys=use.Link.nunique(), n_events=f"{len(events):,}".replace(",", " "),
             L_min=f(use.length_km.min()), L_max=f(use.length_km.max(), 1), L_med=f(use.length_km.median()),
             f_min=f(use.freq_ghz.min(), 1), f_max=f(use.freq_ghz.max(), 1),
             n_short=int((use.length_km <= 1).sum()), ev_med_min=int(events.minutes.median()),
             ref_med_min=int(events.n_ref.median()), n_vertical=int((use.pol == "Vertical").sum()))
    S["class_counts"] = use.cls.value_counts().sort_index().astype(int).to_dict()

    # ---- constant offset
    S.update(W=f(offset["W"]), W_ci=ci(offset["W_ci"]), beta=f(offset["beta"]), beta_ci=ci(offset["beta_ci"]),
             W_rf=f(offset["W_radar_free"]), W_rf_ci=ci(offset["W_radar_free_ci"]), W_cal=f(offset["W_cal"]["native"]),
             W_resid_sd=f(offset["resid_sd"]))
    for name, b in offset["bands"].items():
        key = name.replace("-", "_")
        S[f"W_{key}"], S[f"W_{key}_ci"], S[f"n_{key}"] = f(b["W"]), ci(b["W_ci"]), b["n"]
    sens = [s["W"] for s in offset["sensitivity"]]
    S["W_sens"] = f"{min(sens):.2f}–{max(sens):.2f}"
    use["share"] = offset["W"] / use.att_mean
    share = use.groupby("cls", observed=True).share.median()
    S["share_05"], S["share_1"], S["share_4"] = pct(share.iloc[0]), pct(share.iloc[1]), pct(share.iloc[-1])

    # ---- rain-dependent offset
    full = curves[(curves.variant == "native") & (curves.period == "full")].sort_values("bin")
    cal = curves[(curves.variant == "native") & (curves.period == "cal")].sort_values("bin")
    S["curve"] = full[["bin", "R", "W", "W_lo", "W_hi", "beta", "resid_sd", "W_radar_free", "links"]].round(4).to_dict("records")
    S["curve_cal"] = cal[["bin", "R", "W", "beta"]].round(4).to_dict("records")
    S["W_low"], S["W_high"] = f(full.W.iloc[0]), f(full.W.iloc[-1])
    S["W_low_ci"], S["W_high_ci"] = ci([full.W_lo.iloc[0], full.W_hi.iloc[0]]), ci([full.W_lo.iloc[-1], full.W_hi.iloc[-1]])
    S["R_low"], S["R_high"] = f(full.R.iloc[0], 2), f(full.R.iloc[-1], 1)
    S["W_rf_maxdiff"] = f(np.max(np.abs(full.W_radar_free - full.W)))
    bins = CFG["rain_bins_mmh"]
    S["bin_low"], S["bin_high"] = f"{bins[1]:g}", f"{bins[-2]:g}"
    Wf = offset_function(full.R.values, full.W.values)
    for R in (1.0, 5.0):
        S[f"W_at_{R:g}"] = f(float(Wf(R)))
        ratio = Wf(R) / (use.kL * R ** use.alpha)
        S[f"W_over_rain_{R:g}_05"] = f(np.median(ratio[use.length_km <= 0.5]), 1)
        S[f"W_over_rain_{R:g}_4"] = f(np.median(ratio[use.length_km > 4]), 1)
    # rain rate at which the rain attenuation of a link equals W(R) (crossover)
    grid = np.logspace(-2, 2.5, 2000)
    cross = []
    for m in use.itertuples():
        diff = m.kL * grid ** m.alpha - Wf(grid)
        cross.append(grid[np.argmax(diff > 0)] if (diff > 0).any() else np.inf)
    use["R_cross"] = cross
    S["frac_never_cross_1km"] = pct(np.mean(use.R_cross[use.length_km <= 1] > 50))
    S["R_cross_med"] = {str(k): (f(v, 1) if np.isfinite(v) else ">50") for k, v in use.groupby("cls", observed=True).R_cross.median().items()}

    # ---- reference drift (dry hours)
    lag = lags.groupby("lag_h").agg(sd=("sd_h", "median"), q1=("sd_h", lambda s: s.quantile(.25)),
                                    q3=("sd_h", lambda s: s.quantile(.75)), rms_min=("rms_min", "median"),
                                    bias=("mean_h", "median"), n=("n", "sum"))
    S["drift"] = lag.round(4).reset_index().to_dict("records")
    S["drift_1h"], S["drift_6h"], S["drift_min"] = f(lag.sd.loc[1]), f(lag.sd.loc[6]), f(lag.rms_min.loc[1])
    S["drift_bias"] = f(abs(lag.bias.loc[1]), 3)
    S["n_dry_hours"] = f"{int(lag.n.loc[1]):,}".replace(",", " ")

    # ---- quantization
    ev = events.dropna(subset=["eps_ref_q1", "pred_ref_q1", "eps_ev_q1"])
    S["n_events_q"] = f"{len(ev):,}".replace(",", " ")
    rms = lambda x: float(np.sqrt(np.mean(np.square(x))))  # noqa: E731
    for s in CFG["quantization_steps_db"]:
        q = f"q{s:g}"
        e_r, p_r, e_e, p_e = ev[f"eps_ref_{q}"], ev[f"pred_ref_{q}"], ev[f"eps_ev_{q}"], ev[f"pred_ev_{q}"]
        slope, icpt = np.polyfit(p_r, e_r, 1)
        S[f"{q}_ref_rms"], S[f"{q}_ref_pred"] = f(rms(e_r)), f(rms(p_r))
        S[f"{q}_ev_rms"], S[f"{q}_ev_pred"] = f(rms(e_e)), f(rms(p_e))
        S[f"{q}_tot_rms"] = f(rms(e_r - e_e))
        S[f"{q}_slope"], S[f"{q}_r2"] = f(slope), f(np.corrcoef(p_r, e_r)[0, 1] ** 2)
        S[f"{q}_uniform"] = f(s / np.sqrt(12))
    S["sig_ref_med"], S["sig_ev_med"] = f(ev.sig_ref.median()), f(ev.sig_ev.median())
    S["rms_b_at_med"] = f(float(dither_bias_rms(ev.sig_ref.median(), 1.0)))

    # ---- pilot reference
    pilot["sd_sim"] = np.sqrt(pilot.S_sim.clip(lower=0))
    pilot["rel"] = pilot.S_sim / pilot.S_complete - 1
    S["pilot_rel_range"] = f"{100 * pilot.rel.min():+.0f}%至{100 * pilot.rel.max():+.0f}%".replace("-", "−")
    S["pilot_in_ci"] = int(((pilot.S_lo <= pilot.S_complete) & (pilot.S_complete <= pilot.S_hi)).sum())
    S["pilot_cells"] = len(pilot)
    J = int(pilot.J.iloc[0])
    S["J"] = J
    sd = lambda snr, ms: K_DB / np.sqrt(2 * J * (ms / 0.125) * 10 ** (snr / 10))  # noqa: E731  (shared part of the bound)
    S["pilot_sd_15_1ms"], S["pilot_sd_15_10ms"], S["pilot_sd_5_10ms"] = f(sd(15, 1), 3), f(sd(15, 10), 4), f(sd(5, 10), 3)
    S["pilot_ms_for_001_5db"] = f(0.125 * K_DB ** 2 / (2 * J * 10 ** 0.5 * 0.01 ** 2), 1)
    S["pilot_table"] = pilot[["snr_db", "M", "reference_ms", "S_sim", "S_lo", "S_hi", "S_complete", "S_crb"]].to_dict("records")

    # ---- evaluation (August)
    e = evaluation.merge(use[["j", "length_km", "cls", "kL"]], on="j")
    e["ratio"] = e.acc_link / e.acc_radar
    tab = e.groupby(["cls", "variant", "correction"], observed=True).agg(ratio=("ratio", "median"), rmse=("rmse", "median"),
                                                                        r=("r", "median"), pod=("pod", "median"),
                                                                        n=("j", "size")).reset_index()
    S["evaluation"] = tab.round(4).to_dict("records")

    def ev_get(cls_idx, var, corr, col):
        c = tab.cls.cat.categories[cls_idx] if hasattr(tab.cls, "cat") else sorted(tab.cls.unique())[cls_idx]
        return float(tab[(tab.cls == c) & (tab.variant == var) & (tab.correction == corr)][col].iloc[0])
    for i, tag in ((0, "05"), (1, "1"), (2, "2"), (4, "4")):
        for corr in ("none", "constant", "rain"):
            S[f"ev_ratio_{corr}_{tag}"] = f(ev_get(i, "native", corr, "ratio"), 1)
            S[f"ev_rmse_{corr}_{tag}"] = f(ev_get(i, "native", corr, "rmse"))
        S[f"ev_red_{tag}"] = pct(1 - ev_get(i, "native", "rain", "rmse") / ev_get(i, "native", "none", "rmse"))
        S[f"ev_rmse_q1_{tag}"] = f(ev_get(i, "q1", "rain", "rmse"))
        S[f"ev_rmse_q2_{tag}"] = f(ev_get(i, "q2", "rain", "rmse"))
        S[f"ev_q1_inc_{tag}"] = pct(ev_get(i, "q1", "rain", "rmse") / ev_get(i, "native", "rain", "rmse") - 1)
        S[f"ev_q2_inc_{tag}"] = pct(ev_get(i, "q2", "rain", "rmse") / ev_get(i, "native", "rain", "rmse") - 1)
        S[f"ev_r_none_{tag}"], S[f"ev_r_rain_{tag}"] = f(ev_get(i, "native", "none", "r")), f(ev_get(i, "native", "rain", "r"))

    # ---- error budget at R = 1 mm/h, medians over the sub-links of each class
    q1 = events.groupby("j").apply(lambda d: rms((d.eps_ref_q1 - d.eps_ev_q1).dropna()), include_groups=False)
    lag1 = lags[lags.lag_h == 1].set_index("j").sd_h
    b = use.set_index("j")
    b["signal"] = b.kL
    b["W1"] = float(Wf(1.0))
    resid_1 = float(full[full.bin == 3].resid_sd.iloc[0])  # residual spread of the offset in the 1-2 mm/h bin
    b["resid"] = resid_1
    b["drift"] = lag1.reindex(b.index)
    b["quant"] = q1.reindex(b.index)
    b["pilot"] = sd(15, 10)
    rows = []
    for cls, g in b.groupby("cls", observed=True):
        row = {"cls": str(cls), "n": len(g)}
        for term in ("signal", "W1", "resid", "drift", "quant", "pilot"):
            row[f"{term}_db"] = float(g[term].median())
            row[f"{term}_mmh"] = float(np.median(req(g[term], g.kL, g.alpha)))
        rows.append(row)
    S["budget"] = rows
    S["resid_1"] = f(resid_1)
    tags = dict(zip([r["cls"] for r in rows], ["L05", "L1", "L2", "L4", "L20"]))
    for r in rows:
        t = tags[r["cls"]]
        S[f"b_{t}_n"] = r["n"]
        S[f"b_{t}_signal_db"] = f(r["signal_db"])
        for term in ("W1", "resid", "drift", "quant", "pilot"):
            v = r[f"{term}_mmh"]
            S[f"b_{t}_{term}_mmh"] = f(v, 3 if v < 0.01 else 2 if v < 0.1 else 1)
            S[f"b_{t}_{term}_db"] = f(r[f"{term}_db"], 4 if term == "pilot" else 2)
    rc = use.groupby("cls", observed=True).R_cross.median()
    for (cls, v), t in zip(rc.items(), ["L05", "L1", "L2", "L4", "L20"]):
        S[f"Rc_{t[1:]}"] = f(v, 1)
    S["runtime_s"] = run["runtime_s"]

    # ---- robustness of the intercepts (Theil-Sen)
    def theil_sen(x, y):
        i, k_ = np.triu_indices(len(x), 1)
        dx = x[k_] - x[i]
        ok = np.abs(dx) > 1e-9
        slope = np.median((y[k_] - y[i])[ok] / dx[ok])
        return float(np.median(y - slope * x)), float(slope)
    S["W_ts"] = f(theil_sen(use.exp_mean.values, use.att_mean.values)[0])
    bins_tab = pd.read_csv(RES / "bins.csv")
    bins_tab = bins_tab[(bins_tab.variant == "native") & (bins_tab.period == "full")]
    ts = np.array([theil_sen(g.e.values, g.a.values)[0] for _, g in bins_tab.groupby("bin")])
    S["W_ts_maxdiff"] = f(np.max(np.abs(ts - full.W.values)))
    S["W_ts_maxdiff_low"] = f(np.max(np.abs(ts[:-1] - full.W.values[:-1])))
    S["W_alt_maxdiff_low"] = f(max(np.max(np.abs(ts[:-1] - full.W.values[:-1])),
                                   np.max(np.abs(full.W_radar_free.values[:-1] - full.W.values[:-1]))))
    S["W_cal_top"] = f(cal.W.iloc[-1])

    # ---- effect of the report resolution after the W(R) correction, over the length classes
    inc1 = [ev_get(i, "q1", "rain", "rmse") / ev_get(i, "native", "rain", "rmse") - 1 for i in range(5)]
    inc2 = [ev_get(i, "q2", "rain", "rmse") / ev_get(i, "native", "rain", "rmse") - 1 for i in range(5)]
    S["q1_inc_range"] = f"{100 * min(inc1):+.0f}%至{100 * max(inc1):+.0f}%".replace("-", "−")
    S["q2_inc_range"] = f"{100 * min(inc2):+.0f}%至{100 * max(inc2):+.0f}%".replace("-", "−")

    # ---- allowed shared attenuation error for a rain-rate error of 1 mm/h at 1 mm/h, 28 GHz vertical
    from nrpilot.rain import p838_vertical
    k28, a28 = p838_vertical(28.0)
    S["k28"], S["a28"] = f(k28, 3), f(a28, 3)
    S["req_05_01"] = f((0.1 / (k28 * 0.5)) ** (1 / a28), 1)
    for L in (0.3, 1.0):
        S[f"sigma_req_{str(L).replace('.', '')}"] = f(a28 * k28 * L, 2 if L >= 1 else 3)

    # ---- the example event of Fig. 1 (written by make_figures.py)
    info = RES / "figures_info.json"
    if info.exists():
        ex = json.loads(info.read_text())["example"]
        S.update(ex_sublink=ex["sublink"], ex_L=f(ex["length_km"]), ex_f=f(ex["freq_ghz"], 1), ex_hours=f(ex["hours"], 1),
                 ex_date=ex["start"][:10], ex_radar_mm=f(ex["radar_mm"], 0), ex_link_none_mm=f(ex["link_none_mm"], 0),
                 ex_link_const_mm=f(ex["link_const_mm"], 0), ex_link_rain_mm=f(ex["link_rain_mm"], 0),
                 ex_link_kralt_mm=f(ex.get("link_kralt_mm", float("nan")), 0), ex_red_med=pct(ex["rule_median_reduction"]))

    # ---- extension: drying signature, correction models, link-only classification, gauge reference
    S.update(extension_summary(use))
    (RES / "summary.json").write_text(json.dumps(S, indent=1, ensure_ascii=False, default=float))
    print(json.dumps({k: v for k, v in S.items() if not isinstance(v, (list, dict))}, indent=0, ensure_ascii=False))


if __name__ == "__main__":
    main()
