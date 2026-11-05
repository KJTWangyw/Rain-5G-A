#!/usr/bin/env python3
"""Figures of the Chinese manuscript (600 dpi PNG for the Word build, PDF for reference).

Usage:  python experiments/budget/make_figures.py [--out paper/figures_zh]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments" / "budget"))
from cmlbudget.offset import invert, offset_function  # noqa: E402
from cmlbudget.retrieval import rain_rate  # noqa: E402
from nrpilot.rain import p838_vertical  # noqa: E402
from nrpilot.uncertainty import complete_shared_independent  # noqa: E402

RES = ROOT / "results" / "budget"
CFG = json.loads((ROOT / "experiments" / "budget" / "config.json").read_text())
plt.rcParams.update({
    "font.family": ["Noto Serif CJK SC", "Times New Roman"], "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman", "mathtext.it": "Times New Roman:italic", "mathtext.bf": "Times New Roman:bold",
    "mathtext.sf": "Times New Roman", "mathtext.tt": "Times New Roman", "mathtext.cal": "Times New Roman:italic",
    "mathtext.fallback": "stix", "font.size": 8.5, "axes.linewidth": 0.6, "axes.titlesize": 9,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.frameon": False, "savefig.dpi": 600,
    "axes.unicode_minus": True,
})
C = {"measured": "#222222", "radar": "#0f6e6e", "offset": "#c26a1b", "const": "#8c8c8c", "rain": "#4b3f8f",
     "q1": "#b2182b", "q2": "#e08a8a", "pilot": "#0f6e6e"}
CLASS_COLORS = ["#b2182b", "#e08a3c", "#7f7f7f", "#4b8fb8", "#0f3e6e"]
CLASS_LABELS = ["≤0.5 km", "0.5–1 km", "1–2 km", "2–4 km", ">4 km"]


def save(fig, out, name):
    fig.savefig(out / f"{name}.png", bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def classes(L):
    return np.digitize(L, CFG["length_classes_km"][1:-1], right=True)


def offset_curves():
    cur = pd.read_csv(RES / "offset_curve.csv")
    full = cur[(cur.variant == "native") & (cur.period == "full")].sort_values("bin")
    cal = cur[(cur.variant == "native") & (cur.period == "cal")].sort_values("bin")
    return full, cal


def fig_example(out):
    import analyze_openmrg as A
    from evaluate_corrections import kralt
    meta, arr = A.load(CFG)
    links = pd.read_csv(RES / "links.csv")
    _, cal = offset_curves()
    Wf = offset_function(cal.R.values, cal.W.values)
    Wk = kralt(CFG["extension"]["kralt"])
    W_const = json.loads((RES / "offset.json").read_text())["W_cal"]["native"]
    month = pd.to_datetime(arr["time"], unit="s").month.values
    # rule: the <=0.5 km sub-link whose August RMSE reduction (no correction -> KR-alt) is closest to the class median,
    # then its August event (<= 8 h) with the largest radar accumulation
    evl = pd.read_csv(RES / "evaluation_ext.csv")
    evl = evl[(evl.ref == "radar") & (evl["mask"] == "radar")].pivot_table(index="j", columns="corr", values="rmse")
    short = links[(links.length_km <= 0.5) & (links.wet_min >= CFG["min_wet_minutes"])].set_index("j")
    red = (1 - evl["kralt"] / evl["none"]).reindex(short.index).dropna()
    j = int((red - red.median()).abs().idxmin())
    ev = A.events_of(arr, CFG, j)
    radar = arr["radar_min"][:, j]
    best = None
    for s0, e0 in zip(ev["starts"], ev["stops"]):
        if month[s0] not in CFG["evaluation_months"] or e0 - s0 > 480:
            continue
        acc = np.nansum(radar[s0:e0]) / 60.0
        if best is None or acc > best[0]:
            best = (acc, s0, e0)
    _, s, e = best
    m = meta.loc[j]
    lo, hi = max(s - 60, 0), min(e + 60, len(radar))
    t = (arr["time"][lo:hi] - arr["time"][s]) / 3600.0
    a = ev["rain_att"][lo:hi]
    R = radar[lo:hi]
    exp_att = m.kL * np.where(R > 0, R, 0) ** m.alpha
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.3, 4.4), sharex=True, gridspec_kw={"hspace": 0.28})
    ax1.plot(t, a, color=C["measured"], lw=0.8, label="实测降雨期衰减")
    ax1.plot(t, exp_att, color=C["radar"], lw=1.0, label=r"雷达对应的雨致衰减 $kLR^{\alpha}$")
    ax1.plot(t, Wk(R) + exp_att, color=C["rain"], lw=1.0, ls="--", label=r"KR-alt湿天线衰减 + $kLR^{\alpha}$")
    ax1.set_ylabel("衰减 (dB)")
    ax1.legend(loc="upper left", fontsize=7, ncol=2)
    ax1.set_title(f"(a) 子链路 {int(m.Sublink)}：{m.Length_km:.2f} km，{m.Frequency_GHz:.1f} GHz", loc="left")
    inside = np.isfinite(a)
    dry = ev["dry"][lo:hi]
    rates = {}
    for name, fn in (("none", lambda x: rain_rate(x, m.k, m.alpha, m.Length_km, 0.0)),
                     ("const", lambda x: rain_rate(x, m.k, m.alpha, m.Length_km, W_const)),
                     ("rain", lambda x: invert(x, m.kL, m.alpha, Wf)),
                     ("kralt", lambda x: invert(x, m.kL, m.alpha, Wk))):
        r = np.where(dry, 0.0, np.nan)
        r[inside] = fn(a[inside])
        rates[name] = r
    ax2.plot(t, R, color=C["radar"], lw=1.2, label="雷达路径平均")
    ax2.plot(t, rates["none"], color=C["const"], lw=0.7, label="不扣除偏移")
    ax2.plot(t, rates["rain"], color=C["offset"], lw=0.7, label=r"扣除分档 $W(R)$")
    ax2.plot(t, rates["kralt"], color=C["rain"], lw=0.9, label="扣除KR-alt湿天线衰减")
    ax2.set_ylabel("雨强 (mm/h)")
    ax2.set_xlabel("相对降雨事件开始的时间 (h)")
    ax2.legend(loc="upper left", fontsize=7, ncol=2)
    ax2.set_title("(b) 雨强反演", loc="left")
    ax1.set_ylim(top=1.45 * np.nanmax(a))
    ax2.set_ylim(top=1.4 * max(np.nanmax(rates["none"]), np.nanmax(R)))
    for ax in (ax1, ax2):
        ax.axvspan(0, (e - s) / 60.0, color="#f0f0f0", zorder=0)
        ax.set_xlim(t[0], t[-1])
    save(fig, out, "fig1_example")
    acc = {k: float(np.nansum(v) / 60.0) for k, v in rates.items()}
    return dict(rule_median_reduction=float(red.median()), picked_reduction=float(red.loc[j]), sublink=int(m.Sublink),
                length_km=float(m.Length_km), freq_ghz=float(m.Frequency_GHz),
                start=str(pd.to_datetime(arr["time"][s], unit="s")), hours=(e - s) / 60.0,
                radar_mm=float(np.nansum(R) / 60.0), link_none_mm=acc["none"], link_const_mm=acc["const"],
                link_rain_mm=acc["rain"], link_kralt_mm=acc["kralt"])


def fig_offset(out):
    from evaluate_corrections import kralt
    links = pd.read_csv(RES / "links.csv")
    links = links[links.wet_min >= CFG["min_wet_minutes"]]
    off = json.loads((RES / "offset.json").read_text())
    full, cal = offset_curves()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(6.6, 2.45), gridspec_kw={"wspace": 0.42})
    cls = classes(links.length_km.values)
    for c in range(5):
        s = cls == c
        ax1.scatter(links.exp_mean[s], links.att_mean[s], s=4, color=CLASS_COLORS[c], label=CLASS_LABELS[c], lw=0)
    x = np.linspace(0, links.exp_mean.max() * 1.05, 10)
    ax1.plot(x, off["W"] + off["beta"] * x, color="k", lw=0.8)
    ax1.annotate(rf"$W$ = {off['W']:.2f} dB", xy=(0.0, off["W"]), xytext=(0.6, 0.22), fontsize=7,
                 arrowprops=dict(arrowstyle="->", lw=0.6, color="k"))
    ax1.set_xlabel("雷达对应雨致衰减 (dB)")
    ax1.set_ylabel("链路湿期平均衰减 (dB)")
    ax1.set_xlim(0, None), ax1.set_ylim(0, None)
    ax1.legend(fontsize=5.8, loc="upper left", handletextpad=0.1, markerscale=2, borderaxespad=0.2)
    ax1.set_title("(a) 跨链路回归", loc="left")
    ax2.fill_between(full.R, full.W_lo, full.W_hi, color=C["offset"], alpha=0.25, lw=0)
    ax2.plot(full.R, full.W, "o-", color=C["offset"], ms=2.5, lw=1, label=r"分档 $W(R)$")
    R = np.logspace(-2, np.log10(30), 200)
    ax2.plot(R, kralt(CFG["extension"]["kralt"])(R), color=C["rain"], lw=1, ls="--", label="KR-alt（文献参数）")
    k28, a28 = p838_vertical(28.0)
    for L, ls in ((0.3, ":"), (1.0, "-.")):
        ax2.plot(R, k28 * L * R ** a28, color=C["radar"], lw=0.7, ls=ls, label=f"雨致衰减，{L:g} km")
    ax2.set_xscale("log")
    ax2.set_ylim(0, 3.6), ax2.set_xlim(5e-3, 35)
    ax2.set_xlabel("雨强 $R$ (mm/h)")
    ax2.set_ylabel("长度无关偏移 (dB)")
    ax2.legend(fontsize=5.8, loc="upper left", borderaxespad=0.2)
    ax2.set_title("(b) 随雨强的变化", loc="left")
    d = pd.read_csv(RES / "drying.csv").merge(links[["j", "length_km"]], on="j")
    d = d[d.n >= 5]
    d["c3"] = pd.cut(d.length_km, [0, 1, 4, 20], labels=["≤1 km", "1–4 km", ">4 km"])
    for (lab, g), col in zip(d.groupby("c3", observed=True), ("#b2182b", "#7f7f7f", "#0f3e6e")):
        med = g.groupby("lag_min")["mean"].median()
        ax3.plot(med.index, med.values, "o-", ms=2, lw=0.9, color=col, label=lab)
    ax3.axhline(0, color="k", lw=0.4)
    ax3.set_xlabel("雷达雨停后的时间 (min)")
    ax3.set_ylabel("残余衰减 (dB)")
    ax3.set_xlim(0, 90)
    ax3.legend(fontsize=6, loc="upper right", borderaxespad=0.2)
    ax3.set_title("(c) 雨停后的消退", loc="left")
    save(fig, out, "fig2_offset")


def fig_drift_quant(out):
    lags = pd.read_csv(RES / "dry_lag.csv")
    ev = pd.read_csv(RES / "events.csv.gz").dropna(subset=["eps_ref_q1", "pred_ref_q1"])
    S = json.loads((RES / "summary.json").read_text()) if (RES / "summary.json").exists() else {}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.3, 2.75), gridspec_kw={"wspace": 0.3})
    g = lags.groupby("lag_h").sd_h
    x = g.median().index.values
    ax1.fill_between(x, g.quantile(0.25), g.quantile(0.75), color=C["radar"], alpha=0.2, lw=0)
    ax1.plot(x, g.median(), "o-", color=C["radar"], ms=3, lw=1, label="小时平均误差标准差（中位数，四分位）")
    ax1.plot(x, lags.groupby("lag_h").rms_min.median(), "s--", color=C["measured"], ms=2.5, lw=0.8, label="逐分钟均方根误差")
    ax1.set_xlabel("参考时效 (h)")
    ax1.set_xticks([1, 2, 4, 6])
    ax1.set_ylabel("误差 (dB)")
    ax1.set_ylim(0, None)
    ax1.legend(fontsize=6.5, loc="lower right")
    ax1.set_title("(a) 干期检验：参考漂移", loc="left")
    hb = ax2.hexbin(ev.pred_ref_q1, ev.eps_ref_q1, gridsize=45, bins="log", cmap="Greys", mincnt=1, linewidths=0)
    lim = 0.55
    ax2.plot([-lim, lim], [-lim, lim], color=C["q1"], lw=0.8)
    ax2.set_xlim(-lim, lim), ax2.set_ylim(-lim, lim)
    ax2.set_xlabel(r"预测的参考量化偏差 $b(c;\sigma)$ (dB)")
    ax2.set_ylabel("实际的参考量化误差 (dB)")
    if S:
        ax2.text(0.04, 0.93, f"斜率 {S['q1_slope']}，$R^2$ = {S['q1_r2']}", transform=ax2.transAxes, fontsize=7)
    ax2.set_title("(b) 1 dB 报告：逐事件检验", loc="left")
    cb = fig.colorbar(hb, ax=ax2, pad=0.02, fraction=0.05)
    cb.set_label("事件数", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    save(fig, out, "fig3_drift_quant")


def fig_pilot(out):
    p = pd.read_csv(RES / "pilot_reference.csv")
    S = json.loads((RES / "summary.json").read_text())
    J = int(p.J.iloc[0])
    a2 = 10 ** (-0.5 / 10)
    fig, ax = plt.subplots(figsize=(4.2, 2.9))
    ms = np.logspace(np.log10(0.125), 2, 100)
    for snr, col in zip((5.0, 15.0, 25.0), ("#4b3f8f", "#0f6e6e", "#c26a1b")):
        g_r = (ms / 0.125) * 10 ** (snr / 10)
        S_c, _ = complete_shared_independent(g_r, a2 * 10 ** (snr / 10), J)
        ax.plot(ms, np.sqrt(S_c), color=col, lw=1, label=f"{snr:g} dB，闭式")
        d = p[p.snr_db == snr]
        lo = np.sqrt(d.S_lo.clip(lower=1e-12))
        hi = np.sqrt(d.S_hi)
        ax.errorbar(d.reference_ms, np.sqrt(d.S_sim), yerr=[np.sqrt(d.S_sim) - lo, hi - np.sqrt(d.S_sim)], fmt="o",
                    color=col, ms=3, lw=0.8, capsize=1.5)
    for val, lab, ls, va in ((float(S["q1_tot_rms"]), "1 dB 报告：事件级量化误差", "--", "top"),
                             (float(S["drift_1h"]), "参考漂移（1 h）", ":", "bottom"),
                             (float(S["W_at_1"]), r"长度无关偏移 $W$（1 mm/h）", "-.", "bottom")):
        ax.axhline(val, color="#555555", lw=0.7, ls=ls)
        ax.text(95, val * (0.88 if va == "top" else 1.1), lab, fontsize=6.5, color="#333333", ha="right", va=va)
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xlim(0.1, 100), ax.set_ylim(5e-4, 3)
    ax.set_xlabel("晴空参考时长 (ms)")
    ax.set_ylabel("共享参考误差标准差 (dB)")
    ax.legend(fontsize=6.5, loc="lower left", title=f"导频，$J$ = {J}", title_fontsize=6.5)
    save(fig, out, "fig4_pilot")


def fig_evaluation(out):
    S = json.loads((RES / "summary.json").read_text())
    tab = pd.DataFrame(S["eval_radar"])
    cls_order = list(dict.fromkeys(tab.cls))
    series = [("radar", "none", "不扣除偏移", C["const"], "o"), ("radar", "constant", "常数偏移", C["offset"], "s"),
              ("radar", "rain", r"分档 $W(R)$", "#d9a441", "D"), ("radar", "kralt", "KR-alt", C["rain"], "v"),
              ("link", "kralt", "KR-alt，链路自判干湿", "#4b8fb8", "^"), ("radar_q1", "kralt", "KR-alt，1 dB 报告", C["q1"], "x")]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.3, 2.75), gridspec_kw={"wspace": 0.32, "width_ratios": [1.25, 1]})
    x = np.arange(len(cls_order))
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    for i, (mask, corr, lab, col, mk) in enumerate(series):
        d = tab[(tab["mask"] == mask) & (tab["corr"] == corr)].set_index("cls").reindex(cls_order)
        ax1.plot(x + (i - 2.5) * 0.07, d.rmse, marker=mk, color=col, ms=3, lw=0.8, label=lab)
    ax1.set_xticks(x, ["≤0.5", "0.5–1", "1–2", "2–4", ">4"], fontsize=7.5)
    ax1.set_xlabel("链路长度 (km)")
    ax1.set_yscale("log")
    ax1.yaxis.set_major_locator(FixedLocator([0.2, 0.3, 0.5, 1, 2, 3, 5]))
    ax1.yaxis.set_minor_locator(NullLocator())
    ax1.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax1.set_ylabel("小时雨强均方根误差 (mm/h)")
    ax1.legend(fontsize=5.8, loc="upper right", borderaxespad=0.2)
    ax1.set_title("(a) 以雷达为参考（8月）", loc="left")
    g = pd.DataFrame(S["eval_gauge"])
    bars = [("radar", "none", "不扣除", C["const"]), ("radar", "constant", "常数", C["offset"]),
            ("radar", "rain", "分档W(R)", "#d9a441"), ("radar", "kralt", "KR-alt", C["rain"]), ("-", "radar", "雷达", C["radar"])]
    width = 0.15
    for i, (mask, corr, lab, col) in enumerate(bars):
        vals = [g[(g["mask"] == mask) & (g["corr"] == corr) & (g.c2 == c)].rmse.iloc[0] for c in ("le1", "gt1")]
        ax2.bar(np.arange(2) + (i - 2) * width, vals, width, color=col, label=lab)
    ax2.set_xticks([0, 1], [f"≤1 km（{S['g_n_le1']}条）", f">1 km（{S['g_n_gt1']}条）"], fontsize=7)
    ax2.set_yscale("log")
    ax2.yaxis.set_major_locator(FixedLocator([0.2, 0.3, 0.5, 1, 2, 3]))
    ax2.yaxis.set_minor_locator(NullLocator())
    ax2.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax2.set_ylabel("小时雨强均方根误差 (mm/h)")
    ax2.legend(fontsize=5.8, loc="upper right", borderaxespad=0.2)
    ax2.set_title("(b) 以雨量计为参考（8月）", loc="left")
    save(fig, out, "fig5_evaluation")


def fig_region2(out):
    from evaluate_corrections import kralt
    S = json.loads((RES / "summary.json").read_text())
    R2 = RES / "openrainer"
    links = pd.read_csv(R2 / "links.csv")
    links = links[links.wet_min >= CFG["min_wet_minutes"]]
    off = json.loads((R2 / "offset.json").read_text())
    fig, axes = plt.subplots(2, 2, figsize=(6.3, 5.0), gridspec_kw={"wspace": 0.32, "hspace": 0.45})
    ax1, ax2, ax3, ax4 = axes.ravel()
    edges = CFG["openrainer"]["length_classes_km"]
    labels = ["≤2 km", "2–5 km", "5–10 km", ">10 km"]
    cls = np.digitize(links.length_km.values, edges[1:-1], right=True)
    for c in range(4):
        s_ = cls == c
        ax1.scatter(links.exp_mean[s_], links.att_mean[s_], s=5, color=CLASS_COLORS[[0, 2, 3, 4][c]], label=labels[c], lw=0)
    x = np.linspace(0, links.exp_mean.max() * 1.05, 10)
    ax1.plot(x, off["W"] + off["beta"] * x, color="k", lw=0.8)
    ax1.annotate(rf"$W$ = {off['W']:.2f} dB", xy=(0.0, off["W"]), xytext=(1.6, 0.12), fontsize=7,
                 arrowprops=dict(arrowstyle="->", lw=0.6, color="k"))
    ax1.set_xlabel("雷达对应雨致衰减 (dB)")
    ax1.set_ylabel("链路湿期平均衰减 (dB)")
    ax1.set_xlim(0, None), ax1.set_ylim(0, None)
    ax1.legend(fontsize=6, loc="upper left", handletextpad=0.1, markerscale=2, borderaxespad=0.2)
    ax1.set_title("(a) OpenRainER跨链路回归", loc="left", fontsize=8.5)
    mrg = pd.DataFrame(S["curve"]).sort_values("bin")
    er = pd.DataFrame(S["r2_curve"]).sort_values("bin")
    ax2.fill_between(mrg.R, mrg.W_lo, mrg.W_hi, color=C["offset"], alpha=0.2, lw=0)
    ax2.plot(mrg.R, mrg.W, "o-", color=C["offset"], ms=2.5, lw=1, label="OpenMRG")
    ax2.fill_between(er.R, er.W_lo, er.W_hi, color=C["radar"], alpha=0.2, lw=0)
    ax2.plot(er.R, er.W, "s-", color=C["radar"], ms=2.5, lw=1, label="OpenRainER")
    R = np.logspace(-2, np.log10(30), 200)
    ax2.plot(R, kralt(CFG["extension"]["kralt"])(R), color=C["rain"], lw=1, ls="--", label="KR-alt（文献参数）")
    ax2.set_xscale("log")
    ax2.set_ylim(0, 4.5), ax2.set_xlim(5e-3, 35)
    ax2.set_xlabel("雷达雨强 $R$ (mm/h)")
    ax2.set_ylabel("长度无关偏移 (dB)")
    ax2.legend(fontsize=6, loc="upper left", borderaxespad=0.2)
    ax2.set_title("(b) 偏移随雨强的变化", loc="left", fontsize=8.5)
    for key, lab, col, mk in (("drying_icpt_mrg", "OpenMRG（5 min分档）", C["offset"], "o"),
                              ("drying_icpt_er", "OpenRainER（15 min分档）", C["radar"], "s")):
        d = pd.DataFrame(S[key])
        ax3.plot(d.lag_min, d.W, marker=mk, ms=2.5, lw=1, color=col, label=lab)
    ax3.axhline(0, color="k", lw=0.4)
    ax3.set_xlim(0, 90)
    ax3.set_xlabel("雷达雨停后的时间 (min)")
    ax3.set_ylabel("长度无关残余衰减 (dB)")
    ax3.legend(fontsize=6, loc="upper right", borderaxespad=0.2)
    ax3.set_title("(c) 雨停后的消退（跨链路截距）", loc="left", fontsize=8.5)
    t = pd.DataFrame(S["r2_eval"])
    t = t[t.ref == "gauge"]
    bars = [("none", "不扣除", C["const"]), ("constant", "常数", C["offset"]), ("rain", "分档W(R)", "#d9a441"),
            ("kralt", "KR-alt", C["rain"]), ("radar", "雷达", C["radar"])]
    width = 0.15
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    for i, (corr, lab, col) in enumerate(bars):
        vals = [t[(t["corr"] == corr) & (t.c2 == c)].rmse.iloc[0] for c in ("le2", "gt2")]
        ax4.bar(np.arange(2) + (i - 2) * width, vals, width, color=col, label=lab)
    ax4.set_xticks([0, 1], [f"≤2 km（{S['r2_g_n_le2']}条）", f">2 km（{S['r2_g_n_gt2']}条）"], fontsize=7)
    ax4.set_yscale("log")
    ax4.set_ylim(0.2, 2.6)
    ax4.yaxis.set_major_locator(FixedLocator([0.2, 0.3, 0.5, 1, 2]))
    ax4.yaxis.set_minor_locator(NullLocator())
    ax4.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax4.set_ylabel("小时雨强均方根误差 (mm/h)")
    ax4.legend(fontsize=6, loc="upper right", borderaxespad=0.2)
    ax4.set_title("(d) 2022年5–9月，雨量计参考", loc="left", fontsize=8.5)
    save(fig, out, "fig6_region2")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "paper" / "figures_zh")
    ap.add_argument("--only", nargs="*")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    todo = {"offset": fig_offset, "drift": fig_drift_quant, "pilot": fig_pilot, "evaluation": fig_evaluation,
            "example": fig_example, "region2": fig_region2}
    info = {}
    for name, fn in todo.items():
        if args.only and name not in args.only:
            continue
        r = fn(args.out)
        if r:
            info[name] = r
        print("figure", name, "done", flush=True)
    if info:
        (RES / "figures_info.json").write_text(json.dumps(info, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
