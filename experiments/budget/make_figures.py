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
    meta, arr = A.load(CFG)
    links = pd.read_csv(RES / "links.csv")
    _, cal = offset_curves()
    Wf = offset_function(cal.R.values, cal.W.values)
    W_const = json.loads((RES / "offset.json").read_text())["W_cal"]["native"]
    month = pd.to_datetime(arr["time"], unit="s").month.values
    # rule: the <=0.5 km sub-link whose August RMSE reduction (no correction -> W(R)) is closest to the class median,
    # then its August event (<= 8 h) with the largest radar accumulation
    evl = pd.read_csv(RES / "evaluation.csv")
    evl = evl[evl.variant == "native"].pivot_table(index="j", columns="correction", values="rmse")
    short = links[(links.length_km <= 0.5) & (links.wet_min >= CFG["min_wet_minutes"])].set_index("j")
    red = (1 - evl["rain"] / evl["none"]).reindex(short.index).dropna()
    j_pick = int((red - red.median()).abs().idxmin())
    best = None
    for j in [j_pick]:
        ev = A.events_of(arr, CFG, j)
        radar = arr["radar_min"][:, j]
        for s, e in zip(ev["starts"], ev["stops"]):
            if month[s] not in CFG["evaluation_months"] or e - s > 480:
                continue
            acc = np.nansum(radar[s:e]) / 60.0
            if best is None or acc > best[0]:
                best = (acc, j, s, e)
    _, j, s, e = best
    m = meta.loc[j]
    ev = A.events_of(arr, CFG, j)
    radar = arr["radar_min"][:, j]
    lo, hi = max(s - 60, 0), min(e + 60, len(radar))
    t = (arr["time"][lo:hi] - arr["time"][s]) / 3600.0
    a = ev["rain_att"][lo:hi]
    att_q = arr["tsl"][:, j].astype(float) - arr["rsl_q1"][:, j].astype(float)
    ref_q = np.nanmean(np.where(ev["dry"][max(s - CFG["reference"]["window_min"], 0):s],
                                att_q[max(s - CFG["reference"]["window_min"], 0):s], np.nan))
    a_q = np.where(np.isfinite(a), att_q[lo:hi] - ref_q, np.nan)
    R = radar[lo:hi]
    exp_att = m.kL * np.where(R > 0, R, 0) ** m.alpha
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.3, 4.4), sharex=True, gridspec_kw={"hspace": 0.28})
    ax1.plot(t, a, color=C["measured"], lw=0.8, label="实测降雨期衰减（原始分辨率）")
    ax1.plot(t, a_q, color=C["q1"], lw=0.7, alpha=0.8, drawstyle="steps-mid", label="1 dB 报告")
    ax1.plot(t, exp_att, color=C["radar"], lw=1.0, label=r"雷达对应的雨致衰减 $kLR^{\alpha}$")
    ax1.plot(t, Wf(R) + exp_att, color=C["offset"], lw=1.0, ls="--", label=r"$W(R)+kLR^{\alpha}$")
    ax1.set_ylabel("衰减 (dB)")
    ax1.legend(loc="upper left", fontsize=7, ncol=2)
    ax1.set_title(f"(a) 子链路 {int(m.Sublink)}：{m.Length_km:.2f} km，{m.Frequency_GHz:.1f} GHz", loc="left")
    inside = np.isfinite(a)
    dry = ev["dry"][lo:hi]
    rates = {}
    for name, fn in (("none", lambda x: rain_rate(x, m.k, m.alpha, m.Length_km, 0.0)),
                     ("const", lambda x: rain_rate(x, m.k, m.alpha, m.Length_km, W_const)),
                     ("rain", lambda x: invert(x, m.kL, m.alpha, Wf))):
        r = np.where(dry, 0.0, np.nan)
        r[inside] = fn(a[inside])
        rates[name] = r
    ax2.plot(t, R, color=C["radar"], lw=1.2, label="雷达路径平均")
    ax2.plot(t, rates["none"], color=C["const"], lw=0.7, label="不扣除偏移")
    ax2.plot(t, rates["const"], color=C["offset"], lw=0.7, label="扣除常数偏移")
    ax2.plot(t, rates["rain"], color=C["rain"], lw=0.9, label=r"扣除 $W(R)$")
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
    return dict(rule_median_reduction=float(red.median()), picked_reduction=float(red.loc[j]), sublink=int(m.Sublink), length_km=float(m.Length_km), freq_ghz=float(m.Frequency_GHz),
                start=str(pd.to_datetime(arr["time"][s], unit="s")), hours=(e - s) / 60.0,
                radar_mm=float(np.nansum(R) / 60.0), link_none_mm=float(np.nansum(rates["none"]) / 60.0),
                link_const_mm=float(np.nansum(rates["const"]) / 60.0), link_rain_mm=float(np.nansum(rates["rain"]) / 60.0))


def fig_offset(out):
    links = pd.read_csv(RES / "links.csv")
    links = links[links.wet_min >= CFG["min_wet_minutes"]]
    off = json.loads((RES / "offset.json").read_text())
    full, cal = offset_curves()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.3, 2.75), gridspec_kw={"wspace": 0.3})
    cls = classes(links.length_km.values)
    for c in range(5):
        s = cls == c
        ax1.scatter(links.exp_mean[s], links.att_mean[s], s=5, color=CLASS_COLORS[c], label=CLASS_LABELS[c], lw=0)
    x = np.linspace(0, links.exp_mean.max() * 1.05, 10)
    ax1.plot(x, off["W"] + off["beta"] * x, color="k", lw=0.8)
    ax1.annotate(rf"截距 $W$ = {off['W']:.2f} dB", xy=(0.0, off["W"]), xytext=(0.55, 0.25), fontsize=7.5,
                 arrowprops=dict(arrowstyle="->", lw=0.6, color="k"))
    ax1.set_xlabel("雷达对应的湿期平均雨致衰减 (dB)")
    ax1.set_ylabel("链路湿期平均衰减 (dB)")
    ax1.set_xlim(0, None), ax1.set_ylim(0, None)
    ax1.legend(fontsize=6.5, loc="upper left", handletextpad=0.2, markerscale=2, bbox_to_anchor=(0.18, 1.0))
    ax1.set_title("(a) 跨链路回归", loc="left")
    ax2.fill_between(full.R, full.W_lo, full.W_hi, color=C["offset"], alpha=0.25, lw=0)
    ax2.plot(full.R, full.W, "o-", color=C["offset"], ms=3, lw=1, label=r"$W(R)$，6–8月")
    ax2.plot(cal.R, cal.W, "s--", color=C["rain"], ms=2.5, lw=0.8, label=r"$W(R)$，6–7月（标定）")
    ax2.plot(full.R, full.W_radar_free, "x", color="k", ms=3.5, label=r"按 $kL$ 回归的截距")
    R = np.logspace(-2, np.log10(30), 200)
    k28, a28 = p838_vertical(28.0)
    for L, ls in ((0.3, ":"), (1.0, "-."), (3.0, "-")):
        y = k28 * L * R ** a28
        ax2.plot(R, y, color=C["radar"], lw=0.7, ls=ls)
        i = min(np.searchsorted(y, 2.9), len(R) - 1)
        ax2.text(R[i] * (0.8 if i < len(R) - 1 else 0.97), min(y[i], 2.9), f"{L:g} km", fontsize=6.5, ha="right",
                 va="top", color=C["radar"], clip_on=True)
    ax2.set_xscale("log")
    ax2.set_ylim(0, 3.2), ax2.set_xlim(5e-3, 35)
    ax2.set_xlabel("雷达路径平均雨强 $R$ (mm/h)")
    ax2.set_ylabel("长度无关偏移 (dB)")
    ax2.legend(fontsize=6.5, loc="upper left")
    ax2.set_title("(b) 偏移随雨强的变化", loc="left")
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
    tab = pd.DataFrame(S["evaluation"])
    cls_order = list(dict.fromkeys(tab.cls))
    series = [("native", "none", "不扣除偏移", C["const"], "o"), ("native", "constant", "常数偏移", C["offset"], "s"),
              ("native", "rain", r"$W(R)$", C["rain"], "D"), ("q1", "rain", r"$W(R)$，1 dB 报告", C["q1"], "v"),
              ("q2", "rain", r"$W(R)$，2 dB 报告", C["q2"], "^")]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.3, 2.75), gridspec_kw={"wspace": 0.3})
    x = np.arange(len(cls_order))
    for i, (var, corr, lab, col, mk) in enumerate(series):
        d = tab[(tab.variant == var) & (tab.correction == corr)].set_index("cls").reindex(cls_order)
        dx = (i - 2) * 0.09
        ax1.plot(x + dx, d.rmse, marker=mk, color=col, ms=3.2, lw=0.8, label=lab)
        ax2.plot(x + dx, d.ratio, marker=mk, color=col, ms=3.2, lw=0.8, label=lab)
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    for ax, ticks in ((ax1, [0.3, 0.5, 1, 2, 3, 5]), (ax2, [1, 2, 3, 5, 10])):
        ax.set_xticks(x, ["≤0.5", "0.5–1", "1–2", "2–4", ">4"], fontsize=7.5)
        ax.set_xlabel("链路长度 (km)")
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(FixedLocator(ticks))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax1.set_ylabel("小时雨强均方根误差 (mm/h)")
    ax2.set_ylabel("链路/雷达 累计雨量比")
    ax2.axhline(1, color="k", lw=0.5)
    ax1.legend(fontsize=6.5, loc="upper right")
    ax1.set_title("(a) 均方根误差（8月，中位数）", loc="left")
    ax2.set_title("(b) 累计雨量比（8月，中位数）", loc="left")
    save(fig, out, "fig5_evaluation")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "paper" / "figures_zh")
    ap.add_argument("--only", nargs="*")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    todo = {"offset": fig_offset, "drift": fig_drift_quant, "pilot": fig_pilot, "evaluation": fig_evaluation,
            "example": fig_example}
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
