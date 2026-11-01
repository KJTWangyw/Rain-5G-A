#!/usr/bin/env python3
"""Figures of the manuscript (600 dpi PNG for the Word build, PDF for reference).

Usage:  python experiments/make_figures.py [--results results] [--out paper/figures] [--lang en|zh]
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
from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402
from matplotlib.ticker import NullFormatter  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nrpilot.uncertainty import complete_shared_independent, crb  # noqa: E402

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 8.5, "axes.linewidth": 0.6, "axes.titlesize": 9,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.frameon": False, "savefig.dpi": 600,
})
COLORS = {"independent": "#7f7f7f", "local": "#c26a1b", "complete": "#0f6e6e", "bootstrap": "#4b3f8f",
          "crb": "#000000", "empirical": "#b2182b"}
SNR_COLORS = {0: "#0f6e6e", 5: "#4b3f8f", 15: "#c26a1b", 25: "#7f7f7f"}


TEXT = {
    "en": {
        "independent": "Independent (wrong)", "local": "Fixed-weight", "complete": "Complete (closed form)",
        "bootstrap": "Parametric bootstrap", "crb": "Information bound", "simulated": "Simulated",
        "f1_x": "OFDM symbol in slot", "f1_y": "Subcarrier (2 of 12 RB shown)", "f1_a": "(a) DM-RS type 1 in one slot",
        "f1_b": "(b) One reference reused for $n$ rain observations", "f1_ref": "$M$ clear-sky slots",
        "f1_rain": "rain slots $t = 1, \\ldots, n$:  $z_{tj} = a_t e^{i\\phi_t} h_j + e_{tj}$",
        "f1_shared": "shared reference error $e_{0j}$ enters every $\\hat{A}_t$",
        "f1_mean": "$\\bar{A}_n = \\frac{1}{n}\\sum_t \\hat{A}_t$:  current noise $e_{tj}$ averages out, $e_{0j}$ does not",
        "f2_x": "Number of averaged observations $n$",
        "f3_local": "Fixed-weight (expected)", "f3_complete": "Complete, {snr} dB",
        "f3_x": "Per-pilot reference SNR $\\gamma_{\\mathrm{r}}$", "f3_y": "Shared variance / bound",
        "f3_a": "(a) Shared component", "f3_main": "main cell", "f3_bx": "Reference slots $M$",
        "f3_by": "Averaged observations $n$", "f3_b": "(b) Shared share of the bound, {att} dB",
        "f4_x": "SNR (dB)", "f4_y": "Reference slots $M$", "f4_cb": "95% coverage (%)",
        "f5_x": "Nominal level (%)", "f5_y": "Empirical coverage (%)",
        "f5_a": "(a) Calibration, {snr} dB, $M$ = {M}, $n$ = {n}",
        "f5_ticks": ["Indep.", "Local", "Complete", "Bootstrap"], "f5_by": "Mean width of 95% interval (dB)",
        "f5_b": "(b) Width and coverage at 95%",
    },
    "zh": {
        "independent": "独立假设（错误）", "local": "固定权重线性化", "complete": "完整闭式传播",
        "bootstrap": "参数化自助法", "crb": "信息下界", "simulated": "仿真",
        "f1_x": "时隙内OFDM符号序号", "f1_y": "子载波（显示12个RB中的2个）", "f1_a": "(a) 单个时隙中的类型1 DM-RS",
        "f1_b": "(b) 同一参考复用于$n$次降雨观测", "f1_ref": "$M$个晴空时隙",
        "f1_rain": "降雨时隙$t = 1, \\ldots, n$：$z_{tj} = a_t e^{i\\phi_t} h_j + e_{tj}$",
        "f1_shared": "共享参考误差$e_{0j}$进入每个$\\hat{A}_t$",
        "f1_mean": "$\\bar{A}_n = \\frac{1}{n}\\sum_t \\hat{A}_t$：当前噪声$e_{tj}$被平均，$e_{0j}$不会",
        "f2_x": "平均观测次数$n$",
        "f3_local": "固定权重线性化（期望）", "f3_complete": "完整传播，{snr} dB",
        "f3_x": "单导频参考信噪比$\\gamma_{\\mathrm{r}}$", "f3_y": "共享方差 / 信息下界",
        "f3_a": "(a) 共享分量", "f3_main": "主条件", "f3_bx": "参考时隙数$M$",
        "f3_by": "平均观测次数$n$", "f3_b": "(b) 共享项在下界中的占比，{att} dB",
        "f4_x": "信噪比 (dB)", "f4_y": "参考时隙数$M$", "f4_cb": "95%区间覆盖率 (%)",
        "f5_x": "名义水平 (%)", "f5_y": "实际覆盖率 (%)",
        "f5_a": "(a) 校准，{snr} dB，$M$ = {M}，$n$ = {n}",
        "f5_ticks": ["独立假设", "固定权重", "完整闭式", "自助法"], "f5_by": "95%区间平均宽度 (dB)",
        "f5_b": "(b) 95%区间的宽度与覆盖率",
    },
}
T = TEXT["en"]


def set_language(lang):
    global T
    T = TEXT[lang]
    if lang == "zh":
        # mathtext draws the non-math part of a label with the first family only, so the CJK face goes first;
        # math symbols stay in Times New Roman
        plt.rcParams.update({
            "font.family": ["Noto Serif CJK SC", "Times New Roman"], "mathtext.fontset": "custom",
            "mathtext.rm": "Times New Roman", "mathtext.it": "Times New Roman:italic",
            "mathtext.bf": "Times New Roman:bold", "mathtext.sf": "Times New Roman", "mathtext.tt": "Times New Roman",
            "mathtext.cal": "Times New Roman:italic", "mathtext.fallback": "stix",
        })


def save(fig, out, name):
    fig.savefig(out / f"{name}.png", bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_schematic(out):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.3, 2.5), gridspec_kw={"width_ratios": [1, 1.9], "wspace": 0.15})
    for l in range(14):
        for k in range(24):
            pilot = l in (2, 11) and k % 2 == 0
            ax1.add_patch(Rectangle((l, k), 1, 1, facecolor=COLORS["complete"] if pilot else "#f2f2f2",
                                    edgecolor="white", linewidth=0.4))
    ax1.set_xlim(0, 14), ax1.set_ylim(0, 24)
    ax1.set_xticks([2.5, 11.5], ["2", "11"])
    ax1.set_yticks([0.5, 12.5, 23.5], ["0", "12", "23"])
    ax1.set_xlabel(T["f1_x"])
    ax1.set_ylabel(T["f1_y"])
    ax1.set_title(T["f1_a"], loc="left")
    ax1.tick_params(length=0)
    for s in ax1.spines.values():
        s.set_visible(False)

    ax2.set_xlim(0, 10), ax2.set_ylim(0, 4.6)
    ax2.axis("off")
    ax2.set_title(T["f1_b"], loc="left")
    for i in range(3):
        ax2.add_patch(Rectangle((0.25 + 0.5 * i, 3.2), 0.42, 0.6, facecolor="#9ecae1", edgecolor="#3a6b8a", lw=0.6))
    ax2.text(1.0, 3.95, T["f1_ref"], ha="center", fontsize=7.5)
    ax2.text(5.8, 3.95, T["f1_rain"], ha="center", fontsize=7.5)
    ax2.add_patch(Rectangle((0.05, 1.9), 1.9, 0.6, facecolor="white", edgecolor=COLORS["empirical"], lw=0.9))
    ax2.text(1.0, 2.2, "$r_j = h_j + e_{0j}$", ha="center", va="center", fontsize=8)
    ax2.add_patch(FancyArrowPatch((1.0, 3.15), (1.0, 2.55), arrowstyle="-|>", mutation_scale=7, lw=0.7))
    xs = [3.0, 4.6, 6.2, 8.6]
    ax2.plot([1.0, 1.0, xs[-1]], [1.85, 1.5, 1.5], color=COLORS["empirical"], lw=0.9)
    for k, x in enumerate(xs):
        ax2.add_patch(Rectangle((x - 0.45, 3.2), 0.9, 0.6, facecolor="#c7e9c0", edgecolor="#3f7f3f", lw=0.6))
        ax2.add_patch(FancyArrowPatch((x, 3.15), (x, 2.55), arrowstyle="-|>", mutation_scale=7, lw=0.7))
        ax2.add_patch(FancyArrowPatch((x, 1.5), (x, 1.88), arrowstyle="-|>", mutation_scale=7, lw=0.8,
                                      color=COLORS["empirical"]))
        ax2.add_patch(Rectangle((x - 0.45, 1.9), 0.9, 0.6, facecolor="white", edgecolor="#3f7f3f", lw=0.6))
        ax2.text(x, 2.2, "$\\hat{A}_n$" if k == 3 else f"$\\hat{{A}}_{k + 1}$", ha="center", va="center", fontsize=8)
    ax2.text(7.4, 2.2, "…", ha="center", va="center", fontsize=10)
    ax2.text(4.8, 1.1, T["f1_shared"], ha="center", fontsize=7.5, color=COLORS["empirical"])
    ax2.text(4.8, 0.35, T["f1_mean"], ha="center", fontsize=7.8)
    save(fig, out, "fig1_schematic")


def fig_variance(res, out, cfg):
    var = res["variance"][res["variance"].scenario == "constant"]
    comp = res["components"]
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.6), sharey=False)
    for ax, M, tag in zip(axes, (2, 64), ("(a)", "(b)")):
        snr = cfg["main_cell"]["snr_db"]
        d = var[(var.snr_db == snr) & (var.M == M)].sort_values("n")
        c = comp[(comp.snr_db == snr) & (comp.M == M)].iloc[0]
        n = np.logspace(0, np.log10(32), 60)
        S, U = complete_shared_independent(c.g_r, c.g_t, 144)
        ax.plot(n, S + U / n, color=COLORS["complete"], lw=1.2, label=T["complete"])
        ax.plot(n, crb(c.g_r, c.g_t, 144, n)[2], color=COLORS["crb"], lw=0.9, ls="--", label=T["crb"])
        ax.plot(d.n, d.mean_v_local, "s-", color=COLORS["local"], ms=3.2, lw=0.9, label=T["local"])
        ax.plot(d.n, d.mean_v_independent, "^:", color=COLORS["independent"], ms=3.2, lw=0.9,
                label=T["independent"])
        ax.errorbar(d.n, d.emp_var, yerr=[d.emp_var - d.emp_var_lo, d.emp_var_hi - d.emp_var], fmt="o",
                    color=COLORS["empirical"], ms=3.5, capsize=2, lw=0.8, label=T["simulated"])
        ax.set_xscale("log"), ax.set_yscale("log")
        ax.set_ylim(0.008, 1.2)
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_xticks([1, 4, 16, 32], ["1", "4", "16", "32"])
        ax.set_xlabel(T["f2_x"])
        ax.set_title(f"{tag} {snr} dB, $M$ = {M}", loc="left")
    axes[0].set_ylabel(r"Var($\bar{A}_n$) (dB$^2$)")
    axes[0].legend(fontsize=6.8, loc="lower left")
    fig.tight_layout()
    save(fig, out, "fig2_variance")


def fig_shared(res, out, cfg):
    comp = res["components"]
    a2 = 10.0 ** (-cfg["attenuation_db"] / 10.0)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.3, 2.7), gridspec_kw={"wspace": 0.38})
    g = np.logspace(np.log10(0.8), np.log10(2.5e4), 300)
    ax1.plot(g, g / (g + 1.0), color=COLORS["local"], lw=1.0, label=T["f3_local"])
    ax1.axhline(1.0, color="black", lw=0.8, ls="--", label=T["crb"])
    for snr, col in SNR_COLORS.items():
        g_t = a2 * 10 ** (snr / 10)
        gg = g[g >= 10 ** (snr / 10) * 2 * 0.8]
        ax1.plot(gg, 1 + 2 * gg / (gg + g_t) ** 2, color=col, lw=1.0, label=T["f3_complete"].format(snr=snr))
        d = comp[comp.snr_db == snr]
        ax1.errorbar(d.g_r, d.S_sim / d.S_crb, yerr=[(d.S_sim - d.S_sim_lo) / d.S_crb, (d.S_sim_hi - d.S_sim) / d.S_crb],
                     fmt="o", color=col, ms=3.2, capsize=1.5, lw=0.7)
        ax1.plot(d.g_r, d.S_local_plugin_mean / d.S_crb, "s", color=COLORS["local"], ms=3.0, mfc="white", mew=0.8)
    ax1.set_xscale("log")
    ax1.set_ylim(0.3, 2.0)
    ax1.set_xlabel(T["f3_x"])
    ax1.set_ylabel(T["f3_y"])
    ax1.set_title(T["f3_a"], loc="left")
    ax1.legend(fontsize=6.2, loc="upper right", ncol=1)

    M = np.logspace(0, 4, 300)
    n = np.logspace(0, 4, 300)
    MM, NN = np.meshgrid(M, n)
    share = NN * a2 / (MM + NN * a2)
    im = ax2.pcolormesh(MM, NN, share, shading="auto", cmap="viridis", vmin=0, vmax=1, rasterized=True)
    levels = [0.1, 0.5, 0.9]
    cs = ax2.contour(MM, NN, share, levels=levels, colors="white", linewidths=0.7)
    at = [(Mv, s * Mv / ((1 - s) * a2)) for s, Mv in zip(levels, (3000, 200, 10))]
    ax2.clabel(cs, fmt=lambda x: f"{x:.1f}", fontsize=6.5, manual=at)
    for Mv, lab, off in ((2, T["f3_main"], (5, 4)), (64, "$M$ = 64", (5, -10))):
        ax2.plot(Mv, 32, "*" if Mv == 2 else "o", color="#f97316" if Mv == 2 else "white", ms=8 if Mv == 2 else 4,
                 mec="white", mew=0.5)
        ax2.annotate(lab, (Mv, 32), textcoords="offset points", xytext=off, color="white", fontsize=6.8)
    ax2.set_xscale("log"), ax2.set_yscale("log")
    ax2.set_xlabel(T["f3_bx"])
    ax2.set_ylabel(T["f3_by"])
    ax2.set_title(T["f3_b"].format(att=f"{cfg['attenuation_db']:g}"), loc="left")
    fig.colorbar(im, ax=ax2, fraction=0.05, pad=0.03)
    save(fig, out, "fig3_shared")


def fig_coverage(res, out, cfg):
    cov = res["coverage"]
    n_max = max(cfg["n_avg"])
    sel = cov[(cov.scenario == "constant") & (cov.n == n_max) & np.isclose(cov.level, 0.95)]
    methods = ("independent", "local", "complete", "bootstrap")
    fig, axes = plt.subplots(1, 4, figsize=(6.6, 2.0), sharey=True, gridspec_kw={"wspace": 0.08})
    for ax, m in zip(axes, methods):
        grid = sel[sel.method == m].pivot(index="M", columns="snr_db", values="coverage") * 100
        im = ax.imshow(grid.to_numpy(), cmap="viridis", vmin=30, vmax=100, aspect="auto")
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                val = grid.iat[i, j]
                ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=6.5,
                        color="white" if val < 75 else "black")
        ax.set_xticks(range(grid.shape[1]), [f"{c:g}" for c in grid.columns])
        ax.set_yticks(range(grid.shape[0]), [str(i) for i in grid.index])
        ax.set_xlabel(T["f4_x"])
        ax.set_title(T[m], fontsize=7.5)
        ax.tick_params(length=0)
    axes[0].set_ylabel(T["f4_y"])
    cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cb.set_label(T["f4_cb"], fontsize=7)
    save(fig, out, "fig4_coverage")


def fig_calibration(res, out, cfg):
    cov = res["coverage"]
    mc = cfg["main_cell"]
    n_max = max(cfg["n_avg"])
    d = cov[(cov.scenario == "constant") & (cov.snr_db == mc["snr_db"]) & (cov.M == mc["M"]) & (cov.n == n_max)]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.3, 2.6), gridspec_kw={"wspace": 0.35})
    levels = np.array(sorted(d.level.unique()))
    ax1.plot([45, 100], [45, 100], color="black", lw=0.6, ls="--")
    for m in ("independent", "local", "complete", "bootstrap"):
        e = d[d.method == m].sort_values("level")
        ax1.errorbar(100 * e.level, 100 * e.coverage, yerr=[100 * (e.coverage - e.wilson_lo), 100 * (e.wilson_hi - e.coverage)],
                     fmt="o-", ms=3, lw=0.9, capsize=1.5, color=COLORS[m], label=T[m])
    ax1.set_xlim(45, 100), ax1.set_ylim(0, 100)
    ax1.set_xticks(100 * levels)
    ax1.set_xlabel(T["f5_x"])
    ax1.set_ylabel(T["f5_y"])
    ax1.set_title(T["f5_a"].format(snr=mc["snr_db"], M=mc["M"], n=n_max), loc="left")
    ax1.legend(fontsize=6.5, loc="upper left")

    w = d[np.isclose(d.level, 0.95)].set_index("method")
    order = ["independent", "local", "complete", "bootstrap"]
    ax2.bar(range(4), [w.loc[m, "mean_width_db"] for m in order], color=[COLORS[m] for m in order], width=0.6)
    for k, m in enumerate(order):
        ax2.text(k, w.loc[m, "mean_width_db"] + 0.03, f"{100 * w.loc[m, 'coverage']:.1f}%", ha="center", fontsize=7)
    ax2.set_xticks(range(4), T["f5_ticks"])
    ax2.set_ylabel(T["f5_by"])
    ax2.set_title(T["f5_b"], loc="left")
    ax2.set_ylim(0, max(w["mean_width_db"]) * 1.2)
    save(fig, out, "fig5_calibration")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=ROOT / "results")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--lang", choices=sorted(TEXT), default="en")
    args = ap.parse_args()
    set_language(args.lang)
    args.out = args.out or ROOT / "paper" / ("figures" if args.lang == "en" else f"figures_{args.lang}")
    args.out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((ROOT / "experiments" / "config.json").read_text())
    res = {k: pd.read_csv(args.results / f"{k}.csv") for k in ("variance", "components", "coverage")}
    fig_schematic(args.out)
    fig_variance(res, args.out, cfg)
    fig_shared(res, args.out, cfg)
    fig_coverage(res, args.out, cfg)
    fig_calibration(res, args.out, cfg)
    print("figures written to", args.out)


if __name__ == "__main__":
    main()
