#!/usr/bin/env python3
"""Collect every number quoted in the manuscript into results/summary.json.

Usage:  python experiments/summarize.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nrpilot.rain import p838_horizontal, rain_rate  # noqa: E402
from nrpilot.uncertainty import complete_components, complete_shared_independent, crb, local_expected  # noqa: E402


def f(x, d=3):
    return f"{x:.{d}f}"


def pct(x, d=1):
    return f"{100 * x:.{d}f}"


def pct_range(lo, hi, d=1):
    """'a%' when both ends round equally, otherwise 'a–b%'."""
    a, b = pct(lo, d), pct(hi, d)
    return f"{a}%" if a == b else f"{a}–{b}%"


def p_value(p):
    """Inline-math p-value statement: '$p$ = 0.012', '$p$ = 3.4 × 10^-8' or '$p$ < 10^-300'."""
    if p >= 1e-3:
        return f"$p$ = {p:.3f}" if p < 0.995 else "$p$ = 1.00"
    if p < 1e-300:
        return "$p < 10^{-300}$"
    mant, exp = f"{p:.1e}".split("e")
    return f"$p$ = ${mant}\\times 10^{{{int(exp)}}}$"


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=ROOT / "results")
    res = ap.parse_args().results
    cfg = json.loads((ROOT / "experiments" / "config.json").read_text())
    cov = pd.read_csv(res / "coverage.csv")
    var = pd.read_csv(res / "variance.csv")
    comp = pd.read_csv(res / "components.csv")
    pair = pd.read_csv(res / "paired.csv")
    timing = json.loads((res / "timing.json").read_text())
    info = json.loads((res / "run_info.json").read_text())
    J = info["J"]
    snr0, M0 = cfg["main_cell"]["snr_db"], cfg["main_cell"]["M"]
    n_max = max(cfg["n_avg"])
    A = cfg["attenuation_db"]
    a2 = 10 ** (-A / 10)
    s = {}

    # Link and rain mapping (ITU-R P.838-3)
    k, alpha = p838_horizontal(cfg["link"]["carrier_ghz"])
    L = cfg["link"]["length_km"]
    R = rain_rate(A, k, alpha, L)
    s.update({"k_h": f(k, 4), "alpha_h": f(alpha, 4), "rain_rate": f(R, 1), "J": str(J),
              "refs": str(info["refs_per_cell"]), "draws": str(info["bootstrap_draws"]),
              "n_max": str(n_max), "A_db": f"{A:g}", "length_km": f"{L:g}"})

    # Analytic components at the main cell
    v = 10 ** (-snr0 / 10)
    g_r, g_t = J / (J * v / M0), a2 * J / (J * v)
    c = complete_components(g_r, g_t, J)
    S_c, U_c = complete_shared_independent(g_r, g_t, J)
    S_b, U_b, V_b = crb(g_r, g_t, J, n_max)
    S_l, U_l = local_expected(g_r, g_t, J)
    s.update({
        "g_r": f(g_r, 0), "g_t": f(g_t, 2),
        "crb_S": f(S_b), "crb_U": f(U_b), "crb_V": f(V_b),
        "loc_S": f(S_l), "loc_U": f(U_l), "loc_V": f(S_l + U_l / n_max),
        "comp_S": f(S_c), "comp_Sx": f(c["s_x"]), "comp_Sq": f(c["s_q"]), "comp_U": f(U_c),
        "comp_V": f(S_c + U_c / n_max),
        "loc_over_crb": f(S_l / S_b, 2), "comp_over_crb": f(S_c / S_b, 2),
        "loc_under_pct": pct(1 - S_l / S_c, 0), "loc_under_total_pct": pct(1 - (S_l + U_l / n_max) / (S_c + U_c / n_max), 0),
        "share_main": pct(n_max * a2 / (M0 + n_max * a2), 0), "share_m64": pct(n_max * a2 / (64 + n_max * a2), 0),
        "n_cross_main": f(M0 / a2, 0), "n_cross_m64": f(64 / a2, 0),
    })

    # Simulated variance at the main cell
    main_var = var[(var.scenario == "constant") & (var.snr_db == snr0) & (var.M == M0) & (var.n == n_max)].iloc[0]
    main_comp = comp[(comp.snr_db == snr0) & (comp.M == M0)].iloc[0]
    s.update({
        "emp_V": f(main_var.emp_var), "emp_V_lo": f(main_var.emp_var_lo), "emp_V_hi": f(main_var.emp_var_hi),
        "mean_v_local": f(main_var.mean_v_local), "mean_v_complete": f(main_var.mean_v_complete),
        "mean_v_indep": f(main_var.mean_v_independent), "mean_v_boot": f(main_var.mean_v_bootstrap_normal),
        "sim_S": f(main_comp.S_sim), "sim_S_lo": f(main_comp.S_sim_lo), "sim_S_hi": f(main_comp.S_sim_hi),
        "sim_U": f(main_comp.U_sim), "bias_main": f(main_var.bias, 3),
    })
    cc = comp.copy()
    cc["rel"] = cc.S_sim / cc.S_complete - 1
    inside = ((cc.S_sim_lo <= cc.S_complete) & (cc.S_complete <= cc.S_sim_hi)).sum()
    vv = var[var.scenario == "constant"].copy()
    vv["rel"] = vv.emp_var / vv.complete_analytic - 1
    v_inside = ((vv.emp_var_lo <= vv.complete_analytic) & (vv.complete_analytic <= vv.emp_var_hi)).sum()
    s.update({"S_cells_inside": str(int(inside)), "n_cells": str(len(cc)),
              "V_combos_inside": str(int(v_inside)), "n_combos": str(len(vv)),
              "V_rel_max": pct(vv.rel.abs().max(), 1)})

    # Coverage
    def cov_row(method, scenario="constant", snr=snr0, M=M0, n=n_max, level=0.95):
        return cov[(cov.scenario == scenario) & (cov.snr_db == snr) & (cov.M == M) & (cov.n == n)
                   & np.isclose(cov.level, level) & (cov.method == method)].iloc[0]

    for m in ("independent", "local", "complete", "bootstrap", "bootstrap_normal"):
        r = cov_row(m)
        key = {"independent": "ind", "local": "loc", "complete": "comp", "bootstrap": "boot",
               "bootstrap_normal": "bootn"}[m]
        s.update({f"cov_{key}": pct(r.coverage), f"cov_{key}_lo": pct(r.wilson_lo), f"cov_{key}_hi": pct(r.wilson_hi),
                  f"cnt_{key}": str(int(r.covered)), f"w_{key}": f(r.mean_width_db, 2)})
    s["w_ratio_boot"] = f(cov_row("bootstrap").mean_width_db / cov_row("complete").mean_width_db, 2)
    s["w_increase_pct"] = pct(cov_row("complete").mean_width_db / cov_row("local").mean_width_db - 1, 0)

    grid = cov[(cov.scenario == "constant") & np.isclose(cov.level, 0.95)]
    g32 = grid[grid.n == n_max]
    for m, key in (("complete", "comp"), ("local", "loc"), ("independent", "ind"), ("bootstrap", "boot")):
        d = g32[g32.method == m]
        dall = grid[grid.method == m]
        s.update({f"grid_{key}_min": pct(d.coverage.min()), f"grid_{key}_max": pct(d.coverage.max()),
                  f"grid_{key}_range": pct_range(d.coverage.min(), d.coverage.max()),
                  f"gridall_{key}_range": pct_range(dall.coverage.min(), dall.coverage.max())})
        worst = d.loc[d.coverage.idxmin()]
        s[f"grid_{key}_worst_cell"] = f"{worst.snr_db:g} dB, M = {int(worst.M)}"
    alll = cov[(cov.scenario == "constant") & (cov.method == "complete")]
    alll = alll.assign(dev=(alll.coverage - alll.level).abs())
    s["comp_all_dev_max"] = pct(alll.dev.max())
    s["comp_wilson_contains"] = str(int(((alll.wilson_lo <= alll.level) & (alll.level <= alll.wilson_hi)).sum()))
    s["comp_all_count"] = str(len(alll))
    low = g32[(g32.M == 2) & (g32.method == "local")]
    s["loc_m2_range"] = pct_range(low.coverage.min(), low.coverage.max())
    s["levels_dev_max"] = pct(alll[(alll.snr_db == snr0) & (alll.M == M0) & (alll.n == n_max)].dev.max())

    p = pair[(pair.snr_db == snr0) & (pair.M == M0) & (pair.comparator == "local")].iloc[0]
    s.update({"mcn_comp_only": str(int(p.method_only)), "mcn_loc_only": str(int(p.comparator_only)),
              "mcn_p": p_value(p.p_exact)})
    pb = pair[(pair.snr_db == snr0) & (pair.M == M0) & (pair.comparator == "bootstrap")].iloc[0]
    s.update({"mcnb_comp_only": str(int(pb.method_only)), "mcnb_boot_only": str(int(pb.comparator_only)),
              "mcnb_p": p_value(pb.p_exact)})

    ramp = cov[(cov.scenario == "ramp") & np.isclose(cov.level, 0.95)]
    rc = ramp[ramp.method == "complete"]
    s.update({"ramp_comp_range": pct_range(rc.coverage.min(), rc.coverage.max()),
              "ramp_loc_32": pct(ramp[(ramp.method == "local") & (ramp.n == n_max)].coverage.iloc[0]),
              "ramp_start": f"{cfg['ramp']['start_db']:g}", "ramp_end": f"{cfg['ramp']['end_db']:g}"})

    # Cost
    s.update({"ms_local": f(timing["ms_local"], 2), "ms_complete": f(timing["ms_complete"], 2),
              "ms_boot": f(timing["ms_bootstrap"], 1),
              "speedup": f"{float(f'{timing['ms_bootstrap'] / timing['ms_complete']:.2g}'):.0f}"})

    # Rain-rate floor: shared standard deviation and relative rain-rate error at the configured attenuation
    rows = []
    for snr in cfg["snr_db"][:2]:
        for M in cfg["ref_slots"]:
            vv_ = 10 ** (-snr / 10)
            gr, gt = J / (J * vv_ / M), a2 * J / (J * vv_)
            S_, _ = complete_shared_independent(gr, gt, J)
            rows.append({"snr_db": snr, "M": M, "sd_db": float(np.sqrt(S_)),
                         "rel_rain": float(np.sqrt(S_) / (alpha * A))})
    rr = pd.DataFrame(rows)
    rr.to_csv(res / "rain_floor.csv", index=False)
    s["floor_main_sd"] = f(rr[(rr.snr_db == snr0) & (rr.M == M0)].sd_db.iloc[0], 2)
    s["floor_main_rel"] = pct(rr[(rr.snr_db == snr0) & (rr.M == M0)].rel_rain.iloc[0], 0)
    s["floor_m64_rel"] = pct(rr[(rr.snr_db == snr0) & (rr.M == 64)].rel_rain.iloc[0], 1)

    (res / "summary.json").write_text(json.dumps(s, indent=2, ensure_ascii=False))
    print(json.dumps(s, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
