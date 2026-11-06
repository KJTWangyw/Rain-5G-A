#!/usr/bin/env python3
"""Monte Carlo study through the NR DM-RS chain.

For every (SNR, M) cell, ``refs_per_cell`` clear-sky references are generated from M coherent slots; each
reference is reused for ``n_obs`` rain observations. For the mean of the first n estimates, five intervals
are formed: wrongly independent, fixed-weight local, complete closed form, parametric bootstrap (basic) and
bootstrap variance with a normal quantile. A ramp scenario (attenuation varying over the window) is run at
one cell. Raw per-trial arrays go to results/raw/ (not versioned); aggregates go to results/.

Usage:  python experiments/run_main.py [--refs N] [--draws B] [--processes P]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nrpilot.bootstrap import bootstrap_means, shrunk_channel  # noqa: E402
from nrpilot.estimator import profile_white  # noqa: E402
from nrpilot.stats import mcnemar_exact, wilson, z_value  # noqa: E402
from nrpilot.uncertainty import (complete_shared_independent, crb, local_expected, local_terms,  # noqa: E402
                                 local_variance, plugin_complete_variance, snr_terms)
from nrpilot.waveform import Multipath, NRConfig, observe_slots, pilot_channel  # noqa: E402

NORMAL_METHODS = ("independent", "local", "complete", "bootstrap_normal")
METHODS = NORMAL_METHODS + ("bootstrap",)


def nr_objects(cfg):
    nr = NRConfig(n_fft=cfg["nr"]["n_fft"], n_cp=cfg["nr"]["n_cp"], scs_hz=cfg["nr"]["scs_hz"],
                  n_rb=cfg["nr"]["n_rb"], dmrs_symbols=tuple(cfg["nr"]["dmrs_symbols"]), n_id=cfg["nr"]["n_id"])
    mp = cfg["multipath"]
    gains = tuple(m * np.exp(1j * p) for m, p in zip(mp["magnitudes"], mp["phases_rad"]))
    return nr, Multipath(delays=tuple(mp["delays"]), gains=gains).taps(nr)


def attenuation_profile(task, n_obs):
    if task["scenario"] == "ramp":
        return np.linspace(task["start_db"], task["end_db"], n_obs)
    return np.full(n_obs, task["attenuation_db"])


def run_cell(task):
    cfg = task["cfg"]
    nr, taps = nr_objects(cfg)
    J = nr.n_pilots
    T, n_avg, levels = cfg["n_obs"], np.array(cfg["n_avg"]), cfg["levels"]
    refs, draws, M = task["refs"], task["draws"], task["M"]
    v = 10.0 ** (-task["snr_db"] / 10.0)
    v0 = v / M
    A_t = attenuation_profile(task, T)
    amps = 10.0 ** (-A_t / 20.0)
    A_true = np.cumsum(A_t)[n_avg - 1] / n_avg
    rng = np.random.default_rng([cfg["seed"], task["cell_id"]])

    est = np.empty((refs, T))
    out = {k: np.empty((refs, len(n_avg))) for k in ("abar", "v_independent", "v_local", "v_complete",
                                                     "v_complete_shared", "v_bootstrap_normal")}
    for L in levels:
        out[f"basic_lo_{L}"] = np.empty((refs, len(n_avg)))
        out[f"basic_hi_{L}"] = np.empty((refs, len(n_avg)))
    local_shared_mean = np.empty(refs)
    for i in range(refs):
        first = i * (M + T)
        r = observe_slots(nr, taps, np.ones(M), np.zeros(M), v, rng, first_slot=first).mean(axis=0)
        z = observe_slots(nr, taps, amps, rng.uniform(0.0, 2.0 * np.pi, T), v, rng, first_slot=first + M)
        A_hat = profile_white(r[None, :], z, v0, v)
        a_hat = 10.0 ** (-A_hat / 20.0)
        est[i] = A_hat
        out["abar"][i] = np.cumsum(A_hat)[n_avg - 1] / n_avg

        nu, C0 = local_terms(r, a_hat, v0, v)
        for k, n in enumerate(n_avg):
            out["v_local"][i, k], out["v_independent"][i, k] = local_variance(nu, C0, n)
        local_shared_mean[i] = C0.mean()
        for k, n in enumerate(n_avg):
            out["v_complete"][i, k], out["v_complete_shared"][i, k] = plugin_complete_variance(r, A_hat, v0, v, n)
        reps = bootstrap_means(rng, shrunk_channel(r, v0), a_hat, v0, v, n_avg, draws)
        out["v_bootstrap_normal"][i] = reps.var(axis=0, ddof=1)
        for L in levels:
            q_lo, q_hi = np.quantile(reps, [0.5 - L / 2.0, 0.5 + L / 2.0], axis=0)
            out[f"basic_lo_{L}"][i] = 2.0 * out["abar"][i] - q_hi
            out[f"basic_hi_{L}"][i] = 2.0 * out["abar"][i] - q_lo
    meta = {k: task[k] for k in ("scenario", "cell_id", "snr_db", "M", "refs", "draws")}
    meta.update({"v": v, "v0": v0, "J": J, "H": float(np.sum(np.abs(pilot_channel(nr, taps)) ** 2)),
                 "A_t": A_t.tolist(), "A_true": A_true.tolist()})
    return meta, est, out, local_shared_mean


def covered_flags(out, A_true, level, method):
    if method == "bootstrap":
        return (out[f"basic_lo_{level}"] <= A_true) & (A_true <= out[f"basic_hi_{level}"])
    half = z_value(level) * np.sqrt(out[f"v_{method}"])
    return np.abs(out["abar"] - A_true) <= half


def widths(out, level, method):
    if method == "bootstrap":
        return out[f"basic_hi_{level}"] - out[f"basic_lo_{level}"]
    return 2.0 * z_value(level) * np.sqrt(out[f"v_{method}"])


def bootstrap_ci(x_fn, est, rng, draws=500):
    idx = rng.integers(0, len(est), size=(draws, len(est)))
    vals = np.array([x_fn(est[j]) for j in idx])
    return np.percentile(vals, [2.5, 97.5])


def shared_cov(est):
    c = np.cov(est, rowvar=False)
    return c[~np.eye(len(c), dtype=bool)].mean(), np.diag(c).mean()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=ROOT / "experiments" / "config.json")
    ap.add_argument("--refs", type=int)
    ap.add_argument("--draws", type=int)
    ap.add_argument("--processes", type=int)
    ap.add_argument("--outdir", type=Path, default=ROOT / "results")
    args = ap.parse_args()
    cfg = json.loads(args.config.read_text())
    refs = args.refs or cfg["refs_per_cell"]
    draws = args.draws or cfg["bootstrap_draws"]
    tasks = []
    for snr in cfg["snr_db"]:
        for M in cfg["ref_slots"]:
            tasks.append({"scenario": "constant", "snr_db": snr, "M": M, "attenuation_db": cfg["attenuation_db"]})
    tasks.append({"scenario": "ramp", "snr_db": cfg["ramp"]["snr_db"], "M": cfg["ramp"]["M"],
                  "start_db": cfg["ramp"]["start_db"], "end_db": cfg["ramp"]["end_db"]})
    for k, t in enumerate(tasks):
        t.update({"cfg": cfg, "cell_id": k, "refs": refs, "draws": draws})

    start = time.time()
    with Pool(args.processes or cfg["processes"]) as pool:
        results = pool.map(run_cell, tasks)
    print(f"simulation finished in {time.time() - start:.0f} s", flush=True)

    outdir = args.outdir
    (outdir / "raw").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(cfg["seed"] + 1)
    n_avg, levels = np.array(cfg["n_avg"]), cfg["levels"]
    cov_rows, var_rows, comp_rows, pair_rows = [], [], [], []
    for meta, est, out, local_shared in results:
        tag = f"{meta['scenario']}_snr{meta['snr_db']}_M{meta['M']}"
        np.savez_compressed(outdir / "raw" / f"{tag}.npz", est=est.astype(np.float32),
                            **{k: v.astype(np.float32) for k, v in out.items()})
        A_true = np.array(meta["A_true"])
        g_r, g_t = snr_terms(meta["H"], meta["J"], meta["v"], meta["v0"], 10.0 ** (-cfg["attenuation_db"] / 20.0))
        key = {"scenario": meta["scenario"], "snr_db": meta["snr_db"], "M": meta["M"]}
        for k, n in enumerate(n_avg):
            ab = out["abar"][:, k]
            emp = ab.var(ddof=1)
            lo, hi = bootstrap_ci(lambda x: x.var(ddof=1), ab, rng)
            row = {**key, "n": int(n), "emp_var": emp, "emp_var_lo": lo, "emp_var_hi": hi,
                   "bias": ab.mean() - A_true[k]}
            for m in NORMAL_METHODS:
                row[f"mean_v_{m}"] = out[f"v_{m}"][:, k].mean()
            if meta["scenario"] == "constant":
                row["complete_analytic"] = sum(x / d for x, d in zip(complete_shared_independent(g_r, g_t, meta["J"]),
                                                                     (1, n)))
                row["crb_analytic"] = crb(g_r, g_t, meta["J"], n)[2]
                ls, li = local_expected(g_r, g_t, meta["J"])
                row["local_analytic"] = ls + li / n
            var_rows.append(row)
            for L in levels:
                for m in METHODS:
                    flags = covered_flags(out, A_true[k], L, m)[:, k]
                    c, tot = int(flags.sum()), len(flags)
                    wlo, whi = wilson(c, tot)
                    cov_rows.append({**key, "n": int(n), "level": L, "method": m, "covered": c, "total": tot,
                                     "coverage": c / tot, "wilson_lo": wlo, "wilson_hi": whi,
                                     "mean_width_db": float(widths(out, L, m)[:, k].mean())})
            if n == n_avg[-1]:
                f = {m: covered_flags(out, A_true[k], 0.95, m)[:, k] for m in METHODS}
                for a_m, b_m in (("complete", "local"), ("complete", "bootstrap")):
                    b, c = int((f[a_m] & ~f[b_m]).sum()), int((~f[a_m] & f[b_m]).sum())
                    pair_rows.append({**key, "n": int(n), "method": a_m, "comparator": b_m,
                                      "method_only": b, "comparator_only": c, "p_exact": mcnemar_exact(b, c)})
        if meta["scenario"] == "constant":
            S_sim, D_sim = shared_cov(est)
            s_lo, s_hi = bootstrap_ci(lambda x: shared_cov(x)[0], est, rng, draws=300)
            S_c, U_c = complete_shared_independent(g_r, g_t, meta["J"])
            S_b, U_b, _ = crb(g_r, g_t, meta["J"], 1)
            S_l, U_l = local_expected(g_r, g_t, meta["J"])
            comp_rows.append({**key, "g_r": g_r, "g_t": g_t, "S_sim": S_sim, "S_sim_lo": s_lo, "S_sim_hi": s_hi,
                              "U_sim": D_sim - S_sim, "S_complete": S_c, "U_complete": U_c, "S_crb": S_b,
                              "U_crb": U_b, "S_local": S_l, "U_local": U_l,
                              "S_local_plugin_mean": float(local_shared.mean())})

    pd.DataFrame(cov_rows).to_csv(outdir / "coverage.csv", index=False)
    pd.DataFrame(var_rows).to_csv(outdir / "variance.csv", index=False)
    pd.DataFrame(comp_rows).to_csv(outdir / "components.csv", index=False)
    pd.DataFrame(pair_rows).to_csv(outdir / "paired.csv", index=False)
    (outdir / "run_info.json").write_text(json.dumps({
        "refs_per_cell": refs, "bootstrap_draws": draws, "seed": cfg["seed"],
        "wall_seconds": round(time.time() - start, 1), "J": results[0][0]["J"], "H": results[0][0]["H"]},
        indent=2))
    print("aggregates written to", outdir)


if __name__ == "__main__":
    main()
