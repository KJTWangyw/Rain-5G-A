#!/usr/bin/env python3
"""Shared clear-sky reference error of DM-RS pilot-based attenuation through the NR OFDM chain.

For every (SNR, reference length M) cell, ``refs`` references are built from M coherent clear-sky slots and each
is reused for ``n_obs`` rain slots with random common phases. The shared component is the covariance of two
estimates that use the same reference, estimated as Var(mean over a reference) - mean within-reference
variance / n_obs, and compared with the complete closed form and the Cramer-Rao bound.

Usage:  python experiments/budget/pilot_reference.py [--refs 400]
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from nrpilot.estimator import profile_white  # noqa: E402
from nrpilot.uncertainty import complete_shared_independent, crb  # noqa: E402
from nrpilot.waveform import Multipath, NRConfig, observe_slots, pilot_channel  # noqa: E402

SETTINGS = {"n_rb": 66, "snr_db": [5.0, 15.0, 25.0], "M": [1, 8, 80], "attenuation_db": 0.5, "n_obs": 128,
            "seed": 20261107, "slot_ms": 0.125}


def nr_objects(n_rb):
    base = json.loads((ROOT / "experiments" / "config.json").read_text())
    nr = NRConfig(n_fft=base["nr"]["n_fft"], n_cp=base["nr"]["n_cp"], scs_hz=base["nr"]["scs_hz"], n_rb=n_rb,
                  dmrs_symbols=tuple(base["nr"]["dmrs_symbols"]), n_id=base["nr"]["n_id"])
    mp = base["multipath"]
    gains = tuple(m * np.exp(1j * p) for m, p in zip(mp["magnitudes"], mp["phases_rad"]))
    return nr, Multipath(delays=tuple(mp["delays"]), gains=gains).taps(nr)


def run_cell(task):
    nr, taps = nr_objects(SETTINGS["n_rb"])
    J, M, n_obs, refs = nr.n_pilots, task["M"], SETTINGS["n_obs"], task["refs"]
    v = 10.0 ** (-task["snr_db"] / 10.0)
    amps = np.full(n_obs, 10.0 ** (-SETTINGS["attenuation_db"] / 20.0))
    rng = np.random.default_rng([SETTINGS["seed"], task["cell"]])
    means, within = np.empty(refs), np.empty(refs)
    for i in range(refs):
        first = i * (M + n_obs)
        r = observe_slots(nr, taps, np.ones(M), np.zeros(M), v, rng, first_slot=first).mean(axis=0)
        z = observe_slots(nr, taps, amps, rng.uniform(0, 2 * np.pi, n_obs), v, rng, first_slot=first + M)
        A_hat = profile_white(r[None, :], z, v / M, v)
        means[i], within[i] = A_hat.mean(), A_hat.var(ddof=1)
    H = float(np.sum(np.abs(pilot_channel(nr, taps)) ** 2))
    g_r, g_t = H / (J * v / M), amps[0] ** 2 * H / (J * v)
    S_c, U_c = complete_shared_independent(g_r, g_t, J)
    S_crb, U_crb, _ = crb(g_r, g_t, J, 1)
    shared = lambda idx: np.var(means[idx], ddof=1) - np.mean(within[idx]) / n_obs  # noqa: E731
    b = np.random.default_rng([SETTINGS["seed"], 99, task["cell"]])
    boot = [shared(b.integers(0, refs, refs)) for _ in range(1000)]
    return dict(snr_db=task["snr_db"], M=M, reference_ms=M * SETTINGS["slot_ms"], J=J, refs=refs, n_obs=n_obs,
                g_r=g_r, g_t=g_t, S_sim=shared(np.arange(refs)), S_lo=np.percentile(boot, 2.5),
                S_hi=np.percentile(boot, 97.5), U_sim=np.mean(within), S_complete=S_c, U_complete=U_c,
                S_crb=S_crb, U_crb=U_crb, bias=np.mean(means) - SETTINGS["attenuation_db"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refs", type=int, default=400)
    ap.add_argument("--processes", type=int, default=4)
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "budget" / "pilot_reference.csv")
    args = ap.parse_args()
    tasks = [dict(cell=c, snr_db=s, M=m, refs=args.refs)
             for c, (s, m) in enumerate((s, m) for s in SETTINGS["snr_db"] for m in SETTINGS["M"])]
    t0 = time.time()
    with Pool(args.processes) as pool:
        rows = pool.map(run_cell, tasks)
    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    (args.out.with_suffix(".json")).write_text(json.dumps({**SETTINGS, "refs": args.refs,
                                                           "runtime_s": round(time.time() - t0, 1)}, indent=1))
    with pd.option_context("display.width", 200):
        print(df[["snr_db", "M", "S_sim", "S_lo", "S_hi", "S_complete", "S_crb", "U_sim", "U_complete", "bias"]])
    print(f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
