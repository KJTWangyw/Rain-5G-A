#!/usr/bin/env python3
"""Single-core cost of each interval at the main cell (one reference reused for 32 observations).

Times the variance/interval computation only (the attenuation estimates are shared by all methods).
Usage:  python experiments/benchmark_timing.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nrpilot.bootstrap import bootstrap_means, shrunk_channel  # noqa: E402
from nrpilot.estimator import profile_white  # noqa: E402
from nrpilot.uncertainty import local_terms, local_variance, plugin_complete_variance  # noqa: E402

sys.path.insert(0, str(ROOT / "experiments"))
from run_main import nr_objects  # noqa: E402
from nrpilot.waveform import observe_slots  # noqa: E402


def main(trials=200):
    cfg = json.loads((ROOT / "experiments" / "config.json").read_text())
    nr, taps = nr_objects(cfg)
    snr, M, T = cfg["main_cell"]["snr_db"], cfg["main_cell"]["M"], cfg["n_obs"]
    v = 10.0 ** (-snr / 10.0)
    v0 = v / M
    rng = np.random.default_rng(7)
    amps = np.full(T, 10.0 ** (-cfg["attenuation_db"] / 20.0))
    cases = []
    for i in range(trials):
        r = observe_slots(nr, taps, np.ones(M), np.zeros(M), v, rng, first_slot=i * (M + T)).mean(axis=0)
        z = observe_slots(nr, taps, amps, rng.uniform(0, 2 * np.pi, T), v, rng, first_slot=i * (M + T) + M)
        cases.append((r, profile_white(r[None, :], z, v0, v)))

    def local(r, A):
        nu, C0 = local_terms(r, 10.0 ** (-A / 20.0), v0, v)
        return local_variance(nu, C0, T)

    def complete(r, A):
        return plugin_complete_variance(r, A, v0, v, T)

    def bootstrap(r, A):
        return bootstrap_means(rng, shrunk_channel(r, v0), 10.0 ** (-A / 20.0), v0, v, [T],
                               cfg["bootstrap_draws"]).var()

    out = {}
    for name, fn in (("local", local), ("complete", complete), ("bootstrap", bootstrap)):
        fn(*cases[0])
        t0 = time.perf_counter()
        for r, A in cases:
            fn(r, A)
        out[f"ms_{name}"] = 1e3 * (time.perf_counter() - t0) / trials
    out.update({"trials": trials, "draws": cfg["bootstrap_draws"], "n_obs": T})
    (ROOT / "results" / "timing.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
