#!/usr/bin/env python3
"""OpenMRG 10-s signal levels -> 1-min arrays for the 24-40 GHz sub-links, plus radar path averages.

Writes the cache named in the configuration (about 2 GB, not versioned). Run get_openmrg.py first.
Usage:  python experiments/budget/prepare_openmrg.py [--config experiments/budget/config.json]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from cmlbudget import openmrg  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "experiments" / "budget" / "config.json")
    cfg = json.loads(ap.parse_args().config.read_text())
    data = ROOT / cfg["data_dir"]
    t0 = time.time()
    meta = openmrg.select_sublinks(openmrg.read_metadata(data), *cfg["frequency_ghz"])
    rng = np.random.default_rng(cfg["quantization_seed"])
    steps = cfg["quantization_steps_db"]
    offsets = {s: rng.uniform(0.0, s, len(meta)) for s in steps}
    arrays = openmrg.minute_aggregates(data, meta.Sublink.values, steps, offsets)
    print(f"signal levels done ({time.time() - t0:.0f} s)", flush=True)
    r = cfg["radar"]
    rt, rR = openmrg.radar_path_average(data, meta, r["zr_a"], r["zr_b"], r["min_dbz"], r["path_points"])
    print(f"radar done ({time.time() - t0:.0f} s)", flush=True)
    out = ROOT / cfg["cache"]
    np.savez(out, sublink=meta.Sublink.values, radar_time=rt, radar_R=rR,
             **{f"offset_q{s:g}": offsets[s] for s in steps}, **arrays)
    print("wrote", out, f"({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
