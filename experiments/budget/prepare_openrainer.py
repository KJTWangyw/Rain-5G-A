#!/usr/bin/env python3
"""OpenRainER (CC BY 4.0) -> 1-min CML arrays, 15-min radar path averages and 15-min gauge data for the
calibration and evaluation seasons. Download CML.tar, AWS.tar and RADrain.tar from Zenodo record 22829808
into the configured data directory first (``--download`` does it and checks the MD5 sums).

Usage:  python experiments/budget/prepare_openrainer.py [--download] [--config experiments/budget/config.json]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from cmlbudget import openrainer as O  # noqa: E402


def md5(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "experiments" / "budget" / "config.json")
    ap.add_argument("--download", action="store_true")
    args = ap.parse_args()
    cfg = json.loads(args.config.read_text())["openrainer"]
    data = ROOT / cfg["data_dir"]
    data.mkdir(parents=True, exist_ok=True)
    for name, digest in cfg["archives"].items():
        f = data / name
        if not f.exists() and args.download:
            url = f"https://zenodo.org/records/{cfg['zenodo_record']}/files/{name}?download=1"
            print("downloading", url, flush=True)
            urllib.request.urlretrieve(url, f)
        if md5(f) != digest:
            raise SystemExit(f"MD5 mismatch for {f}")
    for name, sub in (("CML.tar", "cml"), ("AWS.tar", "aws")):
        if not (data / sub).exists():
            with tarfile.open(data / name) as tar:
                tar.extractall(data / sub)
    for gz in (data / "cml").glob("*.nc.gz"):
        nc = gz.with_suffix("")
        if not nc.exists():
            import gzip
            nc.write_bytes(gzip.decompress(gz.read_bytes()))
    t0 = time.time()
    months = O.month_list(*cfg["calibration_months"]) + O.month_list(*cfg["evaluation_months"])
    meta, arr = O.read_cml(data, months, cfg["round_db"])
    print(f"CML {arr['rsl'].shape} ({time.time() - t0:.0f} s)", flush=True)
    rt, rR = O.radar_path_average(data, meta, months, cfg["path_points"], data / "tmp")
    print(f"radar {rR.shape} ({time.time() - t0:.0f} s)", flush=True)
    gt, g, gmeta = O.read_aws(data, months)
    out = ROOT / cfg["cache"]
    np.savez(out, radar_time=rt, radar_R=rR, gauge_time=gt, gauge=g.astype(np.float32), **arr)
    meta.to_csv(out.with_name("openrainer_sublinks.csv"), index=False)
    gmeta.to_csv(out.with_name("openrainer_gauges.csv"), index=False)
    print("wrote", out, f"({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
