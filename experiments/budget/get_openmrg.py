#!/usr/bin/env python3
"""Download OpenMRG (CC BY-SA 4.0) from Zenodo, check its MD5 and extract the files used here.

Usage:  python experiments/budget/get_openmrg.py [--config experiments/budget/config.json]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FILES = ["cml/cml.nc", "cml/cml_metadata.csv", "radar/radar.nc", "gauges/city/CityGauges-2015JJA.csv",
         "gauges/city/CityGauges-metadata.csv", "readme.txt"]


def md5(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "experiments" / "budget" / "config.json")
    cfg = json.loads(ap.parse_args().config.read_text())
    data = ROOT / cfg["data_dir"]
    if all((data / f).exists() for f in FILES):
        print("OpenMRG already present in", data)
        return
    data.mkdir(parents=True, exist_ok=True)
    archive = data / cfg["archive"]
    if not archive.exists():
        url = f"https://zenodo.org/records/{cfg['zenodo_record']}/files/{cfg['archive']}?download=1"
        print("downloading", url)
        urllib.request.urlretrieve(url, archive)
    if md5(archive) != cfg["archive_md5"]:
        raise SystemExit(f"MD5 mismatch for {archive}")
    with zipfile.ZipFile(archive) as z:
        for f in FILES:
            z.extract(f, data)
    print("extracted", ", ".join(FILES), "to", data)


if __name__ == "__main__":
    main()
