#!/usr/bin/env python3
"""Audit archived tables and recompute the physical-link bootstrap.

Run from the repository root with ``python experiments/budget/audit_summary.py``.
The output is ``results/budget/audit_summary.json``. A numerical mismatch exits
with status 1. This audit uses the archived CSV tables, without loading raw
OpenMRG observations or events.csv.gz. Its scope is recorded in the output.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "budget"
sys.path.insert(0, str(ROOT / "src"))
from cmlbudget.offset import cluster_bootstrap, fit_offset, offset_function  # noqa: E402

# These intervals are the ones used by analyze_openmrg.py, including its labels.
BANDS = {"26-30": (24.0, 30.5), "32-34": (30.5, 35.0), "37-40": (35.0, 40.5)}
SOURCE_COMMIT = "3174afd04993fc829858cf5fca0228dd6a351ed1"


def check_equal(checks, name, computed, recorded):
    checks.append({
        "name": name,
        "status": "pass" if computed == recorded else "mismatch",
        "computed": computed,
        "recorded": recorded,
    })


def check_close(checks, name, computed, recorded, tolerance=1e-12):
    actual = np.asarray(computed, dtype=float)
    expected = np.asarray(recorded, dtype=float)
    passed = actual.shape == expected.shape and np.allclose(actual, expected, rtol=0, atol=tolerance)
    checks.append({
        "name": name,
        "status": "pass" if passed else "mismatch",
        "computed": actual.tolist(),
        "recorded": expected.tolist(),
        "absolute_tolerance": tolerance,
    })


def input_fingerprints(paths):
    records = {}
    for path in paths:
        data = path.read_bytes()
        blob = b"blob " + str(len(data)).encode("ascii") + b"\0" + data
        records[str(path.relative_to(ROOT))] = {
            "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
            "computed_local_git_blob_sha1": hashlib.sha1(blob).hexdigest(),
        }
    return records


def main():
    config_path = ROOT / "experiments" / "budget" / "config.json"
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    recorded = json.loads((RES / "summary.json").read_text(encoding="utf-8"))
    offset = json.loads((RES / "offset.json").read_text(encoding="utf-8"))
    run_info = json.loads((RES / "run_info.json").read_text(encoding="utf-8"))
    links = pd.read_csv(RES / "links.csv")
    evaluation = pd.read_csv(RES / "evaluation.csv")
    curves = pd.read_csv(RES / "offset_curve.csv")
    lags = pd.read_csv(RES / "dry_lag.csv")
    pilot = pd.read_csv(RES / "pilot_reference.csv")
    checks = []

    use = links.loc[links.wet_min >= cfg["min_wet_minutes"]].copy()
    cal = links.loc[links.wet_min_cal >= cfg["min_wet_minutes"]].copy()
    edges = cfg["length_classes_km"]
    labels = [f"{left:g}–{right:g}" for left, right in zip(edges[:-1], edges[1:])]
    use["cls"] = pd.cut(use.length_km, edges, labels=labels, right=True)
    counts = use.cls.value_counts().sort_index().astype(int).to_dict()
    selection = {
        "configured_frequency_ghz": cfg["frequency_ghz"],
        "selected_sublinks": len(links),
        "selected_physical_paths": int(links.Link.nunique()),
        "usable_sublinks": len(use),
        "usable_physical_paths": int(use.Link.nunique()),
        "observed_frequency_ghz": [float(use.freq_ghz.min()), float(use.freq_ghz.max())],
        "observed_length_km": [float(use.length_km.min()), float(use.length_km.max())],
        "polarization_counts": use.pol.value_counts().astype(int).to_dict(),
        "length_class_counts": counts,
        "minimum_usable_wet_minutes": cfg["min_wet_minutes"],
        "sub_link_wet_segments_all_selected": int(links.events.sum()),
        "sub_link_wet_segments_usable": int(use.events.sum()),
        "independent_weather_event_count": None,
        "calibration_sublinks": len(cal),
        "calibration_ids_equal_full_period_ids": set(cal.j) == set(use.j),
        "excluded": links.loc[~links.j.isin(use.j), ["j", "Sublink", "Link", "wet_min", "events"]].to_dict("records"),
    }
    values = {
        "n_sel": len(links), "n_use": len(use), "n_phys": int(use.Link.nunique()),
        "n_events": f"{int(links.events.sum()):,}".replace(",", " "),
        "L_min": f"{use.length_km.min():.2f}", "L_max": f"{use.length_km.max():.1f}",
        "L_med": f"{use.length_km.median():.2f}",
        "f_min": f"{use.freq_ghz.min():.1f}", "f_max": f"{use.freq_ghz.max():.1f}",
        "n_short": int((use.length_km <= 1).sum()),
        "n_vertical": int((use.pol == "Vertical").sum()), "class_counts": counts,
    }
    for name, value in values.items():
        check_equal(checks, f"summary.{name}", value, recorded[name])
    for name, value in {"links": len(links), "usable_links": len(use), "events": int(links.events.sum())}.items():
        check_equal(checks, f"run_info.{name}", value, run_info[name])

    def fit(frame):
        return fit_offset(frame.att_mean.to_numpy(), frame.exp_mean.to_numpy())

    def fit_without_radar_amplitude(frame):
        return fit_offset(frame.att_mean.to_numpy(), frame.kL.to_numpy())

    draws, seed = cfg["bootstrap"]["draws"], cfg["bootstrap"]["seed"]
    w, beta = fit(use)
    bootstrap = cluster_bootstrap(use, fit, draws, seed)
    intervals = np.percentile(bootstrap, [2.5, 97.5], axis=0).T
    w_rf, slope_rf = fit_without_radar_amplitude(use)
    bootstrap_rf = cluster_bootstrap(use, fit_without_radar_amplitude, draws, seed + 1)
    interval_rf = np.percentile(bootstrap_rf[:, 0], [2.5, 97.5])
    residual_sd = float(np.std(use.att_mean - w - beta * use.exp_mean))
    for name, value in {
        "W": w, "beta": beta, "W_ci": intervals[0], "beta_ci": intervals[1],
        "W_radar_free": w_rf, "slope_radar_free": slope_rf,
        "W_radar_free_ci": interval_rf, "resid_sd": residual_sd,
    }.items():
        check_close(checks, f"offset.{name}", value, offset[name])
    calibration_offsets = {}
    for variant, column in {"native": "att_mean_cal", "q1": "att_mean_cal_q1", "q2": "att_mean_cal_q2"}.items():
        attenuation = cal[column].to_numpy()
        if variant != "native":
            # analyze_openmrg.py computes these means from its float32 rain_att
            # arrays. Restore that dtype after the CSV round trip, before the
            # estimator's float64 conversion, to recover the original values.
            attenuation = attenuation.astype(np.float32)
        calibration_offsets[variant] = fit_offset(attenuation, cal.exp_mean_cal.to_numpy())[0]
        check_close(checks, f"offset.W_cal.{variant}", calibration_offsets[variant], offset["W_cal"][variant])
    offset_result = {
        "W": w, "beta": beta, "W_ci95": intervals[0].tolist(), "beta_ci95": intervals[1].tolist(),
        "radar_amplitude_free_W": w_rf, "radar_amplitude_free_slope": slope_rf,
        "radar_amplitude_free_W_ci95": interval_rf.tolist(), "residual_sd": residual_sd,
        "calibration_offsets": calibration_offsets,
        "bootstrap_cluster": "physical Link, with both directions kept together",
        "bootstrap_draws": draws, "bootstrap_seed": seed, "weather_event_resampling": False,
    }
    for name, value in {"W": w, "beta": beta, "W_rf": w_rf, "W_resid_sd": residual_sd, "W_cal": calibration_offsets["native"]}.items():
        check_equal(checks, f"summary.{name}", f"{value:.2f}", recorded[name])
    for name, value in {"W_ci": intervals[0], "beta_ci": intervals[1], "W_rf_ci": interval_rf}.items():
        check_equal(checks, f"summary.{name}", f"{value[0]:.2f}–{value[1]:.2f}", recorded[name])

    bands = {}
    for index, (name, (lower, upper)) in enumerate(BANDS.items()):
        group = use.loc[(use.freq_ghz >= lower) & (use.freq_ghz < upper)]
        band_w, band_beta = fit(group)
        band_bootstrap = cluster_bootstrap(group, fit, draws, seed + 2 + index)
        band_interval = np.percentile(band_bootstrap[:, 0], [2.5, 97.5])
        bands[name] = {
            "configured_frequency_interval": [lower, upper], "upper_endpoint_included": False,
            "observed_frequency_ghz": [float(group.freq_ghz.min()), float(group.freq_ghz.max())],
            "sublinks": len(group), "physical_paths": int(group.Link.nunique()),
            "W": band_w, "beta": band_beta, "W_ci95": band_interval.tolist(),
        }
        check_equal(checks, f"offset.bands.{name}.n", len(group), offset["bands"][name]["n"])
        for field, value in {"W": band_w, "beta": band_beta, "W_ci": band_interval}.items():
            check_close(checks, f"offset.bands.{name}.{field}", value, offset["bands"][name][field])
        suffix = name.replace("-", "_")
        check_equal(checks, f"summary.W_{suffix}", f"{band_w:.2f}", recorded[f"W_{suffix}"])
        check_equal(checks, f"summary.W_{suffix}_ci", f"{band_interval[0]:.2f}–{band_interval[1]:.2f}", recorded[f"W_{suffix}_ci"])
        check_equal(checks, f"summary.n_{suffix}", len(group), recorded[f"n_{suffix}"])

    use["share"] = w / use.att_mean
    shares = []
    for label, group in use.groupby("cls", observed=True):
        shares.append({
            "class": str(label), "sublinks": len(group), "physical_paths": int(group.Link.nunique()),
            "median_of_W_over_att_mean": float(group.share.median()),
            "ratio_of_sums_for_comparison": float(w * len(group) / group.att_mean.sum()),
        })
    for position, name in ((0, "share_05"), (1, "share_1"), (4, "share_4")):
        check_equal(checks, f"summary.{name}", f"{100 * shares[position]['median_of_W_over_att_mean']:.0f}", recorded[name])

    full_curve = curves.loc[(curves.variant == "native") & (curves.period == "full")].sort_values("bin")
    curve_columns = ["bin", "R", "W", "W_lo", "W_hi", "beta", "resid_sd", "W_radar_free", "links"]
    check_equal(checks, "summary.curve", full_curve[curve_columns].round(4).to_dict("records"), recorded["curve"])
    W_function = offset_function(full_curve.R.to_numpy(), full_curve.W.to_numpy())
    for rain in (1.0, 5.0):
        check_equal(checks, f"summary.W_at_{rain:g}", f"{float(W_function(rain)):.2f}", recorded[f"W_at_{rain:g}"])

    dry_results = []
    for lag, group in lags.groupby("lag_h"):
        dry_results.append({
            "lag_h": int(lag), "median_sd_h": float(group.sd_h.median()),
            "median_rms_min": float(group.rms_min.median()), "median_mean_h": float(group.mean_h.median()),
            "sub_link_hours": int(group.n.sum()), "sublinks_with_finite_sd": int(group.sd_h.notna().sum()),
        })
    check_equal(checks, "summary.drift_1h", f"{dry_results[0]['median_sd_h']:.2f}", recorded["drift_1h"])
    check_equal(checks, "summary.drift_6h", f"{dry_results[-1]['median_sd_h']:.2f}", recorded["drift_6h"])
    check_equal(checks, "summary.n_dry_hours", f"{dry_results[0]['sub_link_hours']:,}".replace(",", " "), recorded["n_dry_hours"])

    evaluated = evaluation.merge(use[["j", "Link", "cls"]], on="j", validate="many_to_one")
    evaluated["ratio"] = evaluated.acc_link / evaluated.acc_radar
    evaluation_rows = []
    for (label, variant, correction), group in evaluated.groupby(["cls", "variant", "correction"], observed=True):
        evaluation_rows.append({
            "cls": str(label), "variant": variant, "correction": correction,
            "ratio": float(group.ratio.median()), "rmse": float(group.rmse.median()),
            "r": float(group.r.median()), "pod": float(group.pod.median()), "n": len(group),
        })
    check_equal(checks, "summary.evaluation", pd.DataFrame(evaluation_rows).round(4).to_dict("records"), recorded["evaluation"])
    evaluation_table = pd.DataFrame(evaluation_rows)
    scalar_evaluation = {}
    for position, tag in ((0, "05"), (1, "1"), (2, "2"), (4, "4")):
        group = evaluation_table.loc[evaluation_table.cls == labels[position]].set_index(["variant", "correction"])
        for correction in ("none", "constant", "rain"):
            for metric, decimals in (("ratio", 1), ("rmse", 2)):
                key = f"ev_{metric}_{correction}_{tag}"
                scalar_evaluation[key] = float(group.loc[("native", correction), metric])
                check_equal(checks, f"summary.{key}", f"{scalar_evaluation[key]:.{decimals}f}", recorded[key])
        reduction = 1 - group.loc[("native", "rain"), "rmse"] / group.loc[("native", "none"), "rmse"]
        scalar_evaluation[f"reduction_of_group_medians_{tag}"] = float(reduction)
        check_equal(checks, f"summary.ev_red_{tag}", f"{100 * reduction:.0f}", recorded[f"ev_red_{tag}"])
    evaluation_result = {
        "calibration_months": cfg["calibration_months"], "evaluation_months": cfg["evaluation_months"],
        "minimum_valid_minutes_per_hour": cfg["hour_min_minutes"],
        "raw_rows": len(evaluation), "usable_rows": len(evaluated),
        "hours_counts_equal_across_variants": bool(evaluated.groupby("j").hours.nunique().eq(1).all()),
        "radar_accumulations_equal_across_variants": bool(evaluated.groupby("j").acc_radar.nunique().eq(1).all()),
        "rows": evaluation_rows, "selected_scalar_values": scalar_evaluation,
    }

    k_db = 20.0 / np.log(10.0)
    shared_bound = k_db**2 / (2 * pilot.J * pilot.g_r)
    shared_complete = shared_bound + k_db**2 / (pilot.J * (pilot.g_r + pilot.g_t)**2)
    check_close(checks, "pilot.S_crb", shared_bound, pilot.S_crb)
    check_close(checks, "pilot.S_complete", shared_complete, pilot.S_complete)
    in_ci = int(((pilot.S_lo <= shared_complete) & (shared_complete <= pilot.S_hi)).sum())
    relative = pilot.S_sim / shared_complete - 1
    check_equal(checks, "summary.pilot_cells", len(pilot), recorded["pilot_cells"])
    check_equal(checks, "summary.pilot_in_ci", in_ci, recorded["pilot_in_ci"])
    check_equal(checks, "summary.pilot_rel_range", f"{100 * relative.min():+.0f}%至{100 * relative.max():+.0f}%".replace("-", "−"), recorded["pilot_rel_range"])
    J = int(pilot.J.iloc[0])
    reference_sd = float(k_db / np.sqrt(2 * J * 80 * 10**1.5))
    check_equal(checks, "summary.pilot_sd_15_10ms", f"{reference_sd:.4f}", recorded["pilot_sd_15_10ms"])

    inputs = [config_path] + [RES / name for name in (
        "summary.json", "offset.json", "run_info.json", "links.csv", "evaluation.csv",
        "offset_curve.csv", "dry_lag.csv", "pilot_reference.csv",
    )]
    inputs.extend([ROOT / "src/cmlbudget/offset.py", Path(__file__).resolve()])
    mismatches = [check["name"] for check in checks if check["status"] == "mismatch"]
    report = {
        "status": "mismatch" if mismatches else "pass",
        "source_commit_for_archived_inputs": SOURCE_COMMIT,
        "scope": "Archived-table aggregation and physical-link cluster-bootstrap recomputation",
        "checks_run": len(checks), "mismatch_count": len(mismatches), "mismatches": mismatches,
        "versions": {"python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__},
        "definitions": {
            "constant_offset": "Unweighted OLS across one full-period mean attenuation and radar-expected attenuation per sublink",
            "radar_amplitude_free": "OLS on kL, with wet periods and clear-sky reference selection retaining radar input",
            "short_path_share": "Median across sublinks of W / att_mean, using right-closed length intervals",
            "event_count": "Sum of retained wet segments counted separately for each sublink, with independent weather-process count unavailable",
            "evaluation_rmse": "Median of per-sublink hourly RMSE within each length class",
            "evaluation_cumulative_ratio": "Median of per-sublink accumulated-link / accumulated-radar ratios",
            "evaluation_label_usage": "Radar-centered wet masks and radar-confirmed dry outputs of zero are also used in August",
            "support_checks": "Counts and radar accumulation totals are compared across methods, while timestamp checks require hourly inputs",
            "calibration_csv_dtype": "The q1 and q2 calibration attenuation means are restored to their original float32 dtype before fitting",
            "budget_rain_rate_conversion": "Standalone equivalent rain rate (attenuation / kL) ** (1 / alpha)",
        },
        "outside_this_audit": [
            "Raw observations, minute aggregation, radar sampling, wet/dry masks, and individual reference windows",
            "Per-event quantization RMS and dither predictions stored in events.csv.gz",
            "Regression and bootstrap of rain-rate-bin curves from bins.csv",
            "Sensitivity reruns from the original minute cache",
            "NR waveform generation and estimation samples underlying pilot_reference.csv",
        ],
        "selection": selection, "offset": offset_result, "frequency_bands": bands,
        "offset_share_by_length": shares, "dry_reference_statistics": dry_results,
        "evaluation": evaluation_result, "pilot_reference_sd_15db_10ms": reference_sd,
        "input_fingerprints": input_fingerprints(inputs), "checks": checks,
    }
    output = RES / "audit_summary.json"
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"], "checks_run": len(checks), "mismatches": mismatches,
        "output": str(output.relative_to(ROOT)), "scope": report["scope"],
    }, ensure_ascii=False))
    return 1 if mismatches else 0


if __name__ == "__main__":
    raise SystemExit(main())
