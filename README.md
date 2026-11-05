# Rain-5G-A

这是基于5G-A的降雨监测算法 — rainfall sensing on short millimetre-wave links for 5G-Advanced.

This branch holds the code, results and Chinese manuscript of

> **面向5G-A的短距毫米波链路降雨感知误差预算：湿天线偏移的识别与校正**
> (Error budget of rainfall sensing on short millimetre-wave links for 5G-A: identifying and correcting the
> wet-antenna offset)

On a short link the rain-induced attenuation is a few tenths of a decibel, so every attenuation error that does
not scale with the path is amplified as 1/(kL) in rain rate. The study separates four such errors — the
wet-period length-independent offset, drift of the clear-sky reference, quantization of power reports and
pilot estimation noise — and quantifies them on 646 OpenMRG sub-links at 28–39 GHz (0.14–11 km) with radar and gauges:

- the length-independent offset grows with rain rate from about 0.2 dB to 2.4 dB, decays within about half an hour
  after the rain on the path has stopped, identically for short and long links (wet antennas), and dominates links
  below 1 km;
- reference drift (≈0.2–0.3 dB) and 1 dB report quantization (≈0.14 dB per event, mostly in the shared clear-sky
  reference, predicted per event by a closed form) come next;
- the shared reference error of DM-RS pilot estimation is ≈0.004 dB with a 10 ms clear-sky reference (NR chain
  simulation, Cramér–Rao bound);
- with a rain-dependent wet-antenna model (KR-alt, literature parameters) links below 1 km reach an hourly RMSE of
  0.22 mm/h against city gauges in August, better than the weather radar (0.33 mm/h); a constant offset does not.

A second region repeats the checks: OpenRainER (Emilia-Romagna, Italy; 285 sub-links at 24.6–25.6 GHz, 0.16–25 km,
levels logged at 1 dB; May–September 2021 for calibration, 2022 for evaluation). The cross-link regression gives the
same offset (0.54 dB), the length-independent residual again decays after the rain stops, and KR-alt brings links
below 2 km from 1.97 to 0.34 mm/h against gauges, below the radar (0.77 mm/h).

## Layout

| Path | Content |
|---|---|
| `src/cmlbudget/openmrg.py` | OpenMRG reader: 1-min signal levels, re-quantized reports, radar path averages |
| `src/cmlbudget/openrainer.py` | OpenRainER reader: 1-min levels (rounded to 1 dB), 15-min radar path averages, gauges |
| `src/cmlbudget/retrieval.py` | Wet/dry masks (radar or rolling standard deviation), event references, rain rate |
| `src/cmlbudget/quantization.py` | Rounding quantizer, mean rounding error under fluctuation and its RMS |
| `src/cmlbudget/offset.py` | Length-independent offset by regression across links, rain-dependent W(R), inversion |
| `src/nrpilot/` | NR DM-RS chain, profile estimator, bound and closed-form propagation, ITU-R P.838-3 (H and V) |
| `experiments/budget/` | Configuration, data download, preparation, analysis, pilot simulation, summary, figures |
| `results/budget/` | Aggregated results and `summary.json` (every number quoted in the manuscript); `openrainer/` for the second region |
| `paper/` | Manuscript source `content_zh.md`, Word build script, figures, `Manuscript_zh.docx` and `.pdf` |
| `tests/` | Unit tests |

The earlier pilot coverage study (pull request #1) remains in `experiments/run_main.py`, `summarize.py`,
`benchmark_timing.py`, `make_figures.py` and `results/*.csv`; it is not used by the current manuscript.

## Reproduce

```bash
pip install -r requirements.txt
python -m pytest -q
python experiments/budget/get_openmrg.py        # downloads OpenMRG (318 MB) into data/OpenMRG, checks MD5
python experiments/budget/prepare_openmrg.py    # 1-min cache (about 2 GB, 5 min)
python experiments/budget/analyze_openmrg.py    # results/budget/*.csv, offset.json (about 4 min)
python experiments/budget/evaluate_corrections.py  # correction models vs radar and gauges, drying (about 2 min)
python experiments/budget/prepare_openrainer.py --download  # OpenRainER CML, AWS and RADrain (3 GB), MD5, cache
python experiments/budget/analyze_openrainer.py # second region (about 3 min)
python experiments/budget/pilot_reference.py    # NR DM-RS simulation (about 5 min on 3 cores)
python experiments/budget/summarize.py
python experiments/budget/make_figures.py       # paper/figures_zh, results/budget/figures_info.json
python experiments/budget/summarize.py          # adds the Fig. 1 example to summary.json
python paper/build_manuscript.py --lang zh      # paper/Manuscript_zh.docx and .pdf
```

Seeds and all processing settings are in `experiments/budget/config.json`.

Building the manuscript needs XeLaTeX with `unicode-math` and TeX Gyre Termes Math (equations), LibreOffice with
the Math component (PDF conversion), Times New Roman and Arial, and the template's 宋体/黑体 slots mapped to
Noto Serif/Sans CJK SC with a fontconfig alias. `paper/template/reference.docx` holds only the template's styles and
page setup.

## Data

OpenMRG (Andersson et al., 2022, Earth Syst. Sci. Data 14, 5411–5426; CML, radar and city gauges) is published under
CC BY-SA 4.0 at https://doi.org/10.5281/zenodo.7107689. OpenRainER (Covi and Roversi, version 2.0.3; CML, radar,
weather stations) is published under CC BY 4.0 at https://doi.org/10.5281/zenodo.10593848. The raw data are not
stored here; `results/budget/` contains derived aggregates only.
