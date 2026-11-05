# Rain-5G-A

这是基于5G-A的降雨监测算法 — rainfall sensing on short millimetre-wave links for 5G-Advanced.

This branch holds the code, results and Chinese manuscript of

> **面向5G-A的短距毫米波链路降雨感知误差预算：长度无关偏移的识别与标定**
> (Error budget of rainfall sensing on short millimetre-wave links for 5G-A: identifying and calibrating the
> length-independent offset)

On a short link the rain-induced attenuation is a few tenths of a decibel, so every attenuation error that does
not scale with the path is amplified as 1/(kL) in rain rate. The study separates four such errors — the
wet-period length-independent offset, drift of the clear-sky reference, quantization of power reports and
pilot estimation noise — and quantifies them on 646 OpenMRG sub-links at 26–40 GHz (0.14–11 km) with radar:

- the length-independent offset grows with rain rate from about 0.2 dB to 2.4 dB and dominates links below 1 km;
- reference drift (≈0.2–0.3 dB) and 1 dB report quantization (≈0.14 dB per event, mostly in the shared clear-sky
  reference, predicted per event by a closed form) come next;
- the shared reference error of DM-RS pilot estimation is ≈0.004 dB with a 10 ms clear-sky reference (NR chain
  simulation, Cramér–Rao bound);
- a rain-dependent offset identified from link-length diversity, calibrated in June–July, lowers the August hourly
  rain-rate RMSE of links below 0.5 km by 70%.

## Layout

| Path | Content |
|---|---|
| `src/cmlbudget/openmrg.py` | OpenMRG reader: 1-min signal levels, re-quantized reports, radar path averages |
| `src/cmlbudget/retrieval.py` | Wet/dry mask, event references, rain-period attenuation, rain rate |
| `src/cmlbudget/quantization.py` | Rounding quantizer, mean rounding error under fluctuation and its RMS |
| `src/cmlbudget/offset.py` | Length-independent offset by regression across links, rain-dependent W(R), inversion |
| `src/nrpilot/` | NR DM-RS chain, profile estimator, bound and closed-form propagation, ITU-R P.838-3 (H and V) |
| `experiments/budget/` | Configuration, data download, preparation, analysis, pilot simulation, summary, figures |
| `results/budget/` | Aggregated results and `summary.json` (every number quoted in the manuscript) |
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

OpenMRG (Andersson et al., 2022, Earth Syst. Sci. Data 14, 5411–5426) is published under CC BY-SA 4.0 at
https://doi.org/10.5281/zenodo.7107689. The raw data are not stored here; `results/budget/` contains derived
aggregates only.
