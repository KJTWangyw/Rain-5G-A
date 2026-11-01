# Rain-5G-A

这是基于5G-A的降雨监测算法 — rain attenuation estimation from 5G NR pilots.

This repository contains the code, results and manuscript for

> **Shared clear-sky reference errors in 5G NR pilot-based rain attenuation estimation: information bound and closed-form uncertainty propagation**

A clear-sky reference estimated from a finite number of NR slots is reused for every rain observation, so its
error is common to all attenuation estimates and is not reduced by averaging. The code derives and checks

- the Cramér–Rao bound of the averaged attenuation with a shared reference and the shared fraction `n a²/(M + n a²)`;
- a complete closed-form variance of the profile amplitude estimator, compared with the fixed-weight
  linearisation, a wrongly independent model and a parametric bootstrap;
- interval coverage through an NR DM-RS OFDM chain (TS 38.211, numerology 3, 28 GHz) over 16 SNR/reference conditions.

## Layout

| Path | Content |
|---|---|
| `src/nrpilot/waveform.py` | NR DM-RS (type 1, symbols 2 and 11) through IFFT, cyclic prefix, multipath, noise, FFT and LS extraction |
| `src/nrpilot/estimator.py` | Profile amplitude estimator (closed-form root and numeric solver) |
| `src/nrpilot/uncertainty.py` | Information bound, fixed-weight linearisation, complete closed-form propagation |
| `src/nrpilot/bootstrap.py` | Parametric bootstrap with a shared reference draw |
| `src/nrpilot/rain.py` | ITU-R P.838-3 coefficients (horizontal polarisation) |
| `experiments/` | Configuration, Monte Carlo study, timing benchmark, summary and figures |
| `results/` | Aggregated results (`coverage.csv`, `variance.csv`, `components.csv`, `paired.csv`, `summary.json`, ...) |
| `paper/` | Manuscript source (`content.md`; Chinese `content_zh.md`), Word build script, figures, `Manuscript.docx`/`.pdf` and `Manuscript_zh.docx`/`.pdf` |
| `tests/` | Unit tests (pilot extraction, estimator, bound, propagation, P.838-3) |

## Reproduce

```bash
pip install -r requirements.txt
python -m pytest -q                      # about 5 s
python experiments/run_main.py           # 17 conditions x 2000 references; about 31 min on 4 cores
python experiments/benchmark_timing.py   # single-core cost of each interval
python experiments/summarize.py          # results/summary.json (every number quoted in the paper)
python experiments/make_figures.py       # paper/figures/*.png and *.pdf
python paper/build_manuscript.py         # paper/Manuscript.docx and Manuscript.pdf
python experiments/make_figures.py --lang zh && python paper/build_manuscript.py --lang zh   # Chinese version
```

`run_main.py` writes per-trial arrays to `results/raw/` (not versioned) and the aggregates to `results/`.
Seeds are fixed in `experiments/config.json`; a rerun reproduces the aggregates exactly.

Building the manuscript needs XeLaTeX with `unicode-math` and the TeX Gyre Termes Math font (equations),
LibreOffice with the Math component (PDF conversion), and Times New Roman and Arial for an exact match of the
journal template. `paper/template/reference.docx` holds only the template's styles and page setup.
The Chinese version uses the template's 宋体/黑体 slots; on Linux they are rendered with Noto Serif/Sans CJK SC
(map them with a fontconfig alias).
