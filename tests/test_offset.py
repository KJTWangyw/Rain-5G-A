import numpy as np
import pandas as pd

from cmlbudget.offset import bin_means, fit_offset, invert, offset_curve, offset_function


def test_intercept_recovers_constant_offset():
    rng = np.random.default_rng(0)
    e = rng.uniform(0.05, 3.0, 300)
    a = 0.55 + 0.85 * e + 0.05 * rng.standard_normal(300)
    w, beta = fit_offset(a, e)
    assert abs(w - 0.55) < 0.02 and abs(beta - 0.85) < 0.02


def test_binned_curve_recovers_rain_dependent_offset():
    rng = np.random.default_rng(1)
    true_W = lambda R: 0.3 + 0.8 * np.log10(1 + 4 * R)  # noqa: E731
    bins = np.array([0, 0.5, 1, 2, 5, 10, 50])
    rows = []
    for i, L in enumerate(rng.uniform(0.2, 8.0, 200)):
        R = rng.gamma(0.6, 3.0, 4000)
        e = 0.2 * L * R ** 0.93
        a = true_W(R) + e + 0.1 * rng.standard_normal(len(R))
        for r in bin_means(a, e, R, bins, 30):
            rows.append(dict(link=i, **r))
    curve = offset_curve(pd.DataFrame(rows))
    assert len(curve) >= 5
    # W of a bin equals the mean offset over the rain rates in that bin
    for _, row in curve.iterrows():
        assert abs(row.W - true_W(row.R)) < 0.12
        assert abs(row.beta - 1.0) < 0.1


def test_offset_function_is_monotone_and_inversion_round_trips():
    W = offset_function([0.01, 0.3, 1.5, 6.0], [0.2, 0.5, 1.2, 1.1])
    R = np.logspace(-3, 2, 50)
    assert np.all(np.diff(W(R)) >= -1e-12)
    kL, alpha = 0.1, 0.93
    a = W(R) + kL * R ** alpha
    assert np.allclose(invert(a, kL, alpha, W), R, rtol=2e-2)
    assert invert(np.array([0.1, np.nan]), kL, alpha, W)[0] == 0.0
    assert np.isnan(invert(np.array([0.1, np.nan]), kL, alpha, W)[1])
