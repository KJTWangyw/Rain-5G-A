import numpy as np
import pytest

from cmlbudget.quantization import dither_bias, dither_bias_rms, quantize


def test_quantize_grid_and_error_bound():
    x = np.linspace(-3, 3, 1001)
    y = quantize(x, 1.0, offset=0.3)
    steps = (y - 0.3) / 1.0
    assert np.allclose(steps, np.round(steps), atol=1e-9)
    assert np.max(np.abs(y - x)) <= 0.5 + 1e-12


def test_bias_without_fluctuation_is_the_sawtooth():
    c = np.array([0.1, 0.3, -0.2, 0.45])
    assert np.allclose(dither_bias(c, 0.0, 1.0, terms=4000), quantize(c, 1.0) - c, atol=2e-3)


@pytest.mark.parametrize("sigma", [0.05, 0.1, 0.2, 0.3])
def test_bias_matches_monte_carlo(sigma):
    rng = np.random.default_rng(3)
    c = np.linspace(0.0, 1.0, 7, endpoint=False) + 0.37
    mc = [np.mean(quantize(ci + sigma * rng.standard_normal(400_000), 1.0, offset=0.2) - ci) for ci in c]
    assert np.max(np.abs(np.array(mc) - dither_bias(c, sigma, 1.0, offset=0.2))) < 3e-3


def test_rms_limits_and_offset_average():
    assert abs(dither_bias_rms(0.0, 1.0, terms=4000) - 1 / np.sqrt(12)) < 1e-3
    assert dither_bias_rms(0.42, 1.0) < 0.01
    u = np.linspace(0, 1, 2000, endpoint=False)
    for s in (0.1, 0.25):
        assert abs(np.sqrt(np.mean(dither_bias(0.0, s, 1.0, offset=u) ** 2)) - dither_bias_rms(s, 1.0)) < 1e-6
    # the step only rescales
    assert abs(dither_bias_rms(0.5, 2.0) - 2 * dither_bias_rms(0.25, 1.0)) < 1e-12
