import numpy as np
import pytest

from nrpilot.rain import p838, p838_horizontal, p838_vertical, rain_rate

# Table 5 of ITU-R P.838-3
TABLE5 = [  # f, kH, alphaH, kV, alphaV
    (10, 0.01217, 1.2571, 0.01129, 1.2156),
    (25, 0.1571, 0.9991, 0.1533, 0.9491),
    (28, 0.2051, 0.9679, 0.1964, 0.9277),
    (30, 0.2403, 0.9485, 0.2291, 0.9129),
    (40, 0.4431, 0.8673, 0.4274, 0.8421),
]


@pytest.mark.parametrize("f, kh, ah, kv, av", TABLE5)
def test_regression_reproduces_table_5(f, kh, ah, kv, av):
    k, a = p838_horizontal(f)
    assert round(k, 4 if k > 0.1 else 5) == kh and round(a, 4) == ah
    k, a = p838_vertical(f)
    assert round(k, 4 if k > 0.1 else 5) == kv and round(a, 4) == av


def test_polarization_selector_on_arrays():
    k, a = p838(np.array([28.0, 28.0]), np.array(["Vertical", "Horizontal"]))
    assert np.allclose(k, [p838_vertical(28)[0], p838_horizontal(28)[0]])
    assert np.allclose(a, [p838_vertical(28)[1], p838_horizontal(28)[1]])
    assert p838(28.0, "V") == p838_vertical(28.0)


def test_rain_rate_inverts_power_law():
    k, alpha = p838_horizontal(28)
    R = rain_rate(2.0, k, alpha, 0.5)
    assert abs(k * 0.5 * R ** alpha - 2.0) < 1e-12
