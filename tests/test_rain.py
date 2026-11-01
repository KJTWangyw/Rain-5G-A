import pytest

from nrpilot.rain import p838_horizontal, rain_rate


@pytest.mark.parametrize("f, k, alpha", [(25, 0.1571, 0.9991), (28, 0.2051, 0.9679), (30, 0.2403, 0.9485)])
def test_regression_reproduces_table_5(f, k, alpha):
    k_calc, a_calc = p838_horizontal(f)
    assert round(k_calc, 4) == k
    assert round(a_calc, 4) == alpha


def test_rain_rate_inverts_power_law():
    k, alpha = p838_horizontal(28)
    R = rain_rate(2.0, k, alpha, 0.5)
    assert abs(k * 0.5 * R ** alpha - 2.0) < 1e-12
