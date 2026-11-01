import numpy as np
import pytest

from nrpilot.bootstrap import crandn
from nrpilot.estimator import profile_numeric, profile_white


@pytest.mark.parametrize("v, M, A", [(1.0, 2, 2.0), (0.03, 16, 2.0), (0.3, 8, 10.0), (1.0, 64, 0.0)])
def test_closed_form_root_matches_numeric_solver(v, M, A):
    rng = np.random.default_rng(2)
    J = 144
    h = np.exp(2j * np.pi * rng.uniform(size=J)) * rng.uniform(0.5, 1.5, J)
    r = h + crandn(rng, J, v / M)
    z = 10 ** (-A / 20) * np.exp(0.4j) * h[None, :] + crandn(rng, (40, J), v)
    closed = profile_white(r[None, :], z, v / M, v)
    numeric = profile_numeric(r[None, :], z, v / M, v)
    assert np.max(np.abs(closed - numeric)) < 1e-5


def test_exact_reference_reduces_to_correlation_ratio():
    rng = np.random.default_rng(3)
    h = np.exp(2j * np.pi * rng.uniform(size=144))
    z = 0.7 * h[None, :] + crandn(rng, (5, 144), 0.1)
    est = profile_white(h[None, :], z, 0.0, 0.1)
    ratio = np.abs((np.conj(h) * z).sum(-1)) / (np.abs(h) ** 2).sum()
    np.testing.assert_allclose(est, -20 * np.log10(ratio), atol=1e-12)
