import numpy as np
import pytest

from nrpilot.bootstrap import crandn
from nrpilot.estimator import profile_white
from nrpilot.uncertainty import (K_DB, complete_shared_independent, crb, local_expected, local_terms,
                                 snr_terms)


def _numeric_crb_amplitude(J, n, a, v, v0, h, phases):
    """Inverse Fisher information of a with unknown per-observation phases and unknown channel."""
    P = 1 + n + 2 * J
    rows, var = [], []
    for j in range(J):
        g = np.zeros(P, complex)
        g[1 + n + j], g[1 + n + J + j] = 1, 1j
        rows.append(g), var.append(v0)
    for t in range(n):
        e = np.exp(1j * phases[t])
        for j in range(J):
            g = np.zeros(P, complex)
            g[0], g[1 + t] = e * h[j], 1j * a * e * h[j]
            g[1 + n + j], g[1 + n + J + j] = a * e, 1j * a * e
            rows.append(g), var.append(v)
    G = np.array(rows)
    F = np.real((G.conj().T * (2 / np.array(var))) @ G)
    return np.linalg.inv(F)[0, 0]


@pytest.mark.parametrize("J, n", [(8, 3), (24, 6)])
def test_crb_matches_numeric_fisher_information(J, n):
    rng = np.random.default_rng(4)
    h = np.exp(2j * np.pi * rng.uniform(size=J)) * rng.uniform(0.6, 1.4, J)
    a, v, v0 = 10 ** -0.1, 1.0, 0.5
    num = _numeric_crb_amplitude(J, n, a, v, v0, h, rng.uniform(0, 2 * np.pi, n))
    H = np.sum(np.abs(h) ** 2)
    g_r, g_t = snr_terms(H, J, v, v0, a)
    closed_db = crb(g_r, g_t, J, n)[2]
    assert abs(closed_db / ((K_DB / a) ** 2 * num) - 1) < 1e-9


def test_complete_propagation_matches_simulation_at_low_snr():
    rng = np.random.default_rng(5)
    J, v, M, a, T, refs = 144, 1.0, 2, 10 ** -0.1, 32, 3000
    v0 = v / M
    h = np.exp(2j * np.pi * rng.uniform(size=J))
    est = np.empty((refs, T))
    for i in range(refs):
        r = h + crandn(rng, J, v0)
        z = a * np.exp(1j * rng.uniform(0, 2 * np.pi, T))[:, None] * h[None, :] + crandn(rng, (T, J), v)
        est[i] = profile_white(r[None, :], z, v0, v)
    cov = np.cov(est, rowvar=False)
    shared_sim = cov[~np.eye(T, dtype=bool)].mean()
    g_r, g_t = snr_terms(J, J, v, v0, a)
    S, U = complete_shared_independent(g_r, g_t, J)
    assert abs(shared_sim / S - 1) < 0.1
    assert abs((np.diag(cov).mean() - shared_sim) / U - 1) < 0.1


def test_local_propagation_expectation():
    rng = np.random.default_rng(6)
    J, v, M, a = 144, 1.0, 2, 10 ** -0.1
    v0 = v / M
    h = np.exp(2j * np.pi * rng.uniform(size=J))
    shared, indep = [], []
    for _ in range(400):
        r = h + crandn(rng, J, v0)
        nu, C0 = local_terms(r, np.array([a]), v0, v)
        shared.append(C0[0, 0]), indep.append(nu[0])
    g_r, g_t = snr_terms(J, J, v, v0, a)
    exp_shared, exp_indep = local_expected(g_r, g_t, J)
    assert abs(np.mean(shared) / exp_shared - 1) < 0.02
    assert abs(np.mean(indep) / exp_indep - 1) < 0.02
