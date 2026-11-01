"""Variance of the averaged attenuation estimate when one clear-sky reference is reused.

Notation (white noise): J pilots, H = sum_j |h_j|^2, reference noise variance v_0 = v / M, current noise
variance v, per-pilot reference SNR g_r = H / (J v_0) and per-pilot current SNR g_t = a^2 H / (J v),
K = 20 / ln 10. For the mean of n estimates A_bar_n:

  information bound    Var >= K^2/(2 J) [1/g_r + 1/(n g_t)]
  complete propagation Var ~= S + U/n,
                       S = K^2/(2 J g_r) + K^2/(J (g_r+g_t)^2)
                       U = K^2/(2 J g_t) + K^2/(2 J g_t g_r) ((g_r-g_t)/(g_r+g_t))^2 + K^2/(J (g_r+g_t)^2)
  fixed-weight local   linearisation with weights and phase frozen at the estimate; its expectation equals
                       the bound divided by (1 + 1/g_r).
"""
from __future__ import annotations

import numpy as np

K_DB = 20.0 / np.log(10.0)


def snr_terms(H, J, v, v0, a):
    """Per-pilot reference and current SNRs (g_r, g_t)."""
    return H / (J * v0), np.asarray(a) ** 2 * H / (J * v)


def crb(g_r, g_t, J, n):
    """Cramer-Rao bound of the averaged attenuation (dB^2): (shared, independent per observation, total)."""
    shared = K_DB ** 2 / (2 * J * g_r)
    indep = K_DB ** 2 / (2 * J * g_t)
    return shared, indep, shared + indep / n


def complete_components(g_r, g_t, J):
    """Terms of the complete propagation (dB^2). Shared: S = s_x + s_q; independent: U = u_y + u_w + u_q."""
    k2 = K_DB ** 2
    return {
        "s_x": k2 / (2 * J * g_r),
        "s_q": k2 / (J * (g_r + g_t) ** 2),
        "u_y": k2 / (2 * J * g_t),
        "u_w": k2 / (2 * J * g_t * g_r) * ((g_r - g_t) / (g_r + g_t)) ** 2,
        "u_q": k2 / (J * (g_r + g_t) ** 2),
    }


def complete_shared_independent(g_r, g_t, J):
    c = complete_components(g_r, g_t, J)
    return c["s_x"] + c["s_q"], c["u_y"] + c["u_w"] + c["u_q"]


def complete_variance(g_r, g_t, J, n):
    S, U = complete_shared_independent(g_r, g_t, J)
    return S + U / n


def local_expected(g_r, g_t, J):
    """Expected (shared, independent) terms of the fixed-weight local propagation (white noise)."""
    return K_DB ** 2 / (2 * J * (g_r + 1.0)), K_DB ** 2 / (2 * J * g_t * (1.0 + 1.0 / g_r))


def local_terms(r, a_hat, v0, vt):
    """Fixed-weight local propagation: per-observation variance nu_t and shared covariance C_0 (dB^2)."""
    J = r.shape[-1]
    v0 = np.broadcast_to(np.asarray(v0, dtype=float), (J,))
    vt = np.broadcast_to(np.asarray(vt, dtype=float), (J,))
    absr2 = np.abs(r) ** 2
    w = 1.0 / (vt[None, :] + (a_hat[:, None] ** 2) * v0[None, :])
    D = np.sum(w * absr2[None, :], axis=1)
    nu = K_DB ** 2 / (2.0 * a_hat ** 2 * D ** 2) * np.sum(w ** 2 * absr2[None, :] * vt[None, :], axis=1)
    ell = K_DB * w * np.conj(r)[None, :] / D[:, None]
    C0 = 0.5 * np.real((ell * v0[None, :]) @ np.conj(ell).T)
    return nu, C0


def local_variance(nu, C0, n):
    """Variance of the mean of the first n estimates: (local with shared covariance, wrongly independent)."""
    shared = C0[:n, :n].sum()
    return (nu[:n].sum() + shared) / n ** 2, (nu[:n].sum() + np.trace(C0[:n, :n])) / n ** 2


def plugin_complete_variance(r, A_hat, v0, v, n):
    """Complete propagation evaluated with observable quantities only (white noise).

    H is replaced by max(R - J v_0, 0.05 R) and a_t by the estimates; when a_t varies over the window the
    shared noise-energy coefficient is averaged over the window and the independent terms are summed.
    Returns (variance, shared part).
    """
    J = r.shape[-1]
    R = float(np.sum(np.abs(r) ** 2))
    H = max(R - J * v0, 0.05 * R)
    a2 = 10.0 ** (-np.asarray(A_hat[:n], dtype=float) / 10.0)
    g_r, g_t = H / (J * v0), a2 * H / (J * v)
    c = complete_components(g_r, g_t, J)
    shared = K_DB ** 2 * v0 / (2 * H) + (np.mean(K_DB * v / (H * (v + a2 * v0)))) ** 2 * J * v0 ** 2
    indep = np.sum(c["u_y"] + c["u_w"] + c["u_q"]) / n ** 2
    return shared + indep, shared
