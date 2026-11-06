"""NR downlink DM-RS pilot observations through an OFDM chain.

Numerology 3 (120 kHz subcarrier spacing, 1024-point FFT, 122.88 MHz sampling, normal cyclic prefix of
72 samples). The allocation spans ``n_rb`` resource blocks; PDSCH DM-RS configuration type 1, mapping
type A, ``dmrs-TypeA-Position = pos2`` and one additional position place the pilots on OFDM symbols 2
and 11, on every other subcarrier of CDM group 0 (TS 38.211, Sec. 7.4.1.1). The pilot sequence is the
QPSK-mapped length-31 Gold sequence of TS 38.211 Sec. 5.2.1 with the c_init of Sec. 7.4.1.1.1.

Only the two DM-RS symbols of each slot are simulated, because the receiver uses nothing else. The channel
is a fixed tapped delay line whose taps lie inside the cyclic prefix, scaled by the rain amplitude
``a_t = 10^(-A_t/20)`` and rotated by a common phase per slot. White Gaussian noise is added in the time
domain; with the orthonormal FFT its variance per resource element equals the per-sample variance.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

SLOTS_PER_FRAME = 80  # numerology 3
SYMBOLS_PER_SLOT = 14


@dataclass(frozen=True)
class NRConfig:
    n_fft: int = 1024
    n_cp: int = 72
    scs_hz: float = 120e3
    n_rb: int = 12
    dmrs_symbols: tuple = (2, 11)
    n_id: int = 1
    n_scid: int = 0

    @property
    def n_sc(self) -> int:
        return 12 * self.n_rb

    @property
    def pilot_sc(self) -> np.ndarray:
        """Allocated subcarriers carrying DM-RS (configuration type 1, CDM group 0)."""
        return np.arange(0, self.n_sc, 2)

    @property
    def n_pilots(self) -> int:
        return len(self.pilot_sc) * len(self.dmrs_symbols)

    @property
    def fs_hz(self) -> float:
        return self.n_fft * self.scs_hz

    def fft_bins(self, sc) -> np.ndarray:
        """FFT bin of each allocated subcarrier, with the allocation centred on DC."""
        return (np.asarray(sc) - self.n_sc // 2) % self.n_fft


def gold_sequence(c_init: int, length: int) -> np.ndarray:
    """Pseudo-random sequence c(n) of TS 38.211 Sec. 5.2.1."""
    nc = 1600
    total = nc + length
    x1 = np.zeros(total + 31, dtype=np.uint8)
    x2 = np.zeros(total + 31, dtype=np.uint8)
    x1[0] = 1
    x2[:31] = [(c_init >> i) & 1 for i in range(31)]
    for k in range(total):
        x1[k + 31] = x1[k + 3] ^ x1[k]
        x2[k + 31] = x2[k + 3] ^ x2[k + 2] ^ x2[k + 1] ^ x2[k]
    return x1[nc:total] ^ x2[nc:total]


def dmrs_c_init(slot: int, symbol: int, n_id: int, n_scid: int) -> int:
    """c_init of TS 38.211 Sec. 7.4.1.1.1 (single-layer, without the CDM-group offset term)."""
    return (2 ** 17 * (SYMBOLS_PER_SLOT * slot + symbol + 1) * (2 * n_id + 1) + 2 * n_id + n_scid) % 2 ** 31


@lru_cache(maxsize=None)
def dmrs_symbols(cfg: NRConfig, slot: int) -> np.ndarray:
    """QPSK DM-RS values on the pilot subcarriers, shape (len(dmrs_symbols), n_sc // 2)."""
    m = cfg.n_sc // 2
    out = np.empty((len(cfg.dmrs_symbols), m), dtype=complex)
    for i, sym in enumerate(cfg.dmrs_symbols):
        c = gold_sequence(dmrs_c_init(slot % SLOTS_PER_FRAME, sym, cfg.n_id, cfg.n_scid), 2 * m).astype(float)
        out[i] = ((1 - 2 * c[0::2]) + 1j * (1 - 2 * c[1::2])) / np.sqrt(2)
    return out


@dataclass(frozen=True)
class Multipath:
    """Line-of-sight path plus two weak reflections, all inside the 72-sample cyclic prefix."""
    delays: tuple = (0, 9, 23)
    gains: tuple = (1.0, 0.25 * np.exp(1j * np.pi / 3), 0.12 * np.exp(-2j * np.pi / 5))

    def taps(self, cfg: NRConfig) -> np.ndarray:
        """Tap vector normalised so that the mean |H|^2 over the pilot subcarriers is 1."""
        if max(self.delays) >= cfg.n_cp:
            raise ValueError("all taps must lie inside the cyclic prefix")
        h = np.zeros(max(self.delays) + 1, dtype=complex)
        for d, g in zip(self.delays, self.gains):
            h[d] += g
        return h / np.sqrt(np.mean(np.abs(frequency_response(cfg, h)) ** 2))


def frequency_response(cfg: NRConfig, taps: np.ndarray) -> np.ndarray:
    """Channel on the pilot subcarriers, shape (n_sc // 2,)."""
    k = cfg.fft_bins(cfg.pilot_sc)
    d = np.arange(len(taps))
    return np.exp(-2j * np.pi * np.outer(k, d) / cfg.n_fft) @ taps


def pilot_channel(cfg: NRConfig, taps: np.ndarray) -> np.ndarray:
    """True clear-sky channel h_j for all pilots of a slot (both DM-RS symbols), shape (n_pilots,)."""
    return np.tile(frequency_response(cfg, taps), len(cfg.dmrs_symbols))


def observe_slots(cfg: NRConfig, taps: np.ndarray, amp, phase, noise_var: float, rng, first_slot: int = 0):
    """LS pilot estimates of consecutive slots, shape (n_slots, n_pilots).

    Each slot is generated in the time domain: DM-RS mapping, IFFT, cyclic prefix, tapped-delay-line
    channel scaled by ``amp`` and rotated by ``phase``, white Gaussian noise of per-sample variance
    ``noise_var``; the receiver removes the cyclic prefix, applies the FFT and divides by the known DM-RS.
    """
    amp = np.atleast_1d(np.asarray(amp, dtype=float))
    phase = np.atleast_1d(np.asarray(phase, dtype=float))
    n_slots, n_sym = len(amp), len(cfg.dmrs_symbols)
    bins = cfg.fft_bins(cfg.pilot_sc)
    known = np.stack([dmrs_symbols(cfg, first_slot + s) for s in range(n_slots)])  # (S, n_sym, m)

    grid = np.zeros((n_slots, n_sym, cfg.n_fft), dtype=complex)
    grid[:, :, bins] = known
    x = np.fft.ifft(grid, axis=-1, norm="ortho")
    x = np.concatenate([x[..., -cfg.n_cp:], x], axis=-1)
    y = np.zeros_like(x)
    for d, g in enumerate(taps):
        if g != 0:
            y[..., d:] += g * x[..., :x.shape[-1] - d]
    y *= (amp * np.exp(1j * phase))[:, None, None]
    y += np.sqrt(noise_var / 2.0) * (rng.standard_normal(y.shape) + 1j * rng.standard_normal(y.shape))
    received = np.fft.fft(y[..., cfg.n_cp:], axis=-1, norm="ortho")[..., bins]
    return (received / known).reshape(n_slots, n_sym * len(bins))
