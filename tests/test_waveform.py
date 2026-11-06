import numpy as np

from nrpilot.waveform import Multipath, NRConfig, dmrs_symbols, frequency_response, observe_slots, pilot_channel

CFG = NRConfig()
TAPS = Multipath().taps(CFG)


def test_allocation_has_144_pilots():
    assert CFG.n_pilots == 144
    assert len(CFG.pilot_sc) == 72


def test_dmrs_is_unit_modulus_qpsk():
    s = dmrs_symbols(CFG, 5)
    assert s.shape == (2, 72)
    np.testing.assert_allclose(np.abs(s), 1.0, atol=1e-12)
    assert set(np.round(s.real * np.sqrt(2)).astype(int).ravel()) <= {-1, 1}


def test_channel_normalised_to_unit_mean_gain():
    np.testing.assert_allclose(np.mean(np.abs(frequency_response(CFG, TAPS)) ** 2), 1.0, rtol=1e-12)


def test_noise_free_extraction_recovers_scaled_channel():
    rng = np.random.default_rng(0)
    amp, phase = np.array([1.0, 0.7]), np.array([0.0, 1.3])
    z = observe_slots(CFG, TAPS, amp, phase, 0.0, rng, first_slot=3)
    h = pilot_channel(CFG, TAPS)
    np.testing.assert_allclose(z, (amp * np.exp(1j * phase))[:, None] * h[None, :], atol=1e-12)


def test_per_pilot_noise_variance_equals_sample_variance():
    rng = np.random.default_rng(1)
    z = observe_slots(CFG, TAPS, np.zeros(400), np.zeros(400), 0.5, rng)
    assert abs(np.var(z) / 0.5 - 1.0) < 0.02
    assert abs(np.mean(z)) < 0.01
