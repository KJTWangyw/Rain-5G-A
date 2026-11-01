"""Pilot-based rain attenuation estimation with a reused clear-sky reference (5G NR DM-RS).

Modules
-------
waveform     NR DM-RS pilot observations through an OFDM chain (TS 38.211, numerology 3)
estimator    profile amplitude estimator (closed-form root and numeric solver)
uncertainty  information bound, fixed-weight local propagation, complete closed-form propagation
bootstrap    parametric bootstrap with a shared reference draw
rain         ITU-R P.838-3 specific-attenuation coefficients
stats        Wilson interval and exact McNemar test
"""

__version__ = "1.0.0"
