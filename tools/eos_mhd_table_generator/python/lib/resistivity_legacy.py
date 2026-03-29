#!/usr/bin/env python3
"""
Legacy-style resistivity model (eta_and_As-inspired) for EOS table helper.

This module ports the core structure from 1_make/resistivity/eta_and_As.f90:
- legacy species charge / mass table input (59 gas species)
- dust charge bins appended as [-2, -1, 0, +1, +2]
- sigma_pinto-like collision times
- eta_OHM / eta_AD / eta_HALL construction

Scope:
- Input abundances are y64 (legacy-compatible 64 species view).
- Numerical guards are added to keep values finite in Python runs.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np


N_CSP = 59
N_DSP = 5
N_SP = N_CSP + N_DSP

PI = 3.14159265358979
C0 = 2.99792458e10
RMASU = 1.660538782e-24
BOLTZ = 1.380662e-16
QE = 4.80653e-10

G_DEN = 3.0
A_MIN = 5e-7
A_MID = 1e-4
A_MAX = 5e-4
P1 = 3.5
P2 = 5.5

EPS = 1.0e-300
LD = np.longdouble


@dataclass(frozen=True)
class LegacyTables:
    q: np.ndarray  # shape (64,)
    m: np.ndarray  # shape (64,), [g]


def _read_flat_floats(path: Path) -> np.ndarray:
    vals = []
    for raw in path.read_text(encoding="ascii").splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        vals.extend(float(x) for x in s.split())
    return np.asarray(vals, dtype=float)


def dustmass() -> float:
    num = (
        (1.0 - (A_MID / A_MIN) ** (P1 - 4.0)) / (P1 - 4.0)
        + ((A_MID / A_MIN) ** (P2 - 4.0)) / (P2 - 4.0)
        - 1.0
    )
    den = (
        (1.0 - (A_MID / A_MIN) ** (P1 - 1.0)) / (P1 - 1.0)
        + ((A_MID / A_MIN) ** (P2 - 1.0)) / (P2 - 1.0)
        - 1.0
    )
    return 4.0 * PI / 3.0 * G_DEN * A_MID**3 * num / max(den, EPS)


def load_legacy_tables(charge_table_path: Path, mass_table_path: Path) -> LegacyTables:
    q59 = _read_flat_floats(charge_table_path)
    mx59 = _read_flat_floats(mass_table_path)
    if q59.size != N_CSP:
        raise ValueError(f"charge table must have {N_CSP} values, got {q59.size}")
    if mx59.size != N_CSP:
        raise ValueError(f"mass table must have {N_CSP} values, got {mx59.size}")

    q = np.zeros(N_SP, dtype=float)
    m = np.zeros(N_SP, dtype=float)
    q[:N_CSP] = q59
    # qext(i-3), i=1..5 => [-2,-1,0,1,2]
    q[N_CSP:] = np.asarray([-2.0, -1.0, 0.0, 1.0, 2.0], dtype=float)

    m[:N_CSP] = mx59 * RMASU
    md = dustmass()
    m[N_CSP:] = md
    return LegacyTables(q=q, m=m)


def _sigma_pinto_tau(xnH: float, T: float, y: np.ndarray, m: np.ndarray, q: np.ndarray) -> np.ndarray:
    # 1-based comments below are preserved from Fortran for traceability.
    xnH = LD(xnH)
    T = LD(T)
    y = np.asarray(y, dtype=LD)
    m = np.asarray(m, dtype=LD)
    q = np.asarray(q, dtype=LD)
    tau = np.full(N_SP, LD(1.0e70), dtype=LD)
    theta = np.log10(max(T, LD(1.0)))

    # Polarizabilities (only H,H2,He used)
    alp1 = 6.67e-25
    alp2 = 8.08e-25
    alp8 = 2.05e-25

    # i = 1..59
    for i in range(N_CSP):
        if q[i] <= 0.0:
            continue
        # rates_pol for H, H2, He
        mu1 = m[0] * m[i] / max(m[0] + m[i], EPS)
        mu2 = m[1] * m[i] / max(m[1] + m[i], EPS)
        mu8 = m[7] * m[i] / max(m[7] + m[i], EPS)

        rate1 = 2.41 * PI * np.sqrt(QE * QE * alp1 / max(mu1, EPS))
        rate2 = 2.41 * PI * np.sqrt(QE * QE * alp2 / max(mu2, EPS))
        rate8 = 2.41 * PI * np.sqrt(QE * QE * alp8 / max(mu8, EPS))

        y1 = max(y[0], EPS)
        y2 = max(y[1], EPS)
        y8 = max(y[7], EPS)
        t1 = (m[0] + m[i]) / max(m[0], EPS) / max(rate1 * xnH * y1, EPS)
        t2 = (m[1] + m[i]) / max(m[1], EPS) / max(rate2 * xnH * y2, EPS)
        t8 = (m[7] + m[i]) / max(m[7], EPS) / max(rate8 * xnH * y8, EPS)
        tau[i] = 1.0 / max(1.0 / max(t2, EPS) + 1.0 / max(t1, EPS) + 1.0 / max(t8, EPS), EPS)

    # Explicit rates in sigma_pinto
    rates = np.zeros(12, dtype=LD)
    rates[0] = np.sqrt(max(T, 0.0)) * (1.476 - 1.409 * theta + 0.555 * theta**2 - 0.0775 * theta**3)  # HCO+,H2
    rates[1] = 2.693 - 1.238 * theta + 0.664 * theta**2 - 0.089 * theta**3                              # H3+,H2
    rates[2] = 1.003 + 0.050 * theta + 0.136 * theta**2 - 0.014 * theta**3                              # H+,H2
    rates[3] = np.sqrt(max(T, 0.0)) * (0.535 + 0.203 * theta - 0.163 * theta**2 + 0.050 * theta**3)    # e,H2
    rates[4] = 1.983 + 0.425 * theta - 0.431 * theta**2 + 0.114 * theta**3                              # C+,H
    rates[5] = 0.649 * max(T, 0.0) ** 0.375                                                              # H+,H
    rates[6] = np.sqrt(max(T, 0.0)) * (2.841 + 0.093 * theta + 0.245 * theta**2 - 0.089 * theta**3)    # e,H
    rates[7] = 1.424 + 7.438e-6 * T - 6.734e-9 * T**2                                                    # H+,He
    rates[8] = 0.428 * np.sqrt(max(T, 0.0))                                                              # e,He

    core = (
        1.3 * 4.0 * PI / 3.0
        * np.sqrt(8.0 * BOLTZ * max(T, 0.0) / (PI * max(m[0], EPS)))
        * A_MID**2
    )
    num = ((1.0 - (A_MID / A_MIN) ** (P1 - 3.0)) / (P1 - 3.0) + ((A_MID / A_MIN) ** (P2 - 3.0)) / (P2 - 3.0) - 1.0)
    den = ((1.0 - (A_MID / A_MIN) ** (P1 - 1.0)) / (P1 - 1.0) + ((A_MID / A_MIN) ** (P2 - 1.0)) / (P2 - 1.0) - 1.0)
    rates[9] = core * num / max(den, EPS)                 # dust,H
    rates[10] = rates[9] * np.sqrt(max(m[0] / max(m[1], EPS), EPS))  # dust,H2
    rates[11] = rates[9] * np.sqrt(max(m[0] / max(m[7], EPS), EPS))  # dust,He

    rates[:9] *= 1.0e-9

    y1 = max(y[0], EPS)
    y2 = max(y[1], EPS)
    y8 = max(y[7], EPS)

    # tau(45), tau(6), tau(4), tau(3), tau(23) in 1-based index
    # 0-based: 44,5,3,2,22
    tau[44] = (m[1] + m[44]) / max(m[1], EPS) / max(rates[0] * xnH * y2, EPS)
    tau[5] = (m[1] + m[5]) / max(m[1], EPS) / max(rates[1] * xnH * y2, EPS)

    t4a = (m[1] + m[3]) / max(m[1], EPS) / max(rates[2] * xnH * y2, EPS)
    t4b = (m[0] + m[3]) / max(m[0], EPS) / max(rates[5] * xnH * y1, EPS)
    t4c = (m[7] + m[3]) / max(m[7], EPS) / max(rates[7] * xnH * y8, EPS)
    tau[3] = 1.0 / max(1.0 / max(t4a, EPS) + 1.0 / max(t4b, EPS) + 1.0 / max(t4c, EPS), EPS)

    t3a = (m[1] + m[2]) / max(m[1], EPS) / max(rates[3] * xnH * y2, EPS)
    t3b = (m[0] + m[2]) / max(m[0], EPS) / max(rates[6] * xnH * y1, EPS)
    t3c = (m[7] + m[2]) / max(m[7], EPS) / max(rates[8] * xnH * y8, EPS)
    tau[2] = 1.0 / max(1.0 / max(t3a, EPS) + 1.0 / max(t3b, EPS) + 1.0 / max(t3c, EPS), EPS)

    tau[22] = (m[1] + m[22]) / max(m[1], EPS) / max(rates[4] * xnH * y2, EPS)

    # charged dust bins 60..64 (1-based) => 59..63 (0-based)
    for i in range(N_CSP, N_CSP + N_DSP):
        ta = (m[1] + m[i]) / max(m[1], EPS) / max(rates[9] * xnH * y2, EPS)
        tb = (m[0] + m[i]) / max(m[0], EPS) / max(rates[10] * xnH * y1, EPS)
        tc = (m[7] + m[i]) / max(m[7], EPS) / max(rates[11] * xnH * y8, EPS)
        tau[i] = 1.0 / max(1.0 / max(ta, EPS) + 1.0 / max(tb, EPS) + 1.0 / max(tc, EPS), EPS)

    # neutral dust override: tau(N_csp+3)=1e70 -> idx 61
    tau[N_CSP + 2] = LD(1.0e70)
    return tau


def eta_and_as_single(xnH: float, T: float, y64: np.ndarray, B: float, tables: LegacyTables) -> Tuple[float, float, float]:
    q = np.asarray(tables.q, dtype=LD)
    m = np.asarray(tables.m, dtype=LD)
    y = np.asarray(y64, dtype=LD)
    xnH = LD(xnH)
    T = LD(T)
    B = LD(B)

    omega = q * LD(QE) * B / (m * LD(C0))
    tau = _sigma_pinto_tau(xnH, T, y, m, q)

    sigma_nu = np.zeros(N_SP, dtype=LD)
    valid = tau < LD(1.0e60)
    sigma_nu[valid] = y[valid] * xnH * (q[valid] * LD(QE)) ** 2 * tau[valid] / np.maximum(m[valid], LD(EPS))
    tau_eff = tau.copy()
    tau_eff[~valid] = LD(1.0)

    sigma_ohm = max(np.sum(sigma_nu), LD(EPS))
    eta_ohm = LD(C0) ** 2 / (LD(4.0) * LD(PI)) / sigma_ohm

    tw = tau_eff * omega
    den = LD(1.0) + tw**2
    pref = xnH * m * y

    A1 = np.sum(pref * tau_eff * omega**2 / den)
    B1 = np.sum(pref * omega * (tw**3) / den)
    B2 = np.sum(pref * omega * (tw**2) / den)
    A2 = -np.sum(pref * omega * (tw**2) / den)

    sigma_p = A1 * (LD(C0) / max(B, LD(EPS))) ** 2
    sigma_h = A2 * (LD(C0) / max(B, LD(EPS))) ** 2
    den_hall = max(sigma_h**2 + sigma_p**2, LD(EPS))
    eta_hall = LD(C0) ** 2 / (LD(4.0) * LD(PI)) * sigma_h / den_hall

    den_A = max(A1**2 + A2**2, LD(EPS))
    Dratio = (A1 * B1 + B2**2 + LD(2.0) * A2 * B2) / den_A
    eta_ad = LD(C0) ** 2 / (LD(4.0) * LD(PI)) * Dratio / sigma_ohm
    return float(eta_ohm), float(eta_ad), float(eta_hall)


def eta_grid(logn: np.ndarray, T: np.ndarray, y64: np.ndarray, B: float, tables: LegacyTables) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = logn.size
    eta_o = np.empty(n, dtype=float)
    eta_a = np.empty(n, dtype=float)
    eta_h = np.empty(n, dtype=float)
    for i in range(n):
        xnH = 10.0 ** float(logn[i])
        e0, e1, e2 = eta_and_as_single(xnH, float(T[i]), y64[i, :], float(B), tables)
        eta_o[i] = max(abs(e0), EPS)
        eta_a[i] = max(abs(e1), EPS)
        eta_h[i] = max(abs(e2), EPS)
    return eta_o, eta_a, eta_h
