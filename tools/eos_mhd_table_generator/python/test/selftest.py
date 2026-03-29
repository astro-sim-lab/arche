#!/usr/bin/env python3
"""
selftest.py - lightweight checks for eos_table helper modules.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from lib.mapping import map89_to_64
from lib.resistivity_legacy import eta_and_as_single, load_legacy_tables


def test_mapping() -> None:
    y89 = np.zeros((1, 89), dtype=float)
    # direct region
    y89[0, 0] = 1.0
    y89[0, 56] = 2.0
    # Mg/Mg+
    y89[0, 61] = 3.0
    y89[0, 62] = 4.0
    # dust bins Gr2-, Gr-, Gr, Gr+, Gr2+
    y89[0, 67] = 5.0
    y89[0, 66] = 6.0
    y89[0, 63] = 7.0
    y89[0, 64] = 8.0
    y89[0, 65] = 9.0

    y64 = map89_to_64(y89)
    assert y64.shape == (1, 64)
    assert y64[0, 0] == 1.0
    assert y64[0, 56] == 2.0
    assert y64[0, 57] == 3.0
    assert y64[0, 58] == 4.0
    assert y64[0, 59] == 5.0
    assert y64[0, 60] == 6.0
    assert y64[0, 61] == 7.0
    assert y64[0, 62] == 8.0
    assert y64[0, 63] == 9.0


def test_legacy_eta_onepoint() -> None:
    charge_path = Path("tools/eos_mhd_table_generator/include/legacy_species_charge.dat")
    mass_path = Path("tools/eos_mhd_table_generator/include/legacy_species_mass_amu.dat")
    tbl = load_legacy_tables(charge_path, mass_path)
    assert tbl.q.shape[0] == 64
    assert tbl.m.shape[0] == 64

    y = np.zeros(64, dtype=float)
    # simple physically plausible composition
    y[0] = 0.9    # H
    y[1] = 0.05   # H2
    y[2] = 1e-7   # e
    y[3] = 1e-7   # H+
    y[7] = 8.33e-2  # He
    y[22] = 1e-9    # C+
    y[61] = 1e-12   # neutral grain

    e0, e1, e2 = eta_and_as_single(1.0e4, 100.0, y, 1.0e-10, tbl)
    assert np.isfinite(e0) and np.isfinite(e1) and np.isfinite(e2)


def main() -> int:
    test_mapping()
    test_legacy_eta_onepoint()
    print("selftest: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
