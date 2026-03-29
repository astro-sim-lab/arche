"""
Species mapping helpers for EOS-table tools.
"""

from __future__ import annotations

import numpy as np


def map89_to_64(y89: np.ndarray) -> np.ndarray:
    """
    Map arche-dev metal network (89 species) to legacy 64-species view.
    Input shape: (N, 89). Output shape: (N, 64).
    """
    if y89.ndim != 2:
        raise ValueError("y89 must be 2D")
    if y89.shape[1] < 89:
        raise ValueError("expected y with 89 species for metal_grain mapping")

    n = y89.shape[0]
    y64 = np.zeros((n, 64), dtype=float)

    # Fortran 1..57 <= arche 1..57  (0-based: 0..56)
    y64[:, 0:57] = y89[:, 0:57]

    # Fortran 58..59 (M, M+) <= arche 62..63 (Mg, Mg+)
    y64[:, 57] = y89[:, 61]
    y64[:, 58] = y89[:, 62]

    # Fortran 60..64 (g--,g-,g,g+,g++) <= arche (Gr2-,Gr-,Gr,Gr+,Gr2+)
    y64[:, 59] = y89[:, 67]  # g--
    y64[:, 60] = y89[:, 66]  # g-
    y64[:, 61] = y89[:, 63]  # g
    y64[:, 62] = y89[:, 64]  # g+
    y64[:, 63] = y89[:, 65]  # g++
    return y64

