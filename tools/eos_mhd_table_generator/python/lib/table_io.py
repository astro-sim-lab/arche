"""
Common EOS-table (.d) I/O and indexing helpers.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

import numpy as np


EOS_COLS = {
    "ii": 0,
    "jj": 1,
    "nt": 2,
    "tt": 3,
    "mut": 4,
    "pt": 5,
    "bt": 6,
    "etaO": 7,
    "etaA": 8,
    "etaH": 9,
}


def read_table(path: Path, *, strict_cols: bool) -> np.ndarray:
    rows = []
    with path.open("r", encoding="ascii") as f:
        for ln, line in enumerate(f, start=1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            cols = s.split()
            if len(cols) != 10:
                if strict_cols:
                    raise ValueError(f"line {ln}: expected 10 columns, got {len(cols)}")
                continue
            try:
                rows.append([float(x) for x in cols])
            except ValueError as exc:
                token_preview = " ".join(cols[:3])
                if len(token_preview) > 120:
                    token_preview = token_preview[:117] + "..."
                raise ValueError(
                    f"line {ln}: non-float token(s) detected in {path.name}: '{token_preview}'"
                ) from exc
    if not rows:
        raise ValueError(f"no valid rows in {path}")
    return np.asarray(rows, dtype=float)


def to_keyed_nt_bt(arr: np.ndarray) -> Dict[Tuple[float, float], np.ndarray]:
    """Key rows by (nt, bt), keeping the first hit per key."""
    out: Dict[Tuple[float, float], np.ndarray] = {}
    nt_idx = EOS_COLS["nt"]
    bt_idx = EOS_COLS["bt"]
    for r in arr:
        key = (r[nt_idx], r[bt_idx])
        if key not in out:
            out[key] = r
    return out


def monotonic_non_decreasing(x: np.ndarray) -> bool:
    return np.all(np.diff(x) >= -1e-12)


def split_blocks_by_ii(ii: np.ndarray) -> list[tuple[int, int]]:
    starts = [0]
    for k in range(1, ii.size):
        if ii[k] == 1 and ii[k - 1] != 1:
            starts.append(k)
    spans: list[tuple[int, int]] = []
    for s, e in zip(starts, starts[1:] + [ii.size]):
        spans.append((s, e))
    return spans


def extract_series_by_bt(arr: np.ndarray, *, bt: float, col_name: str, atol: float = 1.0e-8) -> tuple[np.ndarray, np.ndarray]:
    if col_name not in EOS_COLS:
        raise ValueError(f"unknown column name: {col_name}")
    bt_idx = EOS_COLS["bt"]
    nt_idx = EOS_COLS["nt"]
    val_idx = EOS_COLS[col_name]
    m = np.isclose(arr[:, bt_idx], bt, atol=atol)
    s = arr[m]
    if s.size == 0:
        raise ValueError(f"no rows for bt={bt}")
    order = np.argsort(s[:, nt_idx])
    x = s[order, nt_idx]
    y = s[order, val_idx]
    return x, y
