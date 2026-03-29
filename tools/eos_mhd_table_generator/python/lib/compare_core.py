"""
Core comparison logic for EOS-table files.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .table_io import EOS_COLS, read_table, to_keyed_nt_bt


USE_COLS = {
    "tt": EOS_COLS["tt"],
    "mut": EOS_COLS["mut"],
    "pt": EOS_COLS["pt"],
    "etaO": EOS_COLS["etaO"],
    "etaA": EOS_COLS["etaA"],
    "etaH": EOS_COLS["etaH"],
}


def _corr(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2:
        return float("nan")
    sx = np.std(x)
    sy = np.std(y)
    if sx == 0.0 or sy == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def calc_metrics(base: Path, target: Path) -> tuple[int, dict[str, dict[str, float]]]:
    a = to_keyed_nt_bt(read_table(base, strict_cols=False))
    b = to_keyed_nt_bt(read_table(target, strict_cols=False))
    keys = sorted(set(a.keys()) & set(b.keys()))
    if not keys:
        raise ValueError("no overlapping (nt, bt) points")

    out: dict[str, dict[str, float]] = {}
    for name, idx in USE_COLS.items():
        x = np.asarray([a[k][idx] for k in keys], dtype=float)
        y = np.asarray([b[k][idx] for k in keys], dtype=float)
        d = y - x
        out[name] = {
            "corr": _corr(x, y),
            "med_abs": float(np.median(np.abs(d))),
            "q05": float(np.percentile(d, 5)),
            "q50": float(np.percentile(d, 50)),
            "q95": float(np.percentile(d, 95)),
        }
    return len(keys), out


def compare_tables(base: Path, target: Path, corr_threshold: float | None = None, med_abs_threshold: float | None = None) -> int:
    overlap_points, metrics = calc_metrics(base, target)
    print(f"overlap_points={overlap_points}")
    failed = False
    for name in USE_COLS:
        corr = metrics[name]["corr"]
        med_abs = metrics[name]["med_abs"]
        print(
            f"{name:5s} "
            f"corr={corr: .4f} "
            f"med_abs={med_abs: .4e} "
            f"q05={metrics[name]['q05']: .4e} "
            f"q50={metrics[name]['q50']: .4e} "
            f"q95={metrics[name]['q95']: .4e}"
        )
        if corr_threshold is not None and (not np.isfinite(corr) or corr < corr_threshold):
            failed = True
            print(f"FAIL: {name} corr {corr:.4f} < threshold {corr_threshold:.4f}")
        if med_abs_threshold is not None and med_abs > med_abs_threshold:
            failed = True
            print(f"FAIL: {name} med_abs {med_abs:.4e} > threshold {med_abs_threshold:.4e}")
    return 2 if failed else 0
