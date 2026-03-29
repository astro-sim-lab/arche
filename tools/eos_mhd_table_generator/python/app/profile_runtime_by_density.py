#!/usr/bin/env python3
"""
profile_runtime_by_density.py - Approximate runtime concentration by density bands.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import numpy as np


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Profile matrix runtime by density bands (row-weighted approximation).")
    p.add_argument("--runtime-csv", default="tools/eos_mhd_table_generator/results/matrix_runtime.csv")
    p.add_argument("--out-csv", default="tools/eos_mhd_table_generator/results/matrix_runtime_density_profile.csv")
    p.add_argument("--bands", default="14,16,18,20,22",
                   help="Comma-separated upper boundaries for log10(nH) bands; lower bound starts at -inf.")
    p.add_argument(
        "--include-nonzero",
        action="store_true",
        help="Include rows even when exit_code != 0 (e.g., timeout-limited profiling runs).",
    )
    return p


def _band_edges(bands_str: str) -> list[float]:
    vals = [float(x.strip()) for x in bands_str.split(",") if x.strip()]
    return sorted(vals)


def _band_label(lo: float, hi: float) -> str:
    if np.isneginf(lo):
        return f"(-inf,{hi}]"
    return f"({lo},{hi}]"


def main() -> int:
    args = _parser().parse_args()
    edges = _band_edges(args.bands)
    rows = list(csv.DictReader(Path(args.runtime_csv).open("r", encoding="ascii", newline="")))
    out_rows: list[list[str]] = []

    for r in rows:
        if (not args.include_nonzero) and r.get("exit_code", "") != "0":
            continue
        h5 = r.get("input_h5", "").strip()
        if not h5:
            continue
        h5p = Path(h5)
        if not h5p.is_file():
            continue
        elapsed = float(r["elapsed_sec"])
        with h5py.File(h5p, "r") as f:
            xnh = np.asarray(f["xnH"][()], dtype=float)
        logn = np.log10(np.maximum(xnh, 1.0e-300))
        n_total = logn.size
        if n_total == 0:
            continue

        bounds = [-np.inf] + edges
        for lo, hi in zip(bounds[:-1], bounds[1:]):
            m = (logn > lo) & (logn <= hi)
            n = int(np.count_nonzero(m))
            frac = n / n_total
            elapsed_alloc = elapsed * frac
            out_rows.append(
                [
                    r["case"],
                    _band_label(lo, hi),
                    str(n),
                    str(n_total),
                    f"{frac:.6f}",
                    f"{elapsed_alloc:.6f}",
                    f"{elapsed:.6f}",
                ]
            )
        lo = edges[-1]
        m = logn > lo
        n = int(np.count_nonzero(m))
        frac = n / n_total
        elapsed_alloc = elapsed * frac
        out_rows.append(
            [
                r["case"],
                f"({lo},+inf)",
                str(n),
                str(n_total),
                f"{frac:.6f}",
                f"{elapsed_alloc:.6f}",
                f"{elapsed:.6f}",
            ]
        )

    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="ascii", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "logn_band", "n_rows_band", "n_rows_total", "row_fraction", "elapsed_sec_alloc", "elapsed_sec_case"])
        w.writerows(out_rows)
    print(f"[OK] density_runtime_csv={out}")
    print(f"[OK] rows={len(out_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
