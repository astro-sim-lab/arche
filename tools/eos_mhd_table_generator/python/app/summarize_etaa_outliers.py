#!/usr/bin/env python3
"""
summarize_etaa_outliers.py - Summarize etaA outlier CSV by case.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Summarize etaA diagnostics by case. "
            "Supports both etaa_outliers.csv and etaa_onset.csv schemas."
        )
    )
    p.add_argument("--in-csv", default="tools/eos_mhd_table_generator/results/etaa_outliers.csv")
    p.add_argument("--out-csv", default="tools/eos_mhd_table_generator/results/etaa_outliers_summary.csv")
    return p


def main() -> int:
    args = _parser().parse_args()
    src = Path(args.in_csv)
    dst = Path(args.out_csv)
    with src.open("r", encoding="ascii", newline="") as f:
        rows = list(csv.DictReader(f))
    by_case: dict[str, list[dict[str, str]]] = {}
    for r in rows:
        k = r["case_request"]
        by_case.setdefault(k, []).append(r)

    if not rows:
        raise SystemExit(f"[ERROR] empty csv: {src}")

    has_outlier_schema = "abs_diff" in rows[0]
    has_onset_schema = "onset_nt" in rows[0] and "max_abs_diff" in rows[0]
    if not (has_outlier_schema or has_onset_schema):
        raise SystemExit(
            "[ERROR] unsupported csv schema. "
            "Expected columns for etaa_outliers.csv or etaa_onset.csv."
        )

    out_rows: list[list[str]] = []
    for case, rs in sorted(by_case.items()):
        if has_outlier_schema:
            abs_vals = [float(r["abs_diff"]) for r in rs]
            nt_vals = [float(r["nt"]) for r in rs]
            bt_vals = [float(r["bt"]) for r in rs]
            rmax = max(rs, key=lambda r: float(r["abs_diff"]))
            n_outliers = len(rs)
        else:
            # etaa_onset.csv: each row is already one case-level summary
            abs_vals = [float(r["max_abs_diff"]) for r in rs]
            nt_vals = [float(r["onset_nt"]) for r in rs if r["onset_nt"]]
            bt_vals = [float(r["onset_bt"]) for r in rs if r["onset_bt"]]
            rmax = max(rs, key=lambda r: float(r["max_abs_diff"]))
            n_outliers = 1
            if not nt_vals:
                nt_vals = [float(rmax["max_abs_diff_nt"])]
            if not bt_vals:
                bt_vals = [float(rmax["max_abs_diff_bt"])]

        out_rows.append(
            [
                case,
                str(n_outliers),
                f"{max(abs_vals):.6e}",
                rmax["max_abs_diff_nt"] if has_onset_schema else rmax["nt"],
                rmax["max_abs_diff_bt"] if has_onset_schema else rmax["bt"],
                f"{min(nt_vals):.6f}",
                f"{max(nt_vals):.6f}",
                f"{min(bt_vals):.6f}",
                f"{max(bt_vals):.6f}",
            ]
        )

    dst.parent.mkdir(parents=True, exist_ok=True)
    with dst.open("w", encoding="ascii", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "case_request",
                "n_outliers",
                "max_abs_diff",
                "max_abs_diff_nt",
                "max_abs_diff_bt",
                "nt_min",
                "nt_max",
                "bt_min",
                "bt_max",
            ]
        )
        w.writerows(out_rows)
    print(f"[OK] summary_csv={dst}")
    print(f"[OK] cases={len(out_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
