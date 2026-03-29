#!/usr/bin/env python3
"""
compare_case_matrix.py - Compare many EOS .d cases against reference tables.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from lib.compare_core import USE_COLS, calc_metrics
from lib.case_paths import resolve_case_file


def _read_cases(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    with path.open("r", encoding="utf-8") as f:
        filtered = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    if not filtered:
        raise ValueError(f"no cases in {path}")
    reader = csv.DictReader(filtered)
    required = ("cr_scale", "sh_scale", "z_metal")
    for key in required:
        if key not in (reader.fieldnames or []):
            raise ValueError(f"{path}: missing required column '{key}'")
    for r in reader:
        rows.append((r["cr_scale"].strip(), r["sh_scale"].strip(), r["z_metal"].strip()))
    if not rows:
        raise ValueError(f"no case rows in {path}")
    return rows


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Batch-compare generated EOS tables against reference tables.")
    p.add_argument("--cases", required=True, help="CSV file with columns: cr_scale,sh_scale,z_metal")
    p.add_argument("--base-dir", default=".dev/low-metal-eos/ref/Lmetal/data/table", help="Directory containing reference .d files")
    p.add_argument("--target-dir", default="tools/eos_mhd_table_generator/results", help="Directory containing generated .d files")
    p.add_argument("--out-csv", default="tools/eos_mhd_table_generator/results/compare_matrix_summary.csv")
    p.add_argument("--strict", action="store_true", help="Return non-zero if any case is missing or fails thresholds.")
    p.add_argument("--corr-threshold", type=float, help="Optional threshold applied to all tracked quantities.")
    p.add_argument("--med-abs-threshold", type=float, help="Optional threshold applied to all tracked quantities.")
    return p


def main() -> int:
    args = _parser().parse_args()
    cases = _read_cases(Path(args.cases))
    base_dir = Path(args.base_dir)
    target_dir = Path(args.target_dir)
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    header = ["case_request", "case_base", "case_target", "status", "overlap_points"]
    for q in USE_COLS:
        header.extend([f"{q}_corr", f"{q}_med_abs", f"{q}_q05", f"{q}_q50", f"{q}_q95"])

    fail_any = False
    rows_out: list[list[str]] = []
    for cr, sh, z in cases:
        req_case = f"CR{cr}_sh{sh}_Z{z}.d"
        base, base_name = resolve_case_file(base_dir, cr, sh, z)
        target, target_name = resolve_case_file(target_dir, cr, sh, z)
        if base is None or target is None:
            status = "missing_base" if base is None else "missing_target"
            print(f"[MISS] {req_case}: {status}")
            row = [req_case, base_name if base is not None else "", target_name if target is not None else "", status, "0"] + [""] * (5 * len(USE_COLS))
            rows_out.append(row)
            fail_any = True
            continue

        overlap, metrics = calc_metrics(base, target)
        status = "ok"
        if args.corr_threshold is not None:
            for q in USE_COLS:
                c = metrics[q]["corr"]
                if (not (c == c)) or c < args.corr_threshold:  # NaN check via c==c
                    status = "fail_corr"
        if args.med_abs_threshold is not None:
            for q in USE_COLS:
                if metrics[q]["med_abs"] > args.med_abs_threshold:
                    status = "fail_med_abs"
        if status != "ok":
            fail_any = True
        print(f"[{status.upper()}] {req_case}: overlap={overlap}")
        row = [req_case, base_name, target_name, status, str(overlap)]
        for q in USE_COLS:
            row.extend(
                [
                    f"{metrics[q]['corr']:.6f}",
                    f"{metrics[q]['med_abs']:.6e}",
                    f"{metrics[q]['q05']:.6e}",
                    f"{metrics[q]['q50']:.6e}",
                    f"{metrics[q]['q95']:.6e}",
                ]
            )
        rows_out.append(row)

    with out_csv.open("w", encoding="ascii", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows_out)
    print(f"[OK] summary_csv={out_csv}")
    if args.strict and fail_any:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
