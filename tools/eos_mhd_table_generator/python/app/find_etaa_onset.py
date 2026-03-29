#!/usr/bin/env python3
"""
find_etaa_onset.py - Detect onset density where etaA absolute difference exceeds threshold.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from lib.case_paths import resolve_case_file
from lib.table_io import EOS_COLS, read_table, to_keyed_nt_bt


def _read_cases(path: Path) -> list[tuple[str, str, str]]:
    out: list[tuple[str, str, str]] = []
    with path.open("r", encoding="utf-8") as f:
        filtered = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    reader = csv.DictReader(filtered)
    for r in reader:
        out.append((r["cr_scale"].strip(), r["sh_scale"].strip(), r["z_metal"].strip()))
    return out


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Find etaA discrepancy onset density per case.")
    p.add_argument("--cases", required=True)
    p.add_argument("--base-dir", default=".dev/low-metal-eos/ref/Lmetal/data/table")
    p.add_argument("--target-dir", default="tools/eos_mhd_table_generator/results")
    p.add_argument("--out-csv", default="tools/eos_mhd_table_generator/results/etaa_onset.csv")
    p.add_argument("--abs-diff-threshold", type=float, default=3.0)
    p.add_argument("--nt-min", type=float, default=14.0)
    return p


def main() -> int:
    args = _parser().parse_args()
    rows_out: list[list[str]] = []
    for cr, sh, z in _read_cases(Path(args.cases)):
        req_case = f"CR{cr}_sh{sh}_Z{z}.d"
        base, base_name = resolve_case_file(Path(args.base_dir), cr, sh, z)
        target, target_name = resolve_case_file(Path(args.target_dir), cr, sh, z)
        if base is None or target is None:
            rows_out.append([req_case, base_name if base else "", target_name if target else "", "missing", "", "", "", ""])
            continue

        a = to_keyed_nt_bt(read_table(base, strict_cols=False))
        b = to_keyed_nt_bt(read_table(target, strict_cols=False))
        keys = sorted(set(a.keys()) & set(b.keys()))
        onset_nt = None
        onset_bt = None
        onset_abs = None
        max_abs = -1.0
        max_nt = None
        max_bt = None
        for nt, bt in keys:
            if nt < args.nt_min:
                continue
            d = float(b[(nt, bt)][EOS_COLS["etaA"]] - a[(nt, bt)][EOS_COLS["etaA"]])
            ad = abs(d)
            if ad > max_abs:
                max_abs = ad
                max_nt = nt
                max_bt = bt
            if onset_nt is None and ad >= args.abs_diff_threshold:
                onset_nt = nt
                onset_bt = bt
                onset_abs = ad
        status = "ok"
        if onset_nt is None:
            status = "no_onset"
        rows_out.append(
            [
                req_case,
                base_name,
                target_name,
                status,
                "" if onset_nt is None else f"{onset_nt:.6f}",
                "" if onset_bt is None else f"{onset_bt:.6f}",
                "" if onset_abs is None else f"{onset_abs:.6e}",
                "" if max_abs < 0 else f"{max_abs:.6e}",
                "" if max_nt is None else f"{max_nt:.6f}",
                "" if max_bt is None else f"{max_bt:.6f}",
            ]
        )

    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="ascii", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "case_request",
                "case_base",
                "case_target",
                "status",
                "onset_nt",
                "onset_bt",
                "onset_abs_diff",
                "max_abs_diff",
                "max_abs_diff_nt",
                "max_abs_diff_bt",
            ]
        )
        w.writerows(rows_out)
    print(f"[OK] etaa_onset_csv={out}")
    print(f"[OK] cases={len(rows_out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

