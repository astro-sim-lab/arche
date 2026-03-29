#!/usr/bin/env python3
"""
scan_etaa_outliers.py - Extract high-density etaA outliers across case matrix.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from lib.case_paths import resolve_case_file
from lib.table_io import EOS_COLS, read_table, to_keyed_nt_bt


def _read_cases(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    with path.open("r", encoding="utf-8") as f:
        filtered = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    reader = csv.DictReader(filtered)
    for r in reader:
        rows.append((r["cr_scale"].strip(), r["sh_scale"].strip(), r["z_metal"].strip()))
    return rows


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Scan etaA outliers at high density.")
    p.add_argument("--cases", required=True)
    p.add_argument("--base-dir", default=".dev/low-metal-eos/ref/Lmetal/data/table")
    p.add_argument("--target-dir", default="tools/eos_mhd_table_generator/results")
    p.add_argument("--out-csv", default="tools/eos_mhd_table_generator/results/etaa_outliers.csv")
    p.add_argument("--nt-min", type=float, default=19.0)
    p.add_argument("--abs-diff-min", type=float, default=3.0)
    p.add_argument("--top-k-per-case", type=int, default=20)
    return p


def main() -> int:
    args = _parser().parse_args()
    rows_out: list[list[str]] = []
    cols = ["case_request", "case_base", "case_target", "nt", "bt", "base_etaA", "target_etaA", "diff", "abs_diff"]

    for cr, sh, z in _read_cases(Path(args.cases)):
        req_case = f"CR{cr}_sh{sh}_Z{z}.d"
        base, base_name = resolve_case_file(Path(args.base_dir), cr, sh, z)
        target, target_name = resolve_case_file(Path(args.target_dir), cr, sh, z)
        if base is None or target is None:
            continue

        a = to_keyed_nt_bt(read_table(base, strict_cols=False))
        b = to_keyed_nt_bt(read_table(target, strict_cols=False))
        keys = sorted(set(a.keys()) & set(b.keys()))
        hit: list[list[str]] = []
        for nt, bt in keys:
            if nt < args.nt_min:
                continue
            ba = float(a[(nt, bt)][EOS_COLS["etaA"]])
            ta = float(b[(nt, bt)][EOS_COLS["etaA"]])
            d = ta - ba
            ad = abs(d)
            if ad >= args.abs_diff_min:
                hit.append(
                    [
                        req_case,
                        base_name,
                        target_name,
                        f"{nt:.6f}",
                        f"{bt:.6f}",
                        f"{ba:.6e}",
                        f"{ta:.6e}",
                        f"{d:.6e}",
                        f"{ad:.6e}",
                    ]
                )
        if hit:
            hit.sort(key=lambda r: float(r[-1]), reverse=True)
            rows_out.extend(hit[: max(args.top_k_per_case, 0)])

    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="ascii", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        w.writerows(rows_out)
    print(f"[OK] outliers_csv={out}")
    print(f"[OK] rows={len(rows_out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
