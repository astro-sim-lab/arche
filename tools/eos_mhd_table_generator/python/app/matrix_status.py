#!/usr/bin/env python3
"""
matrix_status.py - Show completion status for case matrix.
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


def _read_cases(path: Path) -> list[tuple[str, str, str]]:
    out: list[tuple[str, str, str]] = []
    with path.open("r", encoding="utf-8") as f:
        filtered = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    reader = csv.DictReader(filtered)
    for r in reader:
        out.append((r["cr_scale"].strip(), r["sh_scale"].strip(), r["z_metal"].strip()))
    return out


def _find_h5(h5_root: Path, cr: str, sh: str, z: str) -> Path | None:
    cand = sorted(h5_root.glob(f"cr{cr}_sh{sh}_*_z{z}/collapse_*.h5"))
    if not cand and z in ("0", "0.0"):
        cand = sorted(h5_root.glob(f"cr{cr}_sh{sh}_*_z1e-inf/collapse_*.h5"))
    if not cand:
        return None
    return max(cand, key=lambda p: p.stat().st_mtime)


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Show matrix completion status.")
    p.add_argument("--cases", required=True)
    p.add_argument("--target-dir", default="tools/eos_mhd_table_generator/results")
    p.add_argument("--base-dir", default=".dev/low-metal-eos/ref/Lmetal/data/table")
    p.add_argument("--h5-root", default="tools/eos_mhd_table_generator/results/h5")
    p.add_argument("--out-csv", default="tools/eos_mhd_table_generator/results/matrix_status.csv")
    return p


def main() -> int:
    args = _parser().parse_args()
    rows: list[list[str]] = []
    n_done = 0
    n_missing_d = 0
    n_missing_h5 = 0
    n_missing_base = 0
    for cr, sh, z in _read_cases(Path(args.cases)):
        req = f"CR{cr}_sh{sh}_Z{z}.d"
        target, target_name = resolve_case_file(Path(args.target_dir), cr, sh, z)
        base, base_name = resolve_case_file(Path(args.base_dir), cr, sh, z)
        h5 = _find_h5(Path(args.h5_root), cr, sh, z)
        status = "ready_compare"
        if target is None:
            status = "missing_d"
            n_missing_d += 1
        if h5 is None:
            status = "missing_h5" if status == "ready_compare" else f"{status}+missing_h5"
            n_missing_h5 += 1
        if base is None:
            status = "missing_base" if status == "ready_compare" else f"{status}+missing_base"
            n_missing_base += 1
        if status == "ready_compare":
            n_done += 1

        rows.append(
            [
                req,
                target_name if target else "",
                base_name if base else "",
                str(h5) if h5 else "",
                status,
            ]
        )
        print(f"[{status}] {req}")

    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="ascii", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case_request", "case_target", "case_base", "h5_path", "status"])
        w.writerows(rows)
    print(f"[OK] out_csv={out}")
    print(f"[OK] counts: ready_compare={n_done} missing_d={n_missing_d} missing_h5={n_missing_h5} missing_base={n_missing_base}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

