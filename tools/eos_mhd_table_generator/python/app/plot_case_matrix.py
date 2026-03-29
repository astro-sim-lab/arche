#!/usr/bin/env python3
"""
plot_case_matrix.py - Batch plot compare profiles for case matrix.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from datetime import datetime
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from lib.case_paths import resolve_case_file


def _read_cases(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    with path.open("r", encoding="utf-8") as f:
        filtered = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    reader = csv.DictReader(filtered)
    for r in reader:
        rows.append((r["cr_scale"].strip(), r["sh_scale"].strip(), r["z_metal"].strip()))
    return rows


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Batch-plot compare profiles for case matrix.")
    p.add_argument("--cases", required=True)
    p.add_argument("--base-dir", default=".dev/low-metal-eos/ref/Lmetal/data/table")
    p.add_argument("--target-dir", default="tools/eos_mhd_table_generator/results")
    p.add_argument("--out-dir", default="tools/eos_mhd_table_generator/results/plots/matrix")
    p.add_argument("--quantities", default="tt,mut,pt,etaO,etaA,etaH")
    p.add_argument("--bt-list", default="-20,-10,0,2")
    p.add_argument("--run-tag", help="Optional subdirectory name under --out-dir.")
    p.add_argument("--timestamp-tag", action="store_true", help="Use YYYYmmdd-HHMMSS as run tag.")
    return p


def main() -> int:
    args = _parser().parse_args()
    out_root = Path(args.out_dir)
    if args.timestamp_tag and args.run_tag:
        raise SystemExit("ERROR: use either --run-tag or --timestamp-tag")
    if args.timestamp_tag:
        out_root = out_root / datetime.now().strftime("%Y%m%d-%H%M%S")
    elif args.run_tag:
        out_root = out_root / args.run_tag
    out_root.mkdir(parents=True, exist_ok=True)
    plot_script = Path("tools/eos_mhd_table_generator/python/app/plot_compare_eos_table.py")

    n_ok = 0
    n_skip = 0
    for cr, sh, z in _read_cases(Path(args.cases)):
        req_case = f"CR{cr}_sh{sh}_Z{z}.d"
        base, _ = resolve_case_file(Path(args.base_dir), cr, sh, z)
        target, _ = resolve_case_file(Path(args.target_dir), cr, sh, z)
        if base is None or target is None:
            print(f"[SKIP] {req_case}: missing base/target")
            n_skip += 1
            continue
        case_out = out_root / req_case.removesuffix(".d")
        cmd = [
            sys.executable,
            str(plot_script),
            "--base",
            str(base),
            "--target",
            str(target),
            "--out-dir",
            str(case_out),
            f"--quantities={args.quantities}",
            f"--bt-list={args.bt_list}",
        ]
        ret = subprocess.run(cmd).returncode
        if ret == 0:
            print(f"[OK] {req_case}: {case_out}")
            n_ok += 1
        else:
            print(f"[FAIL] {req_case}: plot exit={ret}")
    print(f"[OK] plotted={n_ok} skipped={n_skip}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
