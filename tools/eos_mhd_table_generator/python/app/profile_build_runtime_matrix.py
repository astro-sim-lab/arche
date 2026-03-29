#!/usr/bin/env python3
"""
profile_build_runtime_matrix.py - Measure per-case build_eos_table runtime from existing HDF5.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path


BUILD_CLI = Path("tools/eos_mhd_table_generator/python/app/build_eos_table.py")


def _read_cases(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    with path.open("r", encoding="utf-8") as f:
        filtered = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    reader = csv.DictReader(filtered)
    for r in reader:
        rows.append((r["cr_scale"].strip(), r["sh_scale"].strip(), r["z_metal"].strip()))
    return rows


def _find_h5_for_case(h5_root: Path, cr: str, sh: str, z: str) -> Path | None:
    cand = sorted(h5_root.glob(f"cr{cr}_sh{sh}_*_z{z}/collapse_*.h5"))
    if not cand and z in ("0", "0.0"):
        cand = sorted(h5_root.glob(f"cr{cr}_sh{sh}_*_z1e-inf/collapse_*.h5"))
    if not cand:
        return None
    return max(cand, key=lambda p: p.stat().st_mtime)


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Profile HDF5->.d build runtime for case matrix.")
    p.add_argument("--cases", required=True)
    p.add_argument("--config", default="tools/eos_mhd_table_generator/params/default.conf")
    p.add_argument("--h5-root", default="tools/eos_mhd_table_generator/results/h5")
    p.add_argument("--work-dir", default="tools/eos_mhd_table_generator/results/runtime_build")
    p.add_argument("--eta-model", default="legacy", choices=("placeholder", "legacy"))
    p.add_argument("--eta-resample", default="native", choices=("state", "native"))
    p.add_argument("--timeout-per-case", type=float)
    p.add_argument("--out-csv", default="tools/eos_mhd_table_generator/results/matrix_runtime_build.csv")
    return p


def main() -> int:
    args = _parser().parse_args()
    cases = _read_cases(Path(args.cases))
    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)

    rows: list[list[str]] = []
    failed = 0
    for i, (cr, sh, z) in enumerate(cases, start=1):
        case_name = f"CR{cr}_sh{sh}_Z{z}.d"
        h5 = _find_h5_for_case(Path(args.h5_root), cr, sh, z)
        if h5 is None:
            rows.append([case_name, cr, sh, z, "", "2", "", "missing_h5"])
            print(f"[MISS] {case_name}: h5 not found")
            failed += 1
            continue

        out_d = work / f"{case_name.removesuffix('.d')}.{args.eta_resample}.d"
        cmd = [
            sys.executable,
            str(BUILD_CLI),
            "--config",
            args.config,
            "--input-h5",
            str(h5),
            "--output",
            str(out_d),
            "--eta-model",
            args.eta_model,
            "--eta-resample",
            args.eta_resample,
        ]
        print(f"[{i}/{len(cases)}] {case_name}")
        t0 = time.perf_counter()
        try:
            proc = subprocess.run(cmd, timeout=args.timeout_per_case, capture_output=True, text=True)
            ret = proc.returncode
        except subprocess.TimeoutExpired:
            ret = 124
        elapsed = time.perf_counter() - t0
        rows.append([case_name, cr, sh, z, f"{elapsed:.6f}", str(ret), str(h5), str(out_d)])
        if ret != 0:
            failed += 1
            print(f"  -> FAIL (exit={ret})")
        else:
            print(f"  -> OK ({elapsed:.3f}s)")

    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="ascii", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "cr_scale", "sh_scale", "z_metal", "elapsed_sec", "exit_code", "input_h5", "output_d"])
        w.writerows(rows)
    print(f"[OK] runtime_csv={out}")
    if failed > 0:
        print(f"[WARN] failures={failed}/{len(cases)}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
