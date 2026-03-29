#!/usr/bin/env python3
"""
run_case_matrix.py - Batch runner for EOS table case generation.

Input format (CSV):
  cr_scale,sh_scale,z_metal
  1e-2,1e-2,1e-3
  1e0,1e0,1e-3

Lines starting with '#' are ignored.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path


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
    p = argparse.ArgumentParser(description="Run multiple CR/SH/Z cases via make_case_table.py.")
    p.add_argument("--cases", required=True, help="CSV file with columns: cr_scale,sh_scale,z_metal")
    p.add_argument("--output-dir", default="tools/eos_mhd_table_generator/results")
    p.add_argument("--build-dir", default="build-codex")
    p.add_argument("--mode", default="compat", choices=("compat", "normalized"))
    p.add_argument("--eta-model", default="legacy", choices=("placeholder", "legacy"))
    p.add_argument("--eta-resample", default="native", choices=("state", "native"))
    p.add_argument("--metal-max-iter", type=int)
    p.add_argument("--metal-xnh-stop")
    p.add_argument("--metal-cr-second-frac")
    p.add_argument("--metal-cr-metal-bkgnd")
    p.add_argument("--metal-abundance-set")
    p.add_argument("--metal-log-dir", help="Directory to store per-case metal_collapse logs.")
    p.add_argument("--skip-existing", action="store_true", help="Skip case if output .d already exists.")
    p.add_argument("--dry-run", action="store_true", help="Print commands without executing.")
    p.add_argument("--timeout-per-case", type=float, help="Timeout in seconds for each case process.")
    p.add_argument("--runtime-csv", default="tools/eos_mhd_table_generator/results/matrix_runtime.csv",
                   help="Per-case runtime summary CSV output path.")
    p.add_argument("--stop-on-error", action="store_true", help="Stop at first failure.")
    return p


def main() -> int:
    args = _parser().parse_args()
    cases = _read_cases(Path(args.cases))
    runner = Path("tools/eos_mhd_table_generator/python/app/make_case_table.py")
    if not runner.is_file():
        raise SystemExit(f"ERROR: not found: {runner}")

    failed = 0
    runtime_rows: list[list[str]] = []
    for i, (cr, sh, z) in enumerate(cases, start=1):
        out_name = f"CR{cr}_sh{sh}_Z{z}.d"
        out_path = Path(args.output_dir) / out_name
        if args.skip_existing and out_path.is_file():
            print(f"[{i}/{len(cases)}] {out_name} -> SKIP (exists)")
            continue
        cmd = [
            sys.executable,
            str(runner),
            "--build-dir",
            args.build_dir,
            "--cr-scale",
            cr,
            "--sh-scale",
            sh,
            "--z-metal",
            z,
            "--mode",
            args.mode,
            "--eta-model",
            args.eta_model,
            "--eta-resample",
            args.eta_resample,
            "--output-dir",
            args.output_dir,
            "--output-name",
            out_name,
        ]
        if args.metal_max_iter is not None:
            cmd.extend(["--metal-max-iter", str(args.metal_max_iter)])
        if args.metal_xnh_stop is not None:
            cmd.extend(["--metal-xnh-stop", str(args.metal_xnh_stop)])
        if args.metal_cr_second_frac is not None:
            cmd.extend(["--metal-cr-second-frac", str(args.metal_cr_second_frac)])
        if args.metal_cr_metal_bkgnd is not None:
            cmd.extend(["--metal-cr-metal-bkgnd", str(args.metal_cr_metal_bkgnd)])
        if args.metal_abundance_set is not None:
            cmd.extend(["--metal-abundance-set", str(args.metal_abundance_set)])
        if args.metal_log_dir is not None:
            log_dir = Path(args.metal_log_dir)
            log_dir.mkdir(parents=True, exist_ok=True)
            cmd.extend(["--metal-log", str(log_dir / f"{Path(out_name).stem}.log")])

        print(f"[{i}/{len(cases)}] {out_name}")
        if args.dry_run:
            print("  -> DRYRUN:", " ".join(cmd))
            continue
        t0 = time.perf_counter()
        proc = None
        try:
            proc = subprocess.run(cmd, timeout=args.timeout_per_case, capture_output=True, text=True)
            ret = proc.returncode
        except subprocess.TimeoutExpired:
            ret = 124
            print(f"  -> TIMEOUT ({args.timeout_per_case}s)")
        elapsed = time.perf_counter() - t0
        input_h5 = ""
        if proc is not None:
            if proc.stdout:
                print(proc.stdout, end="" if proc.stdout.endswith("\n") else "\n")
            if proc.stderr:
                print(proc.stderr, end="" if proc.stderr.endswith("\n") else "\n", file=sys.stderr)
            for ln in proc.stdout.splitlines():
                if ln.startswith("[ok] input_h5="):
                    input_h5 = ln.split("=", 1)[1].strip()
                    break
        runtime_rows.append(
            [
                out_name,
                cr,
                sh,
                z,
                f"{elapsed:.6f}",
                str(ret),
                input_h5,
            ]
        )
        if ret != 0:
            failed += 1
            print(f"  -> FAILED (exit={ret})")
            if args.stop_on_error:
                out_csv = Path(args.runtime_csv)
                out_csv.parent.mkdir(parents=True, exist_ok=True)
                with out_csv.open("w", encoding="ascii", newline="") as f:
                    w = csv.writer(f)
                    w.writerow(["case", "cr_scale", "sh_scale", "z_metal", "elapsed_sec", "exit_code", "input_h5"])
                    w.writerows(runtime_rows)
                return ret

    out_csv = Path(args.runtime_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="ascii", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "cr_scale", "sh_scale", "z_metal", "elapsed_sec", "exit_code", "input_h5"])
        w.writerows(runtime_rows)
    print(f"[OK] runtime_csv={out_csv}")

    if failed > 0:
        print(f"done with failures: {failed}/{len(cases)}")
        return 2
    print(f"done: {len(cases)} cases")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
