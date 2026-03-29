#!/usr/bin/env python3
"""
compare_state_native_matrix.py - Compare eta-resample=state vs native on same HDF5 per case.
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
from lib.compare_core import calc_metrics
from lib.config_kv import read_kv_config

# reuse existing builder CLI to keep behavior identical to normal pipeline
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


def _run_build(input_h5: Path, out_d: Path, eta_resample: str, cfg: Path) -> int:
    import subprocess

    cmd = [
        sys.executable,
        str(BUILD_CLI),
        "--config",
        str(cfg),
        "--input-h5",
        str(input_h5),
        "--output",
        str(out_d),
        "--eta-model",
        "legacy",
        "--eta-resample",
        eta_resample,
    ]
    return subprocess.run(cmd).returncode


def _has_nul_bytes(path: Path, *, probe_bytes: int = 2_000_000) -> bool:
    with path.open("rb") as f:
        chunk = f.read(probe_bytes)
    return b"\x00" in chunk


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Batch compare state/native resampling behavior by case.")
    p.add_argument("--cases", required=True)
    p.add_argument("--h5-root", default="tools/eos_mhd_table_generator/results/h5")
    p.add_argument("--base-dir", default=".dev/low-metal-eos/ref/Lmetal/data/table")
    p.add_argument("--config", default="tools/eos_mhd_table_generator/params/default.conf")
    p.add_argument("--work-dir", default="tools/eos_mhd_table_generator/results/state_native")
    p.add_argument("--out-csv", default="tools/eos_mhd_table_generator/results/state_native_compare.csv")
    p.add_argument(
        "--skip-build",
        action="store_true",
        help="Do not rebuild state/native .d files; compare files already in --work-dir.",
    )
    return p


def main() -> int:
    args = _parser().parse_args()
    cfg = Path(args.config)
    _ = read_kv_config(cfg)  # validate file format early
    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    out_rows: list[list[str]] = []

    for cr, sh, z in _read_cases(Path(args.cases)):
        req = f"CR{cr}_sh{sh}_Z{z}.d"
        h5 = _find_h5_for_case(Path(args.h5_root), cr, sh, z)
        base, base_name = resolve_case_file(Path(args.base_dir), cr, sh, z)
        if h5 is None or base is None:
            out_rows.append([req, base_name if base else "", "", "missing_input", "", "", "", "", "", ""])
            print(f"[MISS] {req}: h5/base missing")
            continue

        state_d = work / f"{req.removesuffix('.d')}.state.d"
        native_d = work / f"{req.removesuffix('.d')}.native.d"
        if args.skip_build:
            if not state_d.exists() or not native_d.exists():
                out_rows.append([req, base_name, str(h5), "missing_state_native_d", "", "", "", "", "", ""])
                print(f"[MISS] {req}: state/native .d missing in {work}")
                continue
            r1 = 0
            r2 = 0
        else:
            r1 = _run_build(h5, state_d, "state", cfg)
            r2 = _run_build(h5, native_d, "native", cfg)
            if r1 != 0 or r2 != 0:
                out_rows.append([req, base_name, str(h5), "build_fail", str(r1), str(r2), "", "", "", ""])
                print(f"[FAIL] {req}: build state/native failed ({r1},{r2})")
                continue

        if _has_nul_bytes(state_d) or _has_nul_bytes(native_d):
            out_rows.append([req, base_name, str(h5), "corrupt_output", str(r1), str(r2), "", "", "", ""])
            print(f"[FAIL] {req}: NUL bytes detected in state/native output")
            continue

        try:
            _, m_state = calc_metrics(base, state_d)
            _, m_native = calc_metrics(base, native_d)
        except Exception as exc:
            out_rows.append([req, base_name, str(h5), f"compare_fail:{type(exc).__name__}", str(r1), str(r2), "", "", "", ""])
            print(f"[FAIL] {req}: compare failed ({type(exc).__name__})")
            continue
        # focus on etaA and thermals
        out_rows.append(
            [
                req,
                base_name,
                str(h5),
                "ok",
                f"{m_state['etaA']['corr']:.6f}",
                f"{m_native['etaA']['corr']:.6f}",
                f"{m_state['etaA']['med_abs']:.6e}",
                f"{m_native['etaA']['med_abs']:.6e}",
                f"{(m_native['etaA']['med_abs'] - m_state['etaA']['med_abs']):.6e}",
                f"{(m_native['tt']['med_abs'] - m_state['tt']['med_abs']):.6e}",
            ]
        )
        print(f"[OK] {req}: etaA med_abs state={m_state['etaA']['med_abs']:.3e} native={m_native['etaA']['med_abs']:.3e}")

    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="ascii", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "case_request",
                "case_base",
                "input_h5",
                "status",
                "etaA_corr_state",
                "etaA_corr_native",
                "etaA_med_abs_state",
                "etaA_med_abs_native",
                "etaA_med_abs_native_minus_state",
                "tt_med_abs_native_minus_state",
            ]
        )
        w.writerows(out_rows)
    print(f"[OK] out_csv={out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
