#!/usr/bin/env python3
# Copyright (C) 2026 Shingo Hirano and Sho Higashi
# Licensed under the MIT found in the
# https://github.com/astro-sim-lab/arche/blob/main/LICENSE
"""Run collapse test cells and record what happened (no pass/fail decision).

Each cell runs in its own output directory under --out-dir, from the repository
root (the applications read data/*.h5 by relative path), with
OMP_NUM_THREADS=1.  A cell that exceeds its wall-time gate is killed.

The gate is defined on a reference machine (--gate-sec, default 300 s) and
scaled by the speed of the machine running the tests: before the cells, the
calibration cell (pf, zeta 0, f_ret 1, z 20) runs --jobs times concurrently and
its median wall time is divided by CALIB_REF_SEC, the time of the same
measurement on the reference machine.  The scale is never below 1.

Output: <out-dir>/results.json (per-cell record) and, per cell,
<out-dir>/<cell>/trace.npz (log10 n_H and log10 T of every HDF5 row) and a
trimmed stdout log.  check_runs.py applies the criteria.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import platform
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import h5py  # type: ignore[import-untyped]  # h5py ships no type stubs
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cells import BINARIES, Cell, grid, select  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]

# Reference machine: Intel Core i5-12600KF (6 P-cores), WSL2, g++ 11.4.0 Release
# -march=x86-64-v3, calibration cell run 1-4 times concurrently: median 4.2-4.4 s
# (2026-11-03).  At 8 and 16 concurrent runs the median grows to 5.2 s and 7.4 s;
# calibrating at the same concurrency as the cells carries that into the scale.
CALIB_REF_SEC = 4.3

EXIT_RE = re.compile(r"=== EXIT code=(\d+): (.*) ===")
CONS_RE = re.compile(r"=== CONSERVATION rows=\d+ projected=(\d+)/(\d+)")


def _calibration_cell() -> Cell:
    for c in grid():
        if c.network == "pf" and c.group == "G1" and c.zeta == 0.0 and c.f_ret == 1.0:
            return c
    raise RuntimeError("calibration cell not found in the grid")


def _trim_log(path: Path, head: int = 20, tail: int = 40) -> None:
    """Keep the first and last lines and every '===', WARNING or ERROR line."""
    lines = path.read_text(errors="replace").splitlines()
    if len(lines) <= head + tail:
        return
    keep_mid = [ln for ln in lines[head:-tail] if ln.startswith("===") or "WARNING" in ln or "ERROR" in ln]
    out = lines[:head] + [f"... ({len(lines) - head - tail} lines trimmed) ..."] + keep_mid + lines[-tail:]
    path.write_text("\n".join(out) + "\n")


def _read_h5(cell_dir: Path) -> dict[str, Any]:
    files = sorted(cell_dir.glob("*.h5"))
    if len(files) != 1:
        return {"h5": None, "h5_error": f"{len(files)} HDF5 files"}
    rec: dict[str, Any] = {"h5": str(files[0])}
    try:
        with h5py.File(files[0], "r") as h:
            nh = np.asarray(h["nH"][:], dtype=np.float64)
            t = np.asarray(h["T_K"][:], dtype=np.float64)
            ci = np.asarray(h["charge_imbal"][:], dtype=np.float64)
            msg = h.attrs.get("exit_message", b"")
            rec.update(
                exit_code_h5=int(h.attrs["exit_code"]) if "exit_code" in h.attrs else None,
                exit_message=msg.decode() if isinstance(msg, bytes) else str(msg),
                n_rows=int(nh.size),
                last_nH=float(nh[-1]) if nh.size else None,
                max_abs_charge_imbal=float(np.max(np.abs(ci))) if ci.size else None,
                finite=bool(np.all(np.isfinite(nh)) and np.all(np.isfinite(t))),
            )
            with np.errstate(divide="ignore", invalid="ignore"):
                np.savez_compressed(cell_dir / "trace.npz", log_nH=np.log10(nh), log_T=np.log10(t))
    except (OSError, KeyError) as exc:
        rec["h5_error"] = f"{type(exc).__name__}: {exc}"
    return rec


def run_cell(cell: Cell, build_dir: Path, out_dir: Path, timeout: float, nh_stop: float) -> dict[str, Any]:
    cell_dir = out_dir / cell.name
    cell_dir.mkdir(parents=True, exist_ok=True)
    for old in cell_dir.glob("*.h5"):
        old.unlink()  # output of an earlier run of this script
    exe = build_dir / "src" / "apps" / BINARIES[cell.network]
    env = dict(os.environ)
    env.update(cell.env(str(cell_dir), nh_stop))
    env["OMP_NUM_THREADS"] = "1"
    log = cell_dir / "stdout.log"
    rec: dict[str, Any] = {"name": cell.name, **cell.params(), "timeout_sec": timeout}
    t0 = time.monotonic()
    with log.open("w") as fout:
        try:
            proc = subprocess.run(
                [str(exe)], cwd=REPO_ROOT, env=env, stdout=fout, stderr=subprocess.STDOUT, timeout=timeout
            )
            rec["status"] = "finished"
            rec["returncode"] = proc.returncode
        except subprocess.TimeoutExpired:
            rec["status"] = "timeout"
            rec["returncode"] = None
    rec["wall_sec"] = time.monotonic() - t0
    text = log.read_text(errors="replace")
    m = EXIT_RE.findall(text)
    rec["exit_code_stdout"] = int(m[-1][0]) if m else None
    c = CONS_RE.findall(text)
    rec["projected"] = [int(c[-1][0]), int(c[-1][1])] if c else None
    _trim_log(log)
    rec["log"] = str(log)
    rec.update(_read_h5(cell_dir))
    return rec


def calibrate(build_dir: Path, out_dir: Path, jobs: int, nh_stop: float) -> list[float]:
    cell = _calibration_cell()
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = [
            ex.submit(run_cell, cell, build_dir, out_dir / f"_calibration_{i}", 600.0, nh_stop) for i in range(jobs)
        ]
        return [float(f.result()["wall_sec"]) for f in futs]


def _cmake_cache(build_dir: Path) -> dict[str, str]:
    keys = ("CMAKE_CXX_COMPILER", "CMAKE_BUILD_TYPE", "ARCHE_MARCH", "ARCHE_XRAY")
    out: dict[str, str] = {}
    cache = build_dir / "CMakeCache.txt"
    if cache.is_file():
        for ln in cache.read_text(errors="replace").splitlines():
            for k in keys:
                if ln.startswith(k + ":"):
                    out[k] = ln.split("=", 1)[1]
    return out


def _cpu_model() -> str:
    try:
        for ln in Path("/proc/cpuinfo").read_text().splitlines():
            if ln.startswith("model name"):
                return ln.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor()


def _git_commit() -> str:
    try:
        r = subprocess.run(
            ["git", "describe", "--always", "--dirty"], cwd=REPO_ROOT, capture_output=True, text=True, check=True
        )
        return r.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--subset", choices=("G", "S"), default="G")
    ap.add_argument("--shard", default="1/1", help="i/n: run every n-th cell starting at the i-th (1-based)")
    ap.add_argument("--only", default=None, help="regular expression on cell names")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--gate-sec", type=float, default=300.0, help="wall-time gate on the reference machine")
    ap.add_argument("--warn-sec", type=float, default=120.0, help="report threshold on the reference machine")
    ap.add_argument("--nh-stop", type=float, default=1e23)
    ap.add_argument("--no-calibrate", action="store_true", help="use scale 1 (no calibration run)")
    ap.add_argument("--list", action="store_true", help="print the selected cells and exit")
    args = ap.parse_args()

    i, n = (int(v) for v in args.shard.split("/"))
    if not 1 <= i <= n:
        ap.error(f"bad --shard {args.shard}")
    cells = select(args.subset)[i - 1 :: n]
    if args.only:
        cells = [c for c in cells if re.search(args.only, c.name)]
    if args.list:
        for c in cells:
            print(c.name)
        return 0

    build_dir = args.build_dir.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    missing = sorted(
        {str(build_dir / "src" / "apps" / BINARIES[c.network]) for c in cells}
        - {str(p) for p in (build_dir / "src" / "apps").glob("*/*") if p.is_file()}
    )
    if missing:
        print("missing executables:\n  " + "\n  ".join(missing), file=sys.stderr)
        return 2

    calib: list[float] = []
    scale = 1.0
    if not args.no_calibrate:
        calib = calibrate(build_dir, out_dir, args.jobs, args.nh_stop)
        scale = max(1.0, statistics.median(calib) / CALIB_REF_SEC)
        print(
            f"calibration: {len(calib)} concurrent runs, median {statistics.median(calib):.2f} s "
            f"(reference {CALIB_REF_SEC} s) -> scale {scale:.2f}",
            flush=True,
        )
    gate = args.gate_sec * scale
    warn = args.warn_sec * scale

    meta = {
        "commit": _git_commit(),
        "cmake": _cmake_cache(build_dir),
        "cpu": _cpu_model(),
        "host": platform.node(),
        "subset": args.subset,
        "shard": args.shard,
        "jobs": args.jobs,
        "nh_stop": args.nh_stop,
        "calibration_sec": calib,
        "calibration_ref_sec": CALIB_REF_SEC,
        "scale": scale,
        "gate_sec": gate,
        "warn_sec": warn,
        "started": time.strftime("%Y-%m-%d %H:%M:%S %z"),
    }
    print(f"{len(cells)} cells, {args.jobs} concurrent, gate {gate:.0f} s, report above {warn:.0f} s", flush=True)

    records: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(run_cell, c, build_dir, out_dir, gate, args.nh_stop): c for c in cells}
        for k, f in enumerate(concurrent.futures.as_completed(futs), 1):
            r = f.result()
            records.append(r)
            print(
                f"[{k:3d}/{len(cells)}] {r['name']:<40s} {r['status']:<8s} "
                f"exit={r.get('exit_code_h5')} {r['wall_sec']:7.1f} s",
                flush=True,
            )
    order = {c.name: k for k, c in enumerate(cells)}
    records.sort(key=lambda r: order[str(r["name"])])
    meta["finished"] = time.strftime("%Y-%m-%d %H:%M:%S %z")
    (out_dir / "results.json").write_text(json.dumps({"meta": meta, "cells": records}, indent=1) + "\n")
    print(f"wrote {out_dir / 'results.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
