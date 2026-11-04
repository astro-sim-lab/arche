#!/usr/bin/env bash
# Copyright (C) 2026 Shingo Hirano and Sho Higashi
# Licensed under the MIT found in the
# https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#
# Build ARCHE's collapse applications and run the collapse test cells
# (tests/collapse/), as the CI does, on this machine.
#
#   tools/ci/run_local.sh [--compiler gcc|icx] [--march VALUE] [--subset G|S]
#                         [--shard I/N] [--jobs N] [--build-dir DIR] [--out-dir DIR]
#                         [--no-build] [--title TEXT]
#
# Defaults: --compiler gcc --march x86-64-v3 --subset G --jobs $(nproc)
#           --build-dir build_test_<compiler>_<march> --out-dir results/tests/<compiler>_<march>_<subset>
# --march '' builds without -march (FMA off on x86-64).  For icx, source
# /opt/intel/oneapi/setvars.sh first.  Run from anywhere; paths are relative to
# the repository root.  Exit status 1 if any cell fails.  When
# GITHUB_STEP_SUMMARY is set (GitHub Actions), the report is appended to it.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

COMPILER=gcc
MARCH=x86-64-v3
SUBSET=G
SHARD=1/1
JOBS="$(nproc)"
TITLE=""
BUILD_DIR=""
OUT_DIR=""
NO_BUILD=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --compiler) COMPILER="$2"; shift 2 ;;
    --march) MARCH="$2"; shift 2 ;;
    --subset) SUBSET="$2"; shift 2 ;;
    --shard) SHARD="$2"; shift 2 ;;
    --title) TITLE="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --build-dir) BUILD_DIR="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --no-build) NO_BUILD=1; shift ;;
    -h|--help) sed -n '5,19p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

case "$COMPILER" in
  gcc) CXX=g++; CC=gcc ;;
  icx) CXX=icpx; CC=icx ;;
  *) echo "--compiler must be gcc or icx" >&2; exit 2 ;;
esac
TAG="${COMPILER}_${MARCH:-nomarch}"
BUILD_DIR="${BUILD_DIR:-build_test_${TAG}}"
OUT_DIR="${OUT_DIR:-results/tests/${TAG}_${SUBSET}}"

if [[ "$NO_BUILD" -eq 0 ]]; then
  cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_C_COMPILER="$CC" -DARCHE_MARCH="$MARCH"
  cmake --build "$BUILD_DIR" -j "$JOBS" --target prim_collapse arche_collapse_prim_minimal \
    metal_collapse arche_collapse_metal_minimal
fi

python3 tests/collapse/run_cells.py --build-dir "$BUILD_DIR" --out-dir "$OUT_DIR" \
  --subset "$SUBSET" --shard "$SHARD" --jobs "$JOBS"
SUMMARY=()
if [[ -n "${GITHUB_STEP_SUMMARY:-}" ]]; then SUMMARY=(--summary "$GITHUB_STEP_SUMMARY"); fi
python3 tests/collapse/check_runs.py "$OUT_DIR" "${SUMMARY[@]}" \
  --title "${TITLE:-collapse cells ($TAG, $SUBSET, shard $SHARD)}"
