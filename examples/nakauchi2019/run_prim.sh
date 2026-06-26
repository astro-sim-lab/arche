#!/usr/bin/env bash
# Copyright (C) 2026 Shingo Hirano and Sho Higashi
# Licensed under the MIT found in the
# https://github.com/astro-sim-lab/arche/blob/main/LICENSE
# run_prim.sh — primordial (zero-metal) collapse grid for the reduced-vs-full
# comparison figure, in the spirit of Nakauchi, Omukai & Susa (2019,
# MNRAS 488, 1846).
#
# Grid: 3 cosmic-ray ionization rates x {reduced, full}
#   zeta0 = 0, 1e-17, 1e-15  [s^-1]
#   net   = reduced Nakauchi2019_Minimal (arche_collapse_prim_minimal, 15 sp)
#           full   Nakauchi2019         (prim_collapse, 23 sp)
# Each case is a one-zone gravitational collapse to nH = 1e23 cm^-3.
# (Primordial chemistry has no metallicity axis, so unlike the metal grid the
#  curves are split by chemical species rather than by Z.)
#
# Usage (from project root):
#   bash examples/nakauchi2019/run_prim.sh              # build + run all
#   bash examples/nakauchi2019/run_prim.sh --no-build   # reuse binaries
#
# Output: results/nakauchi2019/collapse_CR<zeta>[_min].h5
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${ROOT_DIR}"

BUILD_DIR="${BUILD_DIR:-build}"
OUT="${OUT:-results/nakauchi2019}"
BIN="${BUILD_DIR}/src/apps/collapse_primordial"
MAXJOBS="${MAXJOBS:-6}"

ZETAS="0 1e-17 1e-15"

# ── build ─────────────────────────────────────────────────────────────────
if [[ "${1-}" != "--no-build" ]]; then
  cmake -S . -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17
  cmake --build "${BUILD_DIR}" \
    --target prim_collapse arche_collapse_prim_minimal -j"$(nproc)"
fi

mkdir -p "${OUT}"
export PRIM_XNH_STOP=1e23 PRIM_OUTPUT_STRIDE=20 PRIM_OUTDIR="${OUT}"

run_one() {
  local zeta="$1" net="$2" bin
  if [[ "${net}" == min ]]; then
    bin=arche_collapse_prim_minimal
  else
    bin=prim_collapse
  fi
  PRIM_ZETA0="${zeta}" "${BIN}/${bin}" \
    > "${OUT}/log_${net}_cr${zeta}.log" 2>&1
}

for zeta in ${ZETAS}; do
  for net in min full; do
    run_one "${zeta}" "${net}" &
    while [[ "$(jobs -r | wc -l)" -ge "${MAXJOBS}" ]]; do wait -n; done
  done
done
wait

echo "GRID DONE -> ${OUT}"
ls "${OUT}"/collapse_CR*.h5 | wc -l
