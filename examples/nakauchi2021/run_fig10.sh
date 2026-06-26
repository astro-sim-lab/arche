#!/usr/bin/env bash
# Copyright (C) 2026 Shingo Hirano and Sho Higashi
# Licensed under the MIT found in the
# https://github.com/astro-sim-lab/arche/blob/main/LICENSE
# run_fig10.sh — reproduce the simulation grid behind Nakauchi et al. (2021)
# Figure 10 with the arche-dev metal_grain network.
#
# Grid: 7 metallicities x 3 cosmic-ray ionization rates x {reduced, full}
#   Z/Zsun = 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6
#   zeta0  = 0, 1e-17, 1e-15  [s^-1]
#   net    = reduced metal Minimal (arche_collapse_metal_minimal)
#            full   Nakauchi2021  (metal_collapse)
# Each case is a one-zone gravitational collapse to nH = 1e23 cm^-3.
#
# Usage (from project root):
#   bash examples/nakauchi2021/run_fig10.sh              # build + run all
#   bash examples/nakauchi2021/run_fig10.sh --no-build   # reuse binaries
#
# Output: results/nakauchi2021/collapse_CR<zeta>_Z<z>[_min].h5
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${ROOT_DIR}"

BUILD_DIR="${BUILD_DIR:-build}"
OUT="${OUT:-results/nakauchi2021}"
BIN="${BUILD_DIR}/src/apps/collapse_metal_grain"
MAXJOBS="${MAXJOBS:-6}"

ZS="1 1e-1 1e-2 1e-3 1e-4 1e-5 1e-6"
ZETAS="0 1e-17 1e-15"

# ── build ─────────────────────────────────────────────────────────────────
if [[ "${1-}" != "--no-build" ]]; then
  cmake -S . -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17
  cmake --build "${BUILD_DIR}" \
    --target metal_collapse arche_collapse_metal_minimal -j"$(nproc)"
fi

mkdir -p "${OUT}"
export METAL_XNH_STOP=1e23 METAL_OUTPUT_STRIDE=20 METAL_OUTDIR="${OUT}"

run_one() {
  local z="$1" zeta="$2" net="$3" bin
  if [[ "${net}" == min ]]; then
    bin=arche_collapse_metal_minimal
  else
    bin=metal_collapse
  fi
  METAL_Z_METAL="${z}" METAL_ZETA0="${zeta}" "${BIN}/${bin}" \
    > "${OUT}/log_${net}_z${z}_cr${zeta}.log" 2>&1
}

for z in ${ZS}; do
  for zeta in ${ZETAS}; do
    for net in min full; do
      run_one "${z}" "${zeta}" "${net}" &
      while [[ "$(jobs -r | wc -l)" -ge "${MAXJOBS}" ]]; do wait -n; done
    done
  done
done
wait

echo "GRID DONE -> ${OUT}"
ls "${OUT}"/collapse_CR*_Z*.h5 | wc -l
