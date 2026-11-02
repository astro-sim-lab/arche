#!/usr/bin/env bash
# Copyright (C) 2026 Shingo Hirano and Sho Higashi
# Licensed under the MIT found in the
# https://github.com/astro-sim-lab/arche/blob/main/LICENSE
# run_fig11.sh — magnetic Reynolds number maps of Nakauchi et al. (2021)
# Figure 11 with the full arche `metal_grain` network.
#
# Grid: 4 metallicities x 3 cosmic-ray ionization rates
#   Z/Zsun = 1e-6, 1e-4, 1e-2, 1
#   zeta0  = 0, 1e-17, 1e-15  [s^-1]
# Each case is a one-zone gravitational collapse to nH = 1e23 cm^-3
# (metal_collapse, the same settings as run_fig10.sh), followed by
# magnetic_reynolds, which evaluates Rm(lambda_J) on log10(B/G) = -10 ... 10.
#
# Usage (from project root):
#   bash examples/nakauchi2021/run_fig11.sh                 # build + run all
#   bash examples/nakauchi2021/run_fig11.sh --no-build      # reuse binaries
#   bash examples/nakauchi2021/run_fig11.sh --no-collapse   # reuse binaries and
#        the collapse files already in OUT (e.g. from run_fig10.sh)
#
# Output: results/nakauchi2021/collapse_CR<zeta>_Z<z>.h5  (collapse track)
#         results/nakauchi2021/rm_CR<zeta>_Z<z>.h5        (Rm map)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${ROOT_DIR}"

BUILD_DIR="${BUILD_DIR:-build}"
OUT="${OUT:-results/nakauchi2021}"
MAXJOBS="${MAXJOBS:-6}"
COLLAPSE="${BUILD_DIR}/src/apps/collapse_metal_grain/metal_collapse"
MAGRE="${BUILD_DIR}/src/apps/magnetic_reynolds/magnetic_reynolds"

# Written with an exponent: the collapse application then names its output
# by the value as given (Z1e-6 ... Z1e+0).
ZS="1e-6 1e-4 1e-2 1e+0"
ZETAS="0 1e-17 1e-15"

mode="${1-}"

# ── build ─────────────────────────────────────────────────────────────────
if [[ "${mode}" != "--no-build" && "${mode}" != "--no-collapse" ]]; then
  cmake -S . -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17
  cmake --build "${BUILD_DIR}" \
    --target metal_collapse magnetic_reynolds -j"$(nproc)"
fi

mkdir -p "${OUT}"
export METAL_XNH_STOP=1e23 METAL_OUTPUT_STRIDE=20 METAL_OUTDIR="${OUT}"

run_one() {
  local z="$1" zeta="$2" tag
  tag="CR${zeta}_Z${z}"
  if [[ "${mode}" != "--no-collapse" ]]; then
    METAL_Z_METAL="${z}" METAL_ZETA0="${zeta}" "${COLLAPSE}" \
      > "${OUT}/log_full_z${z}_cr${zeta}.log" 2>&1
  fi
  "${MAGRE}" "${OUT}/collapse_${tag}.h5" "${OUT}/rm_${tag}.h5" \
    > "${OUT}/log_rm_${tag}.log" 2>&1
}

for z in ${ZS}; do
  for zeta in ${ZETAS}; do
    run_one "${z}" "${zeta}" &
    while [[ "$(jobs -r | wc -l)" -ge "${MAXJOBS}" ]]; do wait -n; done
  done
done
wait

echo "GRID DONE -> ${OUT}"
ls "${OUT}"/rm_CR*_Z*.h5 | wc -l
