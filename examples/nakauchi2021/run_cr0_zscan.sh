#!/usr/bin/env bash
# run_cr0_zscan.sh
#
# Run metal_grain collapse for CR=0 across 7 metallicities.
# Results go to results/cr0_zscan/metal_grain/collapse_CR0_Z<z>.h5
#
# Usage (from project root):
#   bash examples/nakauchi2021/run_cr0_zscan.sh
#   bash examples/nakauchi2021/run_cr0_zscan.sh --no-build

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${ROOT_DIR}"

NO_BUILD=""
[[ "${1-}" == "--no-build" ]] && NO_BUILD="--no-build"

OUT="results/cr0_zscan"

Z_VALUES=(1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1)

for Z in "${Z_VALUES[@]}"; do
    echo ">>> Z_metal = ${Z}"
    bash run_collapse.sh \
        ${NO_BUILD} \
        --no-prim \
        --metal-zeta0 0 \
        --metal-z-metal "${Z}" \
        --no-plot \
        --no-resample \
        --out-dir  "${OUT}" \
        --save-dir "${OUT}"
    NO_BUILD="--no-build"
done

echo ""
echo "Done. HDF5 files:"
ls "${OUT}/metal_grain/"collapse_CR0_Z*.h5
