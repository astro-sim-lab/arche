#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DEFAULT_CONFIG="${SCRIPT_DIR}/params/default.conf"
DEFAULT_REF_TABLE=".dev/low-metal-eos/ref/Lmetal/data/table/CR1e-2_sh1e-2_Z1e-3.d"
DEFAULT_CASE_TABLE="tools/eos_mhd_table_generator/results/CR1e-2_sh1e-2_Z1e-3.d"
DEFAULT_PLOT_DIR="tools/eos_mhd_table_generator/results/plots/quickstart"

cmd="${1:-build}"
shift || true

case "${cmd}" in
  build)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/build_eos_table.py" --config "tools/eos_mhd_table_generator/params/default.conf" "$@")
    ;;
  case)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/make_case_table.py" "$@")
    ;;
  matrix)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/run_case_matrix.py" "$@")
    ;;
  matrix-compare)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/compare_case_matrix.py" "$@")
    ;;
  matrix-state-native)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/compare_state_native_matrix.py" "$@")
    ;;
  matrix-status)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/matrix_status.py" "$@")
    ;;
  matrix-report)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/report_compare_matrix.py" "$@")
    ;;
  matrix-sync)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/update_todo_from_eval.py" "$@")
    ;;
  matrix-plot)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/plot_case_matrix.py" "$@")
    ;;
  matrix-runtime-profile)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/profile_runtime_by_density.py" "$@")
    ;;
  matrix-build-runtime)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/profile_build_runtime_matrix.py" "$@")
    ;;
  matrix-etaa)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/scan_etaa_outliers.py" "$@")
    ;;
  matrix-etaa-onset)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/find_etaa_onset.py" "$@")
    ;;
  matrix-etaa-summary)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/summarize_etaa_outliers.py" "$@")
    ;;
  matrix-md)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/generate_matrix_report_md.py" "$@")
    ;;
  validate)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/validate_eos_table.py" --input "${DEFAULT_CASE_TABLE}" "$@")
    ;;
  selftest)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/test/selftest.py" "$@")
    ;;
  compare)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/compare_eos_table.py" \
      --base "${DEFAULT_REF_TABLE}" \
      --target "${DEFAULT_CASE_TABLE}" "$@")
    ;;
  plot)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/plot_compare_eos_table.py" \
      --base "${DEFAULT_REF_TABLE}" \
      --target "${DEFAULT_CASE_TABLE}" \
      --out-dir "${DEFAULT_PLOT_DIR}" "$@")
    ;;
  quickstart)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/test/selftest.py")
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/make_case_table.py" \
      --cr-scale 1e-2 --sh-scale 1e-2 --z-metal 1e-3 --eta-model legacy --eta-resample native)
    (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/validate_eos_table.py" --input "${DEFAULT_CASE_TABLE}")
    if [[ -f "${ROOT_DIR}/${DEFAULT_REF_TABLE}" ]]; then
      (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/compare_eos_table.py" \
        --base "${DEFAULT_REF_TABLE}" \
        --target "${DEFAULT_CASE_TABLE}")
      (cd "${ROOT_DIR}" && python3 "tools/eos_mhd_table_generator/python/app/plot_compare_eos_table.py" \
        --base "${DEFAULT_REF_TABLE}" \
        --target "${DEFAULT_CASE_TABLE}" \
        --out-dir "${DEFAULT_PLOT_DIR}")
      echo "[quickstart] plots: ${DEFAULT_PLOT_DIR}"
    else
      echo "[quickstart] reference table not found: ${DEFAULT_REF_TABLE}"
      echo "[quickstart] skip compare/plot (case + validate completed)"
    fi
    ;;
  *)
    cat <<EOF
Usage: $(basename "$0") [build|case|matrix|matrix-status|matrix-compare|matrix-state-native|matrix-report|matrix-sync|matrix-plot|matrix-runtime-profile|matrix-build-runtime|matrix-etaa|matrix-etaa-onset|matrix-etaa-summary|matrix-md|validate|compare|plot|selftest|quickstart] [extra args]

Examples:
  $(basename "$0") quickstart
  $(basename "$0") build
  $(basename "$0") case --cr-scale 1e-2 --sh-scale 1e-2 --z-metal 1e-3
  $(basename "$0") matrix --cases tools/eos_mhd_table_generator/params/cases_core.csv
  $(basename "$0") matrix --cases tools/eos_mhd_table_generator/params/cases_extended.csv --skip-existing --metal-log-dir tools/eos_mhd_table_generator/results/logs --timeout-per-case 7200
  $(basename "$0") matrix-status --cases tools/eos_mhd_table_generator/params/cases_extended.csv
  $(basename "$0") matrix-compare --cases tools/eos_mhd_table_generator/params/cases_core.csv
  $(basename "$0") matrix-state-native --cases tools/eos_mhd_table_generator/params/cases_core.csv
  $(basename "$0") matrix-report --summary-csv tools/eos_mhd_table_generator/results/compare_matrix_summary.csv --profile todo
  $(basename "$0") matrix-sync --eval-csv tools/eos_mhd_table_generator/results/compare_matrix_eval.csv
  $(basename "$0") matrix-plot --cases tools/eos_mhd_table_generator/params/cases_core.csv --timestamp-tag
  $(basename "$0") matrix-runtime-profile --runtime-csv tools/eos_mhd_table_generator/results/matrix_runtime.csv
  $(basename "$0") matrix-build-runtime --cases tools/eos_mhd_table_generator/params/cases_extended.csv --out-csv tools/eos_mhd_table_generator/results/matrix_runtime_build.csv
  $(basename "$0") matrix-etaa --cases tools/eos_mhd_table_generator/params/cases_core.csv
  $(basename "$0") matrix-etaa-onset --cases tools/eos_mhd_table_generator/params/cases_core.csv
  $(basename "$0") matrix-etaa-summary --in-csv tools/eos_mhd_table_generator/results/etaa_outliers.csv
  $(basename "$0") matrix-md --out-md .dev/low-metal-eos/matrix_report.md
  $(basename "$0") validate
  $(basename "$0") plot
  $(basename "$0") selftest
  $(basename "$0") compare
EOF
    exit 1
    ;;
esac
