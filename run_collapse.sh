#!/usr/bin/env bash
# run_collapse.sh — build, simulate, visualize, and resample in one step
#
# Usage:
#   bash run_collapse.sh [options]
#
# Options:
#   Shared parameters (apply to both primordial and metal_grain):
#     --common-zeta0         VALUE   CR ionization rate [s^-1]
#     --common-ff-ret        VALUE   scalar free-fall retardation factor (>0)
#     --common-fret-table    FILE    2-col step-function table: nH[cm^-3] f_ret
#     --common-jlw21         VALUE   Lyman-Werner intensity J_21 [10^-21 erg/s/cm^2/Hz/sr]
#     --common-redshift      VALUE   cosmological redshift z (T_rad = 2.725 * (1+z) [K])
#     --common-tk0           VALUE   initial gas temperature [K]
#     --common-ye0           VALUE   initial electron / H+ fraction
#     --common-yh2           VALUE   initial H2 fraction
#     --common-yhd           VALUE   initial HD fraction
#     --common-abundance-set VALUE   abundance preset name (solar/default/alpha-enhanced)
#     --common-xnh0          VALUE   initial H number density [cm^-3]
#     --common-output-stride VALUE   HDF5 output stride
#     --common-max-iter      VALUE   maximum integration steps
#     --common-dt-factor     VALUE   timestep factor after initial steps (>0)
#     --common-dt-factor-init VALUE  timestep factor during initial steps (>0)
#     --common-n-init-steps  VALUE   number of initial short-timestep steps (>=0)
#     --common-xnh-stop      VALUE   stop threshold for nH [cm^-3] (>0)
#     --common-cr-col-dens   VALUE   CR attenuation column density [g cm^-2] (>0)
#     --common-data-dir      DIR     reaction data directory for both runs
#   Per-network override:
#     Replace `common` with `prim` or `metal` to override one side only.
#     Example: --common-zeta0 1e-17 --metal-zeta0 1e-16
#   Network-only parameter:
#     --metal-z-metal        VALUE   metallicity for metal_grain [Z_sun]
#     --metal-sra-rate       VALUE   short-lived radionuclide scale (>=0)
#     --metal-lra-rate       VALUE   long-lived radionuclide scale (>=0)
#   --build-dir  DIR          CMake build directory (default: build)
#   --out-dir    DIR          root directory for HDF5 output (default: results/h5)
#   --save-dir   DIR          root directory for figures (default: results/fig)
#   --config     FILE         parameter file path (default: params/default.conf; required unless --no-config)
#   --no-config               disable parameter file loading
#   --no-build                skip cmake/build step (use existing binaries)
#   --no-prim                 skip primordial_collapse
#   --no-metal                skip metal_grain_collapse
#   --no-plot                 skip analyze_collapse.py
#   --fig-combo               also output fig_summary_<tag>.png
#   --no-resample             skip resample_collapse.py (step 5)
#   --resample                run resample_collapse.py even when --no-plot is set (default: on)
#
# Config/env variable details: docs/parameters.md and params/default.conf
#
# Examples:
#   bash run_collapse.sh --prim-zeta0 0 --metal-zeta0 1e-17 --metal-z-metal 1e-3
#   bash run_collapse.sh --no-build --no-prim --metal-zeta0 1e-16 --metal-z-metal 1.2e-4
#   bash run_collapse.sh --no-metal --prim-zeta0 1.5e-18
#   bash run_collapse.sh --prim-zeta0 1e-17 --prim-ff-ret 3.0 \
#                              --metal-zeta0 1e-17 --metal-z-metal 1e-3 --metal-ff-ret 3.0
#   bash run_collapse.sh --prim-zeta0 1e-17 --prim-redshift 20 \
#                              --metal-zeta0 1e-17 --metal-z-metal 1e-3 --metal-redshift 20
#   bash run_collapse.sh --prim-zeta0 1e-17 --prim-fret-table data/fret_table/fret_step_sample.dat \
#                              --metal-zeta0 1e-17 --metal-z-metal 1e-3
#   bash run_collapse.sh --prim-zeta0 1e-17 --no-plot   # simulate + resample only

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# ── Environment snapshot (for precedence merge) ──────────────────────────────
ENV_PRIM_ZETA0="${PRIM_ZETA0-}"
ENV_METAL_ZETA0="${METAL_ZETA0-}"
ENV_METAL_Z_METAL="${METAL_Z_METAL-}"
ENV_PRIM_FF_RET="${PRIM_FF_RET-}"
ENV_PRIM_FRET_TABLE="${PRIM_FRET_TABLE-}"
ENV_METAL_FF_RET="${METAL_FF_RET-}"
ENV_METAL_FRET_TABLE="${METAL_FRET_TABLE-}"
ENV_PRIM_JLW21="${PRIM_JLW21-}"
ENV_METAL_JLW21="${METAL_JLW21-}"
ENV_PRIM_REDSHIFT="${PRIM_REDSHIFT-}"
ENV_METAL_REDSHIFT="${METAL_REDSHIFT-}"
ENV_PRIM_TK0="${PRIM_TK0-}"
ENV_PRIM_YE0="${PRIM_YE0-}"
ENV_PRIM_YH2="${PRIM_YH2-}"
ENV_PRIM_YHD="${PRIM_YHD-}"
ENV_PRIM_ABUNDANCE_SET="${PRIM_ABUNDANCE_SET-}"
ENV_METAL_TK0="${METAL_TK0-}"
ENV_METAL_YE0="${METAL_YE0-}"
ENV_METAL_YH2="${METAL_YH2-}"
ENV_METAL_YHD="${METAL_YHD-}"
ENV_METAL_ABUNDANCE_SET="${METAL_ABUNDANCE_SET-}"
ENV_PRIM_XNH0="${PRIM_XNH0-}"
ENV_PRIM_OUTPUT_STRIDE="${PRIM_OUTPUT_STRIDE-}"
ENV_PRIM_MAX_ITER="${PRIM_MAX_ITER-}"
ENV_PRIM_DT_FACTOR="${PRIM_DT_FACTOR-}"
ENV_PRIM_DT_FACTOR_INIT="${PRIM_DT_FACTOR_INIT-}"
ENV_PRIM_N_INIT_STEPS="${PRIM_N_INIT_STEPS-}"
ENV_PRIM_XNH_STOP="${PRIM_XNH_STOP-}"
ENV_PRIM_CR_ATTEN_COL_DENS="${PRIM_CR_ATTEN_COL_DENS-}"
ENV_PRIM_DATA_DIR="${PRIM_DATA_DIR-}"
ENV_METAL_XNH0="${METAL_XNH0-}"
ENV_METAL_OUTPUT_STRIDE="${METAL_OUTPUT_STRIDE-}"
ENV_METAL_MAX_ITER="${METAL_MAX_ITER-}"
ENV_METAL_DT_FACTOR="${METAL_DT_FACTOR-}"
ENV_METAL_DT_FACTOR_INIT="${METAL_DT_FACTOR_INIT-}"
ENV_METAL_N_INIT_STEPS="${METAL_N_INIT_STEPS-}"
ENV_METAL_XNH_STOP="${METAL_XNH_STOP-}"
ENV_METAL_CR_ATTEN_COL_DENS="${METAL_CR_ATTEN_COL_DENS-}"
ENV_METAL_CR_ATTEN_SECOND_FRAC="${METAL_CR_ATTEN_SECOND_FRAC-}"
ENV_METAL_CR_METAL_BKGND="${METAL_CR_METAL_BKGND-}"
ENV_METAL_SRA_RATE="${METAL_SRA_RATE-}"
ENV_METAL_LRA_RATE="${METAL_LRA_RATE-}"
ENV_METAL_T_CR_DES="${METAL_T_CR_DES-}"
ENV_METAL_C_GAS_FRAC="${METAL_C_GAS_FRAC-}"
ENV_METAL_O_GAS_FRAC="${METAL_O_GAS_FRAC-}"
ENV_METAL_MG_GAS_FRAC="${METAL_MG_GAS_FRAC-}"
ENV_METAL_DATA_DIR="${METAL_DATA_DIR-}"
ENV_BUILD_DIR="${BUILD_DIR-}"
ENV_OUT_DIR="${OUT_DIR-}"
ENV_SAVE_DIR="${SAVE_DIR-}"

# ── Defaults ─────────────────────────────────────────────────────────────────
PRIM_ZETA0=""
METAL_ZETA0=""
METAL_Z_METAL=""
COMMON_ZETA0=""
COMMON_FF_RET=""
COMMON_FRET_TABLE=""
COMMON_JLW21=""
COMMON_REDSHIFT=""
COMMON_TK0=""
COMMON_YE0=""
COMMON_YH2=""
COMMON_YHD=""
COMMON_ABUNDANCE_SET=""
COMMON_XNH0=""
COMMON_OUTPUT_STRIDE=""
COMMON_MAX_ITER=""
COMMON_DT_FACTOR=""
COMMON_DT_FACTOR_INIT=""
COMMON_N_INIT_STEPS=""
COMMON_XNH_STOP=""
COMMON_CR_ATTEN_COL_DENS=""
COMMON_DATA_DIR=""
PRIM_FF_RET=""
PRIM_FRET_TABLE=""
METAL_FF_RET=""
METAL_FRET_TABLE=""
PRIM_JLW21=""
METAL_JLW21=""
PRIM_REDSHIFT=""
METAL_REDSHIFT=""
PRIM_TK0=""
PRIM_YE0=""
PRIM_YH2=""
PRIM_YHD=""
PRIM_ABUNDANCE_SET=""
METAL_TK0=""
METAL_YE0=""
METAL_YH2=""
METAL_YHD=""
METAL_ABUNDANCE_SET=""
PRIM_XNH0=""
PRIM_OUTPUT_STRIDE=""
PRIM_MAX_ITER=""
PRIM_DT_FACTOR=""
PRIM_DT_FACTOR_INIT=""
PRIM_N_INIT_STEPS=""
PRIM_XNH_STOP=""
PRIM_CR_ATTEN_COL_DENS=""
PRIM_DATA_DIR=""
METAL_XNH0=""
METAL_OUTPUT_STRIDE=""
METAL_MAX_ITER=""
METAL_DT_FACTOR=""
METAL_DT_FACTOR_INIT=""
METAL_N_INIT_STEPS=""
METAL_XNH_STOP=""
METAL_CR_ATTEN_COL_DENS=""
METAL_CR_ATTEN_SECOND_FRAC=""
METAL_CR_METAL_BKGND=""
METAL_SRA_RATE=""
METAL_LRA_RATE=""
METAL_T_CR_DES=""
METAL_C_GAS_FRAC=""
METAL_O_GAS_FRAC=""
METAL_MG_GAS_FRAC=""
METAL_DATA_DIR=""
BUILD_DIR=""
OUT_DIR=""
SAVE_DIR=""
CONFIG_PATH=""
USE_CONFIG=1
CONFIG_SOURCE=""
DEFAULT_CONFIG_PATH="${SCRIPT_DIR}/params/default.conf"
DO_BUILD=1
DO_PRIM=1
DO_METAL=1
DO_PLOT=1
DO_FIG_COMBO=0
DO_RESAMPLE=1

# ── Argument parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --common-zeta0)  COMMON_ZETA0="$2";  shift 2 ;;
        --common-ff-ret) COMMON_FF_RET="$2"; shift 2 ;;
        --common-fret-table) COMMON_FRET_TABLE="$2"; shift 2 ;;
        --common-jlw21)  COMMON_JLW21="$2";  shift 2 ;;
        --common-redshift) COMMON_REDSHIFT="$2"; shift 2 ;;
        --common-tk0)    COMMON_TK0="$2";    shift 2 ;;
        --common-ye0)    COMMON_YE0="$2";    shift 2 ;;
        --common-yh2)    COMMON_YH2="$2";    shift 2 ;;
        --common-yhd)    COMMON_YHD="$2";    shift 2 ;;
        --common-abundance-set) COMMON_ABUNDANCE_SET="$2"; shift 2 ;;
        --common-xnh0)   COMMON_XNH0="$2";   shift 2 ;;
        --common-output-stride) COMMON_OUTPUT_STRIDE="$2"; shift 2 ;;
        --common-max-iter) COMMON_MAX_ITER="$2"; shift 2 ;;
        --common-dt-factor) COMMON_DT_FACTOR="$2"; shift 2 ;;
        --common-dt-factor-init) COMMON_DT_FACTOR_INIT="$2"; shift 2 ;;
        --common-n-init-steps) COMMON_N_INIT_STEPS="$2"; shift 2 ;;
        --common-xnh-stop) COMMON_XNH_STOP="$2"; shift 2 ;;
        --common-cr-col-dens) COMMON_CR_ATTEN_COL_DENS="$2"; shift 2 ;;
        --common-data-dir) COMMON_DATA_DIR="$2"; shift 2 ;;
        --prim-zeta0)    PRIM_ZETA0="$2";    shift 2 ;;
        --metal-zeta0)   METAL_ZETA0="$2";   shift 2 ;;
        --metal-z-metal) METAL_Z_METAL="$2"; shift 2 ;;
        --prim-ff-ret)      PRIM_FF_RET="$2";      shift 2 ;;
        --prim-fret-table)  PRIM_FRET_TABLE="$2";  shift 2 ;;
        --metal-ff-ret)     METAL_FF_RET="$2";     shift 2 ;;
        --metal-fret-table) METAL_FRET_TABLE="$2"; shift 2 ;;
        --prim-jlw21)       PRIM_JLW21="$2";       shift 2 ;;
        --metal-jlw21)      METAL_JLW21="$2";      shift 2 ;;
        --prim-redshift)    PRIM_REDSHIFT="$2";    shift 2 ;;
        --metal-redshift)   METAL_REDSHIFT="$2";   shift 2 ;;
        --prim-tk0)         PRIM_TK0="$2";         shift 2 ;;
        --prim-ye0)         PRIM_YE0="$2";         shift 2 ;;
        --prim-yh2)         PRIM_YH2="$2";         shift 2 ;;
        --prim-yhd)         PRIM_YHD="$2";         shift 2 ;;
        --prim-abundance-set) PRIM_ABUNDANCE_SET="$2"; shift 2 ;;
        --prim-xnh0)        PRIM_XNH0="$2";        shift 2 ;;
        --prim-output-stride) PRIM_OUTPUT_STRIDE="$2"; shift 2 ;;
        --prim-max-iter)    PRIM_MAX_ITER="$2";    shift 2 ;;
        --prim-dt-factor)   PRIM_DT_FACTOR="$2"; shift 2 ;;
        --prim-dt-factor-init) PRIM_DT_FACTOR_INIT="$2"; shift 2 ;;
        --prim-n-init-steps) PRIM_N_INIT_STEPS="$2"; shift 2 ;;
        --prim-xnh-stop)    PRIM_XNH_STOP="$2"; shift 2 ;;
        --prim-cr-col-dens) PRIM_CR_ATTEN_COL_DENS="$2"; shift 2 ;;
        --prim-data-dir)    PRIM_DATA_DIR="$2";    shift 2 ;;
        --metal-tk0)        METAL_TK0="$2";        shift 2 ;;
        --metal-ye0)        METAL_YE0="$2";        shift 2 ;;
        --metal-yh2)        METAL_YH2="$2";        shift 2 ;;
        --metal-yhd)        METAL_YHD="$2";        shift 2 ;;
        --metal-abundance-set) METAL_ABUNDANCE_SET="$2"; shift 2 ;;
        --metal-xnh0)       METAL_XNH0="$2";       shift 2 ;;
        --metal-output-stride) METAL_OUTPUT_STRIDE="$2"; shift 2 ;;
        --metal-max-iter)   METAL_MAX_ITER="$2";   shift 2 ;;
        --metal-dt-factor)  METAL_DT_FACTOR="$2"; shift 2 ;;
        --metal-dt-factor-init) METAL_DT_FACTOR_INIT="$2"; shift 2 ;;
        --metal-n-init-steps) METAL_N_INIT_STEPS="$2"; shift 2 ;;
        --metal-xnh-stop)   METAL_XNH_STOP="$2"; shift 2 ;;
        --metal-cr-col-dens) METAL_CR_ATTEN_COL_DENS="$2"; shift 2 ;;
        --metal-cr-second-frac) METAL_CR_ATTEN_SECOND_FRAC="$2"; shift 2 ;;
        --metal-cr-metal-bkgnd) METAL_CR_METAL_BKGND="$2"; shift 2 ;;
        --metal-sra-rate) METAL_SRA_RATE="$2"; shift 2 ;;
        --metal-lra-rate) METAL_LRA_RATE="$2"; shift 2 ;;
        --metal-t-cr-des) METAL_T_CR_DES="$2"; shift 2 ;;
        --metal-c-gas-frac) METAL_C_GAS_FRAC="$2"; shift 2 ;;
        --metal-o-gas-frac) METAL_O_GAS_FRAC="$2"; shift 2 ;;
        --metal-mg-gas-frac) METAL_MG_GAS_FRAC="$2"; shift 2 ;;
        --metal-data-dir)   METAL_DATA_DIR="$2";   shift 2 ;;
        --build-dir)  BUILD_DIR="$2"; shift 2 ;;
        --out-dir)    OUT_DIR="$2";   shift 2 ;;
        --save-dir)   SAVE_DIR="$2";  shift 2 ;;
        --config)     CONFIG_PATH="$2"; USE_CONFIG=1; shift 2 ;;
        --no-config)  USE_CONFIG=0; shift ;;
        --no-build)    DO_BUILD=0;     shift ;;
        --no-prim)     DO_PRIM=0;      shift ;;
        --no-metal)    DO_METAL=0;     shift ;;
        --no-plot)     DO_PLOT=0;      shift ;;
        --fig-combo)   DO_FIG_COMBO=1; shift ;;
        --resample)    DO_RESAMPLE=1;  shift ;;
        --no-resample) DO_RESAMPLE=0;  shift ;;
        -h|--help)
            sed -n '2,/^set -/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Expand shared CLI parameters into per-network values unless explicitly overridden.
apply_common_cli() {
    local common_val="$1" prim_val="$2" metal_val="$3"
    if [[ -z "$prim_val" && -n "$common_val" ]]; then
        prim_val="$common_val"
    fi
    if [[ -z "$metal_val" && -n "$common_val" ]]; then
        metal_val="$common_val"
    fi
    printf '%s\n%s\n' "$prim_val" "$metal_val"
}

mapfile -t _pair < <(apply_common_cli "$COMMON_ZETA0" "$PRIM_ZETA0" "$METAL_ZETA0")
PRIM_ZETA0="${_pair[0]}"; METAL_ZETA0="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_FF_RET" "$PRIM_FF_RET" "$METAL_FF_RET")
PRIM_FF_RET="${_pair[0]}"; METAL_FF_RET="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_FRET_TABLE" "$PRIM_FRET_TABLE" "$METAL_FRET_TABLE")
PRIM_FRET_TABLE="${_pair[0]}"; METAL_FRET_TABLE="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_JLW21" "$PRIM_JLW21" "$METAL_JLW21")
PRIM_JLW21="${_pair[0]}"; METAL_JLW21="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_REDSHIFT" "$PRIM_REDSHIFT" "$METAL_REDSHIFT")
PRIM_REDSHIFT="${_pair[0]}"; METAL_REDSHIFT="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_TK0" "$PRIM_TK0" "$METAL_TK0")
PRIM_TK0="${_pair[0]}"; METAL_TK0="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_YE0" "$PRIM_YE0" "$METAL_YE0")
PRIM_YE0="${_pair[0]}"; METAL_YE0="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_YH2" "$PRIM_YH2" "$METAL_YH2")
PRIM_YH2="${_pair[0]}"; METAL_YH2="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_YHD" "$PRIM_YHD" "$METAL_YHD")
PRIM_YHD="${_pair[0]}"; METAL_YHD="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_ABUNDANCE_SET" "$PRIM_ABUNDANCE_SET" "$METAL_ABUNDANCE_SET")
PRIM_ABUNDANCE_SET="${_pair[0]}"; METAL_ABUNDANCE_SET="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_XNH0" "$PRIM_XNH0" "$METAL_XNH0")
PRIM_XNH0="${_pair[0]}"; METAL_XNH0="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_OUTPUT_STRIDE" "$PRIM_OUTPUT_STRIDE" "$METAL_OUTPUT_STRIDE")
PRIM_OUTPUT_STRIDE="${_pair[0]}"; METAL_OUTPUT_STRIDE="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_MAX_ITER" "$PRIM_MAX_ITER" "$METAL_MAX_ITER")
PRIM_MAX_ITER="${_pair[0]}"; METAL_MAX_ITER="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_DT_FACTOR" "$PRIM_DT_FACTOR" "$METAL_DT_FACTOR")
PRIM_DT_FACTOR="${_pair[0]}"; METAL_DT_FACTOR="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_DT_FACTOR_INIT" "$PRIM_DT_FACTOR_INIT" "$METAL_DT_FACTOR_INIT")
PRIM_DT_FACTOR_INIT="${_pair[0]}"; METAL_DT_FACTOR_INIT="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_N_INIT_STEPS" "$PRIM_N_INIT_STEPS" "$METAL_N_INIT_STEPS")
PRIM_N_INIT_STEPS="${_pair[0]}"; METAL_N_INIT_STEPS="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_XNH_STOP" "$PRIM_XNH_STOP" "$METAL_XNH_STOP")
PRIM_XNH_STOP="${_pair[0]}"; METAL_XNH_STOP="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_CR_ATTEN_COL_DENS" "$PRIM_CR_ATTEN_COL_DENS" "$METAL_CR_ATTEN_COL_DENS")
PRIM_CR_ATTEN_COL_DENS="${_pair[0]}"; METAL_CR_ATTEN_COL_DENS="${_pair[1]}"
mapfile -t _pair < <(apply_common_cli "$COMMON_DATA_DIR" "$PRIM_DATA_DIR" "$METAL_DATA_DIR")
PRIM_DATA_DIR="${_pair[0]}"; METAL_DATA_DIR="${_pair[1]}"

# ── Config loading and precedence merge ──────────────────────────────────────
declare -A CFG

is_allowed_config_key() {
    case "$1" in
        PRIM_ZETA0|METAL_ZETA0|METAL_Z_METAL|PRIM_FF_RET|PRIM_FRET_TABLE|METAL_FF_RET|METAL_FRET_TABLE|PRIM_JLW21|METAL_JLW21|PRIM_REDSHIFT|METAL_REDSHIFT|PRIM_TK0|PRIM_YE0|PRIM_YH2|PRIM_YHD|PRIM_ABUNDANCE_SET|METAL_TK0|METAL_YE0|METAL_YH2|METAL_YHD|METAL_ABUNDANCE_SET|PRIM_XNH0|PRIM_OUTPUT_STRIDE|PRIM_MAX_ITER|PRIM_DT_FACTOR|PRIM_DT_FACTOR_INIT|PRIM_N_INIT_STEPS|PRIM_XNH_STOP|PRIM_CR_ATTEN_COL_DENS|PRIM_DATA_DIR|METAL_XNH0|METAL_OUTPUT_STRIDE|METAL_MAX_ITER|METAL_DT_FACTOR|METAL_DT_FACTOR_INIT|METAL_N_INIT_STEPS|METAL_XNH_STOP|METAL_CR_ATTEN_COL_DENS|METAL_CR_ATTEN_SECOND_FRAC|METAL_CR_METAL_BKGND|METAL_SRA_RATE|METAL_LRA_RATE|METAL_T_CR_DES|METAL_C_GAS_FRAC|METAL_O_GAS_FRAC|METAL_MG_GAS_FRAC|METAL_DATA_DIR|BUILD_DIR|OUT_DIR|SAVE_DIR|COMMON_ZETA0|COMMON_FF_RET|COMMON_FRET_TABLE|COMMON_JLW21|COMMON_REDSHIFT|COMMON_TK0|COMMON_YE0|COMMON_YH2|COMMON_YHD|COMMON_ABUNDANCE_SET|COMMON_XNH0|COMMON_OUTPUT_STRIDE|COMMON_MAX_ITER|COMMON_DT_FACTOR|COMMON_DT_FACTOR_INIT|COMMON_N_INIT_STEPS|COMMON_XNH_STOP|COMMON_CR_ATTEN_COL_DENS|COMMON_DATA_DIR)
            return 0 ;;
        *)
            return 1 ;;
    esac
}

trim_ws() {
    local s="$1"
    s="${s#"${s%%[![:space:]]*}"}"
    s="${s%"${s##*[![:space:]]}"}"
    printf '%s' "$s"
}

load_config_file() {
    local path="$1"
    [[ -f "$path" ]] || { echo "ERROR: config file not found: $path" >&2; exit 1; }

    local line_no=0 raw line key val
    while IFS= read -r raw || [[ -n "$raw" ]]; do
        line_no=$((line_no + 1))
        line="$(trim_ws "$raw")"
        [[ -z "$line" || "${line:0:1}" == "#" ]] && continue

        if [[ "$line" != *=* ]]; then
            echo "ERROR: ${path}:${line_no}: expected KEY=VALUE format" >&2
            exit 1
        fi

        key="$(trim_ws "${line%%=*}")"
        val="$(trim_ws "${line#*=}")"

        if ! [[ "$key" =~ ^[A-Z0-9_]+$ ]]; then
            echo "ERROR: ${path}:${line_no}: invalid key '$key'" >&2
            exit 1
        fi
        if ! is_allowed_config_key "$key"; then
            echo "ERROR: ${path}:${line_no}: unknown key '$key'" >&2
            exit 1
        fi

        if [[ "$val" == \"*\" && "$val" == *\" && ${#val} -ge 2 ]]; then
            val="${val:1:${#val}-2}"
        elif [[ "$val" == \'*\' && "$val" == *\' && ${#val} -ge 2 ]]; then
            val="${val:1:${#val}-2}"
        fi
        CFG["$key"]="$val"
    done < "$path"
}

pick_value() {
    local cli="$1" env="$2" cfg="$3"
    if [[ -n "$cli" ]]; then
        printf '%s' "$cli"
    elif [[ -n "$env" ]]; then
        printf '%s' "$env"
    else
        printf '%s' "$cfg"
    fi
}

cfg_or_common() {
    local specific="$1" common="$2"
    if [[ -n "${CFG[$specific]-}" ]]; then
        printf '%s' "${CFG[$specific]}"
    else
        printf '%s' "${CFG[$common]-}"
    fi
}

if [[ $USE_CONFIG -eq 1 ]]; then
    if [[ -n "$CONFIG_PATH" ]]; then
        CONFIG_SOURCE="$CONFIG_PATH"
        load_config_file "$CONFIG_SOURCE"
    else
        CONFIG_SOURCE="$DEFAULT_CONFIG_PATH"
        load_config_file "$CONFIG_SOURCE"
    fi
fi

PRIM_ZETA0="$(pick_value "$PRIM_ZETA0" "$ENV_PRIM_ZETA0" "$(cfg_or_common PRIM_ZETA0 COMMON_ZETA0)")"
METAL_ZETA0="$(pick_value "$METAL_ZETA0" "$ENV_METAL_ZETA0" "$(cfg_or_common METAL_ZETA0 COMMON_ZETA0)")"
METAL_Z_METAL="$(pick_value "$METAL_Z_METAL" "$ENV_METAL_Z_METAL" "${CFG[METAL_Z_METAL]-}")"
PRIM_FF_RET="$(pick_value "$PRIM_FF_RET" "$ENV_PRIM_FF_RET" "$(cfg_or_common PRIM_FF_RET COMMON_FF_RET)")"
PRIM_FRET_TABLE="$(pick_value "$PRIM_FRET_TABLE" "$ENV_PRIM_FRET_TABLE" "$(cfg_or_common PRIM_FRET_TABLE COMMON_FRET_TABLE)")"
METAL_FF_RET="$(pick_value "$METAL_FF_RET" "$ENV_METAL_FF_RET" "$(cfg_or_common METAL_FF_RET COMMON_FF_RET)")"
METAL_FRET_TABLE="$(pick_value "$METAL_FRET_TABLE" "$ENV_METAL_FRET_TABLE" "$(cfg_or_common METAL_FRET_TABLE COMMON_FRET_TABLE)")"
PRIM_JLW21="$(pick_value "$PRIM_JLW21" "$ENV_PRIM_JLW21" "$(cfg_or_common PRIM_JLW21 COMMON_JLW21)")"
METAL_JLW21="$(pick_value "$METAL_JLW21" "$ENV_METAL_JLW21" "$(cfg_or_common METAL_JLW21 COMMON_JLW21)")"
PRIM_REDSHIFT="$(pick_value "$PRIM_REDSHIFT" "$ENV_PRIM_REDSHIFT" "$(cfg_or_common PRIM_REDSHIFT COMMON_REDSHIFT)")"
METAL_REDSHIFT="$(pick_value "$METAL_REDSHIFT" "$ENV_METAL_REDSHIFT" "$(cfg_or_common METAL_REDSHIFT COMMON_REDSHIFT)")"
PRIM_TK0="$(pick_value "$PRIM_TK0" "$ENV_PRIM_TK0" "$(cfg_or_common PRIM_TK0 COMMON_TK0)")"
PRIM_YE0="$(pick_value "$PRIM_YE0" "$ENV_PRIM_YE0" "$(cfg_or_common PRIM_YE0 COMMON_YE0)")"
PRIM_YH2="$(pick_value "$PRIM_YH2" "$ENV_PRIM_YH2" "$(cfg_or_common PRIM_YH2 COMMON_YH2)")"
PRIM_YHD="$(pick_value "$PRIM_YHD" "$ENV_PRIM_YHD" "$(cfg_or_common PRIM_YHD COMMON_YHD)")"
PRIM_ABUNDANCE_SET="$(pick_value "$PRIM_ABUNDANCE_SET" "$ENV_PRIM_ABUNDANCE_SET" "$(cfg_or_common PRIM_ABUNDANCE_SET COMMON_ABUNDANCE_SET)")"
METAL_TK0="$(pick_value "$METAL_TK0" "$ENV_METAL_TK0" "$(cfg_or_common METAL_TK0 COMMON_TK0)")"
METAL_YE0="$(pick_value "$METAL_YE0" "$ENV_METAL_YE0" "$(cfg_or_common METAL_YE0 COMMON_YE0)")"
METAL_YH2="$(pick_value "$METAL_YH2" "$ENV_METAL_YH2" "$(cfg_or_common METAL_YH2 COMMON_YH2)")"
METAL_YHD="$(pick_value "$METAL_YHD" "$ENV_METAL_YHD" "$(cfg_or_common METAL_YHD COMMON_YHD)")"
METAL_ABUNDANCE_SET="$(pick_value "$METAL_ABUNDANCE_SET" "$ENV_METAL_ABUNDANCE_SET" "$(cfg_or_common METAL_ABUNDANCE_SET COMMON_ABUNDANCE_SET)")"
PRIM_XNH0="$(pick_value "$PRIM_XNH0" "$ENV_PRIM_XNH0" "$(cfg_or_common PRIM_XNH0 COMMON_XNH0)")"
PRIM_OUTPUT_STRIDE="$(pick_value "$PRIM_OUTPUT_STRIDE" "$ENV_PRIM_OUTPUT_STRIDE" "$(cfg_or_common PRIM_OUTPUT_STRIDE COMMON_OUTPUT_STRIDE)")"
PRIM_MAX_ITER="$(pick_value "$PRIM_MAX_ITER" "$ENV_PRIM_MAX_ITER" "$(cfg_or_common PRIM_MAX_ITER COMMON_MAX_ITER)")"
PRIM_DT_FACTOR="$(pick_value "$PRIM_DT_FACTOR" "$ENV_PRIM_DT_FACTOR" "$(cfg_or_common PRIM_DT_FACTOR COMMON_DT_FACTOR)")"
PRIM_DT_FACTOR_INIT="$(pick_value "$PRIM_DT_FACTOR_INIT" "$ENV_PRIM_DT_FACTOR_INIT" "$(cfg_or_common PRIM_DT_FACTOR_INIT COMMON_DT_FACTOR_INIT)")"
PRIM_N_INIT_STEPS="$(pick_value "$PRIM_N_INIT_STEPS" "$ENV_PRIM_N_INIT_STEPS" "$(cfg_or_common PRIM_N_INIT_STEPS COMMON_N_INIT_STEPS)")"
PRIM_XNH_STOP="$(pick_value "$PRIM_XNH_STOP" "$ENV_PRIM_XNH_STOP" "$(cfg_or_common PRIM_XNH_STOP COMMON_XNH_STOP)")"
PRIM_CR_ATTEN_COL_DENS="$(pick_value "$PRIM_CR_ATTEN_COL_DENS" "$ENV_PRIM_CR_ATTEN_COL_DENS" "$(cfg_or_common PRIM_CR_ATTEN_COL_DENS COMMON_CR_ATTEN_COL_DENS)")"
PRIM_DATA_DIR="$(pick_value "$PRIM_DATA_DIR" "$ENV_PRIM_DATA_DIR" "$(cfg_or_common PRIM_DATA_DIR COMMON_DATA_DIR)")"
METAL_XNH0="$(pick_value "$METAL_XNH0" "$ENV_METAL_XNH0" "$(cfg_or_common METAL_XNH0 COMMON_XNH0)")"
METAL_OUTPUT_STRIDE="$(pick_value "$METAL_OUTPUT_STRIDE" "$ENV_METAL_OUTPUT_STRIDE" "$(cfg_or_common METAL_OUTPUT_STRIDE COMMON_OUTPUT_STRIDE)")"
METAL_MAX_ITER="$(pick_value "$METAL_MAX_ITER" "$ENV_METAL_MAX_ITER" "$(cfg_or_common METAL_MAX_ITER COMMON_MAX_ITER)")"
METAL_DT_FACTOR="$(pick_value "$METAL_DT_FACTOR" "$ENV_METAL_DT_FACTOR" "$(cfg_or_common METAL_DT_FACTOR COMMON_DT_FACTOR)")"
METAL_DT_FACTOR_INIT="$(pick_value "$METAL_DT_FACTOR_INIT" "$ENV_METAL_DT_FACTOR_INIT" "$(cfg_or_common METAL_DT_FACTOR_INIT COMMON_DT_FACTOR_INIT)")"
METAL_N_INIT_STEPS="$(pick_value "$METAL_N_INIT_STEPS" "$ENV_METAL_N_INIT_STEPS" "$(cfg_or_common METAL_N_INIT_STEPS COMMON_N_INIT_STEPS)")"
METAL_XNH_STOP="$(pick_value "$METAL_XNH_STOP" "$ENV_METAL_XNH_STOP" "$(cfg_or_common METAL_XNH_STOP COMMON_XNH_STOP)")"
METAL_CR_ATTEN_COL_DENS="$(pick_value "$METAL_CR_ATTEN_COL_DENS" "$ENV_METAL_CR_ATTEN_COL_DENS" "$(cfg_or_common METAL_CR_ATTEN_COL_DENS COMMON_CR_ATTEN_COL_DENS)")"
METAL_CR_ATTEN_SECOND_FRAC="$(pick_value "$METAL_CR_ATTEN_SECOND_FRAC" "$ENV_METAL_CR_ATTEN_SECOND_FRAC" "${CFG[METAL_CR_ATTEN_SECOND_FRAC]-}")"
METAL_CR_METAL_BKGND="$(pick_value "$METAL_CR_METAL_BKGND" "$ENV_METAL_CR_METAL_BKGND" "${CFG[METAL_CR_METAL_BKGND]-}")"
METAL_SRA_RATE="$(pick_value "$METAL_SRA_RATE" "$ENV_METAL_SRA_RATE" "${CFG[METAL_SRA_RATE]-}")"
METAL_LRA_RATE="$(pick_value "$METAL_LRA_RATE" "$ENV_METAL_LRA_RATE" "${CFG[METAL_LRA_RATE]-}")"
METAL_T_CR_DES="$(pick_value "$METAL_T_CR_DES" "$ENV_METAL_T_CR_DES" "${CFG[METAL_T_CR_DES]-}")"
METAL_C_GAS_FRAC="$(pick_value "$METAL_C_GAS_FRAC" "$ENV_METAL_C_GAS_FRAC" "${CFG[METAL_C_GAS_FRAC]-}")"
METAL_O_GAS_FRAC="$(pick_value "$METAL_O_GAS_FRAC" "$ENV_METAL_O_GAS_FRAC" "${CFG[METAL_O_GAS_FRAC]-}")"
METAL_MG_GAS_FRAC="$(pick_value "$METAL_MG_GAS_FRAC" "$ENV_METAL_MG_GAS_FRAC" "${CFG[METAL_MG_GAS_FRAC]-}")"
METAL_DATA_DIR="$(pick_value "$METAL_DATA_DIR" "$ENV_METAL_DATA_DIR" "$(cfg_or_common METAL_DATA_DIR COMMON_DATA_DIR)")"
BUILD_DIR="$(pick_value "$BUILD_DIR" "$ENV_BUILD_DIR" "${CFG[BUILD_DIR]-}")"
OUT_DIR="$(pick_value "$OUT_DIR" "$ENV_OUT_DIR" "${CFG[OUT_DIR]-}")"
SAVE_DIR="$(pick_value "$SAVE_DIR" "$ENV_SAVE_DIR" "${CFG[SAVE_DIR]-}")"

[[ -z "$BUILD_DIR" ]] && BUILD_DIR="build"
[[ -z "$OUT_DIR"   ]] && OUT_DIR="results/"
[[ -z "$SAVE_DIR"  ]] && SAVE_DIR="results/"

cd "$SCRIPT_DIR"

# ── Validation ───────────────────────────────────────────────────────────────
# Check that value is a valid number and >= 0 (awk only, no python3 required)
validate_nonneg() {
    local name="$1" val="$2"
    awk -v n="$name" -v v="$val" 'BEGIN {
        if (v !~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/) {
            print "ERROR: " n "=" v " is not a valid number" > "/dev/stderr"; exit 1
        }
        if (v+0 < 0) {
            print "ERROR: " n " must be >= 0, got " v > "/dev/stderr"; exit 1
        }
    }' || exit 1
}

# Check that value is a valid number and > 0 (used for free-fall retardation factor)
validate_pos() {
    local name="$1" val="$2"
    awk -v n="$name" -v v="$val" 'BEGIN {
        if (v !~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/) {
            print "ERROR: " n "=" v " is not a valid number" > "/dev/stderr"; exit 1
        }
        if (v+0 <= 0) {
            print "ERROR: " n " must be > 0, got " v > "/dev/stderr"; exit 1
        }
    }' || exit 1
}

validate_abundance_set() {
    local name="$1" val="$2"
    local lc="${val,,}"
    case "$lc" in
        solar|default|alpha-enhanced) ;;
        *)
            echo "ERROR: $name must be one of: solar, default, alpha-enhanced (got '$val')" >&2
            exit 1
            ;;
    esac
}

if [[ $DO_PRIM -eq 1 ]]; then
    if [[ -z "$PRIM_ZETA0" ]]; then
        echo "ERROR: --prim-zeta0 is required (use --no-prim to skip)" >&2; exit 1
    fi
    validate_nonneg "--prim-zeta0" "$PRIM_ZETA0"
    if [[ -n "$PRIM_FRET_TABLE" ]]; then
        [[ -f "$PRIM_FRET_TABLE" ]] || {
            echo "ERROR: --prim-fret-table file not found: $PRIM_FRET_TABLE" >&2; exit 1
        }
    else
        [[ -n "$PRIM_FF_RET" ]] && validate_pos "--prim-ff-ret" "$PRIM_FF_RET"
    fi
    [[ -n "$PRIM_JLW21"    ]] && validate_nonneg "--prim-jlw21"    "$PRIM_JLW21"
    [[ -n "$PRIM_REDSHIFT" ]] && validate_nonneg "--prim-redshift" "$PRIM_REDSHIFT"
    [[ -n "$PRIM_TK0"      ]] && validate_pos    "--prim-tk0"      "$PRIM_TK0"
    [[ -n "$PRIM_YE0"      ]] && validate_nonneg "--prim-ye0"      "$PRIM_YE0"
    [[ -n "$PRIM_YH2"      ]] && validate_nonneg "--prim-yh2"      "$PRIM_YH2"
    [[ -n "$PRIM_YHD"      ]] && validate_nonneg "--prim-yhd"      "$PRIM_YHD"
    [[ -n "$PRIM_ABUNDANCE_SET" ]] && validate_abundance_set "--prim-abundance-set" "$PRIM_ABUNDANCE_SET"
    [[ -n "$PRIM_DT_FACTOR" ]] && validate_pos "--prim-dt-factor" "$PRIM_DT_FACTOR"
    [[ -n "$PRIM_DT_FACTOR_INIT" ]] && validate_pos "--prim-dt-factor-init" "$PRIM_DT_FACTOR_INIT"
    [[ -n "$PRIM_XNH_STOP" ]] && validate_pos "--prim-xnh-stop" "$PRIM_XNH_STOP"
    [[ -n "$PRIM_CR_ATTEN_COL_DENS" ]] && validate_pos "--prim-cr-col-dens" "$PRIM_CR_ATTEN_COL_DENS"
    if [[ -n "$PRIM_N_INIT_STEPS" && ! "$PRIM_N_INIT_STEPS" =~ ^[0-9]+$ ]]; then
        echo "ERROR: --prim-n-init-steps must be an integer >= 0, got $PRIM_N_INIT_STEPS" >&2; exit 1
    fi
fi

if [[ $DO_METAL -eq 1 ]]; then
    if [[ -z "$METAL_ZETA0" ]]; then
        echo "ERROR: --metal-zeta0 is required (use --no-metal to skip)" >&2; exit 1
    fi
    if [[ -z "$METAL_Z_METAL" ]]; then
        echo "ERROR: --metal-z-metal is required (use --no-metal to skip)" >&2; exit 1
    fi
    validate_nonneg "--metal-zeta0"   "$METAL_ZETA0"
    validate_nonneg "--metal-z-metal" "$METAL_Z_METAL"
    if [[ -n "$METAL_FRET_TABLE" ]]; then
        [[ -f "$METAL_FRET_TABLE" ]] || {
            echo "ERROR: --metal-fret-table file not found: $METAL_FRET_TABLE" >&2; exit 1
        }
    else
        [[ -n "$METAL_FF_RET" ]] && validate_pos "--metal-ff-ret" "$METAL_FF_RET"
    fi
    [[ -n "$METAL_JLW21"    ]] && validate_nonneg "--metal-jlw21"    "$METAL_JLW21"
    [[ -n "$METAL_REDSHIFT" ]] && validate_nonneg "--metal-redshift" "$METAL_REDSHIFT"
    [[ -n "$METAL_TK0"      ]] && validate_pos    "--metal-tk0"      "$METAL_TK0"
    [[ -n "$METAL_YE0"      ]] && validate_nonneg "--metal-ye0"      "$METAL_YE0"
    [[ -n "$METAL_YH2"      ]] && validate_nonneg "--metal-yh2"      "$METAL_YH2"
    [[ -n "$METAL_YHD"      ]] && validate_nonneg "--metal-yhd"      "$METAL_YHD"
    [[ -n "$METAL_ABUNDANCE_SET" ]] && validate_abundance_set "--metal-abundance-set" "$METAL_ABUNDANCE_SET"
    [[ -n "$METAL_DT_FACTOR" ]] && validate_pos "--metal-dt-factor" "$METAL_DT_FACTOR"
    [[ -n "$METAL_DT_FACTOR_INIT" ]] && validate_pos "--metal-dt-factor-init" "$METAL_DT_FACTOR_INIT"
    [[ -n "$METAL_XNH_STOP" ]] && validate_pos "--metal-xnh-stop" "$METAL_XNH_STOP"
    [[ -n "$METAL_CR_ATTEN_COL_DENS" ]] && validate_pos "--metal-cr-col-dens" "$METAL_CR_ATTEN_COL_DENS"
    [[ -n "$METAL_CR_ATTEN_SECOND_FRAC" ]] && validate_nonneg "--metal-cr-second-frac" "$METAL_CR_ATTEN_SECOND_FRAC"
    [[ -n "$METAL_CR_METAL_BKGND" ]] && validate_nonneg "--metal-cr-metal-bkgnd" "$METAL_CR_METAL_BKGND"
    [[ -n "$METAL_SRA_RATE" ]] && validate_nonneg "--metal-sra-rate" "$METAL_SRA_RATE"
    [[ -n "$METAL_LRA_RATE" ]] && validate_nonneg "--metal-lra-rate" "$METAL_LRA_RATE"
    [[ -n "$METAL_T_CR_DES" ]] && validate_pos "--metal-t-cr-des" "$METAL_T_CR_DES"
    [[ -n "$METAL_C_GAS_FRAC" ]] && validate_nonneg "--metal-c-gas-frac" "$METAL_C_GAS_FRAC"
    [[ -n "$METAL_O_GAS_FRAC" ]] && validate_nonneg "--metal-o-gas-frac" "$METAL_O_GAS_FRAC"
    [[ -n "$METAL_MG_GAS_FRAC" ]] && validate_nonneg "--metal-mg-gas-frac" "$METAL_MG_GAS_FRAC"
    if [[ -n "$METAL_N_INIT_STEPS" && ! "$METAL_N_INIT_STEPS" =~ ^[0-9]+$ ]]; then
        echo "ERROR: --metal-n-init-steps must be an integer >= 0, got $METAL_N_INIT_STEPS" >&2; exit 1
    fi
    awk -v c="${METAL_C_GAS_FRAC:-0}" -v o="${METAL_O_GAS_FRAC:-0}" -v m="${METAL_MG_GAS_FRAC:-0}" 'BEGIN{
        if (c>1 || o>1 || m>1) { exit 1 }
    }' || { echo "ERROR: gas fractions must be <= 1" >&2; exit 1; }
fi

# ── Tag generation ───────────────────────────────────────────────────────────
# Rule: replace '.' with 'p' in the input string; zero values (0, 0.0, 0e-0, ...) become "0"
make_tag() {
    awk -v val="$1" 'BEGIN {
        if (val+0 == 0.0) { print "0" }
        else { gsub(/\./, "p", val); print val }
    }'
}

echo "============================================================"
echo " run_collapse.sh"
if [[ $USE_CONFIG -eq 0 ]]; then
    echo "  config     : disabled (--no-config)"
elif [[ -n "$CONFIG_SOURCE" ]]; then
    echo "  config     : ${CONFIG_SOURCE}"
else
    echo "  config     : (none)"
fi
if [[ $DO_PRIM -eq 1 ]]; then
    echo "  primordial : zeta0=${PRIM_ZETA0} s^-1"
    if [[ -n "$PRIM_FRET_TABLE" ]]; then
        echo "               f_ret=step-table (${PRIM_FRET_TABLE})"
    elif [[ -n "$PRIM_FF_RET" ]]; then
        echo "               f_ret=${PRIM_FF_RET}"
    else
        echo "               f_ret=1.0 (free-fall)"
    fi
    [[ -n "$PRIM_JLW21"    ]] && echo "               J_LW21=${PRIM_JLW21} J_21" \
                              || echo "               J_LW21=0.0 (no LW field)"
    [[ -n "$PRIM_REDSHIFT" ]] && echo "               redshift=${PRIM_REDSHIFT}" \
                              || echo "               redshift=0.0 (z=0)"
    [[ -n "$PRIM_TK0"      ]] && echo "               T_K0=${PRIM_TK0} K"    || echo "               T_K0=100.0 K (default)"
    [[ -n "$PRIM_YE0"      ]] && echo "               y_e0=${PRIM_YE0}"      || echo "               y_e0=1e-4 (default)"
    [[ -n "$PRIM_YH2"      ]] && echo "               y_H2=${PRIM_YH2}"      || echo "               y_H2=6e-7 (default)"
    [[ -n "$PRIM_YHD"      ]] && echo "               y_HD=${PRIM_YHD}"      || echo "               y_HD=4e-10 (default)"
    [[ -n "$PRIM_ABUNDANCE_SET" ]] && echo "               abundance_set=${PRIM_ABUNDANCE_SET}" \
                                  || echo "               abundance_set=solar (default)"
fi
if [[ $DO_METAL -eq 1 ]]; then
    echo "  metal_grain: zeta0=${METAL_ZETA0} s^-1"
    echo "               Z_metal=${METAL_Z_METAL} Z☉"
    if [[ -n "$METAL_FRET_TABLE" ]]; then
        echo "               f_ret=step-table (${METAL_FRET_TABLE})"
    elif [[ -n "$METAL_FF_RET" ]]; then
        echo "               f_ret=${METAL_FF_RET}"
    else
        echo "               f_ret=1.0 (free-fall)"
    fi
    [[ -n "$METAL_JLW21"    ]] && echo "               J_LW21=${METAL_JLW21} J_21" \
                               || echo "               J_LW21=0.0 (no LW field)"
    [[ -n "$METAL_REDSHIFT" ]] && echo "               redshift=${METAL_REDSHIFT}" \
                               || echo "               redshift=0.0 (z=0)"
    [[ -n "$METAL_TK0"      ]] && echo "               T_K0=${METAL_TK0} K"   || echo "               T_K0=100.0 K (default)"
    [[ -n "$METAL_YE0"      ]] && echo "               y_e0=${METAL_YE0}"     || echo "               y_e0=1e-4 (default)"
    [[ -n "$METAL_YH2"      ]] && echo "               y_H2=${METAL_YH2}"     || echo "               y_H2=6e-7 (default)"
    [[ -n "$METAL_YHD"      ]] && echo "               y_HD=${METAL_YHD}"     || echo "               y_HD=4e-10 (default)"
    [[ -n "$METAL_ABUNDANCE_SET" ]] && echo "               abundance_set=${METAL_ABUNDANCE_SET}" \
                                   || echo "               abundance_set=solar (default)"
    [[ -n "$METAL_SRA_RATE" ]] && echo "               sra_rate=${METAL_SRA_RATE}" \
                               || echo "               sra_rate=0.0 (default)"
    [[ -n "$METAL_LRA_RATE" ]] && echo "               lra_rate=${METAL_LRA_RATE}" \
                               || echo "               lra_rate=0.0 (default)"
fi
echo "  build_dir  : ${BUILD_DIR}"
echo "  out_dir    : ${OUT_DIR}"
echo "  save_dir   : ${SAVE_DIR}"
echo "============================================================"

# ── [1] Build ────────────────────────────────────────────────────────────────
if [[ $DO_BUILD -eq 1 ]]; then
    echo ""
    echo ">>> [1/5] Build"
    mkdir -p "$BUILD_DIR"
    cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17
    cmake --build "$BUILD_DIR" --target prim_collapse metal_collapse -j"$(nproc)"
else
    echo ">>> [1/5] Build skipped"
fi

# ── [2] primordial_collapse ──────────────────────────────────────────────────
PRIM_OUT="${OUT_DIR}/primordial"
PRIM_TAG=""
PRIM_FRET_TAG=""
PRIM_JLW_TAG=""
PRIM_ZRED_TAG=""
if [[ $DO_PRIM -eq 1 ]]; then
    echo ""
    echo ">>> [2/5] primordial_collapse  zeta0=${PRIM_ZETA0} s^-1"
    mkdir -p "$PRIM_OUT"
    PRIM_OUTDIR="$PRIM_OUT" PRIM_ZETA0="$PRIM_ZETA0" \
        PRIM_FF_RET="$PRIM_FF_RET" PRIM_FRET_TABLE="$PRIM_FRET_TABLE" \
        PRIM_JLW21="$PRIM_JLW21" \
        PRIM_REDSHIFT="$PRIM_REDSHIFT" \
        PRIM_ABUNDANCE_SET="$PRIM_ABUNDANCE_SET" \
        PRIM_TK0="$PRIM_TK0" PRIM_YE0="$PRIM_YE0" \
        PRIM_YH2="$PRIM_YH2" PRIM_YHD="$PRIM_YHD" \
        PRIM_XNH0="$PRIM_XNH0" PRIM_OUTPUT_STRIDE="$PRIM_OUTPUT_STRIDE" \
        PRIM_DT_FACTOR="$PRIM_DT_FACTOR" PRIM_DT_FACTOR_INIT="$PRIM_DT_FACTOR_INIT" \
        PRIM_N_INIT_STEPS="$PRIM_N_INIT_STEPS" PRIM_XNH_STOP="$PRIM_XNH_STOP" \
        PRIM_CR_ATTEN_COL_DENS="$PRIM_CR_ATTEN_COL_DENS" \
        PRIM_MAX_ITER="$PRIM_MAX_ITER" PRIM_DATA_DIR="$PRIM_DATA_DIR" \
        "${BUILD_DIR}/src/apps/collapse_primordial/prim_collapse"
    PRIM_TAG=$(make_tag "$PRIM_ZETA0")
    # fret tag: table mode → fixed "-step"; scalar → numeric tag (only when != 1.0)
    if [[ -n "$PRIM_FRET_TABLE" ]]; then
        PRIM_FRET_TAG="-step"
    elif [[ -n "$PRIM_FF_RET" ]]; then
        # append fret tag only when f_ret != 1.0
        if awk -v v="$PRIM_FF_RET" 'BEGIN { exit (v+0 == 1.0) ? 1 : 0 }'; then
            PRIM_FRET_TAG=$(make_tag "$PRIM_FF_RET")
        fi
    fi
    if [[ -n "$PRIM_JLW21" ]]; then
        # append JLW tag only when J_LW21 != 0.0
        if awk -v v="$PRIM_JLW21" 'BEGIN { exit (v+0 == 0.0) ? 1 : 0 }'; then
            PRIM_JLW_TAG=$(make_tag "$PRIM_JLW21")
        fi
    fi
    if [[ -n "$PRIM_REDSHIFT" ]]; then
        # append redshift tag only when z != 0.0
        if awk -v v="$PRIM_REDSHIFT" 'BEGIN { exit (v+0 == 0.0) ? 1 : 0 }'; then
            PRIM_ZRED_TAG=$(make_tag "$PRIM_REDSHIFT")
        fi
    fi
    _prim_h5="${PRIM_OUT}/collapse_CR${PRIM_TAG}${PRIM_FRET_TAG:+_fret${PRIM_FRET_TAG}}${PRIM_JLW_TAG:+_JLW${PRIM_JLW_TAG}}${PRIM_ZRED_TAG:+_z${PRIM_ZRED_TAG}}.h5"
    echo "    HDF5: ${_prim_h5}"
else
    echo ">>> [2/5] primordial_collapse skipped"
fi

# ── [3] metal_grain_collapse ─────────────────────────────────────────────────
METAL_OUT="${OUT_DIR}/metal_grain"
METAL_CR_TAG=""
METAL_Z_TAG=""
METAL_FRET_TAG=""
METAL_JLW_TAG=""
METAL_ZRED_TAG=""
if [[ $DO_METAL -eq 1 ]]; then
    echo ""
    echo ">>> [3/5] metal_grain_collapse  zeta0=${METAL_ZETA0} s^-1  Z_metal=${METAL_Z_METAL} Z☉"
    mkdir -p "$METAL_OUT"
    METAL_OUTDIR="$METAL_OUT" METAL_ZETA0="$METAL_ZETA0" METAL_Z_METAL="$METAL_Z_METAL" \
        METAL_FF_RET="$METAL_FF_RET" METAL_FRET_TABLE="$METAL_FRET_TABLE" \
        METAL_JLW21="$METAL_JLW21" \
        METAL_REDSHIFT="$METAL_REDSHIFT" \
        METAL_ABUNDANCE_SET="$METAL_ABUNDANCE_SET" \
        METAL_TK0="$METAL_TK0" METAL_YE0="$METAL_YE0" \
        METAL_YH2="$METAL_YH2" METAL_YHD="$METAL_YHD" \
        METAL_XNH0="$METAL_XNH0" METAL_OUTPUT_STRIDE="$METAL_OUTPUT_STRIDE" \
        METAL_DT_FACTOR="$METAL_DT_FACTOR" METAL_DT_FACTOR_INIT="$METAL_DT_FACTOR_INIT" \
        METAL_N_INIT_STEPS="$METAL_N_INIT_STEPS" METAL_XNH_STOP="$METAL_XNH_STOP" \
        METAL_CR_ATTEN_COL_DENS="$METAL_CR_ATTEN_COL_DENS" \
        METAL_CR_ATTEN_SECOND_FRAC="$METAL_CR_ATTEN_SECOND_FRAC" \
        METAL_CR_METAL_BKGND="$METAL_CR_METAL_BKGND" \
        METAL_SRA_RATE="$METAL_SRA_RATE" \
        METAL_LRA_RATE="$METAL_LRA_RATE" \
        METAL_T_CR_DES="$METAL_T_CR_DES" \
        METAL_C_GAS_FRAC="$METAL_C_GAS_FRAC" METAL_O_GAS_FRAC="$METAL_O_GAS_FRAC" \
        METAL_MG_GAS_FRAC="$METAL_MG_GAS_FRAC" \
        METAL_MAX_ITER="$METAL_MAX_ITER" METAL_DATA_DIR="$METAL_DATA_DIR" \
        PRIM_DATA_DIR="$PRIM_DATA_DIR" \
        "${BUILD_DIR}/src/apps/collapse_metal_grain/metal_collapse"
    METAL_CR_TAG=$(make_tag "$METAL_ZETA0")
    METAL_Z_TAG=$(make_tag "$METAL_Z_METAL")
    # fret tag: table mode → fixed "-step"; scalar → numeric tag (only when != 1.0)
    if [[ -n "$METAL_FRET_TABLE" ]]; then
        METAL_FRET_TAG="-step"
    elif [[ -n "$METAL_FF_RET" ]]; then
        # append fret tag only when f_ret != 1.0
        if awk -v v="$METAL_FF_RET" 'BEGIN { exit (v+0 == 1.0) ? 1 : 0 }'; then
            METAL_FRET_TAG=$(make_tag "$METAL_FF_RET")
        fi
    fi
    if [[ -n "$METAL_JLW21" ]]; then
        # append JLW tag only when J_LW21 != 0.0
        if awk -v v="$METAL_JLW21" 'BEGIN { exit (v+0 == 0.0) ? 1 : 0 }'; then
            METAL_JLW_TAG=$(make_tag "$METAL_JLW21")
        fi
    fi
    if [[ -n "$METAL_REDSHIFT" ]]; then
        # append redshift tag only when z != 0.0
        if awk -v v="$METAL_REDSHIFT" 'BEGIN { exit (v+0 == 0.0) ? 1 : 0 }'; then
            METAL_ZRED_TAG=$(make_tag "$METAL_REDSHIFT")
        fi
    fi
    _metal_h5="${METAL_OUT}/collapse_CR${METAL_CR_TAG}_Z${METAL_Z_TAG}${METAL_FRET_TAG:+_fret${METAL_FRET_TAG}}${METAL_JLW_TAG:+_JLW${METAL_JLW_TAG}}${METAL_ZRED_TAG:+_z${METAL_ZRED_TAG}}.h5"
    echo "    HDF5: ${_metal_h5}"
else
    echo ">>> [3/5] metal_grain_collapse skipped"
fi

# ── [4] Visualisation ────────────────────────────────────────────────────────
if [[ $DO_PLOT -eq 1 ]]; then
    echo ""
    echo ">>> [4/5] analyze_collapse.py"
    PRIM_SAVE="${SAVE_DIR}/primordial"
    METAL_SAVE="${SAVE_DIR}/metal_grain"
    mkdir -p "$PRIM_SAVE" "$METAL_SAVE"

    PLOT_ARGS=()
    if [[ $DO_PRIM -eq 1 ]]; then
        PLOT_ARGS+=(--h5dir "$PRIM_OUT" --cr-tag "$PRIM_TAG" --save "$PRIM_SAVE")
        [[ -n "$PRIM_FRET_TAG"  ]] && PLOT_ARGS+=(--fret-tag "$PRIM_FRET_TAG")
        [[ -n "$PRIM_JLW_TAG"   ]] && PLOT_ARGS+=(--jlw-tag  "$PRIM_JLW_TAG")
        [[ -n "$PRIM_ZRED_TAG"  ]] && PLOT_ARGS+=(--zred-tag "$PRIM_ZRED_TAG")
    fi
    if [[ $DO_METAL -eq 1 ]]; then
        PLOT_ARGS+=(--metal-h5dir "$METAL_OUT"
                    --metal-cr "$METAL_CR_TAG"
                    --metal-z  "$METAL_Z_TAG"
                    --metal-save "$METAL_SAVE")
        [[ -n "$METAL_FRET_TAG" ]] && PLOT_ARGS+=(--metal-fret "$METAL_FRET_TAG")
        [[ -n "$METAL_JLW_TAG"  ]] && PLOT_ARGS+=(--metal-jlw  "$METAL_JLW_TAG")
        [[ -n "$METAL_ZRED_TAG" ]] && PLOT_ARGS+=(--metal-zred "$METAL_ZRED_TAG")
    fi
    [[ $DO_FIG_COMBO -eq 1 ]] && PLOT_ARGS+=(--fig-combo)

    python3 tools/python/analyze_collapse.py "${PLOT_ARGS[@]}"

    echo ""
    echo "    Figures:"
    [[ $DO_PRIM  -eq 1 ]] && ls "${PRIM_SAVE}"/*.png  2>/dev/null | sed 's/^/      /'
    [[ $DO_METAL -eq 1 ]] && ls "${METAL_SAVE}"/*.png 2>/dev/null | sed 's/^/      /'
else
    echo ">>> [4/5] Plot skipped"
fi

# ── [5] Resample ─────────────────────────────────────────────────────────────
# Convert internal fret tag (-step or numeric) to resample_collapse.py CLI form.
# "-step" → "step"  (avoids argparse leading-dash problem)
# "3p0"   → "3p0"   (unchanged)
_resample_fret_arg() {
    if [[ "$1" == "-step" ]]; then echo "step"; else echo "$1"; fi
}

if [[ $DO_RESAMPLE -eq 1 ]]; then
    if [[ $DO_PRIM -eq 0 && $DO_METAL -eq 0 ]]; then
        echo ">>> [5/5] Resample skipped (both --no-prim and --no-metal)"
    else
        echo ""
        echo ">>> [5/5] resample_collapse.py"

        RESAMPLE_ARGS=()
        if [[ $DO_PRIM -eq 1 && -n "$PRIM_TAG" ]]; then
            RESAMPLE_ARGS+=(--h5dir "$PRIM_OUT" --cr-tag "$PRIM_TAG" --save "$PRIM_OUT")
            if [[ -n "$PRIM_FRET_TAG" ]]; then
                RESAMPLE_ARGS+=(--fret-tag "$(_resample_fret_arg "$PRIM_FRET_TAG")")
            fi
            [[ -n "$PRIM_JLW_TAG"  ]] && RESAMPLE_ARGS+=(--jlw-tag  "$PRIM_JLW_TAG")
            [[ -n "$PRIM_ZRED_TAG" ]] && RESAMPLE_ARGS+=(--zred-tag "$PRIM_ZRED_TAG")
        fi
        if [[ $DO_METAL -eq 1 && -n "$METAL_CR_TAG" ]]; then
            RESAMPLE_ARGS+=(--metal-h5dir "$METAL_OUT"
                            --metal-cr    "$METAL_CR_TAG"
                            --metal-z     "$METAL_Z_TAG"
                            --metal-save  "$METAL_OUT")
            if [[ -n "$METAL_FRET_TAG" ]]; then
                RESAMPLE_ARGS+=(--metal-fret "$(_resample_fret_arg "$METAL_FRET_TAG")")
            fi
            [[ -n "$METAL_JLW_TAG"  ]] && RESAMPLE_ARGS+=(--metal-jlw  "$METAL_JLW_TAG")
            [[ -n "$METAL_ZRED_TAG" ]] && RESAMPLE_ARGS+=(--metal-zred "$METAL_ZRED_TAG")
        fi

        if [[ ${#RESAMPLE_ARGS[@]} -gt 0 ]]; then
            python3 tools/python/resample_collapse.py "${RESAMPLE_ARGS[@]}"
            echo ""
            echo "    CSV files:"
            [[ $DO_PRIM  -eq 1 ]] && ls "${PRIM_OUT}"/resample_*.csv  2>/dev/null | sed 's/^/      /'
            [[ $DO_METAL -eq 1 ]] && ls "${METAL_OUT}"/resample_*.csv 2>/dev/null | sed 's/^/      /'
        fi
    fi
else
    echo ">>> [5/5] Resample skipped"
fi

echo ""
echo "============================================================"
echo " Done."
echo "============================================================"
