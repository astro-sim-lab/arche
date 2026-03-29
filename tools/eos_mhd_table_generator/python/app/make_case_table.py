#!/usr/bin/env python3
"""
make_case_table.py - End-to-end EOS table case runner for tools/eos_mhd_table_generator.

Flow:
1) Run metal_collapse with selected CR/Z/SRA/LRA parameters.
2) Build 10-column EOS table (.d) from produced HDF5.

This helper is intentionally kept under tools/eos_mhd_table_generator (not core pipeline).
"""

from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from lib.build_core import build_table
from lib.config_kv import coalesce_cli_cfg, parse_bool, read_kv_config


def _token_to_float(token: str) -> float:
    t = token.strip().lower().replace("q", "e")
    if t in ("1e-inf", "1e-infinity"):
        return 0.0
    return float(t)


def _norm_label(token: str) -> str:
    return token.strip().lower().replace("q", "e")


def _fmt_label(val: float) -> str:
    if val == 0.0:
        return "0"
    s = f"{val:.6e}".lower()
    m, e = s.split("e")
    m = m.rstrip("0").rstrip(".")
    if m == "":
        m = "0"
    exp = int(e)
    return f"{m}e{exp}"


def _resolve_metal_binary(cli_path: str | None, build_dir: str | None) -> Path:
    candidates = []
    if cli_path:
        candidates.append(Path(cli_path))
    if build_dir:
        candidates.append(Path(build_dir) / "src/apps/collapse_metal_grain/metal_collapse")
    candidates.append(Path("build") / "src/apps/collapse_metal_grain/metal_collapse")
    candidates.append(Path("build-codex") / "src/apps/collapse_metal_grain/metal_collapse")
    for p in candidates:
        if p.is_file():
            return p
    raise FileNotFoundError(
        "metal_collapse binary not found. Build first (cmake --build <build_dir> --target metal_collapse), "
        "or set --metal-binary/--build-dir."
    )


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run one CR/SH/Z case and emit EOS .d table.")
    p.add_argument("--config", help="Optional KEY=VALUE config file. CLI options override config.")
    p.add_argument("--metal-binary", help="Path to metal_collapse binary.")
    p.add_argument("--build-dir", help="Build directory containing metal_collapse.")
    p.add_argument("--output-dir", default="tools/eos_mhd_table_generator/results", help="Output directory for final .d file.")
    p.add_argument("--h5-dir", default="tools/eos_mhd_table_generator/results/h5", help="Directory to store intermediate HDF5.")
    p.add_argument("--input-h5", help="Skip collapse run and use existing HDF5 directly.")

    p.add_argument("--cr-ref", type=float, default=1.0e-17, help="Reference CR rate used for CR scale conversion.")
    p.add_argument("--cr-scale", help="Cosmic-ray scale factor (zeta0 = cr_scale * cr_ref).")
    p.add_argument("--zeta0", help="Absolute CR rate [s^-1]. Overrides --cr-scale.")
    p.add_argument("--z-metal", help="Metallicity [Z_sun]. Supports token '1e-inf' -> 0.")
    p.add_argument("--sh-scale", help="Legacy SH scale. If set, applies to both SRA/LRA unless explicitly overridden.")
    p.add_argument("--sra-rate", help="Short-lived radionuclide rate scale.")
    p.add_argument("--lra-rate", help="Long-lived radionuclide rate scale.")
    p.add_argument("--metal-max-iter", type=int, help="Optional METAL_MAX_ITER override.")
    p.add_argument("--metal-xnh-stop", help="Optional METAL_XNH_STOP override.")
    p.add_argument("--metal-cr-second-frac", help="Optional METAL_CR_ATTEN_SECOND_FRAC override.")
    p.add_argument("--metal-cr-metal-bkgnd", help="Optional METAL_CR_METAL_BKGND override.")
    p.add_argument("--metal-abundance-set", help="Optional METAL_ABUNDANCE_SET override.")
    p.add_argument("--metal-log", help="Optional path to capture metal_collapse stdout/stderr.")

    p.add_argument("--mode", choices=("compat", "normalized"), default="compat")
    p.add_argument("--eta-model", choices=("placeholder", "legacy"), default="legacy")
    p.add_argument("--eta-resample", choices=("state", "native"),
                   help="eta resampling strategy passed to build_eos_table.py")
    p.add_argument("--etaa-compat-mode", choices=("none", "drop_c2_high_n"),
                   help="Optional etaA compatibility post-correction.")
    p.add_argument("--etaa-drop-c2-logn-min", type=float,
                   help="Apply drop_c2_high_n only for log10(nH) >= this threshold.")
    p.add_argument("--charge-table", default="tools/eos_mhd_table_generator/include/legacy_species_charge.dat")
    p.add_argument("--mass-table", default="tools/eos_mhd_table_generator/include/legacy_species_mass_amu.dat")
    p.add_argument("--with-header", action="store_true")
    p.add_argument("--n-log-min", type=float, default=0.1)
    p.add_argument("--n-log-max", type=float, default=22.0)
    p.add_argument("--n-log-step", type=float, default=0.1)
    p.add_argument("--b-log-min", type=float, default=-20.0)
    p.add_argument("--b-log-max", type=float, default=4.0)
    p.add_argument("--b-log-step", type=float, default=0.2)

    p.add_argument("--output-name", help="Optional explicit output filename (e.g. CR1e-2_sh1e-2_Z1e-3.d).")
    return p


def main() -> int:
    args = _build_parser().parse_args()
    cfg: dict[str, str] = {}
    if args.config:
        cfg = read_kv_config(Path(args.config))

    z_metal_raw = coalesce_cli_cfg(args.z_metal, cfg, "EOS_CASE_Z_METAL")
    if z_metal_raw is None:
        raise SystemExit("ERROR: --z-metal or EOS_CASE_Z_METAL is required")
    z_metal = _token_to_float(str(z_metal_raw))
    if z_metal < 0.0:
        raise SystemExit("ERROR: --z-metal must be >= 0")

    cr_ref = float(coalesce_cli_cfg(args.cr_ref, cfg, "EOS_CASE_CR_REF", 1.0e-17))
    zeta0_raw = coalesce_cli_cfg(args.zeta0, cfg, "EOS_CASE_ZETA0")
    cr_scale_raw = coalesce_cli_cfg(args.cr_scale, cfg, "EOS_CASE_CR_SCALE")
    if zeta0_raw is not None:
        zeta0 = _token_to_float(str(zeta0_raw))
        cr_scale = zeta0 / cr_ref if cr_ref != 0.0 else math.nan
    elif cr_scale_raw is not None:
        cr_scale = _token_to_float(str(cr_scale_raw))
        zeta0 = cr_scale * cr_ref
    else:
        raise SystemExit("ERROR: specify either --zeta0/--cr-scale (or EOS_CASE_ZETA0/EOS_CASE_CR_SCALE)")

    if not math.isfinite(zeta0) or zeta0 < 0.0:
        raise SystemExit("ERROR: zeta0 must be finite and >= 0")

    sh_scale_raw = coalesce_cli_cfg(args.sh_scale, cfg, "EOS_CASE_SH_SCALE")
    sra_raw = coalesce_cli_cfg(args.sra_rate, cfg, "EOS_CASE_SRA_RATE")
    lra_raw = coalesce_cli_cfg(args.lra_rate, cfg, "EOS_CASE_LRA_RATE")
    if sra_raw is not None:
        sra_rate = _token_to_float(str(sra_raw))
    elif sh_scale_raw is not None:
        sra_rate = _token_to_float(str(sh_scale_raw))
    else:
        sra_rate = 0.0

    if lra_raw is not None:
        lra_rate = _token_to_float(str(lra_raw))
    elif sh_scale_raw is not None:
        lra_rate = _token_to_float(str(sh_scale_raw))
    else:
        lra_rate = 0.0

    if sra_rate < 0.0 or lra_rate < 0.0:
        raise SystemExit("ERROR: --sra-rate/--lra-rate must be >= 0")

    cr_label = _norm_label(str(cr_scale_raw)) if cr_scale_raw is not None else _fmt_label(cr_scale)
    z_label = _norm_label(str(z_metal_raw)) if z_metal_raw is not None else _fmt_label(z_metal)
    if z_metal == 0.0 and z_label in ("0", "0.0"):
        z_label = "1e-inf"

    if sh_scale_raw is not None:
        sh_label = _norm_label(str(sh_scale_raw))
    elif sra_rate == lra_rate:
        sh_label = _fmt_label(sra_rate)
    else:
        sh_label = ""

    output_name_raw = coalesce_cli_cfg(args.output_name, cfg, "EOS_CASE_OUTPUT_NAME")
    if output_name_raw:
        output_name = str(output_name_raw)
    elif sh_label:
        output_name = f"CR{cr_label}_sh{sh_label}_Z{z_label}.d"
    else:
        output_name = f"CR{cr_label}_sra{_fmt_label(sra_rate)}_lra{_fmt_label(lra_rate)}_Z{z_label}.d"

    output_dir = Path(str(coalesce_cli_cfg(args.output_dir, cfg, "EOS_CASE_OUTPUT_DIR", args.output_dir)))
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / output_name

    input_h5_raw = coalesce_cli_cfg(args.input_h5, cfg, "EOS_CASE_INPUT_H5")
    if input_h5_raw:
        input_h5 = Path(str(input_h5_raw))
        if not input_h5.is_file():
            raise SystemExit(f"ERROR: --input-h5 not found: {input_h5}")
    else:
        metal_binary_raw = coalesce_cli_cfg(args.metal_binary, cfg, "EOS_CASE_METAL_BINARY")
        build_dir_raw = coalesce_cli_cfg(args.build_dir, cfg, "EOS_CASE_BUILD_DIR")
        metal_bin = _resolve_metal_binary(str(metal_binary_raw) if metal_binary_raw else None,
                                          str(build_dir_raw) if build_dir_raw else None)
        h5_dir = Path(str(coalesce_cli_cfg(args.h5_dir, cfg, "EOS_CASE_H5_DIR", args.h5_dir)))
        h5_dir.mkdir(parents=True, exist_ok=True)

        case_label = f"cr{cr_label}_sh{_fmt_label(sra_rate)}_{_fmt_label(lra_rate)}_z{z_label}"
        run_out_dir = h5_dir / case_label
        run_out_dir.mkdir(parents=True, exist_ok=True)

        env = os.environ.copy()
        env["METAL_OUTDIR"] = str(run_out_dir)
        env["METAL_ZETA0"] = f"{zeta0:.16g}"
        env["METAL_Z_METAL"] = f"{z_metal:.16g}"
        env["METAL_SRA_RATE"] = f"{sra_rate:.16g}"
        env["METAL_LRA_RATE"] = f"{lra_rate:.16g}"
        metal_max_iter_raw = coalesce_cli_cfg(args.metal_max_iter, cfg, "EOS_CASE_METAL_MAX_ITER")
        if metal_max_iter_raw is not None:
            env["METAL_MAX_ITER"] = str(metal_max_iter_raw)
        metal_xnh_stop_raw = coalesce_cli_cfg(args.metal_xnh_stop, cfg, "EOS_CASE_METAL_XNH_STOP")
        if metal_xnh_stop_raw is not None:
            env["METAL_XNH_STOP"] = str(metal_xnh_stop_raw)
        metal_cr_second_raw = coalesce_cli_cfg(args.metal_cr_second_frac, cfg, "EOS_CASE_METAL_CR_ATTEN_SECOND_FRAC")
        if metal_cr_second_raw is not None:
            env["METAL_CR_ATTEN_SECOND_FRAC"] = str(metal_cr_second_raw)
        metal_cr_bkgnd_raw = coalesce_cli_cfg(args.metal_cr_metal_bkgnd, cfg, "EOS_CASE_METAL_CR_METAL_BKGND")
        if metal_cr_bkgnd_raw is not None:
            env["METAL_CR_METAL_BKGND"] = str(metal_cr_bkgnd_raw)
        metal_abund_raw = coalesce_cli_cfg(args.metal_abundance_set, cfg, "EOS_CASE_METAL_ABUNDANCE_SET")
        if metal_abund_raw is not None:
            env["METAL_ABUNDANCE_SET"] = str(metal_abund_raw)

        metal_log_raw = coalesce_cli_cfg(args.metal_log, cfg, "EOS_CASE_METAL_LOG")
        print(f"[run] {' '.join([str(metal_bin)])}")
        if metal_log_raw:
            metal_log = Path(str(metal_log_raw))
            metal_log.parent.mkdir(parents=True, exist_ok=True)
            print(f"[run] metal_log={metal_log}")
            with metal_log.open("w", encoding="ascii") as lf:
                subprocess.run([str(metal_bin)], env=env, check=True, stdout=lf, stderr=subprocess.STDOUT)
        else:
            subprocess.run([str(metal_bin)], env=env, check=True)

        h5_files = sorted(run_out_dir.glob("collapse_*.h5"), key=lambda p: p.stat().st_mtime)
        if not h5_files:
            raise SystemExit(f"ERROR: metal_collapse did not produce HDF5 under {run_out_dir}")
        input_h5 = h5_files[-1]

    build_table(
        input_h5=input_h5,
        output_path=output_path,
        mode=str(coalesce_cli_cfg(args.mode, cfg, "EOS_CASE_MODE", args.mode)),
        eta_model=str(coalesce_cli_cfg(args.eta_model, cfg, "EOS_CASE_ETA_MODEL", args.eta_model)),
        eta_resample=str(coalesce_cli_cfg(args.eta_resample, cfg, "EOS_CASE_ETA_RESAMPLE", "native")),
        charge_table_path=Path(str(coalesce_cli_cfg(args.charge_table, cfg, "EOS_CASE_CHARGE_TABLE", args.charge_table)))
        if coalesce_cli_cfg(args.charge_table, cfg, "EOS_CASE_CHARGE_TABLE", args.charge_table) else None,
        mass_table_path=Path(str(coalesce_cli_cfg(args.mass_table, cfg, "EOS_CASE_MASS_TABLE", args.mass_table)))
        if coalesce_cli_cfg(args.mass_table, cfg, "EOS_CASE_MASS_TABLE", args.mass_table) else None,
        with_header=parse_bool(coalesce_cli_cfg(args.with_header, cfg, "EOS_CASE_WITH_HEADER", args.with_header)),
        n_log_min=float(coalesce_cli_cfg(args.n_log_min, cfg, "EOS_CASE_N_LOG_MIN", args.n_log_min)),
        n_log_max=float(coalesce_cli_cfg(args.n_log_max, cfg, "EOS_CASE_N_LOG_MAX", args.n_log_max)),
        n_log_step=float(coalesce_cli_cfg(args.n_log_step, cfg, "EOS_CASE_N_LOG_STEP", args.n_log_step)),
        b_log_min=float(coalesce_cli_cfg(args.b_log_min, cfg, "EOS_CASE_B_LOG_MIN", args.b_log_min)),
        b_log_max=float(coalesce_cli_cfg(args.b_log_max, cfg, "EOS_CASE_B_LOG_MAX", args.b_log_max)),
        b_log_step=float(coalesce_cli_cfg(args.b_log_step, cfg, "EOS_CASE_B_LOG_STEP", args.b_log_step)),
        etaa_compat_mode=str(coalesce_cli_cfg(args.etaa_compat_mode, cfg, "EOS_CASE_ETAA_COMPAT_MODE", "none")),
        etaa_drop_c2_logn_min=float(
            coalesce_cli_cfg(args.etaa_drop_c2_logn_min, cfg, "EOS_CASE_ETAA_DROP_C2_LOGN_MIN", 19.0)
        ),
    )

    print(f"[ok] input_h5={input_h5}")
    print(f"[ok] output_d={output_path}")
    print(f"[ok] params: zeta0={zeta0:.6e}, Z={z_metal:.6e}, sra={sra_rate:.6e}, lra={lra_rate:.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
