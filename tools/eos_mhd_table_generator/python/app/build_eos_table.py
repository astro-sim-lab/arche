#!/usr/bin/env python3
"""
CLI wrapper: build EOS-MHD 10-column .d table from collapse HDF5.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from lib.build_core import build_table
from lib.config_kv import read_kv_config


def _coalesce_typed(args: argparse.Namespace, cfg: dict[str, str], arg_name: str, cfg_key: str, fallback: Any) -> Any:
    v = getattr(args, arg_name)
    if v is not None:
        return v
    if cfg_key in cfg:
        raw = cfg[cfg_key]
        if isinstance(fallback, float):
            return float(raw)
        if isinstance(fallback, bool):
            s = raw.strip().lower()
            if s in ("1", "true", "yes", "on"):
                return True
            if s in ("0", "false", "no", "off"):
                return False
            raise ValueError(f"{cfg_key}: invalid bool value '{raw}'")
        return raw
    return fallback


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build EOS-MHD 10-column .d table from arche-dev HDF5.")
    p.add_argument("--config", help="Optional config file (KEY=VALUE). CLI overrides config.")
    p.add_argument("--input-h5", help="Path to collapse HDF5 (metal_grain; y shape must include 89 species).")
    p.add_argument("--output", help="Output .d path.")
    p.add_argument("--mode", choices=("compat", "normalized"), help="Output indexing mode.")
    p.add_argument("--eta-model", choices=("placeholder", "legacy"), help="Resistivity model.")
    p.add_argument(
        "--eta-resample",
        choices=("state", "native"),
        help="eta resampling strategy: state=interpolate state then eta; native=eta on native rows then interpolate eta.",
    )
    p.add_argument("--charge-table", help="Path to legacy species charge table (needed for --eta-model legacy).")
    p.add_argument("--mass-table", help="Path to legacy species mass(amu) table (needed for --eta-model legacy).")
    p.add_argument("--q-path", dest="charge_table", help=argparse.SUPPRESS)
    p.add_argument("--mx-path", dest="mass_table", help=argparse.SUPPRESS)
    p.add_argument("--with-header", action="store_true", help="Write comment header lines at file top.")
    p.add_argument("--etaa-compat-mode", choices=("none", "drop_c2_high_n"), help="Optional etaA compatibility post-correction.")
    p.add_argument("--etaa-drop-c2-logn-min", type=float, help="Apply drop_c2_high_n only for log10(nH) >= this threshold.")
    p.add_argument("--n-log-min", type=float)
    p.add_argument("--n-log-max", type=float)
    p.add_argument("--n-log-step", type=float)
    p.add_argument("--b-log-min", type=float)
    p.add_argument("--b-log-max", type=float)
    p.add_argument("--b-log-step", type=float)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    try:
        cfg: dict[str, str] = {}
        if args.config:
            cfg = read_kv_config(Path(args.config))

        input_h5 = _coalesce_typed(args, cfg, "input_h5", "EOS_INPUT_H5", None)
        output = _coalesce_typed(args, cfg, "output", "EOS_OUTPUT", None)
        mode = _coalesce_typed(args, cfg, "mode", "EOS_MODE", "compat")
        eta_model = _coalesce_typed(args, cfg, "eta_model", "EOS_ETA_MODEL", "placeholder")
        eta_resample = _coalesce_typed(args, cfg, "eta_resample", "EOS_ETA_RESAMPLE", None)
        if eta_resample is None:
            eta_resample = "native" if str(eta_model) == "legacy" else "state"

        charge_table = _coalesce_typed(args, cfg, "charge_table", "EOS_CHARGE_TABLE_PATH", None)
        mass_table = _coalesce_typed(args, cfg, "mass_table", "EOS_MASS_TABLE_PATH", None)
        if charge_table is None:
            charge_table = _coalesce_typed(args, cfg, "charge_table", "EOS_Q_PATH", None)
        if mass_table is None:
            mass_table = _coalesce_typed(args, cfg, "mass_table", "EOS_MX_PATH", None)

        with_header = bool(_coalesce_typed(args, cfg, "with_header", "EOS_WITH_HEADER", False))
        n_log_min = _coalesce_typed(args, cfg, "n_log_min", "EOS_N_LOG_MIN", 0.1)
        n_log_max = _coalesce_typed(args, cfg, "n_log_max", "EOS_N_LOG_MAX", 22.0)
        n_log_step = _coalesce_typed(args, cfg, "n_log_step", "EOS_N_LOG_STEP", 0.1)
        b_log_min = _coalesce_typed(args, cfg, "b_log_min", "EOS_B_LOG_MIN", -20.0)
        b_log_max = _coalesce_typed(args, cfg, "b_log_max", "EOS_B_LOG_MAX", 4.0)
        b_log_step = _coalesce_typed(args, cfg, "b_log_step", "EOS_B_LOG_STEP", 0.2)
        etaa_compat_mode = _coalesce_typed(args, cfg, "etaa_compat_mode", "EOS_ETAA_COMPAT_MODE", "none")
        etaa_drop_c2_logn_min = _coalesce_typed(args, cfg, "etaa_drop_c2_logn_min", "EOS_ETAA_DROP_C2_LOGN_MIN", 19.0)

        if input_h5 is None or output is None:
            raise ValueError("input/output is required (use --input-h5/--output or EOS_INPUT_H5/EOS_OUTPUT in --config)")

        build_table(
            input_h5=Path(str(input_h5)),
            output_path=Path(str(output)),
            mode=str(mode),
            eta_model=str(eta_model),
            eta_resample=str(eta_resample),
            charge_table_path=Path(str(charge_table)) if charge_table is not None else None,
            mass_table_path=Path(str(mass_table)) if mass_table is not None else None,
            with_header=bool(with_header),
            n_log_min=float(n_log_min),
            n_log_max=float(n_log_max),
            n_log_step=float(n_log_step),
            b_log_min=float(b_log_min),
            b_log_max=float(b_log_max),
            b_log_step=float(b_log_step),
            etaa_compat_mode=str(etaa_compat_mode),
            etaa_drop_c2_logn_min=float(etaa_drop_c2_logn_min),
        )
    except Exception as e:  # noqa: BLE001
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
