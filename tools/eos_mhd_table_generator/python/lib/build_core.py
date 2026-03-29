"""
Core EOS-table build implementation (no CLI parsing).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np

from .mapping import map89_to_64
from .resistivity_legacy import LegacyTables, eta_grid as eta_grid_legacy, load_legacy_tables

try:
    import h5py
except ImportError as e:  # pragma: no cover - import guard for CLI use
    raise RuntimeError("h5py is required for EOS table building") from e


KB_CGS = 1.380662e-16
LOG_FLOOR = 1.0e-300
C2_OVER_4PI_LOG10 = math.log10((2.99792458e10 ** 2) / (4.0 * math.pi))


def safe_log10(x: np.ndarray | float) -> np.ndarray | float:
    if isinstance(x, np.ndarray):
        return np.log10(np.maximum(x, LOG_FLOOR))
    return math.log10(max(float(x), LOG_FLOOR))


def _build_axis(vmin: float, vmax: float, step: float) -> np.ndarray:
    if step <= 0.0:
        raise ValueError("step must be > 0")
    n = int(round((vmax - vmin) / step))
    if n < 0:
        raise ValueError("vmax must be >= vmin")
    axis = vmin + step * np.arange(n + 1, dtype=float)
    if axis[-1] < vmax - 1e-12:
        axis = np.append(axis, vmax)
    return axis


def _load_h5(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"HDF5 not found: {path}")
    with h5py.File(path, "r") as f:
        if "xnH" not in f or "T_K" not in f or "y" not in f:
            raise KeyError("HDF5 must include datasets: xnH, T_K, y")
        data = {
            "xnH": np.asarray(f["xnH"][()], dtype=float),
            "T_K": np.asarray(f["T_K"][()], dtype=float),
            "y": np.asarray(f["y"][()], dtype=float),
        }
    if data["y"].ndim != 2:
        raise ValueError("dataset y must be 2D")
    if data["y"].shape[0] != data["xnH"].shape[0]:
        raise ValueError("row count mismatch between y and xnH")
    if data["xnH"].size == 0:
        raise ValueError(
            "input HDF5 has zero rows. Increase integration length "
            "(e.g., METAL_MAX_ITER) so collapse output contains samples."
        )
    return data


def _unique_sorted_xy(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    order = np.argsort(x)
    xs = x[order]
    ys = y[order]
    ux, idx = np.unique(xs, return_index=True)
    return ux, ys[idx]


def _interp_logy(x_src: np.ndarray, y_src: np.ndarray, x_tgt: np.ndarray) -> np.ndarray:
    xs, ys = _unique_sorted_xy(x_src, y_src)
    ly = safe_log10(np.maximum(ys, LOG_FLOOR))
    return np.power(10.0, np.interp(x_tgt, xs, ly, left=ly[0], right=ly[-1]))


def _resample_state_to_log_n(data: dict, logn_axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    x_src = safe_log10(np.maximum(data["xnH"], LOG_FLOOR))
    t_res = _interp_logy(x_src, data["T_K"], logn_axis)

    y_src = data["y"]
    y_res = np.empty((logn_axis.size, y_src.shape[1]), dtype=float)
    for j in range(y_src.shape[1]):
        y_res[:, j] = _interp_logy(x_src, y_src[:, j], logn_axis)
    return t_res, y_res


def _compute_mu_and_p(logn: np.ndarray, t: np.ndarray, y64: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    y_he = y64[:, 7]
    denom = y64[:, 0] + y64[:, 1] + y64[:, 2] + y64[:, 3] + y64[:, 7] + y64[:, 8] + y64[:, 9]
    denom = np.maximum(denom, LOG_FLOOR)
    mu = (1.0 + 4.0 * y_he) / denom
    n_h = np.power(10.0, logn)
    p = n_h * KB_CGS * t * (1.0 + 4.0 * y_he) / np.maximum(mu, LOG_FLOOR)
    return mu, p


def _eta_placeholder(logn: np.ndarray, t: np.ndarray, y64: np.ndarray, logb: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_h = np.power(10.0, logn)
    b = np.power(10.0, logb)
    y_e = np.maximum(y64[:, 2], 1.0e-20)

    eta_ohm = 1.0e18 * np.sqrt(np.maximum(t, 1.0)) / np.maximum(n_h * y_e, 1.0e-30)
    eta_amb = 1.0e19 * np.abs(b) / np.maximum(np.sqrt(n_h), 1.0e-20)
    eta_hal = 1.0e17 * np.abs(b) / np.maximum(np.power(n_h, 0.2), 1.0e-20)
    return eta_ohm, eta_amb, eta_hal


def _apply_etaa_compat(logn: np.ndarray, eta_a: np.ndarray, mode: str, drop_c2_logn_min: float) -> np.ndarray:
    if mode == "none":
        return eta_a
    if mode == "drop_c2_high_n":
        out = np.array(eta_a, copy=True)
        m = logn >= float(drop_c2_logn_min)
        out[m] = out[m] / (10.0 ** C2_OVER_4PI_LOG10)
        return out
    raise ValueError(f"unknown etaa_compat_mode: {mode}")


def _jj_label(mode: str, block_index_zero_based: int) -> int:
    if mode == "normalized":
        return block_index_zero_based + 1
    if mode == "compat":
        return max(1, block_index_zero_based)
    raise ValueError(f"unknown mode: {mode}")


def _iter_rows(
    mode: str,
    logn_axis: np.ndarray,
    logb_axis: np.ndarray,
    t_res: np.ndarray,
    y64_res: np.ndarray,
    x_src_native: np.ndarray | None,
    t_native: np.ndarray | None,
    y64_native: np.ndarray | None,
    eta_model: str,
    eta_resample: str,
    legacy_tables: LegacyTables | None,
    etaa_compat_mode: str,
    etaa_drop_c2_logn_min: float,
) -> Iterable[str]:
    mu, p = _compute_mu_and_p(logn_axis, t_res, y64_res)
    n_n = logn_axis.size
    for kb, logb in enumerate(logb_axis):
        jj = _jj_label(mode, kb)
        if eta_model == "placeholder":
            eta_o, eta_a, eta_h = _eta_placeholder(logn_axis, t_res, y64_res, float(logb))
        elif eta_model == "legacy":
            if legacy_tables is None:
                raise ValueError("legacy eta model requires loaded tables")
            if eta_resample == "native":
                if x_src_native is None or t_native is None or y64_native is None:
                    raise ValueError("eta_resample=native requires native state arrays")
                eo_raw, ea_raw, eh_raw = eta_grid_legacy(
                    x_src_native, t_native, y64_native, 10.0 ** float(logb), legacy_tables
                )
                eta_o = _interp_logy(x_src_native, eo_raw, logn_axis)
                eta_a = _interp_logy(x_src_native, ea_raw, logn_axis)
                eta_h = _interp_logy(x_src_native, eh_raw, logn_axis)
            elif eta_resample == "state":
                eta_o, eta_a, eta_h = eta_grid_legacy(logn_axis, t_res, y64_res, 10.0 ** float(logb), legacy_tables)
            else:
                raise ValueError(f"unknown eta_resample: {eta_resample}")
        else:
            raise ValueError(f"unknown eta model: {eta_model}")
        eta_a = _apply_etaa_compat(logn_axis, eta_a, etaa_compat_mode, etaa_drop_c2_logn_min)
        for i in range(n_n):
            ii = i + 1
            row = (
                f"{ii:d} {jj:d} "
                f"{logn_axis[i]:.7e} {safe_log10(t_res[i]):.7e} {safe_log10(mu[i]):.7e} {safe_log10(p[i]):.7e} "
                f"{logb:.7e} {safe_log10(eta_o[i]):.7e} {safe_log10(eta_a[i]):.7e} {safe_log10(eta_h[i]):.7e}"
            )
            yield row
        yield ""


def build_table(
    input_h5: Path,
    output_path: Path,
    mode: str,
    eta_model: str,
    eta_resample: str,
    charge_table_path: Path | None,
    mass_table_path: Path | None,
    with_header: bool,
    n_log_min: float,
    n_log_max: float,
    n_log_step: float,
    b_log_min: float,
    b_log_max: float,
    b_log_step: float,
    etaa_compat_mode: str,
    etaa_drop_c2_logn_min: float,
) -> None:
    data = _load_h5(input_h5)
    x_src_native = safe_log10(np.maximum(data["xnH"], LOG_FLOOR))
    logn_axis = _build_axis(n_log_min, n_log_max, n_log_step)
    logb_axis = _build_axis(b_log_min, b_log_max, b_log_step)

    t_res, y89_res = _resample_state_to_log_n(data, logn_axis)
    y64_res = map89_to_64(y89_res)
    t_native = np.asarray(data["T_K"], dtype=float)
    y64_native = map89_to_64(np.asarray(data["y"], dtype=float))
    legacy_tables: LegacyTables | None = None
    if eta_model == "legacy":
        if charge_table_path is None or mass_table_path is None:
            raise ValueError("legacy eta model needs charge_table_path and mass_table_path")
        legacy_tables = load_legacy_tables(charge_table_path, mass_table_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="ascii") as f:
        if with_header:
            f.write(f"# mode={mode}\n")
            f.write(f"# eta_model={eta_model}\n")
            f.write(f"# eta_resample={eta_resample}\n")
            f.write(f"# etaa_compat_mode={etaa_compat_mode}\n")
            f.write(f"# etaa_drop_c2_logn_min={etaa_drop_c2_logn_min}\n")
            f.write("# columns: ii jj nt tt mut pt bt etaO etaA etaH\n")
            if eta_model == "placeholder":
                f.write("# note: eta columns are placeholder values in Phase-A\n")
        for line in _iter_rows(
            mode,
            logn_axis,
            logb_axis,
            t_res,
            y64_res,
            x_src_native,
            t_native,
            y64_native,
            eta_model,
            eta_resample,
            legacy_tables,
            etaa_compat_mode,
            etaa_drop_c2_logn_min,
        ):
            f.write(line + "\n")

