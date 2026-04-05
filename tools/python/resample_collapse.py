#!/usr/bin/env python3
"""
resample_collapse.py — HDF5 collapse output → uniform log10(nH) bin CSV table

Reads one HDF5 file and resamples all physical quantities onto a uniform
log10(nH) grid by averaging data points within each bin.
Empty bins are filled with NaN.

Usage
-----
    python3 tools/resample_collapse.py [options]

Prim options
    --h5dir  DIR    Directory containing collapse_CR<tag>.h5
    --cr-tag TAG    CR tag (e.g. "1p5e-17", "0")
    --fret-tag TAG  free-fall retardation tag
                      scalar PRIM_FF_RET    → e.g. "3p0"   → suffix _fret3p0
                      table  PRIM_FRET_TABLE → use "step"  → suffix _fret-step
                      gamma  PRIM_FF_GAMMA   → use "gamma" → suffix _fret-gamma
                      (omit or '' for f_ret=1, no suffix)
    --jlw-tag  TAG  Lyman-Werner tag (e.g. "1p5" → suffix _JLW1p5; omit for J_LW21=0)
    --zred-tag TAG  redshift tag (omit or '' for z=0)
    --save   DIR    Output directory for CSV (default: same as --h5dir)

Metal options
    --metal-h5dir DIR   Directory containing collapse_CR<cr>_Z<z>.h5
    --metal-cr    TAG   CR tag for metal (e.g. "1p5e-17" or "0")
    --metal-z     TAG   Z tag for metal  (e.g. "1p2e-4" or "0")
    --metal-fret  TAG   free-fall retardation tag for metal ("step" for table mode)
    --metal-jlw   TAG   Lyman-Werner tag for metal (omit for J_LW21=0)
    --metal-zred  TAG   redshift tag for metal
    --metal-save  DIR   Output directory for CSV (default: same as --metal-h5dir)

Grid options
    --log-nH-min  FLOAT  log10(nH) lower edge (default: -2.0)
    --log-nH-max  FLOAT  log10(nH) upper edge (default: 23.0)
    --log-nH-step FLOAT  bin width in log10(nH) (default: 0.1)

Output columns (CSV)
    log10_nH      — bin centre in log10(nH [cm^-3])
    T_K           — gas temperature  [K]            (log-averaged)
    T_gr_K        — grain temperature [K]  (metal only, log-averaged)
    xLmbd_net     — net cooling [erg g^-1 s^-1]     (linear-averaged)
    y_0 … y_N     — species abundances               (log-averaged, NaN if all zero in bin)

Averaging conventions (pluggable — edit AVG_FUNCS dict to customise)
    log_avg  : 10^(mean(log10(x)))   for x > 0; NaN if all values are <= 0
    lin_avg  : mean(x)               for signed/zero quantities (e.g. xLmbd_net)
"""

import argparse
import os
import sys
from typing import Callable, Dict

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("ERROR: h5py is not installed.  Run: pip3 install h5py")

# ─────────────────────────────────────────────────────────────────────────────
# Averaging functions
# ─────────────────────────────────────────────────────────────────────────────

def log_avg(arr: np.ndarray) -> float:
    """Geometric mean: 10^mean(log10(x)) for positive entries only.
    Returns NaN if no positive values exist."""
    pos = arr[arr > 0.0]
    if pos.size == 0:
        return np.nan
    return 10.0 ** np.mean(np.log10(pos))


def lin_avg(arr: np.ndarray) -> float:
    """Arithmetic mean.  Returns NaN on empty input."""
    if arr.size == 0:
        return np.nan
    return float(np.mean(arr))


# ── Assign averaging functions to output columns ──────────────────────────────
# Keys: exact column names as they appear in the output CSV header.
# Fallback (not in dict) → log_avg.
AVG_FUNCS: Dict[str, Callable[[np.ndarray], float]] = {
    "xLmbd_net": lin_avg,
    # Add overrides here, e.g.:
    #   "xGam_cmp": lin_avg,
}


def _avg(col_name: str, arr: np.ndarray) -> float:
    fn = AVG_FUNCS.get(col_name, log_avg)
    return fn(arr)


# ─────────────────────────────────────────────────────────────────────────────
# HDF5 loading
# ─────────────────────────────────────────────────────────────────────────────

def load_h5(path: str) -> dict:
    """Load one HDF5 file into a dict.  Arrays are numpy; attrs are scalars/str."""
    if not os.path.isfile(path):
        sys.exit(
            f"ERROR: HDF5 file not found: {path}\n"
            f"       Check --h5dir / --cr-tag / --fret-tag / --zred-tag options.\n"
            f"       (Relative paths are resolved from the current working directory.)"
        )
    with h5py.File(path, "r") as f:
        d: dict = {}
        for k in f.attrs:
            v = f.attrs[k]
            d[k] = v.decode() if isinstance(v, bytes) else v
        for k in f:
            arr = f[k][()]
            if isinstance(arr, bytes):
                arr = arr.decode()
            d[k] = arr
        if "y" in f:
            sp_raw = f["y"].attrs.get("species", b"")
            if isinstance(sp_raw, bytes):
                sp_raw = sp_raw.decode()
            d["species"] = [s.strip() for s in sp_raw.split(",")]
    return d


# ─────────────────────────────────────────────────────────────────────────────
# Grid definition
# ─────────────────────────────────────────────────────────────────────────────

def make_grid(log_min: float, log_max: float, step: float):
    """Return (bin_edges, bin_centers) arrays.

    Edges span [log_min, log_max] with spacing *step*.
    Centers are midpoints: edges[i] + step/2.
    """
    edges = np.arange(log_min, log_max + 0.5 * step, step)
    centers = edges[:-1] + 0.5 * step
    return edges, centers


# ─────────────────────────────────────────────────────────────────────────────
# Core resampling
# ─────────────────────────────────────────────────────────────────────────────

def resample(data: dict,
             log_min: float = -2.0,
             log_max: float = 23.0,
             step: float    = 0.5,
             has_grain: bool = False) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Resample *data* onto a uniform log10(nH) grid.

    Returns
    -------
    centers  : (N_bins,) — log10(nH) bin centres
    colnames : list of column names (excluding log10_nH)
    table    : (N_bins, len(colnames)) — float64 array, NaN for empty bins
    """
    edges, centers = make_grid(log_min, log_max, step)
    N_bins = len(centers)

    log_nH_data = np.log10(data["xnH"])
    # Assign each data point to a bin
    bin_idx = np.digitize(log_nH_data, edges) - 1   # 0-based; -1 = below, N = above

    # ── Build ordered column list ────────────────────────────────────────────
    # Fixed columns (in order)
    scalar_cols = ["T_K"]
    if has_grain:
        scalar_cols.append("T_gr_K")
    scalar_cols.append("xLmbd_net")

    species_list = data.get("species", [])
    N_sp = data["y"].shape[1] if "y" in data else 0
    sp_cols = [f"y_{i}" for i in range(N_sp)]   # y_0, y_1, …

    colnames = scalar_cols + sp_cols
    N_col = len(colnames)
    table = np.full((N_bins, N_col), np.nan)

    # ── Iterate over bins ────────────────────────────────────────────────────
    for b in range(N_bins):
        mask = (bin_idx == b)
        if not np.any(mask):
            continue   # leave as NaN

        col_i = 0
        for col in scalar_cols:
            if col in data:
                table[b, col_i] = _avg(col, data[col][mask])
            col_i += 1

        if "y" in data:
            y_block = data["y"]   # shape (N_rows, N_sp)
            for j in range(N_sp):
                table[b, col_i] = _avg(sp_cols[j], y_block[mask, j])
                col_i += 1

    return centers, colnames, table


# ─────────────────────────────────────────────────────────────────────────────
# I/O helpers
# ─────────────────────────────────────────────────────────────────────────────

def _h5_stem(h5dir: str, stem: str) -> str:
    return os.path.join(h5dir, stem + ".h5")


def _fret_suffix(fret_tag: str) -> str:
    """Return the _fret<…> suffix for a given fret_tag.

    Named modes are passed without their leading dash (stripped in run_collapse.sh
    to avoid argparse treating leading-dash values as flags):
      "step"  → "_fret-step"   (PRIM/METAL_FRET_TABLE mode)
      "gamma" → "_fret-gamma"  (PRIM/METAL_FF_GAMMA mode; Higuchi+2018 Eq.5-7)
    Numeric tags are passed as-is:
      "3p0"   → "_fret3p0"
    """
    if not fret_tag:
        return ""
    _NAMED = {"step", "gamma"}
    return f"_fret-{fret_tag}" if fret_tag in _NAMED else f"_fret{fret_tag}"


def _jlw_suffix(jlw_tag: str) -> str:
    """Return the _JLW<…> suffix for a given jlw_tag (e.g. "1p5" → "_JLW1p5").
    Empty string → no suffix (J_LW21 = 0 case)."""
    return f"_JLW{jlw_tag}" if jlw_tag else ""


def _prim_stem(cr_tag: str, fret_tag: str = "", zred_tag: str = "",
               jlw_tag: str = "") -> str:
    s = f"collapse_CR{cr_tag}"
    s += _fret_suffix(fret_tag)
    s += _jlw_suffix(jlw_tag)
    if zred_tag:
        s += f"_z{zred_tag}"
    return s


def _metal_stem(cr_tag: str, z_tag: str,
                fret_tag: str = "", zred_tag: str = "",
                jlw_tag: str = "") -> str:
    s = f"collapse_CR{cr_tag}_Z{z_tag}"
    s += _fret_suffix(fret_tag)
    s += _jlw_suffix(jlw_tag)
    if zred_tag:
        s += f"_z{zred_tag}"
    return s


def write_csv(out_path: str, centers: np.ndarray,
              colnames: list, table: np.ndarray,
              species_list: list) -> None:
    """Write resampled table to CSV.

    Column order: log10_nH, <scalar_cols>, y_0, y_1, …
    Species names are annotated in the header comment.
    """
    header_lines = ["# resample_collapse.py output"]
    if species_list:
        sp_str = ",".join(species_list)
        header_lines.append(f"# species: {sp_str}")
    header_lines.append("# NaN = empty bin (no data points)")

    all_colnames = ["log10_nH"] + colnames
    header_lines.append(",".join(all_colnames))

    with open(out_path, "w", encoding="utf-8") as fout:
        for line in header_lines:
            fout.write(line + "\n")

        for b in range(len(centers)):
            row_vals = [f"{centers[b]:.6f}"]
            for j in range(len(colnames)):
                v = table[b, j]
                row_vals.append("NaN" if np.isnan(v) else f"{v:.6e}")
            fout.write(",".join(row_vals) + "\n")

    print(f"  saved: {out_path}  ({len(centers)} bins, "
          f"{int(np.sum(~np.isnan(table[:, 0])))}/{len(centers)} non-empty)")


# ─────────────────────────────────────────────────────────────────────────────
# Per-case processing
# ─────────────────────────────────────────────────────────────────────────────

def process_prim(args) -> None:
    if not (args.h5dir and args.cr_tag):
        return

    stem    = _prim_stem(args.cr_tag, args.fret_tag or "", args.zred_tag or "",
                         args.jlw_tag or "")
    h5path  = _h5_stem(args.h5dir, stem)
    save_dir = args.save if args.save else args.h5dir
    os.makedirs(save_dir, exist_ok=True)

    print(f"[prim] loading {h5path}")
    data = load_h5(h5path)
    centers, colnames, table = resample(
        data, args.log_nH_min, args.log_nH_max, args.log_nH_step,
        has_grain=False)

    out_path = os.path.join(save_dir, f"resample_{stem}.csv")
    write_csv(out_path, centers, colnames, table,
              data.get("species", []))


def process_metal(args) -> None:
    if not (args.metal_h5dir and args.metal_cr and args.metal_z):
        return

    stem     = _metal_stem(args.metal_cr, args.metal_z,
                           args.metal_fret or "", args.metal_zred or "",
                           args.metal_jlw or "")
    h5path   = _h5_stem(args.metal_h5dir, stem)
    save_dir = args.metal_save if args.metal_save else args.metal_h5dir
    os.makedirs(save_dir, exist_ok=True)

    print(f"[metal] loading {h5path}")
    data = load_h5(h5path)
    centers, colnames, table = resample(
        data, args.log_nH_min, args.log_nH_max, args.log_nH_step,
        has_grain=True)

    out_path = os.path.join(save_dir, f"resample_{stem}.csv")
    write_csv(out_path, centers, colnames, table,
              data.get("species", []))


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Resample HDF5 collapse output onto uniform log10(nH) bins.")

    # ── prim ────────────────────────────────────────────────────────────────
    pg = p.add_argument_group("Prim")
    pg.add_argument("--h5dir",    metavar="DIR", help="Directory with prim HDF5 files")
    pg.add_argument("--cr-tag",   metavar="TAG", help="CR tag  (e.g. '1p5e-17')")
    pg.add_argument("--fret-tag", metavar="TAG", default="",
                    help="f_ret tag: scalar → e.g. '3p0'; table → 'step' (→ _fret-step)")
    pg.add_argument("--jlw-tag",  metavar="TAG", default="",
                    help="Lyman-Werner tag (e.g. '1p5' → _JLW1p5; omit for J_LW21=0)")
    pg.add_argument("--zred-tag", metavar="TAG", default="",
                    help="Redshift tag (omit for z=0)")
    pg.add_argument("--save",     metavar="DIR",
                    help="Output directory for prim CSV (default: same as --h5dir)")

    # ── metal ────────────────────────────────────────────────────────────────
    mg = p.add_argument_group("Metal")
    mg.add_argument("--metal-h5dir", metavar="DIR", help="Directory with metal HDF5 files")
    mg.add_argument("--metal-cr",    metavar="TAG", help="CR tag for metal")
    mg.add_argument("--metal-z",     metavar="TAG", help="Z tag for metal (e.g. '1p2e-4')")
    mg.add_argument("--metal-fret",  metavar="TAG", default="",
                    help="f_ret tag for metal: scalar → e.g. '3p0'; table → 'step'")
    mg.add_argument("--metal-jlw",   metavar="TAG", default="",
                    help="Lyman-Werner tag for metal (e.g. '1p5' → _JLW1p5)")
    mg.add_argument("--metal-zred",  metavar="TAG", default="",
                    help="Redshift tag for metal")
    mg.add_argument("--metal-save",  metavar="DIR",
                    help="Output directory for metal CSV (default: same as --metal-h5dir)")

    # ── grid ─────────────────────────────────────────────────────────────────
    gg = p.add_argument_group("Grid")
    gg.add_argument("--log-nH-min",  type=float, default=-2.0,
                    metavar="FLOAT", help="log10(nH) lower edge (default: -2.0)")
    gg.add_argument("--log-nH-max",  type=float, default=23.0,
                    metavar="FLOAT", help="log10(nH) upper edge (default: 23.0)")
    gg.add_argument("--log-nH-step", type=float, default=0.1,
                    metavar="FLOAT", help="Bin width in log10(nH) (default: 0.5)")

    return p


def main() -> None:
    parser = _build_parser()
    args   = parser.parse_args()

    ran = False

    if args.h5dir and args.cr_tag:
        process_prim(args)
        ran = True

    if args.metal_h5dir and args.metal_cr and args.metal_z:
        process_metal(args)
        ran = True

    if not ran:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
