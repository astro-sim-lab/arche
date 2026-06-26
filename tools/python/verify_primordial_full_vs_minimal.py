#!/usr/bin/env python3
# Copyright (C) 2026 Shingo Hirano and Sho Higashi
# Licensed under the MIT found in the
# https://github.com/astro-sim-lab/arche/blob/main/LICENSE
"""Run and compare Nakauchi2019 and Nakauchi2019_Minimal collapse outputs.

This script is a verification helper for the current codebase. It can:

1. Build the two primordial collapse binaries if needed.
2. Run both models for three parameter cases (a baseline plus two extra cases,
   e.g. f_ret=3/cr=1e-17 and f_ret=1/cr=1e-15), six runs in total.
3. Load the resulting HDF5 outputs.
4. Plot a single comparison figure containing:
   - Temperature
   - Ionization fraction
   - Net cooling function
   - Abundances of all species included in the minimal model

The figure overlays all runs per panel:
  color    encodes the model  (full = blue, minimal = red)
  linestyle encodes the case  (baseline = solid, 2nd = dashed, 3rd = dotted)

Example
-------
python3 tools/python/verify_primordial_full_vs_minimal.py \
    --build-dir build \
    --out-dir results/verify/full_vs_minimal \
    --zeta0 0 --fret 1 \
    --zeta0-2 1e-17 --fret-2 3 \
    --zeta0-3 1e-15 --fret-3 1
"""

from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Sequence

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


FULL_SPECIES_23 = [
    "H",
    "H2",
    "e-",
    "H+",
    "H2+",
    "H3+",
    "H-",
    "He",
    "He+",
    "He++",
    "HeH+",
    "D",
    "HD",
    "D+",
    "HD+",
    "D-",
    "Li",
    "LiH",
    "Li+",
    "Li-",
    "LiH+",
    "Li++",
    "Li+++",
]

MINIMAL_SPECIES_15 = [
    "H",
    "H2",
    "e-",
    "H+",
    "H2+",
    "H3+",
    "H-",
    "He",
    "He+",
    "D",
    "HD",
    "D+",
    "Li",
    "LiH",
    "Li+",
]

Y_D_REF = 2.58e-5
Y_LI_REF = 4.65e-10


@dataclass
class CaseData:
    """Container for one collapse output."""

    path: Path
    nH: np.ndarray
    t_k: np.ndarray
    y: np.ndarray
    species: List[str]
    lambda_net: np.ndarray
    y_plus: np.ndarray
    network: str


@dataclass
class PlotCase:
    """One full+minimal pair plotted as a single overlaid case.

    color encodes the model (full vs minimal); linestyle encodes the case
    (baseline vs second parameter set), so each panel carries 4 curves.
    """

    full: CaseData
    mini: CaseData
    label: str  # e.g. "cr=0, f_ret=1"
    linestyle: str


def parse_args() -> argparse.Namespace:
    """Parse command-line options."""
    parser = argparse.ArgumentParser(
        description=(
            "Run Nakauchi2019 and Nakauchi2019_Minimal and generate "
            "a single multi-panel comparison figure."
        )
    )
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Project root. Relative paths are resolved from this path.",
    )
    parser.add_argument(
        "--build-dir",
        default="build",
        help="CMake build directory.",
    )
    parser.add_argument(
        "--out-dir",
        default="results/verify/full_vs_minimal",
        help="Directory where run outputs and logs are stored.",
    )
    parser.add_argument(
        "--figure",
        default="",
        help=("Output PNG path. Default: <out-dir>/fig/compare_full_vs_minimal.png"),
    )
    parser.add_argument(
        "--zeta0",
        type=float,
        default=0.0,
        help="Baseline-case CR ionization rate [s^-1], applied to full+minimal.",
    )
    parser.add_argument(
        "--fret",
        type=float,
        default=1.0,
        help="Baseline-case free-fall retardation factor (PRIM_FF_RET).",
    )
    parser.add_argument(
        "--zeta0-2",
        type=float,
        default=1.0e-17,
        help="Second-case CR ionization rate [s^-1] (overlaid on the figure).",
    )
    parser.add_argument(
        "--fret-2",
        type=float,
        default=3.0,
        help="Second-case free-fall retardation factor (PRIM_FF_RET).",
    )
    parser.add_argument(
        "--zeta0-3",
        type=float,
        default=1.0e-15,
        help="Third-case CR ionization rate [s^-1] (overlaid on the figure).",
    )
    parser.add_argument(
        "--fret-3",
        type=float,
        default=1.0,
        help="Third-case free-fall retardation factor (PRIM_FF_RET).",
    )
    parser.add_argument(
        "--xnh-stop",
        type=float,
        default=1.0e23,
        help="Stop density nH [cm^-3], applied to both runs.",
    )
    parser.add_argument(
        "--output-stride",
        type=int,
        default=10,
        help="HDF5 output stride, applied to both runs.",
    )
    parser.add_argument(
        "--chem-table",
        default="data/primordial.h5",
        help="Path to primordial chemistry HDF5 table.",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Skip CMake configure/build and use existing binaries.",
    )
    parser.add_argument(
        "--skip-run",
        action="store_true",
        help="Skip simulation and only plot from existing HDF5 files.",
    )
    parser.add_argument(
        "--full-h5",
        default="",
        help="Existing baseline full-model HDF5. Required with --skip-run.",
    )
    parser.add_argument(
        "--minimal-h5",
        default="",
        help="Existing baseline minimal-model HDF5. Required with --skip-run.",
    )
    parser.add_argument(
        "--full-h5-2",
        default="",
        help="Existing second-case full-model HDF5. Required with --skip-run.",
    )
    parser.add_argument(
        "--minimal-h5-2",
        default="",
        help="Existing second-case minimal-model HDF5. Required with --skip-run.",
    )
    parser.add_argument(
        "--full-h5-3",
        default="",
        help="Existing third-case full-model HDF5. Required with --skip-run.",
    )
    parser.add_argument(
        "--minimal-h5-3",
        default="",
        help="Existing third-case minimal-model HDF5. Required with --skip-run.",
    )
    return parser.parse_args()


def resolve_path(repo_root: Path, raw_path: str) -> Path:
    """Resolve raw_path against repo_root unless already absolute."""
    path = Path(raw_path)
    return path if path.is_absolute() else (repo_root / path)


def ensure_inside_repo(repo_root: Path, path: Path) -> None:
    """Reject output paths that point outside the repository."""
    try:
        path.resolve().relative_to(repo_root.resolve())
    except ValueError as exc:
        raise ValueError(f"Path must stay inside the repository: {path}") from exc


def run_and_log(
    cmd: Sequence[str],
    cwd: Path,
    env: Mapping[str, str],
    log_path: Path,
) -> None:
    """Run one command and write stdout/stderr to log_path."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as logf:
        logf.write("$ " + " ".join(cmd) + "\n")
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            env=dict(env),
            stdout=logf,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed. See log: {log_path}")


def configure_and_build(repo_root: Path, build_dir: Path, log_dir: Path) -> None:
    """Configure CMake if needed, then build both primordial targets."""
    env = os.environ.copy()
    if not (build_dir / "CMakeCache.txt").is_file():
        run_and_log(
            [
                "cmake",
                "-S",
                str(repo_root),
                "-B",
                str(build_dir),
            ],
            cwd=repo_root,
            env=env,
            log_path=log_dir / "cmake_configure.log",
        )

    run_and_log(
        [
            "cmake",
            "--build",
            str(build_dir),
            "--target",
            "prim_collapse",
            "arche_collapse_prim_minimal",
            "-j",
        ],
        cwd=repo_root,
        env=env,
        log_path=log_dir / "cmake_build.log",
    )


def newest_h5(out_dir: Path) -> Path:
    """Return newest collapse HDF5 path under out_dir."""
    cands = sorted(out_dir.glob("collapse_CR*.h5"), key=lambda p: p.stat().st_mtime)
    if not cands:
        raise FileNotFoundError(f"No collapse HDF5 found in: {out_dir}")
    return cands[-1]


def run_one_model(
    binary: Path,
    repo_root: Path,
    out_dir: Path,
    log_path: Path,
    zeta0: float,
    f_ret: float,
    xnh_stop: float,
    output_stride: int,
    chem_table: Path,
) -> Path:
    """Run one primordial binary and return its output HDF5 path."""
    out_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(
        {
            "PRIM_ZETA0": f"{zeta0:.16g}",
            "PRIM_FF_RET": f"{f_ret:.16g}",
            "PRIM_OUTDIR": str(out_dir),
            "PRIM_XNH_STOP": f"{xnh_stop:.16g}",
            "PRIM_OUTPUT_STRIDE": str(output_stride),
            "PRIM_CHEM_TABLE": str(chem_table),
        }
    )
    run_and_log([str(binary)], cwd=repo_root, env=env, log_path=log_path)
    return newest_h5(out_dir)


def decode_text(value: object) -> str:
    """Decode HDF5 attribute values to text."""
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def normalize_species(species: List[str], n_cols: int, network: str) -> List[str]:
    """Repair species names if HDF5 metadata length does not match y-shape."""
    if len(species) == n_cols:
        return species
    if n_cols == len(MINIMAL_SPECIES_15):
        return MINIMAL_SPECIES_15.copy()
    if n_cols == len(FULL_SPECIES_23):
        return FULL_SPECIES_23.copy()
    raise ValueError(
        "Could not infer species list for output "
        f"(network={network}, y columns={n_cols}, species attr len={len(species)})"
    )


def load_case(path: Path) -> CaseData:
    """Load one HDF5 collapse file into a CaseData object."""
    with h5py.File(path, "r") as f:
        n_h = np.asarray(f["nH"])  # type: ignore[index]
        t_k = np.asarray(f["T_K"])  # type: ignore[index]
        y = np.asarray(f["y"])  # type: ignore[index]
        lam = np.asarray(f["Lambda_net"])  # type: ignore[index]
        y_plus = np.asarray(f["y_plus"]) if "y_plus" in f else np.zeros_like(n_h)
        network = decode_text(f.attrs.get("network", ""))

        sp_raw = f["y"].attrs.get("species", b"")
        sp_text = decode_text(sp_raw)
        species = [s.strip() for s in sp_text.split(",") if s.strip()]

    species = normalize_species(species, y.shape[1], network)
    return CaseData(
        path=path,
        nH=n_h,
        t_k=t_k,
        y=y,
        species=species,
        lambda_net=lam,
        y_plus=y_plus,
        network=network,
    )


def species_series(case: CaseData, name: str) -> np.ndarray:
    """Fetch abundance series y(name) from one case."""
    try:
        idx = case.species.index(name)
    except ValueError as exc:
        raise KeyError(f"Species '{name}' not found in {case.path}") from exc
    return case.y[:, idx]


def sum_available_species(case: CaseData, names: Sequence[str]) -> np.ndarray:
    """Sum abundances for species that exist in the case output."""
    acc = np.zeros(case.nH.shape[0], dtype=np.float64)
    for name in names:
        if name in case.species:
            acc += species_series(case, name)
    return acc


def at_lognh_idx(case: CaseData, target_log_nh: float) -> int:
    """Index of closest point to target log10(nH)."""
    return int(np.argmin(np.abs(np.log10(case.nH) - target_log_nh)))


def positive(arr: np.ndarray, floor: float = 1.0e-300) -> np.ndarray:
    """Return strictly positive series for log plotting."""
    return np.clip(np.asarray(arr, dtype=np.float64), floor, np.inf)


FULL_COLOR = "#1f77b4"  # model encoded by color: full = blue
MINI_COLOR = "#d62728"  # model encoded by color: minimal = red


def panel_series(case: CaseData, kind: str, sp: str = "") -> np.ndarray:
    """Return the y-series for one panel kind from one case."""
    if kind == "T":
        return case.t_k
    if kind == "ion":
        return case.y_plus if np.any(case.y_plus > 0.0) else species_series(case, "e-")
    if kind == "lam":
        return np.abs(case.lambda_net)
    return species_series(case, sp)


def plot_comparison(
    cases: Sequence[PlotCase],
    out_png: Path,
) -> List[Dict[str, float]]:
    """Overlay every PlotCase (full+minimal) on one multi-panel figure."""
    panels = [
        ("Temperature", "T", "", r"$T$ [K]"),
        ("Ionization Fraction", "ion", "", r"$x_{\mathrm{ion}}$"),
        (
            "Net Cooling Function",
            "lam",
            "",
            r"$|\Lambda_{\mathrm{net}}|$ [erg g$^{-1}$ s$^{-1}$]",
        ),
    ]
    for sp in MINIMAL_SPECIES_15:
        panels.append((f"Species: {sp}", "sp", sp, f"y({sp})"))

    n_cols = 5
    n_rows = int(math.ceil(len(panels) / n_cols))
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(3.3 * n_cols, 2.6 * n_rows),
        sharex=True,
    )
    ax_list = np.atleast_1d(axes).ravel()

    for idx, (title, kind, sp, ylabel) in enumerate(panels):
        ax = ax_list[idx]
        for case in cases:
            ax.plot(
                case.full.nH,
                positive(panel_series(case.full, kind, sp)),
                color=FULL_COLOR,
                lw=1.4,
                ls=case.linestyle,
            )
            ax.plot(
                case.mini.nH,
                positive(panel_series(case.mini, kind, sp)),
                color=MINI_COLOR,
                lw=1.2,
                ls=case.linestyle,
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1.0e-1, 1.0e23)
        ax.grid(True, which="both", alpha=0.3)
        ax.set_title(title, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=8)

    empty_idx = list(range(len(panels), len(ax_list)))
    for idx in empty_idx:
        ax_list[idx].axis("off")

    for ax in ax_list[-n_cols:]:
        ax.set_xlabel(r"$n_{\mathrm{H}}$ [cm$^{-3}$]", fontsize=9)

    # Legend: color = model, linestyle = case.  Render it in the empty panel
    # slots at the lower-right instead of over the first panel's curves.
    from matplotlib.lines import Line2D

    handles: List[Line2D] = []
    for case in cases:
        handles.append(
            Line2D([], [], color=FULL_COLOR, ls=case.linestyle,
                   label=f"full ({case.label})")
        )
        handles.append(
            Line2D([], [], color=MINI_COLOR, ls=case.linestyle,
                   label=f"minimal ({case.label})")
        )
    if empty_idx:
        ax_list[empty_idx[0]].legend(
            handles=handles, fontsize=10, loc="center", frameon=True
        )
    else:
        fig.legend(handles=handles, fontsize=10, loc="center right")

    metrics: List[Dict[str, float]] = []
    title_lines = ["Primordial collapse: Nakauchi2019 vs Nakauchi2019_Minimal"]
    for case in cases:
        full, mini = case.full, case.mini
        jf = at_lognh_idx(full, 14.0)
        jm = at_lognh_idx(mini, 14.0)
        ye_full = species_series(full, "e-")[jf]
        ye_mini = species_series(mini, "e-")[jm]
        ye_dex14 = float(np.log10(ye_mini / ye_full))

        d_full = sum_available_species(full, ["D", "HD", "D+", "HD+", "D-"])
        d_mini = sum_available_species(mini, ["D", "HD", "D+"])
        li_full = sum_available_species(
            full, ["Li", "LiH", "Li+", "Li-", "LiH+", "Li++", "Li+++"]
        )
        li_mini = sum_available_species(mini, ["Li", "LiH", "Li+"])

        title_lines.append(
            f"[{case.label}]  y_e residual @ log nH=14: {ye_dex14:+.4f} dex,  "
            f"D cons full/min = {d_full[jf] / Y_D_REF:.4f}/"
            f"{d_mini[jm] / Y_D_REF:.4f},  "
            f"Li cons full/min = {li_full[jf] / Y_LI_REF:.4f}/"
            f"{li_mini[jm] / Y_LI_REF:.4f}"
        )
        metrics.append(
            {
                "label": case.label,
                "ye_full_lognh14": float(ye_full),
                "ye_mini_lognh14": float(ye_mini),
                "ye_dex_lognh14": ye_dex14,
                "d_cons_full": float(d_full[jf] / Y_D_REF),
                "d_cons_mini": float(d_mini[jm] / Y_D_REF),
                "li_cons_full": float(li_full[jf] / Y_LI_REF),
                "li_cons_mini": float(li_mini[jm] / Y_LI_REF),
            }
        )

    fig.suptitle("\n".join(title_lines), fontsize=9)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return metrics


def main() -> int:
    """Entry point."""
    args = parse_args()

    repo_root = resolve_path(Path.cwd(), args.repo_root).resolve()
    build_dir = resolve_path(repo_root, args.build_dir).resolve()
    out_dir = resolve_path(repo_root, args.out_dir).resolve()
    chem_table = resolve_path(repo_root, args.chem_table).resolve()

    figure_path = (
        resolve_path(repo_root, args.figure).resolve()
        if args.figure
        else (out_dir / "fig" / "compare_full_vs_minimal.png").resolve()
    )

    ensure_inside_repo(repo_root, out_dir)
    ensure_inside_repo(repo_root, figure_path)

    log_dir = out_dir / "logs"

    def case_label(zeta0: float, fret: float) -> str:
        cr_str = "0" if zeta0 == 0.0 else f"{zeta0:.0e}"
        return f"cr={cr_str}, f_ret={fret:g}"

    # Overlaid cases: baseline (--zeta0/--fret), second (--zeta0-2/--fret-2),
    # third (--zeta0-3/--fret-3).
    # (subdir, zeta0, f_ret, linestyle, full_h5, minimal_h5)
    case_specs = [
        ("baseline", args.zeta0, args.fret, "-", args.full_h5, args.minimal_h5),
        ("case2", args.zeta0_2, args.fret_2, "--", args.full_h5_2, args.minimal_h5_2),
        ("case3", args.zeta0_3, args.fret_3, ":", args.full_h5_3, args.minimal_h5_3),
    ]

    if args.skip_run:
        missing = [
            n
            for n, v in (
                ("--full-h5", args.full_h5),
                ("--minimal-h5", args.minimal_h5),
                ("--full-h5-2", args.full_h5_2),
                ("--minimal-h5-2", args.minimal_h5_2),
                ("--full-h5-3", args.full_h5_3),
                ("--minimal-h5-3", args.minimal_h5_3),
            )
            if not v
        ]
        if missing:
            raise ValueError("--skip-run requires: " + ", ".join(missing))
    else:
        if not args.skip_build:
            configure_and_build(repo_root, build_dir, log_dir)

        full_bin = build_dir / "src/apps/collapse_primordial/prim_collapse"
        mini_bin = (
            build_dir / "src/apps/collapse_primordial/arche_collapse_prim_minimal"
        )
        if not full_bin.is_file():
            raise FileNotFoundError(f"Binary not found: {full_bin}")
        if not mini_bin.is_file():
            raise FileNotFoundError(f"Binary not found: {mini_bin}")
        if not chem_table.is_file():
            raise FileNotFoundError(f"Chem table not found: {chem_table}")

    cases: List[PlotCase] = []
    for sub, zeta0, fret, linestyle, full_h5, minimal_h5 in case_specs:
        if args.skip_run:
            full_h5_path = resolve_path(repo_root, full_h5).resolve()
            mini_h5_path = resolve_path(repo_root, minimal_h5).resolve()
        else:
            full_h5_path = run_one_model(
                binary=full_bin,
                repo_root=repo_root,
                out_dir=out_dir / f"full_{sub}",
                log_path=log_dir / f"run_full_{sub}.log",
                zeta0=zeta0,
                f_ret=fret,
                xnh_stop=args.xnh_stop,
                output_stride=args.output_stride,
                chem_table=chem_table,
            )
            mini_h5_path = run_one_model(
                binary=mini_bin,
                repo_root=repo_root,
                out_dir=out_dir / f"minimal_{sub}",
                log_path=log_dir / f"run_minimal_{sub}.log",
                zeta0=zeta0,
                f_ret=fret,
                xnh_stop=args.xnh_stop,
                output_stride=args.output_stride,
                chem_table=chem_table,
            )
        cases.append(
            PlotCase(
                full=load_case(full_h5_path),
                mini=load_case(mini_h5_path),
                label=case_label(zeta0, fret),
                linestyle=linestyle,
            )
        )

    metrics = plot_comparison(cases, figure_path)

    print("=== Verification Summary ===")
    print(f"figure      : {figure_path}")
    for case, m in zip(cases, metrics):
        print(f"--- [{m['label']}] ---")
        print(f"  full h5     : {case.full.path}")
        print(f"  minimal h5  : {case.mini.path}")
        print(
            "  y_e @ lognH=14: "
            f"full={m['ye_full_lognh14']:.4e}, "
            f"minimal={m['ye_mini_lognh14']:.4e}, "
            f"delta={m['ye_dex_lognh14']:+.4f} dex"
        )
        print(
            "  D conservation full/min: "
            f"{m['d_cons_full']:.4f}/{m['d_cons_mini']:.4f}"
        )
        print(
            "  Li conservation full/min: "
            f"{m['li_cons_full']:.4f}/{m['li_cons_mini']:.4f}"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
