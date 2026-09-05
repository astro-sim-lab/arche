#!/usr/bin/env python3
# Copyright (C) 2026 Shingo Hirano and Sho Higashi
# Licensed under the MIT found in the
# https://github.com/astro-sim-lab/arche/blob/main/LICENSE
"""Reduced-vs-full primordial collapse comparison (Nakauchi et al. 2019 style).

4 rows x 3 columns:
  rows    : (a) T_gas
            (b) y(H2) & y(HD)
            (c) Major species carried by the Minimal network
            (d) Minor species carried by the Minimal network
  columns : ionization rate zeta_ion = 0, 1e-17, 1e-15 s^-1
  style   : reduced Nakauchi2019_Minimal (thin solid line, 15 sp) over full
            Nakauchi2019 (thick dashed line, 23 sp) — overlap = agreement

Unlike the metal Fig.10 reproduction (which colours curves by metallicity),
primordial chemistry has no metallicity axis, so rows (c)/(d) instead colour
curves by chemical species — the arche analogue of replacing the paper's
log y(e-) / log y(gr-) panels with a major/minor species breakdown.  Every
species the Minimal network carries is shown, so the figure doubles as a
species-by-species reduced-vs-full validation.

Usage (from project root, after run_prim.sh):
  python3 examples/nakauchi2019/plot_prim.py [--h5dir DIR] [--save DIR] [--show]
"""

import argparse
import os

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ZETAS = [
    ("0", r"$\zeta_{\rm ion}=0$"),
    ("1e-17", r"$\zeta_{\rm ion}=10^{-17}\,{\rm s^{-1}}$"),
    ("1e-15", r"$\zeta_{\rm ion}=10^{-15}\,{\rm s^{-1}}$"),
]

# Species the Minimal network carries (15), split by typical abundance.
# H2 and HD live in row (b); the remaining 13 are split major/minor here.
MAJOR = ["H", "He", "e-", "H+", "D"]
MINOR = ["H-", "H2+", "H3+", "He+", "D+", "Li", "LiH", "Li+"]

# Stable per-species colours (tab10 / tab20 picks).
SPCOL = {
    "H": "#1f77b4", "He": "#ff7f0e", "e-": "#2ca02c",
    "H+": "#d62728", "D": "#9467bd",
    "H-": "#8c564b", "H2+": "#e377c2", "H3+": "#7f7f7f",
    "He+": "#bcbd22", "D+": "#17becf", "Li": "#aec7e8",
    "LiH": "#ffbb78", "Li+": "#c5b0d5",
}
# Display labels (matplotlib mathtext).
SPLAB = {
    "H": r"H", "He": r"He", "e-": r"$e^-$", "H+": r"H$^+$", "D": r"D",
    "H-": r"H$^-$", "H2+": r"H$_2^+$", "H3+": r"H$_3^+$", "He+": r"He$^+$",
    "D+": r"D$^+$", "Li": r"Li", "LiH": r"LiH", "Li+": r"Li$^+$",
}
H2COL, HDCOL = "#1f77b4", "#d62728"


def load(h5dir, zeta, net):
    suf = "_min" if net == "min" else ""
    path = f"{h5dir}/collapse_CR{zeta}{suf}.h5"
    if not os.path.exists(path):
        return None
    f = h5py.File(path, "r")
    sp = f["y"].attrs["species"].decode().split(",")
    idx = {s: i for i, s in enumerate(sp)}
    nH = f["nH"][:] if "nH" in f else f["xnH"][:]
    return dict(nH=nH, T=f["T_K"][:], y=f["y"][:], idx=idx)


def frac(d, name):
    return d["y"][:, d["idx"][name]] if name in d["idx"] else None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5dir", default="results/nakauchi2019",
                    help="Directory containing collapse_CR*.h5")
    ap.add_argument("--save", default="results/nakauchi2019",
                    help="Output directory for the PNG")
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    fig, axes = plt.subplots(4, 3, figsize=(13.5, 14), sharex="col")

    for col, (zeta, ztitle) in enumerate(ZETAS):
        dmin = load(args.h5dir, zeta, "min")
        dful = load(args.h5dir, zeta, "full")
        # Full = thick dashed line (drawn first); reduced = thin solid line on
        # top. Agreement shows as the thin solid line tracking the dashes;
        # divergence shows as the solid line leaving the dashes.
        for d, lw, ls, alpha, z in ((dful, 2.6, "--", 0.9, 1),
                                    (dmin, 1.4, "-", 1.0, 2)):
            if d is None:
                continue
            nH = d["nH"]
            kw = dict(lw=lw, ls=ls, alpha=alpha, zorder=z,
                      solid_capstyle="round")
            # (a) temperature — single black curve per network
            axes[0, col].plot(nH, d["T"], color="k", **kw)
            # (b) H2 / HD
            axes[1, col].plot(nH, frac(d, "H2"), color=H2COL, **kw)
            hd = frac(d, "HD")
            if hd is not None:
                axes[1, col].plot(nH, hd, color=HDCOL, **kw)
            # (c) major species  /  (d) minor species
            for sp in MAJOR:
                yv = frac(d, sp)
                if yv is not None:
                    axes[2, col].plot(nH, yv, color=SPCOL[sp], **kw)
            for sp in MINOR:
                yv = frac(d, sp)
                if yv is not None:
                    axes[3, col].plot(nH, yv, color=SPCOL[sp], **kw)

        axes[0, col].set_title(f"(a{col + 1}) {ztitle}", fontsize=11)

    # Row scales / labels
    axes[0, 0].set_ylabel(r"$T_{\rm gas}\ [{\rm K}]$")
    axes[1, 0].set_ylabel(r"$\log\,y({\rm H_2}),\ y({\rm HD})$")
    axes[2, 0].set_ylabel(r"$\log\,y\ ({\rm major})$")
    axes[3, 0].set_ylabel(r"$\log\,y\ ({\rm minor})$")

    for c in range(3):
        axes[0, c].set_yscale("log")
        axes[0, c].set_ylim(1, 1e4)
    for r in (1, 2, 3):
        for c in range(3):
            axes[r, c].set_yscale("log")
    for c in range(3):
        axes[1, c].set_ylim(1e-8, 1.0)
        axes[2, c].set_ylim(1e-13, 3.0)
        axes[3, c].set_ylim(1e-20, 1e-4)

    for c in range(3):
        axes[3, c].set_xlabel(r"$\log\,n_{\rm H}\ [{\rm cm^{-3}}]$")
    for row in range(4):
        for c in range(3):
            ax = axes[row, c]
            ax.set_xscale("log")
            ax.set_xlim(1.0, 1e23)
            ax.grid(alpha=0.25, which="major")

    # H2 / HD annotations
    for c in range(3):
        axes[1, c].text(0.62, 0.80, r"H$_2$", color=H2COL,
                        transform=axes[1, c].transAxes, fontsize=10)
        axes[1, c].text(0.62, 0.20, "HD", color=HDCOL,
                        transform=axes[1, c].transAxes, fontsize=10)

    # Legends
    nethandles = [
        Line2D([0], [0], color="k", lw=1.4, ls="-", alpha=1.0,
               label="reduced (Minimal, 15 sp)"),
        Line2D([0], [0], color="k", lw=2.6, ls="--", alpha=0.9,
               label="full (Nakauchi2019, 23 sp)"),
    ]
    axes[0, 0].legend(handles=nethandles, fontsize=8, loc="lower right",
                      framealpha=0.9)

    major_h = [Line2D([0], [0], color=SPCOL[s], lw=2, label=SPLAB[s])
               for s in MAJOR]
    minor_h = [Line2D([0], [0], color=SPCOL[s], lw=2, label=SPLAB[s])
               for s in MINOR]
    axes[2, 2].legend(handles=major_h, fontsize=8, ncol=2, loc="lower right",
                      framealpha=0.9)
    axes[3, 2].legend(handles=minor_h, fontsize=8, ncol=2, loc="lower right",
                      framealpha=0.9)

    fig.suptitle(
        "Primordial reduced-vs-full collapse — arche zero-metal network "
        "(Nakauchi et al. 2019 style)\nreduced Nakauchi2019_Minimal (thin solid "
        "line, 15 sp) vs full Nakauchi2019 (thick dashed line, 23 sp); "
        r"one-zone collapse to $n_{\rm H}=10^{23}\,{\rm cm^{-3}}$",
        fontsize=13,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    if args.show:
        plt.show()
    else:
        os.makedirs(args.save, exist_ok=True)
        out = os.path.join(args.save, "fig_prim_reduced_vs_full.png")
        fig.savefig(out, dpi=130)
        print("wrote", out)


if __name__ == "__main__":
    main()
