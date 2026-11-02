#!/usr/bin/env python3
# Copyright (C) 2026 Shingo Hirano and Sho Higashi
# Licensed under the MIT found in the
# https://github.com/astro-sim-lab/arche/blob/main/LICENSE
"""Reproduction of Nakauchi et al. (2021, MNRAS 502, 3394) Figure 11.

4 rows x 3 columns of log10 Rm(lambda_J) in the (log n_H, log B) plane:
  rows    : metallicity Z/Zsun = 1e-6, 1e-4, 1e-2, 1
  columns : ionization rate zeta_ion = 0, 1e-17, 1e-15 s^-1
  colour  : log10 Rm, clipped to [-2, 10]
  contours: log10 Rm = 0 (thick) and 2, 4, 6, 8, 10 (thin)
  white   : B = B_crit = (4 pi G M_J rho / lambda_J)^(1/2)

Rm(lambda_J) = lambda_J^2 / (3 t_ff (eta_O + eta_A)) along the full
Nakauchi2021 one-zone collapse tracks.

Usage (from project root, after run_fig11.sh):
  python3 examples/nakauchi2021/plot_fig11.py [--h5dir DIR] [--save DIR] [--show]
"""

import argparse
import os

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ZTAGS = [
    ("1e-6", r"Z=10^{-6}Z_\odot"),
    ("1e-4", r"Z=10^{-4}Z_\odot"),
    ("1e-2", r"Z=10^{-2}Z_\odot"),
    ("1e+0", r"Z=1Z_\odot"),
]
ZETAS = [
    ("0", r"\zeta_{\rm ion}=0"),
    ("1e-17", r"\zeta_{\rm ion}=10^{-17}"),
    ("1e-15", r"\zeta_{\rm ion}=10^{-15}"),
]
ROWS = "abcd"

LOGRM_MIN, LOGRM_MAX = -2.0, 10.0
# Red (dissipating) through yellow and green to blue (well coupled).
CMAP = LinearSegmentedColormap.from_list(
    "rm",
    [(0.0, "#ff0000"), (2.0 / 12, "#ff9900"), (4.0 / 12, "#ffff00"),
     (6.0 / 12, "#33dd33"), (8.0 / 12, "#00a060"), (10.0 / 12, "#0060c0"),
     (1.0, "#0000ff")])


def load(h5dir, zeta, ztag):
    path = f"{h5dir}/rm_CR{zeta}_Z{ztag}.h5"
    if not os.path.exists(path):
        return None
    with h5py.File(path, "r") as f:
        return dict(nH=f["nH"][:], logB=f["log10_B"][:], Rm=f["Rm"][:],
                    B_cr=f["B_cr"][:])


def panel(ax, d, title):
    ax.set_xlim(0, 20)
    ax.set_ylim(-10, 10)
    ax.set_title(title, loc="left", fontsize=9)
    ax.set_xlabel(r"$\log\,n_{\rm H}\ ({\rm cm^{-3}})$")
    ax.set_ylabel(r"$\log\,B\ ({\rm G})$")
    if d is None:
        ax.text(10, 0, "missing", ha="center", va="center")
        return None
    x = np.log10(d["nH"])
    keep = (x >= -0.5) & (x <= 20.5)
    x = x[keep]
    with np.errstate(divide="ignore", invalid="ignore"):
        logrm = np.log10(d["Rm"][keep]).T  # (N_B, N)
    img = ax.pcolormesh(x, d["logB"], np.clip(logrm, LOGRM_MIN, LOGRM_MAX),
                        cmap=CMAP, vmin=LOGRM_MIN, vmax=LOGRM_MAX,
                        shading="nearest", rasterized=True)
    ax.contour(x, d["logB"], logrm, levels=[0.0], colors="k",
               linewidths=1.8)
    ax.contour(x, d["logB"], logrm, levels=[2, 4, 6, 8, 10], colors="k",
               linewidths=0.6)
    ax.plot(x, np.log10(d["B_cr"][keep]), color="white", lw=1.5)
    return img


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5dir", default="results/nakauchi2021",
                    help="Directory containing rm_CR*_Z*.h5")
    ap.add_argument("--save", default="results/nakauchi2021",
                    help="Directory for fig11_reproduction.png")
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    fig, axes = plt.subplots(4, 3, figsize=(13, 14.5),
                             constrained_layout=True)
    for r, (ztag, zlab) in enumerate(ZTAGS):
        for c, (zeta, clab) in enumerate(ZETAS):
            ax = axes[r, c]
            title = f"({ROWS[r]}{c + 1}) ${clab},\\ {zlab}$"
            img = panel(ax, load(args.h5dir, zeta, ztag), title)
            if img is not None:
                cb = fig.colorbar(img, ax=ax, pad=0.01)
                cb.set_label(r"$\log\,{\rm Rm}$")
    os.makedirs(args.save, exist_ok=True)
    out = os.path.join(args.save, "fig11_reproduction.png")
    fig.savefig(out, dpi=130)
    print(f"saved -> {out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
