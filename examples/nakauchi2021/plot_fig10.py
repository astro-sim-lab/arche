#!/usr/bin/env python3
# Copyright (C) 2026 Shingo Hirano and Sho Higashi
# Licensed under the MIT found in the
# https://github.com/astro-sim-lab/arche/blob/main/LICENSE
"""Reproduction of Nakauchi et al. (2021, MNRAS 502, 3394) Figure 10.

4 rows x 3 columns:
  rows    : (a) T_gas, (b) y(H2) & y(HD), (c) y(e-), (d) y(gr-)
  columns : ionization rate zeta_ion = 0, 1e-17, 1e-15 s^-1
  curves  : metallicity Z/Zsun = 1 ... 1e-6 (colour)
  style   : reduced metal Minimal (thin solid line) over full Nakauchi2021
            (thick dashed line) — overlap = agreement at a glance

This is the arche analogue of the paper's reduced-vs-full validation
figure, run as one-zone gravitational collapse up to n_H = 1e23 cm^-3.

Usage (from project root, after run_fig10.sh):
  python3 examples/nakauchi2021/plot_fig10.py [--h5dir DIR] [--save DIR] [--show]
"""

import argparse
import os

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Metallicity tag -> (Z value, colour).  Ordered high -> low like the paper.
# Colours read off Nakauchi+2021 Fig.10 legend (gnuplot-style palette).
ZTAGS = [
    ("1e+0", 1.0),
    ("1e-1", 1e-1),
    ("1e-2", 1e-2),
    ("1e-3", 1e-3),
    ("1e-4", 1e-4),
    ("1e-5", 1e-5),
    ("1e-6", 1e-6),
]
ZCOL = {
    "1e+0": "#000000",  # black
    "1e-1": "#808080",  # gray
    "1e-2": "#0000ff",  # blue
    "1e-3": "#006400",  # dark green
    "1e-4": "#00e0e0",  # cyan
    "1e-5": "#00cc00",  # green
    "1e-6": "#ff9900",  # orange
}

ZETAS = [
    ("0", r"$\zeta_{\rm ion}=0$"),
    ("1e-17", r"$\zeta_{\rm ion}=10^{-17}\,{\rm s^{-1}}$"),
    ("1e-15", r"$\zeta_{\rm ion}=10^{-15}\,{\rm s^{-1}}$"),
]


def load(h5dir, zeta, ztag, net):
    suf = "_min" if net == "min" else ""
    path = f"{h5dir}/collapse_CR{zeta}_Z{ztag}{suf}.h5"
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
    ap.add_argument("--h5dir", default="results/nakauchi2021",
                    help="Directory containing collapse_CR*_Z*.h5")
    ap.add_argument("--save", default="results/nakauchi2021",
                    help="Output directory for the PNG")
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    fig, axes = plt.subplots(4, 3, figsize=(13.5, 14), sharex="col")

    # x-axis upper limit per row: extend to the library's validated density
    # ceiling n_H = 1e23 cm^-3 across all panels.
    XMAX = [23, 23, 23, 23]

    for col, (zeta, ztitle) in enumerate(ZETAS):
        for ztag, _zv in ZTAGS:
            c = ZCOL[ztag]
            dmin = load(args.h5dir, zeta, ztag, "min")
            dful = load(args.h5dir, zeta, ztag, "full")
            # Draw full as a thick dashed line first, then the reduced model as
            # a thin solid line on top: where they agree the solid line tracks
            # the dashes; where they diverge the solid line leaves the dashes.
            for d, lw, ls, alpha, z in ((dful, 2.6, "--", 0.9, 1),
                                        (dmin, 1.3, "-", 1.0, 2)):
                if d is None:
                    continue
                nH = d["nH"]
                kw = dict(color=c, lw=lw, ls=ls, alpha=alpha, zorder=z,
                          solid_capstyle="round")
                # (a) temperature
                axes[0, col].plot(nH, d["T"], **kw)
                # (b) H2 (upper) and HD (lower)
                axes[1, col].plot(nH, frac(d, "H2"), **kw)
                hd = frac(d, "HD")
                if hd is not None:
                    axes[1, col].plot(nH, hd, **kw)
                # (c) electron fraction
                axes[2, col].plot(nH, frac(d, "e-"), **kw)
                # (d) gr- fraction
                gr = frac(d, "Gr-")
                if gr is not None:
                    axes[3, col].plot(nH, gr, **kw)

        axes[0, col].set_title(f"(a{col + 1}) {ztitle}", fontsize=11)

    # Row scales / labels
    axes[0, 0].set_ylabel(r"$T_{\rm gas}\ [{\rm K}]$")
    for c in range(3):
        axes[0, c].set_yscale("log")
        axes[0, c].set_ylim(1, 1e4)
    axes[1, 0].set_ylabel(r"$\log\,y({\rm H_2}),\ y({\rm HD})$")
    axes[2, 0].set_ylabel(r"$\log\,y(e^-)$")
    axes[3, 0].set_ylabel(r"$\log\,y({\rm gr^-})$")

    for r in (1, 2, 3):
        for c in range(3):
            axes[r, c].set_yscale("log")
    for c in range(3):
        axes[1, c].set_ylim(1e-8, 1.0)
        axes[2, c].set_ylim(1e-14, 1e-4)
        axes[3, c].set_ylim(1e-18, 1e-8)

    for c in range(3):
        axes[3, c].set_xlabel(r"$\log\,n_{\rm H}\ [{\rm cm^{-3}}]$")
    for row in range(4):
        for c in range(3):
            ax = axes[row, c]
            ax.set_xscale("log")
            ax.set_xlim(1.0, float(10.0 ** XMAX[row]))
            ax.grid(alpha=0.25, which="major")

    # H2 / HD annotations
    for c in range(3):
        axes[1, c].text(0.6, 0.78, r"H$_2$", transform=axes[1, c].transAxes,
                        fontsize=10)
        axes[1, c].text(0.6, 0.22, "HD", transform=axes[1, c].transAxes,
                        fontsize=10)

    # Legends: metallicity (colour) + network (style)
    zhandles = [
        Line2D([0], [0], color=ZCOL[t], lw=2,
               label=(r"$1\,Z_\odot$" if t == "1e+0"
                      else rf"$10^{{{int(np.log10(zv))}}}\,Z_\odot$"))
        for t, zv in ZTAGS
    ]
    nethandles = [
        Line2D([0], [0], color="k", lw=1.3, ls="-", alpha=1.0,
               label="reduced (Minimal, 40 sp)"),
        Line2D([0], [0], color="k", lw=2.6, ls="--", alpha=0.9,
               label="full (Nakauchi2021, 89 sp)"),
    ]
    axes[0, 0].legend(handles=zhandles, fontsize=7.5, ncol=2,
                      loc="lower right", framealpha=0.9)
    axes[0, 2].legend(handles=nethandles, fontsize=8, loc="lower right",
                      framealpha=0.9)

    fig.suptitle(
        "Reproduction of Nakauchi et al. (2021) Figure 10 — arche metal "
        "network\nreduced metal Minimal (thin solid line) vs full Nakauchi2021 "
        r"(thick dashed line); one-zone collapse to "
        r"$n_{\rm H}=10^{23}\,{\rm cm^{-3}}$",
        fontsize=13,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    if args.show:
        plt.show()
    else:
        os.makedirs(args.save, exist_ok=True)
        out = os.path.join(args.save, "fig10_reproduction.png")
        fig.savefig(out, dpi=130)
        print("wrote", out)


if __name__ == "__main__":
    main()
