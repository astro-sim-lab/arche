#!/usr/bin/env python3
"""
plot_cr0_zscan.py

Plot temperature and abundance profiles for CR=0 metal_grain collapse
across multiple metallicities.

Panels:
  [0] T vs nH          (y: 1e0–1e4, legend lower right)
  [1] y(e-)   vs nH    (y: 1e-15–1e-3)
  [2] y(H2) + y(HD) vs nH  (merged; H2=solid, HD=dashed)
  [3] y(Gr)   vs nH    (y: 1e-20–auto)

Species indices (0-based, metal_grain N_sp=89):
  e-  = 2
  H2  = 1
  HD  = 12
  Gr- = 66

Usage:
  python3 tools/python/plot_cr0_zscan.py [--h5dir DIR] [--save DIR] [--show]
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import h5py

# ── Species indices (0-based) ────────────────────────────────────────────────
SP = {
    "e-": 2,
    "H2": 1,
    "HD": 12,
    "Gr": 66,   # Gr^- (grain-)
}

# ── Z values and labels ──────────────────────────────────────────────────────
Z_VALUES = [
    ("1e-6", r"$Z = 10^{-6}\,Z_\odot$"),
    ("1e-5", r"$Z = 10^{-5}\,Z_\odot$"),
    ("1e-4", r"$Z = 10^{-4}\,Z_\odot$"),
    ("1e-3", r"$Z = 10^{-3}\,Z_\odot$"),
    ("1e-2", r"$Z = 10^{-2}\,Z_\odot$"),
    ("1e-1", r"$Z = 10^{-1}\,Z_\odot$"),
    ("1e+0", r"$Z = 1\,Z_\odot$"),
]

# viridis: high value = yellow (Z=1e-6), low value = dark purple/black (Z=1)
COLORS = plt.cm.viridis(np.linspace(0.95, 0.05, len(Z_VALUES)))


def load_h5(path):
    with h5py.File(path, "r") as f:
        xnH = f["xnH"][:]
        T_K = f["T_K"][:]
        y   = f["y"][:]   # (N, 89)
    return xnH, T_K, y


def main():
    parser = argparse.ArgumentParser(
        description="Plot CR=0 metal_grain collapse profiles across metallicities"
    )
    parser.add_argument("--h5dir", default="results/cr0_zscan/metal_grain",
                        help="Directory containing collapse_CR0_Z*.h5")
    parser.add_argument("--save", default="results/cr0_zscan",
                        help="Output directory for PNG")
    parser.add_argument("--show", action="store_true",
                        help="Show figure interactively")
    args = parser.parse_args()

    fig, axes = plt.subplots(4, 1, figsize=(6, 12), sharex=True)
    ax_T, ax_e, ax_mol, ax_Gr = axes

    for (z_tag, z_label), color in zip(Z_VALUES, COLORS):
        fname = os.path.join(args.h5dir, f"collapse_CR0_Z{z_tag}.h5")
        if not os.path.exists(fname):
            print(f"  [skip] not found: {fname}")
            continue

        xnH, T_K, y = load_h5(fname)

        ax_T.plot(  xnH, T_K,             color=color, lw=1.2, label=z_label)
        ax_e.plot(  xnH, y[:, SP["e-"]], color=color, lw=1.2)
        ax_mol.plot(xnH, y[:, SP["H2"]], color=color, lw=1.2, ls="-")
        ax_mol.plot(xnH, y[:, SP["HD"]], color=color, lw=1.2, ls="--")
        ax_Gr.plot( xnH, y[:, SP["Gr"]], color=color, lw=1.2)

    # ── Axes formatting ──────────────────────────────────────────────────────
    ax_T.set_ylabel(r"$T\,[\mathrm{K}]$")
    ax_T.set_yscale("log")
    ax_T.set_ylim(1e0, 1e4)

    ax_e.set_ylabel(r"$y(e^-)$")
    ax_e.set_yscale("log")
    ax_e.set_ylim(1e-15, 1e-3)

    ax_mol.set_ylabel(r"$y$")
    ax_mol.set_yscale("log")
    ax_mol.set_ylim(1e-8, 1e0)

    ax_Gr.set_ylabel(r"$y(\mathrm{Gr}^-)$")
    ax_Gr.set_yscale("log")
    ax_Gr.set_ylim(1e-18, 1e-8)

    ax_Gr.set_xlabel(r"$n_\mathrm{H}\,[\mathrm{cm}^{-3}]$")
    ax_Gr.set_xscale("log")
    ax_Gr.set_xlim(1e-1, 1e20)

    for ax in axes:
        ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

    # Top panel: Z legend (lower right)
    ax_T.legend(loc="lower right", fontsize=7, ncol=1,
                title=r"$\zeta_0 = 0$", title_fontsize=7)
    ax_T.set_title(r"metal\_grain collapse, CR = 0", fontsize=10)

    # Molecule panel: linestyle legend for H2 / HD
    h2_line = mlines.Line2D([], [], color="gray", ls="-",  lw=1.2, label=r"$\mathrm{H_2}$")
    hd_line = mlines.Line2D([], [], color="gray", ls="--", lw=1.2, label=r"$\mathrm{HD}$")
    ax_mol.legend(handles=[h2_line, hd_line], loc="lower right", fontsize=8)

    fig.tight_layout()

    if args.show:
        plt.show()
    else:
        os.makedirs(args.save, exist_ok=True)
        out_path = os.path.join(args.save, "fig_cr0_zscan.png")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
