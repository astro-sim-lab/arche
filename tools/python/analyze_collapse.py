#!/usr/bin/env python3
"""
analyze_collapse.py — single-case analysis and visualization (HDF5 input)

Reads one HDF5 file per run and generates figures.
Figure filenames include the parameter tags, e.g.:
    fig1_phase_cooling_CR1p5e-17.png
    fig3_jeans_CR1p5e-17_Z1p2e-4.png

Usage
-----
    python3 tools/analyze_collapse.py [options]

Prim options
    --h5dir  DIR    Directory containing collapse_CR<tag>.h5
    --cr-tag TAG    CR tag (e.g. "1p5e-17", "0")
    --fret-tag TAG  free-fall retardation tag (omit or '' for f_ret=1)
    --jlw-tag  TAG  Lyman-Werner tag (omit or '' for J_LW21=0)
    --zred-tag TAG  redshift tag for prim (omit or '' for z=0)
    --save   DIR    Output directory for PNG files (default: results/prim/fig)

Metal options
    --metal-h5dir DIR   Directory containing collapse_CR<cr>_Z<z>.h5
    --metal-cr    TAG   CR tag for metal (e.g. "1p5e-17" or "0")
    --metal-z     TAG   Z tag for metal  (e.g. "1p2e-4" or "0")
    --metal-fret  TAG   free-fall retardation tag for metal (omit or '' for f_ret=1)
    --metal-jlw   TAG   Lyman-Werner tag for metal (omit or '' for J_LW21=0)
    --metal-zred  TAG   redshift tag for metal (omit or '' for z=0)
    --metal-save  DIR   Output directory for PNG files (default: results/metal/fig)

Common
    --show              Interactive display instead of saving
    --fig-combo         Also output a single-panel summary figure (fig_summary_<tag>.png)

Output figures (prim, saved to --save)
    fig1_phase_cooling_CR<tag>[_JLW<jlw>].png
    fig2_species_CR<tag>[_JLW<jlw>].png
    fig3_jeans_CR<tag>[_JLW<jlw>].png
    fig4_thermal_balance_CR<tag>[_JLW<jlw>].png
    fig_summary_CR<tag>[_JLW<jlw>].png          (only with --fig-combo)

Output figures (metal, saved to --metal-save)
    fig1_phase_cooling_CR<cr>_Z<z>[_JLW<jlw>].png
    fig2_species_CR<cr>_Z<z>[_JLW<jlw>].png
    fig3_jeans_CR<cr>_Z<z>[_JLW<jlw>].png
    fig4_thermal_balance_CR<cr>_Z<z>[_JLW<jlw>].png
    fig_summary_CR<cr>_Z<z>[_JLW<jlw>].png      (only with --fig-combo)
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import h5py
except ImportError:
    sys.exit("ERROR: h5py is not installed.  Run: pip3 install h5py")

# ─── Physical constants ───────────────────────────────────────────────────────
M_SUN = 1.989e33   # g
YR    = 3.156e7    # s

# ─────────────────────────────────────────────────────────────────────────────
# I/O helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_case(h5path: str) -> dict:
    """Load one HDF5 file into a dict of numpy arrays + scalar attributes."""
    with h5py.File(h5path, "r") as f:
        d = {}
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


def get_species_idx(data: dict, name: str):
    try:
        return data["species"].index(name)
    except (ValueError, KeyError):
        return None


_SHOW = False

def _save(fig, save_dir: str, fname: str):
    if _SHOW:
        plt.show()
    else:
        path = os.path.join(save_dir, fname)
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"  saved: {path}")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Label helpers  (param label for titles — single math expression)
# ─────────────────────────────────────────────────────────────────────────────

def _zred_part(d: dict) -> str:
    """Return LaTeX redshift term if zred_tag is set, e.g. ',\\ z=20'."""
    tag = d.get("zred_tag", "")
    if not tag:
        return ""
    val = tag.replace("p", ".")   # "1p5" → "1.5", "20" → "20"
    return fr",\ z = {val}"


def _jlw_part(d: dict) -> str:
    """Return LaTeX J_LW term if J_LW21 attribute is present and non-zero."""
    jlw = float(d.get("J_LW21", 0.0))
    if jlw == 0.0:
        return ""
    return fr",\ J_{{LW}} = {jlw:.3g}\ J_{{21}}"


def _prim_label(d: dict) -> str:
    """Return LaTeX param string for prim titles, e.g. '$\\zeta_0=...$'."""
    zeta0 = d.get("zeta0_cgs", 0.0)
    f_ret = d.get("f_ret", 1.0)
    zs = r"\zeta_0 = 0" if zeta0 == 0.0 \
         else fr"\zeta_0 = {zeta0:.3e}\ \mathrm{{s}}^{{-1}}"
    fr_part = fr",\ f_{{\mathrm{{ret}}}} = {f_ret:g}" if f_ret != 1.0 else ""
    return f"${zs}{fr_part}{_jlw_part(d)}{_zred_part(d)}$"


def _metal_label(d: dict) -> str:
    """Return LaTeX param string for metal titles."""
    zeta0   = d.get("zeta0_cgs", 0.0)
    Z_metal = d.get("Z_metal",   0.0)
    f_ret   = d.get("f_ret", 1.0)
    zs = r"0" if zeta0   == 0.0 else fr"{zeta0:.3e}"
    Zs = r"0" if Z_metal == 0.0 else fr"{Z_metal:.3e}"
    base = fr"\zeta_0={zs}\ \mathrm{{s}}^{{-1}},\ Z/Z_{{\odot}}={Zs}"
    if f_ret != 1.0:
        base += fr",\ f_{{\mathrm{{ret}}}} = {f_ret:g}"
    return f"${base}{_jlw_part(d)}{_zred_part(d)}$"


# ─────────────────────────────────────────────────────────────────────────────
# Cooling component lists  (key, legend label, color, linestyle)
# ─────────────────────────────────────────────────────────────────────────────

_PRIM_COOL = [
    ("xLmbd_H2",  "H2 line",   "#1f77b4", "-"),
    ("xLmbd_HD",  "HD line",   "#2ca02c", "--"),
    ("xLmbd_Lya", "Ly-alpha",  "#9467bd", "-."),
    ("xLmbd_cnt", "continuum", "#7f7f7f", "-"),
    ("xLmbd_ch",  "chemical",  "#e377c2", (0, (3, 1, 1, 1))),
    ("xGam_CR",   "CR heat",   "#d62728", "-"),
    ("xLmbd_net", "net",       "k",       "-"),
]

_METAL_COOL = [
    ("xLmbd_H2",  "H2",             "#1f77b4", "-"),
    ("xLmbd_HD",  "HD",             "#2ca02c", "--"),
    ("xLmbd_CO",  "CO",             "#d62728", "-."),
    ("xLmbd_OH",  "OH",             "#9467bd", ":"),
    ("xLmbd_H2O", "H2O",            "#8c564b", (0, (5, 1))),
    ("xLmbd_CII", "C II",           "#e377c2", (0, (3, 1, 1, 1))),
    ("xLmbd_CI",  "C I",            "#17becf", (0, (3, 1))),
    ("xLmbd_OI",  "O I",            "#bcbd22", (0, (1, 1))),
    ("xLmbd_cnt", "cnt (total)",    "#7f7f7f", "-"),
    ("xLmbd_gr",  "cnt (grain)",    "#a0a000", "--"),
    ("xLmbd_gas", "cnt (gas)",      "#b0b0b0", ":"),
    ("xGam_CR",   "CR heat",        "#ff7f0e", "-"),
    ("xLmbd_net", "net",            "k",       "-"),
]

# Species color mapping for metal summary (aligned with _METAL_COOL where possible)
_METAL_SP_COLOR = {
    "H":   "tab:gray",
    "H2":  "#1f77b4",
    "e-":  "tab:purple",
    "C+":  "#e377c2",
    "C":   "#17becf",
    "CO":  "#d62728",
    "O":   "#bcbd22",
    "H2O": "#8c564b",
    "Gr":  "tab:olive",
}


# ─────────────────────────────────────────────────────────────────────────────
# Shared draw helpers
# ─────────────────────────────────────────────────────────────────────────────

def _draw_phase(ax, d: dict, is_metal: bool):
    """Draw nH–T phase diagram into *ax*.  Curve labels are descriptive only."""
    ax.plot(d["xnH"], d["T_K"], lw=1.6, label="gas $T$")
    if is_metal and "T_gr_K" in d:
        ax.plot(d["xnH"], d["T_gr_K"], lw=1.2, ls="--",
                color="tab:orange", label="grain $T$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylabel(r"$T$  [K]", fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)
    ax.set_xlim(1e-1, 1e23)
    ax.set_ylim(3, 3e4) if is_metal else ax.set_ylim(1e1, 2e4)


def _draw_cooling(ax, d: dict, is_metal: bool):
    """Draw cooling/heating rate curves into *ax*."""
    components = _METAL_COOL if is_metal else _PRIM_COOL
    nH = d["xnH"]
    for key, label, color, ls in components:
        if key not in d:
            continue
        vals = d[key]
        pos = np.where(vals > 0, vals, np.nan)
        neg = np.where(vals < 0, -vals, np.nan)
        lw = 1.5 if key in ("xLmbd_cnt", "xLmbd_net") else 1.2
        ax.plot(nH, pos, color=color, ls=ls, lw=lw, label=f"{label} (+)")
        if np.any(np.isfinite(neg)):
            ax.plot(nH, neg, color=color, ls="--", lw=0.7, alpha=0.5,
                    label=f"{label} (−)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylabel(r"$|\Lambda|$  [erg g$^{-1}$ s$^{-1}$]", fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim(1e-1, 1e23)
    ax.set_ylim(1e-8, 1e6) if is_metal else ax.set_ylim(1e-5, 1e5)
    h, l = ax.get_legend_handles_labels()
    seen = {}
    for hh, ll in zip(h, l):
        if ll and ll not in seen:
            seen[ll] = hh
    ncol = 3 if is_metal else 2
    ax.legend(seen.values(), seen.keys(), fontsize=7, ncol=ncol, loc="upper left")


# ─────────────────────────────────────────────────────────────────────────────
# fig1: phase + cooling (2-panel, sharex)
# ─────────────────────────────────────────────────────────────────────────────

def _phase_cooling(d: dict, label: str, tag_sfx: str, save_dir: str,
                   is_metal: bool):
    fig, (ax_phase, ax_cool) = plt.subplots(
        2, 1, figsize=(7, 10), sharex=True,
        gridspec_kw={"hspace": 0.05})

    _draw_phase(ax_phase, d, is_metal)
    _draw_cooling(ax_cool, d, is_metal)

    prefix = "Metal-grain" if is_metal else "Zero-metal"
    ax_phase.set_title(
        fr"{prefix} collapse: $n_\mathrm{{H}}$–$T$ and cooling/heating"
        f"\n[{label}]",
        fontsize=11)
    ax_cool.set_xlabel(r"$n_\mathrm{H}$  [cm$^{-3}$]", fontsize=11)

    fig.tight_layout()
    _save(fig, save_dir, f"fig1_phase_cooling_{tag_sfx}.png")


def prim_phase_cooling(d: dict, tag_sfx: str, save_dir: str):
    _phase_cooling(d, _prim_label(d), tag_sfx, save_dir, is_metal=False)


def metal_phase_cooling(d: dict, tag_sfx: str, save_dir: str):
    _phase_cooling(d, _metal_label(d), tag_sfx, save_dir, is_metal=True)


# ─────────────────────────────────────────────────────────────────────────────
# fig2: species
# ─────────────────────────────────────────────────────────────────────────────

def prim_species(d: dict, tag_sfx: str, save_dir: str):
    spec_list = [
        ("H",  r"$y_\mathrm{H}$"),
        ("H2", r"$y_\mathrm{H_2}$"),
        ("e-", r"$y_{e^-}$"),
        ("D",  r"$y_\mathrm{D}$"),
        ("HD", r"$y_\mathrm{HD}$"),
        ("H+", r"$y_\mathrm{H^+}$"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True)
    axes = axes.flatten()
    for ax, (sp, ylabel) in zip(axes, spec_list):
        idx = get_species_idx(d, sp)
        if idx is not None:
            ax.plot(d["xnH"], d["y"][:, idx], lw=1.4)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_ylabel(ylabel, fontsize=11)
        ax.grid(True, which="both", alpha=0.3)
        ax.set_xlim(1e-1, 1e23)
    for ax in axes[3:]:
        ax.set_xlabel(r"$n_\mathrm{H}$  [cm$^{-3}$]", fontsize=11)
    fig.suptitle(f"Zero-metal collapse: species\n[{_prim_label(d)}]",
                 fontsize=11)
    fig.tight_layout()
    _save(fig, save_dir, f"fig2_species_{tag_sfx}.png")


_METAL_SPECIES = [
    ("H",   r"$y_\mathrm{H}$"),
    ("H2",  r"$y_\mathrm{H_2}$"),
    ("e-",  r"$y_{e^-}$"),
    ("C+",  r"$y_\mathrm{C^+}$"),
    ("C",   r"$y_\mathrm{C}$"),
    ("CO",  r"$y_\mathrm{CO}$"),
    ("O",   r"$y_\mathrm{O}$"),
    ("H2O", r"$y_\mathrm{H_2O}$"),
    ("Gr",  r"$y_\mathrm{Gr^0}$"),
]


def metal_species(d: dict, tag_sfx: str, save_dir: str):
    nsp = len(_METAL_SPECIES)
    fig, axes = plt.subplots(3, 3, figsize=(12, 9), sharex=True)
    axes = axes.flatten()
    for ax, (sp, ylabel) in zip(axes, _METAL_SPECIES):
        idx = get_species_idx(d, sp)
        if idx is not None:
            ax.plot(d["xnH"], d["y"][:, idx], lw=1.3)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_ylabel(ylabel, fontsize=11)
        ax.grid(True, which="both", alpha=0.3)
        ax.set_xlim(1e-1, 1e23); ax.set_ylim(1e-20, 2.0)
    for ax in axes[6:]:
        ax.set_xlabel(r"$n_\mathrm{H}$  [cm$^{-3}$]", fontsize=10)
    for ax in axes[nsp:]:
        ax.set_visible(False)
    fig.suptitle(f"Metal-grain collapse: species\n[{_metal_label(d)}]",
                 fontsize=10)
    fig.tight_layout()
    _save(fig, save_dir, f"fig2_species_{tag_sfx}.png")


# ─────────────────────────────────────────────────────────────────────────────
# fig3: Jeans mass + free-fall time + γ_eff (3 panels)
# ─────────────────────────────────────────────────────────────────────────────

def _plot_jeans(d: dict, label: str, tag_sfx: str, save_dir: str):
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7, 11), sharex=True,
                                         gridspec_kw={"hspace": 0.05})
    nH  = d["xnH"]
    MJ  = d["xMJ"]  / M_SUN
    tff = d["t_ff"] / YR

    # γ_eff = d(ln P)/d(ln ρ) ≈ 1 + d(ln T)/d(ln ρ)
    lnrho     = np.log(d["rho"])
    lnT       = np.log(d["T_K"])
    gamma_eff = 1.0 + np.gradient(lnT, lnrho)
    gamma_eff = np.clip(gamma_eff, -1.0, 5.0)

    # ── panel 1: Jeans mass ───────────────────────────────────────────────────
    ax1.plot(nH, MJ, lw=1.6)
    ax1.axhline(1.0, color="gray", lw=0.8, ls="--", alpha=0.6,
                label=r"$1\,M_\odot$")
    ax1.set_yscale("log")
    ax1.set_ylabel(r"$M_\mathrm{J}$  [$M_\odot$]", fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, which="both", alpha=0.3)
    ax1.set_title(
        r"Jeans mass, free-fall time, and $\gamma_\mathrm{eff}$"
        f"\n[{label}]",
        fontsize=11)

    # ── panel 2: free-fall time ───────────────────────────────────────────────
    ax2.plot(nH, tff, lw=1.6)
    ax2.set_yscale("log")
    ax2.set_ylabel(r"$t_\mathrm{ff}$  [yr]", fontsize=12)
    ax2.grid(True, which="both", alpha=0.3)

    # ── panel 3: γ_eff ────────────────────────────────────────────────────────
    ax3.plot(nH, gamma_eff, lw=1.6)   # default color cycle (blue)
    ax3.axhline(1.0,     color="gray",   lw=1.4, ls="--", alpha=0.8,
                label=r"$\gamma_\mathrm{eff}=1$ (isothermal)")
    ax3.axhline(4.0/3.0, color="salmon", lw=1.4, ls=":",  alpha=0.8,
                label=r"$\gamma_\mathrm{eff}=4/3$")
    ax3.set_ylabel(r"$\gamma_\mathrm{eff}$", fontsize=12)
    ax3.set_xlabel(r"$n_\mathrm{H}$  [cm$^{-3}$]", fontsize=12)
    ax3.legend(fontsize=9)
    ax3.grid(True, which="both", alpha=0.3)
    ax3.set_ylim(-0.5, 4.5)

    for ax in (ax1, ax2, ax3):
        ax.set_xscale("log"); ax.set_xlim(1e-1, 1e23)

    fig.tight_layout()
    _save(fig, save_dir, f"fig3_jeans_{tag_sfx}.png")


def prim_jeans(d: dict, tag_sfx: str, save_dir: str):
    _plot_jeans(d, _prim_label(d), tag_sfx, save_dir)


def metal_jeans(d: dict, tag_sfx: str, save_dir: str):
    _plot_jeans(d, _metal_label(d), tag_sfx, save_dir)


# ─────────────────────────────────────────────────────────────────────────────
# fig4: thermal balance (3 panels)
#   panel 1: |Λ_net|  (+ Λ_gr / Λ_gas breakdown if present)
#   panel 2: Γ_CR / |Λ_net|
#   panel 3: t_cool/t_ff  and  t_chem/t_ff
# ─────────────────────────────────────────────────────────────────────────────

def _plot_thermal_balance(d: dict, label: str, tag_sfx: str, save_dir: str):
    fig, axes = plt.subplots(3, 1, figsize=(7, 10), sharex=True,
                             gridspec_kw={"hspace": 0.05})
    ax_net, ax_cr, ax_ratio = axes
    nH    = d["xnH"]
    Lnet  = d["xLmbd_net"]
    GCR   = d["xGam_CR"]
    tcool = d["t_cool"]
    tff   = d["t_ff"]

    # ── panel 1: |Λ_net| ─────────────────────────────────────────────────────
    pos_net = np.where(Lnet > 0, Lnet,  np.nan)
    neg_net = np.where(Lnet < 0, -Lnet, np.nan)
    ax_net.plot(nH, pos_net, lw=1.6)
    ax_net.plot(nH, neg_net, ls="--", lw=1.0, alpha=0.5)

    # grain / gas continuum breakdown (metal only — data-driven)
    for key, lbl, color, ls in [
        ("xLmbd_gr",  r"$\Lambda_\mathrm{gr}$",  "#a0a000", "--"),
        ("xLmbd_gas", r"$\Lambda_\mathrm{gas}$", "#888888", ":"),
    ]:
        if key in d:
            vals = d[key]
            ax_net.plot(nH, np.where(vals > 0, vals, np.nan),
                        color=color, ls=ls, lw=1.4, alpha=1.0, label=lbl)

    ax_net.set_yscale("log")
    ax_net.set_ylabel(r"$|\Lambda_\mathrm{net}|$  [erg g$^{-1}$ s$^{-1}$]",
                      fontsize=11)
    ax_net.set_title(f"Thermal balance\n[{label}]", fontsize=11)
    ax_net.legend(fontsize=9)
    ax_net.grid(True, which="both", alpha=0.3)

    # ── panel 2: CR fraction ──────────────────────────────────────────────────
    cr_frac = np.abs(GCR) / np.maximum(np.abs(Lnet), 1e-40)
    ax_cr.plot(nH, cr_frac, lw=1.4)
    ax_cr.set_yscale("log")
    ax_cr.set_ylabel(r"$\Gamma_\mathrm{CR}/|\Lambda_\mathrm{net}|$",
                     fontsize=11)
    ax_cr.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.5)
    ax_cr.grid(True, which="both", alpha=0.3)
    ax_cr.set_ylim(1e-6, 1e4)

    # ── panel 3: t_cool/t_ff ─────────────────────────────────────────────────
    tff_safe   = np.maximum(tff, 1.0)
    ratio_cool = np.minimum(tcool / tff_safe, 1e10)
    ax_ratio.plot(nH, ratio_cool, lw=1.4,
                  label=r"$t_\mathrm{cool}/t_\mathrm{ff}$")
    if "t_chem" in d:
        ratio_chem = np.minimum(d["t_chem"] / tff_safe, 1e10)
        ax_ratio.plot(nH, ratio_chem, lw=1.4, ls="--", color="tab:orange",
                      label=r"$t_\mathrm{chem}/t_\mathrm{ff}$")
    ax_ratio.set_yscale("log")
    ax_ratio.set_xlabel(r"$n_\mathrm{H}$  [cm$^{-3}$]", fontsize=12)
    ax_ratio.set_ylabel(r"$t/t_\mathrm{ff}$", fontsize=11)
    ax_ratio.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.5)
    ax_ratio.legend(fontsize=9)
    ax_ratio.grid(True, which="both", alpha=0.3)
    ax_ratio.set_ylim(1e-3, 1e7)

    for ax in axes:
        ax.set_xscale("log"); ax.set_xlim(1e-1, 1e23)
    fig.tight_layout()
    _save(fig, save_dir, f"fig4_thermal_balance_{tag_sfx}.png")


def prim_thermal_balance(d: dict, tag_sfx: str, save_dir: str):
    _plot_thermal_balance(d, _prim_label(d), tag_sfx, save_dir)


def metal_thermal_balance(d: dict, tag_sfx: str, save_dir: str):
    _plot_thermal_balance(d, _metal_label(d), tag_sfx, save_dir)


# ─────────────────────────────────────────────────────────────────────────────
# fig_summary: 2×2 combo figure  (--fig-combo)
#   [0,0] phase diagram      [0,1] key species
#   [1,0] cooling rates      [1,1] t_cool/t_ff + t_chem/t_ff (left) + M_J (right, twin axis)
# ─────────────────────────────────────────────────────────────────────────────

def _plot_summary(d: dict, tag_sfx: str, save_dir: str, is_metal: bool):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9),
                             gridspec_kw={"hspace": 0.30, "wspace": 0.32})
    ax_phase, ax_sp   = axes[0]
    ax_cool,  ax_diag = axes[1]

    # ── top-left: phase ───────────────────────────────────────────────────────
    _draw_phase(ax_phase, d, is_metal)
    prefix = "Metal-grain" if is_metal else "Zero-metal"
    ax_phase.set_title(fr"{prefix}: $n_\mathrm{{H}}$–$T$", fontsize=10)
    ax_phase.set_xlabel(r"$n_\mathrm{H}$  [cm$^{-3}$]", fontsize=10)

    # ── bottom-left: cooling ──────────────────────────────────────────────────
    _draw_cooling(ax_cool, d, is_metal)
    ax_cool.set_title("Cooling/heating", fontsize=10)
    ax_cool.set_xlabel(r"$n_\mathrm{H}$  [cm$^{-3}$]", fontsize=10)

    # ── top-right: species ────────────────────────────────────────────────────
    if is_metal:
        # All 9 species from _METAL_SPECIES with colors aligned to _METAL_COOL
        for sp, ylabel in _METAL_SPECIES:
            idx = get_species_idx(d, sp)
            if idx is not None:
                color = _METAL_SP_COLOR.get(sp)
                ax_sp.plot(d["xnH"], d["y"][:, idx],
                           lw=1.3, color=color, label=sp)
        ax_sp.legend(fontsize=8, ncol=3)
    else:
        for sp, lbl in [("H",  r"H"), ("H2", r"H$_2$"),
                        ("e-", r"$e^-$"), ("HD", r"HD")]:
            idx = get_species_idx(d, sp)
            if idx is not None:
                ax_sp.plot(d["xnH"], d["y"][:, idx], lw=1.4, label=lbl)
        ax_sp.legend(fontsize=9, ncol=2)

    ax_sp.set_xscale("log"); ax_sp.set_yscale("log")
    ax_sp.set_xlabel(r"$n_\mathrm{H}$  [cm$^{-3}$]", fontsize=10)
    ax_sp.set_ylabel("abundance", fontsize=10)
    ax_sp.set_title("Key species", fontsize=10)
    ax_sp.set_xlim(1e-1, 1e23)
    ax_sp.grid(True, which="both", alpha=0.3)

    # ── bottom-right: t_cool/t_ff + t_chem/t_ff (left) + M_J (right, twin) ───
    nH  = d["xnH"]
    tff = d["t_ff"]
    MJ  = d["xMJ"] / M_SUN
    tff_safe   = np.maximum(tff, 1.0)
    ratio_cool = np.minimum(d["t_cool"]  / tff_safe, 1e10)
    ratio_chem = np.minimum(d["t_chem"]  / tff_safe, 1e10)

    c_cool = "tab:blue"
    c_chem = "tab:green"
    c_mj   = "tab:orange"

    ax_diag.plot(nH, ratio_cool, lw=1.4, color=c_cool,
                 label=r"$t_\mathrm{cool}/t_\mathrm{ff}$")
    ax_diag.plot(nH, ratio_chem, lw=1.4, color=c_chem, ls="--",
                 label=r"$t_\mathrm{chem}/t_\mathrm{ff}$")
    ax_diag.axhline(1.0, color="gray", lw=0.8, ls=":", alpha=0.6)
    ax_diag.set_xscale("log"); ax_diag.set_yscale("log")
    ax_diag.set_xlabel(r"$n_\mathrm{H}$  [cm$^{-3}$]", fontsize=10)
    ax_diag.set_ylabel(r"$t / t_\mathrm{ff}$", fontsize=10)
    ax_diag.set_xlim(1e-1, 1e23); ax_diag.set_ylim(1e-3, 1e7)
    ax_diag.grid(True, which="both", alpha=0.3)
    ax_diag.set_title(
        r"Characteristic timescales / $t_\mathrm{ff}$  and  $M_\mathrm{J}$",
        fontsize=10)

    ax_mj = ax_diag.twinx()
    ax_mj.plot(nH, MJ, lw=1.4, color=c_mj, ls="-.",
               label=r"$M_\mathrm{J}$")
    ax_mj.axhline(1.0, color=c_mj, lw=0.7, ls=":", alpha=0.5)
    ax_mj.set_yscale("log")
    ax_mj.set_ylabel(r"$M_\mathrm{J}$  [$M_\odot$]",
                     fontsize=10, color=c_mj)
    ax_mj.tick_params(axis="y", labelcolor=c_mj)
    ax_mj.set_ylim(1e-3, 1e8)

    h1, l1 = ax_diag.get_legend_handles_labels()
    h2, l2 = ax_mj.get_legend_handles_labels()
    ax_diag.legend(h1 + h2, l1 + l2, fontsize=9, loc="upper right")

    label = _metal_label(d) if is_metal else _prim_label(d)
    fig.suptitle(f"{prefix} collapse summary\n[{label}]", fontsize=10)
    fig.tight_layout()
    _save(fig, save_dir, f"fig_summary_{tag_sfx}.png")


def prim_summary(d: dict, tag_sfx: str, save_dir: str):
    _plot_summary(d, tag_sfx, save_dir, is_metal=False)


def metal_summary(d: dict, tag_sfx: str, save_dir: str):
    _plot_summary(d, tag_sfx, save_dir, is_metal=True)


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    global _SHOW

    parser = argparse.ArgumentParser(
        description="Analyze and visualize one-zone collapse (single case, HDF5 input)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument("--h5dir",       default=None,
                        help="Directory with collapse_CR<tag>.h5 (prim)")
    parser.add_argument("--cr-tag",      default=None, dest="cr_tag",
                        help="CR tag for prim, e.g. '1p5e-17' or '0'")
    parser.add_argument("--fret-tag",    default="",   dest="fret_tag",
                        help="free-fall retardation tag for prim (omit or '' for f_ret=1)")
    parser.add_argument("--jlw-tag",     default="",   dest="jlw_tag",
                        help="Lyman-Werner tag for prim (omit or '' for J_LW21=0)")
    parser.add_argument("--zred-tag",    default="",   dest="zred_tag",
                        help="redshift tag for prim (omit or '' for z=0)")
    parser.add_argument("--save",        default="results/prim/fig",
                        help="Output directory for prim PNG files")
    parser.add_argument("--metal-h5dir", default=None, dest="metal_h5dir",
                        help="Directory with collapse_CR<cr>_Z<z>.h5 (metal)")
    parser.add_argument("--metal-cr",    default=None, dest="metal_cr",
                        help="CR tag for metal, e.g. '1p5e-17' or '0'")
    parser.add_argument("--metal-z",     default=None, dest="metal_z",
                        help="Z tag for metal, e.g. '1p2e-4' or '0'")
    parser.add_argument("--metal-fret",  default="",   dest="metal_fret",
                        help="free-fall retardation tag for metal (omit or '' for f_ret=1)")
    parser.add_argument("--metal-jlw",   default="",   dest="metal_jlw",
                        help="Lyman-Werner tag for metal (omit or '' for J_LW21=0)")
    parser.add_argument("--metal-zred",  default="",   dest="metal_zred",
                        help="redshift tag for metal (omit or '' for z=0)")
    parser.add_argument("--metal-save",  default="results/metal/fig", dest="metal_save",
                        help="Output directory for metal PNG files")
    parser.add_argument("--show", action="store_true",
                        help="Show figures interactively instead of saving")
    parser.add_argument("--fig-combo", action="store_true", dest="fig_combo",
                        help="Also output fig_summary_<tag>.png (2×2 overview)")
    args = parser.parse_args()

    _SHOW = args.show
    if _SHOW:
        matplotlib.use("TkAgg")

    # ── Prim figures ───────────────────────────────────────────────────────────
    if args.h5dir is not None:
        if not args.cr_tag:
            sys.exit("ERROR: --cr-tag is required with --h5dir")
        cr_tag   = args.cr_tag
        fret_tag = args.fret_tag
        jlw_tag  = args.jlw_tag
        zred_tag = args.zred_tag
        _stem = (f"collapse_CR{cr_tag}"
                 + (f"_fret{fret_tag}" if fret_tag else "")
                 + (f"_JLW{jlw_tag}"   if jlw_tag  else "")
                 + (f"_z{zred_tag}"    if zred_tag  else ""))
        h5path = os.path.join(args.h5dir, f"{_stem}.h5")
        if not os.path.isfile(h5path):
            sys.exit(f"ERROR: file not found: {h5path}")
        os.makedirs(args.save, exist_ok=True)
        print(f"Loading {h5path}")
        d = load_case(h5path)
        d["zred_tag"] = zred_tag
        zeta0 = d.get("zeta0_cgs", float("nan"))
        f_ret = d.get("f_ret", 1.0)
        jlw21 = float(d.get("J_LW21", 0.0))
        print(f"  zeta0={zeta0:.3e} s^-1  f_ret={f_ret}"
              + (f"  J_LW21={jlw21:.3g}" if jlw21 > 0.0 else "")
              + (f"  z={zred_tag}" if zred_tag else "")
              + f"  ({len(d['xnH'])} rows)")
        prim_sfx = (f"CR{cr_tag}"
                    + (f"_fret{fret_tag}" if fret_tag else "")
                    + (f"_JLW{jlw_tag}"   if jlw_tag  else "")
                    + (f"_z{zred_tag}"    if zred_tag  else ""))
        print("Generating prim figures (fig1–fig4)...")
        prim_phase_cooling   (d, prim_sfx, args.save)
        prim_species         (d, prim_sfx, args.save)
        prim_jeans           (d, prim_sfx, args.save)
        prim_thermal_balance (d, prim_sfx, args.save)
        if args.fig_combo:
            print("Generating prim summary figure...")
            prim_summary(d, prim_sfx, args.save)

    # ── Metal figures ──────────────────────────────────────────────────────────
    if args.metal_h5dir is not None:
        if not args.metal_cr or not args.metal_z:
            sys.exit("ERROR: --metal-cr and --metal-z are required with --metal-h5dir")
        cr_tag   = args.metal_cr
        z_tag    = args.metal_z
        fret_tag = args.metal_fret
        jlw_tag  = args.metal_jlw
        zred_tag = args.metal_zred
        _stem = (f"collapse_CR{cr_tag}_Z{z_tag}"
                 + (f"_fret{fret_tag}" if fret_tag else "")
                 + (f"_JLW{jlw_tag}"   if jlw_tag  else "")
                 + (f"_z{zred_tag}"    if zred_tag  else ""))
        h5path = os.path.join(args.metal_h5dir, f"{_stem}.h5")
        if not os.path.isfile(h5path):
            sys.exit(f"ERROR: file not found: {h5path}")
        os.makedirs(args.metal_save, exist_ok=True)
        print(f"Loading {h5path}")
        d = load_case(h5path)
        d["zred_tag"] = zred_tag
        zeta0   = d.get("zeta0_cgs", float("nan"))
        Z_metal = d.get("Z_metal",   float("nan"))
        f_ret   = d.get("f_ret", 1.0)
        jlw21   = float(d.get("J_LW21", 0.0))
        print(f"  zeta0={zeta0:.3e} s^-1  Z_metal={Z_metal:.3e}  f_ret={f_ret}"
              + (f"  J_LW21={jlw21:.3g}" if jlw21 > 0.0 else "")
              + (f"  z={zred_tag}" if zred_tag else "")
              + f"  ({len(d['xnH'])} rows)")
        metal_sfx = (f"CR{cr_tag}_Z{z_tag}"
                     + (f"_fret{fret_tag}" if fret_tag else "")
                     + (f"_JLW{jlw_tag}"   if jlw_tag  else "")
                     + (f"_z{zred_tag}"    if zred_tag  else ""))
        print("Generating metal figures (fig1–fig4)...")
        metal_phase_cooling   (d, metal_sfx, args.metal_save)
        metal_species         (d, metal_sfx, args.metal_save)
        metal_jeans           (d, metal_sfx, args.metal_save)
        metal_thermal_balance (d, metal_sfx, args.metal_save)
        if args.fig_combo:
            print("Generating metal summary figure...")
            metal_summary(d, metal_sfx, args.metal_save)

    print("Done.")


if __name__ == "__main__":
    main()
