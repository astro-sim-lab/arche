"""
Core plotting logic for EOS-table comparisons.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .table_io import EOS_COLS, extract_series_by_bt


PLOT_COLS = {
    "tt": EOS_COLS["tt"],
    "mut": EOS_COLS["mut"],
    "pt": EOS_COLS["pt"],
    "etaO": EOS_COLS["etaO"],
    "etaA": EOS_COLS["etaA"],
    "etaH": EOS_COLS["etaH"],
}

YLABEL = {
    "tt": r"$\log_{10} T$",
    "mut": r"$\log_{10} \mu$",
    "pt": r"$\log_{10} P$",
    "etaO": r"$\log_{10}\eta_{\mathrm{O}}$",
    "etaA": r"$\log_{10}\eta_{\mathrm{A}}$",
    "etaH": r"$\log_{10}\eta_{\mathrm{H}}$",
}


def common_interp(xb: np.ndarray, yb: np.ndarray, xt: np.ndarray, yt: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    xmin = max(float(np.min(xb)), float(np.min(xt)))
    xmax = min(float(np.max(xb)), float(np.max(xt)))
    if xmax <= xmin:
        raise ValueError("no overlap in nt axis")
    x = xb[(xb >= xmin) & (xb <= xmax)]
    if x.size < 2:
        x = xt[(xt >= xmin) & (xt <= xmax)]
    yb_i = np.interp(x, xb, yb)
    yt_i = np.interp(x, xt, yt)
    return x, yb_i, yt_i


def rel_err_from_log(yb_log: np.ndarray, yt_log: np.ndarray) -> np.ndarray:
    b = np.power(10.0, yb_log)
    t = np.power(10.0, yt_log)
    den = np.maximum(np.abs(b), 1.0e-300)
    return (t - b) / den


def plot_quantity(base: np.ndarray, target: np.ndarray, qname: str, bt_list: List[float], out_png: Path) -> None:
    if qname not in PLOT_COLS:
        raise ValueError(f"unknown quantity: {qname}")

    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(8.2, 6.2), sharex=True, gridspec_kw={"height_ratios": [2.0, 1.0]}
    )
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(bt_list), 1)))

    for i, bt in enumerate(bt_list):
        xb, yb = extract_series_by_bt(base, bt=bt, col_name=qname)
        xt, yt = extract_series_by_bt(target, bt=bt, col_name=qname)
        x, yb_i, yt_i = common_interp(xb, yb, xt, yt)
        rel = rel_err_from_log(yb_i, yt_i)

        c = colors[i]
        ax0.plot(x, yb_i, color=c, lw=1.6, ls="-", label=f"base bt={bt:g}")
        ax0.plot(x, yt_i, color=c, lw=1.2, ls="--", label=f"target bt={bt:g}")
        ax1.plot(x, rel, color=c, lw=1.2, label=f"bt={bt:g}")

    ax0.set_ylabel(YLABEL[qname])
    ax0.grid(True, alpha=0.25)
    ax0.legend(ncol=2, fontsize=8)

    ax1.axhline(0.0, color="k", lw=0.8, alpha=0.7)
    ax1.set_xlabel(r"$\log_{10} n_{\mathrm{H}}$")
    ax1.set_ylabel("rel. err")
    ax1.grid(True, alpha=0.25)

    fig.suptitle(f"{qname}: base vs target and relative error")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

