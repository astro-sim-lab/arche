#!/usr/bin/env python3
"""
plot_compare_eos_table.py

Plot base/target EOS-table profiles (x: nt=log10 nH) and relative error
in stacked panels for selected quantities.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from lib.plot_core import PLOT_COLS, plot_quantity
from lib.table_io import read_table


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot EOS table comparison (stacked value + relative error).")
    p.add_argument("--base", required=True, help="Reference .d")
    p.add_argument("--target", required=True, help="Generated .d")
    p.add_argument("--out-dir", required=True, help="Output directory for PNGs")
    p.add_argument("--quantities", default="tt,mut,pt,etaO,etaA,etaH",
                   help="Comma-separated quantities from tt,mut,pt,etaO,etaA,etaH")
    p.add_argument("--bt-list", default="-20,-10,0,2",
                   help="Comma-separated log10(B) slices to plot")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    base = read_table(Path(args.base), strict_cols=False)
    target = read_table(Path(args.target), strict_cols=False)
    out_dir = Path(args.out_dir)

    qnames = [q.strip() for q in args.quantities.split(",") if q.strip()]
    for q in qnames:
        if q not in PLOT_COLS:
            raise ValueError(f"unknown quantity: {q}")

    bt_list = [float(x.strip()) for x in args.bt_list.split(",") if x.strip()]
    for q in qnames:
        plot_quantity(base, target, q, bt_list, out_dir / f"compare_{q}.png")

    print(f"OK: wrote {len(qnames)} plot(s) to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
