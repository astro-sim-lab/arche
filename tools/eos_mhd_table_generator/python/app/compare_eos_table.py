#!/usr/bin/env python3
"""
compare_eos_table.py - Quick qualitative comparison for 10-column EOS tables.

Compares two .d files on overlapping (nt, bt) points and reports:
- per-column Pearson correlation
- median absolute difference
- 5/50/95 percentile of signed difference
for tt, mut, pt, etaO, etaA, etaH.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from lib.compare_core import compare_tables


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare two EOS-MHD .d tables.")
    p.add_argument("--base", required=True, help="Reference .d file")
    p.add_argument("--target", required=True, help="Generated .d file")
    p.add_argument("--corr-threshold", type=float,
                   help="Optional pass/fail threshold: require corr >= value for all tracked columns.")
    p.add_argument("--med-abs-threshold", type=float,
                   help="Optional pass/fail threshold: require median absolute diff <= value for all tracked columns.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    return compare_tables(
        Path(args.base),
        Path(args.target),
        corr_threshold=args.corr_threshold,
        med_abs_threshold=args.med_abs_threshold,
    )


if __name__ == "__main__":
    raise SystemExit(main())
