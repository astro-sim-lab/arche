#!/usr/bin/env python3
"""
validate_eos_table.py - Minimal format/axis validator for EOS-MHD .d files.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from lib.validate_core import validate_table


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Validate EOS-MHD .d table format.")
    p.add_argument("--input", required=True, help="Path to .d file")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    try:
        validate_table(Path(args.input))
    except Exception as e:  # noqa: BLE001
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
