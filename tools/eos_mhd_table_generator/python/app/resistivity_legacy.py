#!/usr/bin/env python3
"""
Compatibility shim for legacy imports.

Implementation moved to tools/eos_mhd_table_generator/python/lib/resistivity_legacy.py.
"""

import sys
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from lib.resistivity_legacy import *  # noqa: F401,F403
