"""
Library modules for tools/eos_mhd_table_generator.

This package hosts reusable implementation modules.
CLI scripts in tools/eos_mhd_table_generator should stay thin wrappers.
"""

from .build_core import build_table
from .compare_core import compare_tables
from .config_kv import coalesce_cli_cfg, parse_bool, read_kv_config
from .mapping import map89_to_64
from .resistivity_legacy import LegacyTables, eta_and_as_single, eta_grid, load_legacy_tables
from .table_io import EOS_COLS, extract_series_by_bt, monotonic_non_decreasing, read_table, split_blocks_by_ii, to_keyed_nt_bt
from .validate_core import validate_table

__all__ = [
    "build_table",
    "compare_tables",
    "coalesce_cli_cfg",
    "parse_bool",
    "read_kv_config",
    "validate_table",
    "EOS_COLS",
    "extract_series_by_bt",
    "monotonic_non_decreasing",
    "read_table",
    "split_blocks_by_ii",
    "to_keyed_nt_bt",
    "LegacyTables",
    "eta_and_as_single",
    "eta_grid",
    "load_legacy_tables",
    "map89_to_64",
]
