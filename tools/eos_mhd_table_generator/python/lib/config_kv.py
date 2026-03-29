"""
Shared KEY=VALUE config helpers for eos_table tools.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def read_kv_config(path: Path) -> dict[str, str]:
    cfg: dict[str, str] = {}
    if not path.is_file():
        raise FileNotFoundError(f"config file not found: {path}")
    for ln, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        if "=" not in s:
            raise ValueError(f"{path}:{ln}: expected KEY=VALUE")
        k, v = s.split("=", 1)
        cfg[k.strip()] = v.strip()
    return cfg


def coalesce_cli_cfg(cli_val: Any, cfg: dict[str, str], cfg_key: str, fallback: Any = None) -> Any:
    if cli_val is not None:
        return cli_val
    if cfg_key in cfg:
        return cfg[cfg_key]
    return fallback


def parse_bool(v: Any) -> bool:
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in ("1", "true", "yes", "on"):
        return True
    if s in ("0", "false", "no", "off"):
        return False
    raise ValueError(f"invalid bool value: {v}")

