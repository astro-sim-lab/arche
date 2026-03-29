"""
Case-name helpers for EOS-MHD matrix operations.
"""

from __future__ import annotations

from pathlib import Path


def build_case_name(cr: str, sh: str, z: str) -> str:
    return f"CR{cr}_sh{sh}_Z{z}.d"


def z_aliases(z: str) -> list[str]:
    z0 = z.strip().lower()
    aliases = [z.strip()]
    if z0 in ("0", "0.0", "0e0", "0.0e0", "1e-inf", "1e-infinity"):
        for c in ("0", "0.0", "1e-inf"):
            if c not in aliases:
                aliases.append(c)
    return aliases


def case_name_candidates(cr: str, sh: str, z: str) -> list[str]:
    names: list[str] = []
    for zc in z_aliases(z):
        n = build_case_name(cr, sh, zc)
        if n not in names:
            names.append(n)
    return names


def resolve_case_file(directory: Path, cr: str, sh: str, z: str) -> tuple[Path | None, str]:
    for name in case_name_candidates(cr, sh, z):
        p = directory / name
        if p.is_file():
            return p, name
    return None, build_case_name(cr, sh, z)

