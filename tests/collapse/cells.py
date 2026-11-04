# Copyright (C) 2026 Shingo Hirano and Sho Higashi
# Licensed under the MIT found in the
# https://github.com/astro-sim-lab/arche/blob/main/LICENSE
"""Collapse test cells: the grid G (160 cells) and its subset S (16 cells).

Every cell is a one-zone collapse that must reach n_H = 1e23 cm^-3 normally.
The cells cover the code paths of the four networks with initial parameters
that can occur physically:

  G1  primordial  {pf, pm} x zeta {0, 1e-17, 1e-15} x f_ret {1, 3, 10}, z = 20
  G2  metal       {mf, mm} x Z {1, ..., 1e-6} x zeta {0, 1e-17, 1e-15} x f_ret {1, 3, 10}, z = 0
  G3  metal, Z=0  {mf, mm} x zeta 0, z = 20 (compared with pf, zeta 0, f_ret 1, z = 20)
  G4  prim. LW    {pf, pm} x J_LW21 {1e3, 1e4, 1e5} x f_ret {1, 3}, zeta 0, z = 20
  G5  metal LW    {mf, mm} x Z 1e-3 x J_LW21 1e5, zeta 0

Primordial gas and Z = 0 use z = 20 (T_CMB = 57.2 K): primordial collapse
happens only at high redshift.  S is the part of G run on the additional builds
(cross-build comparison): the four networks at zeta {0, 1e-17}, f_ret 1, with
Z {1e-6, 1e-3, 1} for the metal networks.

Units: zeta [s^-1], Z [Z_sun], J_LW21 [1e-21 erg s^-1 cm^-2 Hz^-1 sr^-1].
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

Network = Literal["pf", "pm", "mf", "mm"]

# Executable of each network, relative to <build-dir>/src/apps.
BINARIES: dict[Network, str] = {
    "pf": "collapse_primordial/prim_collapse",
    "pm": "collapse_primordial/arche_collapse_prim_minimal",
    "mf": "collapse_metal_grain/metal_collapse",
    "mm": "collapse_metal_grain/arche_collapse_metal_minimal",
}

ZETAS: tuple[float, ...] = (0.0, 1e-17, 1e-15)
F_RETS: tuple[float, ...] = (1.0, 3.0, 10.0)
METALLICITIES: tuple[float, ...] = (1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
JLW21S: tuple[float, ...] = (1e3, 1e4, 1e5)
Z_HIGH = 20.0  # redshift of primordial and Z = 0 cells

# Cells that are known to fail, mapped to the issue that tracks them.  The
# checker reports them as XFAIL (not a failure) and as XPASS once they pass,
# so that the entry can be removed.
KNOWN_FAILURES: dict[str, str] = {}


def _tag(x: float) -> str:
    """Value as passed to the application (any form strtod accepts)."""
    return f"{x:g}"


def _name_tag(x: float) -> str:
    """Short value for cell names: powers of ten as 1e<n> (1e-3, 1e5), else %g."""
    if x > 0.0:
        e = round(math.log10(x))
        if abs(x - 10.0**e) <= 1e-12 * x and e not in (0, 1):
            return f"1e{e}"
    return f"{x:g}"


@dataclass(frozen=True)
class Cell:
    group: str
    network: Network
    zeta: float
    f_ret: float
    z_metal: float | None  # None for the primordial networks
    jlw21: float
    redshift: float
    subset_s: bool

    @property
    def is_metal(self) -> bool:
        return self.network in ("mf", "mm")

    @property
    def name(self) -> str:
        parts = [self.group, self.network]
        if self.z_metal is not None:
            parts.append(f"Z{_name_tag(self.z_metal)}")
        parts += [f"zeta{_name_tag(self.zeta)}", f"fret{_name_tag(self.f_ret)}"]
        if self.jlw21 > 0.0:
            parts.append(f"J{_name_tag(self.jlw21)}")
        parts.append(f"z{_name_tag(self.redshift)}")
        return "_".join(parts)

    def env(self, outdir: str, nh_stop: float) -> dict[str, str]:
        """Environment variables that run this cell (see docs/parameters.md)."""
        p = "METAL" if self.is_metal else "PRIM"
        env = {
            f"{p}_ZETA0": _tag(self.zeta),
            f"{p}_FF_RET": _tag(self.f_ret),
            f"{p}_JLW21": _tag(self.jlw21),
            f"{p}_REDSHIFT": _tag(self.redshift),
            f"{p}_XNH_STOP": _tag(nh_stop),
            f"{p}_OUTDIR": outdir,
        }
        if self.is_metal:
            assert self.z_metal is not None
            env["METAL_Z_METAL"] = _tag(self.z_metal)
        return env

    def params(self) -> dict[str, object]:
        return {
            "group": self.group,
            "network": self.network,
            "zeta": self.zeta,
            "f_ret": self.f_ret,
            "z_metal": self.z_metal,
            "jlw21": self.jlw21,
            "redshift": self.redshift,
            "subset_s": self.subset_s,
        }


def _prim(group: str, net: Network, zeta: float, f_ret: float, jlw21: float, s: bool) -> Cell:
    return Cell(group, net, zeta, f_ret, None, jlw21, Z_HIGH, s)


def grid() -> list[Cell]:
    """The grid G, in a fixed order."""
    cells: list[Cell] = []
    prim: tuple[Network, Network] = ("pf", "pm")
    metal: tuple[Network, Network] = ("mf", "mm")
    for net in prim:  # G1
        for zeta in ZETAS:
            for f_ret in F_RETS:
                s = zeta in (0.0, 1e-17) and f_ret == 1.0
                cells.append(_prim("G1", net, zeta, f_ret, 0.0, s))
    for net in metal:  # G2
        for z_metal in METALLICITIES:
            for zeta in ZETAS:
                for f_ret in F_RETS:
                    s = z_metal in (1e-6, 1e-3, 1.0) and zeta in (0.0, 1e-17) and f_ret == 1.0
                    cells.append(Cell("G2", net, zeta, f_ret, z_metal, 0.0, 0.0, s))
    for net in metal:  # G3
        cells.append(Cell("G3", net, 0.0, 1.0, 0.0, 0.0, Z_HIGH, False))
    for net in prim:  # G4
        for jlw21 in JLW21S:
            for f_ret in (1.0, 3.0):
                cells.append(_prim("G4", net, 0.0, f_ret, jlw21, False))
    for net in metal:  # G5
        cells.append(Cell("G5", net, 0.0, 1.0, 1e-3, 1e5, 0.0, False))
    return cells


def select(subset: str) -> list[Cell]:
    """subset 'G' (all 160 cells) or 'S' (16 cells)."""
    cells = grid()
    if subset == "G":
        return cells
    if subset == "S":
        return [c for c in cells if c.subset_s]
    raise ValueError(f"unknown subset {subset!r} (expected 'G' or 'S')")


# The G3 cells are compared with this primordial cell (same physics at Z = 0).
G3_REFERENCE = _prim("G1", "pf", 0.0, 1.0, 0.0, True).name

if __name__ == "__main__":
    g = grid()
    print(f"G: {len(g)} cells, S: {len(select('S'))} cells")
    for c in g:
        print(("S " if c.subset_s else "  ") + c.name)
