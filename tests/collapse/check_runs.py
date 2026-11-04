#!/usr/bin/env python3
# Copyright (C) 2026 Shingo Hirano and Sho Higashi
# Licensed under the MIT found in the
# https://github.com/astro-sim-lab/arche/blob/main/LICENSE
"""Check collapse test runs written by run_cells.py.

Per cell (FAIL unless stated):
  J1  normal end: finished before the gate, HDF5 exit_code 0, stdout EXIT code 0,
      all n_H and T finite.  (The applications return 0 to the shell even on an
      abnormal end, so the HDF5 attribute and stdout decide.)
  J2  conservation: projected steps = all steps, max |charge_imbal| <= 1e-14.
  J4  wall time: above the gate -> FAIL (J1, killed); above the report
      threshold -> WARN.
Across cells of one build:
  G3  metal full (mf) at Z = 0 vs primordial full (zeta 0, f_ret 1, z 20):
      max |dlog10 T| <= 1e-4 dex on the common log n_H range.  mf contains the
      primordial full network, so the two must agree at Z = 0.  The metal
      minimal network (mm) is a different reduction (neither pf nor pm is a
      subset of it); its difference is reported, not judged.
Across builds (--compare):
  J3  each cell present in the reference and in another build: |dlog10 T| <= 1e-3 dex
      at equal n_H on the common log n_H range.  A point above 1e-3 dex is accepted
      only if shifting the other track by at most 0.01 dex in log n_H brings it within
      3e-3 dex (checked both ways), and the accepted points span at most 0.1 dex of
      log n_H in total.  This admits a steep feature that two builds place slightly
      apart in n_H: a rapid cooling (T drops 0.095 dex) that g++ 11.4 and g++ 13.3
      start 0.0064 dex apart, after a rise that the later one continues 2.5e-3 dex
      higher.  An offset over more than 0.1 dex of log n_H still fails.

Cells in cells.KNOWN_FAILURES are reported as XFAIL instead of FAIL, and as
XPASS when they pass.  Exit status 1 if any FAIL.

  check_runs.py DIR [DIR ...]                       # results of one build (shards merged)
  check_runs.py --compare REF_DIR OTHER_DIR [...]   # J3; each argument may be a comma-separated list of shard dirs
  check_runs.py --self-test                         # J3 on synthetic tracks
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cells import G3_REFERENCE, KNOWN_FAILURES  # noqa: E402

CHARGE_IMBAL_MAX = 1e-14
G3_TOL_DEX = 1e-4
CROSS_TOL_DEX = 1e-3
CROSS_SHIFT_DEX = 1e-2  # largest shift in log n_H that may explain a difference
CROSS_SHIFT_TOL_DEX = 3e-3  # |dlog T| allowed after that shift
CROSS_SHIFT_SPAN_DEX = 1e-1  # largest total log n_H span explained only by a shift
CROSS_GRID_STEP_DEX = 5e-4
N_GRID = 1000


@dataclass
class Report:
    fails: list[str] = field(default_factory=list)
    warns: list[str] = field(default_factory=list)
    xfails: list[str] = field(default_factory=list)
    xpasses: list[str] = field(default_factory=list)
    lines: list[str] = field(default_factory=list)


def load(dirs: list[Path]) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]], list[str]]:
    """Merge results.json of several shard directories; also return the missing ones."""
    metas: list[dict[str, Any]] = []
    cells: dict[str, dict[str, Any]] = {}
    missing: list[str] = []
    for d in dirs:
        if not (d / "results.json").is_file():
            missing.append(str(d))
            continue
        data = json.loads((d / "results.json").read_text())
        metas.append(data["meta"])
        for rec in data["cells"]:
            rec["_dir"] = str(d)
            cells[rec["name"]] = rec
    return metas, cells, missing


def trace(rec: dict[str, Any]) -> tuple[NDArray[np.float64], NDArray[np.float64]] | None:
    p = Path(rec["_dir"]) / rec["name"] / "trace.npz"
    if not p.is_file():
        return None
    with np.load(p) as z:
        x = np.asarray(z["log_nH"], dtype=np.float64)
        y = np.asarray(z["log_T"], dtype=np.float64)
    ok = np.isfinite(x) & np.isfinite(y)
    return x[ok], y[ok]


def max_dlogt(a: dict[str, Any], b: dict[str, Any]) -> float | None:
    """max |log T_a - log T_b| on a common log n_H grid (both tracks interpolated)."""
    ta, tb = trace(a), trace(b)
    if ta is None or tb is None or ta[0].size < 2 or tb[0].size < 2:
        return None
    lo = max(ta[0][0], tb[0][0])
    hi = min(ta[0][-1], tb[0][-1])
    if not hi > lo:
        return None
    x = np.linspace(lo, hi, N_GRID)
    # n_H increases monotonically along a collapse; np.maximum.accumulate guards
    # against equal neighbours from the output rounding.
    ya = np.interp(x, np.maximum.accumulate(ta[0]), ta[1])
    yb = np.interp(x, np.maximum.accumulate(tb[0]), tb[1])
    return float(np.max(np.abs(ya - yb)))


@dataclass(frozen=True)
class CrossDiff:
    max_equal_nh: float  # max |dlog T| at equal n_H
    max_unexplained: float  # max |dlog T| at equal n_H over points a shift does not explain
    shifted_span: float  # log n_H span of points explained only by a shift

    @property
    def ok(self) -> bool:
        return self.max_unexplained <= CROSS_TOL_DEX and self.shifted_span <= CROSS_SHIFT_SPAN_DEX


def _window_min_max(y: NDArray[np.float64], w: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """min and max of y over [i - w, i + w] for every i (edge values repeated)."""
    win = np.lib.stride_tricks.sliding_window_view(np.pad(y, w, mode="edge"), 2 * w + 1)
    return win.min(axis=1), win.max(axis=1)


def _shift_distance(ya: NDArray[np.float64], yb: NDArray[np.float64], w: int) -> NDArray[np.float64]:
    """Distance from ya[i] to the range of yb over i +- w grid points."""
    lo, hi = _window_min_max(yb, w)
    return np.where(ya < lo, lo - ya, np.where(ya > hi, ya - hi, 0.0))


def cross_diff_arrays(
    xa: NDArray[np.float64], ya: NDArray[np.float64], xb: NDArray[np.float64], yb: NDArray[np.float64]
) -> CrossDiff | None:
    """J3 on two tracks (log n_H increasing, log T)."""
    if xa.size < 2 or xb.size < 2:
        return None
    lo, hi = max(xa[0], xb[0]), min(xa[-1], xb[-1])
    if not hi > lo:
        return None
    n = int(np.ceil((hi - lo) / CROSS_GRID_STEP_DEX)) + 1
    x = np.linspace(lo, hi, n)
    step = (hi - lo) / (n - 1)
    ga = np.interp(x, np.maximum.accumulate(xa), ya)
    gb = np.interp(x, np.maximum.accumulate(xb), yb)
    d = np.abs(ga - gb)
    w = max(1, int(round(CROSS_SHIFT_DEX / step)))
    shift = np.maximum(_shift_distance(ga, gb, w), _shift_distance(gb, ga, w))
    over = d > CROSS_TOL_DEX
    explained = over & (shift <= CROSS_SHIFT_TOL_DEX)
    rest = d[~explained]
    return CrossDiff(
        max_equal_nh=float(d.max()),
        max_unexplained=float(rest.max()) if rest.size else 0.0,
        shifted_span=float(np.count_nonzero(explained) * step),
    )


def cross_diff(a: dict[str, Any], b: dict[str, Any]) -> CrossDiff | None:
    ta, tb = trace(a), trace(b)
    if ta is None or tb is None:
        return None
    return cross_diff_arrays(ta[0], ta[1], tb[0], tb[1])


def self_test() -> int:
    """J3 on synthetic tracks: a smooth track with a steep 0.025 dex drop."""
    x = np.linspace(-1.0, 23.0, 6000)

    def track(x0: float) -> NDArray[np.float64]:
        return 2.0 + 0.1 * x - 0.0125 * (1.0 + np.tanh((x - x0) / 4e-5))

    def rise_drop(x0: float) -> NDArray[np.float64]:
        # rises with slope 0.23 over 0.02 dex before the drop, as in the observed track
        rise = 0.23 * np.clip(x - (x0 - 0.02), 0.0, None) * (x < x0)
        return track(x0) + rise + 0.23 * 0.02 * (x >= x0)

    y = track(11.37)
    cases = [
        ("identical", x, y, True),
        ("drop shifted by 0.0066 dex in n_H", x, track(11.3766), True),
        ("drop shifted by 0.02 dex in n_H", x, track(11.39), False),
        ("offset 5e-4 dex everywhere", x, y + 5e-4, True),
        ("offset 2e-3 dex everywhere", x, y + 2e-3, False),
        ("bump 2e-3 dex over 0.2 dex of n_H", x, y + 2e-3 * (np.abs(x - 5.0) < 0.1), False),
        ("bump 5e-3 dex over 0.02 dex of n_H", x, y + 5e-3 * (np.abs(x - 5.0) < 0.01), False),
        ("rise 0.23 then drop, onset 0.0064 dex later", x, rise_drop(11.3764), True),
        ("whole track shifted by 0.03 dex in n_H", x + 0.03, y, False),
    ]
    bad = 0
    for label, xb, yb, expect in cases:
        ya = rise_drop(11.37) if label.startswith("rise") else y
        r = cross_diff_arrays(x, ya, xb, yb)
        assert r is not None
        status = "ok" if r.ok == expect else "WRONG"
        bad += status != "ok"
        print(
            f"{status:5s} {label:42s} pass={r.ok!s:5s} (expected {expect!s:5s}) "
            f"equal-n_H {r.max_equal_nh:.2e}, unexplained {r.max_unexplained:.2e}, shifted span {r.shifted_span:.4f}"
        )
    return 1 if bad else 0


def cell_problems(rec: dict[str, Any], warn_sec: float) -> tuple[list[str], list[str]]:
    fails: list[str] = []
    warns: list[str] = []
    if rec["status"] == "timeout":
        return [f"J4 killed at the gate ({rec['timeout_sec']:.0f} s)"], warns
    if rec.get("h5_error"):
        fails.append(f"J1 HDF5: {rec['h5_error']}")
    if rec.get("exit_code_h5") != 0:
        fails.append(f"J1 exit_code {rec.get('exit_code_h5')} ({rec.get('exit_message', '')})")
    if rec.get("exit_code_stdout") != 0:
        fails.append(f"J1 stdout EXIT code {rec.get('exit_code_stdout')}")
    if rec.get("finite") is False:
        fails.append("J1 non-finite n_H or T")
    proj = rec.get("projected")
    if proj is None or proj[0] != proj[1]:
        fails.append(f"J2 projected {proj}")
    ci = rec.get("max_abs_charge_imbal")
    if ci is None or not ci <= CHARGE_IMBAL_MAX:
        fails.append(f"J2 max|charge_imbal| {ci}")
    if rec["wall_sec"] > warn_sec:
        warns.append(f"J4 slow: {rec['wall_sec']:.0f} s (report threshold {warn_sec:.0f} s)")
    return fails, warns


def classify(rep: Report, name: str, fails: list[str]) -> None:
    if name in KNOWN_FAILURES:
        if fails:
            rep.xfails.append(f"{name}: {'; '.join(fails)} [{KNOWN_FAILURES[name]}]")
        else:
            rep.xpasses.append(f"{name}: passes; remove it from KNOWN_FAILURES [{KNOWN_FAILURES[name]}]")
    elif fails:
        rep.fails.append(f"{name}: {'; '.join(fails)}")


def check_build(dirs: list[Path]) -> Report:
    rep = Report()
    metas, cells, missing = load(dirs)
    rep.fails += [f"{d}: results.json missing (the run did not finish)" for d in missing]
    if not metas:
        return rep
    warn_sec = min(float(m["warn_sec"]) for m in metas)
    m0 = metas[0]
    rep.lines.append(
        f"build: {m0['cmake'].get('CMAKE_CXX_COMPILER', '?')} ARCHE_MARCH={m0['cmake'].get('ARCHE_MARCH', '?')!r}, "
        f"commit {m0['commit']}, CPU {m0['cpu']}, scale {m0['scale']:.2f}, {len(cells)} cells"
    )
    for name, rec in cells.items():
        fails, warns = cell_problems(rec, warn_sec)
        classify(rep, name, fails)
        rep.warns += [f"{name}: {w}" for w in warns]

    g3 = [r for r in cells.values() if r["group"] == "G3"]
    if g3:
        ref = cells.get(G3_REFERENCE)
        for r in g3:
            if ref is None:
                rep.warns.append(f"{r['name']}: G3 not compared ({G3_REFERENCE} not in these results)")
                continue
            d = max_dlogt(r, ref)
            judged = r["network"] == "mf"
            msg = f"G3 vs {G3_REFERENCE}: max|dlogT| = {d}" + ("" if judged else " (reported only)")
            rep.lines.append(f"{r['name']}: {msg}")
            if judged:
                classify(rep, r["name"] + " (G3)", [] if d is not None and d <= G3_TOL_DEX else [msg])

    walls = sorted(((float(r["wall_sec"]), str(r["name"])) for r in cells.values()), reverse=True)
    rep.lines.append("slowest: " + ", ".join(f"{n} {w:.0f} s" for w, n in walls[:5]))
    rep.lines.append(f"total CPU time of the cells: {sum(w for w, _ in walls) / 3600:.2f} h")
    return rep


def check_compare(groups: list[list[Path]]) -> Report:
    rep = Report()
    _, ref, missing = load(groups[0])
    rep.fails += [f"{d}: results.json missing (the run did not finish)" for d in missing]
    for g in groups[1:]:
        metas, other, missing = load(g)
        rep.fails += [f"{d}: results.json missing (the run did not finish)" for d in missing]
        if not metas:
            continue
        label = f"{metas[0]['cmake'].get('CMAKE_CXX_COMPILER', '?')} ARCHE_MARCH={metas[0]['cmake'].get('ARCHE_MARCH', '?')!r}"
        common = [n for n in other if n in ref]
        worst = 0.0
        span = 0.0
        for n in common:
            r = cross_diff(ref[n], other[n])
            if r is None:
                if n not in KNOWN_FAILURES:
                    rep.warns.append(f"{n}: J3 not compared ({label}; a trace is missing or too short)")
                continue
            worst = max(worst, r.max_equal_nh)
            span = max(span, r.shifted_span)
            msg = (
                f"J3 max|dlogT| at equal n_H {r.max_equal_nh:.3e} dex, not explained by a shift <= "
                f"{CROSS_SHIFT_DEX:g} dex (within {CROSS_SHIFT_TOL_DEX:g} dex): {r.max_unexplained:.3e} dex, "
                f"shifted span {r.shifted_span:.4f} dex"
            )
            classify(rep, f"{n} (J3 {label})", [] if r.ok else [msg])
            if r.ok and r.max_equal_nh > CROSS_TOL_DEX:
                rep.lines.append(f"{n} (J3 {label}): accepted as a shift; {msg}")
        rep.lines.append(
            f"J3 {label}: {len(common)} common cells, max|dlogT| at equal n_H = {worst:.3e} dex, "
            f"largest shifted span {span:.4f} dex (tolerance {CROSS_TOL_DEX:g} dex; a shift <= {CROSS_SHIFT_DEX:g} dex "
            f"within {CROSS_SHIFT_TOL_DEX:g} dex over <= {CROSS_SHIFT_SPAN_DEX:g} dex)"
        )
    return rep


def render(rep: Report, title: str) -> str:
    status = "FAIL" if rep.fails else "PASS"
    out = [f"## {title}: {status}", ""]
    out += [f"- {ln}" for ln in rep.lines]
    for head, items in (("FAIL", rep.fails), ("XFAIL", rep.xfails), ("XPASS", rep.xpasses), ("WARN", rep.warns)):
        if items:
            out += ["", f"**{head}** ({len(items)})", ""] + [f"- {s}" for s in items]
    return "\n".join(out) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="*", help="result directories (with --compare: comma-separated groups)")
    ap.add_argument("--self-test", action="store_true", help="check J3 on synthetic tracks and exit")
    ap.add_argument("--compare", action="store_true", help="J3: compare each group with the first")
    ap.add_argument("--summary", type=Path, default=None, help="append the Markdown report to this file")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()
    if args.self_test:
        return self_test()
    if not args.dirs:
        ap.error("no result directories given")
    if args.compare:
        groups = [[Path(p) for p in g.split(",")] for g in args.dirs]
        if len(groups) < 2:
            ap.error("--compare needs a reference and at least one other group")
        rep = check_compare(groups)
        text = render(rep, args.title or "cross-build comparison")
    else:
        rep = check_build([Path(p) for p in args.dirs])
        text = render(rep, args.title or "collapse cells")
    print(text)
    if args.summary is not None:
        with args.summary.open("a") as f:
            f.write(text + "\n")
    return 1 if rep.fails else 0


if __name__ == "__main__":
    sys.exit(main())
