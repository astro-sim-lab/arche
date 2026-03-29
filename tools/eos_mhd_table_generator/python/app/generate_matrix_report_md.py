#!/usr/bin/env python3
"""
generate_matrix_report_md.py - Generate consolidated markdown report from matrix artifacts.
"""

from __future__ import annotations

import argparse
import csv
import math
from datetime import datetime
from pathlib import Path


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="ascii", newline="") as f:
        return list(csv.DictReader(f))


def _top(rows: list[dict[str, str]], key: str, n: int = 5) -> list[dict[str, str]]:
    vals = []
    for r in rows:
        try:
            vals.append((float(r.get(key, "")), r))
        except ValueError:
            continue
    vals.sort(key=lambda x: x[0], reverse=True)
    return [r for _, r in vals[:n]]


def _fval(row: dict[str, str], key: str) -> float | None:
    s = row.get(key, "")
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Generate EOS-MHD matrix markdown report.")
    p.add_argument("--status-csv", default="tools/eos_mhd_table_generator/results/matrix_status.csv")
    p.add_argument("--eval-csv", default="tools/eos_mhd_table_generator/results/compare_matrix_eval.csv")
    p.add_argument("--onset-csv", default="tools/eos_mhd_table_generator/results/etaa_onset.csv")
    p.add_argument("--outlier-summary-csv", default="tools/eos_mhd_table_generator/results/etaa_outliers_summary.csv")
    p.add_argument("--state-native-csv", default="tools/eos_mhd_table_generator/results/state_native_compare.csv")
    p.add_argument("--runtime-csv", default="tools/eos_mhd_table_generator/results/matrix_runtime.csv")
    p.add_argument("--runtime-density-csv", default="tools/eos_mhd_table_generator/results/matrix_runtime_density_profile.csv")
    p.add_argument("--out-md", default=".dev/low-metal-eos/matrix_report.md")
    return p


def main() -> int:
    args = _parser().parse_args()
    status = _read_csv(Path(args.status_csv))
    ev = _read_csv(Path(args.eval_csv))
    onset = _read_csv(Path(args.onset_csv))
    outsum = _read_csv(Path(args.outlier_summary_csv))
    sn = _read_csv(Path(args.state_native_csv))
    rt = _read_csv(Path(args.runtime_csv))
    rt_prof = _read_csv(Path(args.runtime_density_csv))

    n_ready = sum(1 for r in status if r.get("status") == "ready_compare")
    n_eval_ok = sum(1 for r in ev if r.get("eval_status") == "ok")
    n_eval_warn = sum(1 for r in ev if r.get("eval_status") == "warn")
    n_eval_fail = sum(1 for r in ev if r.get("eval_status") == "fail")
    n_eval_missing = sum(1 for r in ev if r.get("eval_status") == "missing")
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines = [
        "# EOS-MHD Matrix Report",
        "",
        f"- Updated: `{ts}`",
        f"- Ready for compare: `{n_ready}`",
        f"- Eval counts: `ok={n_eval_ok}` `warn={n_eval_warn}` `fail={n_eval_fail}` `missing={n_eval_missing}`",
        "",
        "## Eval Warn/Fail/Missing",
    ]
    wf = [r for r in ev if r.get("eval_status") in ("warn", "fail", "missing")]
    if not wf:
        lines.append("- none")
    else:
        for r in wf:
            lines.append(f"- `{r.get('case_request', r.get('case','?'))}`: `{r.get('eval_status')}` ({r.get('eval_reason','')})")

    lines += ["", "## etaA Onset (Top max_abs_diff)"]
    for r in _top(onset, "max_abs_diff", n=5):
        lines.append(
            f"- `{r.get('case_request','?')}` max=`{r.get('max_abs_diff','')}` onset_nt=`{r.get('onset_nt','')}` onset_bt=`{r.get('onset_bt','')}`"
        )
    if len(lines) and lines[-1] == "## etaA Onset (Top max_abs_diff)":
        lines.append("- none")

    lines += ["", "## etaA Outlier Summary"]
    if not outsum:
        lines.append("- none")
    else:
        for r in outsum:
            lines.append(
                f"- `{r.get('case_request','?')}` n=`{r.get('n_outliers','')}` max=`{r.get('max_abs_diff','')}` at nt=`{r.get('max_abs_diff_nt','')}` bt=`{r.get('max_abs_diff_bt','')}`"
            )

    lines += ["", "## State vs Native (eta-resample)"]
    if not sn:
        lines.append("- none")
    else:
        improve = 0
        worse = 0
        tie = 0
        nan = 0
        for r in sn:
            d = _fval(r, "etaA_med_abs_native_minus_state")
            if d is None or math.isnan(d):
                nan += 1
            elif d < 0.0:
                improve += 1
            elif d > 0.0:
                worse += 1
            else:
                tie += 1
        lines.append(f"- Source: `{args.state_native_csv}`")
        lines.append(f"- Coverage: `{len(sn)}` cases")
        lines.append(f"- `etaA_med_abs_native_minus_state`: 改善={improve}, 悪化={worse}, 同値={tie}, NaN={nan}")

    lines += ["", "## Phase-1 Gate (Thermal + etaO/H)"]
    phase1_pass = []
    phase1_hold = []
    for r in ev:
        case = r.get("case_request", r.get("case", "?"))
        tt_corr = _fval(r, "tt_corr")
        mut_corr = _fval(r, "mut_corr")
        pt_corr = _fval(r, "pt_corr")
        tt_med = _fval(r, "tt_med_abs")
        mut_med = _fval(r, "mut_med_abs")
        pt_med = _fval(r, "pt_med_abs")
        etao_corr = _fval(r, "etaO_corr")
        etah_corr = _fval(r, "etaH_corr")
        ok = (
            tt_corr is not None and tt_corr >= 0.97
            and mut_corr is not None and mut_corr >= 0.97
            and pt_corr is not None and pt_corr >= 0.97
            and tt_med is not None and tt_med <= 0.15
            and mut_med is not None and mut_med <= 0.15
            and pt_med is not None and pt_med <= 0.15
            and etao_corr is not None and etao_corr >= 0.90
            and etah_corr is not None and etah_corr >= 0.90
        )
        (phase1_pass if ok else phase1_hold).append(case)
    lines.append(f"- pass: `{len(phase1_pass)}`")
    if phase1_pass:
        for c in phase1_pass:
            lines.append(f"  - `{c}`")
    lines.append(f"- hold: `{len(phase1_hold)}`")
    if phase1_hold:
        for c in phase1_hold:
            lines.append(f"  - `{c}`")

    lines += ["", "## Runtime Coverage"]
    ok_rt = [r for r in rt if r.get("exit_code") == "0"]
    lines.append(f"- runtime rows (exit=0): `{len(ok_rt)}`")
    lines.append(f"- runtime profile rows: `{len(rt_prof)}`")
    if ok_rt:
        top_rt = sorted(ok_rt, key=lambda r: float(r.get("elapsed_sec", "0") or 0.0), reverse=True)[:3]
        lines.append("- top elapsed cases:")
        for r in top_rt:
            lines.append(f"  - `{r.get('case','?')}`: `{r.get('elapsed_sec','?')} s`")

    lines += ["", "## Artifacts", ""]
    for p in (
        args.status_csv,
        args.eval_csv,
        args.onset_csv,
        args.outlier_summary_csv,
        args.state_native_csv,
        args.runtime_csv,
        args.runtime_density_csv,
    ):
        lines.append(f"- `{p}`")

    out = Path(args.out_md)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[OK] report_md={out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
