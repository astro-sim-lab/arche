#!/usr/bin/env python3
"""
report_compare_matrix.py - Evaluate compare_matrix_summary.csv with pass/fail rules.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _f(row: dict[str, str], key: str) -> float:
    v = row.get(key, "").strip()
    return float(v) if v else float("nan")


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Evaluate matrix comparison summary CSV.")
    p.add_argument("--summary-csv", default="tools/eos_mhd_table_generator/results/compare_matrix_summary.csv")
    p.add_argument("--out-csv", default="tools/eos_mhd_table_generator/results/compare_matrix_eval.csv")
    p.add_argument("--strict", action="store_true", help="Return non-zero if any case is fail or missing.")
    p.add_argument("--profile", choices=("todo", "screening", "strict"), default="todo",
                   help="Threshold preset profile.")
    p.add_argument("--tt-corr-min", type=float, default=None)
    p.add_argument("--mut-corr-min", type=float, default=None)
    p.add_argument("--pt-corr-min", type=float, default=None)
    p.add_argument("--tt-med-abs-max", type=float, default=None)
    p.add_argument("--mut-med-abs-max", type=float, default=None)
    p.add_argument("--pt-med-abs-max", type=float, default=None)
    p.add_argument("--etaoh-corr-min", type=float, default=None, help="Minimum corr for etaO and etaH.")
    p.add_argument("--etaa-q95-abs-max", type=float, default=None, help="Warn threshold for |etaA_q95|.")
    return p


def main() -> int:
    args = _parser().parse_args()
    src = Path(args.summary_csv)
    dst = Path(args.out_csv)
    rows: list[dict[str, str]] = []
    with src.open("r", encoding="ascii", newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit(f"ERROR: no rows in {src}")

    presets = {
        "todo": dict(tt_corr_min=0.97, mut_corr_min=0.97, pt_corr_min=0.97,
                     tt_med_abs_max=0.15, mut_med_abs_max=0.15, pt_med_abs_max=0.15,
                     etaoh_corr_min=0.90, etaa_q95_abs_max=15.0),
        "screening": dict(tt_corr_min=0.95, mut_corr_min=0.95, pt_corr_min=0.95,
                          tt_med_abs_max=0.20, mut_med_abs_max=0.20, pt_med_abs_max=0.20,
                          etaoh_corr_min=0.85, etaa_q95_abs_max=20.0),
        "strict": dict(tt_corr_min=0.98, mut_corr_min=0.98, pt_corr_min=0.98,
                       tt_med_abs_max=0.10, mut_med_abs_max=0.10, pt_med_abs_max=0.10,
                       etaoh_corr_min=0.93, etaa_q95_abs_max=10.0),
    }
    th = presets[args.profile]
    for k in th:
        v = getattr(args, k)
        if v is not None:
            th[k] = float(v)

    out_fields = list(rows[0].keys()) + ["eval_status", "eval_reason"]
    out_rows: list[dict[str, str]] = []

    n_ok = 0
    n_warn = 0
    n_fail = 0
    n_missing = 0
    for r in rows:
        st = r.get("status", "")
        if st != "ok":
            if st.startswith("missing"):
                eval_status = "missing"
                eval_reason = st
                n_missing += 1
            else:
                eval_status = "fail"
                eval_reason = st
                n_fail += 1
        else:
            reasons: list[str] = []
            if _f(r, "tt_corr") < th["tt_corr_min"]:
                reasons.append("tt_corr")
            if _f(r, "mut_corr") < th["mut_corr_min"]:
                reasons.append("mut_corr")
            if _f(r, "pt_corr") < th["pt_corr_min"]:
                reasons.append("pt_corr")
            if _f(r, "tt_med_abs") > th["tt_med_abs_max"]:
                reasons.append("tt_med_abs")
            if _f(r, "mut_med_abs") > th["mut_med_abs_max"]:
                reasons.append("mut_med_abs")
            if _f(r, "pt_med_abs") > th["pt_med_abs_max"]:
                reasons.append("pt_med_abs")
            if _f(r, "etaO_corr") < th["etaoh_corr_min"]:
                reasons.append("etaO_corr")
            if _f(r, "etaH_corr") < th["etaoh_corr_min"]:
                reasons.append("etaH_corr")

            q95a = abs(_f(r, "etaA_q95"))
            if reasons:
                eval_status = "fail"
                eval_reason = ",".join(reasons)
                n_fail += 1
            elif q95a > th["etaa_q95_abs_max"]:
                eval_status = "warn"
                eval_reason = "etaA_q95_outlier"
                n_warn += 1
            else:
                eval_status = "ok"
                eval_reason = ""
                n_ok += 1

        ro = dict(r)
        ro["eval_status"] = eval_status
        ro["eval_reason"] = eval_reason
        out_rows.append(ro)

    dst.parent.mkdir(parents=True, exist_ok=True)
    with dst.open("w", encoding="ascii", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_fields)
        w.writeheader()
        w.writerows(out_rows)

    print(f"[OK] eval_csv={dst}")
    print(f"[OK] counts: ok={n_ok} warn={n_warn} fail={n_fail} missing={n_missing}")
    print(f"[OK] profile={args.profile} thresholds={th}")
    if args.strict and (n_warn > 0 or n_fail > 0 or n_missing > 0):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
