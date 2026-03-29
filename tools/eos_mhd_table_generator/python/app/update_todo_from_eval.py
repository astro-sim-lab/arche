#!/usr/bin/env python3
"""
update_todo_from_eval.py - Sync matrix evaluation summary into TODO_EOS-MHD.md.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path


START_MARK = "<!-- AUTO: MATRIX_EVAL START -->"
END_MARK = "<!-- AUTO: MATRIX_EVAL END -->"


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Update TODO_EOS-MHD.md from compare_matrix_eval.csv.")
    p.add_argument("--eval-csv", default="tools/eos_mhd_table_generator/results/compare_matrix_eval.csv")
    p.add_argument("--todo", default=".dev/TODO_EOS-MHD.md")
    return p


def _fmt_block(rows: list[dict[str, str]]) -> str:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    n_ok = sum(1 for r in rows if r.get("eval_status") == "ok")
    n_warn = sum(1 for r in rows if r.get("eval_status") == "warn")
    n_fail = sum(1 for r in rows if r.get("eval_status") == "fail")
    n_missing = sum(1 for r in rows if r.get("eval_status") == "missing")
    if n_missing > 0:
        phase = "収集中（missingケースあり）"
    elif n_fail > 0:
        phase = "解析中（metric fail あり）"
    elif n_warn > 0:
        phase = "解析中（warnのみ）"
    else:
        phase = "完了（全ケース ok）"
    lines = [
        START_MARK,
        "## Matrix Eval (Auto)",
        "",
        f"- Updated: `{ts}`",
        f"- Phase: `{phase}`",
        f"- Counts: `ok={n_ok}` `warn={n_warn}` `fail={n_fail}` `missing={n_missing}`",
        "",
        "### Warn / Fail Cases",
    ]
    wf = [r for r in rows if r.get("eval_status") in ("warn", "fail", "missing")]
    if not wf:
        lines.append("- none")
    else:
        for r in wf:
            case = r.get("case_request") or r.get("case") or "unknown"
            st = r.get("eval_status", "")
            reason = r.get("eval_reason", "")
            lines.append(f"- `{case}`: `{st}` ({reason})")
    lines.extend(
        [
            "",
            "### Source",
            "",
            "- `tools/eos_mhd_table_generator/results/compare_matrix_eval.csv`",
            END_MARK,
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    args = _parser().parse_args()
    eval_csv = Path(args.eval_csv)
    todo = Path(args.todo)
    with eval_csv.open("r", encoding="ascii", newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit(f"ERROR: no rows in {eval_csv}")
    content = todo.read_text(encoding="utf-8")
    block = _fmt_block(rows)
    if START_MARK in content and END_MARK in content:
        s = content.index(START_MARK)
        e = content.index(END_MARK) + len(END_MARK)
        new_content = content[:s] + block + content[e:]
    else:
        sep = "" if content.endswith("\n") else "\n"
        new_content = content + sep + "\n" + block
    todo.write_text(new_content, encoding="utf-8")
    print(f"[OK] updated {todo}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
