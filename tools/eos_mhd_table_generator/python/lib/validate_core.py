"""
Core validation logic for EOS-table files.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .table_io import EOS_COLS, monotonic_non_decreasing, read_table, split_blocks_by_ii


def validate_table(path: Path) -> None:
    arr = read_table(path, strict_cols=True)
    finite = np.isfinite(arr).all()
    if not finite:
        raise ValueError("table has non-finite values")

    ii = arr[:, EOS_COLS["ii"]].astype(int)
    jj = arr[:, EOS_COLS["jj"]].astype(int)
    nt = arr[:, EOS_COLS["nt"]]
    bt = arr[:, EOS_COLS["bt"]]

    if np.min(ii) < 1 or np.min(jj) < 1:
        raise ValueError("ii/jj must be >= 1")

    spans = split_blocks_by_ii(ii)
    if not spans:
        raise ValueError("could not detect row blocks")

    for bi, (s, e) in enumerate(spans, start=1):
        i_blk = ii[s:e]
        nt_blk = nt[s:e]
        if i_blk[0] != 1:
            raise ValueError(f"block {bi}: ii must start at 1")
        if not monotonic_non_decreasing(i_blk):
            raise ValueError(f"block {bi}: ii not monotonic")
        if not monotonic_non_decreasing(nt_blk):
            raise ValueError(f"block {bi}: nt not monotonic")

    bt_block = []
    jj_block = []
    for s, _e in spans:
        bt_block.append(bt[s])
        jj_block.append(jj[s])
    bt_block = np.asarray(bt_block, dtype=float)
    jj_block = np.asarray(jj_block, dtype=int)
    if not monotonic_non_decreasing(bt_block):
        raise ValueError("bt block sequence is not monotonic")
    if np.min(jj_block) < 1:
        raise ValueError("jj block labels must be >= 1")

    print("OK")
    print(f"  rows           : {arr.shape[0]}")
    print(f"  unique ii      : {len(np.unique(ii))} ({np.min(ii)}..{np.max(ii)})")
    print(f"  unique jj      : {len(np.unique(jj))} ({np.min(jj)}..{np.max(jj)})")
    print(f"  blocks         : {len(spans)}")
    print(f"  nt range       : {np.min(nt):.4f} .. {np.max(nt):.4f}")
    print(f"  bt block range : {np.min(bt_block):.4f} .. {np.max(bt_block):.4f}")

