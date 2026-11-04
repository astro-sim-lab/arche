# Collapse tests

One-zone collapses with initial parameters that can occur physically, each of
which must reach n_H = 1e23 cm⁻³ normally and within a wall-time gate.

## Cells

`cells.py` defines the grid **G** (160 cells) and its subset **S** (16 cells).
Networks: `pf` / `pm` = primordial full / minimal, `mf` / `mm` = metal full /
minimal. Unlisted parameters keep their defaults
([`docs/parameters.md`](../../docs/parameters.md)).

| Group | Cells | Count |
|---|---|---|
| G1 | {pf, pm} × ζ {0, 1e-17, 1e-15} s⁻¹ × f_ret {1, 3, 10}, z = 20 | 18 |
| G2 | {mf, mm} × Z {1, 1e-1, …, 1e-6} Z☉ × ζ {0, 1e-17, 1e-15} s⁻¹ × f_ret {1, 3, 10}, z = 0 | 126 |
| G3 | {mf, mm} × Z = 0, ζ = 0, z = 20 | 2 |
| G4 | {pf, pm} × J_LW21 {1e3, 1e4, 1e5} × f_ret {1, 3}, ζ = 0, z = 20 | 12 |
| G5 | {mf, mm} × Z = 1e-3 × J_LW21 = 1e5, ζ = 0 | 2 |

Primordial gas and Z = 0 use z = 20 (T_CMB = 57.2 K). S: the four networks at
ζ {0, 1e-17} s⁻¹, f_ret = 1, with Z {1e-6, 1e-3, 1} for the metal networks.

## Criteria (`check_runs.py`)

| ID | Criterion |
|---|---|
| J1 | Finished before the gate; HDF5 `exit_code` = 0 and stdout `EXIT code=0`; n_H and T finite. The applications return 0 to the shell even on an abnormal end, so the process status is not used. |
| J2 | Every step projected (`CONSERVATION ... projected`), max \|`charge_imbal`\| ≤ 1e-14. |
| J3 | Across builds, per cell: \|Δlog₁₀ T\| ≤ 1e-3 dex at equal n_H on the common log n_H range. A point above it is accepted only if shifting the other track by ≤ 0.01 dex in log n_H brings it within 3e-3 dex (both ways), and such points span ≤ 0.1 dex of log n_H in total. This admits a steep feature placed slightly apart by two builds (a 0.095 dex rapid cooling that g++ 11.4 and 13.3 start 0.0064 dex apart in n_H); a difference ≤ 3e-3 dex confined to ≤ 0.1 dex of log n_H is therefore not detected. `check_runs.py --self-test` checks the criterion on synthetic tracks. |
| J4 | Wall time ≤ 300 s on the reference machine (gate, fail); above 120 s it is reported. |
| G3 | `mf` at Z = 0 vs `pf` (ζ = 0, f_ret = 1, z = 20): max \|Δlog₁₀ T\| ≤ 1e-4 dex (`mf` contains the primordial full network). For `mm`, a different reduction, the difference is only reported. |

The gate and the report threshold are scaled by the speed of the machine: before
the cells, `run_cells.py` runs the calibration cell (`pf`, ζ = 0, f_ret = 1,
z = 20) as many times concurrently as cells will run, and divides the median
wall time by its value on the reference machine (4.3 s on an Intel Core
i5-12600KF). The scale is never below 1.

Cells listed in `cells.KNOWN_FAILURES` are reported as XFAIL, and as XPASS once
they pass.

## Running

```bash
tools/ci/run_local.sh                                   # G on g++ -march=x86-64-v3
tools/ci/run_local.sh --compiler icx --subset S         # S on icpx (source setvars.sh first)
tools/ci/run_local.sh --march "" --subset S             # S on g++ without -march (FMA off)
tools/ci/run_local.sh --march x86-64-v4 --subset S      # S on AVX-512 (needs an AVX-512 CPU)
python3 tests/collapse/check_runs.py --compare results/tests/gcc_x86-64-v3_G \
    results/tests/icx_x86-64-v3_S results/tests/gcc_nomarch_S   # J3
```

Outputs go to `results/tests/<compiler>_<march>_<subset>/`: `results.json`, and
per cell the HDF5 output, `trace.npz` (log n_H, log T of every row) and a
trimmed stdout log. The pull-request workflow
(`.github/workflows/collapse-tests.yml`) runs G on g++ `x86-64-v3` in 8 shards,
S on icpx `x86-64-v3` and on g++ without `-march`, and J3 between them.
