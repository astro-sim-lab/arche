# EOS-MHD Table Generator (WIP)

`tools/eos_mhd_table_generator` は、`arche-dev` の `metal_collapse` 出力（HDF5）から
MHD 用 EOS テーブル（10列 `.d`）を生成する補助ツールです。

## CAUTION

- 本ツールは開発中です。参考元である Higuchi et al. (2019) コードとの厳密一致は未達成です。
- 現在の目標は「定性的整合」です。厳密比較が必要な用途では必ず比較・可視化を実施してください。
- 参照コード側に最終 `.d` 組み立て手順の欠損がある可能性があり、再現限界があります。

## 関連論文（ADS）

- Higuchi et al. 2018 (MNRAS 475, 3331): https://ui.adsabs.harvard.edu/abs/2018MNRAS.475.3331H/abstract
- Higuchi et al. 2019 (MNRAS 486, 3741): https://ui.adsabs.harvard.edu/#abs/2019MNRAS.486.3741H/abstract

## Quickstart (5 min)

前提:
- プロジェクトルートで実行
- `metal_collapse` をビルド済み

```bash
cmake -S . -B build-codex
cmake --build build-codex --target metal_collapse -j
```

最短実行（計算 -> 検証 -> 比較/可視化）:

```bash
bash tools/eos_mhd_table_generator/run_eos_table.sh quickstart
```

既定出力:
- 生成テーブル: `tools/eos_mhd_table_generator/results/CR1e-2_sh1e-2_Z1e-3.d`
- 比較プロット: `tools/eos_mhd_table_generator/results/plots/quickstart/`

## よく使うコマンド

```bash
bash tools/eos_mhd_table_generator/run_eos_table.sh case --cr-scale 1e-2 --sh-scale 1e-2 --z-metal 1e-3
bash tools/eos_mhd_table_generator/run_eos_table.sh validate
bash tools/eos_mhd_table_generator/run_eos_table.sh compare
bash tools/eos_mhd_table_generator/run_eos_table.sh plot
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix --cases tools/eos_mhd_table_generator/params/cases_core.csv
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix --cases tools/eos_mhd_table_generator/params/cases_extended.csv --skip-existing --metal-log-dir tools/eos_mhd_table_generator/results/logs --timeout-per-case 7200
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-status --cases tools/eos_mhd_table_generator/params/cases_extended.csv
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-compare --cases tools/eos_mhd_table_generator/params/cases_core.csv
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-state-native --cases tools/eos_mhd_table_generator/params/cases_core.csv
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-state-native --cases tools/eos_mhd_table_generator/params/cases_extended.csv --skip-build
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-report --summary-csv tools/eos_mhd_table_generator/results/compare_matrix_summary.csv --profile todo
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-sync --eval-csv tools/eos_mhd_table_generator/results/compare_matrix_eval.csv
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-plot --cases tools/eos_mhd_table_generator/params/cases_core.csv --timestamp-tag
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-runtime-profile --runtime-csv tools/eos_mhd_table_generator/results/matrix_runtime.csv
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-build-runtime --cases tools/eos_mhd_table_generator/params/cases_extended.csv --out-csv tools/eos_mhd_table_generator/results/matrix_runtime_build.csv --timeout-per-case 60
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-etaa --cases tools/eos_mhd_table_generator/params/cases_core.csv --nt-min 19 --abs-diff-min 3
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-etaa-onset --cases tools/eos_mhd_table_generator/params/cases_core.csv --abs-diff-threshold 3 --nt-min 14
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-etaa-summary --in-csv tools/eos_mhd_table_generator/results/etaa_outliers.csv
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-md --out-md .dev/low-metal-eos/matrix_report.md
bash tools/eos_mhd_table_generator/run_eos_table.sh selftest
```

## 長時間計算時の運用フロー

```bash
# 1) 計算 (再開しやすい設定)
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix \
  --cases tools/eos_mhd_table_generator/params/cases_extended.csv \
  --skip-existing \
  --metal-log-dir tools/eos_mhd_table_generator/results/logs \
  --timeout-per-case 7200

# 2) 比較と判定
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-status --cases tools/eos_mhd_table_generator/params/cases_extended.csv
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-compare --cases tools/eos_mhd_table_generator/params/cases_extended.csv
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-state-native --cases tools/eos_mhd_table_generator/params/cases_core.csv
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-report --summary-csv tools/eos_mhd_table_generator/results/compare_matrix_summary.csv --profile todo

# 3) 低速化とetaA外れ値を確認
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-runtime-profile --runtime-csv tools/eos_mhd_table_generator/results/matrix_runtime.csv
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-build-runtime \
  --cases tools/eos_mhd_table_generator/params/cases_extended.csv \
  --out-csv tools/eos_mhd_table_generator/results/matrix_runtime_build.csv \
  --timeout-per-case 60
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-etaa --cases tools/eos_mhd_table_generator/params/cases_extended.csv --nt-min 19 --abs-diff-min 3
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-etaa-onset --cases tools/eos_mhd_table_generator/params/cases_extended.csv --abs-diff-threshold 3 --nt-min 14
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-etaa-summary --in-csv tools/eos_mhd_table_generator/results/etaa_outliers.csv

# 4) 可視化とTODO反映
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-plot --cases tools/eos_mhd_table_generator/params/cases_extended.csv --timestamp-tag
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-sync --eval-csv tools/eos_mhd_table_generator/results/compare_matrix_eval.csv
bash tools/eos_mhd_table_generator/run_eos_table.sh matrix-md --out-md .dev/low-metal-eos/matrix_report.md
```

ケース定義:
- 主要確認用: `tools/eos_mhd_table_generator/params/cases_core.csv`
- 拡張確認用: `tools/eos_mhd_table_generator/params/cases_extended.csv`

## 主要パラメータ（どこで変えるか）

| 目的 | CLI（推奨） | `params/default.conf` | 既定 |
|---|---|---|---|
| 宇宙線基準電離率 | `--cr-scale` / `--zeta0` | `EOS_CASE_CR_SCALE` / `EOS_CASE_ZETA0` | quickstartでは `1e-2` |
| 金属量 | `--z-metal` | `EOS_CASE_Z_METAL` | quickstartでは `1e-3` |
| SRA/LRA係数 | `--sh-scale` または `--sra-rate --lra-rate` | `EOS_CASE_SH_SCALE` / `EOS_CASE_SRA_RATE` / `EOS_CASE_LRA_RATE` | `0.0` |
| etaモデル | `--eta-model` | `EOS_CASE_ETA_MODEL` | `legacy` |
| eta再サンプル | `--eta-resample` | `EOS_CASE_ETA_RESAMPLE` | `native` |
| 出力先 `.d` | `--output-dir --output-name` | `EOS_CASE_OUTPUT_DIR` / `EOS_CASE_OUTPUT_NAME` | `tools/eos_mhd_table_generator/results/...` |
| 軸範囲（上級） | `--n-log-* --b-log-*` | `EOS_N_LOG_*` / `EOS_B_LOG_*` | `n:[0.1,22], b:[-20,4]` |

## ディレクトリ構成

- 実行: `python/app/`
- ライブラリ: `python/lib/`
- テスト: `python/test/`
- 設定: `params/`
- 参照テーブル（電荷・質量）: `include/`

## トラブルシュート

- `metal_collapse binary not found`:
  - `cmake --build build-codex --target metal_collapse -j` を実行
- 比較が失敗する:
  - まず `validate` を通して形式エラーを除外
  - 次に `plot` で `tt/mut/pt/eta*` のどこがずれているか確認
- `matrix-state-native` が長い:
  - 先に `matrix-state-native` を通常実行して `results/state_native/*.d` を作成
  - その後 `--skip-build` を付けて比較CSVのみ再生成

## 開発ドキュメント

開発再開・設計判断・未解決課題は以下を参照:
- `.dev/TODO_EOS-MHD.md`
- `.dev/low-metal-eos/eos_mhd_spec.md`
