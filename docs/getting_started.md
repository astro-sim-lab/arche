# Getting Started

This document contains setup and detailed execution guidance moved out of `README.md`.

## Prerequisites

| Component | Requirement |
|---|---|
| C++ compiler | C++17 or later (GCC >= 9 / Clang >= 9) |
| CMake | >= 3.15 |
| Eigen3 | Required (CMake package: `Eigen3`) |
| HDF5 (C library) | Required by `prim_collapse` and `metal_collapse` |
| Python | >= 3.11 |
| Python packages | `numpy >= 2.0`, `h5py`; `matplotlib` required for plotting |

## Build

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icpx -DCMAKE_C_COMPILER=icx # Intel Compiler
cmake --build build -j$(nproc)
```

Binaries are generated under `build/src/apps/...`.

## Python setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r tools/requirements.txt
```

## Core workflow (`run_collapse.sh`)

Default pipeline:

1. Build
2. Collapse simulation
3. Figure generation (`tools/python/analyze_collapse.py`)
4. Resampling (`tools/python/resample_collapse.py`)

Basic run:

```bash
bash run_collapse.sh --prim-zeta0 1e-17 --metal-zeta0 1e-17 --metal-z-metal 1e-3
```

Useful variants:

```bash
# Skip build
bash run_collapse.sh --no-build --prim-zeta0 1e-17 --metal-zeta0 1e-17 --metal-z-metal 1e-3

# Plot only (reuse existing HDF5)
bash run_collapse.sh --no-build --no-prim --no-metal --prim-zeta0 1e-17 --metal-zeta0 1e-17 --metal-z-metal 1e-3

# Disable config and use CLI/env only
bash run_collapse.sh --no-config --prim-zeta0 1e-17 --metal-zeta0 1e-17 --metal-z-metal 1e-3
```

## Parameter precedence

In `run_collapse.sh`, precedence is:

1. CLI options
2. Environment variables
3. Config file (`params/default.conf` or `--config`)

## More details

- Full parameter reference: `docs/parameters.md`
- Plot tool reference: `tools/analyze_collapse.md`
- Resample tool reference: `tools/resample_collapse.md`
- API reference: `docs/api_reference.md`
- Output schema: `docs/output_schema.md`
