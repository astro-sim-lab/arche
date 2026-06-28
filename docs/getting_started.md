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
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

### Intel compiler (optional)

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icpx -DCMAKE_C_COMPILER=icx
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

Both networks use the full model by default. Add `--minimal` to run the compact
primordial network (`arche_collapse_prim_minimal`), and/or `--metal-minimal` to
run the compact metal-grain network (`Nakauchi2021_Minimal`, 40 species, the
`arche_collapse_metal_minimal` binary). The minimal variants take the same
CLI/env parameters as their full counterparts and reuse the same HDF5 tables;
their output stem gains a `_min` tag (with a non-`_min` symlink alias for the
plot/resample tools). The compact networks are also reachable from the library
(`prim_minimal_*` / `metal_minimal_*` entry points, or the `"…_Minimal"` registry
names — see `docs/api_reference.md`).

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

## Using the library only (external application)

To drive the chemistry network from your own application — without the bundled
`*_collapse` / `*_chemistry` executables — link **`libarche.a`** and include one
stable header (`arche_api.h` for C++, `arche_capi.h` for C). No Eigen or HDF5
headers appear in your own sources.

Build, install (`find_package(arche)`), and consume (`find_package` /
`add_subdirectory` / `FetchContent`) instructions, plus standalone runnable
examples, live in **[`integration/`](../integration/README.md)**.

## More details

- Full parameter reference: `docs/parameters.md`
- Plot tool reference: `tools/analyze_collapse.md`
- Resample tool reference: `tools/resample_collapse.md`
- API reference: `docs/api_reference.md`
- Output schema: `docs/output_schema.md`
