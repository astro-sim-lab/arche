# Getting Started

This document contains setup and detailed execution guidance moved out of `README.md`.

## Prerequisites

| Component | Requirement |
|---|---|
| C++ compiler | C++17 or later (GCC >= 9 / Clang >= 9) |
| CMake | >= 3.18 |
| Eigen3 | Required (CMake package: `Eigen3`) |
| HDF5 (C library) | Required (`libarche` reads its reaction tables from HDF5) |
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

With `icpx`, ARCHE also compiles with `-fp-model=precise`. The classic Intel
compiler (`icpc`) is not supported.

### Target CPU (`ARCHE_MARCH`)

ARCHE compiles libarche and the bundled applications with `-march=<ARCHE_MARCH>`
in every build type. The default `native` targets the CPU of the build machine.
If the machines that run ARCHE have a different CPU (for example, a login node
and compute nodes of different generations or vendors), set the instruction set
of the run machines:

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DARCHE_MARCH=x86-64-v4   # AVX-512
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DARCHE_MARCH=x86-64-v3   # AVX2
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DARCHE_MARCH=            # no -march flag
```

The level names `x86-64-v3` and `x86-64-v4` require GCC >= 11, Clang >= 12 or
`icpx`. A CPU name accepted by the compiler (e.g. `znver4`, `sapphirerapids`)
also works. These flags apply only to ARCHE's own targets: a project that links
`arche::arche` keeps its own CPU and floating-point flags. When
`CMAKE_BUILD_TYPE` is not given, ARCHE builds as `Release`.

Binaries are generated under `build/src/apps/...`: the collapse applications
(`prim_collapse`, `arche_collapse_prim_minimal`, `metal_collapse`,
`arche_collapse_metal_minimal`), the standalone chemistry demos (`prim_chem`,
`metal_chem`), and `magnetic_reynolds`, which evaluates the Ohmic and ambipolar
diffusivities and the magnetic Reynolds number along a collapse output (see
`docs/parameters.md` and `docs/magnetic_resistivity.md`).

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
run the compact metal-grain network (`Nakauchi2021_Minimal`, 44 species, the
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
- Changes between releases: `docs/release_notes.md`
