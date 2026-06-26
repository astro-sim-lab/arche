# Using ARCHE as a library

This directory shows how to drive the ARCHE chemistry network from an **external
application** — without the bundled `*_collapse` / `*_chemistry` executables or
`run_collapse.sh`. You link **`libarche.a`** and include one stable header; no
Eigen or HDF5 headers appear in your own sources (the kernel's linear solver and
the HDF5 table loader are compiled inside the library and reached only through
opaque handles).

| Public header | Language | Entry points |
|---|---|---|
| `arche_api.h` | C++ (`namespace arche`) | RAII handles, `chem_full_step_*`, runtime `Model` dispatch. Errors via C++ exceptions. |
| `arche_capi.h` | C ABI (`extern "C"`) | Mirror C structs, `arche_*` functions, status codes + `arche_last_error()`. For C / Fortran / `ctypes`. |

The CMake target is `arche` (artefact `libarche.a`), exposed to consumers as the
namespaced alias **`arche::arche`**. Its usage requirements (include path,
transitive HDF5 link, Eigen) are `PUBLIC`, so linking the target is all you need.

## Contents

- [`cpp_minimal/`](cpp_minimal) — standalone CMake project, C++ via `arche_api.h`.
- [`c_minimal/`](c_minimal) — standalone CMake project, C via `arche_capi.h`.

Both `find_package(arche)` and link `arche::arche`; neither references Eigen or
HDF5 in its own build description.

## 1. Build and install ARCHE

From the repository root, build the library and install it to a prefix of your
choice (here a throwaway `_prefix` next to the build tree):

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
cmake --install build --prefix "$PWD/_prefix"
```

> **Build your consumer with the same compiler/toolchain as `libarche.a`.** The
> archive carries compiled object code, so it is not compiler-agnostic. If ARCHE
> was built with Intel oneAPI (`icpx`/`icx` — required for the bit-for-bit verify
> flow), pass `-DCMAKE_CXX_COMPILER=icpx -DCMAKE_C_COMPILER=icx` when configuring
> the consumer too; otherwise a GNU `g++` link fails with undefined references to
> Intel runtime helpers (`__cpu_core_type`, `__detect_cpu_core_type`, …). Build
> ARCHE with GNU if you need GNU consumers.

This installs:

```
_prefix/
├─ include/arche/api/{arche_api.h, arche_capi.h}
│  └─ core/{state.h, species_index.h, species_catalog.h}
└─ lib/
   ├─ libarche.a
   └─ cmake/arche/{archeConfig.cmake, archeConfigVersion.cmake, archeTargets.cmake}
```

Only the lightweight facade headers are installed — the closure reachable from
`api/arche_api.h` (`-> core/state.h -> core/species_index.h -> core/species_catalog.h`)
plus the standalone `api/arche_capi.h`. The
Eigen/HDF5-dependent kernel headers (`solve/chemistry.h`, `core/newton.h`,
`kinetics/topology.h`, …) are intentionally **not** shipped, preserving the
"no Eigen/HDF5 in your own sources" contract for installed consumers.

Set `-DARCHE_INSTALL=OFF` to skip generating these rules.

## 2. Consume it

### A. `find_package` (against an installed ARCHE)

```cmake
cmake_minimum_required(VERSION 3.18)
project(my_app LANGUAGES CXX)

find_package(arche REQUIRED)            # pulls in Eigen3 + HDF5 transitively

add_executable(my_app my_app.cc)
target_link_libraries(my_app PRIVATE arche::arche)
```

Point the consumer at the install prefix with `CMAKE_PREFIX_PATH`:

```bash
cmake -B build -S . -DCMAKE_PREFIX_PATH=/path/to/_prefix
cmake --build build
```

> `archeConfig.cmake` re-runs `find_dependency(Eigen3)` and
> `find_dependency(HDF5 COMPONENTS C)`, so those must be discoverable on the
> consuming machine. The HDF5 library is recorded by absolute path, so the
> install is not relocatable across machines with different HDF5 layouts —
> re-install on the target machine if needed.

### B. `add_subdirectory` (vendored source, no install)

```cmake
add_subdirectory(path/to/arche arche-build)
target_link_libraries(my_app PRIVATE arche::arche)
```

### C. `FetchContent`

```cmake
include(FetchContent)
FetchContent_Declare(arche
    GIT_REPOSITORY https://github.com/astro-sim-lab/arche.git
    GIT_TAG        main)
FetchContent_MakeAvailable(arche)
target_link_libraries(my_app PRIVATE arche::arche)
```

For B and C you may want `-DARCHE_INSTALL=OFF` so the vendored copy does not add
its install rules to your project.

## 3. Build and run the examples here

After step 1 (install to `$PWD/_prefix`), from the repository root:

```bash
# C++ example
cmake -B integration/cpp_minimal/build -S integration/cpp_minimal \
      -DCMAKE_PREFIX_PATH="$PWD/_prefix" -DCMAKE_BUILD_TYPE=Release
cmake --build integration/cpp_minimal/build
./integration/cpp_minimal/build/cpp_minimal data/primordial.h5

# C example
cmake -B integration/c_minimal/build -S integration/c_minimal \
      -DCMAKE_PREFIX_PATH="$PWD/_prefix" -DCMAKE_BUILD_TYPE=Release
cmake --build integration/c_minimal/build
./integration/c_minimal/build/c_minimal data/primordial.h5
```

Each prints one `chem_full_step` result with `solver_failed=0`.

## Coupling to a fluid / hydro code

The intended use is: keep your own per-cell **species abundance array** plus
`nH`, temperature and internal energy, hand them to ARCHE each step, and read the
updated abundances and the net cooling/heating rate back out. Concretely, per
cell per step you

- **set inputs** on the cell: `nH` `[cm^-3]`, `T_K` `[K]`, `mu`, `gamma`, and the
  abundance vector `y[i] = n(species i) / nH` (dimensionless, relative to H);
  plus a `ChemShielding` (attenuated CR rate, continuum optical depth, optional
  LW/X-ray field) and a `ChemParams` (`T_rad`, and for metal `Z_metal`,
  `T_gr_K`);
- **get outputs** back: ARCHE updates `y`, `mu`, `gamma` in place and returns a
  `ChemFullRates` whose `Lambda_net` `[erg g^-1 s^-1]` is the net cooling for your
  energy equation (`de/dt = -Lambda_net` on the specific internal energy), along
  with the per-process cooling breakdown and the opacities (`k_gas`, `k_gr`) you
  feed into the next step's `tau_cnt`.

**Which slot of `y[]` is which species, and what to initialise it to**, is
documented per model — including a name → index pattern so you never hard-code
literal indices — in [`docs/api_reference.md`](../docs/api_reference.md) (§3
species tables, §4 initialisation, §6 return values). The two index-free ways to
fill `y[]` are: named C++ enums from `core/species_index.h`, or a runtime name→index
map built from `model_species()` (works from C too).

## Minimal API shape

```cpp
#include "api/arche_api.h"   // one header is enough: it already pulls in the
                             // public species metadata (core/species_index.h:
                             // zero_metal::Sp index enums, abundance_ref::).
                             // Include "core/species_index.h" explicitly only if
                             // you prefer include-what-you-use hygiene.
using namespace arche;
namespace zm = zero_metal;

PrimTablePtr tbl  = load_prim_table_owned("data/primordial.h5");
PrimCellPtr  cell = prim_cell_create_owned();

ChemStateZM& s = prim_cell_state(*cell);
s.nH = 1.0e4; s.T_K = 200.0; s.mu = 1.22; s.gamma = 5.0 / 3.0;
s.y.fill(0.0);                     // zero first, then set what you know
s.y[zm::H]  = 1.0 - 1.0e-4;        // atomic H (~all of nH)
s.y[zm::e]  = 1.0e-4;              // electrons
s.y[zm::Hp] = 1.0e-4;             // H+  (balances e-)
s.y[zm::He] = abundance_ref::yHe;  // He, n(He)/nH

ChemParams    params; params.T_rad = 2.725;
ChemShielding shield; shield.zeta  = 1.0e-17;

ChemFullRates r = chem_full_step_prim(*cell, /*dt=*/1.0e8, params, shield, *tbl);
// r.Lambda_net → net cooling; s.y / s.mu / s.gamma now updated.
```

Indexing `y[]` by name (`zm::H`, `zm::Hp`, …) instead of by literal `0`/`3` keeps
the initialisation explicit and model-correct — the index order differs between
models. Swap `prim_*` for `metal_*` (and `ChemStateZM`/23 species + `zm::` for
`ChemStateMG`/89 + `metal_grain::`) for the metal+grain network, or use the
runtime `Model` / `Cell` / `Table` dispatch layer to pick the model at runtime.
The two compact networks have their own per-model entry points in the same shape:
`prim_minimal_*` (`ChemStatePrimMinimal`/15 + `zero_metal_minimal::`) loads
`primordial.h5`, and `metal_minimal_*` (`ChemStateMetalMinimal`/40 +
`metal_grain_minimal::`) loads `metal_grain.h5`. See
[`cpp_minimal/main.cc`](cpp_minimal/main.cc) and
[`c_minimal/main.c`](c_minimal/main.c) for the full, runnable versions.

## Name-based model selection (species metadata)

For a host that chooses the network by name and maps its own columns without
hard-coding the species count, the registry API selects a model by name and
exposes its species count and host-facing names. The handle owns its reaction
table; per-cell stepping is numerically identical to the per-model entries.
Registered names: `Nakauchi2019` (23 species), `Nakauchi2019_Minimal`
(compact 15 species), `Nakauchi2021` (89 species), `Nakauchi2021_Minimal`
(compact 40 species).

`arche_model_species()` returns the species names **in `y[]` order**, so a host
can build a name → index map once and address every column by species instead of
by literal index (recommended over the `y[0]`/`y[3]` form shown below — the order
is model-specific):

```c
/* C ABI (arche_capi.h) */
ArcheModelRuntime* m = arche_model_create("Nakauchi2019_Minimal",
                                          "data/primordial.h5");
int nsp = arche_model_n_species(m);            /* 15 */
const char* const* names = arche_model_species(m);  /* {"H","H2","e","H+",...} */

ArcheModelCell* c = arche_model_cell_create(m);
double* y = arche_model_cell_y(c);             /* length nsp */
for (int i = 0; i < nsp; ++i) y[i] = 0.0;
/* index by name (build a lookup from `names`); literal indices shown for brevity: */
y[0] = 1.0 - 1.0e-4; y[2] = 1.0e-4; y[3] = 1.0e-4; y[7] = 8.33e-2;  /* H, e, H+, He */
arche_model_cell_set_scalars(c, 1.0e4, 200.0, 1.22, 5.0/3.0);

/* arche_chem_*_init set the library defaults (NOT the same as = {0}:
 * esc_cnt defaults to 1.0 = optically thin, E_X_eV to 300, ...). */
ArcheChemParams p;    arche_chem_params_init(&p);    p.T_rad = 2.725;
ArcheChemShielding s; arche_chem_shielding_init(&s); s.zeta = 1.0e-17;
ArcheChemFullRates r;
arche_model_step(m, c, 1.0e8, &p, &s, &r);

arche_model_cell_free(c);
arche_model_destroy(m);
```

The C++ equivalent (`arche_api.h`) is `model_create` / `model_n_species` /
`model_species` / `model_cell_create` / `model_cell_y` / `model_step`, with
`ModelRuntimePtr` / `ModelCellPtr` RAII wrappers. This layer is additive: the
per-model entries and the `Model`-enum dispatch above are unchanged.

## Linking by hand (non-CMake)

`libarche.a` already contains the compiled Eigen kernel; you only need the
installed include directory at compile time and the HDF5 C library (plus the C++
runtime) at link time:

```bash
# C++ consumer
g++ -std=c++17 -I _prefix/include/arche my_app.cc _prefix/lib/libarche.a -lhdf5 -o my_app

# C consumer (compile as C, link with g++ for the C++ runtime)
gcc -std=c11 -I _prefix/include/arche -c my_app.c -o my_app.o
g++ my_app.o _prefix/lib/libarche.a -lhdf5 -o my_app
```

## Runtime data

Both stepping APIs need a reaction-table HDF5 file, loaded at startup via
`load_prim_table` / `load_metal_table` (C++) or `arche_*_table_load` (C). The
prebuilt tables ship in `data/primordial.h5` and `data/metal_grain.h5`.

See also `docs/api_reference.md` (full API) and `docs/output_schema.md`.
