# API Reference — libarche public API

Public headers (installed under `include/arche/`):

| Header | Language | Entry points |
|---|---|---|
| `api/arche_api.h` | C++ (`namespace arche`) | RAII handles, `chem_full_step_*`, runtime `Model` / name-registry dispatch. Errors via C++ exceptions. |
| `api/arche_capi.h` | C ABI (`extern "C"`) | mirror C structs, `arche_*` functions, status codes + `arche_last_error()`. For C / Fortran / `ctypes`. |
| `core/state.h` | C++ | `ChemState`, `ChemParams`, `ChemShielding`, `ChemFullRates`, `phys::` constants. |
| `core/species_index.h`, `core/species_catalog.h` | C++ | species catalog, per-model `Sp` index enums, cosmic-abundance constants. |

> This is the reference for **driving the network from an external application**
> (e.g. coupling ARCHE to a hydro code). For build / link / `find_package`
> instructions see [`integration/README.md`](../integration/README.md). For the
> HDF5 output of the bundled collapse apps see
> [`docs/output_schema.md`](output_schema.md).

The library is linked as a single static archive `libarche.a` (CMake target
`arche::arche`). The Eigen linear solver and the HDF5 table loader are compiled
*inside* the library; your own sources include only the headers above and never
touch Eigen or HDF5.

---

## 1. The mental model

For a hydro/fluid coupling the loop is always the same:

1. **Once at startup** — load a reaction table (an HDF5 file) and create one
   reusable per-cell scratch object (the "cell") per thread.
2. **Per cell, per step** — copy the cell's gas state *into* ARCHE (density,
   temperature, the species abundance vector), call one stepping function, then
   copy the *updated* abundances and the returned cooling/heating rates back into
   your fluid fields.

ARCHE owns no global state: the table is read-only and shareable; the cell holds
the per-cell species vector plus an inter-step rate cache used by the
predictor–corrector solver.

### The abundance vector `y[]`

The chemical state is a single array `y[]` of **dimensionless abundances
relative to hydrogen**: `y[i] = n(species i) / nH`, where `nH` is the hydrogen
nucleus number density `[cm^-3]`. So `y[H] ≈ 1` for mostly-atomic gas, and a
trace ion sits at `y ≈ 10^-4`. Helium is `y[He] = n(He)/nH ≈ 0.083`. **Indices
into `y[]` are model-specific** (see §3); never assume a literal index — resolve
it from a name or a named enum (§4).

---

## 2. Models

ARCHE ships four networks. A model fixes the species set, the species *order*
in `y[]`, the reaction count, and which HDF5 table file it loads.

| Model | Registry name | C++ `Model` enum | `N_sp` | `N_react` | HDF5 file |
|---|---|---|---|---|---|
| Primordial (Nakauchi 2019) | `"Nakauchi2019"` | `Model::Prim` | 23 | 140 | `data/primordial.h5` |
| Primordial minimal | `"Nakauchi2019_Minimal"` | — | **15** | 33 | `data/primordial.h5` |
| Metal + grain (Nakauchi 2021) | `"Nakauchi2021"` | `Model::Metal` | 95 | 1200 | `data/metal_grain.h5` |
| Metal minimal | `"Nakauchi2021_Minimal"` | — | **44** | 113 | `data/metal_grain.h5` |

- **Nakauchi2019** — full primordial H/He/D/Li network (23 species).
- **Nakauchi2019_Minimal** — compact 15-species subset of the primordial
  network (H/D/Li chemistry + cosmic-ray ionization), for embedded use where the
  full network is too heavy. Loads the *same* `primordial.h5` and remaps it
  internally. He⁺ is tracked only under cosmic rays (no collisional He
  ionization), so without a CR field He⁺ stays ≈ 0.
- **Nakauchi2021** — metal + dust-grain network (95 species), adds C/O/Mg/Na/K
  chemistry, grain charge states, ice-mantle surface species, and the dust
  reservoirs of the grain evaporation (below).
- **Nakauchi2021_Minimal** — compact 44-species composed subset of the
  metal+grain network: the H/He/D core, a reduced C/O molecular network, Li, Mg,
  four grain charge states, five ice-mantle (physisorbed) species, and the dust
  reservoirs of the grain evaporation (C, O, Mg). The
  113 reactions split into 103 gas-phase (68 standard + 8 cosmic-ray + 27
  ion–grain charge transfer) and 10 grain-surface freeze-out. A first-class
  composed model — not a runtime mask — derived from the full network's keep-set;
  it loads the *same* `metal_grain.h5` (reusing its partition functions and Saha
  table) and its coupled Newton solves 45 unknowns (44 species and the
  temperature). It drops the species the reduced network omits (He⁺⁺, the
  higher C/O ions, Na, K, Gr²⁺).

**Grain evaporation (both metal models).** The metals locked in the dust are
species of the model (`C_s`, `O_s`, `Mg_s`; the full network also `Na_s`, `K_s`),
counted as their element by the conservation projection, and `Gr_v` is the
smallest grain volume reached so far relative to the volume at T → 0 (0 = not yet
set: the first step sets it from the grain temperature and releases nothing).
Each `chem_full_step` releases, from its input state, the fraction of the locked
metals that the grain volume lost since `Gr_v` (once the grain temperature
T_gr ≥ 0.975 T_pyr), at a constant rate over the interval, and removes the grain
species of every charge state in proportion to the volume, their net charge going
to the electrons (full network, while the metals are released: a net positive
charge to K⁺). The release is irreversible. The caller sets the locked amounts
initially (e.g. element abundance × (1 − gas fraction)) and carries all these
species like the others; see `solve/grain_evaporation.h`.

The two minimal models are reachable from C++ via the `prim_minimal_*` /
`metal_minimal_*` entry points and from both C++ and C via the **name registry**
(§7.3). The fixed-enum C ABI (`arche_*_prim` / `arche_*_metal`) covers only the
23- and 95-species networks.

---

## 3. Species tables

These tables are the authoritative answer to *"which slot of `y[]` is which
species, and what should I initialise it to?"* The **name** column is exactly
the string returned by `model_species()` / `arche_model_species()` (§4); the
**C++ enum** column is the named constant in `core/species_index.h` (except the
minimal models' `zero_metal_minimal::Sp` / `metal_grain_minimal::Sp`, which live
in the in-tree headers `models/primordial/minimal.h` /
`models/metal_grain/minimal.h` and are **not** part of the installed public set —
see §3.2, §3.4).

### 3.1 Nakauchi2019 — 23 species (`zero_metal::Sp`)

| idx | name | C++ enum | | idx | name | C++ enum |
|---:|---|---|---|---:|---|---|
| 0 | `H`    | `zero_metal::H`    | | 12 | `HD`    | `zero_metal::HD`   |
| 1 | `H2`   | `zero_metal::H2`   | | 13 | `D+`    | `zero_metal::Dp`   |
| 2 | `e`    | `zero_metal::e`    | | 14 | `HD+`   | `zero_metal::HDp`  |
| 3 | `H+`   | `zero_metal::Hp`   | | 15 | `D-`    | `zero_metal::Dm`   |
| 4 | `H2+`  | `zero_metal::H2p`  | | 16 | `Li`    | `zero_metal::Li`   |
| 5 | `H3+`  | `zero_metal::H3p`  | | 17 | `LiH`   | `zero_metal::LiH`  |
| 6 | `H-`   | `zero_metal::Hm`   | | 18 | `Li+`   | `zero_metal::Lip`  |
| 7 | `He`   | `zero_metal::He`   | | 19 | `Li-`   | `zero_metal::Lim`  |
| 8 | `He+`  | `zero_metal::Hep`  | | 20 | `LiH+`  | `zero_metal::LiHp` |
| 9 | `He++` | `zero_metal::Hepp` | | 21 | `Li++`  | `zero_metal::Lipp` |
| 10 | `HeH+`| `zero_metal::HeHp` | | 22 | `Li+++` | `zero_metal::Lippp`|
| 11 | `D`   | `zero_metal::D`    | |    |         |                    |

### 3.2 Nakauchi2019_Minimal — 15 species (`zero_metal_minimal::Sp`)

| idx | name | C++ enum | | idx | name | C++ enum |
|---:|---|---|---|---:|---|---|
| 0 | `H`   | `zero_metal_minimal::H`   | | 8  | `He+` | `zero_metal_minimal::Hep` |
| 1 | `H2`  | `zero_metal_minimal::H2`  | | 9  | `D`   | `zero_metal_minimal::D`   |
| 2 | `e`   | `zero_metal_minimal::e`   | | 10 | `HD`  | `zero_metal_minimal::HD`  |
| 3 | `H+`  | `zero_metal_minimal::Hp`  | | 11 | `D+`  | `zero_metal_minimal::Dp`  |
| 4 | `H2+` | `zero_metal_minimal::H2p` | | 12 | `Li`  | `zero_metal_minimal::Li`  |
| 5 | `H3+` | `zero_metal_minimal::H3p` | | 13 | `LiH` | `zero_metal_minimal::LiH` |
| 6 | `H-`  | `zero_metal_minimal::Hm`  | | 14 | `Li+` | `zero_metal_minimal::Lip` |
| 7 | `He`  | `zero_metal_minimal::He`  | |    |       |                           |

Note the first 9 slots (0–8, through `He+`) match Nakauchi2019, but slot 9 onward
differs (the minimal set drops He++, HeH+, HD+, D-, Li-, LiH+, Li++, Li+++). This
is exactly why you must not reuse one model's literal indices for another —
resolve by name (§4).

The `zero_metal_minimal::Sp` enum above is defined in the in-tree header
`models/primordial/minimal.h` (used by the bundled apps), **not** in the
installed public header `core/species_index.h`. An installed external consumer
therefore resolves the minimal model's columns by name through the registry
(§4 Pattern B, §7.3) rather than by this enum.

### 3.3 Nakauchi2021 — 95 species (`metal_grain::Sp`)

Slots 0–15 are the H/He/D core (same order as Nakauchi2019's 0–15). The C++
enum is `metal_grain::<Name>` (e.g. `metal_grain::CO`).

| idx range | family | names (in order) |
|---|---|---|
| 0–15 | H/He/D core | `H` `H2` `e` `H+` `H2+` `H3+` `H-` `He` `He+` `He++` `HeH+` `D` `HD` `D+` `HD+` `D-` |
| 16–28 | carbon | `C` `C2` `CH` `CH2` `CH3` `CH4` `C+` `C2+` `CH+` `CH2+` `CH3+` `CH4+` `CH5+` |
| 29–38 | oxygen (neutral) | `O` `O2` `OH` `CO` `H2O` `HCO` `O2H` `CO2` `H2CO` `H2O2` |
| 39–49 | oxygen (ions) | `O+` `O2+` `OH+` `CO+` `H2O+` `HCO+` `O2H+` `H3O+` `H2CO+` `HOCO+` `H2COH+` |
| 50–56 | lithium | `Li` `LiH` `Li+` `Li-` `LiH+` `Li++` `Li+++` |
| 57–62 | K / Na / Mg | `K` `K+` `Na` `Na+` `Mg` `Mg+` |
| 63–67 | grain charge | `Gr` `Gr+` `Gr2+` `Gr-` `Gr2-` |
| 68–88 | grain-surface / ice | `H_p` `H_c` `H2_p` `D_p` `D_c` `HD_p` `O_p` `O2_p` `OH_p` `CO_p` `CO2_p` `H2O_p` `HO2_p` `H2O2_p` `HCO_p` `H2CO_p` `C_p` `CH_p` `CH2_p` `CH3_p` `CH4_p` |
| 89–94 | dust reservoirs | `C_s` `O_s` `Mg_s` `Na_s` `K_s` `Gr_v` |

(`_p` = physisorbed on grain, `_c` = chemisorbed, `_s` = locked in the refractory
dust; the gas-phase parent is named without the suffix. `Gr_v` is the
grain-volume tracer of the evaporation, §2.)

### 3.4 Nakauchi2021_Minimal — 44 species (`metal_grain_minimal::Sp`)

Slots 0–11 are the H/He/D core (same order, and same first-9 layout, as the other
networks). The C++ enum is `metal_grain_minimal::<Name>` (e.g.
`metal_grain_minimal::CO`); the five ice-mantle slots use a `_p` suffix in the
enum (`O_p`, `OH_p`, `CO_p`, `H2O_p`, `C_p`) but the host-facing `name` column is
`O(gr)` … `C(gr)`.

| idx range | family | names (in order) |
|---|---|---|
| 0–11 | H/He/D core | `H` `H2` `e` `H+` `H2+` `H3+` `H-` `He` `He+` `D` `HD` `D+` |
| 12–26 | carbon + oxygen | `C` `CH` `CH2` `C+` `O` `O2` `OH` `CO` `H2O` `HCO` `O+` `OH+` `H2O+` `HCO+` `H3O+` |
| 27–30 | Li / Mg | `Li` `Li+` `Mg` `Mg+` |
| 31–34 | grain charge | `Gr` `Gr+` `Gr-` `Gr2-` |
| 35–39 | ice mantle (physisorbed) | `O(gr)` `OH(gr)` `CO(gr)` `H2O(gr)` `C(gr)` |
| 40–43 | dust reservoirs | `C_s` `O_s` `Mg_s` `Gr_v` |

Like the primordial minimal set, the order is model-specific — do **not** reuse
the 95-species `metal_grain::Sp` indices here; resolve by name (§4) or via the
`metal_grain_minimal::Sp` enum. That enum lives in the in-tree header
`models/metal_grain/minimal.h` (used by the bundled apps), **not** in the
installed public `core/species_index.h`; an installed external consumer resolves
this model's columns by name through the registry (§4 Pattern B, §7.3).

---

## 4. Initialising `y[]` without hard-coded indices

Hard-coding `y[7] = 0.083` is fragile: the index 7 means `He` in
Nakauchi2019/Minimal but the *order is model-specific*, and a comment is the only
thing tying the number to the species. ARCHE gives you two index-free ways to set
the vector explicitly.

### Pattern A — named C++ enums (compile-time, recommended for C++)

`api/arche_api.h` already pulls in `core/species_index.h` (via `core/state.h`),
so the one facade include is enough to index `y[]` by name — no second include
is needed:

```cpp
#include "api/arche_api.h"   // brings in zero_metal::Sp, abundance_ref:: too
using namespace arche;
namespace zm = zero_metal;

ChemStateZM& s = prim_cell_state(*cell);
s.y.fill(0.0);
s.y[zm::H]  = 1.0 - 1.0e-4;         // atomic hydrogen, ~all of nH
s.y[zm::e]  = 1.0e-4;               // electrons
s.y[zm::Hp] = 1.0e-4;              // H+  (charge balances e-)
s.y[zm::He] = abundance_ref::yHe;  // helium, n(He)/nH
s.y[zm::D]  = abundance_ref::yD;   // deuterium
s.y[zm::Li] = abundance_ref::yLi;  // lithium
```

The enum for each model lives in its namespace: `zero_metal::`,
`zero_metal_minimal::`, `metal_grain::`, `metal_grain_minimal::` (table §3).

### Pattern B — name → index map at runtime (C and C++)

The name-registry API reports each model's species names *in `y[]` order*, so the
host builds its own column map without knowing `N_sp` at compile time. This is the
right pattern when the model is chosen from a config file, or for a C / Fortran
host.

```c
/* C — arche_capi.h */
ArcheModelRuntime* m = arche_model_create("Nakauchi2019", "data/primordial.h5");
int nsp = arche_model_n_species(m);                 /* 23 */
const char* const* names = arche_model_species(m);  /* {"H","H2","e","H+",...} */

/* Build name -> index once, then use it for every cell. */
int iH = -1, ie = -1, iHp = -1, iHe = -1;
for (int i = 0; i < nsp; ++i) {
    if      (!strcmp(names[i], "H"))  iH  = i;
    else if (!strcmp(names[i], "e"))  ie  = i;
    else if (!strcmp(names[i], "H+")) iHp = i;
    else if (!strcmp(names[i], "He")) iHe = i;
}

ArcheModelCell* c = arche_model_cell_create(m);
double* y = arche_model_cell_y(c);                  /* length nsp */
for (int i = 0; i < nsp; ++i) y[i] = 0.0;
y[iH]  = 1.0 - 1.0e-4;
y[ie]  = 1.0e-4;
y[iHp] = 1.0e-4;
y[iHe] = 8.33e-2;
arche_model_cell_set_scalars(c, /*nH=*/1.0e4, /*T_K=*/200.0, 1.22, 5.0/3.0);
```

Use whichever name strings your host needs; the catalog labels are listed in §3
and `arche_model_species()` returns them verbatim. `arche_model_registry_name(i)`
/ `arche_model_count()` enumerate the registered model names.

### Recommended initial values

For pristine (un-ionized, mostly atomic) gas, the conventional starting point is:

| species | initial `y` | meaning |
|---|---|---|
| `H`  | `1 − x_e` | atomic H carries (almost) all of nH |
| `e`  | `x_e` | free electron fraction (e.g. `1e-4`) |
| `H+` | `x_e` | matches `e` so net charge ≈ 0 |
| `He` | `abundance_ref::yHe` = `8.33e-2` | helium relative to H |
| `D`  | `abundance_ref::yD` = `2.58e-5` | deuterium |
| `Li` | `abundance_ref::yLi` = `4.65e-10` | lithium |
| all others | `0` | molecules/ions build up during the step |

For the metal network also seed `C`, `O`, `Mg`, `Na`, `K` from
`abundance_ref::yC`, `yO`, `yMg`, `yNa`, `yK` (scaled by your metallicity), and
set `params.Z_metal`. Constants live in `core/species_index.h`
(`arche::abundance_ref::`); the app layer's `arche::abundance::get_*_set(preset)`
exposes named presets (`solar`, `alpha-enhanced`).

> **Always zero the whole vector first** (`fill(0.0)` / a loop), then set the
> species you know. ARCHE does not assume a default composition — an
> uninitialised `y[]` slot is read as a real abundance.

Element and charge conservation is **enforced on all four shipped networks**.
After each internal backward-Euler step `solve/conservation.h` projects `y` onto
the affine set that restores the element totals the step was handed and drives
the net charge to zero. A step whose projection would move any abundance by more
than its error weight (`atol + rtol*|y_i|`) is rejected and retried with a smaller
step, so an accepted step is never silently repaired.

| Network | Tracked rows | Rows |
|---|---:|---|
| `Nakauchi2019` | 5 | H / He / D / Li + charge |
| `Nakauchi2019_Minimal` | 5 | H / He / D / Li + charge |
| `Nakauchi2021_Minimal` | 9 | the above + C / O / Mg + grain number + charge (the compact network carries no K or Na) |
| `Nakauchi2021` | 11 | the above + K / Na + grain number + charge |

Read the count back from the table you actually loaded rather than trusting this
table — `prim_table_n_invariants(*tbl)` and its three siblings (§7.1) report how
many rows that handle is **configured** with. A table built through a path that
dropped the rows returns `0` here; every internal step is then rejected
(`SolveFailure::invariants_unavailable`) and the call fails.

⚠ The configured count is not the count enforced on a given step. Weights are
the species' own abundances, so a row whose carriers are all zero drops out of
that step's solve — correctly, since it owes nothing — while the step still
reports `conservation_projected = true`.

The metals locked in the dust (`C_s`, `O_s`, `Mg_s`, `Na_s`, `K_s`) are species
of `y[]` and count toward their elements, so the metal networks conserve the
**gas + dust** total of each element, including over the grain evaporation that
releases those metals inside the call (`solve/grain_evaporation.h`).

Measured element conservation, from the `y` dataset and its `species` attribute
of the collapse applications (icpx 2025.0 Release, ζ = 1e-17 s⁻¹,
`nH` = 0.1 → 1e23 cm⁻³): the largest relative deviation of an element total from
its initial value is 4.4e-16 for `Nakauchi2019` and 7.1e-15 for `Nakauchi2021` at
Z = 1e-3 and Z = 1 (the Mg, Na and K rows; H, He, D, Li, C and O stay ≤ 8.9e-16).
The net charge stays ≤ 8.4e-16 of the positive charge.

Pass a charge-neutral state (e.g. `y[e] ≈ y[H+] + y[He+] + …`). Over any fluid
cell a plasma can hold a net charge of only (λ_D / L)² of its electrons, with
λ_D = 6.90 cm (T/K)^½ (n_e/cm⁻³)^−½ the Debye length, so a charged input is a
numerical defect of the caller's state, not physics. The projection that
neutralises it shifts every species by at most 1e-2 of itself
(`conservation::kMaxRelShift`). A species of charge q moves by q times the
electrons' share, and multiply charged species count even when they are absent, so
it neutralises a net charge up to (1e-2 / max|q|) · Σ qᵢ² yᵢ: 3.33e-3 of it for
`Nakauchi2019` (Li³⁺) and 1e-2 for `Nakauchi2019_Minimal`. `chem_full_step` treats
a charged input in four ways:

- Neutralisable, and neutralising it would move no abundance by more than its
  error weight: the input is integrated as given, and the steps' own projections
  remove the imbalance.
- Neutralisable, with a larger repair: the input is projected once onto zero net
  charge (element totals kept) before the first internal step.
  `status.input_neutralized` (C: `input_neutralized`) is then `true`.
- Not neutralisable, but |net charge| ≤ `atol`: the net charge is moved onto the
  electrons (`y[e] += Σ qᵢ yᵢ`), which carry no element. The change is at most
  `atol`, inside every error weight, and `input_neutralized` stays `false`.
- Otherwise: `std::invalid_argument` (C: `ARCHE_ERR_INVALID_ARGUMENT`), and the
  cell is not modified.

`status.input_charge_imbalance` reports the input's |Σ qᵢ yᵢ| / Σ |qᵢ| yᵢ in every
case. For example, at `nH` = 1e10 cm⁻³ and T = 750 K a net charge of 1e-6 per H
(y(e⁻) = 2e-4, imbalance 2.5e-3) is neutralised on all four models, while a state
without electrons is rejected. The neutralisation changes ionisation states without
accounting for their ionisation energy; its energy error is at most the shifted
abundance times the ionisation potential, e.g. 1.4e-4 of 1.5 kT at 750 K for this
example, so treat a frequent `input_neutralized` as a defect to fix in the caller.
Element totals are not repaired either: the projection targets
the totals *the step was handed*, not any absolute reference, so it preserves a
composition that was wrong on entry.

---

## 5. Inputs you set each step

### `ChemState<N_sp>` — the per-cell gas state (you fill `y[]` + scalars)

```cpp
template<int N_sp>
struct ChemState {
    std::array<double, N_sp> y{};  // abundances n(i)/nH  (dimensionless)
    double nH    = 0.0;            // H nucleus number density [cm^-3]
    double T_K   = 0.0;            // gas temperature [K]
    double mu    = 1.0;            // mean molecular weight [m_p]
    double gamma = 5.0/3.0;        // adiabatic index
};
```

Aliases: `ChemStateZM` (`N_sp=23`), `ChemStateMG` (`N_sp=95`),
`ChemStatePrimMinimal` (`N_sp=15`), `ChemStateMetalMinimal` (`N_sp=44`). Access it
via `prim_cell_state(cell)` / `metal_cell_state(cell)` /
`prim_minimal_cell_state(cell)` / `metal_minimal_cell_state(cell)` (C++), or
`arche_*_cell_set_state` / `arche_model_cell_y` + `…_set_scalars` (C).

`mu` and `gamma` are **updated by the kernel** (see §6) — pass your current
values in and read the new ones out.

### `ChemParams` — read-only environment (one per thread)

| field | unit | meaning |
|---|---|---|
| `zeta` | s⁻¹ | CR ionization rate. **Ignored by `chem_full_step` — it is overwritten by `shield.zeta`.** |
| `T_rad` | K | radiation (CMB) temperature, e.g. `2.725·(1+z)`. Must be finite > 0. |
| `T_gr_K` | K | grain-temperature warm-start seed (metal only); the *solved* value comes back in `ChemFullRates::T_gr_K`. |
| `Z_metal` | Z⊙ | metallicity (metal only). |
| `rtol`, `atol`, `rtol_T` | — | error-control weights of the integrator: `atol + rtol·|y_i|` per abundance (per H nucleus, the same `atol` for every species) and `rtol_T·T`; defaults `1e-4`, `1e-14`, `1e-4`. |
| `internal_max_steps` | — | attempted internal steps per call (accepted + rejected); default `100000`. |
| `internal_max_step` | s | upper bound on the internal step width; default `+inf` (no bound). |
| `saha_nonthermal_tol`, `use_saha`, `saha_threshold_nH` | — | the Saha regime above `saha_threshold_nH` (default `1e18`), with non-thermal ionisation only where it moves `y_e` at the Saha state by at most `saha_nonthermal_tol` (default `1e-3`). |
| `J_H2`, `J_H2O`, `J_tot`, `H`, `H2`, `He` | — | filled by the kernel; leave 0. |

The C mirror `ArcheChemParams` carries the same fields (`use_saha` as `int`,
1 = true); `arche_chem_params_init()` sets the defaults above.

### `ChemShielding` — pre-computed shielding / radiation environment

You fill this every call; it carries no state.

| field | unit | meaning |
|---|---|---|
| `zeta` | s⁻¹ | effective CR ionization rate, **already attenuated** by your shielding. `chem_full_step` uses it directly. |
| `Nc_H`, `Nc_H2`, `Nc_HD` | cm⁻² | column densities for line self-shielding (optically-thin limit: 0). |
| `cnt_column` | g cm⁻² | column mass `ρ·L_shield` of the continuum shielding length; `0` = optically thin (see below). |
| `J_LW21` | 10⁻²¹ erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹ | Lyman–Werner intensity; drives H₂/HD photodissociation + H⁻ photodetachment. `0` = off (default). |
| `zeta_X`, `E_X_eV` | s⁻¹, eV | X-ray HI ionization rate (pre-attenuated) and representative photon energy (`ARCHE_XRAY` builds only; `0` = off). |
| `Nc_CO`, `Nc_OH`, `Nc_H2O`, `Nc_CII`, `Nc_CI`, `Nc_OI` | cm⁻² | metal-line column densities (metal only; 0 otherwise). |

The continuum optical depth is not an input. `chem_full_step` derives it at every
state it evaluates from the opacity of that state and the column you pass, over a local
shielding length `L_shield` (Jeans length, Sobolev length, …):

```cpp
shield.cnt_column = rho * L_shield;
// inside chem_full_step, at each evaluated state (T, and T_gr for the grain models):
//   tau_cnt = (k_gas(T) + k_gr(T_gr)) * cnt_column,  esc_cnt = tau_cnt > 1 ? 1/tau_cnt^2 : 1
```

For the grain models the grain temperature and `tau_cnt` are solved together. A depth
built from a previous step's opacity lags the state by one step and can drive a
period-2 oscillation of the grain temperature, which is why the depth is not an input.

---

## 6. What ARCHE returns

### Updated in place (read these back into your fluid fields)

| object | fields the kernel writes |
|---|---|
| `cell` state | `y[]` (new abundances), `mu`, `gamma` |
| `cell` history | rate cache, conservation carry, step-size hint and (metal) escape-probability state — keep it between steps of the same cell; `model_cell_reset` before reusing the cell for another one |
| `ChemFullRates::T_gr_K` | solved grain temperature (metal only); feed it back as next step's `params.T_gr_K` |

### `ChemFullRates` — full cooling/heating breakdown (return of `chem_full_step_*`)

All rates are **specific powers in `erg g⁻¹ s⁻¹`**, averaged over the requested
interval (the accepted internal steps, weighted by their widths). Metal-only fields are
0 for the primordial networks.

```cpp
struct ChemFullRates {
    // aggregate
    double Lambda_net;   // net cooling = Lambda_line + Lambda_cnt + Lambda_chem − Gamma_CR (− Gamma_X)
    double Lambda_line;  // total line cooling
    double Lambda_cnt;   // total continuum cooling (grain + gas)
    double Lambda_chem;  // chemistry (endothermic-reaction) cooling
    double Gamma_CR;     // cosmic-ray ionization heating
    double Gamma_X;      // X-ray heating (ARCHE_XRAY only)
    // per-line (all networks)
    double Lambda_H2, Lambda_HD, Lambda_Lya;
    // per-line (metal only)
    double Lambda_CO, Lambda_OH, Lambda_H2O, Lambda_CII, Lambda_CI, Lambda_OI;
    double Lambda_gr;    // grain continuum cooling
    double Lambda_gas;   // gas (H ff + H2 CIA) continuum cooling
    // end state of the interval
    double k_gas;        // gas Planck-mean opacity [cm^2/g]
    double k_gr;         // grain opacity × Z_metal [cm^2/g]  (metal only)
    double tau_cnt;      // continuum optical depth (from shield.cnt_column)
    double esc_cnt;      // continuum escape fraction
    double T_gr_K;       // solved grain temperature [K]  (metal only)
    SolveStatus status;  // internal steps, Newton iterations, failure reason
    bool   solver_failed;// true → step rolled back; state unchanged
    bool   conservation_projected;  // true → every accepted step was projected
                         // onto the element/charge invariants (§4).
};
```

**Sign / energy-equation convention.** `Lambda_net` is the *net cooling* rate
(positive cools). Apply it to the **specific internal energy** `e` `[erg/g]`:

```
de/dt |_chem = − Lambda_net          // e = k_B·T / ((γ−1)·μ·m_p)
```

i.e. add `−Lambda_net·dt` to `e` alongside your adiabatic / compression terms,
then recompute `T` from the returned `mu`, `gamma`. (The bundled collapse driver
does exactly this in `update_thermodynamics`.) The kernel integrates its own
temperature inside the interval but returns `T_K` at its input value: `e − dt·Lambda_net`
is the energy of its end state, so your EOS inversion recovers that temperature.

**On solver failure** `solver_failed == true`: the cell state is left rolled back
(unchanged) and the rates are not meaningful. The kernel has already shortened its
internal steps under its own error control before giving up (`status.reason` says
why). The outcome can still depend on `dt`: `no_progress` stops once refused
attempts shorten the step below `8ε·dt`, a floor that grows with `dt`, so a call that
fails can succeed over a shorter `dt` (e.g. `nH` = 1e22 cm⁻³ with `cnt_column` = 0
fails at `dt` = 1e9 s and succeeds at `dt` ≤ 1e8 s; the free-fall time there is
1.4e4 s). Retrying the interval in shorter calls is therefore a valid recovery.

`status.reason` (C: `failure_reason`) is the outcome of the call and
`status.last_failure` (C: `last_failure`) the cause of the last rejected internal
attempt, which is also set on a successful call that rejected some attempts. The C
API names the numeric values `ARCHE_SOLVE_*` (`api/arche_capi.h`):

| code | `SolveFailure` | meaning |
|---:|---|---|
| 0 | `none` | success |
| 1 | `max_iterations` | Newton iteration limit |
| 2 | `stagnation` | not produced by this version |
| 3 | `nonfinite` | NaN/Inf in a residual or update |
| 4 | `singular_jacobian` | Newton matrix not factorisable |
| 5 | `line_search` | damping found no acceptable step |
| 6 | `temperature_bracket` | no bracket for the temperature root |
| 7 | `temperature_iterations` | temperature root iteration limit |
| 8 | `projection_rejected` | conservation projection declined or exceeded the error weights (§4) |
| 9 | `invariants_unavailable` | the table carries no invariant rows |
| 10 | `output_rejected` | local error above tolerance, or the end state failed validation |
| 11 | `thermal_rejected` | not produced by this version |
| 12 | `interval_budget` | `internal_max_steps` exhausted |
| 13 | `no_progress` | a rejection for a cause other than the local error shortened the step below the interval's resolution `8ε·dt` (ε = 2.2e-16), the step can no longer advance the time, or `dt <= 0` |

The error control may take the step below `8ε·dt`: a stiff start can need its
first steps many decades below `dt` (e.g. ~1e-11 s at n_H = 1e19 cm⁻³ and
4000 K). A rejection for any other cause — a failed solve, a refused projection,
a non-finite error estimate — that shortens the step below `8ε·dt` ends the
call. A cause that does not shrink with the step width — a defect of the state
rather than of the step — therefore ends the call after at most log₄(1/(8ε)) ≈ 25
shortenings, and `last_failure` is the cause of the narrowest refused attempt.
`status.last_failure_run` counts the attempts rejected for it in a row since the
last accepted step, and `status.last_failure_run_h` [s] is the width of the first
of them; `status.rejections_by_reason[k]` counts the rejected attempts of the
call by code `k`. A run that starts at a width independent of `dt` and reaches
`8ε·dt` points at the input state. Example (a charged input with the entry
neutralisation of §4 switched off): every model rejects every width from about
2.9e-3 s down for `dt` = 10, 1e3 and 1e5 s, in 24–30 attempts, with
`last_failure = projection_rejected`. The C fields have the same names
(`ARCHE_SOLVE_COUNT` entries in `rejections_by_reason`).

---

## 7. Functions

Three interchangeable layers select the model differently but step through the
same kernel (numerically identical). The C ABI in `arche_capi.h` mirrors each
with `ArcheChem*` POD structs and `int` status codes (`ARCHE_OK` /
`ARCHE_ERR_*`, message via `arche_last_error()`).

### 7.1 Per-model, model fixed at compile time

```cpp
// Tables (load once, share across threads):
PrimTablePtr         tbl   = load_prim_table_owned("data/primordial.h5");
PrimMinimalTablePtr  tblm  = load_minimal_table_owned("data/primordial.h5");
MetalTablePtr        tblM  = load_metal_table_owned("data/metal_grain.h5");
MetalMinimalTablePtr tblMm = load_metal_minimal_table_owned("data/metal_grain.h5");

// Cells (one per thread):
PrimCellPtr         cell   = prim_cell_create_owned();
PrimMinimalCellPtr  cellm  = prim_minimal_cell_create_owned();
MetalCellPtr        cellM  = metal_cell_create_owned();
MetalMinimalCellPtr cellMm = metal_minimal_cell_create_owned();

// Invariant rows the conservation projection will enforce for each table
// (0 = none; see §4). Read from the handle, not from the factory:
int n   = prim_table_n_invariants(*tbl);           // 5
int nm  = prim_minimal_table_n_invariants(*tblm);  // 5
int nM  = metal_table_n_invariants(*tblM);         // 10
int nMm = metal_minimal_table_n_invariants(*tblMm);// 8

ChemStateZM&           s   = prim_cell_state(*cell);
ChemStatePrimMinimal&  sm  = prim_minimal_cell_state(*cellm);
ChemStateMG&           sM  = metal_cell_state(*cellM);
ChemStateMetalMinimal& sMm = metal_minimal_cell_state(*cellMm);

ChemFullRates r  = chem_full_step_prim          (*cell,   dt, params, shield, *tbl);
ChemFullRates m  = chem_full_step_prim_minimal  (*cellm,  dt, params, shield, *tblm);
ChemFullRates M  = chem_full_step_metal         (*cellM,  dt, params, shield, *tblM);
ChemFullRates Mm = chem_full_step_metal_minimal (*cellMm, dt, params, shield, *tblMm);
```

C ABI (Prim/Metal only): `arche_prim_table_load` / `arche_prim_cell_create` /
`arche_prim_cell_set_state` / `arche_chem_full_step_prim` (+ `_metal`).

### 7.2 Runtime `Model`-enum dispatch (length-erased)

Pick `Model::Prim` or `Model::Metal` at runtime; `y` is exposed as a pointer of
length `cell_n_sp()`.

```cpp
TablePtr tbl  = load_table_owned(Model::Prim, "data/primordial.h5");
CellPtr  cell = cell_create_owned(Model::Prim);
double*  y    = cell_y(*cell);                 // length cell_n_sp(*cell)
cell_nH(*cell) = 1e4; cell_T_K(*cell) = 200; /* … */
ChemFullRates r = chem_full_step(*cell, dt, params, shield, *tbl);
```

C ABI: `arche_table_load` / `arche_cell_create` / `arche_cell_y` /
`arche_cell_set_scalars` / `arche_chem_full_step`.

### 7.3 Name registry (species metadata; includes the minimal model)

Choose by name, read the species count and names, map your columns (§4), step.
This is the only runtime (name-erased) path that exposes the minimal models
(`Nakauchi2019_Minimal`, `Nakauchi2021_Minimal`).

```cpp
ModelRuntimePtr m = model_create_owned("Nakauchi2019_Minimal", "data/primordial.h5");
int n = model_n_species(*m);                   // 15
const char* const* names = model_species(*m);  // y-order names
ModelCellPtr c = model_cell_create_owned(*m);
double* y = model_cell_y(*c);                  // length n
model_cell_nH(*c) = 1e4; model_cell_T_K(*c) = 200; /* mu, gamma … */
ChemFullRates r = model_step(*m, *c, dt, params, shield);
```

C ABI: `arche_model_create` / `arche_model_n_species` / `arche_model_species` /
`arche_model_cell_create` / `arche_model_cell_y` / `arche_model_cell_set_scalars`
/ `arche_model_step`.

### 7.4 Registry cell reset and read-only thermodynamic maps

The name-registry layer additionally exposes these C++ functions:

```cpp
void model_cell_reset(ModelCell& c) noexcept;
double model_mu_from_y(const ModelCell& c) noexcept;
double model_gamma_from_y(const ModelCell& c, double T_K) noexcept;
double model_T_from_e(const ModelCell& c, double e_cgs) noexcept;
```

The C ABI mirrors them:

```c
void arche_model_cell_reset(ArcheModelCell* cell);
double arche_model_mu_from_y(const ArcheModelCell* cell);
double arche_model_gamma_from_y(const ArcheModelCell* cell, double T_K);
double arche_model_T_from_e(const ArcheModelCell* cell, double e_cgs);
```

`model_cell_reset()` clears only hidden integration history: the reaction-rate
cache `var`, the conservation-projection remainder `cons_carry`, the step-size
hint `step_hint`, and, for the two `Nakauchi2021*` models, the persisted
line-cooling `EscapeState` warm-start arrays. It does not change `y[]`, `nH`,
`T_K`, stored `mu`, or stored `gamma`.

`step_hint` is the internal step size that the error control proposed at the end
of the cell's previous call. The next call starts with it, rather than with the
whole interval, when `nH`, `T_K`, `zeta`, `J_LW21` and `cnt_column` have each
changed by at most 10% since that call. The hint sets only where the error
control starts: every accepted step still meets `rtol` / `atol` / `rtol_T`. A
stale hint from another physical cell therefore costs only extra or rejected
steps; it does not change the accuracy contract.
Reset a thread-local scratch cell when assigning it to a different physical
cell. Do not reset between consecutive steps of the same physical cell: that
would discard the intended conservation remainder and solver warm starts.

The three thermo functions are read-only and do not write their results back to
the cell. `mu` is dimensionless in proton-mass units, `gamma` is dimensionless,
`T_K` is in K, and `e_cgs` is in erg g⁻¹. They evaluate the same
internal `thermo` EOS implementation that the kernel uses for the same `y` and
`T`, including its species-summation order and `c_H2()` evaluation. The
`thermo` headers and implementation are private library details: external code
uses only this facade (or the C ABI below) and must not include them. This
guarantee concerns identical inputs; a completed step may subsequently alter
`y[]` through the conservation projection or the operator-split photo channels.

`model_gamma_from_y()` requires a positive, finite `T_K`. It returns quiet NaN
for `T_K <= 0`, NaN, and either infinity. Positive finite inputs retain the
kernel EOS summation order and bit identity described above.

`model_T_from_e()` inverts

```text
e(T) = k_B T [1.5 S_atoms(y) + c_H2(T) y_H2]
       / [(1 + 4 yHe) m_p]
```

by secant iteration in `log(T)` inside `[1, 10¹²] K`, stopping when a step
changes `log(T)` by at most `1e-14`; where the iteration leaves the bracket or
does not converge, bracket-preserving bisection in `log(T)` over the bracket
decides, stopping when the bracket endpoints are adjacent at double precision.
A finite energy below or above the bracket maps to `1 K` or `10¹² K`,
respectively. An `e_cgs = NaN` input produces a NaN result.

`c_H2()` evaluates the para/ortho (1:3) rotational partition functions
continuously across `1000 K`; it does not switch to the classical rotational
limit at that temperature.  Through `10^6 K`, the rotational sum retains
enough states that the last included level has a Boltzmann exponent of at
least 60; above that temperature it uses the continuous rigid-rotor
high-temperature expansion.  Thus the EOS does not have a branch-induced
multivalued energy interval near `1000 K`, and `model_T_from_e()` returns the
unique root, to the rounding level of `e(T)`.  The residual
step in `c_H2()` at `1000 K`, left by the change in the number of summed
rotational states, is `2.2e-14` relative.  Measured `T -> e -> T` round trips on
H₂-dominated compositions (`y_H2` = 0.5 and 0.1, which maximise the `c_H2`
contribution) stay within `1.3e-11` relative over the whole `[1, 10¹²] K`
bracket — the worst case sits at the cold end, near `1.3 K` — tightening to
`1.6e-12` for `T ≥ 10 K` and `1.5e-12` for `T ≥ 100 K`.

For C, reset is a no-op when `cell == NULL`; each NULL thermo query returns a
quiet NaN.

### 7.5 Ohmic and ambipolar magnetic diffusivities

```cpp
struct MagneticDiffusivities {
  double eta_O;  // Ohmic      [cm^2 s^-1]
  double eta_A;  // ambipolar  [cm^2 s^-1]
};
MagneticDiffusivities model_resistivity(const ModelCell& c, double T_K,
                                        double B_G) noexcept;
```

```c
typedef struct { double eta_O; double eta_A; } ArcheMagneticDiffusivities;
ArcheMagneticDiffusivities arche_model_resistivity(const ArcheModelCell* cell,
                                                   double T_K, double B_G);
```

`model_resistivity()` evaluates η_O and η_A for the cell's `y[]` and `nH` at
the gas temperature `T_K` [K] and field strength `B_G` [G]. It reads the cell
and changes nothing.

- **Charged particles:** every charged species of the model; the grain charge
  states are spread over the MRN size distribution.
- **Collision partners:** H, H₂ and He.

Both values are quiet NaN unless `nH` and `T_K` are positive and finite and
`B_G` is non-negative and finite. Without charged particles both are `+inf`;
`eta_A` is `0` for `B_G = 0`. The C function returns NaN for a NULL cell.

The formulas, the momentum-transfer rate coefficients and their sources, the
grain size bins, and the differences from Nakauchi, Omukai & Susa (2021) are in
[`magnetic_resistivity.md`](magnetic_resistivity.md).

---

## 8. What `chem_full_step_*` does internally

It integrates the interval `dt` with the error-controlled backward Euler of
`solve/integrator.h`:

1. Copies `shield.zeta` → `params.zeta` (the kernel uses the attenuated rate).
2. Takes internal steps. The first spans the whole interval, or the cell's
   `step_hint` (the step size the error control proposed at the end of the previous
   call) where it is shorter and `nH`, `T_K`, `zeta`, `J_LW21` and `cnt_column` have
   each changed by at most 10% since that call. Each attempt solves one step of width
   `h` and two of `h/2` from the same start (step doubling); the difference, in the
   weighted norm of `rtol`/`atol`/`rtol_T`, is the local error, and `h` is shortened
   or grown from it. The tolerances bound this local error (an RMS over the species
   and T, so one component may exceed its weight), not the accuracy of the call's
   result: the steps are first order, so the error accumulated over a call falls
   only as about `rtol`^½. Over one call of `dt` = 1e9 s at `nH` = 1e14 cm⁻³
   (`Nakauchi2019`), `rtol` = `rtol_T` = 1e-4 leaves abundances off by up to 1.2e-2
   of themselves against `rtol` = 1e-9. In the collapse applications the error from
   the outer time step dominates instead: tightening `rtol` and `rtol_T` from 1e-4
   to 1e-7 moves T by ≤ 8.4e-5 dex (`Nakauchi2019` and `Nakauchi2021_Minimal` at
   Z = 1e-3, ζ = 1e-17 s⁻¹).
3. Each step solves the abundances and the specific internal energy in **one** Newton
   system: `y − y0 − h·r(y, T) = 0` and
   `e(y, T) − e0 + h·(Lambda_chem + (Lambda_rad(start) + Lambda_rad(end))/2 − Gamma_CR) = 0`,
   with the rate coefficients, line cooling (`line_cool`, `line_cool_metal`),
   continuum cooling with its depth from `cnt_column` (`cnt_cool` /
   `cnt_cool_metal_consistent`, which also solves `T_gr`) and CR heating all evaluated
   at the iterate. Where that Newton fails, the step is solved with the species at a
   fixed trial temperature and the energy equation as a bracketed root in `T`.
   With `shield.J_LW21 > 0` the Lyman–Werner channels are first-order sinks in `r` and
   its Jacobian, evaluated at the iterate:
   H₂ + hν(LW) → 2H  (k = 1.38×10⁻¹² · J₂₁ · f_sh(H₂); Wolcott-Green & Haiman 2019
   eqs. 7–8, the exponent α held at its 8000 K value above 8000 K),
   HD + hν(LW) → H + D  (same prefactor, WG2011 self-shielding form),
   H⁻ + hν → H + e⁻  (k = 1.10×10⁻¹⁰ · J₂₁; Tegmark et al. 1997; no shielding).
   They act at fixed specific internal energy: their part of the composition change
   does not enter `Lambda_chem`.
4. Above `saha_threshold_nH` (subject to `saha_nonthermal_tol` with non-thermal
   ionisation) the composition is the model's Saha state at `T` and the step solves the
   same energy equation for `T` alone.
5. Each accepted step is projected onto the element and charge invariants (refused,
   and the step shortened, if the repair exceeds the tolerance weights), then gets the
   operator-split channels:
   - **Lyman–Werner** after a Saha step (if `shield.J_LW21 > 0`) — the same three
     channels as first-order exponential decay over the step width;
   - **X-ray secondary ionization / heating** (if `shield.zeta_X > 0`, `ARCHE_XRAY`).
6. Returns the width-weighted means of the step rates,
   `Lambda_net = Lambda_line + Lambda_cnt + Lambda_chem − Gamma_CR (− Gamma_X)`,
   with `k_gas`, `k_gr`, `tau_cnt`, `esc_cnt` and `T_gr_K` of the end state.

---

## 9. Thread safety

| object | rule |
|---|---|
| cell / `ChemState` | **per-thread** — never share a cell between threads |
| `ChemParams` | **per-thread** (`T_gr_K` warm-start is per-cell) |
| table (`PrimTable` / `MetalTable` / `ModelRuntime`) | **read-only, freely shared** across threads |
| `ChemShielding` | local / per-call, no state |

Cells are not large but are not trivially small either (the metal cell plus the
solver's stack scratch is tens of KB). For OpenMP, allocate one cell per thread
on the heap (the `*_create_owned()` helpers) and reuse it across cells/steps,
calling the stepping function in the parallel region.

When a name-registry cell is a per-thread scratch object, reset it once per new
physical cell before loading that cell's visible state:

```cpp
#pragma omp parallel
{
  ModelCellPtr scratch = model_cell_create_owned(*model);
  #pragma omp for
  for (std::size_t i = 0; i < grid.size(); ++i) {
    model_cell_reset(*scratch);  // grid[i] is unrelated to the prior i
    load_visible_state(grid[i], *scratch);
    for (int substep = 0; substep < n_substeps; ++substep) {
      // No reset here: these are consecutive steps of the same physical cell.
      ChemFullRates r = model_step(*model, *scratch, dt, params, shield);
      /* update the same physical cell */
    }
    store_visible_state(*scratch, grid[i]);
  }
}
```

---

## 10. Fluid-code coupling — end-to-end sketch (C++)

```cpp
#include "api/arche_api.h"   // single include: also provides zero_metal::Sp, abundance_ref::
using namespace arche;
namespace zm = zero_metal;

// ── startup (once) ──────────────────────────────────────────────────────────
PrimTablePtr tbl = load_prim_table_owned("data/primordial.h5");

// ── per thread ──────────────────────────────────────────────────────────────
PrimCellPtr cell = prim_cell_create_owned();
ChemParams  params; params.T_rad = 2.725 * (1.0 + z);

for (int step = 0; step < n_steps; ++step) {
  for (Cell3D& g : grid) {                    // your hydro cells
    ChemStateZM& s = prim_cell_state(*cell);

    // 1) copy fluid → ARCHE
    s.nH = g.nH; s.T_K = g.T; s.mu = g.mu; s.gamma = g.gamma;
    for (int i = 0; i < zm::N_sp; ++i) s.y[i] = g.y[i];   // your stored vector

    // 2) shielding: attenuated CR rate, columns, continuum column
    ChemShielding shield;
    shield.zeta       = g.zeta_attenuated;
    shield.cnt_column = g.rho * g.L_shield;

    // 3) advance one dt (the kernel controls its own internal steps)
    ChemFullRates r = chem_full_step_prim(*cell, g.dt, params, shield, *tbl);
    if (r.solver_failed) { /* cell rolled back: keep its old state for this dt */ continue; }

    // 4) copy ARCHE → fluid
    for (int i = 0; i < zm::N_sp; ++i) g.y[i] = s.y[i];
    g.mu = s.mu; g.gamma = s.gamma;
    g.e += -r.Lambda_net * g.dt;          // specific internal energy update
  }
}
```

For the metal network swap the `prim_*` calls for `metal_*`, use `ChemStateMG` /
`metal_grain::`, set `params.Z_metal` and `params.T_gr_K`, fill the metal column
densities, and feed `r.T_gr_K` back into `params.T_gr_K` each step.

---

## 11. Physical constants (`arche::phys`, via `core/state.h`)

| constant | value | unit |
|---|---|---|
| `k_B` | 1.380662×10⁻¹⁶ | erg K⁻¹ |
| `h_P` | 6.626176×10⁻²⁷ | erg s |
| `m_p` | 1.67262×10⁻²⁴ | g |
| `m_e` | 9.19941×10⁻²⁸ | g |
| `G` | 6.6720×10⁻⁸ | cm³ g⁻¹ s⁻² |
| `sigma_B` | 5.67×10⁻⁵ | erg cm⁻² s⁻¹ K⁻⁴ |
| `eV_to_erg` | 1.6022×10⁻¹² | erg eV⁻¹ |

## Network / abundance constants (`core/species_index.h`)

| constant | value | meaning |
|---|---|---|
| `zero_metal::N_sp` / `N_react` | 23 / 140 | full primordial |
| `zero_metal_minimal::N_sp` / `N_react` | 15 / 33 | compact primordial |
| `metal_grain::N_sp` / `N_react` | 95 / 1200 | metal + grain |
| `metal_grain_minimal::N_sp` / `N_react` | 44 / 113 | compact metal + grain |
| `abundance_ref::yHe` | 8.33×10⁻² | He/H |
| `abundance_ref::yD` | 2.58×10⁻⁵ | D/H |
| `abundance_ref::yLi` | 4.65×10⁻¹⁰ | Li/H |
| `abundance_ref::yC` | 2.69×10⁻⁴ | C/H (solar) |
| `abundance_ref::yO` | 4.90×10⁻⁴ | O/H (solar) |
| `abundance_ref::yMg` | 3.98×10⁻⁵ | Mg/H (solar) |
| `abundance_ref::yNa` | 1.74×10⁻⁶ | Na/H (solar) |
| `abundance_ref::yK` | 1.07×10⁻⁷ | K/H (solar) |

---

## Note for in-tree work

The public facade above wraps an internal template kernel
(`chem_full_step<Model>(ChemCell<N_sp,N_react>&, …)` in `solve/chemistry.h`,
plus `make_*_table()` loaders). Those headers pull in Eigen/HDF5 and are **not
installed** — they are for building the bundled apps inside this repository, not
for external consumers. External couplings use only the installed public headers
listed at the top of this page.
