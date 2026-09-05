# API Reference — libarche public API

Public headers (installed under `include/arche/`):

| Header | Language | Entry points |
|---|---|---|
| `api/arche_api.h` | C++ (`namespace arche`) | RAII handles, `chem_full_step_*`, runtime `Model` / name-registry dispatch. Errors via C++ exceptions. |
| `api/arche_capi.h` | C ABI (`extern "C"`) | mirror C structs, `arche_*` functions, status codes + `arche_last_error()`. For C / Fortran / `ctypes`. |
| `core/state.h` | C++ | `ChemState`, `ChemParams`, `ChemShielding`, `ChemFullRates`, `ChemRates`, `phys::` constants. |
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
| Metal + grain (Nakauchi 2021) | `"Nakauchi2021"` | `Model::Metal` | 89 | 1200 | `data/metal_grain.h5` |
| Metal minimal | `"Nakauchi2021_Minimal"` | — | **40** | 113 | `data/metal_grain.h5` |

- **Nakauchi2019** — full primordial H/He/D/Li network (23 species).
- **Nakauchi2019_Minimal** — compact 15-species subset of the primordial
  network (H/D/Li chemistry + cosmic-ray ionization), for embedded use where the
  full network is too heavy. Loads the *same* `primordial.h5` and remaps it
  internally. He⁺ is tracked only under cosmic rays (no collisional He
  ionization), so without a CR field He⁺ stays ≈ 0.
- **Nakauchi2021** — metal + dust-grain network (89 species), adds C/O/Mg/Na/K
  chemistry, grain charge states, and ice-mantle surface species.
- **Nakauchi2021_Minimal** — compact 40-species composed subset of the
  metal+grain network: the H/He/D core, a reduced C/O molecular network, Li, Mg,
  four grain charge states, and five ice-mantle (physisorbed) species. The
  113 reactions split into 103 gas-phase (68 standard + 8 cosmic-ray + 27
  ion–grain charge transfer) and 10 grain-surface freeze-out. A first-class
  composed model — not a runtime mask — derived from the full network's keep-set;
  it loads the *same* `metal_grain.h5` (reusing its partition functions and Saha
  table) and solves a 42×42 system. It drops the species the reduced network
  omits (He⁺⁺, the higher C/O ions, Na, K, Gr²⁺).

The two minimal models are reachable from C++ via the `prim_minimal_*` /
`metal_minimal_*` entry points and from both C++ and C via the **name registry**
(§7.3). The fixed-enum C ABI (`arche_*_prim` / `arche_*_metal`) covers only the
23- and 89-species networks.

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

### 3.3 Nakauchi2021 — 89 species (`metal_grain::Sp`)

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

(`_p` = physisorbed on grain, `_c` = chemisorbed; the gas-phase parent is named
without the suffix.)

### 3.4 Nakauchi2021_Minimal — 40 species (`metal_grain_minimal::Sp`)

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

Like the primordial minimal set, the order is model-specific — do **not** reuse
the 89-species `metal_grain::Sp` indices here; resolve by name (§4) or via the
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
After each implicit chemistry step `solve/conservation.h` projects `y` back onto
the affine set that restores the element totals the step was handed and drives
the net charge to zero.

| Network | Tracked rows | Elements |
|---|---:|---|
| `Nakauchi2019` | 5 | H / He / D / Li + charge |
| `Nakauchi2019_Minimal` | 5 | H / He / D / Li + charge |
| `Nakauchi2021_Minimal` | 8 | the above + C / O / Mg + charge (the compact network carries no K or Na) |
| `Nakauchi2021` | 10 | the above + K / Na + charge |

Read the count back from the table you actually loaded rather than trusting this
table — `prim_table_n_invariants(*tbl)` and its three siblings (§7.1) report how
many rows that handle is **configured** with. A table built through a path that
dropped the rows returns `0` here, and then `conservation_projected` is `false`
on every step with no other outward sign.

⚠ The configured count is not the count enforced on a given step. Weights are
the species' own abundances, so a row whose carriers are all zero drops out of
that step's solve — correctly, since it owes nothing — while the step still
reports `conservation_projected = true`.

Primordial: measured residuals stay at the float64 floor of each element total
(≤ 4.4e-16 relative; ≤ 1.5e-13 when `J_LW21 > 0`), against ≤ 7.6e-11 without it.

Metal (`Nakauchi2021`, Z = 1, `nH_stop = 1e23`): the residual is held near
1.3e-12 from `nH` = 1e14 to 1e22. Without the projection it grows by roughly two
orders per density decade over that range, reaching 2.3e-3 by 1e22.

> **Metal networks gained the projection in 2026-08; their numerical output
> changed as a result.** Element conservation was previously not enforced there
> and `conservation_projected` was `false` on every metal step. If you pinned
> results against an older ARCHE, re-baseline them.

⚠ The grain-surface reservoirs the `collapse_metal_grain` application keeps
(depleted C / O / Mg / K / Na) live outside `y[]`, so the projection conserves
the **gas-phase** total of those elements across a chemistry step. That is the
correct invariant for the step, because the application moves mass between the
grains and the gas only between steps, never inside one.

In both cases charge neutrality and element totals should hold in your initial
condition (e.g. `y[e] ≈ y[H+] + y[He+] + …`). Note the projection targets the
totals *the step was handed*, not any absolute reference, so it preserves — it
does not repair — a composition that was wrong on entry.

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

Aliases: `ChemStateZM` (`N_sp=23`), `ChemStateMG` (`N_sp=89`),
`ChemStatePrimMinimal` (`N_sp=15`), `ChemStateMetalMinimal` (`N_sp=40`). Access it
via `prim_cell_state(cell)` / `metal_cell_state(cell)` /
`prim_minimal_cell_state(cell)` / `metal_minimal_cell_state(cell)` (C++), or
`arche_*_cell_set_state` / `arche_model_cell_y` + `…_set_scalars` (C).

`mu` and `gamma` are **updated by the kernel** (see §6) — pass your current
values in and read the new ones out.

### `ChemParams` — read-only environment (one per thread)

| field | unit | meaning |
|---|---|---|
| `zeta` | s⁻¹ | CR ionization rate. **Ignored by `chem_full_step` — it is overwritten by `shield.zeta`.** Used only by the low-level `chem_step`. |
| `T_rad` | K | radiation (CMB) temperature, e.g. `2.725·(1+z)`. Must be finite > 0. |
| `T_gr_K` | K | grain-temperature warm-start seed (metal only); the *solved* value comes back in `ChemFullRates::T_gr_K`. |
| `Z_metal` | Z⊙ | metallicity (metal only). |
| `J_H2`, `J_H2O`, `J_tot`, `H`, `H2`, `He` | — | filled by the kernel; leave 0. |

### `ChemShielding` — pre-computed shielding / radiation environment

You fill this every call; it carries no state.

| field | unit | meaning |
|---|---|---|
| `zeta` | s⁻¹ | effective CR ionization rate, **already attenuated** by your shielding. `chem_full_step` uses it directly. |
| `Nc_H`, `Nc_H2`, `Nc_HD` | cm⁻² | column densities for line self-shielding (optically-thin limit: 0). |
| `tau_cnt` | — | continuum optical depth (see formula below). |
| `esc_cnt` | — | continuum escape fraction; `1.0` = optically thin. |
| `J_LW21` | 10⁻²¹ erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹ | Lyman–Werner intensity; drives H₂/HD photodissociation + H⁻ photodetachment. `0` = off (default). |
| `zeta_X`, `E_X_eV` | s⁻¹, eV | X-ray HI ionization rate (pre-attenuated) and representative photon energy (`ARCHE_XRAY` builds only; `0` = off). |
| `Nc_CO`, `Nc_OH`, `Nc_H2O`, `Nc_CII`, `Nc_CI`, `Nc_OI` | cm⁻² | metal-line column densities (metal only; 0 otherwise). |

In a 3-D fluid code the continuum optical depth is built from the opacity the
*previous* step returned, over a local shielding length `L_shield` (Jeans length,
Sobolev length, …):

```cpp
tau_cnt = (k_gr + k_gas) * rho * L_shield;       // k_* from last step's ChemFullRates
esc_cnt = (tau_cnt > 1.0) ? 1.0 / (tau_cnt * tau_cnt) : 1.0;
```

---

## 6. What ARCHE returns

### Updated in place (read these back into your fluid fields)

| object | fields the kernel writes |
|---|---|
| `cell` state | `y[]` (new abundances), `mu`, `gamma` |
| `cell` cache | inter-step rate cache (warm-start for the next step — do not clear between steps for the same cell) |
| `ChemFullRates::T_gr_K` | solved grain temperature (metal only); feed it back as next step's `params.T_gr_K` |

### `ChemFullRates` — full cooling/heating breakdown (return of `chem_full_step_*`)

All rates are **specific powers in `erg g⁻¹ s⁻¹`**. Metal-only fields are 0 for
the primordial networks.

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
    // opacities — feed into next step's tau_cnt
    double k_gas;        // gas opacity   [cm^2/g]
    double k_gr;         // grain opacity × Z_metal [cm^2/g]  (metal only)
    double T_gr_K;       // solved grain temperature [K]  (metal only)
    bool   solver_failed;// true → step rolled back; state unchanged
    bool   conservation_projected;  // true → element/charge projection ran.
                         // All four shipped networks register rows (§4), so on
                         // a converged step false means the projection itself
                         // declined; see core/state.h for the four causes.
};
```

**Sign / energy-equation convention.** `Lambda_net` is the *net cooling* rate
(positive cools). Apply it to the **specific internal energy** `e` `[erg/g]`:

```
de/dt |_chem = − Lambda_net          // e = k_B·T / ((γ−1)·μ·m_p)
```

i.e. add `−Lambda_net·dt` to `e` alongside your adiabatic / compression terms,
then recompute `T` from the returned `mu`, `gamma`. (The bundled collapse driver
does exactly this in `update_thermodynamics`.)

**On solver failure** `solver_failed == true`: the cell state is left rolled back
(unchanged) and the rates are not meaningful — shrink `dt` and retry, or carry the
old state.

### `ChemRates` — return of the low-level `chem_step_*`

```cpp
struct ChemRates {
    double Lambda_chem;  // chemistry cooling   [erg g^-1 s^-1]
    double Gamma_CR;     // CR heating          [erg g^-1 s^-1]
    double Gamma_X;      // X-ray heating       [erg g^-1 s^-1]
    bool   conservation_projected;  // see ChemFullRates
};
```

`chem_step_*` advances the chemistry only; the caller is responsible for line and
continuum cooling. Prefer `chem_full_step_*` for new couplings.

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
// low-level chemistry-only variants:
//   chem_step_prim / _prim_minimal / _metal / _metal_minimal
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
cache `var`, the conservation-projection remainder `cons_carry`, and, for the
two `Nakauchi2021*` models, the persisted line-cooling `EscapeState` warm-start
arrays. It does not change `y[]`, `nH`, `T_K`, stored `mu`, or stored `gamma`.
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
`y[]` through the conservation, LW, or X-ray operator-split paths.

`model_gamma_from_y()` requires a positive, finite `T_K`. It returns quiet NaN
for `T_K <= 0`, NaN, and either infinity. Positive finite inputs retain the
kernel EOS summation order and bit identity described above.

`model_T_from_e()` inverts

```text
e(T) = k_B T [1.5 S_atoms(y) + c_H2(T) y_H2]
       / [(1 + 4 yHe) m_p]
```

by bracket-preserving bisection in `log(T)` over `[1, 10¹²] K`, stopping
when the bracket endpoints are adjacent at double precision. A finite energy
below or above the bracket maps to `1 K` or `10¹² K`, respectively. An
`e_cgs = NaN` input produces a NaN result.

`c_H2()` evaluates the para/ortho (1:3) rotational partition functions
continuously across `1000 K`; it does not switch to the classical rotational
limit at that temperature.  Through `10^6 K`, the rotational sum retains
enough states that the last included level has a Boltzmann exponent of at
least 60; above that temperature it uses the continuous rigid-rotor
high-temperature expansion.  Thus the EOS does not have a branch-induced
multivalued energy interval near `1000 K`, and `model_T_from_e()` returns the
unique root selected by its deterministic log-space bisection.  The residual
step in `c_H2()` at `1000 K`, left by the change in the number of summed
rotational states, is `2.2e-14` relative.  Measured `T -> e -> T` round trips on
H₂-dominated compositions (`y_H2` = 0.5 and 0.1, which maximise the `c_H2`
contribution) stay within `1.1e-11` relative over the whole `[1, 10¹²] K`
bracket — the worst case sits at the cold end, near `1.4 K` — tightening to
`1.5e-12` for `T ≥ 10 K` and `6e-13` for `T ≥ 100 K`.

For C, reset is a no-op when `cell == NULL`; each NULL thermo query returns a
quiet NaN.

---

## 8. What `chem_full_step_*` does internally

In order:

1. Copies `shield.zeta` → `params.zeta` (the kernel uses the attenuated rate).
2. `line_cool()` — H₂, HD, Ly-α line cooling.
3. `line_cool_metal()` — CO, OH, H₂O, C II, C I, O I line cooling (metal only).
4. `cnt_cool()` / `cnt_cool_metal()` — continuum cooling; updates `T_gr_K`,
   returns `k_gas`, `k_gr`.
5. `chemcool()` — advances the network by `dt`; updates `y`, `mu`, `gamma`;
   returns `Lambda_chem`. On non-convergence it sub-cycles, then sets
   `solver_failed` and rolls the cell back.
6. **Lyman–Werner operator split** (if `shield.J_LW21 > 0`) — first-order
   exponential decay applied after the solver:
   - H₂ + hν(LW) → 2H  (k = 1.38×10⁻¹² · J₂₁ · f_sh(H₂); Wolcott-Green & Haiman 2019)
   - HD + hν(LW) → H + D  (same prefactor, WG2011 self-shielding form)
   - H⁻ + hν → H + e⁻  (k = 1.10×10⁻¹⁰ · J₂₁; Tegmark et al. 1997; no shielding)
7. **X-ray secondary ionization / heating** (if `shield.zeta_X > 0`, `ARCHE_XRAY`).
8. Sums CR heating from the rate cache → `Gamma_CR`.
9. `Lambda_net = Lambda_line + Lambda_cnt + Lambda_chem − Gamma_CR (− Gamma_X)`.

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

double k_gas_prev = 0.0;   // opacity carried between steps for tau_cnt

for (int step = 0; step < n_steps; ++step) {
  for (Cell3D& g : grid) {                    // your hydro cells
    ChemStateZM& s = prim_cell_state(*cell);

    // 1) copy fluid → ARCHE
    s.nH = g.nH; s.T_K = g.T; s.mu = g.mu; s.gamma = g.gamma;
    for (int i = 0; i < zm::N_sp; ++i) s.y[i] = g.y[i];   // your stored vector

    // 2) shielding from the previous step's opacity
    ChemShielding shield;
    shield.zeta    = g.zeta_attenuated;
    shield.tau_cnt = k_gas_prev * g.rho * g.L_shield;
    shield.esc_cnt = (shield.tau_cnt > 1.0) ? 1.0/(shield.tau_cnt*shield.tau_cnt) : 1.0;

    // 3) advance one dt
    ChemFullRates r = chem_full_step_prim(*cell, g.dt, params, shield, *tbl);
    if (r.solver_failed) { /* shrink dt and retry */ continue; }

    // 4) copy ARCHE → fluid
    for (int i = 0; i < zm::N_sp; ++i) g.y[i] = s.y[i];
    g.mu = s.mu; g.gamma = s.gamma;
    g.e += -r.Lambda_net * g.dt;          // specific internal energy update
    k_gas_prev = r.k_gas;
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
| `metal_grain::N_sp` / `N_react` | 89 / 1200 | metal + grain |
| `metal_grain_minimal::N_sp` / `N_react` | 40 / 113 | compact metal + grain |
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
