# API Reference — `chem_full_step()` and chemistry module

Header: `src/include/chemistry.h`
Namespace: `chemistry`

---

## Overview

The chemistry module exposes two entry points:

| Function | Use case |
|---|---|
| `chem_full_step()` | **Recommended.** Advances chemistry *and* computes all cooling/heating rates in one call. |
| `chem_step()` | Low-level. Advances chemistry only; returns chemistry cooling + CR heating. Caller handles line and continuum cooling separately. |

---

## Data types

### `ChemState<N_sp>`

Per-cell thermodynamic and chemical state.  Updated in-place by the solver.

```cpp
template<int N_sp>
struct ChemState {
    std::array<double, N_sp> y{};  // species abundances [dimensionless, / xnH]
    double xnH   = 0.0;            // H number density [cm^-3]
    double T_K   = 0.0;            // gas temperature [K]
    double xmu   = 1.0;            // mean molecular weight [m_p]
    double gamma = 5.0/3.0;        // adiabatic index
};
```

Convenience aliases: `ChemStateZM` (N_sp=23), `ChemStateMG` (N_sp=89).

---

### `ChemCell<N_sp, N_react>`

Bundles per-cell state with the inter-step reaction rate cache used by the
predictor-corrector solver.  Each thread / computational cell must own its own
`ChemCell`.

```cpp
template<int N_sp, int N_react>
struct ChemCell {
    ChemState<N_sp>               state{};
    std::array<double, 2*N_react> var{};   // inter-step rate cache; zero-init
    void reset_var() noexcept;             // clear var[] (call when cell is moved)
};
```

The metal-grain specialisation adds an `EscapeState es{}` member that stores
escape-probability fractions for the warm-start Newton-Raphson solver in
`line_cool_metal()`.

**Convenience aliases**

| Alias | N_sp | N_react | Extra member |
|---|---|---|---|
| `ZeroMetalCell` | 23 | 140 | — |
| `MetalGrainCell` | 89 | 1200 | `EscapeState es{}` |

> **Stack vs heap**: `MetalGrainCell` is ~20 KB.  Inside `chem_step` /
> `chem_full_step` the solver also allocates ~126 KB on the stack.
> For OpenMP use per-thread heap allocation (`std::make_unique`).

---

### `ChemParams`

External parameters that are **read-only** for the chemistry kernel, except
`T_gr_K` which is updated in-place (metal-grain only).

```cpp
struct ChemParams {
    double zeta    = 0.0;    // CR ionization rate [s^-1]
                             // NOTE: ignored by chem_full_step(); use shield.zeta
    double T_rad   = 2.725;  // CMB radiation temperature [K]
    double T_gr_K  = 0.0;    // grain temperature [K]   (metal_grain only; updated in-place)
    double Z_metal = 0.0;    // metallicity [Z_sun]     (metal_grain only)
    // (additional internal fields used by the kernel)
};
```

---

### `ChemShielding`

Pre-computed shielding environment.  Must be filled by the caller before each
call to `chem_full_step()`.

```cpp
struct ChemShielding {
    double zeta;            // effective CR ionization rate [s^-1]  (pre-attenuated)
    double xNc_H   = 0.0;  // H  column density [cm^-2]
    double xNc_H2  = 0.0;  // H2 column density [cm^-2]
    double xNc_HD  = 0.0;  // HD column density [cm^-2]
    double tau_cnt = 0.0;  // continuum optical depth
    double esc_cnt = 1.0;  // continuum escape fraction  (1 = optically thin)
    double J_LW21  = 0.0;  // Lyman-Werner intensity [10^-21 erg/s/cm^2/Hz/sr]
                            // Drives H2/HD photodissociation and H- photodetachment
                            // (operator-split step inside chem_full_step).
                            // 0.0 = no LW field (default).
    // metal_grain only (leave 0 for zero_metal network):
    double xNc_CO  = 0.0;
    double xNc_OH  = 0.0;
    double xNc_H2O = 0.0;
    double xNc_CII = 0.0;
    double xNc_CI  = 0.0;
    double xNc_OI  = 0.0;
};
```

**`shield.zeta`** must carry the full attenuation already applied.
The chemistry kernel uses it directly without further reduction.

**`shield.J_LW21`** is the mean specific intensity of the Lyman-Werner (11.2–13.6 eV)
radiation field in units of J₂₁ = 10⁻²¹ erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹.
Set to 0.0 (the default) to disable LW photodissociation entirely.
Column densities `xNc_H2` and `xNc_HD` are used for self-shielding corrections;
they may be 0.0 (optically thin limit).

In a 3-D fluid code, `tau_cnt` and `esc_cnt` are computed from the opacities
returned by the *previous* step:

```cpp
tau_cnt = (xk_gr + xk_gas) * rho * L_shield;
esc_cnt = (tau_cnt > 1.0) ? 1.0 / (tau_cnt * tau_cnt) : 1.0;
```

---

### `ChemFullRates`

All cooling and heating rate components returned by `chem_full_step()`.
All rates are in **erg g⁻¹ s⁻¹** (specific power).
Metal-grain-only fields are 0 for the zero-metal network.

```cpp
struct ChemFullRates {
    // Aggregate
    double xLmbd_net  = 0.0;  // net cooling = line + cnt + ch − CR
    double xLmbd_line = 0.0;  // total line cooling
    double xLmbd_cnt  = 0.0;  // total continuum cooling (grain + gas)
    double xLmbd_ch   = 0.0;  // chemistry (endothermic reaction) cooling
    double xGam_CR    = 0.0;  // CR ionization heating
    // Per-line (both networks)
    double xLmbd_H2   = 0.0;  // H2 ro-vib line cooling
    double xLmbd_HD   = 0.0;  // HD line cooling
    double xLmbd_Lya  = 0.0;  // Lyman-alpha line cooling
    // Per-line (metal_grain only)
    double xLmbd_CO   = 0.0;
    double xLmbd_OH   = 0.0;
    double xLmbd_H2O  = 0.0;
    double xLmbd_CII  = 0.0;
    double xLmbd_CI   = 0.0;
    double xLmbd_OI   = 0.0;
    double xLmbd_gr   = 0.0;  // grain continuum cooling
    double xLmbd_gas  = 0.0;  // gas (H ff + H2 CIA) continuum cooling
    // Opacities for next step's tau_cnt
    double xk_gas     = 0.0;  // gas opacity [cm^2/g]
    double xk_gr      = 0.0;  // grain opacity × Z_metal [cm^2/g]  (metal only)
};
```

---

### `ChemRates`

Simplified rates returned by `chem_step()`.

```cpp
struct ChemRates {
    double xLmbdch = 0.0;  // chemistry cooling [erg g^-1 s^-1]
    double xGam_CR = 0.0;  // CR heating        [erg g^-1 s^-1]
};
```

---

## Functions

### `chem_full_step()`

```cpp
template<int N_sp, int N_react>
ChemFullRates chem_full_step(ChemCell<N_sp, N_react>& cell,
                              double dt,
                              ChemParams& params,
                              const ChemShielding& shield,
                              const ReactionTable<N_sp, N_react>& tbl);
```

**What it does (in order)**

1. Copies `shield.zeta` into `params.zeta` (so the kernel sees the attenuated rate).
2. `line_cool()` — H₂, HD, Ly-α line cooling.
3. `line_cool_metal()` — CO, OH, H₂O, C II, C I, O I line cooling (metal only).
4. `cnt_cool()` / `cnt_cool_metal()` — continuum cooling; updates `params.T_gr_K` and
   returns `xk_gas`, `xk_gr`.
5. `chemcool()` — advances the chemical network by `dt`; updates `cell.state.y`,
   `xmu`, `gamma`; returns `xLmbd_ch`.
6. **LW operator split** (if `shield.J_LW21 > 0`) — applies Lyman-Werner photodissociation
   as first-order exponential decay after the chemistry solver:
   - H₂ + hν(LW) → H + H  (k = 1.38×10⁻¹² J₂₁ f_sh(H₂) s⁻¹;  self-shielding: Wolcott-Green & Haiman 2019 eq. 7–8)
   - HD + hν(LW) → H + D  (k = 1.38×10⁻¹² J₂₁ f_sh(HD) s⁻¹;  WG2011 functional form)
   - H⁻ + hν    → H + e⁻  (k = 1.10×10⁻¹⁰ J₂₁ s⁻¹;  Tegmark et al. 1997; no self-shielding)

   Species indices (both networks): `y[0]`=H, `y[1]`=H₂, `y[2]`=e⁻, `y[6]`=H⁻, `y[11]`=D, `y[12]`=HD.
7. Sums CR heating from the rate cache `cell.var[]`.
8. Computes `xLmbd_net`.

**Updated in-place**

| Object | Fields |
|---|---|
| `cell.state` | `y`, `xmu`, `gamma` |
| `cell.var` | rate cache (warm-start for next step) |
| `params.T_gr_K` | grain temperature (metal only) |

**Return value**: `ChemFullRates` — see above.

**`params.zeta` note**: `chem_full_step()` **ignores** `params.zeta` on entry and
overwrites it with `shield.zeta`.  Always set `shield.zeta`, not `params.zeta`.

---

### `chem_step()`

```cpp
template<int N_sp, int N_react>
ChemRates chem_step(ChemCell<N_sp, N_react>& cell,
                    double dt,
                    const ChemParams& params,
                    const ReactionTable<N_sp, N_react>& tbl);
```

Low-level entry point: runs only `chemcool()` and CR heating.
Line and continuum cooling must be handled by the caller.
Use `chem_full_step()` for new integrations.

**Updated in-place**: `cell.state.y`, `cell.state.xmu`, `cell.state.gamma`,
`cell.var`.

---

### `make_zero_metal_table()`

```cpp
inline ZeroMetalTable make_zero_metal_table(const std::string& data_dir);
```

Loads the zero-metal (23-species, 140-reaction) reaction tables from `data_dir`.

Required files in `data_dir`:

| File | Content |
|---|---|
| `react_prm.dat` | Reaction network (140 entries) |
| `mass.dat` | Species masses [g] |
| `react_prm_saha.dat` | Saha equilibrium entries |
| `pf_H3p.dat` | Partition function for H₃⁺ |
| `pf_HD.dat` | Partition function for HD |
| `pf_HDp.dat` | Partition function for HD⁺ |
| `pf_LiH.dat` | Partition function for LiH |
| `pf_LiHp.dat` | Partition function for LiH⁺ |

---

### `make_metal_grain_table()`

```cpp
inline MetalGrainTable make_metal_grain_table(const std::string& data_dir,
                                               const std::string& metal_data_dir);
```

Loads the metal-grain (89-species, 1200-reaction) reaction tables from two directories.

`data_dir` (shared primordial data, same as zero-metal `data_dir`):

| File | Content |
|---|---|
| `mass_metal.dat` | Species masses [g] (89 species) |
| `pf_H3p.dat` | H₃⁺ |
| `pf_HD.dat` | HD |
| `pf_HDp.dat` | HD⁺ |

`metal_data_dir` (metal-grain-specific data):

| File | Content |
|---|---|
| `react_metal_grain.dat` | Reaction network (1200 entries) |
| `react_grain_surface.dat` | Grain surface reactions |
| `react_metal_saha.dat` | Saha equilibrium entries |
| `pf_CH3.dat` | CH₃ |
| `pf_CH4.dat` | CH₄ |
| `pf_HO2.dat` | HO₂ |
| `pf_CO2.dat` | CO₂ |
| `pf_H2CO.dat` | H₂CO |
| `pf_H2O2.dat` | H₂O₂ |
| `pf_LiHp.dat` | LiH⁺ |

---

## Thread safety

| Object | Thread safety |
|---|---|
| `ChemCell` / `ChemState` | **Per-thread only** — not safe to share |
| `ChemParams` | **Per-thread** (`T_gr_K` is updated in-place) |
| `ZeroMetalTable` / `MetalGrainTable` | **Read-only, fully shared** across threads |
| `ChemShielding` | Local / per-call, no state |

---

## Minimal usage example (zero-metal)

```cpp
#include "chemistry.h"

// ── Startup (once per simulation) ────────────────────────────────────────────
chemistry::ZeroMetalTable tbl = chemistry::make_zero_metal_table("/path/to/data");

// ── Per-cell allocation ───────────────────────────────────────────────────────
auto cell = std::make_unique<chemistry::ZeroMetalCell>();
auto& y   = cell->state.y;   // species abundances [/ n_H]
// ... initialise y[0..22], cell->state.xnH, T_K, xmu, gamma ...

chemistry::ChemParams params;
params.T_rad = 2.725 * (1.0 + z);   // CMB temperature [K]

// ── Per time step ─────────────────────────────────────────────────────────────
chemistry::ChemShielding shield;
shield.zeta    = zeta_attenuated;    // pre-attenuated CR rate [s^-1]
shield.xNc_H2  = col_H2;            // H2 column density [cm^-2]
shield.xNc_HD  = col_HD;
shield.tau_cnt = tau_cnt;
shield.esc_cnt = (tau_cnt > 1.0) ? 1.0 / (tau_cnt * tau_cnt) : 1.0;

cell->state.xnH = nH;  cell->state.T_K = T;
cell->state.xmu = mu;  cell->state.gamma = gamma;

const auto rates = chemistry::chem_full_step(*cell, dt, params, shield, tbl);

// Update fluid state
mu     = cell->state.xmu;
gamma  = cell->state.gamma;
xk_gas = rates.xk_gas;           // opacity for next tau_cnt
// rates.xLmbd_net → net cooling [erg g^-1 s^-1] for energy equation
```

For the metal-grain network, replace `ZeroMetalCell` → `MetalGrainCell`,
`make_zero_metal_table()` → `make_metal_grain_table(data_dir, metal_data_dir)`,
set `params.Z_metal` and `params.T_gr_K`, and fill the additional column density
fields (`xNc_CO`, `xNc_OH`, `xNc_H2O`, `xNc_CII`, `xNc_CI`, `xNc_OI`) in
`ChemShielding`.

---

## Physical constants (`chemistry::phys`)

Accessible without including any extra header (exposed via `state.h`).

| Constant | Value | Unit |
|---|---|---|
| `xk_B` | 1.380662 × 10⁻¹⁶ | erg K⁻¹ |
| `h_P` | 6.626176 × 10⁻²⁷ | erg s |
| `xm_p` | 1.67262 × 10⁻²⁴ | g |
| `pi` | 3.14159265358979 | — |
| `G` | 6.6720 × 10⁻⁸ | cm³ g⁻¹ s⁻² |
| `sigma_B` | 5.67 × 10⁻⁵ | erg cm⁻² s⁻¹ K⁻⁴ |

---

## Network constants

| Constant | Zero-metal | Metal-grain |
|---|---|---|
| `zero_metal::N_sp` | 23 | — |
| `zero_metal::N_react` | 140 | — |
| `metal_grain::N_sp` | — | 89 |
| `metal_grain::N_react` | — | 1200 |
| `metal_grain::yHe` | — | 8.33 × 10⁻² |
| `metal_grain::yD` | — | 2.58 × 10⁻⁵ |
| `metal_grain::yLi` | — | 4.65 × 10⁻¹⁰ |
| `metal_grain::yC` | — | 2.69 × 10⁻⁴ |
| `metal_grain::yO` | — | 4.90 × 10⁻⁴ |
| `metal_grain::yMg` | — | 3.98 × 10⁻⁵ |
