# HDF5 Output Format

Both `prim_collapse` and `metal_collapse` write a single HDF5 file per run.
Each file contains 1-D datasets of length N_rows (one row per recorded step)
plus a few root-level attributes. The compact-network binaries
(`arche_collapse_prim_minimal`, `arche_collapse_metal_minimal`) write the **same**
schema with a reduced `y` column count and a `_min` filename tag (see below).

---

## Schema version

Every file carries a root attribute **`schema_version`** (int).

| Version | Dataset naming |
|---|---|
| **2** *(current)* | Fortran `x`-prefix dropped and abbreviations expanded: `nH`, `Lambda_*` (cooling), `Gamma_*` (heating), `M_J` (Jeans mass), `lambda_J` (Jeans length). |
| **1** *(legacy)* | x-prefixed names: `xnH`, `xLmbd_*`, `xGam_*`, `xMJ`, `xlmbd_J`. Files written before this rename have **no** `schema_version` attribute — treat a missing attribute as version 1. |

The bundled tools (`analyze_collapse.py`, `resample_collapse.py`) read **both**
versions: legacy v1 dataset names are aliased to their v2 names on load
(see `LEGACY_TO_V2`), so the rest of each tool uses the v2 names uniformly.
Only the dataset **names** changed at v1→v2; the stored **values are identical**.

v1→v2 dataset name map: `xnH`→`nH`, `xLmbd_net`→`Lambda_net`,
`xLmbd_line`→`Lambda_line`, `xLmbd_cnt`→`Lambda_cnt`, `xLmbd_ch`→`Lambda_chem`,
`xLmbd_gas`→`Lambda_gas`, `xLmbd_gr`→`Lambda_gr`, `xLmbd_Lya`→`Lambda_Lya`,
`xLmbd_H2`→`Lambda_H2`, `xLmbd_HD`→`Lambda_HD`, `xLmbd_CO`→`Lambda_CO`,
`xLmbd_OH`→`Lambda_OH`, `xLmbd_H2O`→`Lambda_H2O`, `xLmbd_CII`→`Lambda_CII`,
`xLmbd_CI`→`Lambda_CI`, `xLmbd_OI`→`Lambda_OI`, `xGam_CR`→`Gamma_CR`,
`xGam_cmp`→`Gamma_cmp`, `xMJ`→`M_J`, `xlmbd_J`→`lambda_J`.

---

## File naming

### `prim_collapse` (full primordial)

```
collapse_CR<cr_tag>[_fret<fret_tag>][_z<zred_tag>].h5
```

| Component | Value |
|---|---|
| `<cr_tag>` | `PRIM_ZETA0` with `.` replaced by `p` (e.g. `1e-17`, `1p5e-18`, `0`) |
| `_fret<fret_tag>` | Appended only when `PRIM_FF_RET` ≠ 1.0 (scalar) or when `PRIM_FRET_TABLE` is set (always `_fret-step`) |
| `_z<zred_tag>` | Appended only when `PRIM_REDSHIFT` ≠ 0.0 |

Examples:
- Scalar f_ret=3.0 → `collapse_CR1e-17_fret3p0.h5`
- Step-function table → `collapse_CR1e-17_fret-step.h5`

Default output directory: `results/prim/h5/`

### `arche_collapse_prim_minimal` (minimal primordial)

```
collapse_CR<cr_tag>_min[_fret<fret_tag>][_z<zred_tag>].h5
```

This compact model writes the same scalar datasets/attributes as full
primordial, but uses a `_min` file stem and stores a 15-column `y` dataset.

When executed via `run_collapse.sh --minimal`, the wrapper also creates a
compatibility symlink without `_min` so existing plot/resample tools can keep
the same filename pattern.

### `metal_collapse`

```
collapse_CR<cr_tag>_Z<z_tag>[_fret<fret_tag>][_JLW<jlw_tag>][_z<zred_tag>].h5
```

Default output directory: `results/metal/h5/`

### `arche_collapse_metal_minimal` (minimal metal-grain)

```
collapse_CR<cr_tag>_Z<z_tag>[_fret<fret_tag>][_JLW<jlw_tag>][_z<zred_tag>]_min.h5
```

This compact model writes the same scalar datasets/attributes as full
metal-grain, but appends a `_min` tag to the file stem and stores a 40-column
`y` dataset (vs. 89). The `network` attribute reads
`"metal_grain_minimal N_sp=40 N_react=113 (...)"`.

---

## `prim_collapse` — zero-metal network

### Root attributes

| Attribute | Type | Description |
|---|---|---|
| `description` | string | Human-readable label |
| `schema_version` | int | Output schema version (currently `2`; see [Schema version](#schema-version)) |
| `cr_tag` | string | CR tag derived from `PRIM_ZETA0` |
| `zeta0_cgs` | float64 | CR ionization rate [s⁻¹] |
| `f_ret` | float64 | Free-fall retardation factor (1.0 = standard).  When a step-function table is used this is the **initial** (first row) value |
| `f_ret_table` | string | *(optional)* Path to the f_ret step-function table file.  Present only when `PRIM_FRET_TABLE` was set |
| `zred` | float64 | Cosmological redshift z (0.0 if not set) |
| `T_rad` | float64 | CMB radiation temperature [K] = 2.725 × (1+z) |
| `J_LW21` | float64 | Lyman-Werner intensity J₂₁ [10⁻²¹ erg/s/cm²/Hz/sr].  `0.0` = no LW field |
| `zeta_X` | float64 | X-ray primary HI photoionization rate ζ_X [s⁻¹].  `0.0` = no X-ray field.  Present only when compiled with `-DARCHE_XRAY=ON` |
| `E_X_eV` | float64 | Representative X-ray photon energy E_X [eV].  Present only when compiled with `-DARCHE_XRAY=ON` |
| `ic_T_K0` | float64 | Initial gas temperature [K] used for this run (default 100.0) |
| `ic_y_e0` | float64 | Initial electron fraction y(e⁻) = y(H⁺) (default 1e-4) |
| `ic_y_H2` | float64 | Initial H₂ fraction (default 6e-7) |
| `ic_y_HD` | float64 | Initial HD fraction (default 4e-10) |
| `network` | string | Full: `"zero_metal N_sp=23 N_react=140"`; minimal: `"zero_metal_minimal N_sp=15 N_react=33 (...)"` |
| `units_density` | string | `"cm^-3 (nH), g/cm^3 (rho)"` |
| `units_cooling` | string | `"erg g^-1 s^-1"` |
| `units_time` | string | `"s"` |
| `units_length` | string | `"cm"` |
| `units_B` | string | `"G"` |
| `exit_code` | int | Why the run stopped: `0` normal (nH ceiling), `1` max_iter, `2` T ceiling (10⁵ K), `3` internal energy e ≤ 0, `4` NaN/Inf in nH or T, `5` NR subcycle max depth exceeded |
| `exit_message` | string | The same reason in words, e.g. `"Normal: nH reached nH_stop"` |
| `conservation_rows` | int | Element/charge rows the projection is **configured** to enforce (`n_invariants`, read from the loaded table).  `0` = this network registers none.  ⚠ Not the number enforced on any given step — see below |
| `conservation_steps_total` | int | Steps whose chemistry solve succeeded.  ⚠ Counted right after the kernel, so a step whose later thermodynamic update fails (`exit_code` = 3) is still included; `0` means the run solved nothing |
| `conservation_steps_projected` | int | Of those, the steps whose abundances the projection was applied to |
| `conservation_first_unprojected_step` | int | Step index of the first unprojected step; `-1` = every step was projected |
| `conservation_first_unprojected_nH` | float64 | Density [cm⁻³] that same step entered with; `-1.0` = every step was projected |


**Reading the `conservation_*` attributes.** They report whether the element /
charge projection (`solve/conservation.h`) actually ran, which the abundances
themselves do not show: a projection that never fires leaves a file that looks
completely normal.  The metal network spent a development cycle in exactly that
state — `conservation_rows` at its expected 10, most steps unprojected — because
two of the ten rows were empty and an empty row took the whole weighted matrix
down with it.

Read the counts together with `conservation_rows`:

| `conservation_rows` | projected / total | Meaning |
|---|---|---|
| `0` | `0 / N` | This network registers no invariant rows, so the projection declines on every step. Expected, not a fault |
| `> 0` | `N / N` | Every step was projected. Expected |
| `> 0` | anything below `N / N` | ⚠ The projection declined on a step it was configured for. Start at `conservation_first_unprojected_nH` |

The last row is unambiguous in a collapse run: a step that fails to converge
ends the run with `exit_code` = 5 instead of being counted here, so the only
remaining causes are the projection's own decline paths — a row that owes a
residual no surviving species can carry, a weighted matrix that is not positive
definite, or a repair larger than the allowed relative shift.

⚠ **`conservation_rows` counts configured rows, not rows enforced on a step.**
The projection weights each species by its own abundance, so a row whose carrier
species are all zero contributes nothing and is dropped from that step's solve.
Dropping it is correct — there is nothing to conserve and nothing owed — and the
step still counts as projected. But it means `conservation_rows = 10` together
with `projected = total` does **not** assert that ten rows were enforced on every
step.

A measured example: in a metal-grain collapse at Z = 1 with the app holding all
potassium and sodium in its grain reservoir, the K and Na rows carry nothing
until the grains evaporate at n<sub>H</sub> ≈ 1.7 × 10¹⁵ cm⁻³ — **68.9 % of the
run** — so eight of ten rows are solved over that stretch while the attributes
read `rows = 10`, `projected = total`. The remaining eight elements are fully
enforced throughout; no element is silently unconserved. Read the pair as
"the projection ran on every step", not "ten invariants held on every step".

These attributes do not change `schema_version`, which versions **dataset
names** only (see [Schema version](#schema-version)): no dataset name, shape, or
dataset value changes when they are added. The attributes are themselves new
stored values — the guarantee is about the datasets.

### Datasets

All datasets are 1-D with shape `(N_rows,)` and dtype `float64`, except `step`
(int32) and `y` (float64; full: `(N_rows, 23)`, minimal: `(N_rows, 15)`).

One row is written every `PRIM_OUTPUT_STRIDE` steps (default: 10).

| Dataset | Shape | Description |
|---|---|---|
| `step` | (N_rows,) int32 | Integration step number |
| `y` | (N_rows, 23) full / (N_rows, 15) minimal | Species abundances [dimensionless, / nH].  Attribute `species` lists names (see below) |
| `nH` | (N_rows,) | H number density [cm⁻³] |
| `T_K` | (N_rows,) | Gas temperature [K] |
| `rho` | (N_rows,) | Mass density [g cm⁻³] |
| `Lambda_net` | (N_rows,) | Net cooling Λ_line + Λ_cnt + Λ_ch − Γ_CR [erg g⁻¹ s⁻¹] |
| `Lambda_line` | (N_rows,) | Total line cooling (H₂ + HD + Ly-α) |
| `Lambda_cnt` | (N_rows,) | Continuum cooling (dust + H ff + H₂ CIA) |
| `Lambda_chem` | (N_rows,) | Chemical (endothermic) cooling |
| `Lambda_gas` | (N_rows,) | Gas continuum subset (H ff + H₂ CIA) |
| `Lambda_Lya` | (N_rows,) | Lyman-alpha cooling |
| `Lambda_H2` | (N_rows,) | H₂ line cooling |
| `Lambda_HD` | (N_rows,) | HD line cooling |
| `Gamma_CR` | (N_rows,) | CR ionization heating |
| `Gamma_cmp` | (N_rows,) | Compressional heating p/ρ/t_eff |
| `t_ff` | (N_rows,) | True free-fall time [s]  (= t_eff / f_ret) |
| `t_cool` | (N_rows,) | Cooling time e/\|Λ_net\| [s] |
| `t_chem` | (N_rows,) | Chemistry timescale min_i(y_i/\|Δy_i/Δt\|) [s] |
| `tau_cnt` | (N_rows,) | Continuum optical depth |
| `lambda_J` | (N_rows,) | Jeans length [cm] |
| `M_J` | (N_rows,) | Jeans mass [g] |
| `B_cr` | (N_rows,) | Critical (ambipolar) magnetic field [G] |
| `y_plus` | (N_rows,) | Total positive charge fraction |
| `y_minus` | (N_rows,) | Total negative charge fraction |
| `charge_imbal` | (N_rows,) | \|y⁺ − y⁻\| / (y⁺ + y⁻) |

### Species order (`y` dataset, zero-metal network)

For full primordial, `y` stores all 23 species below.
For minimal primordial, `y` stores 15 compact-network species in this order:
H, H2, e-, H+, H2+, H3+, H-, He, He+, D, HD, D+, Li, LiH, Li+.

`y.attrs["species"]` contains a comma-separated list that matches the
actual columns of `y`:

- Full (23 species):

```
H, H2, e-, H+, H2+, H3+, H-, He, He+, He++, HeH+,
D, HD, D+, HD+, D-, Li, LiH, Li+, Li-, LiH+, Li++, Li+++
```

- Minimal (15 species):

```
H, H2, e-, H+, H2+, H3+, H-, He, He+, D, HD, D+, Li, LiH, Li+
```

Index 0 → H, index 1 → H₂, index 2 → e⁻, …

---

## `metal_collapse` — metal-grain network

### Root attributes

| Attribute | Type | Description |
|---|---|---|
| `schema_version` | int | Output schema version (currently `2`; see [Schema version](#schema-version)) |
| `f_ret` | float64 | Free-fall retardation factor.  When a step-function table is used this is the **initial** (first row) value |
| `f_ret_table` | string | *(optional)* Path to the f_ret step-function table file.  Present only when `METAL_FRET_TABLE` was set |
| `zred` | float64 | Cosmological redshift z (0.0 if not set) |
| `T_rad` | float64 | CMB radiation temperature [K] |
| `J_LW21` | float64 | Lyman-Werner intensity J₂₁ [10⁻²¹ erg/s/cm²/Hz/sr].  `0.0` = no LW field |
| `zeta_X` | float64 | X-ray primary HI photoionization rate ζ_X [s⁻¹].  `0.0` = no X-ray field.  Present only when compiled with `-DARCHE_XRAY=ON` |
| `E_X_eV` | float64 | Representative X-ray photon energy E_X [eV] used for Beer-Lambert shielding cross-section.  Present only when compiled with `-DARCHE_XRAY=ON` |
| `ic_T_K0` | float64 | Initial gas temperature [K] used for this run (default 100.0) |
| `ic_y_e0` | float64 | Initial electron fraction y(e⁻) = y(H⁺) (default 1e-4) |
| `ic_y_H2` | float64 | Initial H₂ fraction (default 6e-7) |
| `ic_y_HD` | float64 | Initial HD fraction (default 4e-10) |
| `network` | string | Full: `"metal_grain N_sp=89 N_react=1200"`; minimal: `"metal_grain_minimal N_sp=40 N_react=113 (...)"` |
| `exit_code` | int | Why the run stopped: `0` normal (nH ceiling), `1` max_iter, `2` T ceiling (10⁵ K), `3` internal energy e ≤ 0, `4` NaN/Inf in nH or T, `5` NR subcycle max depth exceeded |
| `exit_message` | string | The same reason in words, e.g. `"Normal: nH reached nH_stop"` |
| `conservation_rows` | int | Element/charge rows the projection is **configured** to enforce (`n_invariants`, read from the loaded table).  `0` = this network registers none.  ⚠ Not the number enforced on any given step — see below |
| `conservation_steps_total` | int | Steps whose chemistry solve succeeded.  ⚠ Counted right after the kernel, so a step whose later thermodynamic update fails (`exit_code` = 3) is still included; `0` means the run solved nothing |
| `conservation_steps_projected` | int | Of those, the steps whose abundances the projection was applied to |
| `conservation_first_unprojected_step` | int | Step index of the first unprojected step; `-1` = every step was projected |
| `conservation_first_unprojected_nH` | float64 | Density [cm⁻³] that same step entered with; `-1.0` = every step was projected |

See [`prim_collapse` root attributes](#root-attributes) for how to read the
`conservation_*` group.

### Datasets

One row is written every `METAL_OUTPUT_STRIDE` steps (default: 10).

| Dataset | Shape | Description |
|---|---|---|
| `step` | (N_rows,) int32 | Integration step number |
| `y` | (N_rows, 89) full / (N_rows, 40) minimal | Species abundances [/ nH].  Attribute `species` lists names (see below) |
| `nH` | (N_rows,) | H number density [cm⁻³] |
| `T_K` | (N_rows,) | Gas temperature [K] |
| `T_gr_K` | (N_rows,) | Grain temperature [K] |
| `rho` | (N_rows,) | Mass density [g cm⁻³] |
| `Lambda_net` | (N_rows,) | Net cooling [erg g⁻¹ s⁻¹] |
| `Lambda_line` | (N_rows,) | Total line cooling |
| `Lambda_cnt` | (N_rows,) | Total continuum cooling (grain + gas) |
| `Lambda_gr` | (N_rows,) | Grain continuum cooling |
| `Lambda_gas` | (N_rows,) | Gas continuum cooling (H ff + H₂ CIA) |
| `Lambda_chem` | (N_rows,) | Chemical cooling |
| `Gamma_cmp` | (N_rows,) | Compressional heating |
| `Lambda_H2` | (N_rows,) | H₂ line cooling |
| `Lambda_HD` | (N_rows,) | HD line cooling |
| `Lambda_Lya` | (N_rows,) | Lyman-alpha cooling |
| `Lambda_CO` | (N_rows,) | CO line cooling |
| `Lambda_OH` | (N_rows,) | OH line cooling |
| `Lambda_H2O` | (N_rows,) | H₂O line cooling |
| `Lambda_CII` | (N_rows,) | C II line cooling |
| `Lambda_CI` | (N_rows,) | C I line cooling |
| `Lambda_OI` | (N_rows,) | O I line cooling |
| `Gamma_CR` | (N_rows,) | CR ionization heating |
| `t_ff` | (N_rows,) | True free-fall time [s] |
| `t_cool` | (N_rows,) | Cooling time [s] |
| `t_chem` | (N_rows,) | Chemistry timescale [s] |
| `tau_cnt` | (N_rows,) | Continuum optical depth |
| `lambda_J` | (N_rows,) | Jeans length [cm] |
| `M_J` | (N_rows,) | Jeans mass [g] |
| `B_cr` | (N_rows,) | Critical magnetic field [G] |
| `y_plus` | (N_rows,) | Total positive charge fraction |
| `y_minus` | (N_rows,) | Total negative charge fraction |
| `charge_imbal` | (N_rows,) | \|y⁺ − y⁻\| / (y⁺ + y⁻) |

### Species order (`y` dataset, metal-grain network)

`y.attrs["species"]` contains a comma-separated list matching the actual columns
of `y`. For the full 89-species order see
[`api_reference.md` §3.3](api_reference.md). The minimal model stores 40
compact-network species in this order:

```
H, H2, e-, H+, H2+, H3+, H-, He, He+, D, HD, D+,
C, CH, CH2, C+, O, O2, OH, CO, H2O, HCO, O+, OH+, H2O+, HCO+, H3O+,
Li, Li+, Mg, Mg+, Gr, Gr+, Gr-, Gr2-,
O(gr), OH(gr), CO(gr), H2O(gr), C(gr)
```

Index 0 → H, index 1 → H₂, index 2 → e⁻, …  (see [`api_reference.md`
§3.4](api_reference.md) for the matching `metal_grain_minimal::Sp` enum.)

---

---

## Resampled CSV output (`resample_collapse.py`)

`tools/resample_collapse.py` reads an HDF5 file and resamples physical quantities
onto a uniform log₁₀(nH) grid, writing a CSV table.
See [`tools/resample_collapse.md`](../tools/resample_collapse.md) for the full CLI reference.

### File naming

```
resample_collapse_CR<cr_tag>[_fret<fret_tag>][_z<zred_tag>].csv          # prim
resample_collapse_CR<cr_tag>_Z<z_tag>[_fret<fret_tag>][_z<zred_tag>].csv # metal
```

The file name mirrors the input HDF5 stem with a `resample_` prefix.

### CSV format

```
# resample_collapse.py output
# species: H,H2,e-,...
# NaN = empty bin (no data points)
log10_nH,T_K,Lambda_net,y_0,y_1,...,y_N
-1.750000,2.726000e+01,1.234567e-20,...
...
```

Lines beginning with `#` are header comments.  The fourth `#` line is the
column-name row (no leading `#`).

### Output columns

| Column | Description | Averaging |
|---|---|---|
| `log10_nH` | Bin centre in log₁₀(nH [cm⁻³]) | — (grid index) |
| `T_K` | Gas temperature [K] | log-average |
| `T_gr_K` | Grain temperature [K] (metal only) | log-average |
| `Lambda_net` | Net cooling rate [erg g⁻¹ s⁻¹] | arithmetic mean |
| `y_0` … `y_N` | Species abundances (prim: N=22, metal: N=88) | log-average |

Empty bins (no data points in range) are written as `NaN`.

### Averaging conventions

| Method | Formula | Used for |
|---|---|---|
| `log_avg` | 10^(mean(log₁₀(x))) for x > 0; NaN if all ≤ 0 | T_K, T_gr_K, species |
| `lin_avg` | mean(x) | Lambda_net (signed/zero-capable) |

The mapping is defined in the `AVG_FUNCS` dict in the script and can be
customised by editing that dict.

---

## Reading HDF5 output with Python

```python
import h5py
import numpy as np

with h5py.File("results/prim/h5/collapse_CR1e-17.h5", "r") as f:
    # Scalar attributes
    zeta0 = f.attrs["zeta0_cgs"]
    f_ret = f.attrs["f_ret"]

    # 1-D arrays
    nH   = f["nH"][:]
    T_K   = f["T_K"][:]
    Lnet  = f["Lambda_net"][:]

    # Species abundances: shape (N_rows, 23)
    y         = f["y"][:]
    species   = f["y"].attrs["species"].split(",")
    idx_H2    = species.index("H2")
    y_H2      = y[:, idx_H2]
```
