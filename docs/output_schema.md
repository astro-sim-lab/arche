# HDF5 Output Format

Both `prim_collapse` and `metal_collapse` write a single HDF5 file per run.
Each file contains 1-D datasets of length N_rows (one row per recorded step)
plus a few root-level attributes.

---

## File naming

### `prim_collapse`

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

### `metal_collapse`

```
collapse_CR<cr_tag>_Z<z_tag>[_fret<fret_tag>][_z<zred_tag>].h5
```

Default output directory: `results/metal/h5/`

---

## `prim_collapse` — zero-metal network

### Root attributes

| Attribute | Type | Description |
|---|---|---|
| `description` | string | Human-readable label |
| `cr_tag` | string | CR tag derived from `PRIM_ZETA0` |
| `zeta0_cgs` | float64 | CR ionization rate [s⁻¹] |
| `f_ret` | float64 | Free-fall retardation factor (1.0 = standard).  When a step-function table is used this is the **initial** (first row) value |
| `f_ret_table` | string | *(optional)* Path to the f_ret step-function table file.  Present only when `PRIM_FRET_TABLE` was set |
| `zred` | float64 | Cosmological redshift z (0.0 if not set) |
| `T_rad` | float64 | CMB radiation temperature [K] = 2.725 × (1+z) |
| `J_LW21` | float64 | Lyman-Werner intensity J₂₁ [10⁻²¹ erg/s/cm²/Hz/sr].  `0.0` = no LW field |
| `ic_T_K0` | float64 | Initial gas temperature [K] used for this run (default 100.0) |
| `ic_y_e0` | float64 | Initial electron fraction y(e⁻) = y(H⁺) (default 1e-4) |
| `ic_y_H2` | float64 | Initial H₂ fraction (default 6e-7) |
| `ic_y_HD` | float64 | Initial HD fraction (default 4e-10) |
| `network` | string | `"zero_metal N_sp=23 N_react=140"` |
| `units_density` | string | `"cm^-3 (xnH), g/cm^3 (rho)"` |
| `units_cooling` | string | `"erg g^-1 s^-1"` |
| `units_time` | string | `"s"` |
| `units_length` | string | `"cm"` |
| `units_B` | string | `"G"` |

### Datasets

All datasets are 1-D with shape `(N_rows,)` and dtype `float64`, except `step`
(int32) and `y` (float64, shape `(N_rows, 23)`).

One row is written every `PRIM_OUTPUT_STRIDE` steps (default: 100).

| Dataset | Shape | Description |
|---|---|---|
| `step` | (N_rows,) int32 | Integration step number |
| `y` | (N_rows, 23) | Species abundances [dimensionless, / xnH].  Attribute `species` lists names (see below) |
| `xnH` | (N_rows,) | H number density [cm⁻³] |
| `T_K` | (N_rows,) | Gas temperature [K] |
| `rho` | (N_rows,) | Mass density [g cm⁻³] |
| `xLmbd_net` | (N_rows,) | Net cooling Λ_line + Λ_cnt + Λ_ch − Γ_CR [erg g⁻¹ s⁻¹] |
| `xLmbd_line` | (N_rows,) | Total line cooling (H₂ + HD + Ly-α) |
| `xLmbd_cnt` | (N_rows,) | Continuum cooling (dust + H ff + H₂ CIA) |
| `xLmbd_ch` | (N_rows,) | Chemical (endothermic) cooling |
| `xLmbd_gas` | (N_rows,) | Gas continuum subset (H ff + H₂ CIA) |
| `xLmbd_Lya` | (N_rows,) | Lyman-alpha cooling |
| `xLmbd_H2` | (N_rows,) | H₂ line cooling |
| `xLmbd_HD` | (N_rows,) | HD line cooling |
| `xGam_CR` | (N_rows,) | CR ionization heating |
| `xGam_cmp` | (N_rows,) | Compressional heating p/ρ/t_eff |
| `t_ff` | (N_rows,) | True free-fall time [s]  (= t_eff / f_ret) |
| `t_cool` | (N_rows,) | Cooling time e/\|Λ_net\| [s] |
| `t_chem` | (N_rows,) | Chemistry timescale min_i(y_i/\|Δy_i/Δt\|) [s] |
| `tau_cnt` | (N_rows,) | Continuum optical depth |
| `xlmbd_J` | (N_rows,) | Jeans length [cm] |
| `xMJ` | (N_rows,) | Jeans mass [g] |
| `B_cr` | (N_rows,) | Critical (ambipolar) magnetic field [G] |
| `y_plus` | (N_rows,) | Total positive charge fraction |
| `y_minus` | (N_rows,) | Total negative charge fraction |
| `charge_imbal` | (N_rows,) | \|y⁺ − y⁻\| / (y⁺ + y⁻) |

### Species order (`y` dataset, zero-metal network)

The `y.attrs["species"]` attribute contains a comma-separated list:

```
H, H2, e-, H+, H2+, H3+, H-, He, He+, He++, HeH+,
D, HD, D+, HD+, D-, Li, LiH, Li+, Li-, LiH+, Li++, Li+++
```

Index 0 → H, index 1 → H₂, index 2 → e⁻, …

---

## `metal_collapse` — metal-grain network

### Root attributes

| Attribute | Type | Description |
|---|---|---|
| `f_ret` | float64 | Free-fall retardation factor.  When a step-function table is used this is the **initial** (first row) value |
| `f_ret_table` | string | *(optional)* Path to the f_ret step-function table file.  Present only when `METAL_FRET_TABLE` was set |
| `zred` | float64 | Cosmological redshift z (0.0 if not set) |
| `T_rad` | float64 | CMB radiation temperature [K] |
| `J_LW21` | float64 | Lyman-Werner intensity J₂₁ [10⁻²¹ erg/s/cm²/Hz/sr].  `0.0` = no LW field |
| `ic_T_K0` | float64 | Initial gas temperature [K] used for this run (default 100.0) |
| `ic_y_e0` | float64 | Initial electron fraction y(e⁻) = y(H⁺) (default 1e-4) |
| `ic_y_H2` | float64 | Initial H₂ fraction (default 6e-7) |
| `ic_y_HD` | float64 | Initial HD fraction (default 4e-10) |

### Datasets

One row is written every `METAL_OUTPUT_STRIDE` steps (default: 10).

| Dataset | Shape | Description |
|---|---|---|
| `step` | (N_rows,) int32 | Integration step number |
| `y` | (N_rows, 89) | Species abundances [/ xnH] |
| `xnH` | (N_rows,) | H number density [cm⁻³] |
| `T_K` | (N_rows,) | Gas temperature [K] |
| `T_gr_K` | (N_rows,) | Grain temperature [K] |
| `rho` | (N_rows,) | Mass density [g cm⁻³] |
| `xLmbd_net` | (N_rows,) | Net cooling [erg g⁻¹ s⁻¹] |
| `xLmbd_line` | (N_rows,) | Total line cooling |
| `xLmbd_cnt` | (N_rows,) | Total continuum cooling (grain + gas) |
| `xLmbd_gr` | (N_rows,) | Grain continuum cooling |
| `xLmbd_gas` | (N_rows,) | Gas continuum cooling (H ff + H₂ CIA) |
| `xLmbd_ch` | (N_rows,) | Chemical cooling |
| `xGam_cmp` | (N_rows,) | Compressional heating |
| `xLmbd_H2` | (N_rows,) | H₂ line cooling |
| `xLmbd_HD` | (N_rows,) | HD line cooling |
| `xLmbd_Lya` | (N_rows,) | Lyman-alpha cooling |
| `xLmbd_CO` | (N_rows,) | CO line cooling |
| `xLmbd_OH` | (N_rows,) | OH line cooling |
| `xLmbd_H2O` | (N_rows,) | H₂O line cooling |
| `xLmbd_CII` | (N_rows,) | C II line cooling |
| `xLmbd_CI` | (N_rows,) | C I line cooling |
| `xLmbd_OI` | (N_rows,) | O I line cooling |
| `xGam_CR` | (N_rows,) | CR ionization heating |
| `t_ff` | (N_rows,) | True free-fall time [s] |
| `t_cool` | (N_rows,) | Cooling time [s] |
| `t_chem` | (N_rows,) | Chemistry timescale [s] |
| `tau_cnt` | (N_rows,) | Continuum optical depth |
| `xlmbd_J` | (N_rows,) | Jeans length [cm] |
| `xMJ` | (N_rows,) | Jeans mass [g] |
| `B_cr` | (N_rows,) | Critical magnetic field [G] |
| `y_plus` | (N_rows,) | Total positive charge fraction |
| `y_minus` | (N_rows,) | Total negative charge fraction |
| `charge_imbal` | (N_rows,) | \|y⁺ − y⁻\| / (y⁺ + y⁻) |

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
log10_nH,T_K,xLmbd_net,y_0,y_1,...,y_N
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
| `xLmbd_net` | Net cooling rate [erg g⁻¹ s⁻¹] | arithmetic mean |
| `y_0` … `y_N` | Species abundances (prim: N=22, metal: N=88) | log-average |

Empty bins (no data points in range) are written as `NaN`.

### Averaging conventions

| Method | Formula | Used for |
|---|---|---|
| `log_avg` | 10^(mean(log₁₀(x))) for x > 0; NaN if all ≤ 0 | T_K, T_gr_K, species |
| `lin_avg` | mean(x) | xLmbd_net (signed/zero-capable) |

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
    xnH   = f["xnH"][:]
    T_K   = f["T_K"][:]
    Lnet  = f["xLmbd_net"][:]

    # Species abundances: shape (N_rows, 23)
    y         = f["y"][:]
    species   = f["y"].attrs["species"].split(",")
    idx_H2    = species.index("H2")
    y_H2      = y[:, idx_H2]
```
