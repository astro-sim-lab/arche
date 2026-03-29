# `resample_collapse.py` — Resampling Reference

`tools/resample_collapse.py` reads one HDF5 collapse output file and resamples
all physical quantities onto a uniform log₁₀(nH) grid, writing a CSV table.

Physical outputs from `prim_collapse` / `metal_collapse` are sampled at
adaptive time steps, so data points cluster in density regions where quantities
change rapidly.  The resampled CSV provides a uniform-spacing table suitable
for:

- Cooling table input to 3-D hydro / semi-analytic codes
- Quantitative comparison across runs with different CR rates, f_ret, or Z
- Data compression / archiving (e.g. thousands of rows → ~100 uniform bins)

---

## Synopsis

```bash
python3 tools/resample_collapse.py [prim options] [metal options] [grid options]
```

At least one of `--h5dir`+`--cr-tag` (prim) or `--metal-h5dir`+`--metal-cr`+`--metal-z`
(metal) must be supplied.

---

## Primordial options

| Option | Argument | Default | Description |
|---|---|---|---|
| `--h5dir` | DIR | — | Directory containing `collapse_CR<tag>.h5` |
| `--cr-tag` | TAG | — | CR tag, e.g. `1e-17`, `0`, `1p5e-18` |
| `--fret-tag` | TAG | `''` | f_ret tag.  For scalar f_ret use e.g. `3p0`; for step-function table use `step` (→ `_fret-step` in filename).  Omit or `''` for f_ret = 1 |
| `--jlw-tag` | TAG | `''` | Lyman-Werner tag (omit or `''` for J_LW21 = 0).  Must match the `_JLW<tag>` suffix in the HDF5 filename |
| `--zred-tag` | TAG | `''` | Redshift tag (omit or `''` for z = 0) |
| `--save` | DIR | same as `--h5dir` | Output directory for CSV |

---

## Metal options

| Option | Argument | Default | Description |
|---|---|---|---|
| `--metal-h5dir` | DIR | — | Directory containing `collapse_CR<cr>_Z<z>.h5` |
| `--metal-cr` | TAG | — | CR tag for metal run |
| `--metal-z` | TAG | — | Metallicity tag, e.g. `1e-3`, `1p2e-4` |
| `--metal-fret` | TAG | `''` | f_ret tag for metal run (same rules as `--fret-tag`) |
| `--metal-jlw` | TAG | `''` | Lyman-Werner tag for metal run (omit or `''` for J_LW21 = 0) |
| `--metal-zred` | TAG | `''` | Redshift tag for metal run |
| `--metal-save` | DIR | same as `--metal-h5dir` | Output directory for CSV |

---

## Grid options

| Option | Type | Default | Description |
|---|---|---|---|
| `--log-nH-min` | float | `-2.0` | log₁₀(nH) lower edge |
| `--log-nH-max` | float | `23.0` | log₁₀(nH) upper edge |
| `--log-nH-step` | float | `0.1` | Bin width in log₁₀(nH) |

The default grid produces 50 bins spanning nH = 10⁻² … 10²³ cm⁻³.

---

## Output file naming

The output file name mirrors the input HDF5 stem with a `resample_` prefix:

```
resample_collapse_CR<cr_tag>[_fret<fret_tag>][_JLW<jlw_tag>][_z<zred_tag>].csv          # prim
resample_collapse_CR<cr_tag>_Z<z_tag>[_fret<fret_tag>][_JLW<jlw_tag>][_z<zred_tag>].csv # metal
```

Suffix order: `_fret` → `_JLW` → `_z`.

---

## Output CSV format

```
# resample_collapse.py output
# species: H,H2,e-,H+,...
# NaN = empty bin (no data points)
log10_nH,T_K,xLmbd_net,y_0,y_1,...,y_N
-1.750000,2.726000e+01,-1.234567e-20,...
```

- Lines 1–3 start with `#` (comments).
- Line 2 lists species names in the same order as the `y_*` columns.
- Line 4 is the CSV header (no `#`).
- Subsequent lines are data rows, one per bin.
- Empty bins are written as `NaN`.

### Columns

| Column | Description | Averaging method |
|---|---|---|
| `log10_nH` | Bin centre in log₁₀(nH [cm⁻³]) | — (grid index variable) |
| `T_K` | Gas temperature [K] | log-average |
| `T_gr_K` | Grain temperature [K] (**metal only**) | log-average |
| `xLmbd_net` | Net cooling rate [erg g⁻¹ s⁻¹] | arithmetic mean |
| `y_0` … `y_22` | Species abundances — prim (23 species) | log-average |
| `y_0` … `y_88` | Species abundances — metal (89 species) | log-average |

---

## Lyman-Werner self-shielding model

The underlying collapse simulation uses the following self-shielding models when `J_LW21 > 0`.
The resampled CSV reproduces these rates as part of the species/cooling columns.

| Species | Model | Reference |
|---|---|---|
| H₂ | WG2019 eq. 7–8: α(n,T)-dependent | Wolcott-Green & Haiman (2019) MNRAS 484, 2467 |
| HD | WG2011 (α = 2 fixed) | Wolcott-Green & Haiman (2011) |
| H⁻ | No self-shielding | Tegmark et al. (1997) |

---

## Averaging conventions

| Method | Formula | Applied to |
|---|---|---|
| `log_avg` | 10^(mean(log₁₀(x))) for x > 0; NaN if all ≤ 0 | T_K, T_gr_K, all `y_*` columns |
| `lin_avg` | mean(x) | `xLmbd_net` (can be negative or zero) |

The mapping is defined in the `AVG_FUNCS` dict near the top of the script:

```python
AVG_FUNCS: Dict[str, Callable[[np.ndarray], float]] = {
    "xLmbd_net": lin_avg,
    # Add overrides here, e.g.:
    #   "xGam_cmp": lin_avg,
}
```

Any column not listed falls back to `log_avg`.

---

## Examples

**Primordial — default grid**

```bash
python3 tools/resample_collapse.py \
    --h5dir results/prim/h5 --cr-tag 1e-17
```

Output: `results/prim/h5/resample_collapse_CR1e-17.csv`

**Primordial — step-function f_ret table**

```bash
python3 tools/resample_collapse.py \
    --h5dir results/prim/h5 --cr-tag 1e-17 --fret-tag step
```

Output: `results/prim/h5/resample_collapse_CR1e-17_fret-step.csv`

> **Note:** Use `--fret-tag step` (no leading dash) for table mode.
> Internally this is converted to the `_fret-step` filename suffix.
> For a scalar f_ret (e.g. f_ret=3.0) use `--fret-tag 3p0`.

**Primordial — scalar f_ret**

```bash
python3 tools/resample_collapse.py \
    --h5dir results/prim/h5 --cr-tag 1e-17 --fret-tag 3p0
```

Output: `results/prim/h5/resample_collapse_CR1e-17_fret3p0.csv`

**Primordial — with redshift**

```bash
python3 tools/resample_collapse.py \
    --h5dir results/prim/h5 --cr-tag 1e-17 --zred-tag 20p0 \
    --save results/prim/csv
```

Output: `results/prim/csv/resample_collapse_CR1e-17_z20p0.csv`

**Primordial — with Lyman-Werner field (J_21 = 1e3)**

```bash
python3 tools/resample_collapse.py \
    --h5dir results/prim/h5 --cr-tag 1e-17 --jlw-tag 1e3
```

Output: `results/prim/h5/resample_collapse_CR1e-17_JLW1e3.csv`

**Metal-grain**

```bash
python3 tools/resample_collapse.py \
    --metal-h5dir results/metal/h5 \
    --metal-cr 1e-17 --metal-z 1e-3
```

Output: `results/metal/h5/resample_collapse_CR1e-17_Z1e-3.csv`

**Metal-grain — with Lyman-Werner field**

```bash
python3 tools/resample_collapse.py \
    --metal-h5dir results/metal/h5 \
    --metal-cr 1e-17 --metal-z 1e-3 --metal-jlw 1e3
```

Output: `results/metal/h5/resample_collapse_CR1e-17_Z1e-3_JLW1e3.csv`

**Metal-grain — step-function f_ret table**

```bash
python3 tools/resample_collapse.py \
    --metal-h5dir results/metal/h5 \
    --metal-cr 1e-17 --metal-z 1e-3 --metal-fret step
```

Output: `results/metal/h5/resample_collapse_CR1e-17_Z1e-3_fret-step.csv`

**Both networks — finer grid, custom save directory**

```bash
python3 tools/resample_collapse.py \
    --h5dir results/prim/h5 --cr-tag 1e-17 \
    --save results/prim/csv \
    --metal-h5dir results/metal/h5 --metal-cr 1e-17 --metal-z 1e-3 \
    --metal-save results/metal/csv \
    --log-nH-min -2.0 --log-nH-max 23.0 --log-nH-step 0.25
```

---

## Reading the CSV with Python

```python
import numpy as np

data = np.genfromtxt(
    "results/prim/h5/resample_collapse_CR1e-17.csv",
    delimiter=",",
    names=True,       # read header line
    comments="#",
)

log_nH  = data["log10_nH"]
T_K     = data["T_K"]
Lnet    = data["xLmbd_net"]
y_H2    = data["y_1"]          # index 1 = H2 in zero-metal network
```

---

## Dependencies

```
python >= 3.11
numpy >= 2.0
h5py
```
