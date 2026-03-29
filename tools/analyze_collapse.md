# `analyze_collapse.py` — Visualisation Reference

`tools/analyze_collapse.py` reads one HDF5 file per collapse run and produces
four standard figures per simulation.  It supports both the zero-metal
(`prim_collapse`) and metal-grain (`metal_collapse`) networks simultaneously.

---

## Synopsis

```bash
python3 tools/analyze_collapse.py [primordial options] [metal options] [common options]
```

---

## Primordial options

| Option | Argument | Description |
|---|---|---|
| `--h5dir` | DIR | Directory containing `collapse_CR<tag>.h5` (output of `prim_collapse`) |
| `--cr-tag` | TAG | CR tag string, e.g. `1e-17`, `0`, `1p5e-18` |
| `--fret-tag` | TAG | Free-fall retardation tag (omit or `''` for f_ret = 1) |
| `--jlw-tag` | TAG | Lyman-Werner tag (omit or `''` for J_LW21 = 0).  Must match the `_JLW<tag>` suffix in the HDF5 filename |
| `--zred-tag` | TAG | Redshift tag for primordial run (omit or `''` for z = 0) |
| `--save` | DIR | Directory for output PNG files (default: `results/prim/fig`) |

---

## Metal options

| Option | Argument | Description |
|---|---|---|
| `--metal-h5dir` | DIR | Directory containing `collapse_CR<cr>_Z<z>.h5` |
| `--metal-cr` | TAG | CR tag for metal run |
| `--metal-z` | TAG | Metallicity tag, e.g. `1e-3`, `1p2e-4` |
| `--metal-fret` | TAG | Free-fall retardation tag for metal (omit or `''` for f_ret = 1) |
| `--metal-jlw` | TAG | Lyman-Werner tag for metal run (omit or `''` for J_LW21 = 0) |
| `--metal-zred` | TAG | Redshift tag for metal run (omit or `''` for z = 0) |
| `--metal-save` | DIR | Directory for output PNG files (default: `results/metal/fig`) |

---

## Common options

| Option | Description |
|---|---|
| `--show` | Display figures interactively instead of saving to PNG |
| `--fig-combo` | Also output a single-panel summary figure (`fig_summary_<tag>.png`) |

---

## Output figures

Each network produces up to 5 figures.  File names encode the parameter tags
so multiple runs can coexist in the same output directory.

### Primordial (saved to `--save`)

| File | Contents |
|---|---|
| `fig1_phase_cooling_CR<tag>[_JLW<jlw>][_z<zred>].png` | Temperature–density phase diagram and net cooling rate |
| `fig2_species_CR<tag>[_JLW<jlw>][_z<zred>].png` | Species abundance evolution vs density |
| `fig3_jeans_CR<tag>[_JLW<jlw>][_z<zred>].png` | Jeans length, Jeans mass, and critical magnetic field vs density |
| `fig4_thermal_balance_CR<tag>[_JLW<jlw>][_z<zred>].png` | Cooling and heating rate breakdown vs density |
| `fig_summary_CR<tag>[_JLW<jlw>][_z<zred>].png` | Single-panel summary (requires `--fig-combo`) |

The `[_fret<fr>]` suffix appears between `CR<tag>` and `[_JLW<jlw>]` when f_ret ≠ 1.

### Metal-grain (saved to `--metal-save`)

| File | Contents |
|---|---|
| `fig1_phase_cooling_CR<cr>_Z<z>[_JLW<jlw>][_z<zred>].png` | T–n phase diagram; includes grain temperature |
| `fig2_species_CR<cr>_Z<z>[_JLW<jlw>][_z<zred>].png` | Species abundance evolution (H, H₂, CO, key metals) |
| `fig3_jeans_CR<cr>_Z<z>[_JLW<jlw>][_z<zred>].png` | Jeans parameters and critical magnetic field |
| `fig4_thermal_balance_CR<cr>_Z<z>[_JLW<jlw>][_z<zred>].png` | Full cooling breakdown including CO, OH, H₂O, CII, CI, OI, grain |
| `fig_summary_CR<cr>_Z<z>[_JLW<jlw>][_z<zred>].png` | Single-panel summary (requires `--fig-combo`) |

---

## Lyman-Werner self-shielding model

When `J_LW21 > 0`, the operator-split LW block in `chemistry.h` applies:

| Species | Self-shielding | Reference |
|---|---|---|
| H₂ | WG2019 eq. 7–8: α(n,T)-dependent exponent | Wolcott-Green & Haiman (2019) MNRAS 484, 2467 |
| HD | WG2011 functional form, α = 2 fixed | Wolcott-Green & Haiman (2011) |
| H⁻ | None (continuum absorption, no band structure) | Tegmark et al. (1997) |

The WG2019 formula uses the **gas** density and temperature (not source temperature) to compute α(n,T), capturing the weakening of self-shielding at high n and T where excited H₂ rovibrational states are significantly populated.
Valid range: n ≤ 10⁷ cm⁻³, T ≤ 8000 K, N(H₂) ≤ 10¹⁷ cm⁻².

---

## Tag format

Tags are formed from parameter values by replacing `.` with `p`.
Zero values are represented as `0`.

| Parameter value | Tag |
|---|---|
| `1e-17` | `1e-17` |
| `1.5e-18` | `1p5e-18` |
| `1e-3` | `1e-3` |
| `1.2e-4` | `1p2e-4` |
| `0` | `0` |

---

## Examples

**Primordial only**

```bash
python3 tools/analyze_collapse.py \
    --h5dir results/prim/h5 --cr-tag 1e-17 \
    --save results/prim/fig
```

**Metal-grain only**

```bash
python3 tools/analyze_collapse.py \
    --metal-h5dir results/metal/h5 --metal-cr 1e-17 --metal-z 1e-3 \
    --metal-save results/metal/fig
```

**Both networks together**

```bash
python3 tools/analyze_collapse.py \
    --h5dir results/prim/h5 --cr-tag 1e-17 \
    --metal-h5dir results/metal/h5 --metal-cr 1e-17 --metal-z 1e-3 \
    --save results/prim/fig --metal-save results/metal/fig
```

**With free-fall retardation (f_ret = 2.0)**

```bash
python3 tools/analyze_collapse.py \
    --h5dir results/prim/h5 --cr-tag 1e-17 --fret-tag 2p0 \
    --save results/prim/fig
```

**With Lyman-Werner radiation field (J_21 = 1e3)**

```bash
python3 tools/analyze_collapse.py \
    --h5dir results/prim/h5 --cr-tag 1e-17 --jlw-tag 1e3 \
    --metal-h5dir results/metal/h5 --metal-cr 1e-17 --metal-z 1e-3 --metal-jlw 1e3 \
    --save results/prim/fig --metal-save results/metal/fig
```

Output (prim): `fig1_phase_cooling_CR1e-17_JLW1e3.png`, …
Output (metal): `fig1_phase_cooling_CR1e-17_Z1e-3_JLW1e3.png`, …

The figure title automatically shows `J_{LW} = 1000 J_{21}` (read from HDF5 `J_LW21` attribute).

**With redshift (z = 20)**

```bash
python3 tools/analyze_collapse.py \
    --h5dir results/prim/h5 --cr-tag 1e-17 --zred-tag 20 \
    --save results/prim/fig
```

**Interactive display**

```bash
python3 tools/analyze_collapse.py \
    --h5dir results/prim/h5 --cr-tag 1e-17 --show
```

**Summary figure**

```bash
python3 tools/analyze_collapse.py \
    --h5dir results/prim/h5 --cr-tag 1e-17 \
    --metal-h5dir results/metal/h5 --metal-cr 1e-17 --metal-z 1e-3 \
    --save results/prim/fig --metal-save results/metal/fig \
    --fig-combo
```

---

## See also

- [`tools/resample_collapse.md`](resample_collapse.md) — `tools/resample_collapse.py`: resample HDF5 onto a
  uniform log₁₀(nH) grid and export a CSV cooling/chemistry table.

---

## Using via `run_collapse.sh`

`run_collapse.sh` calls `analyze_collapse.py` automatically in step 4.
The `--fig-combo` flag is forwarded directly:

```bash
bash run_collapse.sh --prim-zeta0 1e-17 --metal-zeta0 1e-17 --metal-z-metal 1e-3 \
                   --fig-combo
```

To skip the plot step: `--no-plot`.
To run the plot step on previously generated HDF5 files without re-running the
simulation: `--no-build --no-prim --no-metal`.

---

## Dependencies

```
python >= 3.11
numpy >= 2.0
matplotlib
h5py
```

`yt` is not required by this script at present.
