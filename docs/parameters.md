# Environment Variable Reference

For direct binary execution (`prim_collapse`, `metal_collapse`, `prim_chem`, `metal_chem`),
parameters are set as environment variables before invoking a binary.
Unset (or empty) optional variables fall back to their listed defaults.

For `run_collapse.sh`, parameters can be provided by:
- `params/default.conf` (loaded by default; required unless `--no-config`)
- `--config <file>` (custom parameter file)
- environment variables
- CLI options (`--prim-*`, `--metal-*`, ...)

Precedence in `run_collapse.sh` is:
1. CLI options
2. environment variables
3. config file

---

## `prim_collapse` — zero-metal one-zone collapse

| Variable | Required | Default | Description |
|---|---|---|---|
| `PRIM_ZETA0` | **yes** | — | CR ionization rate ζ₀ [s⁻¹].  Example: `1e-17` |
| `PRIM_OUTDIR` | no | `results/prim/h5` | Output directory for HDF5 files |
| `PRIM_FF_RET` | no | `1.0` | Free-fall retardation factor f_ret.  `1.0` = standard free-fall; `> 1` slows collapse.  Ignored when `PRIM_FRET_TABLE` is set |
| `PRIM_FRET_TABLE` | no | *(unset)* | Path to a 2-column ASCII table `nH [cm⁻³]  f_ret` (step-function, ratchet-forward).  When set, overrides `PRIM_FF_RET`; filename gets `_fret-step` suffix |
| `PRIM_XNH0` | no | `0.1` | Initial H number density [cm⁻³] |
| `PRIM_TK0` | no | `100.0` | Initial gas temperature [K] |
| `PRIM_YE0` | no | `1e-4` | Initial electron fraction y(e⁻) = y(H⁺).  Must be in [0, 1) |
| `PRIM_YH2` | no | `6e-7` | Initial H₂ fraction.  Must be in [0, 0.5) |
| `PRIM_YHD` | no | `4e-10` | Initial HD fraction.  Must be ≥ 0.  y(H) = 1 − y(H⁺) − 2·y(H₂) − y(HD) must remain > 0 |
| `PRIM_ABUNDANCE_SET` | no | `solar` | Abundance preset (`solar`, `default`, `alpha-enhanced`).  For primordial species, all presets currently use identical values |
| `PRIM_OUTPUT_STRIDE` | no | `100` | Write one HDF5 row every N integration steps |
| `PRIM_MAX_ITER` | no | `10000000` | Maximum number of integration steps |
| `PRIM_DT_FACTOR` | no | `1e-3` | Timestep factor for regular steps (> 0) |
| `PRIM_DT_FACTOR_INIT` | no | `1e-8` | Timestep factor for the first `PRIM_N_INIT_STEPS` steps (> 0) |
| `PRIM_N_INIT_STEPS` | no | `10` | Number of initial short-timestep steps (integer, ≥ 0) |
| `PRIM_XNH_STOP` | no | `1e23` | Stop threshold for H number density [cm⁻³] (> 0) |
| `PRIM_CR_ATTEN_COL_DENS` | no | `96.0` | CR attenuation column density scale [g cm⁻²] (> 0) |
| `PRIM_DATA_DIR` | no | compile-time `DATA_DIR` | Path to primordial reaction data directory |
| `PRIM_JLW21` | no | `0.0` | Lyman-Werner intensity J₂₁ [10⁻²¹ erg/s/cm²/Hz/sr].  `0.0` = no LW field.  Activates H₂/HD photodissociation and H⁻ photodetachment (operator-split after chemistry).  H₂ self-shielding uses Wolcott-Green & Haiman (2019) eq. 7–8 with density/temperature-dependent exponent α(n,T); HD uses WG2011 functional form (α = 2); H⁻ has no self-shielding |
| `PRIM_REDSHIFT` | no | `0.0` | Cosmological redshift z.  Sets T_rad = 2.725 × (1+z) K |

**Example — scalar f_ret**

```bash
PRIM_ZETA0=1e-17 \
PRIM_FF_RET=2.0 \
PRIM_XNH0=0.1 \
./build/src/apps/collapse_primordial/prim_collapse
```

Output file: `results/prim/h5/collapse_CR1e-17.h5`
(with f_ret=2.0: `collapse_CR1e-17_fret2p0.h5`)

**Example — step-function f_ret table**

```bash
PRIM_ZETA0=1e-17 \
PRIM_FRET_TABLE=data/fret_table/fret_step_sample.dat \
./build/src/apps/collapse_primordial/prim_collapse
```

Output file: `results/prim/h5/collapse_CR1e-17_fret-step.h5`

Table file format (`data/fret_table/fret_step_sample.dat`):

```
# nH [cm^-3]   f_ret
1e-2            3.0
1e6             1.0
```

Lines beginning with `#` are comments.  Rows must be in ascending `nH` order.
`f_ret` is held constant until `nH` reaches the next row's threshold (ratchet-forward, never backward).

---

## `metal_collapse` — metal-grain one-zone collapse

| Variable | Required | Default | Description |
|---|---|---|---|
| `METAL_ZETA0` | **yes** | — | CR ionization rate ζ₀ [s⁻¹] |
| `METAL_Z_METAL` | **yes** | — | Metallicity Z [Z☉].  Example: `1e-3` |
| `METAL_OUTDIR` | no | `results/metal/h5` | Output directory for HDF5 files |
| `METAL_FF_RET` | no | `1.0` | Free-fall retardation factor.  Ignored when `METAL_FRET_TABLE` is set |
| `METAL_FRET_TABLE` | no | *(unset)* | Path to a 2-column ASCII table `nH [cm⁻³]  f_ret` (step-function, ratchet-forward).  When set, overrides `METAL_FF_RET`; filename gets `_fret-step` suffix |
| `METAL_XNH0` | no | `1.0` | Initial H number density [cm⁻³] |
| `METAL_TK0` | no | `100.0` | Initial gas temperature [K] |
| `METAL_YE0` | no | `1e-4` | Initial electron fraction y(e⁻) = y(H⁺).  Must be in [0, 1) |
| `METAL_YH2` | no | `6e-7` | Initial H₂ fraction.  Must be in [0, 0.5) |
| `METAL_YHD` | no | `4e-10` | Initial HD fraction.  Must be ≥ 0.  y(H) = 1 − y(H⁺) − 2·y(H₂) − y(HD) must remain > 0 |
| `METAL_ABUNDANCE_SET` | no | `solar` | Abundance preset (`solar`, `default`, `alpha-enhanced`).  `alpha-enhanced` is a sensitivity-scan sample: O and Mg are scaled by +0.4 dex, others unchanged |
| `METAL_OUTPUT_STRIDE` | no | `10` | Write one HDF5 row every N integration steps |
| `METAL_MAX_ITER` | no | `1000000` | Maximum number of integration steps |
| `METAL_DT_FACTOR` | no | `1e-3` | Timestep factor for regular steps (> 0) |
| `METAL_DT_FACTOR_INIT` | no | `1e-8` | Timestep factor for the first `METAL_N_INIT_STEPS` steps (> 0) |
| `METAL_N_INIT_STEPS` | no | `10` | Number of initial short-timestep steps (integer, ≥ 0) |
| `METAL_XNH_STOP` | no | `1e23` | Stop threshold for H number density [cm⁻³] (> 0) |
| `METAL_CR_ATTEN_COL_DENS` | no | `96.0` | CR attenuation column density scale [g cm⁻²] (> 0) |
| `METAL_CR_ATTEN_SECOND_FRAC` | no | `7.6e-2` | Secondary CR attenuation fraction (≥ 0) |
| `METAL_CR_METAL_BKGND` | no | `1.4e-22` | Metal CR floor coefficient (≥ 0) |
| `METAL_T_CR_DES` | no | `70.0` | Effective CR desorption spike temperature [K] (advanced, > 0) |
| `METAL_C_GAS_FRAC` | no | `0.28` | Initial C gas-phase fraction ([0, 1]) |
| `METAL_O_GAS_FRAC` | no | `0.54` | Initial O gas-phase fraction ([0, 1]) |
| `METAL_MG_GAS_FRAC` | no | `0.02` | Initial Mg gas-phase fraction ([0, 1]) |
| `METAL_DATA_DIR` | no | compile-time value | Path to metal-grain reaction data directory |
| `PRIM_DATA_DIR` | no | compile-time value | Path to shared primordial data directory |
| `METAL_JLW21` | no | `0.0` | Lyman-Werner intensity J₂₁ [10⁻²¹ erg/s/cm²/Hz/sr].  `0.0` = no LW field.  Same physics as `PRIM_JLW21` (WG2019 H₂ self-shielding, WG2011 HD, no H⁻ shielding) |
| `METAL_REDSHIFT` | no | `0.0` | Cosmological redshift z |

**Example**

```bash
METAL_ZETA0=1e-17 METAL_Z_METAL=1e-3 \
./build/src/apps/collapse_metal_grain/metal_collapse
```

---

## `prim_chem` — zero-metal standalone chemistry demo

| Variable | Required | Default | Description |
|---|---|---|---|
| `PRIM_ZETA0` | **yes** | — | CR ionization rate ζ₀ [s⁻¹] |
| `PRIM_DATA_DIR` | no | compile-time `DATA_DIR` | Path to primordial reaction data directory |
| `PRIM_XNH` | no | `1e4` | Fixed H number density [cm⁻³] |
| `PRIM_T_K` | no | `100.0` | Fixed gas temperature [K] |
| `PRIM_YE0` | no | `1e-4` | Initial electron / H⁺ fraction.  Must be in [0, 1) |
| `PRIM_YH2` | no | `6e-7` | Initial H₂ fraction.  Must be in [0, 0.5) |
| `PRIM_YHD` | no | `4e-10` | Initial HD fraction.  Must be ≥ 0.  y(H) = 1 − y(H⁺) − 2·y(H₂) − y(HD) must remain > 0 |
| `PRIM_ABUNDANCE_SET` | no | `solar` | Abundance preset (`solar`, `default`, `alpha-enhanced`) |
| `PRIM_CR_ATTEN_COL_DENS` | no | `96.0` | CR attenuation column density scale [g cm⁻²] used in `zeta = zeta0 * exp(-rho * L / col)` (> 0) |
| `PRIM_REDSHIFT` | no | `0.0` | Cosmological redshift z.  Sets `T_rad = 2.725*(1+z)` K (≥ 0) |
| `PRIM_JLW21` | no | `0.0` | Lyman-Werner intensity J₂₁ for chemistry update (≥ 0) |
| `PRIM_NSTEPS` | no | `200` | Number of integration steps |
| `PRIM_DT` | no | `1e10` | Time step size [s] |

**Example**

```bash
PRIM_ZETA0=1e-17 PRIM_XNH=1e4 PRIM_T_K=100 \
./build/src/apps/chemistry_primordial/prim_chem
```

Output is written to stdout (step, t, species abundances, cooling rates).

> **Note — abundance preset intent**:
> `alpha-enhanced` is provided as a practical sensitivity-test sample,
> not as a uniquely recommended physical standard.

---

## `metal_chem` — metal-grain standalone chemistry demo

| Variable | Required | Default | Description |
|---|---|---|---|
| `METAL_ZETA0` | **yes** | — | CR ionization rate ζ₀ [s⁻¹] |
| `METAL_Z_METAL` | **yes** | — | Metallicity Z [Z☉] |
| `PRIM_DATA_DIR` | no | compile-time value | Path to shared primordial data directory |
| `METAL_DATA_DIR` | no | compile-time value | Path to metal-grain data directory |
| `METAL_XNH` | no | `1e4` | Fixed H number density [cm⁻³] |
| `METAL_T_K` | no | `100.0` | Fixed gas temperature [K] |
| `METAL_YE0` | no | `1e-4` | Initial electron / H⁺ fraction.  Must be in [0, 1) |
| `METAL_YH2` | no | `6e-7` | Initial H₂ fraction.  Must be in [0, 0.5) |
| `METAL_YHD` | no | `4e-10` | Initial HD fraction.  Must be ≥ 0.  y(H) = 1 − y(H⁺) − 2·y(H₂) − y(HD) must remain > 0 |
| `METAL_ABUNDANCE_SET` | no | `solar` | Abundance preset (`solar`, `default`, `alpha-enhanced`) |
| `METAL_CR_ATTEN_COL_DENS` | no | `96.0` | CR attenuation column density scale [g cm⁻²] (> 0) |
| `METAL_CR_ATTEN_SECOND_FRAC` | no | `7.6e-2` | Secondary CR attenuation fraction in shielding model (≥ 0) |
| `METAL_CR_METAL_BKGND` | no | `1.4e-22` | CR metal-floor coefficient in shielding model (≥ 0) |
| `METAL_T_CR_DES` | no | `70.0` | Effective CR desorption spike temperature [K] used in CR desorption terms (advanced, > 0) |
| `METAL_C_GAS_FRAC` | no | `0.28` | Initial C gas-phase fraction ([0,1]) |
| `METAL_O_GAS_FRAC` | no | `0.54` | Initial O gas-phase fraction ([0,1]) |
| `METAL_MG_GAS_FRAC` | no | `0.02` | Initial Mg gas-phase fraction ([0,1]) |
| `METAL_REDSHIFT` | no | `0.0` | Cosmological redshift z.  Sets `T_rad = 2.725*(1+z)` K (≥ 0) |
| `METAL_JLW21` | no | `0.0` | Lyman-Werner intensity J₂₁ for chemistry update (≥ 0) |
| `METAL_NSTEPS` | no | `200` | Number of integration steps |
| `METAL_DT` | no | `1e10` | Time step size [s] |

**Example**

```bash
METAL_ZETA0=1e-17 METAL_Z_METAL=1e-3 METAL_XNH=1e4 \
./build/src/apps/chemistry_metal_grain/metal_chem
```

---

## `run_collapse.sh` — all-in-one wrapper

`run_collapse.sh` builds the binaries, runs both collapse simulations, and calls
`analyze_collapse.py`.

By default, `run_collapse.sh` loads `params/default.conf` and errors out if that
file is missing (unless `--no-config` is specified).

### Synopsis

```bash
bash run_collapse.sh [options]
```

Run from the project root.

### Options

| Option | Argument | Default | Description |
|---|---|---|---|
| `--common-zeta0` | VALUE | from config (`COMMON_ZETA0`) | Shared CR rate for primordial + metal-grain [s⁻¹] |
| `--common-ff-ret` | VALUE | from config (`COMMON_FF_RET`) | Shared free-fall retardation factor (> 0) |
| `--common-fret-table` | FILE | from config (`COMMON_FRET_TABLE`) | Shared f_ret step-function table |
| `--common-jlw21` | VALUE | from config (`COMMON_JLW21`) | Shared Lyman-Werner J₂₁ |
| `--common-redshift` | VALUE | from config (`COMMON_REDSHIFT`) | Shared cosmological redshift |
| `--common-tk0` | VALUE | from config (`COMMON_TK0`) | Shared initial gas temperature [K] |
| `--common-ye0` | VALUE | from config (`COMMON_YE0`) | Shared initial electron fraction |
| `--common-yh2` | VALUE | from config (`COMMON_YH2`) | Shared initial H₂ fraction |
| `--common-yhd` | VALUE | from config (`COMMON_YHD`) | Shared initial HD fraction |
| `--common-abundance-set` | VALUE | from config (`COMMON_ABUNDANCE_SET`) | Shared abundance preset (`solar`, `default`, `alpha-enhanced`) |
| `--common-xnh0` | VALUE | from config (`COMMON_XNH0`) | Shared initial H number density [cm⁻³] |
| `--common-output-stride` | VALUE | from config (`COMMON_OUTPUT_STRIDE`) | Shared HDF5 output stride |
| `--common-max-iter` | VALUE | from config (`COMMON_MAX_ITER`) | Shared max iteration count |
| `--common-dt-factor` | VALUE | from config (`COMMON_DT_FACTOR`) | Shared timestep factor after initial steps (> 0) |
| `--common-dt-factor-init` | VALUE | from config (`COMMON_DT_FACTOR_INIT`) | Shared timestep factor during initial steps (> 0) |
| `--common-n-init-steps` | VALUE | from config (`COMMON_N_INIT_STEPS`) | Shared number of initial short-timestep steps (integer, ≥ 0) |
| `--common-xnh-stop` | VALUE | from config (`COMMON_XNH_STOP`) | Shared stop threshold for nH [cm⁻³] (> 0) |
| `--common-cr-col-dens` | VALUE | from config (`COMMON_CR_ATTEN_COL_DENS`) | Shared CR attenuation column density scale [g cm⁻²] (> 0) |
| `--common-data-dir` | DIR | from config (`COMMON_DATA_DIR`) | Shared reaction data directory |
| `--prim-zeta0` | VALUE | *(required unless `--no-prim`)* | Primordial CR rate [s⁻¹] |
| `--metal-zeta0` | VALUE | *(required unless `--no-metal`)* | Metal-grain CR rate [s⁻¹] |
| `--metal-z-metal` | VALUE | *(required unless `--no-metal`)* | Metallicity [Z☉] |
| `--prim-ff-ret` | VALUE | `1.0` | Primordial free-fall retardation factor (> 0).  Ignored when `--prim-fret-table` is given |
| `--prim-fret-table` | FILE | *(unset)* | Path to prim f_ret step-function table.  Sets `PRIM_FRET_TABLE`; output gets `_fret-step` suffix |
| `--metal-ff-ret` | VALUE | `1.0` | Metal-grain free-fall retardation factor (> 0).  Ignored when `--metal-fret-table` is given |
| `--metal-fret-table` | FILE | *(unset)* | Path to metal f_ret step-function table.  Sets `METAL_FRET_TABLE`; output gets `_fret-step` suffix |
| `--prim-jlw21` | VALUE | `0.0` | Lyman-Werner J₂₁ for primordial run [10⁻²¹ erg/s/cm²/Hz/sr] |
| `--metal-jlw21` | VALUE | `0.0` | Lyman-Werner J₂₁ for metal-grain run |
| `--prim-redshift` | VALUE | `0.0` | Cosmological redshift for primordial run |
| `--metal-redshift` | VALUE | `0.0` | Cosmological redshift for metal-grain run |
| `--prim-tk0` | VALUE | `100.0` | Primordial initial gas temperature [K] |
| `--prim-ye0` | VALUE | `1e-4` | Primordial initial electron fraction y(e⁻) = y(H⁺) |
| `--prim-yh2` | VALUE | `6e-7` | Primordial initial H₂ fraction |
| `--prim-yhd` | VALUE | `4e-10` | Primordial initial HD fraction |
| `--prim-abundance-set` | VALUE | from config (`PRIM_ABUNDANCE_SET` or `COMMON_ABUNDANCE_SET`) | Primordial abundance preset (`solar`, `default`, `alpha-enhanced`) |
| `--prim-xnh0` | VALUE | from config (`PRIM_XNH0` or `COMMON_XNH0`) | Primordial initial H number density [cm⁻³] |
| `--prim-output-stride` | VALUE | from config (`PRIM_OUTPUT_STRIDE` or `COMMON_OUTPUT_STRIDE`) | Primordial HDF5 output stride |
| `--prim-max-iter` | VALUE | from config (`PRIM_MAX_ITER` or `COMMON_MAX_ITER`) | Primordial max iteration count |
| `--prim-dt-factor` | VALUE | from config (`PRIM_DT_FACTOR` or `COMMON_DT_FACTOR`) | Primordial timestep factor after initial steps (> 0) |
| `--prim-dt-factor-init` | VALUE | from config (`PRIM_DT_FACTOR_INIT` or `COMMON_DT_FACTOR_INIT`) | Primordial timestep factor during initial steps (> 0) |
| `--prim-n-init-steps` | VALUE | from config (`PRIM_N_INIT_STEPS` or `COMMON_N_INIT_STEPS`) | Primordial number of initial short-timestep steps (integer, ≥ 0) |
| `--prim-xnh-stop` | VALUE | from config (`PRIM_XNH_STOP` or `COMMON_XNH_STOP`) | Primordial stop threshold for nH [cm⁻³] (> 0) |
| `--prim-cr-col-dens` | VALUE | from config (`PRIM_CR_ATTEN_COL_DENS` or `COMMON_CR_ATTEN_COL_DENS`) | Primordial CR attenuation column density [g cm⁻²] (> 0) |
| `--prim-data-dir` | DIR | from config (`PRIM_DATA_DIR` or `COMMON_DATA_DIR`) | Primordial reaction data directory |
| `--metal-tk0` | VALUE | `100.0` | Metal-grain initial gas temperature [K] |
| `--metal-ye0` | VALUE | `1e-4` | Metal-grain initial electron fraction y(e⁻) = y(H⁺) |
| `--metal-yh2` | VALUE | `6e-7` | Metal-grain initial H₂ fraction |
| `--metal-yhd` | VALUE | `4e-10` | Metal-grain initial HD fraction |
| `--metal-abundance-set` | VALUE | from config (`METAL_ABUNDANCE_SET` or `COMMON_ABUNDANCE_SET`) | Metal-grain abundance preset (`solar`, `default`, `alpha-enhanced`) |
| `--metal-xnh0` | VALUE | from config (`METAL_XNH0` or `COMMON_XNH0`) | Metal-grain initial H number density [cm⁻³] |
| `--metal-output-stride` | VALUE | from config (`METAL_OUTPUT_STRIDE` or `COMMON_OUTPUT_STRIDE`) | Metal-grain HDF5 output stride |
| `--metal-max-iter` | VALUE | from config (`METAL_MAX_ITER` or `COMMON_MAX_ITER`) | Metal-grain max iteration count |
| `--metal-dt-factor` | VALUE | from config (`METAL_DT_FACTOR` or `COMMON_DT_FACTOR`) | Metal-grain timestep factor after initial steps (> 0) |
| `--metal-dt-factor-init` | VALUE | from config (`METAL_DT_FACTOR_INIT` or `COMMON_DT_FACTOR_INIT`) | Metal-grain timestep factor during initial steps (> 0) |
| `--metal-n-init-steps` | VALUE | from config (`METAL_N_INIT_STEPS` or `COMMON_N_INIT_STEPS`) | Metal-grain number of initial short-timestep steps (integer, ≥ 0) |
| `--metal-xnh-stop` | VALUE | from config (`METAL_XNH_STOP` or `COMMON_XNH_STOP`) | Metal-grain stop threshold for nH [cm⁻³] (> 0) |
| `--metal-cr-col-dens` | VALUE | from config (`METAL_CR_ATTEN_COL_DENS` or `COMMON_CR_ATTEN_COL_DENS`) | Metal-grain CR attenuation column density [g cm⁻²] (> 0) |
| `--metal-cr-second-frac` | VALUE | from config (`METAL_CR_ATTEN_SECOND_FRAC`) | Metal-grain secondary CR attenuation fraction (≥ 0) |
| `--metal-cr-metal-bkgnd` | VALUE | from config (`METAL_CR_METAL_BKGND`) | Metal-grain CR floor coefficient (≥ 0) |
| `--metal-t-cr-des` | VALUE | from config (`METAL_T_CR_DES`) | Metal-grain CR desorption spike temperature [K] (advanced, > 0) |
| `--metal-c-gas-frac` | VALUE | from config (`METAL_C_GAS_FRAC`) | Initial C gas-phase fraction ([0, 1]) |
| `--metal-o-gas-frac` | VALUE | from config (`METAL_O_GAS_FRAC`) | Initial O gas-phase fraction ([0, 1]) |
| `--metal-mg-gas-frac` | VALUE | from config (`METAL_MG_GAS_FRAC`) | Initial Mg gas-phase fraction ([0, 1]) |
| `--metal-data-dir` | DIR | from config (`METAL_DATA_DIR` or `COMMON_DATA_DIR`) | Metal-grain reaction data directory |
| `--build-dir` | DIR | from config (`BUILD_DIR`, fallback `build`) | CMake build directory |
| `--out-dir` | DIR | from config (`OUT_DIR`, fallback `results/`) | Root for HDF5 output (`DIR/primordial/`, `DIR/metal_grain/`) |
| `--save-dir` | DIR | from config (`SAVE_DIR`, fallback `results/`) | Root for figure output (`DIR/primordial/`, `DIR/metal_grain/`) |
| `--config` | FILE | `params/default.conf` | Parameter file path |
| `--no-config` | — | — | Disable parameter file loading |
| `--no-build` | — | — | Skip CMake build step |
| `--no-prim` | — | — | Skip primordial collapse |
| `--no-metal` | — | — | Skip metal-grain collapse |
| `--no-plot` | — | — | Skip `analyze_collapse.py` |
| `--fig-combo` | — | — | Also produce summary figure (`fig_summary_<tag>.png`) |
| `--resample` | — | enabled | Force resampling even when `--no-plot` is set |
| `--no-resample` | — | — | Skip `resample_collapse.py` |

### Config keys (`params/default.conf`)

`run_collapse.sh` accepts the same `PRIM_*` / `METAL_*` variables listed above,
plus optional shared keys:

Policy:
- `COMMON_*` is allowed only for keys used by both primordial and metal runs.
- For metal-only keys, use `METAL_*` only (no `COMMON_*` counterpart).

`params/default.conf` is organized into:
- General operation keys (top/middle blocks)
- Advanced keys (expert knobs) collected at file end

| Key | Description |
|---|---|
| `COMMON_ZETA0` | Shared CR ionization rate (`PRIM_ZETA0`, `METAL_ZETA0`) |
| `COMMON_FF_RET` | Shared free-fall retardation (`PRIM_FF_RET`, `METAL_FF_RET`) |
| `COMMON_FRET_TABLE` | Shared f_ret step-function table |
| `COMMON_JLW21` | Shared Lyman-Werner intensity |
| `COMMON_REDSHIFT` | Shared cosmological redshift |
| `COMMON_TK0` | Shared initial gas temperature |
| `COMMON_YE0` | Shared initial electron fraction |
| `COMMON_YH2` | Shared initial H₂ fraction |
| `COMMON_YHD` | Shared initial HD fraction |
| `COMMON_ABUNDANCE_SET` | Shared abundance preset (`solar`, `default`, `alpha-enhanced`) |
| `COMMON_XNH0` | Shared initial H number density |
| `COMMON_OUTPUT_STRIDE` | Shared output stride |
| `COMMON_MAX_ITER` | Shared max iteration count |
| `COMMON_DT_FACTOR` | Shared timestep factor after initial steps |
| `COMMON_DT_FACTOR_INIT` | Shared timestep factor during initial steps |
| `COMMON_N_INIT_STEPS` | Shared number of initial short-timestep steps |
| `COMMON_XNH_STOP` | Shared stop threshold for nH |
| `COMMON_CR_ATTEN_COL_DENS` | Shared CR attenuation column density scale |
| `COMMON_DATA_DIR` | Shared reaction data directory |

If both `COMMON_*` and run-specific keys are present in config, run-specific
keys (`PRIM_*`, `METAL_*`) take precedence for that run.

> **Policy — `METAL_T_CR_DES`**:
> Keep the default (`70.0`) for routine runs.
> Override only for comparison experiments, and record the chosen value in the run log / notes.

> **Note — initial conditions**: `PRIM_TK0 / PRIM_YE0 / PRIM_YH2 / PRIM_YHD` and
> `METAL_TK0 / METAL_YE0 / METAL_YH2 / METAL_YHD` can be set either as CLI options
> (see table above) **or** as environment variables exported before calling the script;
> the CLI option takes precedence.

> **Note — `PRIM_FRET_TABLE` / `METAL_FRET_TABLE`**: these can also be exported
> before calling `run_collapse.sh` as an alternative to `--prim-fret-table` /
> `--metal-fret-table`.

> **Note — parameter precedence**:
> CLI options > environment variables > config file (`default.conf` or `--config`).
>
> For shared CLI parameters, `--prim-*` / `--metal-*` override `--common-*` on each side.

### Examples

```bash
# With default config (`params/default.conf`)
bash run_collapse.sh --prim-zeta0 1e-17 --metal-zeta0 1e-17 --metal-z-metal 1e-3

# With custom config file
bash run_collapse.sh --config params/default.conf \
    --prim-zeta0 1e-17 --metal-zeta0 1e-17 --metal-z-metal 1e-3

# Disable config loading (legacy env/CLI-only mode)
bash run_collapse.sh --no-config \
    --prim-zeta0 1e-17 --metal-zeta0 1e-17 --metal-z-metal 1e-3

# Skip build and primordial; only re-run metal and plot
bash run_collapse.sh --no-build --no-prim \
    --metal-zeta0 1e-16 --metal-z-metal 1.2e-4

# With free-fall retardation
bash run_collapse.sh --prim-zeta0 1e-17 --prim-ff-ret 3.0 \
                   --metal-zeta0 1e-17 --metal-z-metal 1e-3 --metal-ff-ret 3.0

# With cosmological redshift
bash run_collapse.sh --prim-zeta0 1e-17 --prim-redshift 20 \
                   --metal-zeta0 1e-17 --metal-z-metal 1e-3 --metal-redshift 20

# With Lyman-Werner radiation field (J_21 = 1e3, DCBH scenario)
bash run_collapse.sh --prim-zeta0 1e-17 --prim-jlw21 1e3 \
                   --metal-zeta0 1e-17 --metal-z-metal 1e-3 --metal-jlw21 1e3

# Build and simulate only (no plots)
bash run_collapse.sh --prim-zeta0 1e-17 --metal-zeta0 1e-17 --metal-z-metal 1e-3 \
                   --no-plot

# Custom initial conditions (warm, pre-ionised gas)
bash run_collapse.sh --prim-zeta0 0 \
                   --prim-tk0 500 --prim-ye0 1e-3 --prim-yh2 1e-5 \
                   --no-metal
```
