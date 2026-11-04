# Release notes

## 2.0.1

ARCHE 2.0.1 fixes the switch of the grain reactions, the Newton iteration of the
metal network, the HD line cooling above T ≈ 6.8e3 K and the start of a Saha step
of the grain evaporation, and builds with Intel `icpx` on AMD processors. The
library API and ABI are those of 2.0.0. Numerical output is not bit-identical to
2.0.0; each fix below states how much the results change.

### Fixes

- The grain reactions of the metal networks (grain charge transfer and
  grain-surface reactions) are switched off when the grain volume falls below
  1e-10 of the volume of the cold grains, instead of at a grain temperature of
  1.025 max(T_ir, T_pyr, T_ol). The old switch is where an earlier, linear
  grain-volume law reached zero; with the current law about 2% of the grains
  remain there, and the step in the chemical cooling held the gas temperature at
  the switch: a metal Minimal collapse with Z = 1 Z☉, ζ = 1e-15 s⁻¹ and
  f_ret = 10 took 124 s and now takes 21 s. The rates already scale with the
  remaining grains, so the new switch is a negligible step. While grains remain
  they keep capturing electrons, which removes a jump of the electron abundance
  by up to 4 dex at the old switch. T changes by at most 2.8e-3 dex (metal
  Minimal, Z = 1 Z☉) and 2.8e-6 dex (metal full); the abundances of e⁻, Mg⁺ and
  the charged grains change within 0.1–0.7 dex of n_H around the grain
  evaporation (n_H ≈ 1e15–1e17 cm⁻³). The primordial networks are unchanged.
- The Newton iterations of the metal network (`Nakauchi2021`) include how the
  grain-surface rate coefficients depend on the physisorbed H2 and H2O and on the
  chemisorbed H and D (surface coverage). Without it the iteration for the
  physisorbed H2 oscillated where its desorption depends on its own coverage, and
  every step fell back to the slower fixed-temperature solve: a metal collapse with
  Z = 1 Z☉ and ζ = 1e-15 s⁻¹ took up to 24 min and now takes under 2 min. The
  equations are unchanged; the results change within the solver tolerance (T by
  less than 1e-3 dex, abundances above 1e-10 by less than 1e-2 dex, except where a
  steep front moves by at most 0.005 dex in n_H). The minimal networks and the
  primordial networks are unchanged.
- The HD line escape fraction is bounded by 1. Its α(T) fit turns negative above
  T ≈ 6.8e3 K, where the formula gives an escape fraction above 1 and, for a large
  HD column, overflows. Where the continuum is optically thick enough that
  exp(−τ_cnt) underflows to 0, the HD cooling is now 0 instead of 0 · ∞ = NaN. The
  NaN had stopped metal collapses with `solver_failed` (exit code 5) near
  n_H ≈ 2.6e22 cm⁻³ and T ≈ 1.7e4 K, for example Z = 1e-5 Z☉ with f_ret = 10.
  In collapses that did not stop, T changes by at most 9e-13 dex.
- A step of the grain evaporation in its Saha regime (n_H > 1e18 cm⁻³) starts
  from non-negative amounts of the metals still locked in the grains. On the last
  step of a complete evaporation, an amount that the release drains could start a
  rounding error below 0; the step was then rejected and the collapse stopped with
  `solver_failed` (exit code 5). Whether this happened depended on whether the
  build fuses multiply-adds: with `g++ -march=native` it stopped the metal
  collapses with f_ret = 10, Z = 1e-3, 1e-4 and 1e-6 Z☉ and ζ = 1e-19, 1e-17 and
  1e-15 s⁻¹, which now reach n_H = 1e20 cm⁻³. Only the start of that step
  changes; collapses that did not stop change by less than 1e-9 (relative).

### Build

- With Intel `icpx`, ARCHE compiles with `-march=native` instead of `-xHost`.
  The `-x` options embed an Intel-only processor check, so a binary built with
  `-xHost` aborts on AMD processors. On an Intel processor the two flags gave
  bit-identical output in the four collapses compared (primordial and metal,
  full and minimal networks).
- New CMake cache variable `ARCHE_MARCH` (default `native`) sets the
  `-march=<value>` flag; an empty value adds none. See
  [Getting Started](getting_started.md#target-cpu-arche_march).
- The CPU and floating-point flags (`-march`, and `-fp-model=precise` with
  `icpx`) are applied in every build type, not only `Release`. An unset
  `CMAKE_BUILD_TYPE` defaults to `Release` when ARCHE is the top-level project.
- These flags are no longer part of the usage requirements of `arche::arche`.
  A project that links libarche (installed, or through `add_subdirectory`)
  keeps its own CPU and floating-point flags; before, it received `-march=native`
  (GNU) or `-xHost -fp-model=precise` (`icpx`) in its `Release` build.

### Tests

- New collapse tests in [`tests/collapse/`](../tests/collapse/README.md): 160
  one-zone collapses over the four networks, ζ, f_ret, Z and J_LW21, each of which
  must reach n_H = 1e23 cm⁻³ normally within a wall-time gate, and a comparison of
  the temperature tracks of 16 of them across builds (`icpx`, and `g++` without
  `-march`). `tools/ci/run_local.sh` builds ARCHE and runs them on a local machine.
- The GitHub Actions workflow `.github/workflows/collapse-tests.yml` runs them on
  every pull request, together with an `ARCHE_XRAY=ON` build and the
  [`integration/`](../integration/README.md) consumers of an installed ARCHE.

## 2.0.0

ARCHE 2.0.0 solves every model with one error-controlled integrator, moves the
grain evaporation into the library, and adds Ohmic and ambipolar magnetic
diffusivities. The library API and ABI change incompatibly: code written against
1.x must be updated and rebuilt (see [Migrating from 1.x](#migrating-from-1x)).
Numerical output is not bit-identical to 1.1.0.

### Breaking changes

**Library API (C++ and C)**

- `chem_step_prim`, `chem_step_prim_minimal`, `chem_step_metal_minimal`,
  `chem_step_metal` and their return type `ChemRates` are removed. Use
  `chem_full_step_*` or the registry entry point `model_step`.
- `ChemShielding::tau_cnt` and `ChemShielding::esc_cnt` (C: the same fields of
  `ArcheChemShielding`) are removed. Pass the continuum shielding column
  `cnt_column = rho*L` [g cm⁻²] instead (0 = optically thin). The kernel derives
  the continuum depth from it at every state it evaluates and returns the end
  state's `tau_cnt` and `esc_cnt` in `ChemFullRates`.
- `ChemFullRates` changes meaning: the cooling and heating rates are means over
  the requested interval, and the opacities, continuum depth and grain
  temperature are those of the interval's end state. The new field `status`
  (`SolveStatus`) reports the internal steps and the failure reason.
- The metal networks carry the metals locked in the dust as species:
  - `Nakauchi2021`: 89 → 95 species (`C_s`, `O_s`, `Mg_s`, `Na_s`, `K_s`, `Gr_v`
    at indices 89–94; indices 0–88 are unchanged). `ARCHE_METAL_NSP` is 95.
  - `Nakauchi2021_Minimal`: 40 → 44 species (`C_s`, `O_s`, `Mg_s`, `Gr_v` at
    40–43). `ChemStateMetalMinimal` is `ChemState<44>`.
- The C structs `ArcheChemParams`, `ArcheChemShielding` and
  `ArcheChemFullRates` change layout.

**Build**

- The CMake option `ARCHE_SUBCYCLE` is removed.
- The installed package version is 2.0.0 with `SameMajorVersion` compatibility:
  `find_package(arche 1.x)` does not accept it. Request `find_package(arche 2)`
  or no version.

**Collapse applications**

- There is no temperature ceiling (1.x stopped at T > 10⁵ K). `exit_code` `2` is
  unused, and `5` now means that `chem_full_step` returned `solver_failed`.
- The metal output dataset `y` has 95 columns (44 for the minimal network).
  `schema_version` stays `2`: no dataset name changes.

### New features

- **Error-controlled integrator** (`solve/integrator.h`). Every model advances
  through backward-Euler Newton steps on the abundances and the internal energy,
  with the local error controlled by step doubling, a fixed-temperature fallback
  with a bracketed temperature root, and a Saha regime above
  `saha_threshold_nH`. New `ChemParams` fields: `rtol`, `atol`, `rtol_T`
  (defaults `1e-4`, `1e-14`, `1e-4`), `internal_max_steps`, `internal_max_step`,
  `use_saha`, `saha_threshold_nH`, `saha_nonthermal_tol`. The collapse
  applications read them from `PRIM_*` / `METAL_*` environment variables
  ([`parameters.md`](parameters.md)).
- **Step-size hint.** A call starts from the step size the previous call of the
  same cell proposed, when its conditions changed by at most 10%.
  `model_cell_reset` clears it.
- **Grain evaporation in the library** (`solve/grain_evaporation.h`). The metals
  locked in the dust return to the gas inside the call as the grains evaporate,
  and the conservation projection counts them toward their elements.
- **Magnetic diffusivities.** `arche::model_resistivity` / C
  `arche_model_resistivity` return the Ohmic and ambipolar diffusivities of a cell
  for all four models ([`magnetic_resistivity.md`](magnetic_resistivity.md)). The
  application `magnetic_reynolds` evaluates the magnetic Reynolds number along a
  collapse output, and `examples/nakauchi2021` reproduces Figure 11 of Nakauchi,
  Omukai & Susa (2021).
- **C API.** `ArcheChemParams` carries the integrator fields above;
  `ArcheChemFullRates` adds `tau_cnt`, `esc_cnt`, the `SolveStatus` counters,
  `saha_nonthermal_shift`, `input_charge_imbalance` and `input_neutralized`. The outcome codes in `failure_reason` and
  `last_failure` are named `ARCHE_SOLVE_*` ([`api_reference.md`](api_reference.md)
  §6).

### Fixes and behaviour changes

- Lyman-Werner photodissociation of H₂ and HD and photodetachment of H⁻ enter the
  implicit right-hand side and Jacobian of the network steps instead of being
  applied after each step.
- The H₂ self-shielding exponent α of Wolcott-Green & Haiman (2019, eqs. 7–8),
  fitted for T ≤ 8000 K, keeps its 8000 K value above 8000 K, so the shielding
  factor stays ≤ 1.
- Carriers of an element whose total is exactly zero are reset to zero before the
  conservation projection, so a rounding residue no longer blocks it (metal
  networks at Z = 0).
- The temperature inversion `model_T_from_e` uses secant iteration in log T with a
  bisection fallback; `T -> e -> T` round trips stay within 1.6e-12 relative for
  T ≥ 10 K.
- A charged input state is neutralised once, element totals kept, before the
  first internal step when its imbalance exceeds the error weights, and rejected
  with `std::invalid_argument` when that would shift a species by more than 1e-2
  of itself, unless the net charge is at most `atol`: such an imbalance is moved
  onto the electrons instead. `SolveStatus` reports `input_charge_imbalance` and
  `input_neutralized` (C: the same fields of `ArcheChemFullRates`). 1.1.0 repaired
  such a state at every step without reporting it.
- The metals locked in the dust are released when the grain temperature reaches
  0.975 T_pyr; 1.1.0 compared the gas temperature. Where the gas is much hotter
  than the grains the release moves to higher density: with
  `Nakauchi2021_Minimal` at ζ = 1e-17 s⁻¹ and Z = 1e-6, half of the Mg is back in
  the gas at n_H = 2.0e16 cm⁻³ instead of 4.8e11 cm⁻³ (where 1.1.0 had T = 1274 K
  and T_gr = 429 K). At Z = 1e-5 and 1e-3 the density changes by ≤ 1.1%.
- Outside its table the gas continuum opacity takes the table's edge value: the
  row's last tabulated value above the tabulated density range, and the 750 K value
  below 750 K. 1.1.0 returned 1e5 cm² g⁻¹ and 0 there; at n_H = 1e19 cm⁻³ the
  former is 8 decades above the edge value for T < 1310 K. Below log R = −5
  (R = ρ / T₆³) the opacity is still 0.
- The grain temperature solver keeps an initial guess that is already an exact
  root; it set T_gr = T_gas whenever the residual vanished at the guess (as in
  1.1.0). Warm-started from its own previous result, the solver hits an exact root
  most often in builds without FMA contraction (e.g. clang without
  `-march=native`), where the jump in grain cooling stopped metal collapses with
  e ≤ 0 before this fix.
- Collisional ionisation of H, He and He⁺ (`r0`, `r2`, `r4`): below
  T* = 11959 K, 15512 K and 85254 K the fit's prefactor k exp(ΔE/kT) is held at its
  value at T*. The eighth-order fits otherwise grow without bound below T* (by
  1e17, 1e28 and 1e112 at 1000 K), and so do the three-body recombinations derived
  from them by detailed balance.
- H + Li → LiH + γ (Bovino et al. 2011, reaction 3) is held at its value at its
  minimum, T* = 13722 K, above T*. The fit otherwise overflows near 1.2e5 K.
- Primordial reaction table: the detailed-balance constants (delE, Cmass) of 111
  of the 139 rows are the values implied by the Saha formation reactions, to six
  significant digits instead of three, so the network's equilibrium coincides with
  the Saha branch. delE changes by ≤ 1.2% and Cmass by ≤ 0.6%; at 866 K the reverse
  rates of the H, He and H₂ rows change by up to a factor 1.51
  (e + H₃⁺ → H + H₂).
- A call whose attempts keep being rejected for a cause other than the local
  error stops once those rejections shorten the step below the interval's
  resolution 8ε·dt, after at most ~25 shortenings, with
  `SolveFailure::no_progress`; `last_failure` is then the cause of the narrowest
  refused attempt. The error control itself may take the step below 8ε·dt (a
  stiff start). `SolveStatus` adds `rejections_by_reason`,
  `last_failure_run` and `last_failure_run_h` (C: the same fields, sized by
  `ARCHE_SOLVE_COUNT`).
- README Quick Start case 3 uses `J_LW21 = 1e5`.
- The HDF5 `description` attribute of the collapse outputs reads
  `"zero-metal one-zone collapse CR<tag>"` and
  `"metal_grain one-zone collapse CR<tag> Z<tag>"`.

### Known limitations

- Without the 1.x ceiling the temperature can exceed 1e6 K, where the gas
  opacity is extrapolated beyond its table and the rate coefficients and cooling
  functions are used above 1e5 K, the limit to which 1.x was run. With
  `J_LW21 = 1e5` (README case 3) T passes 1e6 K at n_H = 8.1e22 cm⁻³
  (`Nakauchi2019`; 5.5e22 cm⁻³ for `Nakauchi2019_Minimal`) and reaches 1.16e6 K
  (1.49e6 K) at 1e23 cm⁻³.
- `rtol`, `atol` and `rtol_T` bound the local error of each internal step, not the
  accuracy of a call: the error accumulated over a call falls only as about
  `rtol`^½ ([`api_reference.md`](api_reference.md) §8).
- Abundances below `atol` are not error-controlled. Where every charged particle
  lies below it, the ambipolar diffusivity depends on `atol` by several decades
  ([`magnetic_resistivity.md`](magnetic_resistivity.md)).

### Migrating from 1.x

1. Replace calls of `chem_step_*` by `chem_full_step_*` (or `model_step`) and read
   the scalar rates from `ChemFullRates`.
2. Stop computing `tau_cnt` / `esc_cnt` from the previous step's opacities: set
   `ChemShielding::cnt_column = rho*L` every call.
3. Size metal abundance vectors by `ARCHE_METAL_NSP` / `model_n_species()` (95 and
   44), and initialise the dust species as described in
   [`api_reference.md`](api_reference.md) §2 (grain evaporation). Look species up
   by name rather than by literal index (§4).
4. Initialise every C struct with `arche_chem_params_init` /
   `arche_chem_shielding_init` before setting fields, and rebuild against the 2.0
   headers and `libarche.a`.
5. Pass charge-neutral states; check `status.input_neutralized` to find where the
   caller produces charged ones.
6. On `solver_failed` the cell is unchanged and `status.reason` gives the cause.
   The kernel has already shortened its internal steps; a `no_progress` failure
   can still succeed when the interval is split into shorter calls
   ([`api_reference.md`](api_reference.md) §6).
7. Re-generate any stored reference outputs: results change at the level of the
   integrator tolerances.
