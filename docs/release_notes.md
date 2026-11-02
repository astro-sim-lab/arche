# Release notes

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
  of itself. `SolveStatus` reports `input_charge_imbalance` and
  `input_neutralized` (C: the same fields of `ArcheChemFullRates`). 1.1.0 repaired
  such a state at every step without reporting it.
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
6. Treat `solver_failed` as final for that call: the kernel has already shortened
   its internal steps, and `status.reason` gives the cause.
7. Re-generate any stored reference outputs: results change at the level of the
   integrator tolerances.
