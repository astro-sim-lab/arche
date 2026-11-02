# Benchmark: Per-Step Computation Time

This document summarizes the per-step wall-clock timing of the chemistry
kernel (`chem_full_step`) across the primordial and metal-grain networks,
measured as a function of H number density and metallicity.  It also
includes a simple estimate based on the chemistry timescale `t_chem`.

---

## Purpose

Quantify the relative computational cost of each network variant so that
the overhead of the metal-grain solver (95 species) compared to the
zero-metal solver (23 species) can be estimated as a function of
metallicity Z/Z☉ and density nH.

---

## Measurement record

The per-step wall times come from the `PRIM_BENCH` / `METAL_BENCH` output of the
collapse applications ([`parameters.md`](parameters.md)), which writes the time of
every `chem_full_step` call to a CSV file. The scripts that run the benchmark set
and draw the figures below are not distributed with ARCHE, so the figures are a
record of one measurement rather than a procedure to reproduce:

| Item | Value |
|---|---|
| Code | ARCHE 2.0.0 development tree (integrator as released) |
| Compiler | Intel oneAPI icpx/icx 2025.0, `CMAKE_BUILD_TYPE=Release` (`-xHost -fp-model=precise`) |
| Runs | 20 collapses (prim and metal, full and minimal), one at a time, once each |
| Machine | not recorded |

Timings depend on the machine and the load; compare ratios between networks
measured in the same session, not absolute times across machines.

---

## Benchmark conditions

| Parameter | Value |
|---|---|
| CR ionization rate ζ₀ | 0 (no cosmic rays) |
| Free-fall retardation | 1.0 (standard free-fall) |
| Lyman-Werner field J₂₁ | 0 (none) |
| Redshift | 0 (no CMB contribution) |
| All other parameters | default values |

### Models

| Label | Binary | Z/Z☉ |
|---|---|---|
| prim | `prim_collapse` | 0 (zero-metal, 23 species) |
| metal Z=0 | `metal_collapse` | 0 |
| metal Z=1e-6 | `metal_collapse` | 10⁻⁶ |
| metal Z=1e-5 | `metal_collapse` | 10⁻⁵ |
| metal Z=1e-4 | `metal_collapse` | 10⁻⁴ |
| metal Z=1e-3 | `metal_collapse` | 10⁻³ |
| metal Z=1e-2 | `metal_collapse` | 10⁻² |
| metal Z=1e-1 | `metal_collapse` | 10⁻¹ |
| metal Z=1 | `metal_collapse` | 1 (solar) |
| metal Z=3 | `metal_collapse` | 3 |

The `metal Z=0` case isolates the cost of the larger network (95 species)
independent of metallicity-dependent chemistry.

Each metal case is run with both networks: the **full** 95-species network
(`metal_collapse`) and the compact **minimal** 44-species network
(`arche_collapse_metal_minimal`).  In every figure below the two variants share
a colour per metallicity and are distinguished by line/bar style — full is
solid (and a solid bar), minimal is dashed (and a hatched bar).  This lets the
extra cost of the full network be read directly against the minimal model at
matched Z.

---

## Results

### Test-run trajectories

![fig_bench_rho_T](img/bench/fig_bench_rho_T.png)

- As a visual reference for the benchmark set, this figure shows the
  `nH-T` trajectories of the metal-model test runs used below (solid = full,
  dashed = minimal).
- The benchmark timing and cost-estimate figures should be interpreted along
  these thermal tracks, not as density-independent averages.

### Simple estimate based on `t_chem`

![fig_bench_cost_ratio_summary](img/bench/fig_bench_cost_ratio_summary.png)

- Top panel: `wall_time_metal / wall_time_prim` for one chemistry call.
- Middle panel: `t_chem,metal / t_chem,prim`, which controls the relative
  number of chemistry subcycles if the outer timestep `dt` is fixed.
- Bottom panel: a simple local estimate of the chemistry cost ratio,
  `Cost_metal / Cost_prim`, obtained from the product of the two effects above.

In all three panels the dashed curves are the minimal network and the solid
curves the full network, at matched metallicity.

All three panels are smoothed on a shared `log10(nH)` grid with 0.2 dex
spacing and a moving window of ±0.5 dex.

This ratio alone is **not** the total 3D cost ratio.  To estimate the
chemistry-related computational cost in a hydrodynamical simulation, one would
need the product

`(dt / dt_chem) × (wall time per step)`

and compare it between prim and metal.  In practice, however, the actual outer
timestep `dt` depends on the hydrodynamical timescale and other timestep
criteria, so this repository by itself cannot determine the final
`prim / metal` cost ratio for a coupled 3D code.

Therefore, we use the figure above only as a **simple order-of-magnitude
estimate** of the local chemistry cost ratio.

![fig_bench_cost_ratio_integrated_bar](img/bench/fig_bench_cost_ratio_integrated_bar.png)

As a compact summary, we also show the cumulative cost ratio obtained by
integrating the local estimate over the shared `log10(nH)` range and taking
the final value for each metal model.  This bar chart is the closest proxy in
this repository to the total chemistry cost ratio `metal / prim`, although it
still excludes hydro, gravity, I/O, and any timestep criteria outside the
chemistry module.  Each metallicity shows a pair of bars — the solid bar is the
full network and the hatched bar the minimal network — so the cost saved by the
compact network is read off directly.

The cumulative value integrates the local estimate over the full benchmark
density range, up to the `n_H = 10²³ cm⁻³` ceiling.  It is 4.4–6.7 for the full
network and 0.52–0.93 for the minimal network.  The local estimate (bottom panel
above) is largest at `n_H ~ 10⁸–10¹¹ cm⁻³` and near `10¹⁴ cm⁻³`; integrating
only up to `n_H = 10²⁰ cm⁻³` gives 5.9–13.2 (full) and 0.89–2.3 (minimal).  A
value below 1 means that the minimal network's longer `t_chem` more than offsets
its per-step cost relative to prim.  The number is a property of the chosen
density range as much as of the network, and should be read together with the
`nH-T` trajectories above.

### Full vs minimal: per-step cost within each network

![fig_bench_time_full_vs_min](img/bench/fig_bench_time_full_vs_min.png)

The figures above all reference the prim network.  This last one instead
compares, **within each network**, the full and minimal variants directly:
primordial (`prim_collapse`, 23 species, vs `arche_collapse_prim_minimal`,
15 species) on the left, and metal-grain (`metal_collapse`, 95 species, vs
`arche_collapse_metal_minimal`, 44 species) on the right.

- The bar height is the **10th-percentile per-step wall time** [µs].  A low
  percentile is used deliberately: it estimates the per-call compute cost with
  the least scheduler/load contention.  The run-wide mean or median is far
  noisier — a single contended segment can shift it by an order of magnitude,
  occasionally even making the smaller network appear slower.
- `×N` above each pair is the full/minimal speed-up.  The minimal network is
  ≈2.3× faster for primordial (42 µs vs 18.5 µs) and ≈2.9–3.9× faster for
  metal-grain (330–450 µs vs 100–125 µs).  The minimal per-step time is nearly
  flat across metallicity, while the full one rises with Z (334 µs at Z = 0 to
  449 µs at Z = 3), so the speed-up is largest at high metallicity.
- The speed-up is below the naive species-count ratio (e.g. `(95/44)³`): the
  implicit solve is only one part of `chem_full_step`, and the per-call fixed
  costs (rate evaluation, cooling, table lookups) are shared.

---

## Notes

- Timing is measured around the `chem_full_step` call only; I/O and HDF5
  output are excluded.
- The `metal Z=0` model shows the overhead of the 95-species network
  itself, irrespective of metal chemistry.
- At high density (nH ≳ 10¹⁰ cm⁻³), per-step time can increase due to
  grain surface reactions and additional line-cooling terms activated in
  the metal network.
- The `t_chem`-based cost figures are simplified diagnostics, not a substitute
  for timestep control implemented in a specific hydro code.
