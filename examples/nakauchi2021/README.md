# Example: Nakauchi et al. 2021

Reproduces **Figures 10 and 11** of

> Nakauchi, Omukai & Susa (2021) MNRAS, 502, 3394
> [[ADS](https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.3394N/abstract)]

with the arche `metal_grain` network, comparing the **reduced metal
Minimal** model (44 species) against the **full Nakauchi2021** model (95
species) over a one-zone gravitational collapse to nₕ = 10²³ cm⁻³. Both counts
include the metals locked in the dust (C_s, O_s, Mg_s; the full model also
Na_s, K_s) and the grain-volume tracer Gr_v.

## Grid — Figure 10

`run_fig10.sh` runs 7 metallicities × 3 cosmic-ray ionization rates × {reduced,
full} = 42 collapses:

| Parameter | Values |
|---|---|
| Z/Z☉ | 1, 10⁻¹, 10⁻², 10⁻³, 10⁻⁴, 10⁻⁵, 10⁻⁶ |
| ζ_ion | 0, 10⁻¹⁷, 10⁻¹⁵ s⁻¹ |
| network | reduced (`arche_collapse_metal_minimal`) / full (`metal_collapse`) |

## Usage — Figure 10

Run from the **project root**:

```bash
# Step 1: build both binaries + run the 42-case grid
bash examples/nakauchi2021/run_fig10.sh
#   (re-run with --no-build to reuse already-built binaries)

# Step 2: plot
python3 examples/nakauchi2021/plot_fig10.py
```

Output:
- data: `results/nakauchi2021/collapse_CR<ζ>_Z<z>[_min].h5`
- figure: `results/nakauchi2021/fig10_reproduction.png`

## Panels (paper layout)

| Row | Content | y-range | nₕ range |
|---|---|---|---|
| (a) | T_gas | 10⁰–10⁴ K | –10²³ |
| (b) | y(H₂), y(HD) | 10⁻⁸–10⁰ | –10²³ |
| (c) | y(e⁻) | 10⁻¹⁴–10⁻⁴ | –10²³ |
| (d) | y(gr⁻) | 10⁻¹⁸–10⁻⁸ | –10²³ |

Columns are the three ζ_ion values; colour = metallicity (black = 1 Z☉ → orange
= 10⁻⁶ Z☉, matching the paper's legend). To make the near-perfect reduced/full
agreement readable at a glance, the **full** model is drawn as a **thick dashed
line** and the **reduced** model as a **thin solid line** on top:
where they agree the solid line tracks the dashes; where they diverge
(e.g. e⁻ at nₕ ≳ 10¹⁵) the solid line visibly leaves the dashes.

## Agreement with the paper

- **(a) T**: lower Z → hotter (weaker dust cooling). Reduced/full overlap at all
  Z and ζ.
- **(b) H₂/HD**: reduced and full agree across the grid.
- **(c) e⁻**: reduced ≈ full over most of the density range. For nₕ ≳ 10¹⁵ the
  reduced model falls below the full one — the missing K/Na collisional
  ionization channels (the paper's key point).
- **(d) gr⁻**: grain negative charge grows with metallicity; reduced/full agree.

## Output

![Nakauchi2021 Fig.10 reproduction](../../docs/img/examples/nakauchi2021/fig10_reproduction.png)

---

## Companion grid — CR = 0 metallicity scan

`run_cr0_zscan.sh` runs the **full** `metal_grain` network at ζ_ion = 0 over the
same 7 metallicities, without the reduced/full comparison. It isolates the
metallicity dependence when no cosmic-ray ionization is present.

| Parameter | Values |
|---|---|
| Z/Z☉ | 1, 10⁻¹, 10⁻², 10⁻³, 10⁻⁴, 10⁻⁵, 10⁻⁶ |
| ζ_ion | 0 s⁻¹ |
| network | full (`metal_collapse`) |

Run from the **project root**:

```bash
# Step 1: build the binary + run the 7-case grid
bash examples/nakauchi2021/run_cr0_zscan.sh
#   (re-run with --no-build to reuse an already-built binary)

# Step 2: plot
python3 examples/nakauchi2021/plot_cr0_zscan.py
```

Output:
- data: `results/cr0_zscan/metal_grain/collapse_CR0_Z<z>.h5`
- figure: `results/cr0_zscan/fig_cr0_zscan.png`

Panels: T, y(e⁻), y(H₂) + y(HD) (H₂ solid, HD dashed), and y(gr⁻), each against
nₕ with one curve per metallicity.

---

## Figure 11 — magnetic Reynolds number

`run_fig11.sh` runs the **full** `metal_grain` network for the 4 × 3 cases of the
paper's Figure 11 and evaluates the magnetic Reynolds number at the Jeans scale on
each collapse track with the `magnetic_reynolds` application:

  Rm(λ_J) = v_ff λ_J / (η_O + η_A),  v_ff = λ_J / (3 t_ff)

(paper eq. 41). η_O and η_A come from `arche::model_resistivity` for the
composition, nₕ and T of each output row, on log₁₀(B/G) = −10 … 10 (401 points).
The collapse is computed without magnetic fields, so one track serves every B.

| Parameter | Values |
|---|---|
| Z/Z☉ | 10⁻⁶, 10⁻⁴, 10⁻², 1 |
| ζ_ion | 0, 10⁻¹⁷, 10⁻¹⁵ s⁻¹ |
| network | full (`metal_collapse`) |

Run from the **project root**:

```bash
# Step 1: build both binaries, run the 12 collapses and the Rm maps
bash examples/nakauchi2021/run_fig11.sh
#   --no-build     reuse already-built binaries
#   --no-collapse  reuse binaries and the collapse files already in the output
#                  directory (run_fig10.sh writes the same 12 full-network files)

# Step 2: plot
python3 examples/nakauchi2021/plot_fig11.py
```

Output:
- data: `results/nakauchi2021/collapse_CR<ζ>_Z<z>.h5` (collapse track) and
  `results/nakauchi2021/rm_CR<ζ>_Z<z>.h5` (`nH`, `T_K`, `rho`, `lambda_J`,
  `t_ff`, `B_cr`, `log10_B`, `eta_O` (per row), `eta_A` and `Rm` (row × B))
- figure: `results/nakauchi2021/fig11_reproduction.png`

Panels follow the paper: rows are Z/Z☉ = 10⁻⁶, 10⁻⁴, 10⁻², 1, columns are the
three ζ_ion values; colour is log Rm clipped to [−2, 10]; black contours are
log Rm = 0 (thick) and 2, 4, 6, 8, 10 (thin); the white curve is
B_crit = (4πG M_J ρ / λ_J)^{1/2}.

### Agreement with the paper (ζ_ion = 0)

| Z/Z☉ | Rm < 1 at every B (Ohmic), paper | here | recoupling, here | Rm(B_crit) < 1 from, paper | here |
|---|---|---|---|---|---|
| 10⁻⁶ | nₕ ≳ 10¹⁵ | 10^14.9 | 10^16.2 | — | — |
| 10⁻⁴ | ≳ 10¹⁰ | 10^9.8 | 10^16.5 | ≳ 10⁵ | 10^4.9 |
| 10⁻² | ≳ 10⁸ | 10^7.9 | 10^15.9 | ≳ 10³ | 10^2.4 |
| 1 | ≳ 10⁸ | 10^8.2 | 10^15.3 | ≳ 10¹ | 10^0.0 |

The paper places the recoupling at nₕ ~ 10¹⁵–10¹⁷ cm⁻³. With ionization sources
the dissipating band narrows or disappears as in the paper's middle and right
columns (e.g. Z = 10⁻² Z☉: 10^13.0–10^15.9 at 10⁻¹⁷ s⁻¹, 10^15.0–10^15.9 at
10⁻¹⁵ s⁻¹; none at Z = 10⁻⁶ Z☉).

![Nakauchi2021 Fig.11 reproduction](../../docs/img/examples/nakauchi2021/fig11_reproduction.png)
