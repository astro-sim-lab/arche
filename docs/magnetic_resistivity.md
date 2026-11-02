# Ohmic and ambipolar magnetic diffusivities

arche computes the Ohmic and ambipolar magnetic diffusivities η_O and η_A
[cm² s⁻¹] of a cell from its abundances, n_H, the gas temperature T and the
field strength B. The calculation works for all four registered models.

| Entry point | Header |
|---|---|
| `arche::model_resistivity(cell, T_K, B_G)` → `MagneticDiffusivities{eta_O, eta_A}` | `api/arche_api.h` |
| `arche_model_resistivity(cell, T_K, B_G)` → `ArcheMagneticDiffusivities` | `api/arche_capi.h` |
| `arche::magnetic::diffusivities<Model, PairRates>(y, nH, T_K, B_G)` (kernel) | `magnetic/resistivity.h` |
| `arche::magnetic::mt_rate(MtPair, T_K)` (rate tables) | `magnetic/momentum_transfer_rates.h` |

The Hall diffusivity is not computed. The application `magnetic_reynolds`
evaluates the magnetic Reynolds number over a collapse track; see
[`examples/nakauchi2021`](../examples/nakauchi2021/README.md) (Figure 11).

---

## 1. Definitions

Charged particles ν (charge q_ν, mass m_ν, number density n_ν) collide with the
neutrals H, H₂ and He (Wardle & Ng 1999; Nakauchi, Omukai & Susa 2021, eqs. 37–39):

```text
eta_O   = c^2 / (4 pi sigma_O)
eta_A   = c^2 / (4 pi) [sigma_P / (sigma_P^2 + sigma_H^2) - 1 / sigma_O]
sigma_O = sum g_nu,                 g_nu = n_nu q_nu^2 tau_nu / m_nu
sigma_P = sum p_nu,                 p_nu = g_nu / (1 + beta_nu^2)
sigma_H = (c/B) sum q_nu n_nu / (1 + beta_nu^2)
beta_nu = q_nu B tau_nu / (m_nu c)                       (signed Hall parameter)
1/tau_nu = sum_n rho_n <sigma v>_nu,n / (m_nu + m_n),    n in {H, H2, He}
```

### Numerical form

Two exact rearrangements remove the cancellations of these expressions.

- **σ_H.** Charge neutrality gives Σ h_ν = 0 for h_ν = (c/B) q_ν n_ν, so
  σ_H = Σ h_ν (r_ν − w) for any constant w, with r_ν = 1/(1+β_ν²).
  - w = 0 is the expression above. Its terms cancel when every |β| ≪ 1, for
    example at high density.
  - w = 1 gives −Σ p_ν β_ν. Its terms cancel when every |β| ≫ 1, for example at
    low density or in strong fields.
  - arche takes w as the median of r_ν weighted by |h_ν|. This minimises the
    rounding bound Σ |h_ν| |r_ν − w|.
  - r_ν − w is formed as (1 − w) − (1 − r_ν) when w > 1/2, with
    1 − r = β²/(1+β²) evaluated directly, so β² below the double-precision
    epsilon is not lost.
- **Ambipolar numerator.** σ_O σ_P − σ_P² − σ_H² equals the Lagrange-identity
  sum ½ Σ_ij p_i p_j (β_i − β_j)² = P² Var_p(β), where P = σ_P and Var_p is the
  p-weighted variance of β with mean ⟨β⟩_p = −σ_H/P. Hence

  ```text
  eta_A = eta_O Var_p(beta) / (1 + <beta>_p^2)
  ```

  Var_p is summed about ⟨β⟩_p, so every term is non-negative.

### Inputs and special values

- Abundances below zero count as zero.
- Both values are quiet NaN unless n_H and T are positive and finite and B is
  non-negative and finite. They are also NaN when charged particles are present
  but none of H, H₂ and He is, since there is then no collision time.
- With no charged particle (σ_O = 0), η_O = +∞, and η_A = +∞ for B > 0.
- η_A = 0 for B = 0.

---

## 2. Charged particles

**Gas-phase species.** Every species of the model with a non-zero charge
(`core/species_composition.h`), with its catalog mass (`core/species_catalog.h`)
times m_p.

**Grains.** The charge states Gr⁺, Gr²⁺, Gr⁻ and Gr²⁻
(`composition::grain_count == 1`) are spread over the MRN size distribution of
the metal_grain network (`core/species_index.h`):

- n(a) ∝ a^−3.5 on 0.005–1 µm and a^−5.5 on 1–5 µm, in 15 + 5 bins equally
  spaced in ln a.
- Each bin is represented by its centre radius a_k, with mass
  (4/3)π a_k³ ρ_gr and ρ_gr = 3 g cm⁻³.
- The fraction of grains in a bin is the exact integral of n(a) over the bin
  divided by the integral over the whole range, so the fractions sum to 1.
  The fraction is the same for every charge state.
- The size distribution is fixed. It does not follow the grain evaporation that
  the metal_grain models track with `Gr_v`.

---

## 3. Momentum-transfer rate coefficients

| Pair | ⟨σv⟩ |
|---|---|
| e–H, e–H₂, e–He, H⁺–H, H⁺–H₂, H⁺–He, H₃⁺–H₂, HCO⁺–H₂, C⁺–H | tables in `magnetic/momentum_transfer_rates.h` (below) |
| other ion–neutral pairs | Langevin, 2.41π\|Z\| e √(α_n/μ) (Osterbrock 1961); α_H 0.667, α_H₂ 0.804, α_He 0.207 Å³; μ the reduced mass |
| grain–neutral | (4/3)π a² δ √(8kT/(π m_n)), δ = 1.3 (Pinto & Galli 2008, eq. 25) |

The Langevin rate scales as |Z| because the polarisation potential of a charge
Ze scales as Z².

### Definition of the tables

The tabulated rates are the zero-drift limit of Pinto & Galli (2008, A&A 484, 17,
eq. 8):

```text
<sigma v> = (4/3) sqrt(8 k T / (pi mu)) (1/2) (k T)^-3
            Int_0^inf sigma_mt(E) E^2 exp(-E / k T) dE
```

Here σ_mt is the momentum-transfer cross section, E the centre-of-mass energy
and μ the reduced mass.

- **Constants:** k = 1.380649×10⁻¹⁶ erg K⁻¹; 1 eV = 1.602176634×10⁻¹² erg;
  m_e = 9.109383714×10⁻²⁸ g; atomic mass constant 1.66053906892×10⁻²⁴ g;
  a₀ = 5.29177210544×10⁻⁹ cm.
- **Masses** (amu, `species_catalog.h`): H 1.00783, H₂ 2.01565, He 4.0026,
  H⁺ 1.00728, H₃⁺ 3.02293, HCO⁺ 29.00219, C⁺ 11.99945.

| Pair | Cross section | Range | Interpolation; outside the range |
|---|---|---|---|
| e–H | Dalgarno, Yan & Liu 1999, ApJS 125, 237 (Fig. 1b), as the solid curve of Pinto & Galli 2008 Fig. 14 | 10⁻³–10³ eV | log–log; constant below (zero-energy limit), power law of the last interval above |
| e–H₂ | IST-Lisbon database (LXCat), ELASTIC, set of Baghel, Marques & Alves 2026, Plasma Sources Sci. Technol. 35, 085023 | 0–10⁴ eV | linear; zero above the last energy (as the database prescribes) |
| e–He | IST-Lisbon database (LXCat), ELASTIC, set of Alves et al. 2013, J. Phys. D 46, 334002 and Santos et al. 2014, J. Phys. D 47, 265201 | 0–10³ eV | as e–H₂ |
| H⁺–H | Schultz, Krstić, Lee & Raymond 2008, ApJ 678, 950 (Fig. 8 data; below 10² eV Glassgold, Krstić & Schultz 2005, ApJ 621, 808) | 10⁻⁴–10⁶ eV (c.m.) | log–log; σ ∝ E^−1/2 below, joined to the first point |
| H⁺–H₂ | Krstić & Schultz 1999 (APID 8, 1; Phys. Rev. A 60, 2118), as the solid curve of Pinto & Galli 2008 Fig. 6 | 0.10–76 eV | log–log; 1.05×10⁻¹⁵ (E/eV)^−0.531 cm² below (Pinto & Galli 2008), power law of the last interval above |
| H⁺–He | Krstić & Schultz 1999 (APID 8, 1; J. Phys. B 32, 2415), as the solid curve of Pinto & Galli 2008 Fig. 16 | 0.099–90 eV | log–log; Langevin cross section 2.41π √(α e²/(2E)) below, power law of the last interval above |
| H₃⁺–H₂ | Phelps 1990, J. Phys. Chem. Ref. Data 19, 653, as the solid curve of Pinto & Galli 2008 Fig. 4 (c.m. energy) | 0.040–4×10³ eV | log–log; 1.62×10⁻¹⁵ (E/eV)^−0.507 cm² below (Pinto & Galli 2008), power law of the last interval above |
| HCO⁺–H₂ | Flower 2000, MNRAS 313, L19 (erratum 315, 431), the computed points of Pinto & Galli 2008 Fig. 2 | 8.6×10⁻⁴–0.60 eV | σ linear in log E; σ ∝ E^−1/2 below; above, the power law fitted (least squares, log–log) to the points of the last decade |
| C⁺–H | Flower & Pineau des Forêts 1995, MNRAS 275, 1049, as the solid curve of Pinto & Galli 2008 Fig. 10 | 10⁻⁵–10 eV | log–log; σ ∝ E^−1/2 below, power law of the last interval above |

Notes on the sources:

- The curves of Pinto & Galli (2008) are read from the vector paths of the
  figures in arXiv:0804.2856v1.
- The Schultz et al. (2008) table gives energies to 3 digits; rows that round
  to the same energy are averaged.
- Atomic units of area are a₀².

How the tables are built and evaluated:

- **Integration:** trapezoidal rule in ln E from 10⁻⁷ kT to 80 kT, 40001
  points. The difference from an integration 5 times finer is ≤ 7.3×10⁻⁷.
- **Grid:** log₁₀(T/K) = 0, 0.01, …, 6 (601 points); values written to 8
  significant digits.
- **Run time:** `mt_rate` interpolates linearly in (log T, log ⟨σv⟩). The error
  against direct integration at the interval midpoints is ≤ 5.5×10⁻⁵.
  - Below 1 K: ⟨σv⟩ ∝ T^1/2 for electrons, constant for ions.
  - Above 10⁶ K: the power law of the last interval.
  - T ≤ 0 or non-finite T gives NaN.

arche carries only the rate tables derived from these cross sections. The
cross-section data themselves (in particular the LXCat sets, whose terms do not
allow redistribution by third parties) are not included.

### Comparison with Pinto & Galli (2008)

Ratio of the tables to the fits of Pinto & Galli (2008, Table 1), on log₁₀ T =
1–4 (HCO⁺–H₂: 1–3). Pinto & Galli state that the fits are accurate to about
10%. The larger ratios fall where the fits leave the numerical rates they
approximate, e.g. H⁺–H below the 24 K lower limit of its fit, and C⁺–H and
e–H₂ near 10⁴ K.

| Pair | min | max |
|---|---|---|
| e–H | 0.990 | 1.073 |
| e–H₂ | 0.634 | 1.014 |
| e–He | 0.982 | 1.325 |
| H⁺–H | 1.100 | 2.046 |
| H⁺–H₂ | 1.014 | 1.845 |
| H⁺–He | 0.959 | 1.092 |
| H₃⁺–H₂ | 0.774 | 1.061 |
| HCO⁺–H₂ | 0.983 | 1.045 |
| C⁺–H | 0.609 | 1.012 |

Ratio of the tables to the rate curves that Pinto & Galli (2008) computed from
their cross sections, read from the vector paths of their figures:

| Figure | Pair | table / figure |
|---|---|---|
| Fig. 15 | e–H | 0.996–1.001 at 10²–3×10⁴ K |
| Fig. 13 | H⁺–H | 0.894–1.002 at 10²–10⁴ K |
| Fig. 17 | H⁺–He | within 0.3% at 10–3×10⁴ K |
| Fig. 3 | HCO⁺–H₂ | 0.998–1.015 at 30–10³ K |
| Fig. 11 | C⁺–H | 0.999–1.004 at 10²–10⁴ K |
| Fig. 5 | H₃⁺–H₂ | 0.986 at 10² K; 0.76–0.79 at 10³–10⁴ K |
| Fig. 7 | H⁺–H₂ | 1.02–1.08 at 10³–10⁴ K; 1.32 at 10² K, 1.86 at 10 K |
| Fig. 19 | e–He | 1.19–1.33 at 10³–10⁴ K |

Figs. 5, 7 (at and below 10² K) and 19 do not follow Pinto & Galli's own cross
sections; there the tables follow the cross sections.

- The cross section of Fig. 4 (H₃⁺–H₂) stays below 0.97 times the Langevin one,
  so its rate cannot exceed the Langevin rate. Fig. 5 is 1.10–1.21 times the
  Langevin rate at 10³–10⁴ K.
- Integrating the cross section of Fig. 6 (H⁺–H₂) with Pinto & Galli's
  extension below 0.1 eV gives 2.17×10⁻⁹ cm³ s⁻¹ at 10 K; Fig. 7 shows
  1.17×10⁻⁹.
- The curve of Fig. 19 (e–He) equals 0.428 T^1/2, which corresponds to a
  constant cross section of 5.2×10⁻¹⁶ cm²; the cross sections peak at
  7×10⁻¹⁶ cm² near 2 eV.

---

## 4. Differences from Pinto & Galli (2008)

arche uses the definitions of Nakauchi, Omukai & Susa (2021, eqs. 37–39; §1) and
differs from the rate prescriptions of Pinto & Galli (2008) in two choices:

1. **⟨σv⟩** for the nine pairs of §3 comes from the cross sections, not from the
   Table 1 fits. Several of those fits become negative above their range:
   HCO⁺–H₂ at 3.6×10³ K, H⁺–He at 1.5×10⁴ K, e–H at 3.5×10⁴ K.
2. **Langevin rate** scales as |Z| (Pinto & Galli 2008 eqs. A.1–A.3 use |Z|^1/2).
   This affects He²⁺, Li²⁺ and Li³⁺ only.
