// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
#include <array>
#include <cmath>
#include <limits>

#include "core/species_index.h"

namespace arche {

// ---------------------------------------------------------------------------
// ChemState  — single-cell state updated by the chemistry kernel.
//
// Holds only what the kernel needs to read/write per time step.
// Thermodynamic variables that can be derived (rho, p, e) are kept
// in CollapseState (app layer) to avoid redundancy.
// ---------------------------------------------------------------------------
template <int N_sp>
struct ChemState {
  std::array<double, N_sp>
      y{};           // species abundances (dimensionless, relative to nH)
  double nH = 0.0;   // H number density [cm^-3]
  double T_K = 0.0;  // gas temperature [K]
  double mu = 1.0;   // mean molecular weight [m_p]
  double gamma = 5.0 / 3.0;  // adiabatic index
};

// Convenience aliases
using ChemStateZM = ChemState<zero_metal::N_sp>;
using ChemStateMG = ChemState<metal_grain::N_sp>;

// ---------------------------------------------------------------------------
// ChemRates  — scalar rates returned by the kernel to the app layer.
//
// The full var[2*N_react] reaction rate array is NOT exposed.
// CR heating (var[131:136]) is summed inside the kernel.
// ---------------------------------------------------------------------------
struct ChemRates {
  double Lambda_chem = 0.0;  // net chemistry cooling rate [erg g^-1 s^-1]
  double Gamma_CR = 0.0;     // CR heating rate            [erg g^-1 s^-1]
  double Gamma_X = 0.0;      // X-ray heating rate         [erg g^-1 s^-1]
  // Did the element/charge conservation projection (solve/conservation.h) run
  // to completion on this step?  false has FOUR causes:
  //   1. the model registers no invariant rows (n_invariants == 0).
  //      All four shipped models register rows (5 / 5 / 8 / 10), so case 1
  //      does not fire for them; it is the path for any model whose species
  //      are not all registered in core/species_composition.h, or whose
  //      reaction tables do not balance under them (fail-closed).
  //   2. the step did not converge, so no projection was attempted.
  //   3. a weighted invariant matrix that is not positive definite.
  //   4. a repair larger than conservation::kMaxRelShift.
  // In every case the abundances are exactly what the solver produced.  Only
  // 3 and 4 are anomalies; 1 and 2 are ordinary.  Do NOT read this flag as
  // "something failed" — read it as "conservation was not enforced on this
  // step", which is all it asserts.
  bool conservation_projected = false;
};

// ---------------------------------------------------------------------------
// ChemShielding  — shielding environment pre-computed by the calling code.
//
// For a 1-D collapse app: derived from local Jeans / shielding lengths.
// For a 3-D fluid code:   from RT, TreeCol, or Sobolev approximation.
//
// zeta already carries the full attenuation factor; the chemistry kernel
// uses this value directly without further reduction.
//
// Metal-only column densities (Nc_CO, Nc_OH, ...) are zero-initialised
// and ignored by the zero_metal network.
// ---------------------------------------------------------------------------
struct ChemShielding {
  double zeta = 0.0;     // effective CR ionization rate [s^-1]
  double Nc_H = 0.0;     // H  column density [cm^-2]
  double Nc_H2 = 0.0;    // H2 column density [cm^-2]
  double Nc_HD = 0.0;    // HD column density [cm^-2]
  double tau_cnt = 0.0;  // continuum optical depth
  double esc_cnt = 1.0;  // continuum escape fraction  (1 = optically thin)
  double J_LW21 = 0.0;   // Lyman-Werner intensity [10^-21 erg/s/cm^2/Hz/sr]
                         // Photodissociates H2/HD and photodetaches H-.
                         // 0.0 = no LW field (default).
  // metal_grain only (leave at 0 for zero_metal network):
  double Nc_CO = 0.0;
  double Nc_OH = 0.0;
  double Nc_H2O = 0.0;
  double Nc_CII = 0.0;
  double Nc_CI = 0.0;
  double Nc_OI = 0.0;
  double zeta_X =
      0.0;  // HI primary photoionization rate [s^-1]  (0 = no X-ray)
  double E_X_eV = 300.0;  // representative X-ray photon energy [eV]
};

// ---------------------------------------------------------------------------
// ChemFullRates  — complete cooling/heating breakdown returned by
//                  chem_full_step().
//
// All rates in [erg g^-1 s^-1].
// Metal-only fields are 0 for the zero_metal network.
// ---------------------------------------------------------------------------
struct ChemFullRates {
  // Aggregate rates
  double Lambda_net = 0.0;   // net cooling  = line + cnt + ch − CR
  double Lambda_line = 0.0;  // total line cooling
  double Lambda_cnt = 0.0;   // total continuum cooling  (grain + gas)
  double Lambda_chem = 0.0;  // chemistry cooling
  double Gamma_CR = 0.0;     // CR ionization heating
  double Gamma_X = 0.0;      // X-ray heating rate [erg g^-1 s^-1]
  // Per-line rates
  double Lambda_H2 = 0.0;
  double Lambda_HD = 0.0;
  double Lambda_Lya = 0.0;
  // metal_grain only (0 for zero_metal):
  double Lambda_CO = 0.0;
  double Lambda_OH = 0.0;
  double Lambda_H2O = 0.0;
  double Lambda_CII = 0.0;
  double Lambda_CI = 0.0;
  double Lambda_OI = 0.0;
  double Lambda_gr = 0.0;   // grain continuum cooling
  double Lambda_gas = 0.0;  // gas (ff+CIA) continuum cooling
  // Opacities — needed to compute tau_cnt in the next step
  double k_gas = 0.0;  // gas opacity [cm^2/g]
  double k_gr = 0.0;   // grain opacity × Z_metal [cm^2/g]  (metal only)
  // Grain temperature solved this step [K] (metal only; 0 for zero_metal).
  // The caller feeds this back into ChemParams::T_gr_K as the next step's
  // warm-start seed (see chem_full_step / cnt_cool_metal).
  double T_gr_K = 0.0;
  bool solver_failed = false;
  // See ChemRates::conservation_projected.
  bool conservation_projected = false;
};

// ---------------------------------------------------------------------------
// ChemParams  — external parameters that are read-only for the kernel.
//
// Contains both zero_metal and metal_grain parameters.
// For zero_metal runs: T_gr_K = 0, Z_metal = 0 (unused).
// ---------------------------------------------------------------------------
struct ChemParams {
  double zeta = 0.0;  // cosmic-ray ionization rate [s^-1]  (pre-attenuated)
  double zeta_X =
      0.0;  // X-ray secondary ionization rate [s^-1] (not CR-shielded)
  // NaN by default so missing caller assignment is caught explicitly.
  double T_rad = std::numeric_limits<double>::quiet_NaN();  // CMB radiation
                                                            // temperature [K]
  double T_gr_K = 0.0;   // grain temperature [K]  (metal_grain only; input
                         // warm-start seed — the solved value is returned in
                         // ChemFullRates::T_gr_K, never written back here)
  double Z_metal = 0.0;  // metallicity [Z_sun]    (metal_grain only)
  double T_cr_desorp = 70.0;  // effective CR desorption spike temperature [K]

  // critical densities for LTE/low-density interpolation (filled by kernel)
  double H = 0.0;   // H abundance (copy of y[H]; drives the LTE
                    // critical-density mix)
  double H2 = 0.0;  // H2 abundance (2*y[H2])
  double He = 0.0;  // He abundance (copy of y[He])

  // metal_grain extra: UV radiation field strengths (J_H2, J_H2O, J_tot)
  double J_H2 = 0.0;
  double J_H2O = 0.0;
  double J_tot = 0.0;
};

// ---------------------------------------------------------------------------
// Physical constants  (CGS)
//
// CODATA 2022 values converted to CGS, written to 10 significant digits.
// Constants that are exactly defined in the SI (k_B, h, e, and the values
// derived from them) terminate earlier and carry no uncertainty; G is only
// known to 6 significant digits.
// ---------------------------------------------------------------------------
namespace phys {
constexpr double k_B = 1.380649e-16;      // Boltzmann [erg/K] (exact)
constexpr double h_P = 6.62607015e-27;    // Planck [erg·s] (exact)
constexpr double hbar = 1.054571817e-27;  // reduced Planck h/2π [erg·s]
constexpr double m_p = 1.672621926e-24;   // proton mass [g]
constexpr double m_e = 9.109383714e-28;   // electron mass [g]
constexpr double qe = 4.803204713e-10;    // elementary charge [statC]
constexpr double pi = 3.141592653589793;  // π (nearest double)
constexpr double G = 6.67430e-8;          // gravitational [cm^3 g^-1 s^-2]
constexpr double sigma_B =
    5.670374419e-5;  // Stefan-Boltzmann [erg cm^-2 s^-1 K^-4]
constexpr double eV_to_erg = 1.602176634e-12;  // [erg/eV] (exact)
constexpr double k_B_eV = 8.617333262e-5;      // Boltzmann [eV/K]
}  // namespace phys

namespace numerics {
constexpr double eps_gaussj = 1.0e-13;
constexpr double eps_it_prim = 1.0e-13;
constexpr double eps_it_metal = 1.0e-5;
constexpr double eps_y = 1.0e-10;
constexpr double nH_eq = 1.0e18;

}  // namespace numerics

namespace model {
constexpr double cr_photo_albedo = 0.6;
}  // namespace model

}  // namespace arche
