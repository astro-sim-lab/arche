// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
#include <array>
#include <cmath>
#include <limits>

#include "core/species_index.h"

namespace arche {

enum class SolveFailure {
  none,
  max_iterations,
  stagnation,
  nonfinite,
  singular_jacobian,
  line_search,
  temperature_bracket,
  temperature_iterations,
  projection_rejected,
  invariants_unavailable,
  output_rejected,
  thermal_rejected,
  interval_budget,
  no_progress
};
// Number of SolveFailure values (the size of SolveStatus::rejections_by_reason).
inline constexpr int kSolveFailureCount =
    static_cast<int>(SolveFailure::no_progress) + 1;

struct SolveStatus {
  SolveFailure reason = SolveFailure::none;
  int newton_iterations = 0;  // maximum over inner solves
  int temperature_iterations =
      0;  // steps solved by the fixed-temperature fallback
  int attempted_substeps = 0;
  int accepted_substeps = 0;  // final schedule only; zero on interval failure
  int rejected_substeps = 0;
  int newton_total = 0;
  SolveFailure last_failure = SolveFailure::none;
  double completed_time = 0.0;
  double internal_step = 0.0;
  double residual_absolute =
      0.0;  // common max |y-y0-h*f(y)|, before projection
  double projection_max_relative = 0.0;
  double residual_norm = 0.0;
  double update_norm = 0.0;
  double temperature_norm = 0.0;
  // Estimated non-thermal shift of y_e at the Saha state (saha_nonthermal_shift);
  // -1 when not evaluated (the check is off, below nH_eq, or no non-thermal ionization).
  double saha_nonthermal_shift = -1.0;
  // Net charge of the input, |sum_i q_i y_i| / sum_i |q_i| y_i, and whether the input
  // was neutralised before the first internal step (solve/chemistry.h,
  // chem_integrate_adaptive).
  double input_charge_imbalance = 0.0;
  bool input_neutralized = false;
  // Rejected attempts of the call by cause (index = SolveFailure value).
  std::array<int, kSolveFailureCount> rejections_by_reason{};
  // The attempts rejected for last_failure in a row since the last accepted step, and
  // the width [s] of the first of them; 0 when the call ended with an accepted step.
  // On a call that ended with no_progress, last_failure is the cause of the narrowest
  // refused attempt, and this run shows from which width that cause refused every attempt.
  int last_failure_run = 0;
  double last_failure_run_h = 0.0;
};

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
  double zeta = 0.0;   // effective CR ionization rate [s^-1]
  double Nc_H = 0.0;   // H  column density [cm^-2]
  double Nc_H2 = 0.0;  // H2 column density [cm^-2]
  double Nc_HD = 0.0;  // HD column density [cm^-2]
  // Column mass of the continuum shielding length, rho*L [g cm^-2]; 0 (default) =
  // optically thin. chem_full_step derives the continuum depth from it at every state it
  // evaluates, tau = (k_gas + k_gr) * cnt_column, esc = tau > 1 ? 1/tau^2 : 1, solving
  // the grain temperature and tau together (cooling/grain.h grtemp_consistent). A depth
  // built by the caller from a previous step's opacities lags the state by one step and
  // can drive a period-2 oscillation of the grain temperature, so the depth is not an
  // input.
  double cnt_column = 0.0;
  double J_LW21 = 0.0;  // Lyman-Werner intensity [10^-21 erg/s/cm^2/Hz/sr]
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
// All rates in [erg g^-1 s^-1], averaged over the COMPLETE requested interval.
// Only accepted internal steps contribute; the host applies -dt*Lambda_net.
// Opacity and grain temperature are final accepted values, never averages.
// Density, input gas temperature and shielding stay fixed within the interval.
// No fluid energy or compression work is updated by this API.
// Metal-only fields are 0 for the zero_metal network.
// ---------------------------------------------------------------------------
struct ChemFullRates {
  SolveStatus status;
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
  // Opacities, continuum depth and grain temperature of the interval's end state.
  double k_gas = 0.0;    // gas opacity [cm^2/g]
  double k_gr = 0.0;     // grain opacity × Z_metal [cm^2/g]  (metal only)
  double tau_cnt = 0.0;  // continuum optical depth (ChemShielding::cnt_column)
  double esc_cnt = 1.0;  // continuum escape fraction
  // Grain temperature [K] (metal only; the caller's ChemParams::T_gr_K for zero_metal).
  // The caller feeds this back into ChemParams::T_gr_K as the next call's warm-start
  // seed of the grain-temperature solves.
  double T_gr_K = 0.0;
  bool solver_failed = false;
  // true when every accepted step of the interval was projected onto the element and
  // charge invariants (solve/conservation.h) within the tolerance weights; false when
  // the interval failed (the cell is then returned as it came in).
  bool conservation_projected = false;
};

// ---------------------------------------------------------------------------
// ChemParams  — external parameters that are read-only for the kernel.
//
// Contains both zero_metal and metal_grain parameters.
// For zero_metal runs: T_gr_K = 0, Z_metal = 0 (unused).
// ---------------------------------------------------------------------------
struct ChemParams {
  double internal_max_step =
      std::numeric_limits<double>::infinity();  // internal step bound [s]
  // Error-controlled integrator (solve/integrator.h). Per-component weights
  // atol + rtol*|y_i| for the abundances (per H nucleus) and rtol_T*T for the
  // temperature bound both the Newton correction and the local error of a step.
  // atol is the same for every species: no species is singled out as unimportant.
  // A smaller atol makes the first-order steps resolve trace tails far below any
  // abundance that feeds back on the thermal evolution, at a higher cost; a larger one
  // stops following a species that grows from below atol (e.g. the grain ices of a
  // low-metallicity collapse), which shifts the onset of its freeze-out.
  double rtol = 1.0e-4;
  double atol = 1.0e-14;
  double rtol_T = 1.0e-4;
  int internal_max_steps =
      100000;         // attempted steps per call (accepted + rejected)
  double zeta = 0.0;  // cosmic-ray ionization rate [s^-1]  (pre-attenuated)
  double zeta_X =
      0.0;  // X-ray secondary ionization rate [s^-1] (not CR-shielded)
  // Above saha_threshold_nH, a step is taken in the Saha regime only when there is no
  // non-thermal ionization (zeta = zeta_X = 0). With a positive tolerance it also is when
  // the non-thermal ionization, as a perturbation of the Saha state, shifts y_e by at
  // most this fraction (solve/chemreact.h, saha_nonthermal_shift; reported in
  // SolveStatus::saha_nonthermal_shift). 0 = never. Default 1e-3.
  double saha_nonthermal_tol = 1.0e-3;
  // High-density Saha regime (solve/integrator.h saha_step): for nH > saha_threshold_nH when
  // use_saha is true (and, with zeta or zeta_X > 0, only as saha_nonthermal_tol allows).
  // use_saha = false solves the kinetic network by Newton at every density. The
  // default threshold is nH_eq = 1e18. The Saha solvers carry
  // approximations made for nH > 1e18 (e.g. the primordial one keeps He+ out of the
  // charge balance), so a threshold below 1e18 is outside what they were built for.
  bool use_saha = true;
  double saha_threshold_nH = 1.0e18;  // [cm^-3]
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
// The default Saha threshold (ChemParams::saha_threshold_nH).
constexpr double nH_eq = 1.0e18;

}  // namespace numerics

namespace model {
constexpr double cr_photo_albedo = 0.6;
}  // namespace model

}  // namespace arche
