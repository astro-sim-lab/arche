// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// integrator.h — error-controlled backward Euler for the species AND the internal
// energy (the model's equations are unchanged; only their discretisation lives here).
//
// One step solves, for x = (y, tau = T/T_in),
//   F_i = y_i - y0_i - h r_i(y, T)                                        (i < N)
//   F_N = [ e(y,T) - e0 + h (L_chem + (L_rad(x0) + L_rad(x))/2 - Gamma_CR) ] / e0
// with one damped Newton (Deuflhard's natural monotonicity test: the simplified Newton
// correction at the trial point, measured with the current Jacobian, must shrink; a
// residual-norm line search stalls at nH >~ 1e17, where three-body rates make the
// species residual strongly nonlinear in y). e(y,T) is the EOS specific internal energy
// (thermo::specific_internal_energy), so a host inverting e0 - h*Lambda_net with the same
// EOS recovers T. L_rad (line + continuum) enters as the mean of its start- and end-state
// values: at the highest densities the fully implicit and the explicit forms err by a
// couple of percent in opposite directions, the mean by a few tenths of one.
// Everything the right-hand side depends on is evaluated at the iterate: the rate
// coefficients' LTE mix (p_loc), the UV-shielding inputs of the grain models and the
// mean molecular weight come from the iterate's y; the continuum channel from its T.
// The continuum depth follows the iterate's T through the column mass
// ChemShielding::cnt_column (tau = k_gas(T) * column for gas-only models; for grain
// models tau and the grain temperature are solved together, grtemp_consistent); a zero
// column is optically thin. The grain temperature T_gr(T) -- which the grain rate
// coefficients and the grain continuum read -- depends on T and not on y, so the
// Jacobian's temperature column carries it. The metal fine-structure lines' escape
// probabilities are solved at the iterate (cooling/cooling_molecular.h), each evaluation
// warm-started from the escape state at the step's start, so the residual is a function
// of the iterate only.
//
// Convergence: the weighted RMS of the Newton correction,
//   w = sqrt( (sum_i (dy_i/(atol + rtol|y_i|))^2 + (dT/(rtol_T T))^2) / (N+1) ),
// is at most kNewtonKappa. There is no absolute threshold: an absolute one cannot
// resolve abundances spread over 50 decades.
//
// Error control: step doubling. Each attempt solves one step of width h and two of h/2
// from the same start; the pair's end state is kept and its local error is taken as the
// weighted RMS of the difference (backward Euler's local error is C h^2, so the pair's
// error is the difference). Both paths are L-stable and land on the quasi-steady state of
// the fast species, so the estimate stays small where the true error is small -- which
// the Taylor-based (filtered) estimate does not do: at high density it exceeds the true
// error by more than 1e7. h_new = h * clamp(0.9/sqrt(err), .2, 5).
//
// Fallback: where the coupled Newton fails -- in particular across a temperature at
// which the model's right-hand side jumps (rate laws gated at 1000 K, fits switched
// piecewise), so that the coupled system has no root -- the step is solved with the
// species at a fixed trial temperature and the energy equation as a bracketed root in T
// (fixed_T_step). The same step-doubling estimate applies to it. Where the bracket closes
// on such a jump and the right-hand sides on its two sides point at it, the step ends on
// the jump with the composition that mixes the two one-sided solutions so that the energy
// equation holds (bracketed_T_step).
//
// Lyman-Werner photodissociation / photodetachment (ChemShielding::J_LW21 > 0): in the
// network steps, first-order sinks in the rates and the Jacobian (add_lw_sinks), at
// fixed specific internal energy -- their part of dy does not enter the chemistry heat.
//
// Per accepted step: the element/charge projection (conservation::project), then the
// operator-split photo channels: X-ray, and LW after a Saha step (whose composition is
// the Saha state). The split channels change y at fixed specific internal energy (plus
// X-ray heating), so the step's end temperature is re-derived from the EOS afterwards;
// the returned Lambda_net is the interval mean the caller's e update needs.
//
// Grain evaporation (grain models; solve/grain_evaporation.h): the interval's release
// enters each step as constant rates (the part of dy that is not chemistry) and as
// first-order sinks of the grain species (in the rates and the Jacobian).
//
// Saha regime (above ChemParams::saha_threshold_nH; with non-thermal ionisation, only
// as ChemParams::saha_nonthermal_tol allows): the composition is the model's Saha state at T and the step
// solves the same energy equation for T alone (saha_step).
//
// Scope: every model. Continuum: grain + gas (cnt_cool_metal_consistent) for the grain
// models, gas (cnt_cool) otherwise. Lines: H2/HD/Ly-alpha from line_cool, and for the
// metal models the metal lines and Ly-alpha from line_cool_metal.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>

#ifdef ARCHE_USE_EIGEN
#include <Eigen/Dense>
#endif

#include "cooling/cooling.h"
#include "cooling/cooling_grain.h"      // cnt_cool_metal(_consistent)
#include "cooling/cooling_molecular.h"  // line_cool_metal, EscapeState
#include "core/state.h"
#include "kinetics/rates.h"   // detail::eval_opacity
#include "solve/chemreact.h"  // compute_chemistry_cooling, saha_nonthermal_shift
#include "solve/conservation.h"
#include "solve/grain_evaporation.h"
#include "thermo/eos.h"

namespace arche {
namespace integ {

constexpr int kMaxNewton = 50;
constexpr int kMaxBacktracks = 30;
constexpr double kMaxRelativeTStep = 0.2;  // |dT|/T per Newton iteration
constexpr double kNewtonKappa = 0.1;  // weighted RMS of a converged correction
constexpr double kFdT = 1.0e-6;       // relative temperature FD step
constexpr double kFdYRel = 1.0e-7;    // relative abundance FD step (energy row)
constexpr double kFdYFloor =
    1.0e-12;  // absolute floor of the abundance FD scale
constexpr double kSafety = 0.9;
constexpr double kGrowMax = 5.0;
constexpr double kShrinkMin = 0.2;
constexpr double kShrinkFail = 0.25;  // after a failed solve
constexpr double kRootRelTol =
    1.0e-14;  // bracketed temperature root: relative width
constexpr double kJumpResidual =
    1.0e-10;  // energy residual / e0 that marks a jump

struct Tolerance {
  double rtol, atol, rtol_T;
};

template <class Model>
double eos_energy(const std::array<double, Model::N_sp>& y, double T) {
  return thermo::specific_internal_energy(
      thermo::make_eos_composition<Model>(y), T);
}

// CR ionisation heating [erg g^-1 s^-1] from the reaction-rate cache of compute_rates.
template <class Model>
double cr_heating(const std::array<double, 2 * Model::N_react>& var) {
  double rate_CR = 0.0;
  for (int i = Model::cr_heat_var_begin; i < Model::cr_heat_var_end; ++i)
    rate_CR += var[i];
  return 3.4 * phys::eV_to_erg * rate_CR /
         ((1.0 + 4.0 * abundance_ref::yHe) * phys::m_p);
}

// Dense LU with partial pivoting, kept so several right-hand sides can reuse it.
// factor() is false where a pivot is zero or not finite.  With ARCHE_USE_EIGEN it is
// Eigen's PartialPivLU (the same pivoting rule, a blocked factorisation); otherwise the
// loops below, which swap WHOLE rows (multipliers included, as LAPACK's getrf does), so
// solve() applies the row interchanges to b first, in order, and only then the
// unit-lower forward substitution.
template <int M>
class DenseLU {
 public:
  bool factor(const std::array<double, M * M>& A) {
#ifdef ARCHE_USE_EIGEN
    lu_.compute(Eigen::Map<const Eigen::Matrix<double, M, M, Eigen::RowMajor>>(
        A.data()));
    const auto& lu = lu_.matrixLU();
    for (int i = 0; i < M; ++i)
      if (!(std::abs(lu(i, i)) > 0.0) || !std::isfinite(lu(i, i))) return false;
    return true;
#else
    A_ = A;
    for (int c = 0; c < M; ++c) {
      int p = c;
      for (int r = c + 1; r < M; ++r)
        if (std::abs(A_[r * M + c]) > std::abs(A_[p * M + c])) p = r;
      piv_[c] = p;
      if (!(std::abs(A_[p * M + c]) > 0.0) || !std::isfinite(A_[p * M + c]))
        return false;
      if (p != c)
        for (int k = 0; k < M; ++k) std::swap(A_[c * M + k], A_[p * M + k]);
      for (int r = c + 1; r < M; ++r) {
        A_[r * M + c] /= A_[c * M + c];
        const double f = A_[r * M + c];
        if (f == 0.0) continue;
        for (int k = c + 1; k < M; ++k) A_[r * M + k] -= f * A_[c * M + k];
      }
    }
    return true;
#endif
  }

  // `b` is overwritten with the solution.
  void solve(std::array<double, M>& b) const {
#ifdef ARCHE_USE_EIGEN
    Eigen::Map<Eigen::Matrix<double, M, 1>> x(b.data());
    const Eigen::Matrix<double, M, 1> rhs = x;
    x = lu_.solve(rhs);
#else
    for (int c = 0; c < M; ++c)
      if (piv_[c] != c) std::swap(b[c], b[piv_[c]]);
    for (int c = 0; c < M; ++c)
      for (int r = c + 1; r < M; ++r) b[r] -= A_[r * M + c] * b[c];
    for (int r = M - 1; r >= 0; --r) {
      double s = b[r];
      for (int k = r + 1; k < M; ++k) s -= A_[r * M + k] * b[k];
      b[r] = s / A_[r * M + r];
    }
#endif
  }

 private:
#ifdef ARCHE_USE_EIGEN
  Eigen::PartialPivLU<Eigen::Matrix<double, M, M>> lu_;
#else
  std::array<double, M * M> A_{};
  std::array<int, M> piv_{};
#endif
};

// Weighted RMS of b - a (abundances and temperature).
template <int N>
double weighted_rms(const std::array<double, N>& ya, double Ta,
                    const std::array<double, N>& yb, double Tb,
                    const Tolerance& tol) {
  double s = 0.0;
  for (int i = 0; i < N; ++i) {
    const double w =
        tol.atol + tol.rtol * std::max(std::abs(ya[i]), std::abs(yb[i]));
    const double d = (yb[i] - ya[i]) / w;
    s += d * d;
  }
  const double dT =
      (Tb - Ta) / (tol.rtol_T * std::max(std::abs(Ta), std::abs(Tb)));
  s += dT * dT;
  return std::sqrt(s / (N + 1));
}

// Lyman-Werner rate coefficients [s^-1] at (nH, T): H2 and HD photodissociation with
// their self-shielding factors (H2: Wolcott-Green & Haiman 2019 eqs. 7-8, the exponent
// held at its 8000 K value above 8000 K; HD: the WG2011 form) and H- photodetachment
// (Tegmark et al. 1997; no shielding).
struct LwRates {
  double H2 = 0.0, HD = 0.0, Hm = 0.0;
};

inline LwRates lw_rates(double nH, double T_K, const ChemShielding& shield) {
  const double b5_H2 = std::sqrt(phys::k_B * T_K / phys::m_p) / 1.0e5;
  const double b5_HD = std::sqrt(phys::k_B * T_K / (1.5 * phys::m_p)) / 1.0e5;

  const double log10_nH = std::log10(nH > 0.0 ? nH : 1.0);
  // The WG19 exponent is fitted for T <= 8000 K; its linear extrapolation in log T
  // turns negative at high density from T ~ 1.1e4 K, where f_sh would exceed 1 and grow
  // with the column. Above 8000 K the exponent keeps its 8000 K value (alpha >= 0.130 at
  // every density, so f_sh <= 1), continuous in T; b5 follows the actual T.
  const double log10_TK = std::log10(T_K > 0.0 ? std::min(T_K, 8000.0) : 1.0);
  const double A1_H2 = 0.8711 * log10_TK - 1.928;
  const double A2_H2 = -0.9639 * log10_TK + 3.892;
  const double alpha_H2 = A1_H2 * std::pow(10.0, -0.2856 * log10_nH) + A2_H2;
  const double x_H2 = shield.Nc_H2 / 5.0e14;
  const double fsh_H2 =
      0.965 / std::pow(1.0 + x_H2 / b5_H2, alpha_H2) +
      0.035 / std::sqrt(1.0 + x_H2) * std::exp(-8.5e-4 * std::sqrt(1.0 + x_H2));

  const double x_HD = shield.Nc_HD / 5.0e14;
  const double fsh_HD =
      0.965 / std::pow(1.0 + x_HD / b5_HD, 2.0) +
      0.035 / std::sqrt(1.0 + x_HD) * std::exp(-8.5e-4 * std::sqrt(1.0 + x_HD));

  LwRates k;
  k.H2 = 1.38e-12 * shield.J_LW21 * fsh_H2;
  k.HD = 1.38e-12 * shield.J_LW21 * fsh_HD;
  k.Hm = 1.1e-10 * shield.J_LW21;
  return k;
}

// The Lyman-Werner channels' derivatives (rate coefficients k), added to J.
template <class Model>
void add_lw_jacobian(const LwRates& k,
                     std::array<double, Model::N_sp * Model::N_sp>& A) {
  using Sp = typename Model::Sp;
  constexpr int N = Model::N_sp;
  A[Sp::H2 * N + Sp::H2] -= k.H2;
  A[Sp::H * N + Sp::H2] += 2.0 * k.H2;
  A[Sp::HD * N + Sp::HD] -= k.HD;
  A[Sp::H * N + Sp::HD] += k.HD;
  A[Sp::D * N + Sp::HD] += k.HD;
  A[Sp::Hm * N + Sp::Hm] -= k.Hm;
  A[Sp::H * N + Sp::Hm] += k.Hm;
  A[Sp::e * N + Sp::Hm] += k.Hm;
}

// The Lyman-Werner channels as first-order sinks at y with the rate coefficients k:
// their rates are added to rp and, when J is given, their derivatives to J.
//   H2 + hv(LW) -> H + H,  HD + hv(LW) -> H + D,  H- + hv -> H + e-
template <class Model>
void add_lw_sinks(const LwRates& k, const std::array<double, Model::N_sp>& y,
                  std::array<double, Model::N_sp>& rp,
                  std::array<double, Model::N_sp * Model::N_sp>* J) {
  using Sp = typename Model::Sp;
  const double d_H2 = k.H2 * y[Sp::H2];
  const double d_HD = k.HD * y[Sp::HD];
  const double d_Hm = k.Hm * y[Sp::Hm];
  rp[Sp::H2] -= d_H2;
  rp[Sp::HD] -= d_HD;
  rp[Sp::Hm] -= d_Hm;
  rp[Sp::H] += 2.0 * d_H2 + d_HD + d_Hm;
  rp[Sp::D] += d_HD;
  rp[Sp::e] += d_Hm;
  if (J) add_lw_jacobian<Model>(k, *J);
}

// The same at (y, T).
template <class Model>
void add_lw_sinks(double nH, double T, const ChemShielding& shield,
                  const std::array<double, Model::N_sp>& y,
                  std::array<double, Model::N_sp>& rp,
                  std::array<double, Model::N_sp * Model::N_sp>* J) {
  add_lw_sinks<Model>(lw_rates(nH, T, shield), y, rp, J);
}

// The escape-probability state the metal lines carry (none for the other models).
struct NoEscape {};
template <class Model>
using EscOf = std::conditional_t<Model::has_escape, EscapeState, NoEscape>;

// Radiative (line + continuum) channel at one state, with the depth and grain
// temperature it used.
struct Radiative {
  double line = 0.0, H2 = 0.0, HD = 0.0, Lya = 0.0;
  double CO = 0.0, OH = 0.0, H2O = 0.0, CII = 0.0, CI = 0.0,
         OI = 0.0;                        // metal lines
  double gr = 0.0, gas = 0.0, cnt = 0.0;  // continuum: grain, gas, total
  double k_gas = 0.0, k_gr = 0.0;
  double T_gr =
      0.0;  // grain temperature (grain models; else the caller's T_gr_K)
  double tau = 0.0, esc = 1.0;
  // Line cooling per unit abundance of the species whose own abundance enters the line
  // cooling only linearly (as the emitter): its exact derivative with respect to that
  // abundance.  HD for every model; C (CI), C+ (CII), O (OI), OH, CO, H2O for the metal
  // models.
  double per_HD = 0.0;
  double per_CI = 0.0, per_CII = 0.0, per_OI = 0.0;
  double per_OH = 0.0, per_CO = 0.0, per_H2O = 0.0;
  double total() const { return line + cnt; }
};

// The start/end mean the energy equation uses, field by field (rates only).
inline Radiative mean_rates(const Radiative& a, const Radiative& b) {
  auto m = [](double u, double v) { return 0.5 * (u + v); };
  Radiative r;
  r.line = m(a.line, b.line);
  r.H2 = m(a.H2, b.H2);
  r.HD = m(a.HD, b.HD);
  r.Lya = m(a.Lya, b.Lya);
  r.CO = m(a.CO, b.CO);
  r.OH = m(a.OH, b.OH);
  r.H2O = m(a.H2O, b.H2O);
  r.CII = m(a.CII, b.CII);
  r.CI = m(a.CI, b.CI);
  r.OI = m(a.OI, b.OI);
  r.gr = m(a.gr, b.gr);
  r.gas = m(a.gas, b.gas);
  r.cnt = m(a.cnt, b.cnt);
  return r;
}

// The combination (1 - w) a + w b, field by field (rates, depth, grain temperature).
inline Radiative mix_rates(const Radiative& a, const Radiative& b, double w) {
  auto m = [w](double u, double v) { return (1.0 - w) * u + w * v; };
  Radiative r;
  r.line = m(a.line, b.line);
  r.H2 = m(a.H2, b.H2);
  r.HD = m(a.HD, b.HD);
  r.Lya = m(a.Lya, b.Lya);
  r.CO = m(a.CO, b.CO);
  r.OH = m(a.OH, b.OH);
  r.H2O = m(a.H2O, b.H2O);
  r.CII = m(a.CII, b.CII);
  r.CI = m(a.CI, b.CI);
  r.OI = m(a.OI, b.OI);
  r.gr = m(a.gr, b.gr);
  r.gas = m(a.gas, b.gas);
  r.cnt = m(a.cnt, b.cnt);
  r.k_gas = m(a.k_gas, b.k_gas);
  r.k_gr = m(a.k_gr, b.k_gr);
  r.T_gr = m(a.T_gr, b.T_gr);
  r.tau = m(a.tau, b.tau);
  r.esc = m(a.esc, b.esc);
  return r;
}

// Continuum channel at T: the depth from the column mass, and for grain models the grain
// temperature (seeded from p.T_gr_K) solved together with it. None of it depends on the
// abundances.
template <class Model>
void continuum_at(double nH, double T, const ChemParams& p,
                  const ChemShielding& in, Radiative& r) {
  if constexpr (Model::has_grain) {
    double T_gr = p.T_gr_K;
    cnt_cool_metal_consistent(nH, T, p.T_rad, in.cnt_column, p.Z_metal, T_gr,
                              r.k_gr, r.k_gas, r.gr, r.gas, r.cnt, r.tau,
                              r.esc);
    r.T_gr = T_gr;
  } else {
    const double rho = nH * ((1.0 + 4.0 * abundance_ref::yHe) * phys::m_p);
    r.tau = detail::eval_opacity(T, rho) * in.cnt_column;
    r.esc = r.tau > 1.0 ? 1.0 / (r.tau * r.tau) : 1.0;
    cnt_cool(nH, T, p.T_rad, r.esc, r.k_gas, r.gas, r.cnt);
    r.T_gr = p.T_gr_K;
  }
}

// Line cooling at (y, T) on the depth r.tau that continuum_at put into r. For the metal
// models: H2/HD from line_cool, the metal lines and Ly-alpha from line_cool_metal, whose
// escape probabilities are solved from (and written back to) `es`.
// slopes = true also fills the per-abundance line cooling (Radiative::per_*).
template <class Model>
void lines_at(const std::array<double, Model::N_sp>& y, double nH, double T,
              const ChemParams& p, const ChemShielding& in, EscOf<Model>& es,
              Radiative& r, bool slopes = false) {
  using Sp = typename Model::Sp;
  line_cool<Sp, Model::N_sp>(y, in.Nc_H2, in.Nc_HD, nH, T, p.T_rad, r.tau,
                             r.line, r.H2, r.HD, r.Lya);
  if (slopes) {
    const double rho = nH * ((1.0 + 4.0 * abundance_ref::yHe) * phys::m_p);
    r.per_HD = hd_cooling(nH, T, rho, 1.0, in.Nc_HD, r.tau) -
               hd_cooling(nH, p.T_rad, rho, 1.0, in.Nc_HD, r.tau);
  }
  if constexpr (Model::has_escape) {
    LineCoolRates m{};
    m.H2 = r.H2;
    m.HD = r.HD;
    m.Lya = r.Lya;
    line_cool_metal(nH, T, p.T_rad, y[Sp::H], y[Sp::H2], y[Sp::e], y[Sp::Hp],
                    y[Sp::He], y[Sp::C], y[Sp::Cp], y[Sp::O], y[Sp::OH],
                    y[Sp::CO], y[Sp::H2O], r.tau, in.Nc_CO, in.Nc_OH, in.Nc_H2O,
                    in.Nc_CII, in.Nc_CI, in.Nc_OI, es, m);
    r.Lya = m.Lya;
    r.CO = m.CO;
    r.OH = m.OH;
    r.H2O = m.H2O;
    r.CII = m.CII;
    r.CI = m.CI;
    r.OI = m.OI;
    r.line = m.total();
    r.per_CO = m.per_CO;
    r.per_OH = m.per_OH;
    r.per_H2O = m.per_H2O;
    r.per_CII = m.per_CII;
    r.per_CI = m.per_CI;
    r.per_OI = m.per_OI;
  } else {
    (void)es;
  }
}

template <class Model>
Radiative radiative_at(const std::array<double, Model::N_sp>& y, double nH,
                       double T, const ChemParams& p, const ChemShielding& in,
                       EscOf<Model>& es) {
  Radiative r;
  continuum_at<Model>(nH, T, p, in, r);
  lines_at<Model>(y, nH, T, p, in, es, r);
  return r;
}

// The species the line cooling reads (line_cool, and line_cool_metal for the metal
// models).
template <class Model>
constexpr bool enters_line_cooling(int j) {
  using Sp = typename Model::Sp;
  if (j == Sp::H || j == Sp::H2 || j == Sp::e || j == Sp::Hp || j == Sp::He ||
      j == Sp::HD)
    return true;
  if constexpr (Model::has_escape)
    return j == Sp::C || j == Sp::Cp || j == Sp::O || j == Sp::OH ||
           j == Sp::CO || j == Sp::H2O;
  return false;
}

// The line cooling's exact derivative with respect to y[j] where y[j] enters it only as
// an emitter's abundance (linearly); NaN for every other species.
template <class Model>
double line_slope(int j, const Radiative& r) {
  using Sp = typename Model::Sp;
  if (j == Sp::HD) return r.per_HD;
  if constexpr (Model::has_escape) {
    if (j == Sp::C) return r.per_CI;
    if (j == Sp::Cp) return r.per_CII;
    if (j == Sp::O) return r.per_OI;
    if (j == Sp::OH) return r.per_OH;
    if (j == Sp::CO) return r.per_CO;
    if (j == Sp::H2O) return r.per_H2O;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

// The species the energy row (energy_row) reads at a fixed temperature: the EOS (H, e,
// H+, He, He+, He++, H2), the chemistry heat (compute_chemistry_cooling, through y and
// through dy, whose Lyman-Werner part reads H2, HD and H-) and the line cooling. A
// change of any other abundance leaves every input of the row, hence its value, as it
// was.
template <class Model>
constexpr bool enters_energy_row(int j) {
  using Sp = typename Model::Sp;
  if (enters_line_cooling<Model>(j)) return true;
  if (j == Sp::Hm || j == Sp::H2p) return true;
  if constexpr (Model::cooling::has_he_ion)
    if (j == Sp::Hep) return true;
  if constexpr (Model::cooling::has_he_pp)
    if (j == Sp::Hepp) return true;
  if constexpr (Model::has_grain) {
    if constexpr (Model::cooling::grain_H2::has_surface)
      if (j == Sp::H_p || j == Sp::H_c) return true;
  }
  return false;
}

// The part of one residual evaluation that depends on the point (y, T) alone -- not on
// the step's width or start -- shared by the steps of one attempt: the full step and the
// first half step start at the same point (and take the same Jacobian there), and the
// second half step starts where the first one ended. The rate coefficients are the bulk
// of an evaluation (the metal kernels: ~70% of a one-zone run), so these are worth
// keeping. Keyed by the exact bits of (y, T); for grain models an entry may come from a
// step whose grain-temperature and escape-probability solves were seeded differently,
// which changes their results only within those solvers' tolerances.
// The Jacobian J is kept once a step has needed it at the point (has_J).
// (Members are left uninitialised: an entry is written whole before it is read, and
// zeroing ~100 kB per entry for the metal network would cost a visible fraction.)
template <class Model>
struct PointEval {
  static constexpr int N = Model::N_sp;
  std::array<double, N> y;
  double T;
  std::array<double, N> r;
  bool has_J;
  std::array<double, N * N> J;
  std::array<double, 2 * Model::N_react> k;
  std::array<double, 2 * Model::N_react> var;
  double Gamma;
  Radiative rad;
  EscOf<Model> es;
};

template <class Model>
class PointCache {
 public:
  static constexpr int kSize = 6;
  PointCache() {}  // user-provided: make_unique must not zero the entries
  PointEval<Model>* find(const std::array<double, Model::N_sp>& y, double T) {
    for (int i = 0; i < n_; ++i) {
      PointEval<Model>& e = entry_[i];
      if (e.T == T && e.y == y) return &e;
    }
    return nullptr;
  }
  PointEval<Model>& slot() {  // the entry to overwrite next (oldest first)
    PointEval<Model>& e = entry_[next_];
    next_ = (next_ + 1) % kSize;
    n_ = std::min(n_ + 1, kSize);
    return e;
  }
  void clear() { n_ = next_ = 0; }

 private:
  std::array<PointEval<Model>, kSize> entry_;
  int n_ = 0, next_ = 0;
};

// Fixed context of one step. p->T_gr_K is the grain temperature at the step's start
// (the seed of every grain-temperature solve in it); es0 the escape state there.
template <class Model>
struct Context {
  static constexpr int N = Model::N_sp;
  double nH = 0.0, h = 0.0, T_in = 0.0, e0 = 0.0;
  std::array<double, N> y0{};
  const ChemParams* p = nullptr;
  const ChemShielding* shield = nullptr;
  const ReactionTable<Model::N_sp, Model::N_react>* tbl = nullptr;
  Radiative rad_start;
  EscOf<Model> es0{};
  PointCache<Model>* cache =
      nullptr;  // optional (chem_integrate_adaptive's attempt)
  const GrainRelease<Model>* release =
      nullptr;      // grain evaporation of the call, or none
  bool lw = false;  // the Lyman-Werner channels are part of this step's rates
};

// The context of a step from (y0, T0); false when its start energy is not usable.
template <class Model>
bool make_context(double nH, const std::array<double, Model::N_sp>& y0,
                  double T0, double h, const ChemParams& p,
                  const ChemShielding& shield,
                  const ReactionTable<Model::N_sp, Model::N_react>& tbl,
                  const EscOf<Model>& es0, PointCache<Model>* cache,
                  Context<Model>& c,
                  const GrainRelease<Model>* release = nullptr,
                  bool lw = true) {
  c.nH = nH;
  c.h = h;
  c.T_in = T0;
  c.y0 = y0;
  c.e0 = eos_energy<Model>(y0, T0);
  c.p = &p;
  c.shield = &shield;
  c.tbl = &tbl;
  c.es0 = es0;
  c.cache = cache;
  c.release = release;
  c.lw = lw && shield.J_LW21 > 0.0;
  EscOf<Model> es = es0;
  c.rad_start = radiative_at<Model>(y0, nH, T0, p, shield, es);
  return std::isfinite(c.e0) && c.e0 > 0.0;
}

// One residual evaluation. J, the Jacobian of r, is filled only on demand
// (ensure_jacobian; has_J).
template <class Model>
struct Eval {
  static constexpr int N = Model::N_sp;
  std::array<double, N + 1> F{};
  std::array<double, N> r{};
  bool has_J = false;
  std::array<double, N * N> J{};
  LwRates lw;  // the Lyman-Werner rate coefficients at T (when Context::lw)
  std::array<double, 2 * Model::N_react> k{};
  std::array<double, 2 * Model::N_react> var{};
  double L_chem = 0.0, Gamma = 0.0;
  Radiative rad_end;
  EscOf<Model> es_end{};
};

// Rate coefficients at (y, T): LTE mix, UV-shielding inputs and mean molecular weight of
// this y, grain temperature T_gr (grain models).
template <class Model>
void base_rates(const Context<Model>& c,
                const std::array<double, Model::N_sp>& y, double T, double T_gr,
                std::array<double, 2 * Model::N_react>& k) {
  using Sp = typename Model::Sp;
  ChemParams pl = *c.p;
  pl.H = y[Sp::H];
  pl.H2 = 2.0 * y[Sp::H2];
  pl.He = y[Sp::He];
  if constexpr (Model::has_uv_shield) {
    pl.J_H2 = y[Sp::H2_p];
    pl.J_H2O = y[Sp::H2O_p];
    pl.J_tot = y[Sp::H_c] + y[Sp::D_c];
  }
  if constexpr (Model::has_grain) pl.T_gr_K = T_gr;
  (void)T_gr;
  const double mu =
      thermo::mean_molecular_weight(thermo::make_eos_composition<Model>(y));
  compute_base_rates<Model>(c.nH, T, mu, pl, *c.tbl, k);
}

// The sinks' part of the Jacobian (grain evaporation, then Lyman-Werner), added to the
// reactions' part in ev.J.
template <class Model>
void add_sink_jacobians(const Context<Model>& c, Eval<Model>& ev) {
  if (c.release) add_grain_sink_jacobian<Model>(*c.release, ev.J);
  if (c.lw) add_lw_jacobian<Model>(ev.lw, ev.J);
}

// Residual and rates at (y, T) into ev; its Jacobian only where the cache holds it
// (ensure_jacobian provides it when a step needs it).
template <class Model>
void evaluate(const Context<Model>& c, const std::array<double, Model::N_sp>& y,
              double T, Eval<Model>* ev) {
  constexpr int N = Model::N_sp;
  const PointEval<Model>* hit = c.cache ? c.cache->find(y, T) : nullptr;
  if (hit) {
    ev->r = hit->r;
    ev->has_J = hit->has_J;
    if (hit->has_J) ev->J = hit->J;
    ev->k = hit->k;
    ev->var = hit->var;
    ev->Gamma = hit->Gamma;
    ev->rad_end = hit->rad;
    ev->es_end = hit->es;
  } else {
    Radiative rad;
    continuum_at<Model>(c.nH, T, *c.p, *c.shield, rad);
    base_rates<Model>(c, y, T, rad.T_gr, ev->k);
    ev->var.fill(0.0);
    compute_rhs<Model>(ev->k, c.nH, y, *c.tbl, ev->r, ev->var);
    ev->has_J = false;
    ev->Gamma = cr_heating<Model>(ev->var);
    ev->es_end = c.es0;
    lines_at<Model>(y, c.nH, T, *c.p, *c.shield, ev->es_end, rad, true);
    ev->rad_end = rad;
    if (c.cache) {
      PointEval<Model>& e = c.cache->slot();
      e.y = y;
      e.T = T;
      e.r = ev->r;
      e.has_J = false;
      e.k = ev->k;
      e.var = ev->var;
      e.Gamma = ev->Gamma;
      e.rad = ev->rad_end;
      e.es = ev->es_end;
    }
  }
  // Grain evaporation: the grains' first-order sinks join the rates; the constant
  // rates are the part of dy that is not chemistry.
  if (c.release) add_grain_sink_rates<Model>(*c.release, y, ev->r);
  // Lyman-Werner: first-order sinks in the rates and the Jacobian. They change y at
  // fixed specific internal energy, so their part of dy stays out of the chemistry heat.
  std::array<double, N> rp{};
  if (c.lw) {
    ev->lw = lw_rates(c.nH, T, *c.shield);
    add_lw_sinks<Model>(ev->lw, y, rp, nullptr);
    for (int i = 0; i < N; ++i) ev->r[i] += rp[i];
  }
  if (ev->has_J) add_sink_jacobians<Model>(c, *ev);
  std::array<double, N> dy{};
  for (int i = 0; i < N; ++i) {
    dy[i] = y[i] - c.y0[i];
    if (c.release) dy[i] -= c.h * c.release->src[i];
    ev->F[i] = dy[i] - c.h * ev->r[i];
  }
  if (c.lw)
    for (int i = 0; i < N; ++i) dy[i] -= c.h * rp[i];
  compute_chemistry_cooling<Model>(y, dy, ev->k, c.nH, T, c.h, ev->L_chem);
  const double L_rad = 0.5 * (c.rad_start.total() + ev->rad_end.total());
  ev->F[N] = (eos_energy<Model>(y, T) - c.e0 +
              c.h * (ev->L_chem + L_rad - ev->Gamma)) /
             c.e0;
}

// The Jacobian at (y, T) of the evaluation ev made there: the reactions' part from ev's
// rate coefficients -- taken from the cache where the point's entry already has it, and
// kept there otherwise -- plus the sinks'. Most evaluations never need it (the
// temperature difference, rejected trial points, every point of a bracketed root).
template <class Model>
void ensure_jacobian(const Context<Model>& c,
                     const std::array<double, Model::N_sp>& y, double T,
                     Eval<Model>& ev) {
  if (ev.has_J) return;
  PointEval<Model>* e = c.cache ? c.cache->find(y, T) : nullptr;
  if (e && !(e->k == ev.k)) e = nullptr;  // an entry from another evaluation
  if (e && e->has_J) {
    ev.J = e->J;
  } else {
    compute_jacobian<Model>(ev.k, c.nH, y, *c.tbl, ev.J);
    if (e) {
      e->J = ev.J;
      e->has_J = true;
    }
  }
  add_sink_jacobians<Model>(c, ev);
  ev.has_J = true;
}

// What the energy row's columns at one temperature share: the Lyman-Werner rate
// coefficients and the chemistry heat's factors that depend on T alone.
struct RowTerms {
  LwRates lw;
  ChemCoolingTerms cool;
};

// Energy row at abundances y and the temperature of `base`, reusing base's rate
// coefficients, CR heating and continuum channel (Jacobian only: the residual stays
// exact).  The line cooling is base's plus line_shift, or re-evaluated at y when
// with_line.
template <class Model>
double energy_row(const Context<Model>& c,
                  const std::array<double, Model::N_sp>& y, double T,
                  const Eval<Model>& base, const RowTerms& tt, bool with_line,
                  double line_shift = 0.0) {
  constexpr int N = Model::N_sp;
  std::array<double, N> dy{};
  for (int i = 0; i < N; ++i) {
    dy[i] = y[i] - c.y0[i];
    if (c.release) dy[i] -= c.h * c.release->src[i];
  }
  if (c.lw) {  // the Lyman-Werner part of dy, as in evaluate
    std::array<double, N> rp{};
    add_lw_sinks<Model>(tt.lw, y, rp, nullptr);
    for (int i = 0; i < N; ++i) dy[i] -= c.h * rp[i];
  }
  double L_chem = 0.0;
  compute_chemistry_cooling<Model>(y, dy, base.k, c.nH, c.h, tt.cool, L_chem);
  double line = base.rad_end.line + line_shift;
  if (with_line) {
    Radiative r = base.rad_end;
    EscOf<Model> es = base.es_end;  // warm: the base evaluation's solution
    lines_at<Model>(y, c.nH, T, *c.p, *c.shield, es, r);
    line = r.line;
  }
  const double L_rad = 0.5 * (c.rad_start.total() + line + base.rad_end.cnt);
  return (eos_energy<Model>(y, T) - c.e0 +
          c.h * (L_chem + L_rad - base.Gamma)) /
         c.e0;
}

// One backward-Euler step of width h from (y0, T0).
template <class Model>
struct StepResult {
  bool ok = false;
  SolveFailure reason = SolveFailure::max_iterations;
  int newton = 0;
  std::array<double, Model::N_sp> y{};
  double T = 0.0;
  std::array<double, 2 * Model::N_react> var{};
  double L_chem = 0.0, Gamma = 0.0;
  Radiative rad_mean;  // start/end mean the energy equation used
  Radiative rad_end;   // the end state's channel, depth and grain temperature
  EscOf<Model> es_end{};
  bool fallback = false;  // solved by fixed_T_step (the coupled Newton failed)
  // The coupled Newton ended with a factorised matrix in the caller's LU workspace
  // (StepGuess): a following step of the same width may start from it.
  bool lu_valid = false;
};

// What a step may start from instead of its start state, all of it optional: an
// initial Newton iterate (y, T) -- e.g. the full step's end for the second half step,
// which ends at the same time -- the Newton matrix already factorised in the LU
// workspace by a previous step of the same width (lu), and the temperature at which a
// bracketed temperature root starts its search (T_seed; 0 = the start temperature).
// They change where the iteration starts, not the equations it solves.
template <class Model>
struct StepGuess {
  const std::array<double, Model::N_sp>* y = nullptr;
  double T = 0.0;
  bool lu = false;
  double T_seed = 0.0;
};

template <class Model>
StepResult<Model> coupled_step(
    double nH, const std::array<double, Model::N_sp>& y0, double T0, double h,
    const ChemParams& p, const ChemShielding& shield, const EscOf<Model>& es0,
    const ReactionTable<Model::N_sp, Model::N_react>& tbl, const Tolerance& tol,
    PointCache<Model>* cache = nullptr, bool end_rates = true,
    const GrainRelease<Model>* release = nullptr,
    const StepGuess<Model>* guess = nullptr,
    DenseLU<Model::N_sp + 1>* lu_work = nullptr) {
  constexpr int N = Model::N_sp;
  constexpr int M = N + 1;
  StepResult<Model> out;

  Context<Model> c;
  if (!make_context<Model>(nH, y0, T0, h, p, shield, tbl, es0, cache, c,
                           release)) {
    out.reason = SolveFailure::nonfinite;
    return out;
  }

  auto wnorm_step = [&](const std::array<double, N>& xa, double Ta,
                        const std::array<double, N>& xb, double Tb) {
    return weighted_rms<N>(xa, Ta, xb, Tb, tol);
  };
  auto norm2 = [](const std::array<double, M>& v) {
    double s = 0.0;
    for (double a : v) s += a * a;
    return std::sqrt(s);
  };

  std::array<double, N> x = y0;
  double T = T0;
  if (guess && guess->y) {
    x = *guess->y;
    T = guess->T;
  }
  Eval<Model> ev;
  evaluate<Model>(c, x, T, &ev);
  bool converged = false;
  double last_w = std::numeric_limits<double>::infinity();

  std::array<double, M * M> A{};
  // The factorisation lives here and is copied from / to lu_work: Eigen's vectorised
  // kernels for a matrix whose size is not a multiple of the vector width depend on
  // the storage's alignment, so solving in place in the workspace would change the
  // rounding.
  DenseLU<M> lu;
  bool have_jacobian = false, jacobian_fresh = false, have_direction = false;
  // A matrix factorised by an earlier step is used as a chord: its first correction is
  // taken only through the damped trial below, which tests that it contracts.
  bool inherited = guess && guess->lu && lu_work;
  if (inherited) {
    lu = *lu_work;
    have_jacobian = true;
  }
  std::array<double, M> d{};

  for (int it = 0; it < kMaxNewton; ++it) {
    bool finite = std::isfinite(ev.F[N]);
    for (int i = 0; i < N && finite; ++i) finite = std::isfinite(ev.F[i]);
    if (!finite) {
      out.reason = SolveFailure::nonfinite;
      break;
    }
    out.newton = it + 1;
    if (!have_jacobian) {
      ensure_jacobian<Model>(c, x, T, ev);
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
          A[i * M + j] = (i == j ? 1.0 : 0.0) - h * ev.J[i * N + j];
      const double dT = kFdT * T;
      Eval<Model> e2;
      evaluate<Model>(c, x, T + dT, &e2);
      for (int i = 0; i < M; ++i) A[i * M + N] = (e2.F[i] - ev.F[i]) / dT * T0;
      RowTerms tt;
      if (c.lw) tt.lw = ev.lw;
      tt.cool = chemistry_cooling_terms(c.nH, T);
      const double f0 = energy_row<Model>(c, x, T, ev, tt, false);
      // Finite differences, except that an emitter-only species shifts the line
      // cooling by its exact slope times the step instead of re-evaluating it, and
      // that a species the row does not read leaves it at f0.
      for (int j = 0; j < N; ++j) {
        const double dy = kFdYRel * std::max(std::abs(x[j]), kFdYFloor);
        if (!enters_energy_row<Model>(j)) {
          A[N * M + j] = (f0 - f0) / dy;
          continue;
        }
        std::array<double, N> x2 = x;
        x2[j] += dy;
        const double slope = line_slope<Model>(j, ev.rad_end);
        const double row =
            std::isnan(slope)
                ? energy_row<Model>(c, x2, T, ev, tt,
                                    enters_line_cooling<Model>(j))
                : energy_row<Model>(c, x2, T, ev, tt, false,
                                    slope * (x2[j] - x[j]));
        A[N * M + j] = (row - f0) / dy;
      }
      if (!lu.factor(A)) {
        out.reason = SolveFailure::singular_jacobian;
        break;
      }
      have_jacobian = jacobian_fresh = true;
      have_direction = false;
      inherited = false;
    }
    if (!have_direction) {
      for (int i = 0; i < M; ++i) d[i] = -ev.F[i];
      lu.solve(d);
    }
    if (!inherited) {  // the full correction is already inside the tolerance: stop
      std::array<double, N> xf{};
      for (int i = 0; i < N; ++i) {
        const double v = x[i] + d[i];
        xf[i] = v < 0.0 ? 0.1 * x[i] : v;
      }
      const double Tf = T + d[N] * T0;
      if (Tf > 1.0 && wnorm_step(x, T, xf, Tf) <= kNewtonKappa) {
        if (!end_rates) {  // the caller uses only the end state (y, T)
          x = xf;
          T = Tf;
          converged = true;
          break;
        }
        Eval<Model> ef;
        evaluate<Model>(c, xf, Tf, &ef);
        bool fin = std::isfinite(ef.F[N]);
        for (int i = 0; i < N && fin; ++i) fin = std::isfinite(ef.F[i]);
        if (fin) {
          x = xf;
          T = Tf;
          ev = ef;
          converged = true;
          break;
        }
      }
    }
    const double dnorm = norm2(d);
    double alpha = 1.0;
    if (std::abs(d[N]) * T0 > kMaxRelativeTStep * T)
      alpha = kMaxRelativeTStep * T / (std::abs(d[N]) * T0);
    bool accepted = false;
    std::array<double, N> xt{};
    double Tt = T;
    Eval<Model> et;
    std::array<double, M> dbar{};
    double dbar_norm = 0.0;
    const int backtracks = jacobian_fresh ? kMaxBacktracks : 1;
    for (int bt = 0; bt < backtracks; ++bt) {
      for (int i = 0; i < N; ++i) {
        const double v = x[i] + alpha * d[i];
        xt[i] = v < 0.0 ? 0.1 * x[i] : v;  // the library's positivity clamp
      }
      Tt = T + alpha * d[N] * T0;
      if (Tt > 1.0) {
        evaluate<Model>(c, xt, Tt, &et);
        bool fin = std::isfinite(et.F[N]);
        for (int i = 0; i < N && fin; ++i) fin = std::isfinite(et.F[i]);
        if (fin) {
          for (int i = 0; i < M; ++i) dbar[i] = -et.F[i];
          lu.solve(dbar);
          dbar_norm = norm2(dbar);
          if (dbar_norm <= (1.0 - 0.25 * alpha) * dnorm) {
            accepted = true;
            break;
          }
        }
      }
      alpha *= 0.5;
    }
    if (!accepted) {
      if (!jacobian_fresh) {
        have_jacobian = false;
        --it;
        continue;
      }
      if (last_w <=
          1.0)  // damping stalled with the last update inside the tolerance
        converged = true;
      else
        out.reason = SolveFailure::line_search;
      break;
    }
    inherited = false;
    last_w = wnorm_step(x, T, xt, Tt);
    x = xt;
    T = Tt;
    ev = et;
    if (last_w <= kNewtonKappa) {
      converged = true;
      break;
    }
    if (alpha == 1.0 && dbar_norm <= 0.25 * dnorm) {  // chord step next
      d = dbar;
      have_direction = true;
      jacobian_fresh = false;
    } else {
      have_jacobian = false;
      have_direction = false;
    }
  }
  if (converged) {
    for (double v : x)
      if (!std::isfinite(v) || v < 0.0 || v > 1.0e50) converged = false;
    if (!std::isfinite(T) || !(T > 0.0)) converged = false;
    if (!converged) out.reason = SolveFailure::nonfinite;
  }
  if (!converged) return out;

  out.ok = true;
  out.reason = SolveFailure::none;
  out.y = x;
  out.T = T;
  out.var = ev.var;
  out.L_chem = ev.L_chem;
  out.Gamma = ev.Gamma;
  out.rad_mean = mean_rates(c.rad_start, ev.rad_end);
  out.rad_end = ev.rad_end;
  out.es_end = ev.es_end;
  if (lu_work && have_jacobian) {
    *lu_work = lu;
    out.lu_valid = true;
  }
  return out;
}

// Species at a fixed temperature T: backward Euler Newton on y alone (undamped, the
// library's positivity clamp), converged when the weighted RMS of the correction is at
// most kNewtonKappa. `y` is the starting guess and receives the solution.
template <class Model>
bool species_at_T(const Context<Model>& c, double T,
                  std::array<double, Model::N_sp>& y, const Tolerance& tol,
                  int& newton) {
  constexpr int N = Model::N_sp;
  std::array<double, 2 * Model::N_react> k{}, var{};
  std::array<double, N> r{};
  std::array<double, N * N> J{}, A{};
  DenseLU<N> lu;
  Radiative rc;  // the grain temperature at T (it does not depend on y)
  continuum_at<Model>(c.nH, T, *c.p, *c.shield, rc);
  const LwRates lw = c.lw ? lw_rates(c.nH, T, *c.shield) : LwRates{};
  for (int it = 0; it < kMaxNewton; ++it) {
    ++newton;
    base_rates<Model>(c, y, T, rc.T_gr, k);
    compute_rates<Model>(k, c.nH, y, *c.tbl, r, J, var);
    if (c.release) add_grain_sinks<Model>(*c.release, y, r, J);
    if (c.lw) add_lw_sinks<Model>(lw, y, r, &J);
    std::array<double, N> d{};
    for (int i = 0; i < N; ++i) {
      d[i] = -(y[i] - c.y0[i] - c.h * r[i]);
      if (c.release) d[i] += c.h * c.release->src[i];
      if (!std::isfinite(d[i])) return false;
      for (int j = 0; j < N; ++j)
        A[i * N + j] = (i == j ? 1.0 : 0.0) - c.h * J[i * N + j];
    }
    if (!lu.factor(A)) return false;
    lu.solve(d);
    std::array<double, N> yt{};
    for (int i = 0; i < N; ++i) {
      const double v = y[i] + d[i];
      yt[i] = v < 0.0 ? 0.1 * y[i] : v;
    }
    const double w = weighted_rms<N>(y, T, yt, T, tol);
    y = yt;
    if (!std::isfinite(w)) return false;
    if (w <= kNewtonKappa) {
      for (double v : y)
        if (!std::isfinite(v) || v < 0.0 || v > 1.0e50) return false;
      return true;
    }
  }
  return false;
}

// Fallback for a step the coupled Newton cannot solve: the species at a fixed trial
// temperature (species_at_T) and the energy equation G(T) = F_N as a bracketed root in T
// (Illinois false position). Across a temperature at which the model's right-hand side
// jumps, G may change sign without a root; the bracket then closes on the jump with a
// residual of opposite sign at its two ends. The step then ends on the jump with the
// convex combination of the two ends, y = (1 - a) y_a + a y_b with a = G_a / (G_a - G_b),
// and the same combination of their chemistry heat, radiative and CR rates: the solution
// that stays on a discontinuity whose two sides both drive the state towards it, the
// fraction a of the step taken under the right-hand side of end b. The specific internal
// energy is linear in y at fixed T, so the combination satisfies the energy equation up
// to the rounding of the two ends' temperatures, and it keeps every element and the
// charge, which both ends keep. The model's equations are not changed.
// The bracketed temperature root shared by fixed_T_step and saha_step: the species at a
// trial temperature come from `species_at(c, T, y)` (y holds a starting guess on entry and
// the composition on exit; false = no composition at this T). lw: the Lyman-Werner
// channels are part of the step's rates (false for the Saha composition).
template <class Model, class SpeciesAt>
StepResult<Model> bracketed_T_step(
    double nH, const std::array<double, Model::N_sp>& y0, double T0, double h,
    const ChemParams& p, const ChemShielding& shield, const EscOf<Model>& es0,
    const ReactionTable<Model::N_sp, Model::N_react>& tbl, const Tolerance& tol,
    SpeciesAt species_at, PointCache<Model>* cache = nullptr,
    const GrainRelease<Model>* release = nullptr, bool lw = true,
    double T_seed = 0.0) {
  constexpr int N = Model::N_sp;
  StepResult<Model> out;
  Context<Model> c;
  if (!make_context<Model>(nH, y0, T0, h, p, shield, tbl, es0, cache, c,
                           release, lw)) {
    out.reason = SolveFailure::nonfinite;
    return out;
  }
  std::array<double, N> yw =
      y0;  // warm start carried across trial temperatures
  auto G = [&](double T, std::array<double, N>& y_at, Eval<Model>& ev) {
    y_at = yw;
    if (!species_at(c, T, y_at))
      return std::numeric_limits<double>::quiet_NaN();
    yw = y_at;
    evaluate<Model>(c, y_at, T, &ev);
    return ev.F[N];
  };
  std::array<double, N> ya{}, yb{};
  Eval<Model> ea, eb;
  // A seed (an estimate of the root) starts the search there with a first step of
  // rtol_T, growing fourfold up to the unseeded step; the root is the same.
  bool seeded = T_seed > 0.0 && std::isfinite(T_seed) && T_seed != T0;
  double Ta = seeded ? T_seed : T0, Ga = G(Ta, ya, ea);
  if (seeded && !std::isfinite(Ga)) {
    seeded = false;
    Ta = T0;
    Ga = G(Ta, ya, ea);
  }
  if (!std::isfinite(Ga)) {
    out.reason = SolveFailure::temperature_bracket;
    return out;
  }
  // Bracket: G increases with T (energy and cooling both do), so step away from T0 in the
  // direction that changes the sign of G.
  double Tb = Ta, Gb = Ga;
  const double factor = Ga > 0.0 ? 0.8 : 1.25;
  double rel = tol.rtol_T;  // the seeded step
  bool bracketed = (Ga == 0.0);
  for (int k = 0; k < 60 && !bracketed; ++k) {
    if (seeded && rel < 0.25) {
      Tb = Ga > 0.0 ? Ta / (1.0 + rel) : Ta * (1.0 + rel);
      rel *= 4.0;
    } else {
      Tb = Ta * factor;
    }
    Gb = G(Tb, yb, eb);
    if (!std::isfinite(Gb)) {
      out.reason = SolveFailure::temperature_bracket;
      return out;
    }
    if ((Gb > 0.0) != (Ga > 0.0) || Gb == 0.0) {
      bracketed = true;
    } else {
      Ta = Tb;
      Ga = Gb;
      ya = yb;
      ea = eb;
    }
  }
  if (!bracketed) {
    out.reason = SolveFailure::temperature_bracket;
    return out;
  }
  if (Ga == 0.0) {  // T0 is the root: both ends are that evaluation
    Tb = Ta;
    Gb = Ga;
    yb = ya;
    eb = ea;
  }
  // Illinois false position on [Ta, Tb] (either order), to the rounding level of T: the
  // residual left at a smooth root is carried into the chemistry heat, so it must be
  // negligible, and where the composition depends steeply on T (the Saha regime) a
  // bracket of even 1e-5 T leaves a residual that, carried into the heat call after call,
  // drifts the energy.
  const double t_tol = kRootRelTol;
  int side = 0;
  bool closed = (Ga == 0.0 || Gb == 0.0);
  for (int it = 0; it < 400 && !closed; ++it) {
    if (std::abs(Tb - Ta) <= t_tol * std::max(Ta, Tb)) {
      closed = true;
      break;
    }
    double Tc = (Ta * Gb - Tb * Ga) / (Gb - Ga);
    if (!(Tc > std::min(Ta, Tb) && Tc < std::max(Ta, Tb))) Tc = 0.5 * (Ta + Tb);
    std::array<double, N> yc{};
    Eval<Model> ec;
    const double Gc = G(Tc, yc, ec);
    if (!std::isfinite(Gc)) {
      out.reason = SolveFailure::temperature_iterations;
      return out;
    }
    if ((Gc > 0.0) == (Gb > 0.0)) {
      Tb = Tc;
      Gb = Gc;
      yb = yc;
      eb = ec;
      if (side == -1) Ga *= 0.5;
      side = -1;
    } else {
      Ta = Tc;
      Ga = Gc;
      ya = yc;
      ea = ec;
      if (side == +1) Gb *= 0.5;
      side = +1;
    }
  }
  // An unclosed bracket is neither a root nor a jump: refuse the step (it is shortened).
  if (!closed && !(Ga == 0.0 || Gb == 0.0) &&
      !(std::abs(Tb - Ta) <= t_tol * std::max(Ta, Tb))) {
    out.reason = SolveFailure::temperature_iterations;
    return out;
  }
  // (Ga/Gb carry the Illinois down-weighting; the ends' own residuals are ea/eb.F[N].)
  // A closed bracket whose two ends both keep a residual above kJumpResidual, of opposite
  // sign, lies on a jump of the right-hand side: at a smooth root the ends' residuals are
  // at the rounding level of T. The step then ends on the jump (see above).
  const double Fa = ea.F[N], Fb = eb.F[N];
  if (std::abs(Fa) > kJumpResidual && std::abs(Fb) > kJumpResidual &&
      (Fa > 0.0) != (Fb > 0.0)) {
    const double w = Fa / (Fa - Fb);  // weight of end b, in (0, 1)
    auto mix = [w](double u, double v) { return (1.0 - w) * u + w * v; };
    std::array<double, N> ym{};
    for (int i = 0; i < N; ++i) ym[i] = mix(ya[i], yb[i]);
    const Radiative rm = mix_rates(ea.rad_end, eb.rad_end, w);
    const double Tm = mix(Ta, Tb);
    const double L_chem = mix(ea.L_chem, eb.L_chem);
    const double Gamma = mix(ea.Gamma, eb.Gamma);
    // What the combination leaves of the energy equation: the nonlinearity in T across
    // the rounding-level width of the bracket. Carried into the chemistry heat.
    const double Gm =
        (eos_energy<Model>(ym, Tm) - c.e0 +
         h * (L_chem + 0.5 * (c.rad_start.total() + rm.total()) - Gamma)) /
        c.e0;
    const Eval<Model>& major = w < 0.5 ? ea : eb;
    out.ok = true;
    out.reason = SolveFailure::none;
    out.y = ym;
    out.T = Tm;
    out.var = major.var;
    out.L_chem = L_chem - Gm * c.e0 / h;
    out.Gamma = Gamma;
    out.rad_mean = mean_rates(c.rad_start, rm);
    out.rad_end = rm;
    out.es_end = major.es_end;
    return out;
  }
  // Keep the end with the smaller residual; carry what is left into the chemistry heat.
  const bool take_a = std::abs(Fa) <= std::abs(Fb);
  const double T = take_a ? Ta : Tb;
  const Eval<Model>& ev = take_a ? ea : eb;
  const double Gr = ev.F[N];
  out.ok = true;
  out.reason = SolveFailure::none;
  out.y = take_a ? ya : yb;
  out.T = T;
  out.var = ev.var;
  out.L_chem = ev.L_chem - Gr * c.e0 / h;
  out.Gamma = ev.Gamma;
  out.rad_mean = mean_rates(c.rad_start, ev.rad_end);
  out.rad_end = ev.rad_end;
  out.es_end = ev.es_end;
  return out;
}

template <class Model>
StepResult<Model> fixed_T_step(
    double nH, const std::array<double, Model::N_sp>& y0, double T0, double h,
    const ChemParams& p, const ChemShielding& shield, const EscOf<Model>& es0,
    const ReactionTable<Model::N_sp, Model::N_react>& tbl, const Tolerance& tol,
    PointCache<Model>* cache = nullptr,
    const GrainRelease<Model>* release = nullptr, double T_seed = 0.0) {
  int newton = 0;
  auto at = [&](const Context<Model>& c, double T,
                std::array<double, Model::N_sp>& y) {
    return species_at_T<Model>(c, T, y, tol, newton);
  };
  StepResult<Model> out =
      bracketed_T_step<Model>(nH, y0, T0, h, p, shield, es0, tbl, tol, at,
                              cache, release, true, T_seed);
  out.fallback = true;
  out.newton += newton;
  return out;
}

// Saha composition at T with the element and charge totals of y_in: the model's own
// SahaPolicy. false when it does not return a finite, non-negative composition.
template <class Model>
bool saha_at(double nH, double T, const ChemParams& p,
             const ReactionTable<Model::N_sp, Model::N_react>& tbl,
             const std::array<double, Model::N_sp>& y_in,
             std::array<double, Model::N_sp>& y) {
  y = y_in;
  Model::SahaPolicy::template solve<Model::N_sp, Model::N_react>(nH, T, p, y,
                                                                 tbl);
  for (double v : y)
    if (!std::isfinite(v) || v < 0.0) return false;
  return true;
}

// A step in the Saha regime: the composition is the Saha state y_S(T) (element totals of
// the step's start state), the unknown is T alone, and the energy equation is the same
// as the network's (bracketed_T_step). The composition is the model's SahaPolicy; the
// temperature follows the EOS energy equation.
template <class Model>
StepResult<Model> saha_step(
    double nH, const std::array<double, Model::N_sp>& y0, double T0, double h,
    const ChemParams& p, const ChemShielding& shield, const EscOf<Model>& es0,
    const ReactionTable<Model::N_sp, Model::N_react>& tbl, const Tolerance& tol,
    PointCache<Model>* cache = nullptr,
    const GrainRelease<Model>* release = nullptr, double T_seed = 0.0) {
  // The Saha totals are those of the start state after the call's grain release over h.
  std::array<double, Model::N_sp> y_start = y0;
  if (release) apply_grain_release_to_start<Model>(*release, h, y_start);
  auto at = [&](const Context<Model>& c, double T,
                std::array<double, Model::N_sp>& y) {
    return saha_at<Model>(c.nH, T, p, tbl, y_start, y);
  };
  return bracketed_T_step<Model>(nH, y0, T0, h, p, shield, es0, tbl, tol, at,
                                 cache, release, false, T_seed);
}

// Whether a step from (y0, T0) is taken in the Saha regime: above saha_threshold_nH when
// use_saha, and -- with non-thermal ionisation -- only where saha_nonthermal_shift at the
// Saha state of T0 is within saha_nonthermal_tol.
template <class Model>
bool saha_regime(double nH, const std::array<double, Model::N_sp>& y0,
                 double T0, const ChemParams& p,
                 const ReactionTable<Model::N_sp, Model::N_react>& tbl,
                 const std::array<double, 2 * Model::N_react>& var,
                 double* shift) {
  *shift = -1.0;
  if (!(p.use_saha && nH > p.saha_threshold_nH)) return false;
  const bool nonthermal = (p.zeta > 0.0) || (p.zeta_X > 0.0);
  if (!nonthermal) return true;
  if (!(p.saha_nonthermal_tol > 0.0)) return false;
  std::array<double, Model::N_sp> yS{};
  if (!saha_at<Model>(nH, T0, p, tbl, y0, yS)) return false;
  *shift = saha_nonthermal_shift<Model>(nH, T0, yS, tbl, p, var);
  return *shift <= p.saha_nonthermal_tol;
}

// One backward-Euler step: the coupled Newton, and fixed_T_step where it fails.
// end_rates = false: the caller needs only the end state (y, T), not the rates there.
// guess / lu_work: see StepGuess. A coupled Newton that fails from a guess is repeated
// from the start state before the fallback.
template <class Model>
StepResult<Model> be_step(double nH, const std::array<double, Model::N_sp>& y0,
                          double T0, double h, const ChemParams& p,
                          const ChemShielding& shield, const EscOf<Model>& es0,
                          const ReactionTable<Model::N_sp, Model::N_react>& tbl,
                          const Tolerance& tol,
                          PointCache<Model>* cache = nullptr,
                          bool end_rates = true,
                          const GrainRelease<Model>* release = nullptr,
                          const StepGuess<Model>* guess = nullptr,
                          DenseLU<Model::N_sp + 1>* lu_work = nullptr) {
  StepResult<Model> s =
      coupled_step<Model>(nH, y0, T0, h, p, shield, es0, tbl, tol, cache,
                          end_rates, release, guess, lu_work);
  if (s.ok) return s;
  int newton = s.newton;
  if (guess && (guess->y || guess->lu)) {
    s = coupled_step<Model>(nH, y0, T0, h, p, shield, es0, tbl, tol, cache,
                            end_rates, release, nullptr, lu_work);
    s.newton += newton;
    if (s.ok) return s;
    newton = s.newton;
  }
  StepResult<Model> f =
      fixed_T_step<Model>(nH, y0, T0, h, p, shield, es0, tbl, tol, cache,
                          release, guess ? guess->T_seed : 0.0);
  f.newton += newton;
  return f;
}

// Rates of one solved step, in the ChemFullRates fields the interval average carries.
template <class Model>
void add_step_rates(const StepResult<Model>& s, double w, ChemFullRates& acc) {
  acc.Lambda_chem += w * s.L_chem;
  acc.Gamma_CR += w * s.Gamma;
  acc.Lambda_line += w * s.rad_mean.line;
  acc.Lambda_H2 += w * s.rad_mean.H2;
  acc.Lambda_HD += w * s.rad_mean.HD;
  acc.Lambda_Lya += w * s.rad_mean.Lya;
  acc.Lambda_CO += w * s.rad_mean.CO;
  acc.Lambda_OH += w * s.rad_mean.OH;
  acc.Lambda_H2O += w * s.rad_mean.H2O;
  acc.Lambda_CII += w * s.rad_mean.CII;
  acc.Lambda_CI += w * s.rad_mean.CI;
  acc.Lambda_OI += w * s.rad_mean.OI;
  acc.Lambda_gr += w * s.rad_mean.gr;
  acc.Lambda_gas += w * s.rad_mean.gas;
  acc.Lambda_cnt += w * s.rad_mean.cnt;
}

}  // namespace integ
}  // namespace arche
