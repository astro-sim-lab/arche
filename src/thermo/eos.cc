// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE

#include <array>
#include <cmath>
#include <limits>

#include "core/exact_memo.h"
#include "model_traits.h"
#include "thermo/eos.h"
#include "thermo/h2_thermodynamics.h"

// Keep this implementation in a translation unit separate from arche_api.cc
// for numerical reproducibility. Adding it to that already-large unit changes
// the compiler's optimization budget and can alter ULP-level code generation
// for inline-template COMDAT definitions emitted there (including metal line
// cooling), which then changes existing applications linked to libarche.
// LTO/IPO removes this isolation and can reintroduce the same coupling.

namespace arche::thermo {

double c_H2(double T_K) {
  thread_local ExactMemo<1, 1, 1> memo;
  std::array<double, 1> out{};
  if (memo.find({T_K}, out)) return out[0];
  out[0] = detail::c_H2_eval(T_K);
  memo.insert({T_K}, out);
  return out[0];
}

namespace {

// The specific internal energy at T = exp(L).  Every evaluation below goes through this
// one out-of-line function, so the margin bracket, the clamps and both bisections use
// the same machine code for e(T) whatever the compiler inlines elsewhere.
[[gnu::noinline]] double energy_at_log_T(const EosComposition& c, double L) {
  return specific_internal_energy(c, std::exp(L));
}

// Bisection in ln T over [ln 1, ln 1e12] to full double precision.  A midpoint at or
// below La is taken as below e_cgs and one at or above Lb as above it without
// evaluating the energy; every other midpoint is compared exactly as in the plain
// bisection, which (La, Lb) = (-inf, +inf) is.
double bisect_log_T(const EosComposition& c, double e_cgs, double La,
                    double Lb) {
  double log_lo = std::log(1.0);
  double log_hi = std::log(1.0e12);
  for (int it = 0; it < 80; ++it) {
    const double log_mid = 0.5 * (log_lo + log_hi);
    // The bracket is at full double precision once its midpoint rounds to an
    // endpoint. Further c_H2() evaluations cannot refine either bound.
    if (log_mid == log_lo || log_mid == log_hi) break;
    if (log_mid <= La) {
      log_lo = log_mid;
      continue;
    }
    if (log_mid >= Lb) {
      log_hi = log_mid;
      continue;
    }
    const double e_mid = energy_at_log_T(c, log_mid);
    if (e_mid < e_cgs) {
      log_lo = log_mid;
    } else {
      log_hi = log_mid;
    }
  }
  return std::exp(0.5 * (log_lo + log_hi));
}

// A bracket [La, Lb] of ln T inside (ln 1, ln 1e12) with e(exp La) below e_cgs and
// e(exp Lb) above it, both by a relative margin of 1e-12.  e(T) increases with T, so
// a midpoint at or below La lies below e_cgs and one at or above Lb above it wherever
// the rounding error of the computed e stays below that margin; the bisection then
// returns what the plain bisection returns.  The rounding error of c_H2 (and with it
// of e) reaches ~1e-11 relative below 100 K, where the finite-difference rotational
// heat capacity loses digits, so there the result can differ from the plain
// bisection's at that level.  The centre comes from a secant iteration on
// ln e(exp L) - ln e_cgs.  False where no such bracket is found.
bool margin_bracket(const EosComposition& c, double e_cgs, double& La,
                    double& Lb) {
  const double L_min = std::log(1.0), L_max = std::log(1.0e12);
  auto energy = [&](double L) { return energy_at_log_T(c, L); };
  const double target = std::log(e_cgs);
  double L0 = 3.0;
  double f0 = std::log(energy(L0)) - target;
  double L1 = L0 - f0;
  for (int it = 0; it < 12; ++it) {
    if (!(L1 > L_min + 1.0e-3 && L1 < L_max - 1.0e-3)) return false;
    const double f1 = std::log(energy(L1)) - target;
    if (!std::isfinite(f1)) return false;
    if (std::abs(f1) < 1.0e-13 || f1 == f0) break;
    const double L2 = L1 - f1 * (L1 - L0) / (f1 - f0);
    L0 = L1;
    f0 = f1;
    L1 = L2;
  }
  double h = 1.0e-11;
  for (int w = 0; w < 6; ++w, h *= 10.0) {
    La = L1 - h;
    Lb = L1 + h;
    if (!(La > L_min && Lb < L_max)) return false;
    if (energy(La) < e_cgs * (1.0 - 1.0e-12) &&
        energy(Lb) > e_cgs * (1.0 + 1.0e-12))
      return true;
  }
  return false;
}

// The root of ln e(exp L) = ln e_cgs by secant iteration from L = 3 (T = 20 K).  e(T)
// is smooth and increasing, so the iteration converges in a few steps; it stops once
// a step changes L by at most 1e-14, i.e. T to that relative precision, below the
// rounding level of e itself.  False where it leaves (ln 1, ln 1e12) or does not
// converge; the bisection then decides.
bool secant_log_T(const EosComposition& c, double e_cgs, double& L) {
  const double L_min = std::log(1.0), L_max = std::log(1.0e12);
  const double target = std::log(e_cgs);
  double L0 = 3.0;
  double f0 = std::log(energy_at_log_T(c, L0)) - target;
  double L1 = L0 - f0;
  for (int it = 0; it < 30; ++it) {
    if (!(L1 > L_min + 1.0e-3 && L1 < L_max - 1.0e-3)) return false;
    const double f1 = std::log(energy_at_log_T(c, L1)) - target;
    if (!std::isfinite(f1)) return false;
    if (f1 == 0.0) {
      L = L1;
      return true;
    }
    if (f1 == f0) return false;
    const double L2 = L1 - f1 * (L1 - L0) / (f1 - f0);
    if (std::abs(L2 - L1) <= 1.0e-14) {
      L = L2;
      return true;
    }
    L0 = L1;
    f0 = f1;
    L1 = L2;
  }
  return false;
}

}  // namespace

double temperature_from_specific_internal_energy(const EosComposition& c,
                                                 double e_cgs) noexcept {
  if (std::isnan(e_cgs)) return e_cgs;

  // Inside (1 K, 1e12 K) the secant root; e(1 K) < e_cgs < e(1e12 K) there, so neither
  // clamp applies.  Otherwise the bisection, first with a margin bracket.
  double L = 0.0;
  if (std::isfinite(e_cgs) && e_cgs > 0.0 && secant_log_T(c, e_cgs, L))
    return std::exp(L);
  double La = 0.0, Lb = 0.0;
  if (std::isfinite(e_cgs) && e_cgs > 0.0 && margin_bracket(c, e_cgs, La, Lb))
    return bisect_log_T(c, e_cgs, La, Lb);

  return temperature_from_specific_internal_energy_bisection(c, e_cgs);
}

double temperature_from_specific_internal_energy_bisection(
    const EosComposition& c, double e_cgs) noexcept {
  if (std::isnan(e_cgs)) return e_cgs;
  const double e_lo = energy_at_log_T(c, std::log(1.0));
  const double e_hi = energy_at_log_T(c, std::log(1.0e12));
  if (e_cgs <= e_lo) return std::exp(std::log(1.0));
  if (e_cgs >= e_hi) return std::exp(std::log(1.0e12));
  return bisect_log_T(c, e_cgs, -std::numeric_limits<double>::infinity(),
                      std::numeric_limits<double>::infinity());
}

template <class Model>
double temperature_from_specific_internal_energy(
    const std::array<double, Model::N_sp>& y, double e_cgs) noexcept {
  return temperature_from_specific_internal_energy(
      make_eos_composition<Model>(y), e_cgs);
}

template double temperature_from_specific_internal_energy<Nakauchi2019>(
    const std::array<double, Nakauchi2019::N_sp>&, double) noexcept;
template double temperature_from_specific_internal_energy<Nakauchi2019_Minimal>(
    const std::array<double, Nakauchi2019_Minimal::N_sp>&, double) noexcept;
template double temperature_from_specific_internal_energy<Nakauchi2021>(
    const std::array<double, Nakauchi2021::N_sp>&, double) noexcept;
template double temperature_from_specific_internal_energy<Nakauchi2021_Minimal>(
    const std::array<double, Nakauchi2021_Minimal::N_sp>&, double) noexcept;

}  // namespace arche::thermo
