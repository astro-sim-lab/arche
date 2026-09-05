// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
//
// collapse_defaults.h — shared constants and types for one-zone collapse apps
//
// Included by both collapse_primordial and collapse_metal_grain.
// Network-specific constants (kNSp, kOutputStride, kItMax, grain params, etc.)
// remain in each app's main.cc.
//
#pragma once

#include <cmath>
#include <cstdio>
#include <string>

#include "solve/chemistry.h"

namespace collapse_defaults {

// ─── Physical constants (CGS; aliases of arche::phys, see core/state.h) ─────
constexpr double kKB = arche::phys::k_B;  // [erg/K]
constexpr double kMp = arche::phys::m_p;  // [g]
constexpr double kPi = arche::phys::pi;
constexpr double kGGrav = arche::phys::G;         // [cm^3 g^-1 s^-2]
constexpr double kSigmaB = arche::phys::sigma_B;  // [erg cm^-2 s^-1 K^-4]

// ─── CMB ─────────────────────────────────────────────────────────────────────
constexpr double kTCMB0 = 2.725;  // CMB temperature at z=0 [K]

// ─── CR attenuation ──────────────────────────────────────────────────────────
// Primary attenuation column density [g cm^-2]
constexpr double kCrAttenuColDens = 96.0;

// ─── Integration control defaults ────────────────────────────────────────────
// (overridable at runtime via environment variables)
constexpr double kDtFactor = 1.0e-3;  // timestep fraction of min(t_cool, t_eff)
constexpr double kDtFactorChem =
    1.0e50;  // timestep fraction of t_chem (1e50 = disabled by default)
constexpr double kDtFactorInit =
    1.0e-8;                          // timestep fraction during initial phase
constexpr int kNInitSteps = 10;      // number of initial short-timestep steps
constexpr double kXnHStop = 1.0e23;  // density ceiling [cm^-3]
constexpr double kTHighStop =
    1.0e5;  // temperature ceiling (validated T_max) [K]
// HDF5 output stride: write every N-th step.  Shared by both models.
constexpr int kOutputStride = 10;
// Maximum integration steps before MaxIter exit.  Shared by both models.
constexpr int kItMax = 1000000;

// ─── Initial conditions / scenario defaults ──────────────────────────────────
// (overridable at runtime via environment variables)
constexpr double kXnH0 = 0.1;     // initial H number density [cm^-3]
constexpr double kTK0 = 1.0e2;    // initial gas temperature [K]
constexpr double kYe0 = 1.0e-4;   // initial electron / H+ fraction
constexpr double kYH2 = 6.0e-7;   // initial H2 fraction
constexpr double kYHD = 4.0e-10;  // initial HD fraction
constexpr double kZeta0 = 0.0;    // CR ionization rate [s^-1]
constexpr double kFFRet = 1.0;    // free-fall retardation factor
constexpr double kJLW21 =
    0.0;  // Lyman-Werner intensity [10^-21 erg/s/cm^2/Hz/sr]
constexpr double kZred = 0.0;    // cosmological redshift
constexpr double kZmetal = 0.0;  // metallicity [Z_sun] (metal_grain only)

// ─── Subcycle control (ARCHE_SUBCYCLE) ───────────────────────────────────────
// Maximum bisection depth for chemcool subcycling.
// dt is halved up to kSubcycleMaxLevel times (2^6 = 64 sub-steps at most).
constexpr int kSubcycleMaxLevel = 6;

// ─── ExitReason ──────────────────────────────────────────────────────────────
// Categorises why the time integration loop terminated.  exit_code values
// written to stdout and HDF5 attributes match enum order.  SolverFailed is
// always defined (even without ARCHE_SUBCYCLE) to keep enum values stable
// across build configurations.
enum class ExitReason {
  Normal,    // 0: nH >= kXnHStop     (density ceiling reached — expected end)
  MaxIter,   // 1: it  >  max_iter     (step limit hit before density ceiling)
  HighTemp,  // 2: T_K > kTHighStop    (T ceiling reached — expected end)
  NegEnergy,     // 3: e  <= 0             (internal energy went non-positive)
  NonFinite,     // 4: NaN/Inf in nH or T_K
  SolverFailed,  // 5: subcycle max depth exceeded without NR convergence
                 // (ARCHE_SUBCYCLE)
};

// Normalise a numeric string to the canonical tag form used in filenames.
// Rule: 0 → "0"; already in scientific notation (contains 'e'/'E') → '.'→'p';
// otherwise convert to <mantissa>e<±exp> with mantissa in [1,10).
inline std::string make_tag(const char* raw, double val) {
  if (val == 0.0) return "0";
  std::string s(raw);
  if (s.find('e') != std::string::npos || s.find('E') != std::string::npos) {
    for (char& c : s)
      if (c == '.') c = 'p';
    return s;
  }
  double x = val;
  int exp = 0;
  if (x >= 10.0) {
    while (x >= 10.0) {
      x /= 10.0;
      exp++;
    }
  } else if (x < 1.0) {
    while (x < 1.0) {
      x *= 10.0;
      exp--;
    }
  }
  x = std::round(x * 1.0e9) / 1.0e9;
  char buf[64];
  if (x == static_cast<int>(x)) {
    std::snprintf(buf, sizeof(buf), "%de%+d", static_cast<int>(x + 0.5), exp);
  } else {
    char mbuf[32];
    std::snprintf(mbuf, sizeof(mbuf), "%g", x);
    std::string m(mbuf);
    for (char& c : m)
      if (c == '.') c = 'p';
    std::snprintf(buf, sizeof(buf), "%se%+d", m.c_str(), exp);
  }
  return std::string(buf);
}

}  // namespace collapse_defaults
