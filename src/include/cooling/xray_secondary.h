// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// xray_secondary.h — X-ray secondary ionization and heating
//
// Implements the Sugimura (2026) fitting formulas for the fraction of
// primary photoelectron energy deposited into:
//   - HI secondary ionization:  secion_HI2nd()
//   - HeI secondary ionization: secion_HeI2nd()
//   - Gas heating:              secion_heat()
//
// Fitting is based on Shull & van Steenberg (1985), ApJ, 298, 268 (Table 1),
// with improvements following Ricotti, Gnedin & Shull (2002), ApJ, 575, 33.
//
// Photoionization cross sections: Verner et al. (1996), ApJ, 465, 487, Table 1.
//   sigma_HI_Verner96()   — H  I (1s)
//   sigma_HeI_Verner96()  — He I (1s)
//   sigma_HeII_Verner96() — He II (1s, H-like)
//
// Convention for xe (electron fraction):
//   xe = y[e] / (1 + yHe)   where y[e] = ne/nH, yHe = 8.33e-2
//   This differs slightly from Sugimura's xe = ne/(nH+nHe), but the
//   relative error is < ~8% since yHe << 1.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <cassert>
#include <cmath>

namespace arche {
namespace xray {

// ---------------------------------------------------------------------------
// Verner et al. (1996), ApJ, 465, 487 — photoionization cross sections [cm^2]
// Formula: σ(E) = σ₀ · F(y)  [Mb],  1 Mb = 10^-18 cm^2
//   x  = E/E₀ - y₀
//   y  = sqrt(x^2 + y₁^2)
//   F(y) = [(x-1)^2 + yw^2] · y^(0.5P-5.5) · (1 + sqrt(y/ya))^(-P)
// ---------------------------------------------------------------------------

// H I 1s — Verner et al. (1996) Table 1
// y₀=0, y₁=0, yw=0  → simplifies to x=E/E₀, y=x
static inline double sigma_HI_Verner96(double E_eV) {
  if (E_eV < 13.6) return 0.0;
  constexpr double E0 = 4.298e-1, s0 = 5.475e4, ya = 32.88, P = 2.963;
  const double x = E_eV / E0;             // y₀=0
  const double y = x;                     // y₁=0
  const double F = (x - 1.0) * (x - 1.0)  // yw=0
                   * std::pow(y, 0.5 * P - 5.5) *
                   std::pow(1.0 + std::sqrt(y / ya), -P);
  return s0 * F * 1.0e-18;  // [cm^2]
}

// He I 1s — Verner et al. (1996) Table 1
static inline double sigma_HeI_Verner96(double E_eV) {
  if (E_eV < 24.6) return 0.0;
  constexpr double E0 = 13.61, s0 = 949.2, ya = 1.469, P = 3.188;
  constexpr double yw = 2.039, y0 = 0.4434, y1 = 2.136;
  const double x = E_eV / E0 - y0;
  const double y = std::sqrt(x * x + y1 * y1);
  const double F = ((x - 1.0) * (x - 1.0) + yw * yw) *
                   std::pow(y, 0.5 * P - 5.5) *
                   std::pow(1.0 + std::sqrt(y / ya), -P);
  return s0 * F * 1.0e-18;  // [cm^2]
}

// He II 1s (H-like, Z=2) — Verner et al. (1996) Table 1
// Same form as H I with y₀=0, y₁=0, yw=0
static inline double sigma_HeII_Verner96(double E_eV) {
  if (E_eV < 54.42) return 0.0;
  constexpr double E0 = 1.720, s0 = 1.369e4, ya = 32.88, P = 2.963;
  const double x = E_eV / E0;
  const double F = (x - 1.0) * (x - 1.0) * std::pow(x, 0.5 * P - 5.5) *
                   std::pow(1.0 + std::sqrt(x / ya), -P);
  return s0 * F * 1.0e-18;  // [cm^2]
}

// ---------------------------------------------------------------------------
// secion_HI2nd(E0_eV, xe)
//
// Returns Phi_HI: number of secondary HI ionizations per primary
// photoionization event.
//
// Fitting formula: Sugimura (2026), based on Shull & van Steenberg (1985)
//   and improved by Ricotti et al. (2002).
//
// Parameters:
//   E0_eV — primary photoelectron energy [eV]  (= E_X - E_ion of parent)
//   xe    — electron fraction = ne/(nH + nHe)  [0, 1]
//
// Returns 0 if E0_eV < 28 eV (threshold from Shull & van Steenberg 1985).
// ---------------------------------------------------------------------------
static inline double secion_HI2nd(double E0_eV, double xe) {
  assert(E0_eV >= 0.0);
  assert(xe >= 0.0 && xe <= 1.0);

  constexpr double E_ref = 28.0;      // [eV] threshold
  constexpr double E_ion_HI = 13.60;  // [eV]

  if (E0_eV < E_ref) return 0.0;

  constexpr double C1 = -0.5517, C2 = 0.9904, C3 = -0.006490;
  constexpr double a1 = 0.3692, a2 = -0.1220, a3 = -0.9775;
  constexpr double b1 = 1.476, b2 = 2.787, b3 = -0.6476;

  const double u = E0_eV / E_ref;
  const double C = C1 + C2 * std::pow(u, C3);
  const double a = a1 + a2 * std::pow(u, a3);
  const double b = b1 + b2 * std::pow(u, b3);

  const double Efrac = C * std::pow(1.0 - std::pow(xe, a), b);
  return std::max(0.0, Efrac * E0_eV / E_ion_HI);
}

// ---------------------------------------------------------------------------
// secion_HeI2nd(E0_eV, xe)
//
// Returns Phi_HeI: number of secondary HeI ionizations per primary
// photoionization event.
//
// Same fitting form as secion_HI2nd but with C1 = -0.8986 and
// E_ion_HeI = 24.59 eV.  (C2, C3, a*, b* parameters are shared.)
// ---------------------------------------------------------------------------
static inline double secion_HeI2nd(double E0_eV, double xe) {
  assert(E0_eV >= 0.0);
  assert(xe >= 0.0 && xe <= 1.0);

  constexpr double E_ref = 28.0;       // [eV] threshold (same as HI)
  constexpr double E_ion_HeI = 24.59;  // [eV]

  if (E0_eV < E_ref) return 0.0;

  constexpr double C1 = -0.8986, C2 = 0.9904, C3 = -0.006490;
  constexpr double a1 = 0.3692, a2 = -0.1220, a3 = -0.9775;
  constexpr double b1 = 1.476, b2 = 2.787, b3 = -0.6476;

  const double u = E0_eV / E_ref;
  const double C = C1 + C2 * std::pow(u, C3);
  const double a = a1 + a2 * std::pow(u, a3);
  const double b = b1 + b2 * std::pow(u, b3);

  const double Efrac = C * std::pow(1.0 - std::pow(xe, a), b);
  return std::max(0.0, Efrac * E0_eV / E_ion_HeI);
}

// ---------------------------------------------------------------------------
// secion_heat(E0_eV, xe)
//
// Returns Eh: heating energy deposited per primary photoionization event [eV].
//
// Fitting formula: Sugimura (2026), based on Shull & van Steenberg (1985).
//
// Parameters:
//   E0_eV — primary photoelectron energy [eV]
//   xe    — electron fraction [0, 1]
//
// If E0_eV < 11 eV (threshold), all energy is assumed to go into heating:
//   Eh = E0_eV.
// ---------------------------------------------------------------------------
static inline double secion_heat(double E0_eV, double xe) {
  assert(E0_eV >= 0.0);
  assert(xe >= 0.0 && xe <= 1.0);

  constexpr double E_ref = 11.0;  // [eV] threshold

  if (E0_eV < E_ref) return E0_eV;

  constexpr double a1 = 0.2364, a2 = 0.07122, a3 = -0.3183;
  constexpr double b1 = 1.091, b2 = 3.320, b3 = -0.7662;

  const double u = E0_eV / E_ref;
  const double a = a1 + a2 * std::pow(u, a3);
  const double b = b1 + b2 * std::pow(u, b3);

  const double Efrac = 1.0 - std::pow(1.0 - std::pow(xe, a), b);
  return Efrac * E0_eV;  // [eV]
}

}  // namespace xray
}  // namespace arche
