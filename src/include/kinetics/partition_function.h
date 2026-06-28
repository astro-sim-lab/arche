// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

#include "core/state.h"
#include "kinetics/topology.h"

// ---------------------------------------------------------------------------
// partition_function.h — generic partition-function primitives.
//
// Interpolators and H2 internal-energy / partition data shared by every model,
// plus the PfProvider abstraction that expresses a single species' partition
// function uniformly as one of four kinds: a constant, an analytic fit, a
// 32-point BC16 table, or a tabulated (T, Q) numeric table.  Model layers build
// per-species provider tables on top of this primitive set (see
// models/<model>/partition_function_*.h), so the physical forms live in one
// place and every model reads them through the same evaluator.
//
// Sources:
//   P09   = Pagano+09
//   BC16  = Barklem & Collet 2016 (diatomic, 32-point hardcoded table)
//   PJ16  = Popovas & Jorgensen 2016 (H2 polynomial fits)
//   BFM89 = Borysow, Frommhold & Moraldi 1989 (H2 energy levels)
//   HITRAN / ExoMol — tabulated data read as numeric tables.
// ---------------------------------------------------------------------------

namespace arche {
namespace detail {

// Linear interpolation on a sorted (T, value) vector.
// Clamps to boundary values outside the tabulated range.
inline double pf_interp(const std::vector<std::pair<double, double>>& tab,
                        double T_K) {
  if (tab.empty()) return 1.0;
  if (T_K <= tab.front().first) return tab.front().second;
  if (T_K >= tab.back().first) return tab.back().second;
  // Find first entry with T >= T_K
  auto it = std::lower_bound(
      tab.begin(), tab.end(), std::pair<double, double>{T_K, 0.0},
      [](const auto& a, const auto& b) { return a.first < b.first; });
  auto lo = std::prev(it);
  double t = (T_K - lo->first) / (it->first - lo->first);
  return (1.0 - t) * lo->second + t * it->second;
}

// BC16 32-point temperature grid (shared by H2+, HeH+, Li, LiH)
constexpr std::array<double, 32> kTBC = {
    1.0,    1.3,    1.7,    2.0,    3.0,    5.0,    7.0,    10.0,
    15.0,   20.0,   30.0,   50.0,   70.0,   100.0,  130.0,  170.0,
    200.0,  250.0,  300.0,  500.0,  700.0,  1000.0, 1500.0, 2000.0,
    3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0};

// Linear interpolation on the 32-point BC16 table
inline double bc16_interp(const std::array<double, 32>& ytab, double T_K) {
  int ms = 31;
  for (int i = 0; i < 32; ++i) {
    if (T_K <= kTBC[i]) {
      ms = i;
      break;
    }
  }
  if (ms == 0) ms = 1;
  double t = (T_K - kTBC[ms - 1]) / (kTBC[ms] - kTBC[ms - 1]);
  return (1.0 - t) * ytab[ms - 1] + t * ytab[ms];
}

// H2 vibrational+rotational energy levels, BFM 1989 (ApJ 336,495) eq.(A11)
// Returns E[K] for vibrational quantum number iv and rotational quantum J.
inline double E_H2_BFM(int iv, int J) {
  constexpr double Ev0 = 0.38496, Ev1 = -0.04609, Ev2 = 0.00178, Ev3 = -7.0e-5,
                   Ev4 = 2.9511e-6;
  constexpr double Bv0 = 54.438, Bv1 = 4.6063, Bv2 = -2.0050, Bv3 = 0.19260,
                   Bv4 = -6.2953e-3;
  constexpr double Dv0 = -0.72593, Dv1 = 5.9990, Dv2 = -1.5187, Dv3 = 0.12721,
                   Dv4 = -3.3391e-3;
  constexpr double Fv0 = -12.662, Fv1 = 20.047, Fv2 = -4.7873, Fv3 = 0.36900,
                   Fv4 = -9.003e-3;
  constexpr double Gv0 = -24.006, Gv1 = 33.989, Gv2 = -7.6841, Gv3 = 0.52413,
                   Gv4 = -1.1297e-2;
  constexpr double Hv0 = -22.384, Hv1 = 31.150, Hv2 = -6.6139, Hv3 = 0.39226,
                   Hv4 = -9.496e-3;
  constexpr double Ov0 = -10.541, Ov1 = 14.746, Ov2 = -2.9476, Ov3 = 0.16016,
                   Ov4 = -6.5005e-3;
  constexpr double Pv0 = -2.0021, Pv1 = 2.8400, Pv2 = -0.54654, Pv3 = 0.031636,
                   Pv4 = -2.4398e-3;

  double vv = iv + 0.5;
  double vv2 = vv * vv, vv3 = vv2 * vv, vv4 = vv3 * vv;
  double Ev = Ev0 + Ev1 * vv + Ev2 * vv2 + Ev3 * vv3 + Ev4 * vv4;
  double Bv = Bv0 + Bv1 * vv + Bv2 * vv2 + Bv3 * vv3 + Bv4 * vv4;
  double Dv = Dv0 + Dv1 * vv + Dv2 * vv2 + Dv3 * vv3 + Dv4 * vv4;
  double Fv = Fv0 + Fv1 * vv + Fv2 * vv2 + Fv3 * vv3 + Fv4 * vv4;
  double Gv = Gv0 + Gv1 * vv + Gv2 * vv2 + Gv3 * vv3 + Gv4 * vv4;
  double Hv = Hv0 + Hv1 * vv + Hv2 * vv2 + Hv3 * vv3 + Hv4 * vv4;
  double Ov = Ov0 + Ov1 * vv + Ov2 * vv2 + Ov3 * vv3 + Ov4 * vv4;
  double Pv = Pv0 + Pv1 * vv + Pv2 * vv2 + Pv3 * vv3 + Pv4 * vv4;

  double rrot = static_cast<double>(J * (J + 1));
  return 1.43879 *
         (-Ev * 1e5 + Bv * rrot - Dv * 1e-2 * rrot * rrot +
          Fv * 1e-5 * rrot * rrot * rrot -
          Gv * 1e-8 * rrot * rrot * rrot * rrot +
          Hv * 1e-11 * rrot * rrot * rrot * rrot * rrot -
          Ov * 1e-14 * rrot * rrot * rrot * rrot * rrot * rrot +
          Pv * 1e-17 * rrot * rrot * rrot * rrot * rrot * rrot * rrot);
}

// H2 partition function (equilibrium ortho/para 1:3).
inline double H2_pf(double T_K) {
  constexpr int iv_max = 5, j_max = 25;
  double ET[iv_max + 1][j_max + 1];
  for (int iv = 0; iv <= iv_max; ++iv)
    for (int j = 0; j <= j_max; ++j) ET[iv][j] = E_H2_BFM(iv, j);

  const double ET_diss = 5.196e4;
  double z_p = 0.0, z_o = 0.0;
  for (int iv = 0; iv <= iv_max; ++iv)
    for (int j = 0; j <= j_max; ++j)
      if (ET[iv][j] - ET[0][0] <= ET_diss) {
        if (j % 2 == 0)
          z_p += (2 * j + 1) * std::exp(-(ET[iv][j] - ET[0][0]) / T_K);
        else
          z_o += (2 * j + 1) * std::exp(-(ET[iv][j] - ET[0][1]) / T_K);
      }
  return 0.25 * z_p + 0.75 * z_o;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// PfProvider — one species' partition function as a uniform descriptor.
//
// A species' partition function takes one of four shapes:
//   Constant     — temperature-independent value (e.g. e-, H+, He).
//   Analytic     — a closed-form fit f(T_K) (e.g. the H multi-level sum, the
//                  piecewise H2 polynomial).
//   Bc16         — interpolation on a 32-point BC16 table over the kTBC grid.
//   NumericTable — interpolation on a per-species tabulated (T, Q) vector,
//                  scaled by `divisor` (nuclear-spin degeneracy removal) and
//                  optionally floored at 1.0.
//
// pf_eval() evaluates a descriptor; the four factory helpers below build one.
// The operations match the historical inline evaluators exactly (same order),
// so values are bit-for-bit identical to the previous slot-indexed path.
// ---------------------------------------------------------------------------
enum class PfKind { Constant, Analytic, Bc16, NumericTable };

struct PfProvider {
  PfKind kind = PfKind::Constant;
  double constant = 1.0;                         // Constant
  double (*analytic)(double T_K) = nullptr;      // Analytic
  const std::array<double, 32>* bc16 = nullptr;  // Bc16 (uses kTBC grid)
  double divisor = 1.0;                          // NumericTable scaling
  bool floor_one = false;                        // NumericTable floor at 1.0
};

inline PfProvider pf_const(double c) {
  PfProvider p;
  p.kind = PfKind::Constant;
  p.constant = c;
  return p;
}
inline PfProvider pf_analytic(double (*f)(double)) {
  PfProvider p;
  p.kind = PfKind::Analytic;
  p.analytic = f;
  return p;
}
inline PfProvider pf_bc16(const std::array<double, 32>* table) {
  PfProvider p;
  p.kind = PfKind::Bc16;
  p.bc16 = table;
  return p;
}
inline PfProvider pf_numeric(double divisor, bool floor_one = false) {
  PfProvider p;
  p.kind = PfKind::NumericTable;
  p.divisor = divisor;
  p.floor_one = floor_one;
  return p;
}

// Evaluate one species' partition function.  `table` is the per-species
// numeric (T, Q) data; it is used only by the NumericTable kind and ignored
// otherwise.
inline double pf_eval(const PfProvider& p, double T_K,
                      const std::vector<std::pair<double, double>>& table) {
  switch (p.kind) {
    case PfKind::Constant:
      return p.constant;
    case PfKind::Analytic:
      return p.analytic(T_K);
    case PfKind::Bc16:
      return detail::bc16_interp(*p.bc16, T_K);
    case PfKind::NumericTable: {
      double raw = detail::pf_interp(table, T_K) / p.divisor;
      return p.floor_one ? std::max(1.0, raw) : raw;
    }
  }
  return 1.0;
}

}  // namespace arche
