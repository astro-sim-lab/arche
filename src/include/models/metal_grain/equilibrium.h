// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// equilibrium.h — equichem_metal helper for the metal_grain Saha branch
//
// detail::equichem_metal performs a 4D Newton–Raphson on (y_e, y_H, y_C, y_O),
// filling y[0..62] under the input cell's own gas-phase element budgets.
// ---------------------------------------------------------------------------
#include <algorithm>
#include <cmath>
#include <limits>

#include "core/saha_newton.h"  // saha_newton_log
#include "kinetics/rates.h"
#include "models/metal_grain/minimal.h"  // metal_grain_minimal::Sp / MinimalTable
#include "models/metal_grain/partition_function_metal.h"  // eval_metal_partition_functions

namespace arche {

// ─────────────────────────────────────────────────────────────────────────────
// detail::equichem_metal
//
// Saha equilibrium for the metal_grain network.
// 4D Newton–Raphson on (y_e, y_H, y_C, y_O).
// Fills y[0..62]; leaves y[63..94] unchanged.
// ─────────────────────────────────────────────────────────────────────────────
namespace detail {

// metal_saha_solve — the 4D (y_e, y_H, y_C, y_O) Saha root through core/saha_newton.h,
// from the cell's own abundances, then, if that fails, from a start that depends on T
// and n_H and the cell's budgets only (y_H from the H/H2 balance with ions neglected,
// y_H + 2 K_H2 y_H^2 = y_H,tot, y_e from H+, D+ and Li+, y_C and y_O at their element
// budgets). If both fail the active unknowns become NaN,
// so the integrator's composition check (saha_at) rejects the step. The log-space
// iteration keeps y_C and y_O positive: a floored step max(x + d, 1e-30) can pin them at
// the floor, and the full network then breaks the C or O budget by up to 20 orders below
// ~2800 K.
template <class Residual>
void metal_saha_solve(std::array<double, 4>& x,
                      const std::array<bool, 4>& active, double K_Hp,
                      double K_H2, double ion_sum_per_ye, double yH_tot,
                      double XC, double XO, Residual&& residual) {
  if (saha_newton_log<4>(x, active, residual)) return;
  const double h = 2.0 * yH_tot / (1.0 + std::sqrt(1.0 + 8.0 * K_H2 * yH_tot));
  std::array<double, 4> x2 = {std::sqrt(K_Hp * h + ion_sum_per_ye), h,
                              active[2] ? XC : x[2], active[3] ? XO : x[3]};
  if (saha_newton_log<4>(x2, active, residual)) {
    x = x2;
    return;
  }
  for (int j = 0; j < 4; ++j)
    if (active[j]) x[j] = std::numeric_limits<double>::quiet_NaN();
}

//
// The element budgets and the charge target are the input cell's own, over the
// gas-phase species this solver writes (0..62): the grain charge states, the
// ice mantles (63..88) and the dust reservoirs (89..94) are left as they are, so
// their elements and charge stay
// where they were and the step conserves every element and the total charge. The
// solar composition scaled by Z_metal holds only where nothing has moved the
// gas-phase budgets -- not with ice or grain charge present, nor after advection
// in 3D.
inline void equichem_metal(
    double nH, double T_K, std::array<double, metal_grain::N_sp>& y,
    const ReactionTable<metal_grain::N_sp, metal_grain::N_react>& tbl) {
  using Sp = metal_grain::Sp;  // y[Sp::H] etc. (Sp::name == species index)
  constexpr double k_B = phys::k_B;
  constexpr double h_P = phys::h_P;
  constexpr double pi = phys::pi;
  // Gas-phase budgets of the elements with a closed-form sequence (H, C and O,
  // the Newton unknowns, are summed below with their conservation residuals).
  const double yHe = y[Sp::He] + y[Sp::Hep] + y[Sp::Hepp] + y[Sp::HeHp];
  const double yD = y[Sp::D] + y[Sp::HD] + y[Sp::Dp] + y[Sp::HDp] + y[Sp::Dm];
  const double yLi = y[Sp::Li] + y[Sp::LiH] + y[Sp::Lip] + y[Sp::Lim] +
                     y[Sp::LiHp] + y[Sp::Lipp] + y[Sp::Lippp];
  const double XKa = y[Sp::K] + y[Sp::Kp];
  const double XNa = y[Sp::Na] + y[Sp::Nap];
  const double XMg = y[Sp::Mg] + y[Sp::Mgp];
  // Charge the untouched grains carry; the gas phase balances it.
  const double q_grain =
      y[Sp::Grp] + 2.0 * y[Sp::Gr2p] - y[Sp::Grm] - 2.0 * y[Sp::Gr2m];

  // ── Partition functions (0-based species + 3 sentinel slots) ─────────────
  std::array<double, metal_grain::N_sp + 3> pf{};
  eval_metal_partition_functions<metal_grain::N_sp, metal_grain::N_react>(
      T_K, tbl, pf);

  // ── Build Keqb[0..1199] indexed by (reaction_number - 1) ─────────────────
  std::array<double, metal_grain::N_react> Keqb{};
  Keqb.fill(0.0);
  {
    const double lnT32 = 1.5 * std::log(2.0 * pi * k_B * T_K / (h_P * h_P));
    for (int i = 0; i < tbl.n_saha; ++i) {
      int num = tbl.saha[i].num;
      int r1 = tbl.saha[i].reactants[0];
      int r2 = tbl.saha[i].reactants[1];
      int p1 = tbl.saha[i].products[0];
      int p2 = tbl.saha[i].products[1];
      int nr = tbl.saha[i].n_reactants;
      int np = tbl.saha[i].n_products;
      double Cm = tbl.saha[i].Cmass;
      double dE = tbl.saha[i].delE;

      double lnCpf = std::log(pf[r1]) + std::log(pf[r2]) - std::log(pf[p1]);
      if (np == 2 && p2 != ReactionTable<metal_grain::N_sp,
                                         metal_grain::N_react>::IDX_PHOTON)
        lnCpf -= std::log(pf[p2]);

      double lnKeqb = static_cast<double>(nr - np) * lnT32 + std::log(Cm) +
                      lnCpf - dE / (k_B * T_K);

      if (num >= 1 && num <= metal_grain::N_react)
        Keqb[num - 1] = std::exp(lnKeqb);
    }
  }

  // ── Equilibrium ratio constants (1-based Keqb(N) → 0-based Keqb[N-1]) ────
  // H sequence
  double K_Hp = Keqb[1] / nH;
  double K_Hm = nH / Keqb[6];
  double K_H2p = Keqb[1] / Keqb[8];
  double K_H2 = nH / Keqb[6] / Keqb[7];
  double K_H3p = nH * Keqb[1] / (Keqb[6] * Keqb[7] * Keqb[8] * Keqb[25]);
  // He sequence
  double K_Hep = Keqb[3] / nH;
  double K_He2p = Keqb[3] * Keqb[5] / (nH * nH);
  double K_HeHp = nH / Keqb[640];
  // D sequence
  double K_Dp = Keqb[50] / nH;
  double K_HD = nH / Keqb[53];
  double K_HDp = nH / Keqb[59];
  double K_Dm = nH / Keqb[62];
  // Li sequence
  double K_Lip = Keqb[800] / nH;
  double K_Li2p = K_Lip * Keqb[820] / nH;
  double K_Li3p = K_Li2p * Keqb[821] / nH;
  double K_LiH = nH / Keqb[817];
  double K_LiHp = nH / Keqb[812];
  double K_Lim = nH / Keqb[803];

  // ── fill_y: write yy[0..62] from (y_e, y_H, y_C, y_O) ──────────────────
  // Writes to the supplied array (a scratch array during the NR loop, the
  // real output `y` only at the end) instead of mutating a captured `y`.
  auto fill_y = [&](auto& yy, double ye, double yh, double yc, double yo) {
    double yHp = K_Hp * yh / ye;
    double yHm = K_Hm * yh * ye;
    double yH2p = K_H2p * yh * yh / ye;
    double yH2 = K_H2 * yh * yh;
    double yH3p = K_H3p * yh * yh * yh / ye;

    yy[Sp::H] = yh;      // H
    yy[Sp::H2] = yH2;    // H2
    yy[Sp::e] = ye;      // e-
    yy[Sp::Hp] = yHp;    // H+
    yy[Sp::H2p] = yH2p;  // H2+
    yy[Sp::H3p] = yH3p;  // H3+
    yy[Sp::Hm] = yHm;    // H-

    // He sequence
    yy[Sp::He] = yHe / (1.0 + K_Hep / ye + K_He2p / (ye * ye) + K_HeHp * yHp);
    yy[Sp::Hep] = K_Hep * yy[Sp::He] / ye;
    yy[Sp::Hepp] = K_He2p * yy[Sp::He] / (ye * ye);
    yy[Sp::HeHp] = K_HeHp * yHp * yy[Sp::He];

    // D sequence
    yy[Sp::D] = yD / (1.0 + K_HD * yh + K_Dp / ye + K_HDp * yHp + K_Dm * ye);
    yy[Sp::HD] = K_HD * yy[Sp::D] * yh;
    yy[Sp::Dp] = K_Dp * yy[Sp::D] / ye;
    yy[Sp::HDp] = K_HDp * yy[Sp::D] * yHp;
    yy[Sp::Dm] = K_Dm * yy[Sp::D] * ye;

    // Li sequence  (1-based y(51..57) → 0-based yy[50..56])
    yy[Sp::Li] =
        yLi / (1.0 + K_LiH * yh + K_Lim * ye + K_LiHp * yHp + K_Lip / ye +
               K_Li2p / (ye * ye) + K_Li3p / (ye * ye * ye));
    yy[Sp::LiH] = K_LiH * yy[Sp::Li] * yh;
    yy[Sp::Lip] = K_Lip * yy[Sp::Li] / ye;
    yy[Sp::Lim] = K_Lim * yy[Sp::Li] * ye;
    yy[Sp::LiHp] = K_LiHp * yy[Sp::Li] * yHp;
    yy[Sp::Lipp] = K_Li2p * yy[Sp::Li] / (ye * ye);
    yy[Sp::Lippp] = K_Li3p * yy[Sp::Li] / (ye * ye * ye);

    // Alkali ion sequences
    yy[Sp::Kp] = XKa / (1.0 + nH * ye / Keqb[700]);   // K+
    yy[Sp::Nap] = XNa / (1.0 + nH * ye / Keqb[702]);  // Na+
    yy[Sp::Mgp] = XMg / (1.0 + nH * ye / Keqb[717]);  // Mg+

    // Carbon sequence  (1-based y(17..29) → 0-based yy[16..28])
    yy[Sp::C] = yc;
    yy[Sp::C2] = yc * yc * nH / Keqb[185];
    yy[Sp::CH] = yh * yc * nH / Keqb[184];
    yy[Sp::CH2] = yy[Sp::H2] * yc * nH / Keqb[537];
    yy[Sp::CH3] = yy[Sp::H2] * yy[Sp::CH] * nH / Keqb[642];
    yy[Sp::CH4] = yy[Sp::H2] * yy[Sp::CH3] / (yh * Keqb[135]);
    yy[Sp::Cp] = (yc / ye) * Keqb[506] / nH;
    yy[Sp::C2p] = yc * yy[Sp::Cp] * nH / Keqb[643];
    yy[Sp::CHp] = yh * yy[Sp::Cp] * nH / Keqb[294];
    yy[Sp::CH2p] = yy[Sp::H2] * yy[Sp::Cp] * nH / Keqb[297];
    yy[Sp::CH3p] = (yy[Sp::CH3] / ye) * Keqb[513] / nH;
    yy[Sp::CH4p] = yy[Sp::Hp] * yy[Sp::CH4] / (yh * Keqb[212]);
    yy[Sp::CH5p] = yy[Sp::H2] * yy[Sp::CH3p] * nH / Keqb[339];

    // Oxygen sequence  (1-based y(30..50) → 0-based yy[29..49])
    yy[Sp::O] = yo;
    yy[Sp::O2] = yo * yo * nH / Keqb[193];
    yy[Sp::OH] = yh * yo * nH / Keqb[192];
    yy[Sp::CO] = yc * yo * nH / Keqb[186];
    yy[Sp::H2O] = yh * yy[Sp::OH] * nH / Keqb[641];
    yy[Sp::HCO] = Keqb[183] * yy[Sp::H2] * yy[Sp::CO] / yh;
    yy[Sp::O2H] = Keqb[114] * yo * yy[Sp::H2O] / yh;
    yy[Sp::CO2] = yo * yy[Sp::HCO] / (Keqb[200] * yh);
    yy[Sp::H2CO] = Keqb[111] * yy[Sp::H2] * yy[Sp::HCO] / yh;
    yy[Sp::H2O2] = Keqb[601] * yy[Sp::OH] * yy[Sp::H2O] / yh;
    yy[Sp::Op] = (yo / ye) * Keqb[514] / nH;
    yy[Sp::O2p] = Keqb[533] * yo * yo / ye;
    yy[Sp::OHp] = Keqb[517] * yh * yo / ye;
    yy[Sp::COp] = yc * yy[Sp::Op] * nH / Keqb[644];
    yy[Sp::H2Op] = Keqb[520] * yh * yy[Sp::OH] / ye;
    yy[Sp::HCOp] = Keqb[526] * yh * yy[Sp::CO] / ye;
    yy[Sp::O2Hp] = Keqb[534] * yh * yy[Sp::O2] / ye;
    yy[Sp::H3Op] = Keqb[522] * yh * yy[Sp::H2O] / ye;
    yy[Sp::H2COp] = (yy[Sp::H2CO] / ye) * Keqb[529] / nH;
    yy[Sp::HOCOp] = Keqb[535] * yh * yy[Sp::CO2] / ye;
    yy[Sp::H2COHp] = Keqb[532] * yh * yy[Sp::H2CO] / ye;
  };

  // ── Conservation residuals ────────────────────────────────────────────────
  auto F_cha_fn = [&](const auto& yy) -> double {
    return q_grain + yy[Sp::Hp] + yy[Sp::H2p] + yy[Sp::H3p] + yy[Sp::Hep] +
           yy[Sp::HeHp] + yy[Sp::Dp] + yy[Sp::HDp] + yy[Sp::Cp] + yy[Sp::C2p] +
           yy[Sp::CHp] + yy[Sp::CH2p] + yy[Sp::CH3p] + yy[Sp::CH4p] +
           yy[Sp::CH5p] + yy[Sp::Op] + yy[Sp::O2p] + yy[Sp::OHp] + yy[Sp::COp] +
           yy[Sp::H2Op] + yy[Sp::HCOp] + yy[Sp::O2Hp] + yy[Sp::H3Op] +
           yy[Sp::H2COp] + yy[Sp::H2COHp] + yy[Sp::Lip] + yy[Sp::LiHp] +
           yy[Sp::Kp] + yy[Sp::Nap] + yy[Sp::Mgp] +
           2.0 * (yy[Sp::Hepp] + yy[Sp::Lipp]) + 3.0 * yy[Sp::Lippp] -
           (yy[Sp::e] + yy[Sp::Hm] + yy[Sp::Dm] + yy[Sp::Lim]);
  };
  auto H_content = [&](const auto& yy) -> double {
    return yy[Sp::H] + yy[Sp::Hp] + yy[Sp::Hm] + yy[Sp::HeHp] + yy[Sp::HD] +
           yy[Sp::HDp] + yy[Sp::CH] + yy[Sp::CHp] + yy[Sp::OH] + yy[Sp::HCO] +
           yy[Sp::O2H] + yy[Sp::OHp] + yy[Sp::HCOp] + yy[Sp::O2Hp] +
           yy[Sp::HOCOp] + yy[Sp::LiH] + yy[Sp::LiHp] +
           2.0 * (yy[Sp::H2] + yy[Sp::H2p] + yy[Sp::CH2] + yy[Sp::CH2p] +
                  yy[Sp::H2O] + yy[Sp::H2CO] + yy[Sp::H2O2] + yy[Sp::H2Op] +
                  yy[Sp::H2COp]) +
           3.0 * (yy[Sp::H3p] + yy[Sp::CH3] + yy[Sp::CH3p] + yy[Sp::H3Op] +
                  yy[Sp::H2COHp]) +
           4.0 * (yy[Sp::CH4] + yy[Sp::CH4p]) + 5.0 * yy[Sp::CH5p];
  };
  auto C_content = [&](const auto& yy) -> double {
    return yy[Sp::C] + yy[Sp::CH] + yy[Sp::CH2] + yy[Sp::CH3] + yy[Sp::CH4] +
           yy[Sp::Cp] + yy[Sp::CHp] + yy[Sp::CH2p] + yy[Sp::CH3p] +
           yy[Sp::CH4p] + yy[Sp::CH5p] + yy[Sp::CO] + yy[Sp::HCO] +
           yy[Sp::CO2] + yy[Sp::H2CO] + yy[Sp::COp] + yy[Sp::HCOp] +
           yy[Sp::H2COp] + yy[Sp::HOCOp] + yy[Sp::H2COHp] +
           2.0 * (yy[Sp::C2] + yy[Sp::C2p]);
  };
  auto O_content = [&](const auto& yy) -> double {
    return yy[Sp::O] + yy[Sp::OH] + yy[Sp::CO] + yy[Sp::H2O] + yy[Sp::HCO] +
           yy[Sp::H2CO] + yy[Sp::Op] + yy[Sp::OHp] + yy[Sp::COp] +
           yy[Sp::H2Op] + yy[Sp::HCOp] + yy[Sp::H3Op] + yy[Sp::H2COp] +
           yy[Sp::H2COHp] +
           2.0 * (yy[Sp::O2] + yy[Sp::O2H] + yy[Sp::CO2] + yy[Sp::H2O2] +
                  yy[Sp::O2p] + yy[Sp::O2Hp] + yy[Sp::HOCOp]);
  };
  const double yH_tot = H_content(y);
  const double XC = C_content(y);
  const double XO = O_content(y);
  auto F_hyd_fn = [&](const auto& yy) -> double {
    return H_content(yy) - yH_tot;
  };
  auto F_car_fn = [&](const auto& yy) -> double { return C_content(yy) - XC; };
  auto F_oxy_fn = [&](const auto& yy) -> double { return O_content(yy) - XO; };

  // ── Initial guess from current y ──────────────────────────────────────────
  double y_e = y[Sp::e];
  double y_H = y[Sp::H];
  double y_C = y[Sp::C];
  double y_O = y[Sp::O];
  if (y_e <= 0.0) y_e = 1.0e-20;
  if (y_H <= 0.0) y_H = 1.0e-20;
  if (y_C <= 0.0) y_C = (XC > 0.0) ? 1.0e-20 : 0.0;
  if (y_O <= 0.0) y_O = (XO > 0.0) ? 1.0e-20 : 0.0;

  // Scratch array for the per-iteration residual/Jacobian evaluations; the
  // real output `y` is written only by the final fill_y below.
  std::array<double, metal_grain::N_sp> ytmp{};
  std::array<double, 4> x = {y_e, y_H, y_C, y_O};
  const std::array<bool, 4> active = {true, true, (y_C > 0.0 && XC > 0.0),
                                      (y_O > 0.0 && XO > 0.0)};

  // 4D Newton–Raphson on (y_e, y_H, y_C, y_O): charge / hydrogen / carbon /
  // oxygen conservation.  Carbon and oxygen columns are frozen (inactive) when
  // the element is absent (XC / XO == 0).
  metal_saha_solve(
      x, active, K_Hp, K_H2, K_Dp * yD + K_Lip * yLi, yH_tot, XC, XO,
      [&](const std::array<double, 4>& xx, std::array<double, 4>& fv) {
        fill_y(ytmp, xx[0], xx[1], xx[2], xx[3]);
        fv[0] = F_cha_fn(ytmp);
        fv[1] = F_hyd_fn(ytmp);
        fv[2] = F_car_fn(ytmp);
        fv[3] = F_oxy_fn(ytmp);
      });
  y_e = x[0];
  y_H = x[1];
  y_C = x[2];
  y_O = x[3];

  // ── Final fill at converged (y_e, y_H, y_C, y_O) — writes the real output ─
  fill_y(y, y_e, y_H, y_C, y_O);

  // ── Neutral alkali abundances ────────────────────────────────────────────
  y[Sp::K] = (Keqb[700] > 1.0e-300) ? nH * y_e * y[Sp::Kp] / Keqb[700] : 0.0;
  y[Sp::Na] = (Keqb[702] > 1.0e-300) ? nH * y_e * y[Sp::Nap] / Keqb[702] : 0.0;
  y[Sp::Mg] = (Keqb[717] > 1.0e-300) ? nH * y_e * y[Sp::Mgp] / Keqb[717] : 0.0;
}

// ─────────────────────────────────────────────────────────────────────────────
// detail::equichem_metal_minimal
//
// Compact (44-species) Saha equilibrium for the metal_grain Minimal network.
// The same 4D Newton–Raphson on (y_e, y_H, y_C, y_O) as equichem_metal, reduced
// to the 31 gas-phase species the compact network carries (indices 0..30:
// H..Mg+); the grain charge states, ice-mantle species and dust reservoirs
// (compact 31..43) are left unchanged, mirroring how the full solver leaves
// 63..94.
//
// The equilibrium constants come from the compact Saha sub-table
// (kMetalMinimalSahaKeep, 23 reactions renumbered 1..23), so Keqb is indexed by
// the keep-set position rather than the full reaction number; the K_* ratios
// below name the full reaction id each compact slot carries.  Only ions and
// neutrals present in the compact selection appear in the conservation
// residuals, so no removed species (He++, the higher C/O ions, LiH, K, Na, …)
// is ever referenced in a denominator.
//
// The partition functions are supplied by the caller (the compact species'
// partition functions, evaluated in the full layout) — see MinimalMetalSaha.
// ─────────────────────────────────────────────────────────────────────────────
inline void equichem_metal_minimal(
    double nH, double T_K, std::array<double, metal_grain_minimal::N_sp>& y,
    const metal_grain_minimal::MinimalTable& tbl,
    const std::array<double, metal_grain_minimal::N_sp + 3>& pf) {
  using Sp = metal_grain_minimal::Sp;
  constexpr double k_B = phys::k_B;
  constexpr double h_P = phys::h_P;
  constexpr double pi = phys::pi;

  // The input cell's gas-phase budgets and the charge of its untouched grains, as
  // in equichem_metal.
  const double yHe = y[Sp::He] + y[Sp::Hep];
  const double yD = y[Sp::D] + y[Sp::HD] + y[Sp::Dp];
  const double yLi = y[Sp::Li] + y[Sp::Lip];
  const double XMg = y[Sp::Mg] + y[Sp::Mgp];
  const double q_grain = y[Sp::Grp] - y[Sp::Grm] - 2.0 * y[Sp::Gr2m];

  // ── Equilibrium constants from the compact Saha sub-table.  Keqb[i] is the
  //    equilibrium constant of compact Saha reaction i+1 (keep-set position i);
  //    the comments give the full reaction id it carries. ──
  std::array<double, metal_grain_minimal::N_react> Keqb{};
  Keqb.fill(0.0);
  {
    const double lnT32 = 1.5 * std::log(2.0 * pi * k_B * T_K / (h_P * h_P));
    for (int i = 0; i < tbl.n_saha; ++i) {
      int num = tbl.saha[i].num;
      int r1 = tbl.saha[i].reactants[0];
      int r2 = tbl.saha[i].reactants[1];
      int p1 = tbl.saha[i].products[0];
      int p2 = tbl.saha[i].products[1];
      int nr = tbl.saha[i].n_reactants;
      int np = tbl.saha[i].n_products;
      double Cm = tbl.saha[i].Cmass;
      double dE = tbl.saha[i].delE;

      double lnCpf = std::log(pf[r1]) + std::log(pf[r2]) - std::log(pf[p1]);
      if (np == 2 && p2 != metal_grain_minimal::MinimalTable::IDX_PHOTON)
        lnCpf -= std::log(pf[p2]);

      double lnKeqb = static_cast<double>(nr - np) * lnT32 + std::log(Cm) +
                      lnCpf - dE / (k_B * T_K);
      if (num >= 1 && num <= metal_grain_minimal::N_react)
        Keqb[num - 1] = std::exp(lnKeqb);
    }
  }

  // clang-format off
  // Equilibrium ratio constants (compact slot index → full reaction id):
  // H sequence
  double K_Hp = Keqb[0] / nH;                                 // id 2
  double K_Hm = nH / Keqb[2];                                 // id 7
  double K_H2p = Keqb[0] / Keqb[4];                           // id 2 / id 9
  double K_H2 = nH / Keqb[2] / Keqb[3];                       // id 7, id 8
  double K_H3p = nH * Keqb[0] /                               // id 2
                 (Keqb[2] * Keqb[3] * Keqb[4] * Keqb[5]);     // id 7,8,9,26
  // He sequence (He, He+ only)
  double K_Hep = Keqb[1] / nH;                                // id 4
  // D sequence (D, HD, D+ only)
  double K_Dp = Keqb[6] / nH;                                 // id 51
  double K_HD = nH / Keqb[7];                                 // id 54
  // Li sequence (Li, Li+ only)
  double K_Lip = Keqb[22] / nH;                               // id 801

  // clang-format on
  // ── fill_y: write yy[0..30] from (y_e, y_H, y_C, y_O) ──────────────────────
  auto fill_y = [&](auto& yy, double ye, double yh, double yc, double yo) {
    double yHp = K_Hp * yh / ye;
    double yHm = K_Hm * yh * ye;
    double yH2p = K_H2p * yh * yh / ye;
    double yH2 = K_H2 * yh * yh;
    double yH3p = K_H3p * yh * yh * yh / ye;

    yy[Sp::H] = yh;
    yy[Sp::H2] = yH2;
    yy[Sp::e] = ye;
    yy[Sp::Hp] = yHp;
    yy[Sp::H2p] = yH2p;
    yy[Sp::H3p] = yH3p;
    yy[Sp::Hm] = yHm;

    // He sequence
    yy[Sp::He] = yHe / (1.0 + K_Hep / ye);
    yy[Sp::Hep] = K_Hep * yy[Sp::He] / ye;

    // D sequence
    yy[Sp::D] = yD / (1.0 + K_HD * yh + K_Dp / ye);
    yy[Sp::HD] = K_HD * yy[Sp::D] * yh;
    yy[Sp::Dp] = K_Dp * yy[Sp::D] / ye;

    // Li sequence
    yy[Sp::Li] = yLi / (1.0 + K_Lip / ye);
    yy[Sp::Lip] = K_Lip * yy[Sp::Li] / ye;

    // Mg sequence  (Mg+ from the alkali Saha balance, neutral Mg below)
    yy[Sp::Mgp] = XMg / (1.0 + nH * ye / Keqb[21]);  // id 718
    yy[Sp::Mg] = (Keqb[21] > 1.0e-300) ? nH * ye * yy[Sp::Mgp] / Keqb[21] : 0.0;

    // Carbon sequence  (C, CH, CH2, C+)
    yy[Sp::C] = yc;
    yy[Sp::CH] = yh * yc * nH / Keqb[9];     // id 185
    yy[Sp::CH2] = yH2 * yc * nH / Keqb[19];  // id 538
    yy[Sp::Cp] = (yc / ye) * Keqb[13] / nH;  // id 507

    // Oxygen sequence  (O, O2, OH, CO, H2O, HCO and O/C ions)
    yy[Sp::O] = yo;
    yy[Sp::O2] = yo * yo * nH / Keqb[12];             // id 194
    yy[Sp::OH] = yh * yo * nH / Keqb[11];             // id 193
    yy[Sp::CO] = yc * yo * nH / Keqb[10];             // id 187
    yy[Sp::H2O] = yh * yy[Sp::OH] * nH / Keqb[20];    // id 642
    yy[Sp::HCO] = Keqb[8] * yH2 * yy[Sp::CO] / yh;    // id 184
    yy[Sp::Op] = (yo / ye) * Keqb[14] / nH;           // id 515
    yy[Sp::OHp] = Keqb[15] * yh * yo / ye;            // id 518
    yy[Sp::H2Op] = Keqb[16] * yh * yy[Sp::OH] / ye;   // id 521
    yy[Sp::H3Op] = Keqb[17] * yh * yy[Sp::H2O] / ye;  // id 523
    yy[Sp::HCOp] = Keqb[18] * yh * yy[Sp::CO] / ye;   // id 527
  };

  // ── Conservation residuals over the compact species only ───────────────────
  auto F_cha_fn = [&](const auto& yy) -> double {
    return q_grain + yy[Sp::Hp] + yy[Sp::H2p] + yy[Sp::H3p] + yy[Sp::Hep] +
           yy[Sp::Dp] + yy[Sp::Cp] + yy[Sp::Op] + yy[Sp::OHp] + yy[Sp::H2Op] +
           yy[Sp::HCOp] + yy[Sp::H3Op] + yy[Sp::Lip] + yy[Sp::Mgp] -
           (yy[Sp::e] + yy[Sp::Hm]);
  };
  auto H_content = [&](const auto& yy) -> double {
    return yy[Sp::H] + yy[Sp::Hp] + yy[Sp::Hm] + yy[Sp::HD] + yy[Sp::CH] +
           yy[Sp::OH] + yy[Sp::HCO] + yy[Sp::OHp] + yy[Sp::HCOp] +
           2.0 * (yy[Sp::H2] + yy[Sp::H2p] + yy[Sp::CH2] + yy[Sp::H2O] +
                  yy[Sp::H2Op]) +
           3.0 * (yy[Sp::H3p] + yy[Sp::H3Op]);
  };
  auto C_content = [&](const auto& yy) -> double {
    return yy[Sp::C] + yy[Sp::CH] + yy[Sp::CH2] + yy[Sp::Cp] + yy[Sp::CO] +
           yy[Sp::HCO] + yy[Sp::HCOp];
  };
  auto O_content = [&](const auto& yy) -> double {
    return yy[Sp::O] + yy[Sp::OH] + yy[Sp::CO] + yy[Sp::H2O] + yy[Sp::HCO] +
           yy[Sp::Op] + yy[Sp::OHp] + yy[Sp::H2Op] + yy[Sp::HCOp] +
           yy[Sp::H3Op] + 2.0 * yy[Sp::O2];
  };
  const double yH_tot = H_content(y);
  const double XC = C_content(y);
  const double XO = O_content(y);
  auto F_hyd_fn = [&](const auto& yy) -> double {
    return H_content(yy) - yH_tot;
  };
  auto F_car_fn = [&](const auto& yy) -> double { return C_content(yy) - XC; };
  auto F_oxy_fn = [&](const auto& yy) -> double { return O_content(yy) - XO; };

  // ── Initial guess from current y ──────────────────────────────────────────
  double y_e = y[Sp::e];
  double y_H = y[Sp::H];
  double y_C = y[Sp::C];
  double y_O = y[Sp::O];
  if (y_e <= 0.0) y_e = 1.0e-20;
  if (y_H <= 0.0) y_H = 1.0e-20;
  if (y_C <= 0.0) y_C = (XC > 0.0) ? 1.0e-20 : 0.0;
  if (y_O <= 0.0) y_O = (XO > 0.0) ? 1.0e-20 : 0.0;

  std::array<double, metal_grain_minimal::N_sp> ytmp{};
  std::array<double, 4> x = {y_e, y_H, y_C, y_O};
  const std::array<bool, 4> active = {true, true, (y_C > 0.0 && XC > 0.0),
                                      (y_O > 0.0 && XO > 0.0)};

  metal_saha_solve(
      x, active, K_Hp, K_H2, K_Dp * yD + K_Lip * yLi, yH_tot, XC, XO,
      [&](const std::array<double, 4>& xx, std::array<double, 4>& fv) {
        fill_y(ytmp, xx[0], xx[1], xx[2], xx[3]);
        fv[0] = F_cha_fn(ytmp);
        fv[1] = F_hyd_fn(ytmp);
        fv[2] = F_car_fn(ytmp);
        fv[3] = F_oxy_fn(ytmp);
      });

  // ── Final fill at the converged root — writes the real output ──────────────
  fill_y(y, x[0], x[1], x[2], x[3]);
}

}  // namespace detail

// ─────────────────────────────────────────────────────────────────────────────
// Saha-equilibrium policy for the metal_grain model.  The integrator selects
// this via Model::SahaPolicy (forward-declared in models/metal_grain/traits.h)
// and calls solve() in the high-density regime.  The element budgets come from y
// itself, so the 4D equichem_metal needs nothing from params.
// ─────────────────────────────────────────────────────────────────────────────
struct MetalSaha {
  template <int N_sp, int N_react>
  static void solve(double nH, double T_K, const ChemParams& /*params*/,
                    std::array<double, N_sp>& y,
                    const ReactionTable<N_sp, N_react>& tbl) {
    detail::equichem_metal(nH, T_K, y, tbl);
  }
};

// ─────────────────────────────────────────────────────────────────────────────
// Saha-equilibrium policy for the compact metal_grain Minimal model.
//
// The compact 4D Newton (detail::equichem_metal_minimal) needs the partition
// functions of the 44 compact species.  The metal PF evaluator is written in the
// full 95-species layout; it evaluates only the compact species
// (kKeptFullSpecies) on the full PF-loaded table the runtime owns alongside the
// compact table, and they are copied onto the compact species by canonical
// identity.  That table is
// supplied when the compact metal runtime is constructed (table/PF phase);
// until then this policy is guarded so it cannot be silently instantiated with
// no PF source.
// ─────────────────────────────────────────────────────────────────────────────
struct MinimalMetalSaha {
  template <int N_sp, int N_react>
  static void solve(double nH, double T_K, const ChemParams& /*params*/,
                    std::array<double, N_sp>& y,
                    const ReactionTable<N_sp, N_react>& tbl) {
    const auto& full = *tbl.aux_full_metal;
    std::array<double, metal_grain::N_sp + 3> pf_full{};
    eval_metal_partition_functions<metal_grain::N_sp, metal_grain::N_react>(
        T_K, full, pf_full,
        [](int i) { return metal_grain_minimal::kKeptFullSpecies[i]; });

    std::array<double, metal_grain_minimal::N_sp + 3> pf{};
    pf.fill(1.0);
    for (int i = 0; i < metal_grain_minimal::N_sp; ++i)
      pf[i] = pf_full[metal_grain::Species::local(
          metal_grain_minimal::Species::canonical(i))];

    detail::equichem_metal_minimal(nH, T_K, y, tbl, pf);
  }
};

}  // namespace arche
