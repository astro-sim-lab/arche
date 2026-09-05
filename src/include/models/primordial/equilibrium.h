// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// equilibrium.h — equichem / equichem_minimal helpers for the primordial Saha
// branch
//
// equichem<N_sp, N_react> and equichem_minimal<N_sp, N_react> compute Saha
// equilibrium abundances in the high-density branch (nH > 1e18).  chemreact()
// dispatches to these through forward declarations in solve/chemreact.h; the
// including translation unit (chemistry.h) pulls in these definitions so the
// template instantiations resolve.
// ---------------------------------------------------------------------------
#include <array>
#include <cmath>

#include "kinetics/reaction_index.h"
#include "core/newton.h"  // newton_solve, CramerSolve2, NewtonOpts
#include "kinetics/rates.h"
#include "models/primordial/minimal.h"
#include "models/primordial/partition_function_prim.h"

namespace arche {

// ─────────────────────────────────────────────────────────────────────────────
// Shared primordial-Saha building blocks
//
// build_saha_keqb and fill_H_saha are the parts that every primordial Saha
// solver (full and compact) carries verbatim.  Factoring them keeps the
// model-specific element blocks and conservation residuals as the only
// difference between solvers, without changing the floating-point operations
// (so the high-density branch stays bit-for-bit identical).
// ─────────────────────────────────────────────────────────────────────────────

// Per-reaction Saha equilibrium constants Keqb[num-1] = exp(lnKeqb), built from
// the model's Saha sub-table and partition functions.  Index-agnostic: it reads
// the table's stored species indices, so the same builder serves every model.
template <int N_sp, int N_react>
void build_saha_keqb(double T_K, const ReactionTable<N_sp, N_react>& tbl,
                     const std::array<double, N_sp + 3>& pf,
                     std::array<double, N_react>& Keqb) {
  constexpr double k_B = phys::k_B;
  constexpr double h_P = phys::h_P;
  constexpr double pi = phys::pi;
  Keqb.fill(0.0);
  double lnT32 = 1.5 * std::log(2.0 * pi * k_B * T_K / (h_P * h_P));
  for (int i = 0; i < tbl.n_saha; ++i) {
    int num = tbl.saha[i].num;
    int r1 = tbl.saha[i].reactants[0];
    int r2 = tbl.saha[i].reactants[1];
    int p1 = tbl.saha[i].products[0];
    int p2 = tbl.saha[i].products[1];  // IDX_PHOTON = N_sp + 1
    int nr = tbl.saha[i].n_reactants;
    int np = tbl.saha[i].n_products;
    double Cm = tbl.saha[i].Cmass;
    double dE = tbl.saha[i].delE;
    double lnC1 = static_cast<double>(nr - np) * lnT32;
    double lnCpf = std::log(pf[r1]) + std::log(pf[r2]) - std::log(pf[p1]);
    if (np == 2 && p2 != ReactionTable<N_sp, N_react>::IDX_PHOTON)
      lnCpf -= std::log(pf[p2]);
    double lnKeqb = lnC1 + std::log(Cm) + lnCpf - dE / (k_B * T_K);
    if (num >= 1 && num <= N_react) Keqb[num - 1] = std::exp(lnKeqb);
  }
}

// Fill the seven shared H-network species (H, H2, e-, H+, H2+, H3+, H-) from the
// Saha trial (ye, yh) and the H equilibrium ratios; returns y_H+ (used by the
// He/D/Li blocks).  Templated on the model's Sp enum so the same fill serves the
// full and compact primordial index spaces (both lay H..H- out identically).
template <class Sp, std::size_t N>
double fill_H_saha(double ye, double yh, double K_Hp, double K_Hm, double K_H2p,
                   double K_H2, double K_H3p, std::array<double, N>& yy) {
  double y_Hp = K_Hp * yh / ye;
  double y_Hm = K_Hm * yh * ye;
  double y_H2p = K_H2p * yh * yh / ye;
  double y_H2 = K_H2 * yh * yh;
  double y_H3p = K_H3p * yh * yh * yh / ye;

  yy[Sp::H] = yh;       // H
  yy[Sp::H2] = y_H2;    // H2
  yy[Sp::e] = ye;       // e-
  yy[Sp::Hp] = y_Hp;    // H+
  yy[Sp::H2p] = y_H2p;  // H2+
  yy[Sp::H3p] = y_H3p;  // H3+
  yy[Sp::Hm] = y_Hm;    // H-
  return y_Hp;
}

// ─────────────────────────────────────────────────────────────────────────────
// equichem — Saha equilibrium abundances (high-density branch nH > 1e18)
//
//   Computes equilibrium constants K_eq for each Saha reaction from the
//   primordial partition-function provider (eval_pf_set), then iterates 2D
//   Newton-Raphson on charge neutrality (F_cha) and hydrogen conservation
//   (F_hyd).
// ─────────────────────────────────────────────────────────────────────────────
template <int N_sp, int N_react>
void equichem(double nH, double T_K, std::array<double, N_sp>& y,
              const ReactionTable<N_sp, N_react>& tbl) {
  using Sp = zero_metal::Sp;  // primordial Saha (y[Sp::H] etc.)
  constexpr double yHe = abundance_ref::yHe;
  constexpr double yD = abundance_ref::yD;
  constexpr double yLi = abundance_ref::yLi;
  constexpr double eps_it = numerics::eps_it_prim;

  // Compute partition functions (0-based species + 3 sentinel slots) through
  // the shared primordial provider, addressed by canonical identity.
  std::array<double, N_sp + 3> pf{};
  pf_prim::eval_pf_set<zero_metal::Species, N_react>(T_K, tbl, pf);

  // Saha equilibrium constants, indexed by reaction number (Keqb[num-1]).
  std::array<double, N_react> Keqb{};
  build_saha_keqb<N_sp, N_react>(T_K, tbl, pf, Keqb);

  // Convenience: 1-based Keqb(k) → 0-based Keqb[k-1]
  auto K = [&](int k) -> double { return Keqb[k - 1]; };

  // Build abundance ratios (independent of y_e, y_H except where noted)
  // H sequence
  double K_Hp = K(2) / nH;
  double K_Hm = nH / K(7);
  double K_H2p = K(2) / K(9);
  double K_H2 = nH / K(7) / K(8);
  double K_H3p = nH * K(2) / K(7) / K(8) / K(9) / K(26);
  // He sequence
  double K_Hep = K(4) / nH;
  double K_He2p = K(4) * K(6) / nH / nH;
  double K_HeHp = nH / K(44);
  // D sequence
  double K_Dp = K(51) / nH;
  double K_HD = nH / K(54);
  double K_HDp = nH / K(60);
  double K_Dm = nH / K(63);
  // Li sequence
  double K_Lip = K(101) / nH;
  double K_Li2p = K_Lip * K(121) / nH;
  double K_Li3p = K_Li2p * K(122) / nH;
  double K_LiH = nH / K(118);
  double K_LiHp = nH / K(113);
  double K_Lim = nH / K(104);

  // Initial guess from current y
  double y_e = y[Sp::e];  // e-   (0-based: index 2)
  double y_H = y[Sp::H];  // H    (0-based: index 0)

  // Helper: fill all species from y_e, y_H
  auto fill_y = [&](double ye, double yh, std::array<double, N_sp>& yy) {
    double y_Hp = fill_H_saha<Sp>(ye, yh, K_Hp, K_Hm, K_H2p, K_H2, K_H3p, yy);

    // He sequence
    {
      double den = 1.0 + K_Hep / ye + K_He2p / ye / ye + K_HeHp * y_Hp;
      yy[Sp::He] = yHe / den;                        // He
      yy[Sp::Hep] = K_Hep * yy[Sp::He] / ye;         // He+
      yy[Sp::Hepp] = K_He2p * yy[Sp::He] / ye / ye;  // He++
      yy[Sp::HeHp] = K_HeHp * y_Hp * yy[Sp::He];     // HeH+
    }
    // D sequence
    {
      double den = 1.0 + K_HD * yh + K_Dp / ye + K_HDp * y_Hp + K_Dm * ye;
      yy[Sp::D] = yD / den;                    // D
      yy[Sp::HD] = K_HD * yy[Sp::D] * yh;      // HD
      yy[Sp::Dp] = K_Dp * yy[Sp::D] / ye;      // D+
      yy[Sp::HDp] = K_HDp * yy[Sp::D] * y_Hp;  // HD+
      yy[Sp::Dm] = K_Dm * yy[Sp::D] * ye;      // D-
    }
    // Li sequence (only if N_sp >= 23)
    if constexpr (N_sp >= 23) {
      double den = 1.0 + K_LiH * yh + K_Lim * ye + K_LiHp * y_Hp + K_Lip / ye +
                   K_Li2p / ye / ye + K_Li3p / (ye * ye * ye);
      yy[Sp::Li] = yLi / den;                                // Li
      yy[Sp::LiH] = K_LiH * yy[Sp::Li] * yh;                 // LiH
      yy[Sp::Lip] = K_Lip * yy[Sp::Li] / ye;                 // Li+
      yy[Sp::Lim] = K_Lim * yy[Sp::Li] * ye;                 // Li-
      yy[Sp::LiHp] = K_LiHp * yy[Sp::Li] * y_Hp;             // LiH+
      yy[Sp::Lipp] = K_Li2p * yy[Sp::Li] / ye / ye;          // Li++
      yy[Sp::Lippp] = K_Li3p * yy[Sp::Li] / (ye * ye * ye);  // Li3+
    }
  };

  // Charge conservation:  Σ(z_i * y_i) = 0
  auto F_cha = [&](const std::array<double, N_sp>& yy) -> double {
    double val = yy[Sp::Hp] + yy[Sp::H2p] + yy[Sp::H3p] + yy[Sp::Hep] +
                 yy[Sp::HeHp] + yy[Sp::Dp] + yy[Sp::HDp] +
                 2.0 * (yy[Sp::Hepp] + yy[Sp::Lipp]) + 3.0 * yy[Sp::Lippp] -
                 (yy[Sp::e] + yy[Sp::Hm] + yy[Sp::Dm] + yy[Sp::Lim]);
    if constexpr (N_sp >= 23) val += yy[Sp::Lip] + yy[Sp::LiHp];  // Li+, LiH+
    return val;
  };

  // Hydrogen conservation:  Σ(H-content_i * y_i) = 1
  auto F_hyd = [&](const std::array<double, N_sp>& yy) -> double {
    double val = yy[Sp::H] + yy[Sp::Hp] + yy[Sp::Hm] + yy[Sp::HeHp] +
                 yy[Sp::HD] + yy[Sp::HDp] + 2.0 * (yy[Sp::H2] + yy[Sp::H2p]) +
                 3.0 * yy[Sp::H3p] - 1.0;
    if constexpr (N_sp >= 23) val += yy[Sp::LiH] + yy[Sp::LiHp];  // LiH, LiH+
    return val;
  };

  std::array<double, N_sp> ytmp{};
  std::array<double, 2> x = {y_e, y_H};
  const std::array<bool, 2> active = {true, true};

  // 2D Newton–Raphson on (y_e, y_H): charge neutrality and H conservation.
  newton_solve<2, CramerSolve2>(
      x, active,
      [&](const std::array<double, 2>& xx, std::array<double, 2>& fv) {
        fill_y(xx[0], xx[1], ytmp);
        fv[0] = F_cha(ytmp);
        fv[1] = F_hyd(ytmp);
      },
      [](std::array<double, 2>& xx, const std::array<double, 2>& d) {
        xx[0] += d[0];
        xx[1] += d[1];
      },
      NewtonOpts{eps_it, 100, 1.0e-10});

  // Final fill
  fill_y(x[0], x[1], y);
}

// ─────────────────────────────────────────────────────────────────────────────
// Saha-equilibrium policies for the primordial models.  chemreact() selects one
// per model via Model::SahaPolicy (forward-declared in models/primordial/
// traits.h) and calls solve() in the high-density branch.
// ─────────────────────────────────────────────────────────────────────────────
struct FullPrimSaha {
  template <int N_sp, int N_react>
  static void solve(double nH, double T_K, const ChemParams& /*params*/,
                    std::array<double, N_sp>& y,
                    const ReactionTable<N_sp, N_react>& tbl) {
    equichem<N_sp, N_react>(nH, T_K, y, tbl);
  }
};

// ─────────────────────────────────────────────────────────────────────────────
// equichem_minimal — Saha equilibrium for the compact 15-species
// Nakauchi2019_Minimal.  Same retained-ion physics as equichem, but
// addressed in the compact local index space (zero_metal_minimal::Sp) with the
// 9-row compact Saha sub-table (ids 1..9 in keep-set order), and partition
// functions evaluated by canonical identity (pf_for_prim).
// ─────────────────────────────────────────────────────────────────────────────
template <int N_sp, int N_react>
void equichem_minimal(double nH, double T_K, std::array<double, N_sp>& y,
                      const ReactionTable<N_sp, N_react>& tbl) {
  using Sp = zero_metal_minimal::Sp;
  constexpr double yHe = abundance_ref::yHe;
  constexpr double yD = abundance_ref::yD;
  constexpr double yLi = abundance_ref::yLi;
  constexpr double eps_it = numerics::eps_it_prim;

  // Partition functions by canonical identity (compact local order).
  std::array<double, N_sp + 3> pf{};
  pf_prim::eval_pf_set<zero_metal_minimal::Species, N_react>(T_K, tbl, pf);

  // Saha equilibrium constants per entry (compact id 1..9 → Keqb[id-1]).
  std::array<double, N_react> Keqb{};
  build_saha_keqb<N_sp, N_react>(T_K, tbl, pf, Keqb);
  auto K = [&](int k) -> double { return Keqb[k - 1]; };

  // Compact Saha id map (keep-set order): 1:num2 2:num7 3:num8 4:num9 5:num26
  //                                       6:num51 7:num54 8:num101 9:num118
  double K_Hp = K(1) / nH;
  double K_Hm = nH / K(2);
  double K_H2p = K(1) / K(4);
  double K_H2 = nH / K(2) / K(3);
  double K_H3p = nH * K(1) / K(2) / K(3) / K(4) / K(5);
  double K_Dp = K(6) / nH;
  double K_HD = nH / K(7);
  double K_Lip = K(8) / nH;
  double K_LiH = nH / K(9);

  double y_e = y[Sp::e];
  double y_H = y[Sp::H];

  auto fill_y = [&](double ye, double yh, std::array<double, N_sp>& yy) {
    fill_H_saha<Sp>(ye, yh, K_Hp, K_Hm, K_H2p, K_H2, K_H3p, yy);
    yy[Sp::He] = yHe;  // neutral background
    yy[Sp::Hep] =
        0.0;  // He+ negligible at nH>1e18 (kept out of charge balance)
    {
      double den = 1.0 + K_HD * yh + K_Dp / ye;
      yy[Sp::D] = yD / den;
      yy[Sp::HD] = K_HD * yy[Sp::D] * yh;
      yy[Sp::Dp] = K_Dp * yy[Sp::D] / ye;
    }
    {
      double den = 1.0 + K_LiH * yh + K_Lip / ye;
      yy[Sp::Li] = yLi / den;
      yy[Sp::LiH] = K_LiH * yy[Sp::Li] * yh;
      yy[Sp::Lip] = K_Lip * yy[Sp::Li] / ye;
    }
  };

  auto F_cha = [&](const std::array<double, N_sp>& yy) -> double {
    return yy[Sp::Hp] + yy[Sp::H2p] + yy[Sp::H3p] + yy[Sp::Dp] + yy[Sp::Lip] -
           (yy[Sp::e] + yy[Sp::Hm]);
  };
  auto F_hyd = [&](const std::array<double, N_sp>& yy) -> double {
    return yy[Sp::H] + yy[Sp::Hp] + yy[Sp::Hm] + yy[Sp::HD] +
           2.0 * (yy[Sp::H2] + yy[Sp::H2p]) + 3.0 * yy[Sp::H3p] + yy[Sp::LiH] -
           1.0;
  };

  std::array<double, N_sp> ytmp{};
  std::array<double, 2> x = {y_e, y_H};
  const std::array<bool, 2> active = {true, true};
  newton_solve<2, CramerSolve2>(
      x, active,
      [&](const std::array<double, 2>& xx, std::array<double, 2>& fv) {
        fill_y(xx[0], xx[1], ytmp);
        fv[0] = F_cha(ytmp);
        fv[1] = F_hyd(ytmp);
      },
      [](std::array<double, 2>& xx, const std::array<double, 2>& d) {
        xx[0] += d[0];
        xx[1] += d[1];
      },
      NewtonOpts{eps_it, 100, 1.0e-10});
  fill_y(x[0], x[1], y);
}

struct MinimalSaha {
  template <int N_sp, int N_react>
  static void solve(double nH, double T_K, const ChemParams& /*params*/,
                    std::array<double, N_sp>& y,
                    const ReactionTable<N_sp, N_react>& tbl) {
    equichem_minimal<N_sp, N_react>(nH, T_K, y, tbl);
  }
};

}  // namespace arche
