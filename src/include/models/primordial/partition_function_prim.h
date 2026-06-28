// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// partition_function_prim.h — canonical-species-keyed partition functions for
// the primordial family.
//
// prim_provider(SpId) returns the PfProvider descriptor for one canonical
// primordial species; pf_eval (kinetics/partition_function.h) turns it into a
// value.  Addressing by canonical identity rather than a fixed array slot lets
// any primordial-family model — the full 23-species network or a compact
// selection — obtain correct partition functions for exactly the species it
// carries.  Numeric-table species (H3+, HD, HD+, LiH+) read their tabulated
// data from the model's per-slot table (tbl.pf_table[local_index]); the other
// kinds ignore it.
//
// eval_pf_set<Set>() fills a model's pf[0..N-1] in its local index order and
// leaves the three trailing sentinel slots at 1.0, so vacant / photon / CR
// reactants act as multiplicative no-ops in detailed-balance terms.
// ---------------------------------------------------------------------------
#include <array>
#include <cmath>
#include <utility>
#include <vector>

#include "core/species_catalog.h"
#include "core/state.h"
#include "kinetics/partition_function.h"  // PfProvider, pf_eval, detail::H2_pf
#include "kinetics/topology.h"            // ReactionTable

namespace arche {
namespace pf_prim {

// BC16 32-point tables (interpolated over detail::kTBC).
inline constexpr std::array<double, 32> kQH2p = {
    5.00000e-1, 5.00000e-1, 5.00000e-1, 5.00000e-1, 5.00000e-1, 5.00000e-1,
    5.00024e-1, 5.00921e-1, 5.15632e-1, 5.64397e-1, 7.65738e-1, 1.33880e+0,
    1.91022e+0, 2.68536e+0, 3.40918e+0, 4.35256e+0, 5.05717e+0, 6.23172e+0,
    7.40726e+0, 1.21361e+1, 1.70214e+1, 2.50338e+1, 4.11228e+1, 6.13523e+1,
    1.16205e+2, 1.93054e+2, 2.94603e+2, 4.22340e+2, 5.76489e+2, 7.56297e+2,
    9.60417e+2, 1.18728e+3};
inline constexpr std::array<double, 32> kQHeHp = {
    1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0,
    1.00000e+0, 1.00019e+0, 1.00482e+0, 1.02408e+0, 1.12058e+0, 1.45094e+0,
    1.83823e+0, 2.44356e+0, 3.05868e+0, 3.88521e+0, 4.50784e+0, 5.54901e+0,
    6.59355e+0, 1.08010e+1, 1.50831e+1, 2.18432e+1, 3.48084e+1, 5.08408e+1,
    9.58110e+1, 1.63926e+2, 2.59055e+2, 3.80397e+2, 5.24517e+2, 6.87325e+2,
    8.65041e+2, 1.05451e+3};
inline constexpr std::array<double, 32> kQLi = {
    2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
    2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
    2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
    2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00013e+0,
    2.00473e+0, 2.02869e+0, 2.08800e+0, 2.19942e+0, 2.39054e+0, 2.70188e+0,
    3.17983e+0, 3.86752e+0};
inline constexpr std::array<double, 32> kQLiH = {
    1.00000e+0, 1.00000e+0, 1.00001e+0, 1.00007e+0, 1.00247e+0, 1.04234e+0,
    1.14353e+0, 1.36475e+0, 1.79698e+0, 2.25093e+0, 3.17620e+0, 5.04593e+0,
    6.92359e+0, 9.74658e+0, 1.25749e+1, 1.63530e+1, 1.91923e+1, 2.39405e+1,
    2.87248e+1, 4.88043e+1, 7.16412e+1, 1.13107e+2, 2.04843e+2, 3.28787e+2,
    6.95286e+2, 1.26125e+3, 2.06412e+3, 3.13695e+3, 4.53795e+3, 6.36625e+3,
    8.75776e+3, 1.18704e+4};

// H (P09) — ground state 2, plus excited levels above 1000 K.
inline double pf_H_analytic(double T_K) {
  constexpr double k_B = phys::k_B;
  double q = 2.0;
  if (T_K > 1.0e3) {
    for (int i = 2; i <= 5; ++i) {
      double di = static_cast<double>(i);
      q += (2.0 * di * di) * std::exp(-13.598 * phys::eV_to_erg *
                                      (1.0 - 1.0 / (di * di)) / k_B / T_K);
    }
  }
  return q;
}

// H2 (Popovas & Jorgensen 2016 equilibrium fit, capped above 1000 K by the
// BFM energy-level partition function).
inline double pf_H2_analytic(double T_K) {
  double T2 = T_K * T_K, T3 = T2 * T_K, T4 = T3 * T_K, T5 = T4 * T_K;
  double q;
  if (T_K <= 200.0) {
    q = 2.673071615415136e-1 - 3.495586051757211e-3 * T_K +
        1.227901619954258e-4 * T2 - 5.776440695273789e-7 * T3 +
        9.251224490175610e-10 * T4;
  } else if (T_K <= 1000.0) {
    q = 1.410033600294133e-1 + 6.085738724141971e-3 * T_K -
        4.096994909866605e-7 * T2 + 4.220221708082499e-10 * T3 -
        8.790594164685680e-14 * T4;
  } else {
    q = -9.661842638994980e-1 + 7.302127874247883e-3 * T_K -
        6.760893004505151e-7 * T2 + 3.128741080316710e-10 * T3 -
        1.645206030945910e-14 * T4 + 2.788597060472472e-19 * T5;
    q = std::min(q, detail::H2_pf(T_K));
  }
  return q;
}

// PfProvider descriptor for one canonical primordial species.  Numeric-table
// species carry the divisor (nuclear-spin degeneracy removal) and floor; the
// table data itself is supplied at evaluation time from the model's slot.
inline PfProvider prim_provider(SpId id) {
  switch (id) {
    case SpId::H:
      return pf_analytic(&pf_H_analytic);
    case SpId::H2:
      return pf_analytic(&pf_H2_analytic);
    case SpId::e:
      return pf_const(2.0);
    case SpId::Hp:
      return pf_const(1.0);
    case SpId::H2p:
      return pf_bc16(&kQH2p);
    case SpId::H3p:
      return pf_numeric(8.0, /*floor_one=*/true);
    case SpId::Hm:
      return pf_const(1.0);
    case SpId::He:
      return pf_const(1.0);
    case SpId::Hep:
      return pf_const(2.0);
    case SpId::Hepp:
      return pf_const(1.0);
    case SpId::HeHp:
      return pf_bc16(&kQHeHp);
    case SpId::D:
      return pf_const(2.0);
    case SpId::HD:
      return pf_numeric(6.0);
    case SpId::Dp:
      return pf_const(1.0);
    case SpId::HDp:
      return pf_numeric(6.0);
    case SpId::Dm:
      return pf_const(1.0);
    case SpId::Li:
      return pf_bc16(&kQLi);
    case SpId::LiH:
      return pf_bc16(&kQLiH);
    case SpId::Lip:
      return pf_const(1.0);
    case SpId::Lim:
      return pf_const(1.0);
    case SpId::LiHp:
      return pf_numeric(2.0);
    case SpId::Lipp:
      return pf_const(2.0);
    case SpId::Lippp:
      return pf_const(1.0);
    default:
      return pf_const(1.0);
  }
}

// Partition function for one canonical primordial species.
//   numeric_table — the per-species tabulated (T, pf) data for the model's
//   local slot (used only by H3+/HD/HD+/LiH+); empty otherwise.
inline double pf_for_prim(
    SpId id, double T_K,
    const std::vector<std::pair<double, double>>& numeric_table) {
  return pf_eval(prim_provider(id), T_K, numeric_table);
}

// Fill pf[0..Set::N-1] (+3 sentinel slots = 1.0) for a SpeciesSet, addressing
// each species by canonical identity.  `pf_tables[i]` is the model's per-slot
// numeric table (tbl.pf_table[i]).
template <class Set, int N_react>
void eval_pf_set(double T_K, const ReactionTable<Set::N, N_react>& tbl,
                 std::array<double, Set::N + 3>& pf) {
  pf.fill(1.0);
  for (int i = 0; i < Set::N; ++i)
    pf[i] = pf_eval(prim_provider(Set::canonical(i)), T_K, tbl.pf_table[i]);
}

}  // namespace pf_prim
}  // namespace arche
