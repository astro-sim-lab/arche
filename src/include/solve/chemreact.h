// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// chemreact.h — model terms of the chemistry step shared by the integrator
//
// compute_chemistry_cooling (the chemistry heat term of the energy equation) and
// saha_nonthermal_shift (the test of whether the Saha regime may be used with
// non-thermal ionisation). The step itself is solved by solve/integrator.h.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <utility>

#include "core/newton.h"  // gaussj_solve, newton_solve, NewtonOpts, LinSolve
#include "cooling/cooling.h"              // c_H2
#include "kinetics/partition_function.h"  // PfProvider, pf_eval, detail primitives
#include "models/rate_kernel.h"           // compute_base_rates, compute_rates
#include "kinetics/topology.h"
#include "kinetics/reaction_index.h"
#include "thermo/eos.h"

namespace arche {

// ─────────────────────────────────────────────────────────────────────────────
// compute_chemistry_cooling<Model> — chemistry heating/cooling rate [erg g^-1
// s^-1]
//   Handles both zero_metal (primordial) and metal_grain networks via
//   if constexpr on Model traits.
//
//   Plain `inline` (not force-inlined): inlining is left to the compiler.
// ─────────────────────────────────────────────────────────────────────────────
// The factors of the chemistry heat that depend on T alone (the temperature
// dependence of the H2 critical density, used below nH = 1e13 cm^-3), so that
// evaluations at one T can share them.
struct ChemCoolingTerms {
  double exp_H = 0.0, exp_H2 = 0.0, sqrt_T = 0.0;
};

inline ChemCoolingTerms chemistry_cooling_terms(double nH, double T_K) {
  ChemCoolingTerms t;
  if (nH < 1.0e13) {
    t.exp_H = std::exp(-std::pow(400.0 / T_K, 2.0));
    t.exp_H2 = std::exp(-12000.0 / (T_K + 1200.0));
    t.sqrt_T = std::sqrt(T_K);
  }
  return t;
}

template <class Model>
inline void compute_chemistry_cooling(
    const std::array<double, Model::N_sp>& y,
    const std::array<double, Model::N_sp>& dy,
    const std::array<double, 2 * Model::N_react>& k_rxn, double nH, double dt,
    const ChemCoolingTerms& tt, double& Lambda_chem) {
  using Sp = typename Model::Sp;  // model-specific species names: y[Sp::H] etc.
  constexpr int N_react = Model::N_react;
  constexpr double yHe = abundance_ref::yHe;
  constexpr double m_p = phys::m_p;

  double L_H2 = 0.0;
  double dyHp = dy[Sp::Hp];
  // He+/He++ ionization-energy terms are a composable contribution: models that
  // carry He+ provide dyHep (has_he_ion), models that carry He++ provide dyHepp
  // (has_he_pp); a model without them leaves the term zero and never names
  // Sp::Hep / Sp::Hepp.  With both flags true the assignments below are the full
  // models' additive energy sum, in order.
  double dyHep = 0.0;
  double dyHepp = 0.0;
  if constexpr (Model::cooling::has_he_ion) {
    dyHep = dy[Sp::Hep];
  }
  if constexpr (Model::cooling::has_he_pp) {
    dyHepp = dy[Sp::Hepp];
  }

  if (nH < 1.0e13) {
    double denom_cr = 1.6 * y[Sp::H] * tt.exp_H + 1.4 * y[Sp::H2] * tt.exp_H2;
    double n_cr = (denom_cr > 1.0e-30) ? 1.0e6 / tt.sqrt_T / denom_cr : 1.0e50;
    double crit = 1.0 / (1.0 + n_cr / nH);

    // k_rxn[] indices are addressed through the model's cooling rate-index map
    // so a compact model (different reaction numbering) reaches the same
    // reactions; for the full models these equal the shared rxn:: ids.  The
    // H2–He dissociation channel is a composable term (a model without it sets
    // has_h2_he_cool = false); each accumulation runs left to right in the order
    // written.
    using cidx = typename Model::cooling::rate_idx;
    double rtHm = k_rxn[cidx::Hm_H_to_H2_e] * y[Sp::Hm] * y[Sp::H];
    double rtH2p = k_rxn[cidx::H2p_H_to_H2_Hp] * y[Sp::H2p] * y[Sp::H];
    double rt3b_sum =
        k_rxn[cidx::three_H] * y[Sp::H] * y[Sp::H] * y[Sp::H] +
        k_rxn[cidx::H2H2_dis + N_react] * y[Sp::H] * y[Sp::H] * y[Sp::H2];
    if constexpr (Model::cooling::has_h2_he_cool)
      rt3b_sum +=
          k_rxn[cidx::H2_He_dis + N_react] * y[Sp::H] * y[Sp::H] * y[Sp::He];
    double rt3b = rt3b_sum * nH;
    double rtdis = k_rxn[cidx::three_H + N_react] * y[Sp::H2] * y[Sp::H] +
                   k_rxn[cidx::H2H2_dis] * y[Sp::H2] * y[Sp::H2];
    if constexpr (Model::cooling::has_h2_he_cool)
      rtdis += k_rxn[cidx::H2_He_dis] * y[Sp::H2] * y[Sp::He];

    if constexpr (Model::has_grain) {
      // Grain H2 formation cooling.  The surface channels (physisorbed /
      // chemisorbed H) need H_p / H_c, which the compact network drops, so they
      // are composable (gh2::has_surface); with them the full expression is
      // evaluated.
      using gh2 = typename Model::cooling::grain_H2;
      double rtgrain;
      if constexpr (gh2::has_surface)
        rtgrain = k_rxn[gh2::two_H] * y[Sp::H] +
                  k_rxn[gh2::surface_grgr] * y[Sp::H_p] * y[Sp::H_p] +
                  k_rxn[gh2::surface_grdust] * y[Sp::H_p] * y[Sp::H_c] +
                  k_rxn[gh2::H_dust_3body] * y[Sp::H] * y[Sp::H_c] * nH +
                  k_rxn[gh2::H_gr_3body] * y[Sp::H] * y[Sp::H_p] * nH;
      else
        rtgrain = k_rxn[gh2::two_H] * y[Sp::H];

      L_H2 = -(rtgrain * (0.2 + 4.2 * crit) + rtHm * 3.53 * crit +
               rtH2p * 1.83 * crit + (rt3b * crit - rtdis) * 4.48) *
             nH * 1.60219e-12;
    } else {
      L_H2 = -(rtHm * 3.53 * crit + rtH2p * 1.83 * crit +
               (rt3b * crit - rtdis) * 4.48) *
             nH * 1.60219e-12;
    }

    // H+ recombination minus the CR-producer losses; the CR term is composable
    // (a CR-free model omits it).  Each branch is the model's full expression.
    if constexpr (Model::cooling::has_cr_loss) {
      using Hp = typename Model::cooling::Hp_producers;
      // The second H2 CR channel (CR-induced photon) is composable: a model
      // whose keep-set lacks it sets H2_CR_ch2 = -1.  Full models keep both
      // channels, so the sum (ch1 + ch2) stays byte-identical.
      double h2_cr_loss = k_rxn[Hp::H2_CR_ch1];
      if constexpr (Hp::H2_CR_ch2 >= 0) h2_cr_loss += k_rxn[Hp::H2_CR_ch2];
      dyHp =
          dy[Sp::Hp] + (k_rxn[cidx::Hp_rec_caseB] * y[Sp::Hp] * y[Sp::e] * nH -
                        (k_rxn[Hp::H_CR] + k_rxn[Hp::H_CRph]) * y[Sp::H] -
                        h2_cr_loss * y[Sp::H2]) *
                           dt;
    } else {
      dyHp = dy[Sp::Hp] +
             (k_rxn[cidx::Hp_rec_caseB] * y[Sp::Hp] * y[Sp::e] * nH) * dt;
    }
    // He+ recombination minus CR losses — only for models carrying He ions; the
    // CR-producer loss is composable (a He+-but-CR-free model omits it, mirroring
    // the dyHp split above).  Full models keep both terms, byte-identical.
    if constexpr (Model::cooling::has_he_ion) {
      if constexpr (Model::cooling::has_cr_loss) {
        using Hep = typename Model::cooling::Hep_producers;
        dyHep = dy[Sp::Hep] +
                (k_rxn[cidx::Hep_rec] * y[Sp::Hep] * y[Sp::e] * nH -
                 (k_rxn[Hep::He_CR] + k_rxn[Hep::He_CRph]) * y[Sp::He]) *
                    dt;
      } else {
        dyHep = dy[Sp::Hep] +
                (k_rxn[cidx::Hep_rec] * y[Sp::Hep] * y[Sp::e] * nH) * dt;
      }
    }

  } else {
    double dyH2 = dy[Sp::H2];
    L_H2 = -7.18e-12 * dyH2 / dt;
    constexpr double L_H2_max = 1.0e20;
    L_H2 = std::max(-L_H2_max, std::min(L_H2, L_H2_max));
  }

  Lambda_chem =
      ((2.18e-11 * dyHp + 3.94e-11 * dyHep + 12.66e-11 * dyHepp) / dt + L_H2) /
      ((1.0 + 4.0 * yHe) * m_p);
  // Preserve non-finite output so the acceptance layer rejects it; never
  // silently turn invalid chemical heating into a zero contribution.
}

template <class Model>
inline void compute_chemistry_cooling(
    const std::array<double, Model::N_sp>& y,
    const std::array<double, Model::N_sp>& dy,
    const std::array<double, 2 * Model::N_react>& k_rxn, double nH, double T_K,
    double dt, double& Lambda_chem) {
  compute_chemistry_cooling<Model>(
      y, dy, k_rxn, nH, dt, chemistry_cooling_terms(nH, T_K), Lambda_chem);
}

// ─────────────────────────────────────────────────────────────────────────────
// saha_nonthermal_shift — how far non-thermal ionization moves y_e off Saha
//
// At the Saha state y_S the network is evaluated twice: with the caller's zeta /
// zeta_X, and with both set to zero. Their difference dr is the non-thermal
// production alone (the difference, not r itself, because the network's thermal rates
// and the Saha partition functions need not balance exactly). Treated as a perturbation,
// it shifts the electron abundance by dy_e = dr_e / |J_ee| in steady state, where
// |J_ee| y_e is the electrons' gross destruction (recombination) at y_S. The returned
// value is |dy_e| / y_e. For two-body electron-ion recombination the exact linear shift
// is half of this, so the estimate errs on the side of keeping the Newton solve.
template <class Model>
double saha_nonthermal_shift(
    double nH, double T_K, const std::array<double, Model::N_sp>& y_S,
    const ReactionTable<Model::N_sp, Model::N_react>& tbl,
    const ChemParams& params,
    const std::array<double, 2 * Model::N_react>& var) {
  constexpr int N_sp = Model::N_sp;
  constexpr int N_react = Model::N_react;
  using Sp = typename Model::Sp;
  ChemParams p = params;
  p.H = y_S[Sp::H];
  p.H2 = 2.0 * y_S[Sp::H2];
  p.He = y_S[Sp::He];
  if constexpr (Model::has_uv_shield) {
    p.J_H2 = y_S[Sp::H2_p];
    p.J_H2O = y_S[Sp::H2O_p];
    p.J_tot = y_S[Sp::H_c] + y_S[Sp::D_c];
  }
  ChemParams p0 = p;
  p0.zeta = 0.0;
  p0.zeta_X = 0.0;
  const auto eos = thermo::make_eos_composition<Model>(y_S);
  const double mu = thermo::mean_molecular_weight(eos);
  std::array<double, 2 * N_react> k{}, k0{};
  compute_base_rates<Model>(nH, T_K, mu, p, tbl, k);
  compute_base_rates<Model>(nH, T_K, mu, p0, tbl, k0);
  std::array<double, N_sp> r{}, r0{};
  std::array<double, N_sp * N_sp> J{};
  auto v = var, v0 = var;  // compute_rates writes its cache; keep the caller's
  compute_rates<Model>(k, nH, y_S, tbl, r, J, v);
  compute_rhs<Model>(k0, nH, y_S, tbl, r0, v0);
  const double ye = y_S[Sp::e];
  const double loss = std::abs(J[Sp::e * N_sp + Sp::e]) * ye;
  const double dr = std::abs(r[Sp::e] - r0[Sp::e]);
  if (!(ye > 0.0) || !(loss > 0.0) || !std::isfinite(dr)) {
    return std::numeric_limits<
        double>::infinity();  // cannot judge: keep Newton
  }
  return dr / loss;
}

}  // namespace arche
