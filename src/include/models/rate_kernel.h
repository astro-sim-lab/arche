// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
#include <algorithm>
#include <array>
#include <cmath>

#include "kinetics/reaction_index.h"
#include "core/state.h"
#include "kinetics/partition_function.h"
#include "kinetics/rates.h"
#include "kinetics/topology.h"
#include "model_traits.h"
#include "models/metal_grain/partition_function_metal.h"
#include "models/metal_grain/reaction_metal.h"
#include "models/primordial/minimal.h"  // zero_metal_minimal::loop, cr slots
#include "models/primordial/partition_function_prim.h"
#include "models/primordial/rate_laws.h"

// ---------------------------------------------------------------------------
// rate_kernel.h
//
// Model-parameterised reaction-rate kernel.  Assembles each model's forward
// and reverse rate coefficients from the shared rate primitives and the
// per-model forward laws, then turns k_rxn and the abundance vector into
// dy/dt and the Jacobian for the Newton solver.
//
// Provides:
//   compute_base_rates<Model>(...) — forward + reverse rate coefficients
//                                    k_rxn[0..2*N_react-1].
//   compute_rates<Model>(...)      — dy/dt (r_f) and Jacobian (dr_fdy) from
//                                    k_rxn and the abundance vector y.
//
// All rate coefficients in CGS (cm³/s or cm⁶/s for 3-body). k_rxn indexing:
// k_rxn[num-1] = forward, k_rxn[num-1+N_react] = reverse (1-based num).
// ---------------------------------------------------------------------------

namespace arche {

// ---------------------------------------------------------------------------
// compute_base_rates<Model>
//
// Computes all rate coefficients k_rxn[0..2*N_react-1].
//   k_rxn[num-1]            = forward rate for reaction num (1-based)
//   k_rxn[num-1 + N_react]  = reverse rate (computed from detailed balance)
// ---------------------------------------------------------------------------
template <class Model>
void compute_base_rates(double nH, double T_K, double mu,
                        const ChemParams& params,
                        const ReactionTable<Model::N_sp, Model::N_react>& tbl,
                        std::array<double, 2 * Model::N_react>& k_rxn) {
  constexpr int N_sp = Model::N_sp;
  constexpr int N_react = Model::N_react;

  k_rxn.fill(0.0);

  constexpr double k_B = phys::k_B;
  constexpr double h_P = phys::h_P;
  constexpr double pi = phys::pi;
  constexpr double G = phys::G;
  constexpr double m_p = phys::m_p;
  constexpr double yHe = abundance_ref::yHe;

  // ── Compact primordial model (Nakauchi2019_Minimal): a first-class composed
  //    model with its own 15-species / 24-reaction table.  Each kept reaction
  //    binds to the shared ratelaw::prim law; the generic detailed-balance
  //    reverse loop and the three retained LTE/low-density reverse blends
  //    complete k_rxn.  Returns before the full-network path, whose preamble
  //    indexes reaction slots that do not exist in the compact table. ──
  if constexpr (Model::is_compact_prim) {
    namespace rl = ratelaw::prim;
    PrimRateContext cc =
        build_prim_context(T_K, nH, params.H, params.H2, params.He, 0.0);
    double rho = nH * ((1.0 + 4.0 * yHe) * m_p);
    double lmbd_J = std::sqrt(pi * k_B * T_K / (G * mu * m_p * rho));
    double kap = detail::eval_opacity(T_K, rho);
    double tau_cnt = kap * rho * lmbd_J;

    std::array<double, N_sp + 3> pf;
    pf_prim::eval_pf_set<typename Model::Species, N_react>(T_K, tbl, pf);

    // Forward: compact slot -> canonical law (slot order = network keep-set).
    k_rxn[0] = rl::r1_Hp_rec_caseB(cc);       // num 2:   e- + H+  -> H + γ
    k_rxn[1] = rl::r6_H_e_to_Hm(cc);          // num 7:   H + e-   -> H- + γ
    k_rxn[2] = rl::r7_Hm_H_to_H2_e(cc);       // num 8:   H- + H   -> H2 + e-
    k_rxn[3] = rl::r8_H_Hp_to_H2p(cc);        // num 9:   H + H+   -> H2+ + γ
    k_rxn[4] = rl::r9_H2p_H_to_H2_Hp(cc);     // num 10:  H2+ + H  -> H2 + H+
    k_rxn[5] = rl::r18_three_H(cc);           // num 19:  3H       -> H2 + H
    k_rxn[6] = rl::r20_2H2_to_2H_H2(cc);      // num 21:  2H2      -> 2H + H2
    k_rxn[7] = rl::r25_H2p_H2_to_H3p_H(cc);   // num 26:  H2+ + H2 -> H3+ + H
    k_rxn[8] = rl::r33_H3p_e_to_H2_H(cc);     // num 34:  e- + H3+ -> H2 + H
    k_rxn[9] = rl::r34_H3p_e_to_3H(cc);       // num 35:  e- + H3+ -> 3H
    k_rxn[10] = rl::r51_D_Hp_to_Dp_H(cc);     // num 52:  D + H+   -> D+ + H
    k_rxn[11] = rl::r54_D_H2_to_H_HD(cc);     // num 55:  D + H2   -> H + HD
    k_rxn[12] = rl::r56_Dp_H2_to_Hp_HD(cc);   // num 57:  D+ + H2  -> H+ + HD
    k_rxn[13] = rl::r100_e_Lip_to_Li(cc);     // num 101: e- + Li+ -> Li + γ
    k_rxn[14] = rl::r101_Hm_Lip_to_H_Li(cc);  // num 102: H- + Li+ -> H + Li
    k_rxn[15] =
        rl::r105_Hp_Li_to_H_Lip_ph(cc);  // num 106: H+ + Li  -> H + Li+ + γ
    k_rxn[16] = rl::r110_H_LiH_to_H2_Li(cc);  // num 111: H + LiH  -> H2 + Li
    k_rxn[17] = rl::r117_H_Li_to_LiH(cc);     // num 118: H + Li   -> LiH + γ
    k_rxn[18] =
        rl::r129_H2_Li_to_H2_e_Lip(cc, tbl.reactions[18].delE);  // num 130
    // He+ / ion-processing reactions (appended; thermal He+ via id 4 reverse).
    k_rxn[19] = rl::r3_Hep_rec<RateForm::Primordial>(cc);  // num 4:  He+ + e- -> He + γ
    k_rxn[20] =
        rl::r14_Hm_Hp_to_2H<RateForm::Primordial>(cc);  // num 15: H+ + H-  -> 2H
    k_rxn[21] = rl::r16_H2p_e_to_2H(cc);          // num 17: e- + H2+ -> 2H
    k_rxn[22] = rl::r23_Hep_H2_to_Hp_H_He(cc);    // num 24: He+ + H2 -> H+ + H + He
    k_rxn[23] = rl::r37_Hm_H_to_2H_e(cc);         // num 38: H- + H   -> 2H + e-

    // Cosmic-ray channels (compact slots 24..32 = full ids 131..139), emitted
    // at the compact CR-block head; zero when params.zeta == 0.
    compute_CR_rates_prim<N_react>(k_rxn, params.zeta,
                                   zero_metal_minimal::loop::n_std);

    // Reverse rates via detailed balance for the standard reactions only (the
    // CR channels are first-order and carry no detailed-balance reverse).
    std::array<double, N_react> lnKeqb{};
    compute_reverse_loop<RateForm::Primordial, N_sp, N_react>(
        k_rxn, T_K, tau_cnt, zero_metal_minimal::loop::n_std, tbl, pf, lnKeqb);

    // Retained LTE/low-density reverse blends (mirror the full primordial
    // post-reverse corrections, at the compact slots of the blended reactions).
    {  // slot 2: H- + H -> H2 + e  (reverse blended with k40 low-density form)
      double k40LTE = std::exp(std::log(k_rxn[2]) + lnKeqb[2]);
      k_rxn[2 + N_react] = lte_blend_noguard(cc.k40v0, k40LTE, cc.crit_ratio);
    }
    {  // slot 5: 3H -> H2 + H  (reverse blended with k13 low-density form)
      double k13LTE = std::exp(std::log(k_rxn[5]) + lnKeqb[5]);
      k_rxn[5 + N_react] = lte_blend_noguard(cc.k13v0, k13LTE, cc.crit_ratio);
    }
    {  // slot 6: 2H2 -> 2H + H2  (reverse uses the LTE coefficient directly)
      k_rxn[6 + N_react] = std::exp(std::log(cc.k21LTE) + lnKeqb[6]);
    }
    return;
  }

  // ── Compact metal Minimal model (Nakauchi2021_Minimal): strategy-1 rate
  //    path.  The compact forward/reverse coefficients are gathered from a full
  //    metal_grain coefficient run (the metal forward laws write full slots, so
  //    they cannot be re-bound per compact slot the way the primordial compact
  //    path is — see the rate-path design note).  The gather reads the kept
  //    reactions out of the full k_rxn by num and writes them into the compact
  //    slots, using the full PF-loaded table the runtime owns alongside the
  //    compact table.  That table is supplied when the compact metal runtime is
  //    constructed (rate-gather phase); falling through to the full-network path
  //    below would index k_rxn slots far outside the compact stride, so this
  //    path is guarded until the gather is wired. ──
  if constexpr (Model::is_compact_metal) {
    // Run the full metal_grain coefficient kernel on the runtime-owned full
    // table, then gather the kept reactions into the compact slots.  The full
    // forward laws write full slots, so re-using the full path verbatim keeps
    // the kept coefficients bit-for-bit identical to the full model.
    constexpr int Nf = metal_grain::N_react;
    static thread_local std::array<double, 2 * Nf> k_full;
    compute_base_rates<Nakauchi2021>(nH, T_K, mu, params, *tbl.aux_full_metal,
                                     k_full);

    namespace net = metal_grain::net;
    int s = 0;
    // Standard block: forward + detailed-balance reverse.
    for (int full_num : net::kMetalMinimalKeep) {
      k_rxn[s] = k_full[full_num - 1];
      k_rxn[s + N_react] = k_full[full_num - 1 + Nf];
      ++s;
    }
    // Cosmic-ray block: first-order, forward only.
    for (int full_num : net::kMetalMinimalCRKeep) k_rxn[s++] = k_full[full_num - 1];
    // Ion-grain charge transfer: first-order, forward only.
    for (int full_num : net::kMetalMinimalChargeKeep)
      k_rxn[s++] = k_full[full_num - 1];
    // Grain-surface freeze-out: forward rates land in the grain band
    // (k_rxn[N_gas..]); the full slot is grain_surface_begin + (kgr num - 1).
    int g = metal_grain_minimal::N_gas;
    for (int grain_num : net::kMetalMinimalGrainKeep)
      k_rxn[g++] =
          k_full[metal_grain::slot::grain_surface_begin + grain_num - 1];
    // Grain-catalysed H2 / HD formation.  The full network builds these from
    // physisorbed H/D atoms on the grain surface (H_p + H_p -> H2, H_p + D_p ->
    // HD), whose intermediate species the compact network does not carry.  The
    // closed-form Cazaux & Tielens surface-formation coefficients
    // (kgr[127] = H+H->H2, kgr[128] = H+D->HD) reproduce that formation rate
    // directly from the gas-phase H/D abundances, so the compact model routes
    // them into the special grain reactions 2H+grain->H2 and 2D+grain->HD whose
    // forward rate is k * y[H or D] * nH (see compute_rates special branch).
    // Without this the gathered coefficient is the gas-phase k_rxn[22]/[143],
    // which are identically zero, leaving H2 grain formation absent.
    k_rxn[metal_grain_minimal::special::rxn_2H_grain - 1] =
        k_full[metal_grain::slot::grain_surface_begin + 127];
    k_rxn[metal_grain_minimal::special::rxn_2D_grain - 1] =
        k_full[metal_grain::slot::grain_surface_begin + 128];
    return;
  }

  // Generic full-network path (primordial and full metal).  Guarded so it is
  // not instantiated for the compact metal model, whose has_grain branches
  // below would otherwise index full-network slots that do not exist in the
  // compact table; that model returns from the gather branch above.
  if constexpr (!Model::is_compact_metal) {
  // Species-independent intermediates (temperature powers, LTE critical
  // densities, LTE/low-density blend sub-rates) built once and shared by the
  // forward rate laws and the reverse / LTE-correction step below.
  PrimRateContext c = build_prim_context(T_K, nH, params.H, params.H2,
                                         params.He, tbl.reactions[21].delE);
  double T300 = c.T300;
  double crit_ratio = c.crit_ratio;
  double crit_ratio_HD = c.crit_ratio_HD;

  // rho — parenthesization preserved per model for FP bit identity
  double rho;
  if constexpr (Model::has_grain) {
    rho = nH * (1.0 + 4.0 * yHe) * m_p;
  } else {
    rho = nH * ((1.0 + 4.0 * yHe) * m_p);
  }

  double lmbd_J = std::sqrt(pi * k_B * T_K / (G * mu * m_p * rho));
  double kap = detail::eval_opacity(T_K, rho);

  double tau_cnt;
  if constexpr (Model::has_grain) {
    double k_gr = detail::kp_gr(rho, params.T_gr_K) * params.Z_metal;
    tau_cnt = (kap + k_gr) * rho * lmbd_J;
  } else {
    tau_cnt = kap * rho * lmbd_J;
  }

  std::array<double, N_sp + 3> pf;
  if constexpr (Model::has_grain) {
    eval_metal_partition_functions<N_sp, N_react>(T_K, tbl, pf);
  } else {
    pf_prim::eval_pf_set<typename Model::Species, N_react>(T_K, tbl, pf);
  }

  // -----------------------------------------------------------------------
  // H/He/D forward rates (k_rxn[0..99], shared GA08 network)
  // -----------------------------------------------------------------------
  compute_HHeD_rates<Model::rate_form, N_react>(k_rxn, c);

  // Zero-forcing (common H/He)
  k_rxn[10] = 0.0;
  k_rxn[12] = 0.0;
  k_rxn[19] = 0.0;
  if constexpr (Model::has_grain) {
    k_rxn[24] = 0.0;  // metal zeroes k_rxn[24], not k_rxn[30]
  } else {
    k_rxn[30] = 0.0;  // primordial zeroes k_rxn[30], not k_rxn[24]
  }
  k_rxn[39] = 0.0;
  k_rxn[40] = 0.0;
  k_rxn[42] = 0.0;

  // Primordial k_rxn[13] threshold
  if constexpr (!Model::has_grain) {
    if (phys::k_B_eV * T_K <= 0.04) {
      k_rxn[rxn::Hm_e_to_H_2e] = 0.0;
      k_rxn[78] = 0.0;
    }
  }

  // Zero-forcing (D)
  k_rxn[52] = 0.0;
  k_rxn[57] = 0.0;
  k_rxn[58] = 0.0;
  k_rxn[66] = 0.0;
  k_rxn[70] = 0.0;
  k_rxn[73] = 0.0;
  k_rxn[74] = 0.0;
  k_rxn[92] = 0.0;
  k_rxn[93] = 0.0;
  k_rxn[94] = 0.0;
  k_rxn[95] = 0.0;
  k_rxn[96] = 0.0;

  // Blend sub-rates shared with the reverse / LTE-correction step below.
  double k12LTE = c.k12LTE;
  double k13v0 = c.k13v0;
  double k21LTE = c.k21LTE;
  double k37LTE = c.k37LTE;
  double k40v0 = c.k40v0;
  double k75v0 = c.k75v0;

  // -----------------------------------------------------------------------
  // Model-specific forward rates beyond H/He/D
  // -----------------------------------------------------------------------
  if constexpr (Model::has_grain) {
    // Metal reactions k_rxn[100..542] + k_rxn[600..644]
    compute_metal_rates<N_react>(k_rxn, T_K, T300);
    // Li reactions k_rxn[800..829]
    compute_Li_rates<N_react>(k_rxn, T_K, T300, tbl.reactions[634].delE);
    // K, Na, Mg reactions k_rxn[700..729]
    compute_KNaMg_rates<N_react>(k_rxn, T_K, T300, tbl.reactions[576].delE,
                                 tbl.reactions[591].delE);
  } else {
    // Primordial Li reactions k_rxn[100..129] (bind reaction slots to laws)
    namespace rl = ratelaw::prim;
    k_rxn[100] = rl::r100_e_Lip_to_Li(c);
    k_rxn[101] = rl::r101_Hm_Lip_to_H_Li(c);
    k_rxn[102] = rl::r102_Hp_Lim_to_H_Li(c);
    k_rxn[103] = rl::r103_e_Li_to_Lim(c);
    k_rxn[104] = rl::r104_Hp_Li_to_H_Lip(c);
    k_rxn[105] = rl::r105_Hp_Li_to_H_Lip_ph(c);
    k_rxn[106] = rl::r106_Hm_Li_to_e_LiH(c);
    k_rxn[107] = rl::r107_H_Lim_to_e_LiH(c);
    k_rxn[108] = rl::r108_H_LiHp_to_Hp_LiH(c);
    k_rxn[109] = rl::r109_Hp_LiH_to_H_LiHp(c);
    k_rxn[110] = rl::r110_H_LiH_to_H2_Li(c);
    k_rxn[111] = rl::r111_H_Lip_to_LiHp(c);
    k_rxn[112] = rl::r112_Hp_Li_to_LiHp(c);
    k_rxn[113] = rl::r113_Hp_LiH_to_H2_Lip(c);
    k_rxn[114] = rl::r114_e_LiHp_to_H_Li(c);
    k_rxn[115] = rl::r115_H_LiHp_to_H2p_Li(c);
    k_rxn[116] = rl::r116_H_LiHp_to_H2_Lip(c);
    k_rxn[117] = rl::r117_H_Li_to_LiH(c);
    k_rxn[118] = rl::r118_H2p_Li_to_Hp_LiH(c);
    k_rxn[119] = rl::r119_Hp_LiH_to_H2p_Li(c);
    k_rxn[120] = rl::r120_e_Lipp_to_Lip(c);
    k_rxn[121] = rl::r121_e_Lippp_to_Lipp(c);
    k_rxn[122] = rl::r122_Dm_Lip_to_D_Li(c);
    k_rxn[123] = rl::r123_Dp_Lim_to_D_Li(c);
    k_rxn[124] = rl::r124_e_Li_to_2e_Lip(c);
    k_rxn[125] = rl::r125_e_Lip_to_2e_Lipp(c);
    k_rxn[126] = rl::r126_e_Lipp_to_2e_Lippp(c);
    k_rxn[127] = rl::r127_2H_Li_to_H_LiH(c);
    k_rxn[128] = rl::r128_H_H2_Li_to_H2_LiH(c);
    k_rxn[129] = rl::r129_H2_Li_to_H2_e_Lip(c, tbl.reactions[129].delE);

    // Zero out disabled Li reactions
    k_rxn[108] = 0.0;
    k_rxn[119] = 0.0;
  }

  // -----------------------------------------------------------------------
  // Reverse rates via detailed balance
  // -----------------------------------------------------------------------
  std::array<double, N_react> lnKeqb{};
  if constexpr (Model::has_grain) {
    compute_reverse_loop<RateForm::MetalGrain, N_sp, N_react>(
        k_rxn, T_K, tau_cnt, metal_grain::loop::n_std, tbl, pf, lnKeqb);
  } else {
    compute_reverse_loop<RateForm::Primordial, N_sp, N_react>(
        k_rxn, T_K, tau_cnt, zero_metal::loop::n_std, tbl, pf, lnKeqb);
  }

  // -----------------------------------------------------------------------
  // Post-reverse LTE/low-density corrections
  // -----------------------------------------------------------------------
  if constexpr (Model::has_grain) {
    // ─── Metal LTE corrections ───────────────────────────────────────────
    const double lnC1_base = std::log(2.0 * pi * k_B * T_K / (h_P * h_P));

    // Reaction 273 special case (4-body)
    {
      constexpr double Cm273 = 1.01e71;
      constexpr double dE273 = 7.80e-12;
      double lnC1_273 = -3.0 * lnC1_base;
      double lnCpf_273 = std::log(pf[8]) + std::log(pf[21]) - std::log(pf[0]) -
                         std::log(pf[1]) - std::log(pf[7]) - std::log(pf[24]);
      double lnK273 =
          lnC1_273 + std::log(Cm273) + lnCpf_273 - dE273 / (k_B * T_K);
      k_rxn[272 + N_react] = std::exp(std::log(k_rxn[272]) + lnK273);
    }

    // k_rxn(8+N_react): H- + H -> H2 + e reverse
    {
      double k40LTE = std::exp(std::log(k_rxn[rxn::Hm_H_to_H2_e]) +
                               lnKeqb[rxn::Hm_H_to_H2_e]);
      k_rxn[rxn::Hm_H_to_H2_e + N_react] = lte_blend(k40v0, k40LTE, crit_ratio);
    }

    // k_rxn(12+N_react)
    {
      double lnk12rev = std::log(k12LTE) + lnKeqb[11];
      k_rxn[11 + N_react] = std::exp(lnk12rev);
    }

    // k_rxn(19+N_react)
    {
      double k13LTE_rev =
          std::exp(std::log(k_rxn[rxn::three_H]) + lnKeqb[rxn::three_H]);
      k_rxn[rxn::three_H + N_react] = lte_blend(k13v0, k13LTE_rev, crit_ratio);
    }

    // k_rxn(21+N_react)
    {
      double lnk21rev = std::log(k21LTE) + lnKeqb[rxn::H2H2_dis];
      k_rxn[rxn::H2H2_dis + N_react] = std::exp(lnk21rev);
    }

    // k_rxn(37+N_react)
    {
      double lnk37rev = std::log(k37LTE) + lnKeqb[rxn::H2_He_dis];
      k_rxn[rxn::H2_He_dis + N_react] = std::exp(lnk37rev);
    }

    // k_rxn(47+N_react)
    {
      double k47LTE_val = std::exp(std::log(k_rxn[46]) + lnKeqb[46]);
      double k47v0_val =
          5.09e-9 * std::pow(T_K, 0.128) * std::exp(-103258.0 / T_K);
      k_rxn[46 + N_react] = lte_blend(k47v0_val, k47LTE_val, crit_ratio_HD);
    }

    // k_rxn(48+N_react)
    {
      double lnk48rev = std::log(k37LTE) + lnKeqb[47];
      k_rxn[47 + N_react] = std::exp(lnk48rev);
    }

    // k_rxn(49+N_react)
    {
      double lnk49rev = std::log(k21LTE) + lnKeqb[48];
      k_rxn[48 + N_react] = std::exp(lnk49rev);
    }

    // k_rxn(50+N_react)
    {
      double k13LTE_db =
          std::exp(std::log(k_rxn[rxn::three_H]) + lnKeqb[rxn::three_H]);
      double lnk50rev = std::log(k13LTE_db) + lnKeqb[49];
      k_rxn[49 + N_react] = std::exp(lnk50rev);
    }

    // k_rxn(68+N_react)
    {
      double k75LTE_val = std::exp(std::log(k_rxn[67]) + lnKeqb[67]);
      k_rxn[67 + N_react] = lte_blend(k75v0, k75LTE_val, crit_ratio_HD);
    }

    // Suppress He ionization reversal at low crit_ratio
    if (crit_ratio < 1.0) {
      k_rxn[0 + N_react] = 0.0;
      k_rxn[2 + N_react] = 0.0;
    }

    // T <= 300K: zero specific reverse rates
    if (T_K <= 300.0) {
      k_rxn[101 + N_react] = 0.0;
      k_rxn[106 + N_react] = 0.0;
      k_rxn[108 + N_react] = 0.0;
      k_rxn[113 + N_react] = 0.0;
      k_rxn[133 + N_react] = 0.0;
      k_rxn[137 + N_react] = 0.0;
      k_rxn[138 + N_react] = 0.0;
      k_rxn[140 + N_react] = 0.0;
      k_rxn[141 + N_react] = 0.0;
      k_rxn[174 + N_react] = 0.0;
      k_rxn[181 + N_react] = 0.0;
      k_rxn[229 + N_react] = 0.0;
      k_rxn[539 + N_react] = 0.0;
      k_rxn[602 + N_react] = 0.0;
    }

    // Grain rates
    compute_grain_rates<N_react>(k_rxn, T_K, params.T_gr_K, rho, nH,
                                 params.Z_metal, params.J_H2, params.J_H2O,
                                 params.J_tot, params.zeta, params.T_cr_desorp);

    // CR rates
    compute_CR_rates_metal<N_react>(k_rxn, params.zeta, T300);

  } else {
    // ─── Primordial LTE corrections ──────────────────────────────────────
    // k_rxn(8+N_react)
    {
      double lnk40LTE =
          std::log(k_rxn[rxn::Hm_H_to_H2_e]) + lnKeqb[rxn::Hm_H_to_H2_e];
      double k40LTE = std::exp(lnk40LTE);
      k_rxn[rxn::Hm_H_to_H2_e + N_react] =
          lte_blend_noguard(k40v0, k40LTE, crit_ratio);
    }

    // k_rxn(12+N_react)
    {
      double lnk12rev = std::log(k12LTE) + lnKeqb[11];
      k_rxn[11 + N_react] = std::exp(lnk12rev);
    }

    // k_rxn[18+N_react] — k13LTE declared at outer scope for k_rxn[49]
    double k13LTE;
    {
      double lnk13LTE = std::log(k_rxn[rxn::three_H]) + lnKeqb[rxn::three_H];
      k13LTE = std::exp(lnk13LTE);
      k_rxn[rxn::three_H + N_react] =
          lte_blend_noguard(k13v0, k13LTE, crit_ratio);
    }

    // k_rxn(21+N_react)
    {
      double lnk21rev = std::log(k21LTE) + lnKeqb[rxn::H2H2_dis];
      k_rxn[rxn::H2H2_dis + N_react] = std::exp(lnk21rev);
    }

    // k_rxn(37+N_react)
    {
      double lnk37rev = std::log(k37LTE) + lnKeqb[rxn::H2_He_dis];
      k_rxn[rxn::H2_He_dis + N_react] = std::exp(lnk37rev);
    }

    // k_rxn(47+N_react)
    {
      double k47LTE = c.k47LTE;
      double lnk47rev = std::log(k47LTE) + lnKeqb[46];
      k_rxn[46 + N_react] = std::exp(lnk47rev);
    }

    // k_rxn(48+N_react)
    {
      double lnk48rev = std::log(k37LTE) + lnKeqb[47];
      k_rxn[47 + N_react] = std::exp(lnk48rev);
    }

    // k_rxn(49+N_react)
    {
      double lnk49rev = std::log(k21LTE) + lnKeqb[48];
      k_rxn[48 + N_react] = std::exp(lnk49rev);
    }

    // k_rxn(50+N_react)
    {
      double lnk50rev = std::log(k13LTE) + lnKeqb[49];
      k_rxn[49 + N_react] = std::exp(lnk50rev);
    }

    // k_rxn(68+N_react)
    {
      double lnk75LTE = std::log(k_rxn[67]) + lnKeqb[67];
      double k75LTE = std::exp(lnk75LTE);
      k_rxn[67 + N_react] = lte_blend_noguard(k75v0, k75LTE, crit_ratio_HD);
    }

    // Suppress He ionization reversal at low crit_ratio
    if (crit_ratio < 1.0) {
      k_rxn[0 + N_react] = 0.0;
      k_rxn[2 + N_react] = 0.0;
    }

    // CR reactions k_rxn[130..138]
    compute_CR_rates_prim<N_react>(k_rxn, params.zeta);
  }
  }  // if constexpr (!Model::is_compact_metal)
}

// ---------------------------------------------------------------------------
// compute_rates<Model>
//
// Computes dy/dt = r_f[i] and Jacobian dr_fdy[i][j] for the NR solver.
// ---------------------------------------------------------------------------
template <class Model>
void compute_rates(const std::array<double, 2 * Model::N_react>& k_rxn,
                   double nH, const std::array<double, Model::N_sp>& y,
                   const ReactionTable<Model::N_sp, Model::N_react>& tbl,
                   std::array<double, Model::N_sp>& r_f,
                   std::array<double, Model::N_sp * Model::N_sp>& dr_fdy,
                   std::array<double, 2 * Model::N_react>& var) {
  constexpr int N_sp = Model::N_sp;
  constexpr int N_react = Model::N_react;
  using Sp = typename Model::Sp;

  // Special-form grain reaction ids and the grain solver-block sizes, selected
  // per model: the compact metal Minimal network renumbers these reactions and
  // carries smaller blocks than the full metal_grain network.  The ternaries
  // resolve at compile time; for non-grain models the metal branch is unused.
  constexpr int sp_2H_grain = Model::is_compact_metal
                                  ? metal_grain_minimal::special::rxn_2H_grain
                                  : metal_grain::special::rxn_2H_grain;
  constexpr int sp_2D_grain = Model::is_compact_metal
                                  ? metal_grain_minimal::special::rxn_2D_grain
                                  : metal_grain::special::rxn_2D_grain;
  constexpr int sp_4body_H = Model::is_compact_metal
                                 ? metal_grain_minimal::special::rxn_4body_H
                                 : metal_grain::special::rxn_4body_H;
  // D index for the 2D+grain special (full layout 11, compact 9).  Guarded so
  // Sp::D is only required of grain models.
  constexpr int idx_D_grain = [] {
    if constexpr (Model::has_grain)
      return static_cast<int>(Sp::D);
    else
      return 11;
  }();

  r_f.fill(0.0);
  dr_fdy.fill(0.0);
  var.fill(0.0);

  constexpr int N_ext = N_sp + 3;
  // Extended species vector with vacant/photon/CR sentinels at N_sp..N_sp+2
  // (see make_y_ext); sentinels are 1.0 so absent reactants/products factor
  // out.
  std::array<double, N_ext> y_ext = make_y_ext<N_sp>(y);

  std::array<double, N_ext> r_f_dum{};
  r_f_dum.fill(0.0);
  std::array<double, N_ext * N_ext> dr_fdy_dum{};
  dr_fdy_dum.fill(0.0);

  auto dJ = [&](int i, int j) -> double& { return dr_fdy_dum[i * N_ext + j]; };

  // ─── Loop 1: standard reactions ──────────────────────────────────────────
  int n_std_loop;
  if constexpr (Model::is_compact_metal) {
    n_std_loop = metal_grain_minimal::loop::n_std;  // 66 compact standard
  } else if constexpr (Model::has_grain) {
    n_std_loop = metal_grain::loop::n_std;
  } else if constexpr (Model::is_compact_prim) {
    n_std_loop = zero_metal_minimal::loop::n_std;  // 24 standard; CR in loop 2
  } else {
    n_std_loop = zero_metal::loop::n_std;
  }

  for (int ire = 0; ire < n_std_loop; ++ire) {
    int num = tbl.reactions[ire].num;
    if (num < 1 || num > N_react) continue;
    int r1 = tbl.reactions[ire].reactants[0],
        r2 = tbl.reactions[ire].reactants[1],
        r3 = tbl.reactions[ire].reactants[2];
    int p1 = tbl.reactions[ire].products[0],
        p2 = tbl.reactions[ire].products[1],
        p3 = tbl.reactions[ire].products[2];
    int nr = tbl.reactions[ire].n_reactants, np = tbl.reactions[ire].n_products;

    double nHf = pow_int(nH, nr - 1);

    // Forward rate
    double rate_fwd;
    if constexpr (Model::has_grain) {
      if (num == sp_2H_grain) {
        rate_fwd = k_rxn[sp_2H_grain - 1] * y_ext[0] * nH;
      } else if (num == sp_2D_grain) {
        rate_fwd = k_rxn[sp_2D_grain - 1] * y_ext[idx_D_grain] * nH;
      } else {
        rate_fwd = k_rxn[num - 1] * y_ext[r1] * y_ext[r2] * y_ext[r3] * nHf;
      }
    } else {
      rate_fwd = k_rxn[num - 1] * y_ext[r1] * y_ext[r2] * y_ext[r3] * nHf;
    }

    // Reverse rate
    double rate_rev, nHr;
    if constexpr (Model::has_grain) {
      if (num == sp_4body_H) {
        nHr = pow_int(nH, np);
        rate_rev = k_rxn[num - 1 + N_react] * y_ext[p1] * y_ext[p2] *
                   y_ext[p3] * y_ext[0] * nHr;
      } else {
        nHr = pow_int(nH, np - 1);
        rate_rev =
            k_rxn[num - 1 + N_react] * y_ext[p1] * y_ext[p2] * y_ext[p3] * nHr;
      }
    } else {
      nHr = pow_int(nH, np - 1);
      rate_rev =
          k_rxn[num - 1 + N_react] * y_ext[p1] * y_ext[p2] * y_ext[p3] * nHr;
    }

    double rate = rate_fwd - rate_rev;

    r_f_dum[r1] -= rate;
    r_f_dum[r2] -= rate;
    r_f_dum[r3] -= rate;
    r_f_dum[p1] += rate;
    r_f_dum[p2] += rate;
    r_f_dum[p3] += rate;
    if constexpr (Model::has_grain) {
      if (num == sp_4body_H) r_f_dum[0] += rate;
    }

    var[num - 1] = rate_fwd;
    var[num - 1 + N_react] = rate_rev;

    // ── Jacobian: forward contributions ──────────────────────────────────
    if constexpr (Model::has_grain) {
      if (num == sp_2H_grain) {
        double d1 = k_rxn[sp_2H_grain - 1] * nH;
        dJ(r1, 0) -= d1;
        dJ(r2, 0) -= d1;
        dJ(r3, 0) -= d1;
        dJ(p1, 0) += d1;
        dJ(p2, 0) += d1;
        dJ(p3, 0) += d1;
      } else if (num == sp_2D_grain) {
        double d12 = k_rxn[sp_2D_grain - 1] * nH;
        dJ(r1, idx_D_grain) -= d12;
        dJ(r2, idx_D_grain) -= d12;
        dJ(r3, idx_D_grain) -= d12;
        dJ(p1, idx_D_grain) += d12;
        dJ(p2, idx_D_grain) += d12;
        dJ(p3, idx_D_grain) += d12;
      } else {
        double fJ_r1 = k_rxn[num - 1] * y_ext[r2] * y_ext[r3] * nHf;
        dJ(r1, r1) -= fJ_r1;
        dJ(r2, r1) -= fJ_r1;
        dJ(r3, r1) -= fJ_r1;
        dJ(p1, r1) += fJ_r1;
        dJ(p2, r1) += fJ_r1;
        dJ(p3, r1) += fJ_r1;
        if (num == sp_4body_H) dJ(0, r1) += fJ_r1;

        double fJ_r2 = k_rxn[num - 1] * y_ext[r1] * y_ext[r3] * nHf;
        dJ(r1, r2) -= fJ_r2;
        dJ(r2, r2) -= fJ_r2;
        dJ(r3, r2) -= fJ_r2;
        dJ(p1, r2) += fJ_r2;
        dJ(p2, r2) += fJ_r2;
        dJ(p3, r2) += fJ_r2;
        if (num == sp_4body_H) dJ(0, r2) += fJ_r2;

        double fJ_r3 = k_rxn[num - 1] * y_ext[r1] * y_ext[r2] * nHf;
        dJ(r1, r3) -= fJ_r3;
        dJ(r2, r3) -= fJ_r3;
        dJ(r3, r3) -= fJ_r3;
        dJ(p1, r3) += fJ_r3;
        dJ(p2, r3) += fJ_r3;
        dJ(p3, r3) += fJ_r3;
        if (num == sp_4body_H) dJ(0, r3) += fJ_r3;
      }
    } else {
      double dr1v = k_rxn[num - 1] * y_ext[r2] * y_ext[r3] * nHf;
      dJ(r1, r1) -= dr1v;
      dJ(r2, r1) -= dr1v;
      dJ(r3, r1) -= dr1v;
      dJ(p1, r1) += dr1v;
      dJ(p2, r1) += dr1v;
      dJ(p3, r1) += dr1v;

      double dr2v = k_rxn[num - 1] * y_ext[r1] * y_ext[r3] * nHf;
      dJ(r1, r2) -= dr2v;
      dJ(r2, r2) -= dr2v;
      dJ(r3, r2) -= dr2v;
      dJ(p1, r2) += dr2v;
      dJ(p2, r2) += dr2v;
      dJ(p3, r2) += dr2v;

      double dr3v = k_rxn[num - 1] * y_ext[r1] * y_ext[r2] * nHf;
      dJ(r1, r3) -= dr3v;
      dJ(r2, r3) -= dr3v;
      dJ(r3, r3) -= dr3v;
      dJ(p1, r3) += dr3v;
      dJ(p2, r3) += dr3v;
      dJ(p3, r3) += dr3v;
    }

    // ── Jacobian: reverse contributions ──────────────────────────────────
    if constexpr (Model::has_grain) {
      if (num == sp_4body_H) {
        double rJ_p1 =
            k_rxn[num - 1 + N_react] * y_ext[p2] * y_ext[p3] * y_ext[0] * nHr;
        dJ(r1, p1) += rJ_p1;
        dJ(r2, p1) += rJ_p1;
        dJ(r3, p1) += rJ_p1;
        dJ(p1, p1) -= rJ_p1;
        dJ(p2, p1) -= rJ_p1;
        dJ(p3, p1) -= rJ_p1;
        dJ(0, p1) -= rJ_p1;

        double rJ_p2 =
            k_rxn[num - 1 + N_react] * y_ext[p1] * y_ext[p3] * y_ext[0] * nHr;
        dJ(r1, p2) += rJ_p2;
        dJ(r2, p2) += rJ_p2;
        dJ(r3, p2) += rJ_p2;
        dJ(p1, p2) -= rJ_p2;
        dJ(p2, p2) -= rJ_p2;
        dJ(p3, p2) -= rJ_p2;
        dJ(0, p2) -= rJ_p2;

        double rJ_p3 =
            k_rxn[num - 1 + N_react] * y_ext[p1] * y_ext[p2] * y_ext[0] * nHr;
        dJ(r1, p3) += rJ_p3;
        dJ(r2, p3) += rJ_p3;
        dJ(r3, p3) += rJ_p3;
        dJ(p1, p3) -= rJ_p3;
        dJ(p2, p3) -= rJ_p3;
        dJ(p3, p3) -= rJ_p3;
        dJ(0, p3) -= rJ_p3;

        double rJ_1 =
            k_rxn[num - 1 + N_react] * y_ext[p1] * y_ext[p2] * y_ext[p3] * nHr;
        dJ(r1, 0) += rJ_1;
        dJ(r2, 0) += rJ_1;
        dJ(r3, 0) += rJ_1;
        dJ(p1, 0) -= rJ_1;
        dJ(p2, 0) -= rJ_1;
        dJ(p3, 0) -= rJ_1;
        dJ(0, 0) -= rJ_1;
      } else {
        double rJ_p1 = k_rxn[num - 1 + N_react] * y_ext[p2] * y_ext[p3] * nHr;
        dJ(r1, p1) += rJ_p1;
        dJ(r2, p1) += rJ_p1;
        dJ(r3, p1) += rJ_p1;
        dJ(p1, p1) -= rJ_p1;
        dJ(p2, p1) -= rJ_p1;
        dJ(p3, p1) -= rJ_p1;

        double rJ_p2 = k_rxn[num - 1 + N_react] * y_ext[p1] * y_ext[p3] * nHr;
        dJ(r1, p2) += rJ_p2;
        dJ(r2, p2) += rJ_p2;
        dJ(r3, p2) += rJ_p2;
        dJ(p1, p2) -= rJ_p2;
        dJ(p2, p2) -= rJ_p2;
        dJ(p3, p2) -= rJ_p2;

        double rJ_p3 = k_rxn[num - 1 + N_react] * y_ext[p1] * y_ext[p2] * nHr;
        dJ(r1, p3) += rJ_p3;
        dJ(r2, p3) += rJ_p3;
        dJ(r3, p3) += rJ_p3;
        dJ(p1, p3) -= rJ_p3;
        dJ(p2, p3) -= rJ_p3;
        dJ(p3, p3) -= rJ_p3;
      }
    } else {
      double dp1v = k_rxn[num - 1 + N_react] * y_ext[p2] * y_ext[p3] * nHr;
      dJ(r1, p1) += dp1v;
      dJ(r2, p1) += dp1v;
      dJ(r3, p1) += dp1v;
      dJ(p1, p1) -= dp1v;
      dJ(p2, p1) -= dp1v;
      dJ(p3, p1) -= dp1v;

      double dp2v = k_rxn[num - 1 + N_react] * y_ext[p1] * y_ext[p3] * nHr;
      dJ(r1, p2) += dp2v;
      dJ(r2, p2) += dp2v;
      dJ(r3, p2) += dp2v;
      dJ(p1, p2) -= dp2v;
      dJ(p2, p2) -= dp2v;
      dJ(p3, p2) -= dp2v;

      double dp3v = k_rxn[num - 1 + N_react] * y_ext[p1] * y_ext[p2] * nHr;
      dJ(r1, p3) += dp3v;
      dJ(r2, p3) += dp3v;
      dJ(r3, p3) += dp3v;
      dJ(p1, p3) -= dp3v;
      dJ(p2, p3) -= dp3v;
      dJ(p3, p3) -= dp3v;
    }
  }

  // ─── Loop 2: CR first-order reactions ────────────────────────────────────
  int cr_loop_begin, cr_loop_end;
  if constexpr (Model::is_compact_metal) {
    cr_loop_begin = metal_grain_minimal::loop::n_std;
    cr_loop_end = metal_grain_minimal::loop::n_std +
                  metal_grain_minimal::loop::n_cr_react;
  } else if constexpr (Model::has_grain) {
    cr_loop_begin = metal_grain::loop::n_std;
    cr_loop_end = metal_grain::loop::n_std + metal_grain::loop::n_cr_react;
  } else if constexpr (Model::is_compact_prim) {
    cr_loop_begin = zero_metal_minimal::loop::n_std;
    cr_loop_end = zero_metal_minimal::loop::n_std +
                  zero_metal_minimal::loop::n_cr_react;
  } else {
    cr_loop_begin = zero_metal::loop::n_std;
    cr_loop_end = zero_metal::loop::n_std + zero_metal::loop::n_cr_react;
  }

  for (int ire = cr_loop_begin; ire < cr_loop_end; ++ire) {
    int num = tbl.reactions[ire].num;
    if (num < 1 || num > N_react) continue;
    int r2 = tbl.reactions[ire].reactants[1];
    int p1 = tbl.reactions[ire].products[0],
        p2 = tbl.reactions[ire].products[1],
        p3 = tbl.reactions[ire].products[2];

    double rate = k_rxn[num - 1] * y_ext[r2];

    r_f_dum[r2] -= rate;
    r_f_dum[p1] += rate;
    r_f_dum[p2] += rate;
    r_f_dum[p3] += rate;

    var[num - 1] = rate;

    double d2 = k_rxn[num - 1];
    dJ(r2, r2) -= d2;
    dJ(p1, r2) += d2;
    dJ(p2, r2) += d2;
    dJ(p3, r2) += d2;
  }

  // ─── Loop 3 & 4: metal-only (charge + grain surface) ─────────────────────
  if constexpr (Model::has_grain) {
    // Gas-phase solver-block sizes, compact-aware (the Minimal network carries
    // smaller standard / CR / charge blocks than the full metal_grain network).
    constexpr int g_n_std = Model::is_compact_metal
                                ? metal_grain_minimal::loop::n_std
                                : metal_grain::loop::n_std;
    constexpr int g_n_cr = Model::is_compact_metal
                               ? metal_grain_minimal::loop::n_cr_react
                               : metal_grain::loop::n_cr_react;
    constexpr int g_n_charge = Model::is_compact_metal
                                   ? metal_grain_minimal::loop::n_charge_react
                                   : metal_grain::loop::n_charge_react;
    // Fallback anchor for the grain-surface num run: the start of the grain
    // band (full network: 990; compact: the gas-phase block width N_gas).
    constexpr int g_surface_begin = Model::is_compact_metal
                                        ? metal_grain_minimal::N_gas
                                        : metal_grain::slot::grain_surface_begin;

    // Loop 3: gas-grain charge transfer
    for (int ire = g_n_std + g_n_cr; ire < g_n_std + g_n_cr + g_n_charge;
         ++ire) {
      int num = tbl.reactions[ire].num;
      if (num < 1 || num > N_react) continue;
      int r1 = tbl.reactions[ire].reactants[0],
          r2 = tbl.reactions[ire].reactants[1],
          r3 = tbl.reactions[ire].reactants[2];
      int p1 = tbl.reactions[ire].products[0],
          p2 = tbl.reactions[ire].products[1],
          p3 = tbl.reactions[ire].products[2];
      int nr = tbl.reactions[ire].n_reactants;

      double nHf = pow_int(nH, nr - 1);
      double rate = k_rxn[num - 1] * y_ext[r1] * y_ext[r2] * y_ext[r3] * nHf;

      r_f_dum[r1] -= rate;
      r_f_dum[r2] -= rate;
      r_f_dum[r3] -= rate;
      r_f_dum[p1] += rate;
      r_f_dum[p2] += rate;
      r_f_dum[p3] += rate;

      var[num - 1] = rate;

      double fJ_r1 = k_rxn[num - 1] * y_ext[r2] * y_ext[r3] * nHf;
      dJ(r1, r1) -= fJ_r1;
      dJ(r2, r1) -= fJ_r1;
      dJ(r3, r1) -= fJ_r1;
      dJ(p1, r1) += fJ_r1;
      dJ(p2, r1) += fJ_r1;
      dJ(p3, r1) += fJ_r1;

      double fJ_r2 = k_rxn[num - 1] * y_ext[r1] * y_ext[r3] * nHf;
      dJ(r1, r2) -= fJ_r2;
      dJ(r2, r2) -= fJ_r2;
      dJ(r3, r2) -= fJ_r2;
      dJ(p1, r2) += fJ_r2;
      dJ(p2, r2) += fJ_r2;
      dJ(p3, r2) += fJ_r2;

      double fJ_r3 = k_rxn[num - 1] * y_ext[r1] * y_ext[r2] * nHf;
      dJ(r1, r3) -= fJ_r3;
      dJ(r2, r3) -= fJ_r3;
      dJ(r3, r3) -= fJ_r3;
      dJ(p1, r3) += fJ_r3;
      dJ(p2, r3) += fJ_r3;
      dJ(p3, r3) += fJ_r3;
    }

    // Loop 4: grain surface reactions.  Their forward rates live in the grain
    // band above the gas-phase block (compact: k_rxn[N_gas..N_gas+N_grain);
    // full: k_rxn[990..]), addressed by the run of nums following the last
    // gas-phase reaction.
    {
      constexpr int n_pre_surface = g_n_std + g_n_cr + g_n_charge;
      const int num_start = (tbl.n_loaded >= n_pre_surface &&
                             tbl.reactions[n_pre_surface - 1].num > 0)
                                ? tbl.reactions[n_pre_surface - 1].num
                                : g_surface_begin;
      for (int ire = 0; ire < tbl.n_grain; ++ire) {
        int num = num_start + 1 + ire;
        if (num < 1 || num > N_react) continue;
        int r1 = tbl.grain_reactions[ire].reactants[0],
            r2 = tbl.grain_reactions[ire].reactants[1];
        int p1 = tbl.grain_reactions[ire].products[0],
            p2 = tbl.grain_reactions[ire].products[1];
        int nr = tbl.grain_reactions[ire].n_reactants;

        double nHr = pow_int(nH, nr);
        double rate = k_rxn[num - 1] * y_ext[r1] * y_ext[r2] * nHr;

        r_f_dum[r1] -= rate;
        r_f_dum[r2] -= rate;
        r_f_dum[p1] += rate;
        r_f_dum[p2] += rate;

        var[num - 1] = rate;

        double fJ_r1 = k_rxn[num - 1] * y_ext[r2] * nHr;
        dJ(r1, r1) -= fJ_r1;
        dJ(r2, r1) -= fJ_r1;
        dJ(p1, r1) += fJ_r1;
        dJ(p2, r1) += fJ_r1;

        double fJ_r2 = k_rxn[num - 1] * y_ext[r1] * nHr;
        dJ(r1, r2) -= fJ_r2;
        dJ(r2, r2) -= fJ_r2;
        dJ(p1, r2) += fJ_r2;
        dJ(p2, r2) += fJ_r2;
      }
    }
  }

  // ── Copy [0..N_sp-1] from extended buffers ───────────────────────────────
  for (int i = 0; i < N_sp; ++i) r_f[i] = r_f_dum[i];

  for (int i = 0; i < N_sp; ++i)
    for (int j = 0; j < N_sp; ++j)
      dr_fdy[i * N_sp + j] = dr_fdy_dum[i * N_ext + j];
}

}  // namespace arche
