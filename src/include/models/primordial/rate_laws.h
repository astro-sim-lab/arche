// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// rate_laws.h
//
// Forward rate laws for the GA08 H/He/D primordial network plus the
// cosmic-ray reactions.
//
// Provides:
//   ratelaw::prim::*           — per-reaction forward rate laws, each a pure
//                                function of the PrimRateContext (free of any
//                                global k_rxn slot).  A model composes its
//                                network by binding reaction slots to laws.
//   compute_HHeD_rates(...)    — binds the H/He/D forward laws to k_rxn[0..99].
//   compute_CR_rates_prim(...) — CR rates for the primordial network
//                                (k_rxn[130..138]).
//
// The species-independent PrimRateContext, the LTE/low-density blends, and the
// detailed-balance reverse loop are shared primitives in kinetics/rates.h.
//
// Index convention: 0-based k_rxn[N] = 1-based reaction id N+1.
// All rate coefficients in CGS (cm³/s or cm⁶/s for 3-body).
// ---------------------------------------------------------------------------
#include <array>
#include <cmath>

#include "kinetics/reaction_index.h"
#include "core/state.h"
#include "kinetics/rates.h"

namespace arche {

// ---------------------------------------------------------------------------
// ratelaw::prim — forward rate laws for the GA08 H/He/D network.
//
// Each function returns one reaction's forward rate coefficient as a pure
// function of the PrimRateContext, free of any global k_rxn[] slot.  This is
// the "small call unit" a model composes from: a reaction binds to a law by
// name, and several reactions may bind to the same law (e.g. the D+ radiative
// recombination r50 reuses the H+ case-B law r1).  The names carry the 0-based
// k_rxn index for traceability; the laws themselves are index-independent.
//
// The Form template parameter selects between the two models' expression forms
// where they differ at ULP level (r3, r14, r21, r53).  Laws that reuse a
// Form-dependent law are Form-templated in turn (r63, r64, r81 reuse r14).
// ---------------------------------------------------------------------------
namespace ratelaw {
namespace prim {

// k_rxn[0] H + e -> H+ + 2e
inline double r0_H_e_to_Hp_2e(const PrimRateContext& c) {
  if (c.T_K < 1000.0) return 0.0;
  const double lnT_eV = c.lnT_eV;
  return std::exp(
      -32.71396786 +
      lnT_eV *
          (13.536556 +
           lnT_eV *
               (-5.73932875 +
                lnT_eV *
                    (1.56315498 +
                     lnT_eV *
                         (-0.2877056 +
                          lnT_eV *
                              (3.48255977e-2 +
                               lnT_eV *
                                   (-2.63197617e-3 +
                                    lnT_eV * (1.11954395e-4 +
                                              lnT_eV * (-2.03914985e-6)))))))));
}

// k_rxn[1] H+ + e -> H + ph.  (case B)
inline double r1_Hp_rec_caseB(const PrimRateContext& c) {
  return 2.753e-14 * std::pow(315614.0 / c.T_K, 1.500) *
         std::pow(1.0 + std::pow(115188.0 / c.T_K, 0.407), -2.242);
}

// k_rxn[2] He + e -> He+ + 2e
inline double r2_He_e_to_Hep_2e(const PrimRateContext& c) {
  if (c.T_K < 1000.0) return 0.0;
  const double lnT_eV = c.lnT_eV;
  return std::exp(
      -44.09864886 +
      lnT_eV *
          (23.91596563 +
           lnT_eV *
               (-10.7532302 +
                lnT_eV *
                    (3.05803875 +
                     lnT_eV *
                         (-0.56851189 +
                          lnT_eV *
                              (6.79539123e-2 +
                               lnT_eV *
                                   (-5.00905610e-3 +
                                    lnT_eV * (2.06723616e-4 +
                                              lnT_eV * (-3.64916141e-6)))))))));
}

// k_rxn[3] He+ + e -> He + ph.  (optically thick)
template <RateForm Form>
inline double r3_Hep_rec(const PrimRateContext& c) {
  const double T_K = c.T_K, lgT = c.lgT, lgT2 = c.lgT2, lgT3 = c.lgT3;
  double sqrt_factor;
  if constexpr (Form == RateForm::Primordial)
    sqrt_factor = 1.0e-11 / std::sqrt(T_K);
  else
    sqrt_factor = 1.0e-11 * std::pow(T_K, -0.5);
  double krrA =
      sqrt_factor * (12.72 - 1.615 * lgT - 0.3162 * lgT2 + 0.0493 * lgT3);
  double krrB =
      sqrt_factor * (11.19 - 1.676 * lgT - 0.2852 * lgT2 + 0.04433 * lgT3);
  double kdi = 1.9e-3 * std::pow(T_K, -1.5) * std::exp(-473421.0 / T_K) *
               (1.0 + 0.3 * std::exp(-94684.0 / T_K));
  return 0.68 * krrA + 0.32 * krrB + kdi;
}

// k_rxn[4] He+ + e -> He++ + 2e
inline double r4_Hep_e_to_Hepp_2e(const PrimRateContext& c) {
  if (c.T_K < 1000.0) return 0.0;
  const double lnT_eV = c.lnT_eV;
  return std::exp(
      -68.71040990 +
      lnT_eV *
          (43.93347633 +
           lnT_eV *
               (-18.4806699 +
                lnT_eV *
                    (4.70162649 +
                     lnT_eV *
                         (-0.76924663 +
                          lnT_eV *
                              (8.113042e-2 +
                               lnT_eV *
                                   (-5.32402063e-3 +
                                    lnT_eV * (1.97570531e-4 +
                                              lnT_eV * (-3.16558106e-6)))))))));
}

// k_rxn[5] He++ + e -> He+ + ph.  (case B)
inline double r5_Hepp_rec_caseB(const PrimRateContext& c) {
  return 5.506e-14 * std::pow(1262456.0 / c.T_K, 1.500) *
         std::pow(1.0 + std::pow(460752.0 / c.T_K, 0.407), -2.242);
}

// k_rxn[6] H + e -> H- + ph.
inline double r6_H_e_to_Hm(const PrimRateContext& c) {
  const double lgT = c.lgT, lgT2 = c.lgT2, lgT3 = c.lgT3, lgT4 = c.lgT4,
               lgT6 = c.lgT6;
  if (c.T_K < 6000.0)
    return std::pow(10.0,
                    -17.845 + 0.762 * lgT + 0.1523 * lgT2 - 0.03274 * lgT3);
  return std::pow(
      10.0, -16.4199 + 0.1998 * lgT2 - 5.447e-3 * lgT4 + 4.0415e-5 * lgT6);
}

// k_rxn[7] H- + H -> H2 + e  (Kreckel+2010)
inline double r7_Hm_H_to_H2_e(const PrimRateContext& c) {
  const double T_K = c.T_K;
  return 1.35e-9 *
         (std::pow(T_K, 9.8493e-2) + 3.2852e-1 * std::pow(T_K, 5.561e-1) +
          2.771e-7 * std::pow(T_K, 2.1826)) /
         (1.0 + 6.191e-3 * std::pow(T_K, 1.0461) +
          8.9712e-11 * std::pow(T_K, 3.0424) +
          3.2576e-14 * std::pow(T_K, 3.7741));
}

// k_rxn[8] H + H+ -> H2+ + ph.  (Coppola+2011, Glover+2015)
inline double r8_H_Hp_to_H2p(const PrimRateContext& c) {
  return std::pow(10.0,
                  -18.2 - 3.194 * c.lgT + 1.786 * c.lgT2 - 0.2072 * c.lgT3);
}

// k_rxn[9] H2+ + H -> H2 + H+
inline double r9_H2p_H_to_H2_Hp(const PrimRateContext&) { return 6.4e-10; }

// k_rxn[10] H2 + H+ -> H2+ + H  (GA08 7)
inline double r10_H2_Hp_to_H2p_H(const PrimRateContext& c) {
  const double T_K = c.T_K, lnT = c.lnT;
  if (T_K < 100.0) return 0.0;
  return (-3.3232183e-7 + 3.3735382e-7 * lnT - 1.4491368e-7 * lnT * lnT +
          3.4172805e-8 * std::pow(lnT, 3) - 4.7813720e-9 * std::pow(lnT, 4) +
          3.9731542e-10 * std::pow(lnT, 5) - 1.8171411e-11 * std::pow(lnT, 6) +
          3.5311932e-13 * std::pow(lnT, 7)) *
         std::exp(-21237.15 / T_K);
}

// k_rxn[11] H2 + e -> 2H + e  (GA08 8, LTE/low-density blend)
inline double r11_H2_e_to_2H_e(const PrimRateContext& c) {
  return lte_blend(c.k12v0, c.k12LTE, c.crit_ratio);
}

// k_rxn[12] H2 + H -> 3H  (GA08 9, LTE/low-density blend; inputs never 0)
inline double r12_H2_H_to_3H(const PrimRateContext& c) {
  return lte_blend_noguard(c.k13v0, c.k13LTE, c.crit_ratio);
}

// k_rxn[13] H- + e -> H + 2e  (GA08 14)
inline double r13_Hm_e_to_H_2e(const PrimRateContext& c) {
  if (c.T_K <= 40.0) return 0.0;
  const double lnT_eV = c.lnT_eV;
  return std::exp(
      -18.01849334 +
      lnT_eV *
          (2.3608522 +
           lnT_eV *
               (-0.28274430 +
                lnT_eV *
                    (1.62331664e-2 +
                     lnT_eV *
                         (-3.36501203e-2 +
                          lnT_eV *
                              (1.17832978e-2 +
                               lnT_eV *
                                   (-1.65619470e-3 +
                                    lnT_eV * (1.06827520e-4 +
                                              lnT_eV * (-2.63128581e-6)))))))));
}

// k_rxn[14] H- + H+ -> 2H  (GA08 5)
template <RateForm Form>
inline double r14_Hm_Hp_to_2H(const PrimRateContext& c) {
  const double T_K = c.T_K;
  if constexpr (Form == RateForm::Primordial)
    return 2.4e-6 / std::sqrt(T_K) * (1.0 + T_K / 20000.0);
  else
    return 2.4e-6 * std::pow(T_K, -0.5) * (1.0 + T_K / 20000.0);
}

// k_rxn[15] H- + H+ -> H2+ + e  (GA08 16)
inline double r15_Hm_Hp_to_H2p_e(const PrimRateContext& c) {
  return (c.T_K <= 8000.0) ? 6.9e-9 * std::pow(c.T_K, -0.35)
                           : 9.6e-7 * std::pow(c.T_K, -0.90);
}

// k_rxn[16] H2+ + e -> 2H  (GA08 6)
inline double r16_H2p_e_to_2H(const PrimRateContext& c) {
  return (c.T_K < 617.0) ? 1.0e-8 : 1.32e-6 * std::pow(c.T_K, -0.76);
}

// k_rxn[17] H2+ + H- -> H2 + H  (GA08 21)
inline double r17_H2p_Hm_to_H2_H(const PrimRateContext& c) {
  return 1.4e-7 * std::pow(c.T300, -0.5);
}

// k_rxn[18] 3H -> H2 + H  (Forrey+2013, Glover+2015)
inline double r18_three_H(const PrimRateContext& c) {
  return 6.0e-32 * std::pow(c.T_K, -0.25) + 2.0e-31 * std::pow(c.T_K, -0.5);
}

// k_rxn[19] 2H + H2 -> 2H2  (GA08 31)  = r18 / 8
inline double r19_2H_H2_to_2H2(const PrimRateContext& c) {
  return r18_three_H(c) / 8.0;
}

// k_rxn[20] 2H2 -> 2H + H2  (GA08 10, LTE/low-density blend)
inline double r20_2H2_to_2H_H2(const PrimRateContext& c) {
  return lte_blend(c.k21v0, c.k21LTE, c.crit_ratio);
}

// k_rxn[21] H + H -> H + e + H+  (Lenzuni+1991)
template <RateForm Form>
inline double r21_H_H_to_H_e_Hp(const PrimRateContext& c) {
  constexpr double k_B = phys::k_B;
  const double T_K = c.T_K, delE_22 = c.delE_22;
  if constexpr (Form == RateForm::Primordial)
    return 1.2e-17 * std::pow(T_K, 1.2) *
           std::exp(-std::abs(delE_22) / k_B / T_K);
  else
    return 1.2e-17 * std::pow(T_K, 1.2) *
           std::exp(-std::abs(delE_22) / (k_B * T_K));
}

// k_rxn[22] 2H + grain -> H2  (grain surface; set by grain_coef_rates)
inline double r22_2H_grain(const PrimRateContext&) { return 0.0; }

// k_rxn[23] He+ + H2 -> H+ + H + He  (GA08 24)
inline double r23_Hep_H2_to_Hp_H_He(const PrimRateContext& c) {
  return 3.70e-14 * std::exp(-35.0 / c.T_K);
}

// k_rxn[24] H2+ + He -> HeH+ + H  (GS09 TR33)
inline double r24_H2p_He_to_HeHp_H(const PrimRateContext& c) {
  return 3.0e-10 * std::exp(-6717.0 / c.T_K);
}

// k_rxn[25] H2+ + H2 -> H3+ + H  (GS09 TR1)
inline double r25_H2p_H2_to_H3p_H(const PrimRateContext& c) {
  return 2.24e-9 * std::pow(c.T300, 0.042) * std::exp(-c.T_K / 4.66e4);
}

// k_rxn[26] H3+ + H- -> 2H2
inline double r26_H3p_Hm_to_2H2(const PrimRateContext& c) {
  return 2.30e-7 * std::pow(c.T300, -0.50);
}

// k_rxn[27] He+ + H -> H+ + He  (GA08 26)
inline double r27_Hep_H_to_Hp_He(const PrimRateContext& c) {
  return 1.2e-15 * std::pow(c.T300, 0.25);
}

// k_rxn[28] He+ + H- -> H + He  (GA08 28)
inline double r28_Hep_Hm_to_H_He(const PrimRateContext& c) {
  return 2.32e-7 * std::pow(c.T300, -0.52) * std::exp(c.T_K / 22400.0);
}

// k_rxn[29] He+ + H2 -> H2+ + He  (GA08 25)
inline double r29_Hep_H2_to_H2p_He(const PrimRateContext&) { return 7.20e-15; }

// k_rxn[30] HeH+ + H -> H2+ + He  (GS09 TR37)
inline double r30_HeHp_H_to_H2p_He(const PrimRateContext& c) {
  return 1.04e-9 * std::pow(c.T300, 0.13) * std::exp(-c.T_K / 3.31e4);
}

// k_rxn[31] HeH+ + H2 -> H3+ + He  (GS09 TR39)
inline double r31_HeHp_H2_to_H3p_He(const PrimRateContext& c) {
  return 1.53e-9 * std::pow(c.T300, 0.24) * std::exp(-c.T_K / 1.48e4);
}

// k_rxn[32] H2+ + H- -> H3+ + e  (GS09 AD13)
inline double r32_H2p_Hm_to_H3p_e(const PrimRateContext& c) {
  return 2.7e-10 * std::pow(c.T300, -0.485) * std::exp(c.T_K / 3.12e4);
}

// k_rxn[33] H3+ + e -> H2 + H
inline double r33_H3p_e_to_H2_H(const PrimRateContext& c) {
  return 2.34e-8 * std::pow(c.T300, -0.52);
}

// k_rxn[34] H3+ + e -> 3H
inline double r34_H3p_e_to_3H(const PrimRateContext& c) {
  return 4.36e-8 * std::pow(c.T300, -0.52);
}

// k_rxn[35] HeH+ + e -> H + He  (GS09 DR14)
inline double r35_HeHp_e_to_H_He(const PrimRateContext& c) {
  return 3.0e-8 * std::pow(c.T300, -0.47);
}

// k_rxn[36] H2 + He -> 2H + He  (GA08 11, LTE/low-density blend)
inline double r36_H2_He_to_2H_He(const PrimRateContext& c) {
  return lte_blend(c.k37v0, c.k37LTE, c.crit_ratio);
}

// k_rxn[37] H- + H -> 2H + e  (GA08 15, zero below 0.1 eV)
inline double r37_Hm_H_to_2H_e(const PrimRateContext& c) {
  if (c.T_eV < 0.1) return 0.0;
  const double lnT_eV = c.lnT_eV;
  return std::exp(
      -20.372609 +
      lnT_eV *
          (1.13944933 +
           lnT_eV *
               (-1.4210135e-1 +
                lnT_eV *
                    (8.4644554e-3 +
                     lnT_eV *
                         (-1.4327641e-3 +
                          lnT_eV *
                              (2.0122503e-4 +
                               lnT_eV *
                                   (8.6639632e-5 +
                                    lnT_eV *
                                        (-2.5850097e-5 +
                                         lnT_eV *
                                             (2.4555012e-6 +
                                              lnT_eV * (-8.0683825e-8))))))))));
}

// k_rxn[38] H- + H2+ -> 3H  (GA08 22)
inline double r38_Hm_H2p_to_3H(const PrimRateContext& c) {
  return 1.4e-7 * std::pow(c.T300, -0.5);
}

// k_rxn[39] H2 + e -> H- + H  (GA08 23)
inline double r39_H2_e_to_Hm_H(const PrimRateContext& c) { return c.k40v0; }

// k_rxn[40] He + H+ -> He+ + H  (GA08 27)
inline double r40_He_Hp_to_Hep_H(const PrimRateContext& c) {
  return (c.T_K < 1.0e4)
             ? 1.26e-9 * std::pow(c.T_K, -0.75) * std::exp(-127500.0 / c.T_K)
             : 4.0e-37 * std::pow(c.T_K, 4.74);
}

// k_rxn[41] He + H- -> He + H + e  (GA08 29)
inline double r41_He_Hm_to_He_H_e(const PrimRateContext& c) {
  return 4.1e-17 * c.T_K * c.T_K * std::exp(-19870.0 / c.T_K);
}

// k_rxn[42] 2H + He -> H2 + He  (GA08 32)
inline double r42_2H_He_to_H2_He(const PrimRateContext& c) {
  return 6.9e-32 * std::pow(c.T_K, -0.4);
}

// k_rxn[43] H + H2+ -> H3+ + ph.  (GS09 RA6)
inline double r43_H_H2p_to_H3p(const PrimRateContext& c) {
  return 1.5e-17 * std::pow(c.T300, 1.8) * std::exp(20.0 / c.T_K);
}

// k_rxn[44] H2 + H+ -> H3+ + ph.  (Stancil+98)
inline double r44_H2_Hp_to_H3p(const PrimRateContext&) { return 1.0e-20; }

// k_rxn[45] He + H+ -> HeH+ + ph.  (GS09 RA25)
inline double r45_He_Hp_to_HeHp(const PrimRateContext& c) {
  return 8.0e-20 * std::pow(c.T300, -0.24) * std::exp(-c.T_K / 4.0e3);
}

// k_rxn[46] HD + e -> H + D + e  (GA08 111, LTE/low blend with crit_ratio_HD)
inline double r46_HD_e_to_H_D_e(const PrimRateContext& c) {
  return lte_blend(c.k47v0, c.k47LTE, c.crit_ratio_HD);
}

// k_rxn[47] HD + He -> H + D + He  (GA08 110, uses k37 branches)
inline double r47_HD_He_to_H_D_He(const PrimRateContext& c) {
  return lte_blend(c.k37v0, c.k37LTE, c.crit_ratio_HD);
}

// k_rxn[48] HD + H2 -> H + D + H2  (GA08 109, uses k21 branches)
inline double r48_HD_H2_to_H_D_H2(const PrimRateContext& c) {
  return lte_blend(c.k21v0, c.k21LTE, c.crit_ratio_HD);
}

// k_rxn[49] HD + H -> 2H + D  (GA08 108, uses k13 branches)
inline double r49_HD_H_to_2H_D(const PrimRateContext& c) {
  return lte_blend(c.k13v0, c.k13LTE, c.crit_ratio_HD);
}

// k_rxn[50] D+ + e -> D + ph.  = r1
inline double r50_Dp_rec(const PrimRateContext& c) {
  return r1_Hp_rec_caseB(c);
}

// k_rxn[51] D + H+ -> D+ + H  (GA08 34)
inline double r51_D_Hp_to_Dp_H(const PrimRateContext& c) {
  if (c.T_K < 2.0e5)
    return 2.0e-10 * std::pow(c.T_K, 0.402) * std::exp(-37.1 / c.T_K) -
           3.31e-17 * std::pow(c.T_K, 1.48);
  return 3.44e-10 * std::pow(c.T_K, 0.35);
}

// k_rxn[52] D+ + H -> D + H+  (GA08 35)
inline double r52_Dp_H_to_D_Hp(const PrimRateContext& c) {
  return 2.06e-10 * std::pow(c.T_K, 0.396) * std::exp(-33.0 / c.T_K) +
         2.03e-9 * std::pow(c.T_K, -0.332);
}

// k_rxn[53] D + H -> HD + ph.  (GA08 36)
template <RateForm Form>
inline double r53_D_H_to_HD(const PrimRateContext& c) {
  const double T_K = c.T_K, lnT = c.lnT;
  if constexpr (Form == RateForm::Primordial) {
    if (T_K < 200.0)
      return 1.0e-25 *
             (2.80202 - 6.63697 * lnT + 4.75619 * lnT * lnT -
              1.39325 * std::pow(lnT, 3) + 0.178259 * std::pow(lnT, 4) -
              0.00817097 * std::pow(lnT, 5));
    return 1.0e-25 *
           std::exp(507.207 - 370.889 * lnT + 104.854 * lnT * lnT -
                    14.4192 * std::pow(lnT, 3) + 0.971469 * std::pow(lnT, 4) -
                    0.0258076 * std::pow(lnT, 5));
  } else {
    if (T_K > 10.0 && T_K <= 200.0)
      return 1.0e-25 *
             (2.80202 - 6.63697 * lnT + 4.75619 * lnT * lnT -
              1.39325 * std::pow(lnT, 3) + 0.178259 * std::pow(lnT, 4) -
              0.00817097 * std::pow(lnT, 5));
    else if (T_K > 200.0)
      return 1.0e-25 *
             std::exp(507.207 - 370.889 * lnT + 104.854 * lnT * lnT -
                      14.4192 * std::pow(lnT, 3) + 0.971469 * std::pow(lnT, 4) -
                      0.0258076 * std::pow(lnT, 5));
    else
      return 0.0;
  }
}

// k_rxn[54] D + H2 -> H + HD  (GA08 37)
inline double r54_D_H2_to_H_HD(const PrimRateContext& c) {
  if (c.T_K < 2000.0)
    return std::pow(10.0, -56.4737 + 5.88886 * c.lgT + 7.19692 * c.lgT2 +
                              2.25069 * c.lgT3 - 2.16903 * c.lgT4 +
                              0.317887 * c.lgT5);
  return 3.17e-10 * std::exp(-5207.0 / c.T_K);
}

// k_rxn[55] HD+ + H -> H+ + HD  = r9
inline double r55_HDp_H_to_Hp_HD(const PrimRateContext& c) {
  return r9_H2p_H_to_H2_Hp(c);
}

// k_rxn[56] D+ + H2 -> H+ + HD  (GA08 39)
inline double r56_Dp_H2_to_Hp_HD(const PrimRateContext& c) {
  return (0.417 + 0.846 * c.lgT - 0.137 * c.lgT2) * 1.0e-9;
}

// k_rxn[57] HD + H -> H2 + D  (GA08 40)
inline double r57_HD_H_to_H2_D(const PrimRateContext& c) {
  const double T_K = c.T_K;
  return (T_K < 200.0)
             ? 5.25e-11 * std::exp(-4430.0 / T_K)
             : 5.25e-11 * std::exp(-4430.0 / T_K + 173900.0 / (T_K * T_K));
}

// k_rxn[58] HD + H+ -> H2 + D+  (GA08 41)
inline double r58_HD_Hp_to_H2_Dp(const PrimRateContext& c) {
  return 1.1e-9 * std::exp(-488.0 / c.T_K);
}

// k_rxn[59] D + H+ -> HD+ + ph.  (GA08 42)
inline double r59_D_Hp_to_HDp(const PrimRateContext& c) {
  return 3.9e-19 * std::pow(c.T300, 1.8) * std::exp(20.0 / c.T_K);
}

// k_rxn[60] D+ + H -> HD+ + ph.  (GA08 43)
inline double r60_Dp_H_to_HDp(const PrimRateContext& c) {
  return 3.9e-19 * std::pow(c.T300, 1.8) * std::exp(20.0 / c.T_K);
}

// k_rxn[61] HD+ + e -> H + D  (GA08 44)
inline double r61_HDp_e_to_H_D(const PrimRateContext& c) {
  return 7.2e-8 * std::pow(c.T_K, -0.5);
}

// k_rxn[62] D + e -> D- + ph.  = r6
inline double r62_D_e_to_Dm(const PrimRateContext& c) {
  return r6_H_e_to_Hm(c);
}

// k_rxn[63] D+ + D- -> 2D  = r14
template <RateForm Form>
inline double r63_Dp_Dm_to_2D(const PrimRateContext& c) {
  return r14_Hm_Hp_to_2H<Form>(c);
}

// k_rxn[64] H+ + D- -> D + H  = r14
template <RateForm Form>
inline double r64_Hp_Dm_to_D_H(const PrimRateContext& c) {
  return r14_Hm_Hp_to_2H<Form>(c);
}

// k_rxn[65] H- + D -> H + D-  (GA08 53)
inline double r65_Hm_D_to_H_Dm(const PrimRateContext& c) {
  return 6.4e-9 * std::pow(c.T300, 0.41);
}

// k_rxn[66] D- + H -> D + H-  (GA08 52)
inline double r66_Dm_H_to_D_Hm(const PrimRateContext& c) {
  return 6.4e-9 * std::pow(c.T300, 0.41);
}

// k_rxn[67] D- + H -> HD + e  (GA08 55)  = 0.5 * r7
inline double r67_Dm_H_to_HD_e(const PrimRateContext& c) {
  return 0.5 * r7_Hm_H_to_H2_e(c);
}

// k_rxn[68] D + e -> D+ + 2e  = r0
inline double r68_D_e_to_Dp_2e(const PrimRateContext& c) {
  return r0_H_e_to_Hp_2e(c);
}

// k_rxn[69] He+ + D -> D+ + He  (GA08 46)
inline double r69_Hep_D_to_Dp_He(const PrimRateContext& c) {
  return 1.1e-15 * std::pow(c.T300, 0.25);
}

// k_rxn[70] He + D+ -> D + He+  (GA08 47)
inline double r70_He_Dp_to_D_Hep(const PrimRateContext& c) {
  return (c.T_K < 10000.0)
             ? 1.85e-9 * std::pow(c.T_K, -0.75) * std::exp(-127500.0 / c.T_K)
             : 5.9e-37 * std::pow(c.T_K, 4.74);
}

// k_rxn[71] H2+ + D -> HD+ + H  (GA08 48)
inline double r71_H2p_D_to_HDp_H(const PrimRateContext& c) {
  return 1.07e-9 * std::pow(c.T300, 0.062) * std::exp(-c.T_K / 41400.0);
}

// k_rxn[72] HD+ + D -> HD + D+  = r9
inline double r72_HDp_D_to_HD_Dp(const PrimRateContext& c) {
  return r9_H2p_H_to_H2_Hp(c);
}

// k_rxn[73] HD+ + H -> H2+ + D  (GA08 50)
inline double r73_HDp_H_to_H2p_D(const PrimRateContext& c) {
  return 1.0e-9 * std::exp(-154.0 / c.T_K);
}

// k_rxn[74] HD + e -> H + D-  (GA08 57)
inline double r74_HD_e_to_H_Dm(const PrimRateContext& c) { return c.k75v0; }

// k_rxn[75] HD + e -> D + H-  (GA08 58)
inline double r75_HD_e_to_D_Hm(const PrimRateContext& c) {
  return 1.35e-9 * std::pow(c.T_K, -1.27) * std::exp(-43000.0 / c.T_K);
}

// k_rxn[76] H+ + D- -> HD+ + e  (GA08 60)
inline double r76_Hp_Dm_to_HDp_e(const PrimRateContext& c) {
  return 1.1e-9 * std::pow(c.T300, -0.4);
}

// k_rxn[77] D+ + H- -> HD+ + e  (GA08 61)
inline double r77_Dp_Hm_to_HDp_e(const PrimRateContext& c) {
  return 1.1e-9 * std::pow(c.T300, -0.4);
}

// k_rxn[78] D- + e -> D + 2e  = r13
inline double r78_Dm_e_to_D_2e(const PrimRateContext& c) {
  return r13_Hm_e_to_H_2e(c);
}

// k_rxn[79] D- + H -> D + H + e  = r37
inline double r79_Dm_H_to_D_H_e(const PrimRateContext& c) {
  return r37_Hm_H_to_2H_e(c);
}

// k_rxn[80] D- + He -> D + He + e  (GA08 65)
inline double r80_Dm_He_to_D_He_e(const PrimRateContext& c) {
  return 1.5e-17 * c.T_K * c.T_K * std::exp(-19870.0 / c.T_K);
}

// k_rxn[81] D+ + H- -> D + H  = r14
template <RateForm Form>
inline double r81_Dp_Hm_to_D_H(const PrimRateContext& c) {
  return r14_Hm_Hp_to_2H<Form>(c);
}

// k_rxn[82] H2+ + D- -> H2 + D  (GA08 69)
inline double r82_H2p_Dm_to_H2_D(const PrimRateContext& c) {
  return 1.7e-7 * std::pow(c.T300, -0.5);
}

// k_rxn[83] H2+ + D- -> 2H + D  (GA08 70)
inline double r83_H2p_Dm_to_2H_D(const PrimRateContext& c) {
  return 1.7e-7 * std::pow(c.T300, -0.5);
}

// k_rxn[84] HD+ + H- -> HD + H  (GA08 71)
inline double r84_HDp_Hm_to_HD_H(const PrimRateContext& c) {
  return 1.5e-7 * std::pow(c.T300, -0.5);
}

// k_rxn[85] HD+ + H- -> D + 2H  (GA08 72)
inline double r85_HDp_Hm_to_D_2H(const PrimRateContext& c) {
  return 1.5e-7 * std::pow(c.T300, -0.5);
}

// k_rxn[86] HD+ + D- -> HD + D  (GA08 73)
inline double r86_HDp_Dm_to_HD_D(const PrimRateContext& c) {
  return 1.9e-7 * std::pow(c.T300, -0.5);
}

// k_rxn[87] HD+ + D- -> 2D + H  (GA08 74)
inline double r87_HDp_Dm_to_2D_H(const PrimRateContext& c) {
  return 1.9e-7 * std::pow(c.T300, -0.5);
}

// k_rxn[88] He+ + D- -> He + D  (GA08 79)
inline double r88_Hep_Dm_to_He_D(const PrimRateContext& c) {
  return 3.03e-7 * std::pow(c.T300, -0.52) * std::exp(c.T_K / 22400.0);
}

// k_rxn[89] D + H2+ -> H2 + D+  = r9
inline double r89_D_H2p_to_H2_Dp(const PrimRateContext& c) {
  return r9_H2p_H_to_H2_Hp(c);
}

// k_rxn[90] H2+ + D -> HD + H+  (GA08 82)
inline double r90_H2p_D_to_HD_Hp(const PrimRateContext&) { return 1.0e-9; }

// k_rxn[91] HD+ + H -> H2 + D+  (GA08 83)
inline double r91_HDp_H_to_H2_Dp(const PrimRateContext&) { return 1.0e-9; }

// k_rxn[92] H2 + D+ -> H2+ + D  = r10
inline double r92_H2_Dp_to_H2p_D(const PrimRateContext& c) {
  return r10_H2_Hp_to_H2p_H(c);
}

// k_rxn[93] H2 + D+ -> HD+ + H  (GA08 91)
inline double r93_H2_Dp_to_HDp_H(const PrimRateContext& c) {
  const double T_K = c.T_K;
  return (1.04e-9 + 9.52e-9 * (T_K / 10000.0) -
          1.81e-9 * std::pow(T_K / 10000.0, 2)) *
         std::exp(-21000.0 / T_K);
}

// k_rxn[94] HD + H+ -> HD+ + H  = r10
inline double r94_HD_Hp_to_HDp_H(const PrimRateContext& c) {
  return r10_H2_Hp_to_H2p_H(c);
}

// k_rxn[95] HD + H+ -> H2+ + D  (GA08 93)
inline double r95_HD_Hp_to_H2p_D(const PrimRateContext& c) {
  return 1.0e-9 * std::exp(-21600.0 / c.T_K);
}

// k_rxn[96] HD + D+ -> HD+ + D  = r10
inline double r96_HD_Dp_to_HDp_D(const PrimRateContext& c) {
  return r10_H2_Hp_to_H2p_H(c);
}

// k_rxn[97] HD + He+ -> HD+ + He  = r29
inline double r97_HD_Hep_to_HDp_He(const PrimRateContext& c) {
  return r29_Hep_H2_to_H2p_He(c);
}

// k_rxn[98] HD + He+ -> He + H+ + D  (GA08 102)
inline double r98_HD_Hep_to_He_Hp_D(const PrimRateContext& c) {
  return 1.85e-14 * std::exp(-35.0 / c.T_K);
}

// k_rxn[99] HD + He+ -> He + H + D+  (GA08 103)
inline double r99_HD_Hep_to_He_H_Dp(const PrimRateContext& c) {
  return 1.85e-14 * std::exp(-35.0 / c.T_K);
}

// ── Li reactions  k_rxn[100..129]  (1-based reaction nums 101..130) ─────────
// k_rxn[100] e- + Li+ -> Li + γ
inline double r100_e_Lip_to_Li(const PrimRateContext& c) {
  return 1.036e-11 / (std::sqrt(c.T_K / 107.7) *
                      std::pow(1.0 + std::sqrt(c.T_K / 107.7), 0.6612) *
                      std::pow(1.0 + std::sqrt(c.T_K / 1.177e7), 1.3388));
}

// k_rxn[101] H- + Li+ -> H + Li
inline double r101_Hm_Lip_to_H_Li(const PrimRateContext& c) {
  return 6.3e-9 / std::sqrt(c.T_K) * (1.0 + c.T_K / 14000.0);
}

// k_rxn[102] H+ + Li- -> H + Li
inline double r102_Hp_Lim_to_H_Li(const PrimRateContext& c) {
  return 2.3e-6 / std::sqrt(c.T_K);
}

// k_rxn[103] e- + Li -> Li- + γ
inline double r103_e_Li_to_Lim(const PrimRateContext& c) {
  return 6.1e-17 * std::pow(c.T_K, 0.58) * std::exp(-c.T_K / 1.72e4);
}

// k_rxn[104] H+ + Li -> H + Li+
inline double r104_Hp_Li_to_H_Lip(const PrimRateContext& c) {
  return 2.5e-40 * std::pow(c.T_K, 7.9) * std::exp(-c.T_K / 1210.0);
}

// k_rxn[105] H+ + Li -> H + Li+ + γ
inline double r105_Hp_Li_to_H_Lip_ph(const PrimRateContext& c) {
  return 1.7e-13 * std::pow(c.T_K, -0.051) * std::exp(-c.T_K / 282000.0);
}

// k_rxn[106] H- + Li -> e- + LiH
inline double r106_Hm_Li_to_e_LiH(const PrimRateContext&) { return 4.0e-10; }

// k_rxn[107] H + Li- -> e- + LiH
inline double r107_H_Lim_to_e_LiH(const PrimRateContext&) { return 4.0e-10; }

// k_rxn[108] H + LiH+ -> H+ + LiH
inline double r108_H_LiHp_to_Hp_LiH(const PrimRateContext& c) {
  return 1.0e-11 * std::exp(-67900.0 / c.T_K);
}

// k_rxn[109] H+ + LiH -> H + LiH+
inline double r109_Hp_LiH_to_H_LiHp(const PrimRateContext&) { return 1.0e-9; }

// k_rxn[110] H + LiH -> H2 + Li
inline double r110_H_LiH_to_H2_Li(const PrimRateContext& c) {
  return 2.0e-12 * c.T_K * std::exp(-c.T_K / 1200.0);
}

// k_rxn[111] H + Li+ -> LiH+ + γ
inline double r111_H_Lip_to_LiHp(const PrimRateContext& c) {
  return 1.4e-20 * std::pow(c.T_K, -0.9) * std::exp(-c.T_K / 7.0e3);
}

// k_rxn[112] H+ + Li -> LiH+ + γ
inline double r112_Hp_Li_to_LiHp(const PrimRateContext& c) {
  return 5.3e-14 * std::pow(c.T_K, -0.49);
}

// k_rxn[113] H+ + LiH -> H2 + Li+
inline double r113_Hp_LiH_to_H2_Lip(const PrimRateContext&) { return 1.0e-9; }

// k_rxn[114] e- + LiH+ -> H + Li
inline double r114_e_LiHp_to_H_Li(const PrimRateContext& c) {
  return 3.9e-6 * std::pow(c.T_K, -0.70) * std::exp(-c.T_K / 1200.0);
}

// k_rxn[115] H + LiH+ -> H2+ + Li
inline double r115_H_LiHp_to_H2p_Li(const PrimRateContext& c) {
  return 9.0e-10 * std::exp(-66400.0 / c.T_K);
}

// k_rxn[116] H + LiH+ -> H2 + Li+
inline double r116_H_LiHp_to_H2_Lip(const PrimRateContext& c) {
  return 8.7e-10 * std::pow(c.T_K, 0.040) * std::exp(c.T_K / 5.92e8);
}

// k_rxn[117] H + Li -> LiH + γ
inline double r117_H_Li_to_LiH(const PrimRateContext& c) {
  return 4.0e-20 * std::exp(-c.T_K / 4065.0 + std::pow(c.T_K / 13193.0, 3.0));
}

// k_rxn[118] H2+ + Li -> H+ + LiH
inline double r118_H2p_Li_to_Hp_LiH(const PrimRateContext& c) {
  if (c.T_K < 500.0) return 6.3e-10 * std::exp(-2553.0 / c.T_K);
  return 7.2e-14 * std::pow(c.T_K, 1.18) * std::exp(-1470.0 / c.T_K);
}

// k_rxn[119] H+ + LiH -> H2+ + Li
inline double r119_Hp_LiH_to_H2p_Li(const PrimRateContext& c) {
  return 2.9e-10 * std::pow(c.T_K, 0.59) -
         2.6e-10 * std::pow(c.T_K, 0.6) * std::exp(-400.0 / c.T_K);
}

// k_rxn[120] e- + Li++ -> Li+ + γ
inline double r120_e_Lipp_to_Lip(const PrimRateContext& c) {
  return 5.34e-8 * std::pow(c.T300, -1.23) * std::exp(-c.T_K / 9.23e5);
}

// k_rxn[121] e- + Li+++ -> Li++ + γ
inline double r121_e_Lippp_to_Lipp(const PrimRateContext& c) {
  return 4.83e-11 * std::pow(c.T300, -0.621) * std::exp(-c.T_K / 1.67e6);
}

// k_rxn[122] D- + Li+ -> D + Li
inline double r122_Dm_Lip_to_D_Li(const PrimRateContext& c) {
  return 3.71e-7 * std::pow(c.T300, -0.51) * std::exp(c.T_K / 4.41e4);
}

// k_rxn[123] D+ + Li- -> D + Li
inline double r123_Dp_Lim_to_D_Li(const PrimRateContext& c) {
  return 2.28e-7 * std::pow(c.T300, -0.51) * std::exp(c.T_K / 4.41e4);
}

// k_rxn[124] e- + Li -> e- + e- + Li+
inline double r124_e_Li_to_2e_Lip(const PrimRateContext& c) {
  return 3.11e-8 * std::pow(c.T300, 0.163) * std::exp(-6.27e4 / c.T_K);
}

// k_rxn[125] e- + Li+ -> e- + e- + Li++
inline double r125_e_Lip_to_2e_Lipp(const PrimRateContext& c) {
  return 5.67e-12 * std::pow(c.T300, 0.715) * std::exp(-8.77e5 / c.T_K);
}

// k_rxn[126] e- + Li++ -> e- + e- + Li+++
inline double r126_e_Lipp_to_2e_Lippp(const PrimRateContext& c) {
  return 1.70e-12 * std::pow(c.T300, 0.709) * std::exp(-1.42e6 / c.T_K);
}

// k_rxn[127] H + H + Li -> H + LiH
inline double r127_2H_Li_to_H_LiH(const PrimRateContext& c) {
  return 2.5e-29 / c.T_K;
}

// k_rxn[128] H + H2 + Li -> H2 + LiH
inline double r128_H_H2_Li_to_H2_LiH(const PrimRateContext& c) {
  return 4.1e-30 / c.T_K;
}

// k_rxn[129] H2 + Li -> H2 + e- + Li+  (delE from the reaction table)
inline double r129_H2_Li_to_H2_e_Lip(const PrimRateContext& c, double delE) {
  constexpr double k_B = phys::k_B;
  return 9.9e-9 * std::sqrt(c.T_K) * std::exp(delE / k_B / c.T_K);
}

}  // namespace prim
}  // namespace ratelaw

// ---------------------------------------------------------------------------
// compute_HHeD_rates  — bind the GA08 H/He/D forward laws to k_rxn[0..99].
//
// This is the full network's reaction->law map: each slot is assigned the
// value of its ratelaw::prim law (a composed model would assign only the slots
// it carries, reusing the same laws).  All numeric content lives in the laws;
// this assembler only wires slot indices to law names.  The Form template
// parameter selects the per-model expression forms (r3, r14, r21, r53).
// ---------------------------------------------------------------------------
template <RateForm Form, int N_react>
inline void compute_HHeD_rates(std::array<double, 2 * N_react>& k_rxn,
                               const PrimRateContext& c) {
  namespace rl = ratelaw::prim;

  // ── H, He reactions  k_rxn[0..45] ──────────────────────────────────────
  k_rxn[rxn::H_e_to_Hp_2e] = rl::r0_H_e_to_Hp_2e(c);
  k_rxn[rxn::Hp_rec_caseB] = rl::r1_Hp_rec_caseB(c);
  k_rxn[2] = rl::r2_He_e_to_Hep_2e(c);
  k_rxn[rxn::Hep_rec] = rl::r3_Hep_rec<Form>(c);
  k_rxn[4] = rl::r4_Hep_e_to_Hepp_2e(c);
  k_rxn[5] = rl::r5_Hepp_rec_caseB(c);
  k_rxn[rxn::H_e_to_Hm] = rl::r6_H_e_to_Hm(c);
  k_rxn[rxn::Hm_H_to_H2_e] = rl::r7_Hm_H_to_H2_e(c);
  k_rxn[8] = rl::r8_H_Hp_to_H2p(c);
  k_rxn[rxn::H2p_H_to_H2_Hp] = rl::r9_H2p_H_to_H2_Hp(c);
  k_rxn[10] = rl::r10_H2_Hp_to_H2p_H(c);
  k_rxn[11] = rl::r11_H2_e_to_2H_e(c);
  k_rxn[12] = rl::r12_H2_H_to_3H(c);
  k_rxn[rxn::Hm_e_to_H_2e] = rl::r13_Hm_e_to_H_2e(c);
  k_rxn[rxn::Hm_Hp_to_2H] = rl::r14_Hm_Hp_to_2H<Form>(c);
  k_rxn[15] = rl::r15_Hm_Hp_to_H2p_e(c);
  k_rxn[16] = rl::r16_H2p_e_to_2H(c);
  k_rxn[17] = rl::r17_H2p_Hm_to_H2_H(c);
  k_rxn[rxn::three_H] = rl::r18_three_H(c);
  k_rxn[19] = rl::r19_2H_H2_to_2H2(c);
  k_rxn[rxn::H2H2_dis] = rl::r20_2H2_to_2H_H2(c);
  k_rxn[21] = rl::r21_H_H_to_H_e_Hp<Form>(c);
  k_rxn[22] = rl::r22_2H_grain(c);
  k_rxn[23] = rl::r23_Hep_H2_to_Hp_H_He(c);
  k_rxn[24] = rl::r24_H2p_He_to_HeHp_H(c);
  k_rxn[25] = rl::r25_H2p_H2_to_H3p_H(c);
  k_rxn[26] = rl::r26_H3p_Hm_to_2H2(c);
  k_rxn[27] = rl::r27_Hep_H_to_Hp_He(c);
  k_rxn[28] = rl::r28_Hep_Hm_to_H_He(c);
  k_rxn[rxn::Hep_H2_to_H2p_He] = rl::r29_Hep_H2_to_H2p_He(c);
  k_rxn[30] = rl::r30_HeHp_H_to_H2p_He(c);
  k_rxn[31] = rl::r31_HeHp_H2_to_H3p_He(c);
  k_rxn[32] = rl::r32_H2p_Hm_to_H3p_e(c);
  k_rxn[33] = rl::r33_H3p_e_to_H2_H(c);
  k_rxn[34] = rl::r34_H3p_e_to_3H(c);
  k_rxn[35] = rl::r35_HeHp_e_to_H_He(c);
  k_rxn[rxn::H2_He_dis] = rl::r36_H2_He_to_2H_He(c);
  k_rxn[rxn::Hm_H_to_2H_e] = rl::r37_Hm_H_to_2H_e(c);
  k_rxn[38] = rl::r38_Hm_H2p_to_3H(c);
  k_rxn[39] = rl::r39_H2_e_to_Hm_H(c);
  k_rxn[40] = rl::r40_He_Hp_to_Hep_H(c);
  k_rxn[41] = rl::r41_He_Hm_to_He_H_e(c);
  k_rxn[42] = rl::r42_2H_He_to_H2_He(c);
  k_rxn[43] = rl::r43_H_H2p_to_H3p(c);
  k_rxn[44] = rl::r44_H2_Hp_to_H3p(c);
  k_rxn[45] = rl::r45_He_Hp_to_HeHp(c);

  // ── D reactions  k_rxn[46..99] ─────────────────────────────────────────
  k_rxn[46] = rl::r46_HD_e_to_H_D_e(c);
  k_rxn[47] = rl::r47_HD_He_to_H_D_He(c);
  k_rxn[48] = rl::r48_HD_H2_to_H_D_H2(c);
  k_rxn[49] = rl::r49_HD_H_to_2H_D(c);
  k_rxn[50] = rl::r50_Dp_rec(c);
  k_rxn[51] = rl::r51_D_Hp_to_Dp_H(c);
  k_rxn[52] = rl::r52_Dp_H_to_D_Hp(c);
  k_rxn[53] = rl::r53_D_H_to_HD<Form>(c);
  k_rxn[54] = rl::r54_D_H2_to_H_HD(c);
  k_rxn[55] = rl::r55_HDp_H_to_Hp_HD(c);
  k_rxn[56] = rl::r56_Dp_H2_to_Hp_HD(c);
  k_rxn[57] = rl::r57_HD_H_to_H2_D(c);
  k_rxn[58] = rl::r58_HD_Hp_to_H2_Dp(c);
  k_rxn[59] = rl::r59_D_Hp_to_HDp(c);
  k_rxn[60] = rl::r60_Dp_H_to_HDp(c);
  k_rxn[61] = rl::r61_HDp_e_to_H_D(c);
  k_rxn[62] = rl::r62_D_e_to_Dm(c);
  k_rxn[63] = rl::r63_Dp_Dm_to_2D<Form>(c);
  k_rxn[64] = rl::r64_Hp_Dm_to_D_H<Form>(c);
  k_rxn[65] = rl::r65_Hm_D_to_H_Dm(c);
  k_rxn[66] = rl::r66_Dm_H_to_D_Hm(c);
  k_rxn[67] = rl::r67_Dm_H_to_HD_e(c);
  k_rxn[68] = rl::r68_D_e_to_Dp_2e(c);
  k_rxn[69] = rl::r69_Hep_D_to_Dp_He(c);
  k_rxn[70] = rl::r70_He_Dp_to_D_Hep(c);
  k_rxn[71] = rl::r71_H2p_D_to_HDp_H(c);
  k_rxn[72] = rl::r72_HDp_D_to_HD_Dp(c);
  k_rxn[73] = rl::r73_HDp_H_to_H2p_D(c);
  k_rxn[74] = rl::r74_HD_e_to_H_Dm(c);
  k_rxn[75] = rl::r75_HD_e_to_D_Hm(c);
  k_rxn[76] = rl::r76_Hp_Dm_to_HDp_e(c);
  k_rxn[77] = rl::r77_Dp_Hm_to_HDp_e(c);
  k_rxn[78] = rl::r78_Dm_e_to_D_2e(c);
  k_rxn[79] = rl::r79_Dm_H_to_D_H_e(c);
  k_rxn[80] = rl::r80_Dm_He_to_D_He_e(c);
  k_rxn[81] = rl::r81_Dp_Hm_to_D_H<Form>(c);
  k_rxn[82] = rl::r82_H2p_Dm_to_H2_D(c);
  k_rxn[83] = rl::r83_H2p_Dm_to_2H_D(c);
  k_rxn[84] = rl::r84_HDp_Hm_to_HD_H(c);
  k_rxn[85] = rl::r85_HDp_Hm_to_D_2H(c);
  k_rxn[86] = rl::r86_HDp_Dm_to_HD_D(c);
  k_rxn[87] = rl::r87_HDp_Dm_to_2D_H(c);
  k_rxn[88] = rl::r88_Hep_Dm_to_He_D(c);
  k_rxn[89] = rl::r89_D_H2p_to_H2_Dp(c);
  k_rxn[90] = rl::r90_H2p_D_to_HD_Hp(c);
  k_rxn[91] = rl::r91_HDp_H_to_H2_Dp(c);
  k_rxn[92] = rl::r92_H2_Dp_to_H2p_D(c);
  k_rxn[93] = rl::r93_H2_Dp_to_HDp_H(c);
  k_rxn[94] = rl::r94_HD_Hp_to_HDp_H(c);
  k_rxn[95] = rl::r95_HD_Hp_to_H2p_D(c);
  k_rxn[96] = rl::r96_HD_Dp_to_HDp_D(c);
  k_rxn[97] = rl::r97_HD_Hep_to_HDp_He(c);
  k_rxn[98] = rl::r98_HD_Hep_to_He_Hp_D(c);
  k_rxn[99] = rl::r99_HD_Hep_to_He_H_Dp(c);

  // (No zero-forcing here — each caller applies its own set.)
}

// ---------------------------------------------------------------------------
// compute_CR_rates_prim — CR + CR-photo reactions for the primordial network.
//
// Writes the 9 CR channel rates into k_rxn[cr_base + 0..8].  cr_base defaults to
// the full network's CR block head (k_rxn[130..138]); the compact
// Nakauchi2019_Minimal passes its own CR-block head so the same channel set is
// emitted at its compact slots.  The full default reproduces the original
// k_rxn[130..138] writes verbatim (bit-for-bit).
// ---------------------------------------------------------------------------
template <int N_react>
inline void compute_CR_rates_prim(std::array<double, 2 * N_react>& k_rxn,
                                  double zeta,
                                  int cr_base = zero_metal::slot::cr_begin) {
  constexpr double zeta_ref = 1.36e-17;
  double zeta_ratio = zeta / zeta_ref;
  constexpr double omega = model::cr_photo_albedo;

  k_rxn[cr_base + 0] = 5.98e-18 * zeta_ratio;
  k_rxn[cr_base + 1] = 2.86e-19 * zeta_ratio;
  k_rxn[cr_base + 2] = 1.20e-17 * zeta_ratio;
  k_rxn[cr_base + 3] = 1.30e-18 * zeta_ratio;
  k_rxn[cr_base + 4] = 3.90e-21 * zeta_ratio;
  k_rxn[cr_base + 5] = 6.50e-18 * zeta_ratio;
  k_rxn[cr_base + 6] = 1.3e-17 * zeta_ratio * 0.2 / (1.0 - omega);
  k_rxn[cr_base + 7] = 1.3e-17 * zeta_ratio * 250.0 / (1.0 - omega);
  k_rxn[cr_base + 8] = 1.3e-17 * zeta_ratio * 0.2 / (1.0 - omega);
}

}  // namespace arche
