// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// reaction_metal_grain_surface.h
//
// Grain surface reaction rate coefficients for the metal_grain network.
//
// Computes the grain-surface reaction rate coefficients via grain_coef.
//
// Output: grain_coef_rates() fills kgr[0..150]  (151 entries, 0-based).
// These map to tbl.react[820..970] in the main reaction table
// (1-based kgr(1..151) → global reaction ids 991..1141).
//
// NOTE: kgr[129]  (1-based index 130) is absent — never assigned, stays 0.
//
// References:
//   Hocuk et al. 2016, A&A 586, A35          (adsorption / desorption)
//   Esplugues et al. 2016, A&A 591, A91      (surface reactions)
//   Esplugues et al. 2019, arXiv:1904.03420  (surface reactions)
//   Cazaux & Tielens 2010, ApJ 715, 698      (chemisorption H/D)
//   Hasegawa & Herbst 1993, MNRAS 261, 83    (CR desorption)
//
// All units CGS.
// ---------------------------------------------------------------------------
#include <algorithm>
#include <array>
#include <cmath>

#include "core/species_index.h"  // metal_grain::Species, SpId, make_mass_array
#include "cooling/grain.h"       // vol_gr()

namespace arche {

// ---------------------------------------------------------------------------
// grain_coef_rates
//
//   Computes 151 grain surface reaction rate coefficients.
//
//   1-based reaction index → 0-based array index:  kgr(n) → kgr[n-1]
//   Note: index 130 (0-based [129]) is absent/zero.
//
//   Group summary (1-based):
//     [  1.. 21]  adsorption  X(g) → X(p)  or  X(g) → X(c)
//     [ 22.. 42]  thermal desorption  X(p/c) → X(g)
//     [ 43.. 58]  surface reactions, single product stays on grain
//     [ 59.. 74]  surface reactions, single product desorbs
//     [ 75..116]  surface reactions, two products
//     [117..129]  chemisorption reactions (H/D Cazaux&Tielens)
//     [130]       absent (zero)
//     [131..151]  CR + CRUV desorption  (Hasegawa&Herbst 1993)
//
//   Parameters:
//     nH    : hydrogen number density [cm^{-3}]
//     T_K    : gas temperature [K]
//     T_gr_K : grain temperature [K]
//     Z_metal: metallicity relative to solar
//     JH2   : y[H2_p]  = H2 physisorbed abundance (per H)
//     JH2O  : y[H2O_p] = H2O physisorbed abundance (per H)
//     Jtot  : y[H_c] + y[D_c]  = total chemisorbed H+D abundance (per H)
//     zeta   : cosmic-ray ionization rate [s^{-1}]
//     kgr   : output array, length 151
// ---------------------------------------------------------------------------
inline void grain_coef_rates(double nH, double T_K, double T_gr_K,
                             double Z_metal, double JH2, double JH2O,
                             double Jtot, double zeta, double T_cr_eff,
                             std::array<double, 151>& kgr) {
  // ── Physical constants ────────────────────────────────────────────────
  constexpr double k_B = phys::k_B;
  constexpr double pi = phys::pi;
  constexpr double m_p = phys::m_p;
  constexpr double h_b = phys::hbar;  // ℏ [erg·s]
  constexpr double yHe = abundance_ref::yHe;

  // ── Species masses [g] (21 surface species, 0-based), derived from the
  // canonical catalog.  Local order matches metal_grain Sp 68..88:
  //   0:H_p 1:H_c 2:H2_p 3:D_p 4:D_c 5:HD_p 6:O_p 7:O2_p 8:OH_p 9:CO_p
  //   10:CO2_p 11:H2O_p 12:HO2_p 13:H2O2_p 14:HCO_p 15:H2CO_p
  //   16:C_p 17:CH_p 18:CH2_p 19:CH3_p 20:CH4_p
  static constexpr auto mass = make_mass_array<SpeciesSet<
      SpId::H_p, SpId::H_c, SpId::H2_p, SpId::D_p, SpId::D_c, SpId::HD_p,
      SpId::O_p, SpId::O2_p, SpId::OH_p, SpId::CO_p, SpId::CO2_p, SpId::H2O_p,
      SpId::HO2_p, SpId::H2O2_p, SpId::HCO_p, SpId::H2CO_p, SpId::C_p,
      SpId::CH_p, SpId::CH2_p, SpId::CH3_p, SpId::CH4_p>>();

  // ── Derived quantities ────────────────────────────────────────────────
  double rho = nH * (1.0 + 4.0 * yHe) * m_p;
  double stick = 1.0 / (1.0 + 4.0e-2 * std::sqrt(T_K + T_gr_K) + 2.0e-3 * T_K +
                        8.0e-6 * T_K * T_K);
  double nu_0 = 1.0e12;
  double alpha_MRN = 1.0e-21 * Z_metal *
                     (detail::vol_gr(rho, T_gr_K) / detail::vol_gr(rho, 2.725));
  double vel_H = std::sqrt(8.0 * k_B * T_K / (pi * m_p));
  double nd_nsite = 4.0e-5 * nH / 9.0;

  double f_ice = std::min(1.0, JH2O * nH / nd_nsite);
  double f_bare = std::max(0.0, 1.0 - f_ice);
  double f_h2i = std::min(1.0, JH2 * nH / nd_nsite);
  double f_h2b = std::max(0.0, 1.0 - f_h2i);
  double f_chem = Jtot * nH / nd_nsite;

  // Helper: tunneling + activation barrier P_max for given Ea [K] and masses
  auto P_tunnel = [&](double m_red_amu, double Ea_K) -> double {
    double m_red = m_p * m_red_amu;
    return std::max(
        std::exp(-(2.0e-8 / h_b) * std::sqrt(2.0 * m_red * k_B * Ea_K)),
        std::exp(-Ea_K / T_gr_K));
  };

  kgr.fill(0.0);

  // ==========================================================================
  // (1) Adsorption  kgr(1..21)  →  C++ [0..20]
  //     X(g) → X(p):  k = α_MRN · stick · vel_H / √mass
  //     kgr(2): H(g) → H(c)  includes chemisorption probability
  //     kgr(5): D(g) → D(c)  same
  // ==========================================================================

  // kgr[0]: H(g) → H(p)
  kgr[0] = alpha_MRN * stick * vel_H / std::sqrt(mass[0]);

  // kgr[1]: H(g) → H(c)  Esplugues 2016 appendix B
  {
    double m_red = m_p * mass[0];
    double P =
        std::exp(-(2.0e-8 / h_b) * std::sqrt(2.0 * m_red * k_B * 1.0e3)) +
        std::exp(-1.0e3 / T_gr_K);
    kgr[1] =
        alpha_MRN * stick * (vel_H / std::sqrt(mass[0])) * P * (1.0 - f_chem);
  }

  // kgr[2]: H2(g) → H2(p)
  kgr[2] = alpha_MRN * stick * vel_H / std::sqrt(mass[2]);

  // D and HD are not depleted onto grain surfaces: their gas-phase abundance
  // (y ~ 1e-5) far exceeds the number of available surface binding sites per H
  // nucleus (~1e-6), so bulk physisorption of deuterated species is suppressed
  // by site saturation.  The adsorption coefficients below carry no monolayer
  // limit, which would otherwise freeze several monolayers' worth of HD onto
  // cold (T_gr ~ T_CMB) grains and spuriously deplete gas-phase HD.  D/HD are
  // therefore kept in the gas, where HD remains an important coolant.
  // kgr[3]: D(g) → D(p)   (physisorption, disabled)
  kgr[3] = 0.0;
  // kgr[4]: D(g) → D(c)   (chemisorption, disabled)
  kgr[4] = 0.0;
  // kgr[5]: HD(g) → HD(p) (physisorption, disabled)
  kgr[5] = 0.0;

  // 7..16: O, O2, OH, CO, CO2, H2O, HO2, H2O2, HCO, H2CO
  for (int i = 6; i <= 15; ++i)
    kgr[i] = alpha_MRN * stick * vel_H / std::sqrt(mass[i]);

  // 17..21: C, CH, CH2, CH3, CH4
  for (int i = 16; i <= 20; ++i)
    kgr[i] = alpha_MRN * stick * vel_H / std::sqrt(mass[i]);

  // ==========================================================================
  // (2) Thermal desorption  kgr(22..42)  →  C++ [21..41]
  // ==========================================================================

  // kgr[21]: H(p) → H(g)
  kgr[21] = nu_0 * (f_bare * std::exp(-500.0 / T_gr_K) +
                    f_ice * std::exp(-650.0 / T_gr_K));
  // kgr[22]: H(c) → H(g)
  kgr[22] = nu_0 * std::exp(-1.0e4 / T_gr_K);
  // kgr[23]: H2(p) → H2(g)
  kgr[23] = nu_0 * (f_h2b * (f_bare * std::exp(-300.0 / T_gr_K) +
                             f_ice * std::exp(-500.0 / T_gr_K)) +
                    f_h2i * std::exp(-100.0 / T_gr_K));
  // kgr[24]: D(p) → D(g)   [H + 58K: Thi et al.]
  kgr[24] = nu_0 * (f_bare * std::exp(-558.0 / T_gr_K) +
                    f_ice * std::exp(-708.0 / T_gr_K));
  // kgr[25]: D(c) → D(g)
  kgr[25] = nu_0 * std::exp(-1.0e4 / T_gr_K);
  // kgr[26]: HD(p) → HD(g)
  kgr[26] = nu_0 * (f_h2b * (f_bare * std::exp(-358.0 / T_gr_K) +
                             f_ice * std::exp(-558.0 / T_gr_K)) +
                    f_h2i * std::exp(-158.0 / T_gr_K));
  // kgr[27]: O(p) → O(g)
  kgr[27] = nu_0 * (f_bare * std::exp(-1500.0 / T_gr_K) +
                    f_ice * std::exp(-1420.0 / T_gr_K));
  // kgr[28]: O2(p) → O2(g)
  kgr[28] = nu_0 * (f_bare * std::exp(-1250.0 / T_gr_K) +
                    f_ice * std::exp(-1160.0 / T_gr_K));
  // kgr[29]: OH(p) → OH(g)
  kgr[29] = nu_0 * std::exp(-4600.0 / T_gr_K);
  // kgr[30]: CO(p) → CO(g)
  kgr[30] = nu_0 * (f_bare * std::exp(-1200.0 / T_gr_K) +
                    f_ice * std::exp(-1300.0 / T_gr_K));
  // kgr[31]: CO2(p) → CO2(g)
  kgr[31] = nu_0 * (f_bare * std::exp(-3000.0 / T_gr_K) +
                    f_ice * std::exp(-2670.0 / T_gr_K));
  // kgr[32]: H2O(p) → H2O(g)
  kgr[32] = nu_0 * (f_bare * std::exp(-4800.0 / T_gr_K) +
                    f_ice * std::exp(-5700.0 / T_gr_K));
  // kgr[33]: HO2(p) → O(g) + OH(g)
  kgr[33] = nu_0 * std::exp(-4000.0 / T_gr_K);
  // kgr[34]: H2O2(p) → H2O2(g)
  kgr[34] = nu_0 * std::exp(-6000.0 / T_gr_K);
  // kgr[35]: HCO(p) → HCO(g)
  kgr[35] = nu_0 * std::exp(-1600.0 / T_gr_K);
  // kgr[36]: H2CO(p) → H2CO(g)
  kgr[36] = nu_0 * (f_bare * std::exp(-3700.0 / T_gr_K) +
                    f_ice * std::exp(-3250.0 / T_gr_K));
  // kgr[37]: C(p) → C(g)
  kgr[37] = nu_0 * std::exp(-800.0 / T_gr_K);
  // kgr[38]: CH(p) → CH(g)
  kgr[38] = nu_0 * std::exp(-870.0 / T_gr_K);
  // kgr[39]: CH2(p) → CH2(g)
  kgr[39] = nu_0 * std::exp(-945.0 / T_gr_K);
  // kgr[40]: CH3(p) → CH3(g)
  kgr[40] = nu_0 * std::exp(-1017.0 / T_gr_K);
  // kgr[41]: CH4(p) → CH4(g)
  kgr[41] = nu_0 * std::exp(-1090.0 / T_gr_K);

  // ==========================================================================
  // (3) Surface reactions, single product stays on grain  kgr(43..58) →
  // [42..57]
  //   P_bare = f_bare*(exp(-2E1/(3T)) + exp(-2E2/(3T)))
  //   P_ice  = f_ice *(exp(-2E1'/(3T))+ exp(-2E2'/(3T)))
  //   k = (nu_0/N_site)*(P_bare*(1-alpha_b)+P_ice*(1-alpha_i)) [/2 homonuclear]
  // ==========================================================================

  auto Pbr = [&](double E1b, double E2b) -> double {
    return f_bare * (std::exp(-2.0 * E1b / (3.0 * T_gr_K)) +
                     std::exp(-2.0 * E2b / (3.0 * T_gr_K)));
  };
  auto Pic = [&](double E1i, double E2i) -> double {
    return f_ice * (std::exp(-2.0 * E1i / (3.0 * T_gr_K)) +
                    std::exp(-2.0 * E2i / (3.0 * T_gr_K)));
  };

  // kgr[42]: H(p)+H(p) → H2(p)
  kgr[42] =
      (nu_0 / nd_nsite) *
      (Pbr(500., 500.) * (1. - 9.630e-1) + Pic(650., 650.) * (1. - 9.640e-2)) /
      2.0;
  // kgr[43]: H(p)+D(p) → HD(p)
  kgr[43] = (nu_0 / nd_nsite) * (Pbr(500., 500.) * (1. - 9.630e-1) +
                                 Pic(650., 650.) * (1. - 9.640e-2));
  // kgr[44]: H(p)+O(p) → OH(p)
  kgr[44] = (nu_0 / nd_nsite) * (Pbr(500., 1500.) * (1. - 3.875e-1) +
                                 Pic(650., 1420.) * (1. - 3.880e-2));
  // kgr[45]: H(p)+OH(p) → H2O(p)
  kgr[45] = (nu_0 / nd_nsite) * (Pbr(500., 4600.) * (1. - 2.677e-1) +
                                 Pic(650., 4600.) * (1. - 2.680e-2));
  // kgr[46]: H(p)+O2(p) → HO2(p)
  kgr[46] = (nu_0 / nd_nsite) * (Pbr(500., 1250.) * (1. - 1.380e-2) +
                                 Pic(650., 1160.) * (1. - 1.400e-3));
  // kgr[47]: H(p)+CO(p) → HCO(p)   (with tunneling, Ea=2000K)
  {
    double Pb = Pbr(500., 1200.);
    double Pi = Pic(650., 1300.);
    double mu = mass[0] * mass[9] / (mass[0] + mass[9]);
    double Pm = P_tunnel(mu, 2000.0);
    double Pr = Pm / (Pm + Pb + Pi);
    kgr[47] =
        (nu_0 / nd_nsite) * Pr * (Pb * (1. - 6.700e-3) + Pi * (1. - 7.000e-4));
  }
  // kgr[48]: H(p)+HO2(p) → H2O2(p)
  kgr[48] = (nu_0 / nd_nsite) * (Pbr(500., 4000.) * (1. - 4.600e-3) +
                                 Pic(650., 4000.) * (1. - 5.000e-4));
  // kgr[49]: H(p)+HCO(p) → H2CO(p)
  kgr[49] = (nu_0 / nd_nsite) * (Pbr(500., 1600.) * (1. - 6.610e-2) +
                                 Pic(650., 1600.) * (1. - 6.700e-3));
  // kgr[50]: H(p)+C(p) → CH(p)
  kgr[50] = (nu_0 / nd_nsite) * (Pbr(500., 800.) * (1. - 8.212e-1) +
                                 Pic(650., 800.) * (1. - 8.220e-2));
  // kgr[51]: H(p)+CH(p) → CH2(p)
  kgr[51] = (nu_0 / nd_nsite) * (Pbr(500., 870.) * (1. - 7.668e-1) +
                                 Pic(650., 870.) * (1. - 7.670e-2));
  // kgr[52]: H(p)+CH2(p) → CH3(p)
  kgr[52] = (nu_0 / nd_nsite) * (Pbr(500., 945.) * (1. - 6.937e-1) +
                                 Pic(650., 945.) * (1. - 6.940e-2));
  // kgr[53]: H(p)+CH3(p) → CH4(p)
  kgr[53] = (nu_0 / nd_nsite) * (Pbr(500., 1017.) * (1. - 5.886e-1) +
                                 Pic(650., 1017.) * (1. - 5.890e-2));
  // kgr[54]: O(p)+O(p) → O2(p)
  kgr[54] = (nu_0 / nd_nsite) *
            (Pbr(1500., 1500.) * (1. - 6.884e-1) +
             Pic(1420., 1420.) * (1. - 6.890e-2)) /
            2.0;
  // kgr[55]: O(p)+C(p) → CO(p)
  kgr[55] = (nu_0 / nd_nsite) * (Pbr(1500., 800.) * (1. - 8.659e-1) +
                                 Pic(1420., 800.) * (1. - 8.660e-2));
  // kgr[56]: O(p)+CO(p) → CO2(p)  (with tunneling, Ea=650K)
  {
    double Pb = Pbr(1500., 1200.);
    double Pi = Pic(1420., 1300.);
    double mu = mass[6] * mass[9] / (mass[6] + mass[9]);
    // Tunneling factor: empirical form exp(-6.5 + 2.0/T_gr_K), retained for
    // bit-level reproducibility.
    double Pm = std::max(
        std::exp(-(2.0e-8 / h_b) * std::sqrt(2.0 * m_p * mu * k_B * 650.0)),
        std::exp(-6.5 + 2.0 / T_gr_K));
    double Pr = Pm / (Pm + Pb + Pi);
    kgr[56] =
        (nu_0 / nd_nsite) * Pr * (Pb * (1. - 1.403e-1) + Pi * (1. - 1.400e-2));
  }
  // kgr[57]: OH(p)+OH(p) → H2O2(p)
  kgr[57] = (nu_0 / nd_nsite) *
            (Pbr(4600., 4600.) * (1. - 2.000e-4) +
             Pic(4600., 4600.) * (1. - 1.000e-4)) /
            2.0;

  // ==========================================================================
  // (4) Surface reactions, single product desorbs  kgr(59..74) → [58..73]
  //   k = (nu_0/N_site)*(P_bare*alpha_b + P_ice*alpha_i) [/2 homonuclear]
  // ==========================================================================

  // kgr[58]: H(p)+H(p) → H2(g)
  kgr[58] = (nu_0 / nd_nsite) *
            (Pbr(500., 500.) * 9.630e-1 + Pic(650., 650.) * 9.640e-2) / 2.0;
  // kgr[59]: H(p)+D(p) → HD(g)
  kgr[59] = (nu_0 / nd_nsite) *
            (Pbr(500., 500.) * 9.630e-1 + Pic(650., 650.) * 9.640e-2);
  // kgr[60]: H(p)+O(p) → OH(g)
  kgr[60] = (nu_0 / nd_nsite) *
            (Pbr(500., 1500.) * 3.875e-1 + Pic(650., 1420.) * 3.880e-2);
  // kgr[61]: H(p)+OH(p) → H2O(g)
  kgr[61] = (nu_0 / nd_nsite) *
            (Pbr(500., 4600.) * 2.677e-1 + Pic(650., 4600.) * 2.680e-2);
  // kgr[62]: H(p)+O2(p) → HO2(g)
  kgr[62] = (nu_0 / nd_nsite) *
            (Pbr(500., 1250.) * 1.380e-2 + Pic(650., 1160.) * 1.400e-3);
  // kgr[63]: H(p)+CO(p) → HCO(g)  (with tunneling, Ea=2000K)
  {
    double Pb = Pbr(500., 1200.);
    double Pi = Pic(650., 1300.);
    double mu = mass[0] * mass[9] / (mass[0] + mass[9]);
    double Pm = P_tunnel(mu, 2000.0);
    double Pr = Pm / (Pm + Pb + Pi);
    kgr[63] = (nu_0 / nd_nsite) * Pr * (Pb * 6.700e-3 + Pi * 7.000e-4);
  }
  // kgr[64]: H(p)+HO2(p) → H2O2(g)
  kgr[64] = (nu_0 / nd_nsite) *
            (Pbr(500., 4000.) * 4.600e-3 + Pic(650., 4000.) * 5.000e-4);
  // kgr[65]: H(p)+HCO(p) → H2CO(g)
  kgr[65] = (nu_0 / nd_nsite) *
            (Pbr(500., 1600.) * 6.610e-2 + Pic(650., 1600.) * 6.700e-3);
  // kgr[66]: H(p)+C(p) → CH(g)
  kgr[66] = (nu_0 / nd_nsite) *
            (Pbr(500., 800.) * 8.212e-1 + Pic(650., 800.) * 8.220e-2);
  // kgr[67]: H(p)+CH(p) → CH2(g)
  kgr[67] = (nu_0 / nd_nsite) *
            (Pbr(500., 870.) * 7.668e-1 + Pic(650., 870.) * 7.670e-2);
  // kgr[68]: H(p)+CH2(p) → CH3(g)
  kgr[68] = (nu_0 / nd_nsite) *
            (Pbr(500., 945.) * 6.937e-1 + Pic(650., 945.) * 6.940e-2);
  // kgr[69]: H(p)+CH3(p) → CH4(g)
  kgr[69] = (nu_0 / nd_nsite) *
            (Pbr(500., 1017.) * 5.886e-1 + Pic(650., 1017.) * 5.890e-2);
  // kgr[70]: O(p)+O(p) → O2(g)
  kgr[70] = (nu_0 / nd_nsite) *
            (Pbr(1500., 1500.) * 6.884e-1 + Pic(1420., 1420.) * 6.890e-2) / 2.0;
  // kgr[71]: O(p)+C(p) → CO(g)
  kgr[71] = (nu_0 / nd_nsite) *
            (Pbr(1500., 800.) * 8.659e-1 + Pic(1420., 800.) * 8.660e-2);
  // kgr[72]: O(p)+CO(p) → CO2(g)  (with tunneling, Ea=650K)
  // Same empirical tunneling form as kgr[56] (retained for reproducibility).
  {
    double Pb = Pbr(1500., 1200.);
    double Pi = Pic(1420., 1300.);
    double mu = mass[6] * mass[9] / (mass[6] + mass[9]);
    double Pm = std::max(
        std::exp(-(2.0e-8 / h_b) * std::sqrt(2.0 * m_p * mu * k_B * 650.0)),
        std::exp(-6.5 + 2.0 / T_gr_K));
    double Pr = Pm / (Pm + Pb + Pi);
    kgr[72] = (nu_0 / nd_nsite) * Pr * (Pb * 1.403e-1 + Pi * 1.400e-2);
  }
  // kgr[73]: OH(p)+OH(p) → H2O2(g)
  kgr[73] = (nu_0 / nd_nsite) *
            (Pbr(4600., 4600.) * 2.000e-4 + Pic(4600., 4600.) * 1.000e-4) / 2.0;

  // ==========================================================================
  // (5) Surface reactions, two products  kgr(75..116) → [74..115]
  // ==========================================================================

  // Helper: tunneling P_react given (mu_amu, Ea_K, Pb, Pi)
  auto Pr_tun = [&](double mu_amu, double Ea_K, double Pb,
                    double Pi) -> double {
    double Pm = P_tunnel(mu_amu, Ea_K);
    return Pm / (Pm + Pb + Pi);
  };

  // kgr[74]: H(p)+H2O(p) → H2(p)+OH(p)   Ea=9600K
  {
    double Pb = Pbr(500., 4800.);
    double Pi = Pic(650., 5700.);
    double mu = mass[0] * mass[11] / (mass[0] + mass[11]);
    double Pr = Pr_tun(mu, 9600.0, Pb, Pi);
    kgr[74] = (nu_0 / nd_nsite) * Pr * (Pb + Pi);
  }
  // kgr[75]: H(p)+HO2(p) → OH(p)+OH(p)
  kgr[75] = (nu_0 / nd_nsite) * (Pbr(500., 4000.) * (1. - 3.400e-3) +
                                 Pic(650., 4000.) * (1. - 4.000e-4));
  // kgr[76]: H(p)+HO2(p) → OH(g)+OH(g)
  kgr[76] = (nu_0 / nd_nsite) *
            (Pbr(500., 4000.) * 3.400e-3 + Pic(650., 4000.) * 4.000e-4);
  // kgr[77]: H(p)+H2O2(p) → OH(p)+H2O(p)  Ea=1000K
  {
    double Pb = Pbr(500., 6000.);
    double Pi = Pic(650., 6000.);
    double mu = mass[0] * mass[13] / (mass[0] + mass[13]);
    double Pr = Pr_tun(mu, 1000.0, Pb, Pi);
    kgr[77] = (nu_0 / nd_nsite) * Pr * (Pb * 9.716e-1 + Pi * 9.971e-1);
  }
  // kgr[78]: H(p)+H2O2(p) → OH(g)+H2O(p)  Ea=1000K
  {
    double Pb = Pbr(500., 6000.);
    double Pi = Pic(650., 6000.);
    double mu = mass[0] * mass[13] / (mass[0] + mass[13]);
    double Pr = Pr_tun(mu, 1000.0, Pb, Pi);
    kgr[78] = (nu_0 / nd_nsite) * Pr * (Pb * 7.200e-3 + Pi * 8.000e-4);
  }
  // kgr[79]: H(p)+H2O2(p) → OH(g)+H2O(g)  Ea=1000K
  {
    double Pb = Pbr(500., 6000.);
    double Pi = Pic(650., 6000.);
    double mu = mass[0] * mass[13] / (mass[0] + mass[13]);
    double Pr = Pr_tun(mu, 1000.0, Pb, Pi);
    kgr[79] = (nu_0 / nd_nsite) * Pr * (Pb * 2.120e-2 + Pi * 2.100e-3);
  }
  // kgr[80]: H(p)+HCO(p) → H2(p)+CO(p)
  kgr[80] = (nu_0 / nd_nsite) *
            (Pbr(500., 1600.) * 8.260e-2 + Pic(650., 1600.) * 9.082e-1);
  // kgr[81]: H(p)+HCO(p) → H2(g)+CO(p)
  kgr[81] = (nu_0 / nd_nsite) *
            (Pbr(500., 1600.) * 4.827e-1 + Pic(650., 1600.) * 4.820e-2);
  // kgr[82]: H(p)+HCO(p) → H2(g)+CO(g)
  kgr[82] = (nu_0 / nd_nsite) *
            (Pbr(500., 1600.) * 4.347e-1 + Pic(650., 1600.) * 4.360e-2);
  // kgr[83]: H(p)+H2CO(p) → H2(p)+HCO(p)  Ea=2200K
  {
    double Pb = Pbr(500., 3700.);
    double Pi = Pic(650., 3250.);
    double mu = mass[0] * mass[15] / (mass[0] + mass[15]);
    double Pr = Pr_tun(mu, 2200.0, Pb, Pi);
    kgr[83] = (nu_0 / nd_nsite) * Pr * (Pb * 4.948e-1 + Pi * 9.494e-1);
  }
  // kgr[84]: H(p)+H2CO(p) → H2(g)+HCO(p)  Ea=2200K
  {
    double Pb = Pbr(500., 3700.);
    double Pi = Pic(650., 3250.);
    double mu = mass[0] * mass[15] / (mass[0] + mass[15]);
    double Pr = Pr_tun(mu, 2200.0, Pb, Pi);
    kgr[84] = (nu_0 / nd_nsite) * Pr * (Pb * 5.050e-1 + Pi * 5.050e-2);
  }
  // kgr[85]: H(p)+H2CO(p) → H2(g)+HCO(g)  Ea=2200K
  {
    double Pb = Pbr(500., 3700.);
    double Pi = Pic(650., 3250.);
    double mu = mass[0] * mass[15] / (mass[0] + mass[15]);
    double Pr = Pr_tun(mu, 2200.0, Pb, Pi);
    kgr[85] = (nu_0 / nd_nsite) * Pr * (Pb * 2.000e-4 + Pi * 1.000e-4);
  }
  // kgr[86]: H(p)+CO2(p) → OH(p)+CO(p)  Ea=10000K
  {
    double Pb = Pbr(500., 3000.);
    double Pi = Pic(650., 2670.);
    double mu = mass[0] * mass[10] / (mass[0] + mass[10]);
    double Pr = Pr_tun(mu, 10000.0, Pb, Pi);
    kgr[86] = (nu_0 / nd_nsite) * Pr * (Pb + Pi);
  }
  // kgr[87]: H(p)+CH(p) → H2(p)+C(p)
  kgr[87] = (nu_0 / nd_nsite) *
            (Pbr(500., 870.) * 2.226e-1 + Pic(650., 870.) * 9.222e-1);
  // kgr[88]: H(p)+CH(p) → H2(g)+C(p)
  kgr[88] = (nu_0 / nd_nsite) *
            (Pbr(500., 870.) * 3.859e-1 + Pic(650., 870.) * 3.865e-2);
  // kgr[89]: H(p)+CH(p) → H2(g)+C(g)
  kgr[89] = (nu_0 / nd_nsite) *
            (Pbr(500., 870.) * 3.915e-1 + Pic(650., 870.) * 3.915e-2);
  // kgr[90]: H(p)+CH2(p) → H2(p)+CH(p)
  kgr[90] = (nu_0 / nd_nsite) *
            (Pbr(500., 945.) * 9.554e-1 + Pic(650., 945.) * 9.955e-1);
  // kgr[91]: H(p)+CH2(p) → H2(g)+CH(p)
  kgr[91] = (nu_0 / nd_nsite) *
            (Pbr(500., 945.) * 4.452e-2 + Pic(650., 945.) * 4.452e-3);
  // kgr[92]: H(p)+CH2(p) → H2(g)+CH(g)
  kgr[92] = (nu_0 / nd_nsite) *
            (Pbr(500., 945.) * 8.000e-5 + Pic(650., 945.) * 4.800e-5);
  // kgr[93]: H(p)+CH3(p) → H2(p)+CH2(p)
  kgr[93] = (nu_0 / nd_nsite) * (Pbr(500., 1017.) + Pic(650., 1017.));
  // kgr[94]: H(p)+CH4(p) → H2(p)+CH3(p)
  kgr[94] = (nu_0 / nd_nsite) * (Pbr(500., 1090.) + Pic(650., 1090.));
  // kgr[95]: O(p)+OH(p) → H(p)+O2(p)
  kgr[95] = (nu_0 / nd_nsite) *
            (Pbr(1500., 4600.) * 4.547e-1 + Pic(1420., 4600.) * 9.454e-1);
  // kgr[96]: O(p)+OH(p) → H(g)+O2(p)
  kgr[96] = (nu_0 / nd_nsite) *
            (Pbr(1500., 4600.) * 5.264e-1 + Pic(1420., 4600.) * 5.260e-2);
  // kgr[97]: O(p)+OH(p) → H(g)+O2(g)
  kgr[97] = (nu_0 / nd_nsite) *
            (Pbr(1500., 4600.) * 1.890e-2 + Pic(1420., 4600.) * 2.000e-3);
  // kgr[98]: O(p)+HO2(p) → O2(p)+OH(p)
  kgr[98] = (nu_0 / nd_nsite) *
            (Pbr(1500., 4000.) * 8.265e-1 + Pic(1420., 4000.) * 9.826e-1);
  // kgr[99]: O(p)+HO2(p) → O2(g)+OH(p)
  kgr[99] = (nu_0 / nd_nsite) *
            (Pbr(1500., 4000.) * 1.516e-1 + Pic(1420., 4000.) * 1.510e-2);
  // kgr[100]: O(p)+HO2(p) → O2(g)+OH(g)
  kgr[100] = (nu_0 / nd_nsite) *
             (Pbr(1500., 4000.) * 2.190e-2 + Pic(1420., 4000.) * 2.300e-3);
  // kgr[101]: O(p)+HCO(p) → H(p)+CO2(p)
  kgr[101] = (nu_0 / nd_nsite) *
             (Pbr(1500., 1600.) * 1.141e-1 + Pic(1420., 1600.) * 9.114e-1);
  // kgr[102]: O(p)+HCO(p) → H(g)+CO2(p)
  kgr[102] = (nu_0 / nd_nsite) *
             (Pbr(1500., 1600.) * 8.348e-1 + Pic(1420., 1600.) * 8.340e-2);
  // kgr[103]: O(p)+HCO(p) → H(g)+CO2(g)
  kgr[103] = (nu_0 / nd_nsite) *
             (Pbr(1500., 1600.) * 5.110e-2 + Pic(1420., 1600.) * 5.200e-3);
  // kgr[104]: O(p)+H2CO(p) → H2(p)+CO2(p)  Ea=335K
  {
    double Pb = Pbr(1500., 3700.);
    double Pi = Pic(1420., 3250.);
    double mu = mass[6] * mass[15] / (mass[6] + mass[15]);
    double Pr = Pr_tun(mu, 335.0, Pb, Pi);
    kgr[104] = (nu_0 / nd_nsite) * Pr * (Pb * 7.310e-2 + Pi * 9.073e-1);
  }
  // kgr[105]: O(p)+H2CO(p) → H2(g)+CO2(p)  Ea=335K
  {
    double Pb = Pbr(1500., 3700.);
    double Pi = Pic(1420., 3250.);
    double mu = mass[6] * mass[15] / (mass[6] + mass[15]);
    double Pr = Pr_tun(mu, 335.0, Pb, Pi);
    kgr[105] = (nu_0 / nd_nsite) * Pr * (Pb * 8.901e-1 + Pi * 8.900e-2);
  }
  // kgr[106]: O(p)+H2CO(p) → H2(g)+CO2(g)  Ea=335K
  {
    double Pb = Pbr(1500., 3700.);
    double Pi = Pic(1420., 3250.);
    double mu = mass[6] * mass[15] / (mass[6] + mass[15]);
    double Pr = Pr_tun(mu, 335.0, Pb, Pi);
    kgr[106] = (nu_0 / nd_nsite) * Pr * (Pb * 3.680e-2 + Pi * 3.700e-3);
  }
  // kgr[107]: H2(p)+OH(p) → H(p)+H2O(p)  Ea=2100K
  {
    double Pb = Pbr(300., 4600.);
    double Pi = Pic(500., 4600.);
    double mu = mass[2] * mass[8] / (mass[2] + mass[8]);
    double Pr = Pr_tun(mu, 2100.0, Pb, Pi);
    kgr[107] = (nu_0 / nd_nsite) * Pr * (Pb * 5.948e-1 + Pi * 9.594e-1);
  }
  // kgr[108]: H2(p)+OH(p) → H(g)+H2O(p)  Ea=2100K
  {
    double Pb = Pbr(300., 4600.);
    double Pi = Pic(500., 4600.);
    double mu = mass[2] * mass[8] / (mass[2] + mass[8]);
    double Pr = Pr_tun(mu, 2100.0, Pb, Pi);
    kgr[108] = (nu_0 / nd_nsite) * Pr * (Pb * 4.052e-1 + Pi * 4.060e-2);
  }
  // kgr[109]: OH(p)+CO(p) → H(p)+CO2(p)  Ea=400K
  {
    double Pb = Pbr(4600., 1200.);
    double Pi = Pic(4600., 1300.);
    double mu = mass[8] * mass[9] / (mass[8] + mass[9]);
    double Pr = Pr_tun(mu, 400.0, Pb, Pi);
    kgr[109] = (nu_0 / nd_nsite) * Pr * (Pb * 4.205e-1 + Pi * 9.420e-1);
  }
  // kgr[110]: OH(p)+CO(p) → H(g)+CO2(p)  Ea=400K
  {
    double Pb = Pbr(4600., 1200.);
    double Pi = Pic(4600., 1300.);
    double mu = mass[8] * mass[9] / (mass[8] + mass[9]);
    double Pr = Pr_tun(mu, 400.0, Pb, Pi);
    kgr[110] = (nu_0 / nd_nsite) * Pr * (Pb * 5.794e-1 + Pi * 5.794e-2);
  }
  // kgr[111]: OH(p)+CO(p) → H(g)+CO2(g)  Ea=400K
  {
    double Pb = Pbr(4600., 1200.);
    double Pi = Pic(4600., 1300.);
    double mu = mass[8] * mass[9] / (mass[8] + mass[9]);
    double Pr = Pr_tun(mu, 400.0, Pb, Pi);
    kgr[111] = (nu_0 / nd_nsite) * Pr * (Pb * 1.000e-4 + Pi * 6.000e-5);
  }
  // kgr[112]: OH(p)+HCO(p) → H2(p)+CO2(p)
  kgr[112] = (nu_0 / nd_nsite) *
             (Pbr(4600., 1600.) * 8.063e-2 + Pic(4600., 1600.) * 9.080e-1);
  // kgr[113]: OH(p)+HCO(p) → H2(g)+CO2(p)
  kgr[113] = (nu_0 / nd_nsite) *
             (Pbr(4600., 1600.) * 8.936e-1 + Pic(4600., 1600.) * 8.936e-2);
  // kgr[114]: OH(p)+HCO(p) → H2(g)+CO2(g)
  kgr[114] = (nu_0 / nd_nsite) *
             (Pbr(4600., 1600.) * 2.577e-2 + Pic(4600., 1600.) * 2.640e-3);
  // kgr[115]: H2(p)+HO2(p) → H(p)+H2O2(p)  Ea=5000K
  {
    double Pb = Pbr(300., 4000.);
    double Pi = Pic(500., 4000.);
    double mu = mass[2] * mass[12] / (mass[2] + mass[12]);
    double Pr = Pr_tun(mu, 5000.0, Pb, Pi);
    kgr[115] = (nu_0 / nd_nsite) * Pr * (Pb + Pi);
  }

  // ==========================================================================
  // (6) Chemisorption reactions  kgr(117..129) → [116..128]
  // ==========================================================================

  // Cazaux & Tielens 2010 helper for physisorbed → chemisorbed transition
  // Parameters: E_Hc [K], E_Hp [K], E_s [K], mass [amu]
  // Returns: k_transition factor (not yet multiplied by alpha_MRN)
  auto CT_rate = [&](double E_Hc, double E_Hp, double E_s,
                     double m_amu) -> double {
    double m_red = m_p * m_amu;
    double P =
        std::exp(-(6.0e-8 / h_b) * std::sqrt(2.0 * m_red * k_B * (E_Hp - E_s)));
    return nu_0 * (8.0 * std::sqrt(pi * T_gr_K) * P * std::sqrt(E_Hc - E_Hp) /
                       (E_Hc - E_s) +
                   4.0 * std::sqrt((E_Hp - E_s) / (E_Hc - E_s)) *
                       std::exp(-(E_Hp - E_s) / T_gr_K));
  };

  // kgr[116]: H(p) → H(c)
  kgr[116] = CT_rate(1.0e4, 500.0, 200.0, mass[0]) * (1.0 - f_chem);

  // kgr[117]: H(g) + H(c) → H2(g)
  {
    double m_red = m_p * mass[0];
    double P =
        std::exp(-(2.0e-8 / h_b) * std::sqrt(2.0 * m_red * k_B * 1.0e3)) +
        std::exp(-1.0e3 / T_gr_K);
    kgr[117] = alpha_MRN * stick * (vel_H / std::sqrt(mass[0])) * P / nd_nsite;
  }

  // kgr[118]: H(p) + H(c) → H2(g)
  kgr[118] = CT_rate(1.0e4, 500.0, 200.0, mass[0]) / nd_nsite;

  // kgr[119]: D(p) → D(c)
  kgr[119] = CT_rate(1.0e4, 558.0, 200.0, mass[3]) * (1.0 - f_chem);

  // kgr[120]: H(g) + D(c) → HD(g)  [same rate as 118]
  kgr[120] = kgr[117];

  // kgr[121]: D(g) + H(c) → HD(g)
  {
    double m_red = m_p * mass[3];
    double P =
        std::exp(-(2.0e-8 / h_b) * std::sqrt(2.0 * m_red * k_B * 1.0e3)) +
        std::exp(-1.0e3 / T_gr_K);
    kgr[121] = alpha_MRN * stick * (vel_H / std::sqrt(mass[3])) * P / nd_nsite;
  }

  // kgr[122]: H(p) + D(c) → HD(g)  [same rate as 119]
  kgr[122] = kgr[118];

  // kgr[123]: D(p) + H(c) → HD(g)
  kgr[123] = CT_rate(1.0e4, 558.0, 200.0, mass[3]) / nd_nsite;

  // kgr[124]: H(g) + H(p) → H2(g)  = kgr[0] / nd_nsite
  kgr[124] = kgr[0] / nd_nsite;

  // kgr[125]: H(g) + D(p) → HD(g)  = kgr[0] / nd_nsite
  kgr[125] = kgr[0] / nd_nsite;

  // kgr[126]: D(g) + H(p) → HD(g)  = kgr[3] / nd_nsite
  kgr[126] = kgr[3] / nd_nsite;

  // kgr[127]: H(g) + H(g) → H2(g)  (simple formula, Cazaux&Tielens)
  {
    double m_red = m_p * mass[0];
    double Ptun =
        std::exp(-(6.0e-8 / h_b) * std::sqrt(2.0 * m_red * k_B * 300.0));
    // Reuse CT_rate-like expression for effective chemisorption probability
    double P_eff = 8.0 * std::sqrt(pi * T_gr_K) * Ptun *
                       std::sqrt(1.0e4 - 500.0) / (1.0e4 - 200.0) +
                   4.0 * std::sqrt((500.0 - 200.0) / (1.0e4 - 200.0)) *
                       std::exp(-(500.0 - 200.0) / T_gr_K);
    kgr[127] = alpha_MRN * stick * vel_H / std::sqrt(mass[0]) / 2.0 /
               (1.0 + std::exp(-500.0 / T_gr_K) / P_eff);
  }

  // kgr[128]: H(g) + D(g) → HD(g)
  {
    double m_red = m_p * mass[3];
    double Ptun =
        std::exp(-(6.0e-8 / h_b) * std::sqrt(2.0 * m_red * k_B * 358.0));
    double P_eff = 8.0 * std::sqrt(pi * T_gr_K) * Ptun *
                       std::sqrt(1.0e4 - 558.0) / (1.0e4 - 200.0) +
                   4.0 * std::sqrt((558.0 - 200.0) / (1.0e4 - 200.0)) *
                       std::exp(-(558.0 - 200.0) / T_gr_K);
    kgr[128] = alpha_MRN * stick * vel_H / std::sqrt(mass[3]) /
               (1.0 + std::exp(-558.0 / T_gr_K) / P_eff);
  }

  // kgr[129] (1-based index 130) intentionally absent — stays 0.

  // ==========================================================================
  // (7) CR + CRUV desorption  kgr(131..151) → [130..150]
  //   Pattern:  k_base = nu_0 * f(T_eff=T_cr_eff)
  //   k_CR = (3.16e-19 * k_base + 2.07e-16) * (zeta / 1.36e-17)
  //   [Hasegawa & Herbst 1993 eq. 9 + CRUV term]
  // ==========================================================================

  auto CR = [&](double k_base) -> double {
    return (3.16e-19 * k_base + 2.07e-16) * (zeta / 1.36e-17);
  };

  // kgr[130]: H(p) → H(g)
  kgr[130] = CR(nu_0 * (f_bare * std::exp(-500.0 / T_cr_eff) +
                        f_ice * std::exp(-650.0 / T_cr_eff)));
  // kgr[131]: H(c) → H(g)
  kgr[131] = CR(nu_0 * std::exp(-1.0e4 / T_cr_eff));
  // kgr[132]: H2(p) → H2(g)
  kgr[132] = CR(nu_0 * (f_h2b * (f_bare * std::exp(-300.0 / T_cr_eff) +
                                 f_ice * std::exp(-500.0 / T_cr_eff)) +
                        f_h2i * std::exp(-100.0 / T_cr_eff)));
  // kgr[133]: D(p) → D(g)
  kgr[133] = CR(nu_0 * (f_bare * std::exp(-558.0 / T_cr_eff) +
                        f_ice * std::exp(-708.0 / T_cr_eff)));
  // kgr[134]: D(c) → D(g)
  kgr[134] = CR(nu_0 * std::exp(-1.0e4 / T_cr_eff));
  // kgr[135]: HD(p) → HD(g)
  kgr[135] = CR(nu_0 * (f_h2b * (f_bare * std::exp(-358.0 / T_cr_eff) +
                                 f_ice * std::exp(-558.0 / T_cr_eff)) +
                        f_h2i * std::exp(-158.0 / T_cr_eff)));
  // kgr[136]: O(p) → O(g)
  kgr[136] = CR(nu_0 * (f_bare * std::exp(-1500.0 / T_cr_eff) +
                        f_ice * std::exp(-1420.0 / T_cr_eff)));
  // kgr[137]: O2(p) → O2(g)
  kgr[137] = CR(nu_0 * (f_bare * std::exp(-1250.0 / T_cr_eff) +
                        f_ice * std::exp(-1160.0 / T_cr_eff)));
  // kgr[138]: OH(p) → OH(g)
  kgr[138] = CR(nu_0 * std::exp(-4600.0 / T_cr_eff));
  // kgr[139]: CO(p) → CO(g)
  kgr[139] = CR(nu_0 * (f_bare * std::exp(-1200.0 / T_cr_eff) +
                        f_ice * std::exp(-1300.0 / T_cr_eff)));
  // kgr[140]: CO2(p) → CO2(g)
  kgr[140] = CR(nu_0 * (f_bare * std::exp(-3000.0 / T_cr_eff) +
                        f_ice * std::exp(-2670.0 / T_cr_eff)));
  // kgr[141]: H2O(p) → H2O(g)
  kgr[141] = CR(nu_0 * (f_bare * std::exp(-4800.0 / T_cr_eff) +
                        f_ice * std::exp(-5700.0 / T_cr_eff)));
  // kgr[142]: HO2(p) → O(g)+OH(g)
  kgr[142] = CR(nu_0 * std::exp(-4000.0 / T_cr_eff));
  // kgr[143]: H2O2(p) → H2O2(g)
  kgr[143] = CR(nu_0 * std::exp(-6000.0 / T_cr_eff));
  // kgr[144]: HCO(p) → HCO(g)
  kgr[144] = CR(nu_0 * std::exp(-1600.0 / T_cr_eff));
  // kgr[145]: H2CO(p) → H2CO(g)
  kgr[145] = CR(nu_0 * (f_bare * std::exp(-3700.0 / T_cr_eff) +
                        f_ice * std::exp(-3250.0 / T_cr_eff)));
  // kgr[146]: C(p) → C(g)
  kgr[146] = CR(nu_0 * std::exp(-800.0 / T_cr_eff));
  // kgr[147]: CH(p) → CH(g)
  kgr[147] = CR(nu_0 * std::exp(-870.0 / T_cr_eff));
  // kgr[148]: CH2(p) → CH2(g)
  kgr[148] = CR(nu_0 * std::exp(-945.0 / T_cr_eff));
  // kgr[149]: CH3(p) → CH3(g)
  kgr[149] = CR(nu_0 * std::exp(-1017.0 / T_cr_eff));
  // kgr[150]: CH4(p) → CH4(g)
  kgr[150] = CR(nu_0 * std::exp(-1090.0 / T_cr_eff));
}

}  // namespace arche
