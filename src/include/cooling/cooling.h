// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// cooling.h — radiative cooling and heating rates (zero-metal network)
//
// Provides H2 and HD line cooling (H2_cooling, fesc_H2, c_H2, HD_cooling,
// fesc_HD) and the line/continuum cooling drivers (line_cool, cnt_cool).
//
// All Lambda values are specific cooling rates [erg g⁻¹ s⁻¹]

#include <array>
#include <cmath>

#include "kinetics/rates.h"  // detail::eval_opacity (k_prm)
#include "thermo/h2_thermodynamics.h"

namespace arche {

// ─────────────────────────────────────────────────────────────────────────────
// c_H2 — compatibility entry point for the H2 EOS factor.
//
// The implementation belongs to thermo so the EOS layer does not depend on
// cooling.  Keep this wrapper while existing cooling and solver callers use
// arche::c_H2().
// ─────────────────────────────────────────────────────────────────────────────
inline double c_H2(double T_K) { return thermo::c_H2(T_K); }

// ─────────────────────────────────────────────────────────────────────────────
// fesc_H2 — H2 line escape fraction (fitting formula)
// ─────────────────────────────────────────────────────────────────────────────
inline double fesc_H2(double Nc_H2, double T_K) {
  double a0, a1, a2, a3, a4, a5;
  if (T_K < 1.0e3) {
    a0 = 0.978382;
    a1 = 3.99572e-4;
    a2 = -1.73108e-6;
    a3 = 1.15363e-9;
    a4 = 8.24607e-13;
    a5 = -7.65975e-16;
  } else if (T_K < 4.0e3) {
    a0 = 0.827472;
    a1 = 1.07697e-4;
    a2 = -8.25123e-8;
    a3 = 1.92812e-12;
    a4 = 5.7192e-15;
    a5 = -7.869e-19;
  } else {
    a0 = 1.11569;
    a1 = -3.29302e-4;
    a2 = 1.01846e-7;
    a3 = -1.46666e-11;
    a4 = 1.00764e-15;
    a5 = -2.68873e-20;
  }
  double alpha =
      a0 + T_K * (a1 + T_K * (a2 + T_K * (a3 + T_K * (a4 + T_K * a5))));

  constexpr double b0 = 24.0561, b1 = 1.10043e-3, b2 = -2.87224e-7,
                   b3 = 6.11525e-11, b4 = -6.55034e-15, b5 = 2.54997e-19;
  double beta =
      b0 + T_K * (b1 + T_K * (b2 + T_K * (b3 + T_K * (b4 + T_K * b5))));

  double lnf = -alpha * std::log10(1.0 + Nc_H2 / std::pow(10.0, beta));
  return std::pow(10.0, lnf);
}

// ─────────────────────────────────────────────────────────────────────────────
// h2_cooling — H2 line cooling rate [erg g⁻¹ s⁻¹]
//   Low-density rates: Glover & Abel (2008)
//   LTE limit:        Glover (2015) Eq. 30
//
//   y_a  = y[Sp::H]  (H atoms)
//   y_m  = y[Sp::H2]  (H2)
//   y_e  = y[Sp::e]  (e−)
//   y_Hp = y[Sp::Hp]  (H+)
//   y_He = y[Sp::He]  (He)
// ─────────────────────────────────────────────────────────────────────────────
inline double h2_cooling(double nH, double T_K, double rho, double y_a,
                         double y_m, double y_e, double y_Hp, double y_He,
                         double Nc_H2, double tau_cnt) {
  const double T3 = T_K / 1.0e3;
  const double lgT3 = std::log10(T3);

  // Fixed reference log10 values for clamping
  constexpr double lgT1 = -2.0;          // log10(10/1000)
  constexpr double lgT2 = -1.0;          // log10(100/1000)
  constexpr double lgT4 = 1.0;           // log10(10000/1000)
  constexpr double lgT6 = 0.7781512504;  // log10(6000/1000)

  // Horner-form polynomial evaluation
  auto poly5 = [](double x, double a0, double a1, double a2, double a3,
                  double a4, double a5) -> double {
    return a0 + x * (a1 + x * (a2 + x * (a3 + x * (a4 + x * a5))));
  };
  auto poly8 = [](double x, double a0, double a1, double a2, double a3,
                  double a4, double a5, double a6, double a7,
                  double a8) -> double {
    return a0 +
           x * (a1 +
                x * (a2 +
                     x * (a3 +
                          x * (a4 + x * (a5 + x * (a6 + x * (a7 + x * a8)))))));
  };

  // ── H2–H collision (Glover & Abel 2008) ──────────────────────────────────
  double L0_H2_H;
  if (T_K < 1.0e2) {
    double t = (T_K < 1.0e1) ? lgT1 : lgT3;  // T<10: clamp at T=10K
    L0_H2_H = std::pow(10.0, poly5(t, -16.818342, 37.383713, 58.145166,
                                   48.656103, 20.159831, 3.8479610));
  } else if (T_K <= 1.0e3) {
    L0_H2_H = std::pow(10.0, poly5(lgT3, -24.311209, 3.5692468, -11.332860,
                                   -27.850082, -21.328264, -4.2519023));
  } else if (T_K <= 6.0e3) {
    L0_H2_H = std::pow(10.0, poly5(lgT3, -24.311209, 4.6450521, -3.7209846,
                                   5.9369081, -5.5108047, 1.5538288));
  } else {  // T > 6000: clamp at T=6000K
    L0_H2_H = std::pow(10.0, poly5(lgT6, -24.311209, 4.6450521, -3.7209846,
                                   5.9369081, -5.5108047, 1.5538288));
  }

  // ── H2–H2 collision (Glover & Abel 2008) ─────────────────────────────────
  double L0_H2_H2;
  {
    double t = (T_K < 1.0e2) ? lgT2 : (T_K <= 6.0e3) ? lgT3 : lgT6;
    L0_H2_H2 = std::pow(10.0, poly5(t, -23.962112, 2.0943374, -0.77151436,
                                    0.43693353, -0.14913216, -0.033638326));
  }

  // ── H2–He collision (Glover & Abel 2008) ─────────────────────────────────
  double L0_H2_He;
  {
    double t = (T_K < 1.0e1) ? lgT1 : (T_K <= 6.0e3) ? lgT3 : lgT6;
    L0_H2_He = std::pow(10.0, poly5(t, -23.689237, 2.1892372, -0.81520438,
                                    0.29036281, -0.16596184, 0.19191375));
  }

  // ── H2–H+ collision (Glover 2015) ────────────────────────────────────────
  double L0_H2_Hp;
  {
    double t = (T_K < 1.0e1) ? lgT1 : (T_K <= 1.0e4) ? lgT3 : lgT4;
    L0_H2_Hp = std::pow(10.0, poly5(t, -22.089523, 1.5714711, 0.015391166,
                                    -0.23619985, -0.51002221, 0.32168730));
  }

  // ── H2–e collision (Glover & Abel 2008, 9-coefficient polynomial) ────────
  double L0_H2_e;
  if (T_K < 1.5e2) {
    double t = (T_K < 1.0e1) ? lgT1 : lgT3;  // T<10: clamp at T=10K
    L0_H2_e =
        std::pow(10.0, poly8(t, -34.286155, -48.537163, -77.121176, -51.352459,
                             -15.169160, -0.98120322, 0.0, 0.0, 0.0));
  } else if (T_K <= 5.0e2) {
    L0_H2_e = std::pow(
        10.0, poly8(lgT3, -21.928796, 16.815730, 96.743155, 343.19180,
                    734.71651, 983.67576, 801.81247, 364.14446, 70.609154));
  } else if (T_K <= 1.0e4) {
    L0_H2_e = std::pow(
        10.0, poly8(lgT3, -22.921189, 1.6802758, 0.93310622, 4.0406627,
                    -4.7274036, -8.8077017, 8.9167183, 6.4380698, -6.3701156));
  } else {  // T > 10000: clamp at T=10000K
    L0_H2_e = std::pow(
        10.0, poly8(lgT4, -22.921189, 1.6802758, 0.93310622, 4.0406627,
                    -4.7274036, -8.8077017, 8.9167183, 6.4380698, -6.3701156));
  }

  // Combined low-density rate [erg/s per H2 molecule]
  double L0_H2 = (y_a * L0_H2_H + y_m * L0_H2_H2 + y_He * L0_H2_He +
                  y_Hp * L0_H2_Hp + y_e * L0_H2_e) *
                 nH;

  // ── LTE limit (Glover 2015 Eq. 30) ───────────────────────────────────────
  double L_LTE;
  {
    double t = (T_K <= 1.0e2) ? lgT2 : (T_K <= 1.0e4) ? lgT3 : lgT4;
    L_LTE = std::pow(
        10.0, poly8(t, -20.584225, 5.0194035, -1.5738805, -4.7155769, 2.4714161,
                    5.4710750, -3.9467356, -2.2148338, 1.8161874));
  }

  // Critical-density interpolation
  double frac_cr = L_LTE / L0_H2;
  double Lambda_H2 = L_LTE / (1.0 + frac_cr);

  // Line escape fraction
  double fesc = fesc_H2(Nc_H2, T_K);

  // Specific cooling rate [erg g⁻¹ s⁻¹]
  return Lambda_H2 * y_m * nH * fesc * std::exp(-tau_cnt) / rho;
}

// ─────────────────────────────────────────────────────────────────────────────
// fesc_HD — HD line escape fraction (fitting formula)
// ─────────────────────────────────────────────────────────────────────────────
inline double fesc_HD(double Nc_HD, double T_K) {
  double a0, a1, a2, a3, a4, a5;
  if (T_K < 1.0e3) {
    a0 = 1.02339;
    a1 = -9.04450e-4;
    a2 = 3.26141e-6;
    a3 = -5.73674e-9;
    a4 = 4.98247e-12;
    a5 = -1.68078e-15;
  } else {
    a0 = 8.94334e-1;
    a1 = 6.21583e-5;
    a2 = -1.47443e-9;
    a3 = -9.07612e-12;
    a4 = 3.31028e-15;
    a5 = -3.79804e-19;
  }
  double alpha =
      a0 + T_K * (a1 + T_K * (a2 + T_K * (a3 + T_K * (a4 + T_K * a5))));

  // Ncr: polynomial in log10(T_K)
  double lgT = std::log10(T_K);
  constexpr double b0 = 18.8924, b1 = -1.48020, b2 = 0.85905, b3 = 0.43840,
                   b4 = -0.30569, b5 = 0.046412;
  double lgNc =
      b0 + lgT * (b1 + lgT * (b2 + lgT * (b3 + lgT * (b4 + lgT * b5))));
  double Ncr = std::pow(10.0, lgNc);

  return 1.0 / std::pow(1.0 + Nc_HD / Ncr, alpha);
}

// ─────────────────────────────────────────────────────────────────────────────
// hd_cooling — HD line cooling rate [erg g⁻¹ s⁻¹]
//   Lipovka (2005) density-dependent fitting formula
// ─────────────────────────────────────────────────────────────────────────────
inline double hd_cooling(double nH, double T_K, double rho, double y_HD,
                         double Nc_HD, double tau_cnt) {
  double lgT = std::log10(T_K);
  double L0_HD;

  if (nH < 1.0) {
    // Low-density limit: L0_HD ∝ nH
    double lgL = -42.57688 + 21.93385 * lgT - 10.19097 * (lgT * lgT) +
                 2.19906 * (lgT * lgT * lgT) -
                 0.17334 * (lgT * lgT * lgT * lgT);
    L0_HD = std::pow(10.0, lgL) * nH;
  } else {
    double lgn = (nH < 1.0e8) ? std::log10(nH) : std::log10(1.0e8);

    // Lipovka (2005) bivariate polynomial in log10(T) and log10(n)
    double lgL =
        -42.57688 + 0.92433 * lgn + 0.54962 * (lgn * lgn) -
        0.07676 * (lgn * lgn * lgn) + 0.00275 * (lgn * lgn * lgn * lgn) +
        (21.93385 + 0.77952 * lgn - 1.06447 * (lgn * lgn) +
         0.11864 * (lgn * lgn * lgn) - 0.00366 * (lgn * lgn * lgn * lgn)) *
            lgT +
        (-10.19097 - 0.54263 * lgn + 0.62343 * (lgn * lgn) -
         0.07366 * (lgn * lgn * lgn) + 0.002514 * (lgn * lgn * lgn * lgn)) *
            (lgT * lgT) +
        (2.19906 + 0.11711 * lgn - 0.13768 * (lgn * lgn) +
         0.01759 * (lgn * lgn * lgn) - 0.000666317 * (lgn * lgn * lgn * lgn)) *
            (lgT * lgT * lgT) +
        (-0.17334 - 0.00835 * lgn + 0.0106 * (lgn * lgn) -
         0.001482 * (lgn * lgn * lgn) + 0.000061926 * (lgn * lgn * lgn * lgn)) *
            (lgT * lgT * lgT * lgT);

    L0_HD = std::pow(10.0, lgL);
  }

  // Line escape fraction
  double fesc = fesc_HD(Nc_HD, T_K);

  // Specific cooling rate [erg g⁻¹ s⁻¹]
  return L0_HD * y_HD * nH * fesc * std::exp(-tau_cnt) / rho;
}

// ─────────────────────────────────────────────────────────────────────────────
// line_cool — combined line cooling
//   Subtracts CMB background contribution (evaluated at T_rad)
//
//   Species are addressed through the caller-supplied Sp enum (the model's own
//   index space): y[Sp::H], y[Sp::H2], y[Sp::e], y[Sp::Hp], y[Sp::He],
//   y[Sp::HD]. This matters because the compact Nakauchi2019_Minimal lays HD at
//   a different local index (10) than the full/metal networks (12); using the
//   model's Sp keeps every model's HD/H2 cooling correct (full/metal stay
//   bit-for-bit since their indices are unchanged).
// ─────────────────────────────────────────────────────────────────────────────
template <class Sp, int N_sp>
void line_cool(const std::array<double, N_sp>& y, double Nc_H2, double Nc_HD,
               double nH, double T_K, double T_rad, double tau_cnt,
               double& Lambda_line, double& Lambda_H2, double& Lambda_HD,
               double& Lambda_Lya) {
  constexpr double yHe = abundance_ref::yHe;
  constexpr double m_p = phys::m_p;

  double y_a = y[Sp::H];    // H
  double y_m = y[Sp::H2];   // H2
  double y_e = y[Sp::e];    // e-
  double y_Hp = y[Sp::Hp];  // H+
  double y_He = y[Sp::He];  // He
  double y_HD = y[Sp::HD];  // HD

  double rho = nH * ((1.0 + 4.0 * yHe) * m_p);

  // H2 cooling: net = T_K value − T_rad (CMB) value
  Lambda_H2 =
      h2_cooling(nH, T_K, rho, y_a, y_m, y_e, y_Hp, y_He, Nc_H2, tau_cnt) -
      h2_cooling(nH, T_rad, rho, y_a, y_m, y_e, y_Hp, y_He, Nc_H2, tau_cnt);

  // HD cooling: net = T_K value − T_rad (CMB) value
  Lambda_HD = hd_cooling(nH, T_K, rho, y_HD, Nc_HD, tau_cnt) -
              hd_cooling(nH, T_rad, rho, y_HD, Nc_HD, tau_cnt);

  // HI Lya cooling (Osterbrock 2006 fit; only valid below 10^10 cm^-3)
  if (nH <= 1.0e10) {
    double T_5 = 1.0e-5 * T_K;
    Lambda_Lya = y_e * y_a * 7.50e-19 / (1.0 + std::sqrt(T_5)) *
                 std::exp(-1.18348e5 / T_K) * nH / ((1.0 + 4.0 * yHe) * m_p);
  } else {
    Lambda_Lya = 0.0;
  }

  Lambda_line = Lambda_H2 + Lambda_HD + Lambda_Lya;
}

// ─────────────────────────────────────────────────────────────────────────────
// cnt_cool — continuum (dust+gas) radiative cooling
//   Uses Rosseland mean opacity table k_prm (detail::eval_opacity)
// ─────────────────────────────────────────────────────────────────────────────
inline void cnt_cool(double nH, double T_K, double T_rad, double esc_cnt,
                     double& k_gas, double& Lambda_gas, double& Lambda_cnt) {
  constexpr double yHe = abundance_ref::yHe;
  constexpr double m_p = phys::m_p;
  constexpr double sigma_B =
      phys::sigma_B;  // Stefan-Boltzmann [erg cm^-2 s^-1 K^-4]

  double rho = nH * ((1.0 + 4.0 * yHe) * m_p);
  k_gas = detail::eval_opacity(T_K, rho);

  Lambda_gas = 4.0 * sigma_B *
               (T_K * T_K * T_K * T_K - T_rad * T_rad * T_rad * T_rad) * k_gas *
               esc_cnt;
  Lambda_cnt = Lambda_gas;
}

}  // namespace arche
