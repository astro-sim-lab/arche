// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// cooling_grain.h — grain + gas continuum cooling for the metal_grain network
//
// Provides cnt_cool_metal(), extracted from CntCoolMetal() in
// src/apps/collapse_metal_grain/main.cc so that it can be called from
// chemistry.h's chem_full_step().
// ---------------------------------------------------------------------------
#include <cmath>

#include "cooling/grain.h"   // grtemp, detail::kp_gr
#include "kinetics/rates.h"  // detail::eval_opacity

namespace arche {

// ---------------------------------------------------------------------------
// cnt_cool_metal — bisection-solve for grain temperature, then compute
//                  grain + gas continuum cooling rates.
//
// On entry:
//   nH     — H number density [cm^-3]
//   T_K     — gas temperature [K]
//   T_rad   — CMB radiation temperature [K]
//   tau_cnt — continuum optical depth
//   esc_cnt — continuum escape fraction
//   Z_metal — metallicity [Z_sun]
//   T_gr_K  — grain temperature [K]  (initial guess / previous value)
//
// On exit (all output parameters):
//   T_gr_K    — updated grain temperature [K]
//   k_gr     — grain opacity × Z_metal [cm^2/g]
//   k_gas    — gas (ff+CIA) opacity [cm^2/g]
//   Lambda_gr  — grain cooling [erg g^-1 s^-1]
//   Lambda_gas — gas cooling   [erg g^-1 s^-1]
//   Lambda_cnt — total = grain + gas [erg g^-1 s^-1]
// ---------------------------------------------------------------------------
inline void cnt_cool_metal(double nH, double T_K, double T_rad, double tau_cnt,
                           double esc_cnt, double Z_metal, double& T_gr_K,
                           double& k_gr, double& k_gas, double& Lambda_gr,
                           double& Lambda_gas, double& Lambda_cnt) {
  constexpr double kYHe = abundance_ref::yHe;  // He abundance (y_He / y_H)

  double esc_fact;
  if (tau_cnt <= 1.0e-3)
    esc_fact = tau_cnt - 0.5 * tau_cnt * tau_cnt;
  else
    esc_fact = 1.0 - std::exp(-tau_cnt);

  grtemp(nH, T_K, tau_cnt, esc_cnt, esc_fact, T_rad, T_gr_K);

  double rho = nH * ((1.0 + 4.0 * kYHe) * phys::m_p);
  k_gr = detail::kp_gr(rho, T_gr_K) * Z_metal;
  k_gas = detail::eval_opacity(T_K, rho);

  double T_gr4 = T_gr_K * T_gr_K * T_gr_K * T_gr_K;
  double T_r4 = T_rad * T_rad * T_rad * T_rad;
  double T_K4 = T_K * T_K * T_K * T_K;

  Lambda_gr = 4.0 * phys::sigma_B * (T_gr4 - T_r4) * k_gr * esc_cnt;
  Lambda_gas = 4.0 * phys::sigma_B * (T_K4 - T_r4) * k_gas * esc_cnt;
  Lambda_cnt = Lambda_gr + Lambda_gas;
}

// ---------------------------------------------------------------------------
// cnt_cool_metal_consistent — cnt_cool_metal with tau_cnt/esc_cnt solved in-step
//
// Same outputs as cnt_cool_metal, but instead of taking tau_cnt/esc_cnt from the caller
// it takes the column mass rho*L [g cm^-2] and solves the grain temperature and the
// optical depth together (grtemp_consistent), returning the tau/esc it used.
// ---------------------------------------------------------------------------
inline void cnt_cool_metal_consistent_solve(
    double nH, double T_K, double T_rad, double column, double Z_metal,
    double& T_gr_K, double& k_gr, double& k_gas, double& Lambda_gr,
    double& Lambda_gas, double& Lambda_cnt, double& tau_cnt, double& esc_cnt) {
  constexpr double kYHe = abundance_ref::yHe;
  double rho = nH * ((1.0 + 4.0 * kYHe) * phys::m_p);
  k_gas = detail::eval_opacity(T_K, rho);
  double esc_fact = 0.0;
  grtemp_consistent(nH, T_K, column, k_gas, Z_metal, T_rad, T_gr_K, tau_cnt,
                    esc_cnt, esc_fact);
  k_gr = detail::kp_gr(rho, T_gr_K) * Z_metal;

  double T_gr4 = T_gr_K * T_gr_K * T_gr_K * T_gr_K;
  double T_r4 = T_rad * T_rad * T_rad * T_rad;
  double T_K4 = T_K * T_K * T_K * T_K;

  Lambda_gr = 4.0 * phys::sigma_B * (T_gr4 - T_r4) * k_gr * esc_cnt;
  Lambda_gas = 4.0 * phys::sigma_B * (T_K4 - T_r4) * k_gas * esc_cnt;
  Lambda_cnt = Lambda_gr + Lambda_gas;
}

// The continuum channel: cnt_cool_metal_consistent_solve keeping its last results per
// thread.  It is a pure function of its arguments and the grain-temperature guess
// T_gr_K, and the integrator repeats it at identical inputs (the start of a call and of
// its first attempt, and the full and first half steps of one attempt).  Defined out of
// line in cooling/cooling_memo.cc, so every caller shares one compiled evaluation and a
// kept result is exactly what that evaluation returns (core/exact_memo.h).
void cnt_cool_metal_consistent(double nH, double T_K, double T_rad,
                               double column, double Z_metal, double& T_gr_K,
                               double& k_gr, double& k_gas, double& Lambda_gr,
                               double& Lambda_gas, double& Lambda_cnt,
                               double& tau_cnt, double& esc_cnt);

}  // namespace arche
