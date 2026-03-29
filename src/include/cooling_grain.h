#pragma once
// ---------------------------------------------------------------------------
// cooling_grain.h — grain + gas continuum cooling for the metal_grain network
//
// Provides cnt_cool_metal(), extracted from CntCoolMetal() in
// src/apps/collapse_metal_grain/main.cc so that it can be called from
// chemistry.h's chem_full_step().
//
// Port of the continuum cooling block in:
//   fortran/metal_grain/collapse_metal_grain.f
// ---------------------------------------------------------------------------
#include <cmath>
#include "grain.h"    // grtemp, detail::xkp_gr
#include "reaction.h" // detail::eval_opacity

namespace chemistry {

// ---------------------------------------------------------------------------
// cnt_cool_metal — bisection-solve for grain temperature, then compute
//                  grain + gas continuum cooling rates.
//
// On entry:
//   xnH     — H number density [cm^-3]
//   T_K     — gas temperature [K]
//   T_rad   — CMB radiation temperature [K]
//   tau_cnt — continuum optical depth
//   esc_cnt — continuum escape fraction
//   Z_metal — metallicity [Z_sun]
//   T_gr_K  — grain temperature [K]  (initial guess / previous value)
//
// On exit (all output parameters):
//   T_gr_K    — updated grain temperature [K]
//   xk_gr     — grain opacity × Z_metal [cm^2/g]
//   xk_gas    — gas (ff+CIA) opacity [cm^2/g]
//   xLmbd_gr  — grain cooling [erg g^-1 s^-1]
//   xLmbd_gas — gas cooling   [erg g^-1 s^-1]
//   xLmbd_cnt — total = grain + gas [erg g^-1 s^-1]
// ---------------------------------------------------------------------------
inline void cnt_cool_metal(double xnH,   double T_K,   double T_rad,
                           double tau_cnt, double esc_cnt,
                           double Z_metal,
                           double& T_gr_K,
                           double& xk_gr,
                           double& xk_gas,
                           double& xLmbd_gr,
                           double& xLmbd_gas,
                           double& xLmbd_cnt)
{
    constexpr double kYHe = abundance_ref::yHe;  // He abundance (y_He / y_H)

    double esc_fact;
    if (tau_cnt <= 1.0e-3)
        esc_fact = tau_cnt - 0.5*tau_cnt*tau_cnt;
    else
        esc_fact = 1.0 - std::exp(-tau_cnt);

    grtemp(xnH, T_K, tau_cnt, esc_cnt, esc_fact, T_rad, T_gr_K);

    double rho = xnH * ((1.0 + 4.0*kYHe) * phys::xm_p);
    xk_gr  = detail::xkp_gr(rho, T_gr_K) * Z_metal;
    xk_gas = detail::eval_opacity(T_K, rho);

    double T_gr4 = T_gr_K*T_gr_K*T_gr_K*T_gr_K;
    double T_r4  = T_rad*T_rad*T_rad*T_rad;
    double T_K4  = T_K*T_K*T_K*T_K;

    xLmbd_gr  = 4.0 * phys::sigma_B * (T_gr4 - T_r4) * xk_gr  * esc_cnt;
    xLmbd_gas = 4.0 * phys::sigma_B * (T_K4  - T_r4) * xk_gas * esc_cnt;
    xLmbd_cnt = xLmbd_gr + xLmbd_gas;
}

} // namespace chemistry
