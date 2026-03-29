#pragma once
// ---------------------------------------------------------------------------
// grain.h — Dust grain physics for metal_grain chemistry
//
// Port of:
//   fortran/metal_grain/subs/grain_new.f  (grtemp, xkp_gr, vol_gr, vaptemp,
//                                           phelectr)
//
// All functions operate in CGS.
// ---------------------------------------------------------------------------
#include <cmath>
#include <algorithm>
#include <limits>
#include "state.h"

namespace chemistry {

// ---------------------------------------------------------------------------
// GrainState — grain temperature passed between cnt_cool and chemreact.
// Kept outside ChemState to avoid mixing kernel-internal vs app-layer state.
// ---------------------------------------------------------------------------
struct GrainState {
    // NaN default to make missing caller initialization visible at runtime.
    double T_gr_K  = std::numeric_limits<double>::quiet_NaN(); // grain temperature [K]
    double xk_gr   = 0.0;     // Z-weighted Planck opacity [cm^2/g]
    double vol_gr0 = 0.0;     // vol_gr at T_gr→0 (reference, set at init)
};

namespace detail {

// ---------------------------------------------------------------------------
// vaptemp — density-interpolated vaporization temperatures
// Pollack et al. 1994; 8-point log-density table
// ---------------------------------------------------------------------------
inline void vaptemp(double rho_in,
                    double& T_ice, double& T_vo, double& T_ro,
                    double& T_tr,  double& T_ir, double& T_pyr,
                    double& T_ol)
{
    // log-spaced density nodes [g/cm^3]
    static constexpr double ro[8] = {
        1.0e-18, 1.0e-16, 1.0e-14, 1.0e-12,
        1.0e-10, 1.0e-08, 1.0e-06, 1.0e-04
    };
    // tt[col][row]: 4 columns = (water ice, metallic iron, orthopyroxene, olivine)
    // 8 rows = density nodes
    static constexpr double tt[4][8] = {
        // Water ice evaporation temperature
        { 1.090e2, 1.180e2, 1.290e2, 1.430e2,
          1.590e2, 1.800e2, 2.070e2, 2.440e2 },
        // Metallic iron
        { 8.350e2, 9.080e2, 9.940e2, 1.100e3,
          1.230e3, 1.395e3, 1.612e3, 1.908e3 },
        // Orthopyroxene
        { 9.020e2, 9.800e2, 1.049e3, 1.129e3,
          1.222e3, 1.331e3, 1.462e3, 1.621e3 },
        // Olivine
        { 9.290e2, 9.970e2, 1.076e3, 1.168e3,
          1.277e3, 1.408e3, 1.570e3, 1.774e3 }
    };

    // Linear interpolation in log-rho space (clamp at edges)
    auto interp = [&](const double* tab) -> double {
        if (rho_in <= ro[0]) return tab[0];
        if (rho_in >= ro[7]) return tab[7];
        for (int j = 0; j < 7; ++j) {
            if (rho_in >= ro[j] && rho_in <= ro[j+1]) {
                double frac = (rho_in - ro[j]) / (ro[j+1] - ro[j]);
                return tab[j] + frac * (tab[j+1] - tab[j]);
            }
        }
        return tab[7];
    };

    T_ice = interp(tt[0]);
    T_ir  = interp(tt[1]);
    T_pyr = interp(tt[2]);
    T_ol  = interp(tt[3]);
    // Fixed values (Pollack et al. 1994)
    T_vo  = 275.0;
    T_ro  = 425.0;
    T_tr  = 680.0;
}

// ---------------------------------------------------------------------------
// xkp_gr — Planck mean mass absorption coefficient [cm^2/g] for solar-Z dust
// Semenov et al. 2003; 5 polynomial segments blended with tanh
// ---------------------------------------------------------------------------
inline double xkp_gr(double rho, double T_in)
{
    // eP[seg][coeff] — 5 temperature segments, degree-5 polynomial
    // Coefficients ordered eP[seg][0..5] = a0..a5 (ascending T power)
    static constexpr double eP[5][6] = {
        // Region 1 (lowest T — below water-ice evap)
        { -0.664969066656727e-10,
           0.427373694560880e-07,
          -0.108371821805323e-04,
           0.130732588474620e-02,
          -0.191556936842122e-02,
           0.423062838413742e-03 },
        // Region 2
        { -0.371134628305626e-11,
           0.472321713029246e-08,
          -0.210704968520099e-05,
           0.311604311624702e-03,
           0.302186734201724e-01,
          -0.173587657890234e+01 },
        // Region 3
        {  0.233228177925061e-12,
          -0.681138400996940e-09,
           0.784339065020565e-06,
          -0.436299138970822e-03,
           0.120954274022502e+00,
          -0.638046050114383e+01 },
        // Region 4
        {  0.115933380913977e-12,
          -0.349467552126090e-09,
           0.417313950448598e-06,
          -0.240971082612604e-03,
           0.664915248499724e-01,
          -0.345042341508906e+01 },
        // Region 5 (highest T)
        { -0.503398974713655e-15,
           0.417753047583439e-11,
          -0.140583908595144e-07,
           0.238329845360675e-04,
          -0.166789974601131e-01,
           0.667823366719512e+01 }
    };

    double T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol;
    vaptemp(rho, T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol);

    double Tv[5];
    Tv[0] = T_ice;
    Tv[1] = T_vo;
    Tv[2] = T_ro;
    Tv[3] = T_tr;
    double tmax1 = std::max(T_ir, T_pyr);
    double tmax2 = std::max(T_pyr, T_ol);
    Tv[4] = std::min(tmax1, tmax2);

    // Evaluate each polynomial at T_in
    double xk[5];
    for (int s = 0; s < 5; ++s) {
        xk[s] = ((((eP[s][0]*T_in + eP[s][1])*T_in + eP[s][2])*T_in
                    + eP[s][3])*T_in + eP[s][4])*T_in + eP[s][5];
    }

    // Blend with tanh step functions at vaporization temperatures
    auto step = [](double T, double Tv_) {
        return 0.5 * (std::tanh(40.0 * (1.0 - T / Tv_)) + 1.0);
    };

    return (xk[0]-xk[1]) * step(T_in, Tv[0])
         + (xk[1]-xk[2]) * step(T_in, Tv[1])
         + (xk[2]-xk[3]) * step(T_in, Tv[2])
         + (xk[3]-xk[4]) * step(T_in, Tv[3])
         +  xk[4]        * step(T_in, Tv[4]);
}

// ---------------------------------------------------------------------------
// vol_gr — grain volume per unit mass of gas [cm^3/g] at solar metallicity
// ---------------------------------------------------------------------------
inline double vol_gr(double rho, double T_K)
{
    double T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol;
    vaptemp(rho, T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol);

    auto x = [](double T, double Tv_) {
        return 0.5 * (std::tanh(40.0 * (1.0 - T / Tv_)) + 1.0);
    };

    return (0.78e-4 - 0.46e-4 * x(T_K, T_tr)) * x(T_K, T_ir)
         + 7.19e-4 * x(T_K, T_ol)
         + 2.16e-4 * x(T_K, T_pyr)
         + 1.18e-4 * x(T_K, T_tr)
         + 23.53e-4 * x(T_K, T_ro)
         + 6.02e-4  * x(T_K, T_vo)
         + 12.93e-4 * x(T_K, T_ice);
}

} // namespace detail

// ---------------------------------------------------------------------------
// grtemp — solve for grain temperature by bisection
//
// Solves func(T_gr) = 0 where:
//   func(Td) = xkp_gr(rho,Td)*Td^4*esc_cnt
//            - xnH * funcL(Td)
//            - xkp_gr(rho,Td)*(T_gas^4*esc_fact + T_rad^4*exp(-tau))
//   funcL(Td) = 4.9e-2 * vol_gr(rho,Td) * sqrt(T_gas*1e-3)
//               * (1 - 0.8*exp(-75/T_gas)) * (T_gas - Td)
//
// On entry:  T_gr is the initial guess (updated in-place).
// On exit:   T_gr = solution, clamped to min(T_gr, T_gas).
// ---------------------------------------------------------------------------
inline void grtemp(double xnH, double T_gas,
                   double tau_cnt, double esc_cnt, double esc_fact,
                   double T_rad,
                   double& T_gr)
{
    constexpr double yHe = abundance_ref::yHe;
    constexpr double xm_p = phys::legacy::xm_p;
    const double rho = (1.0 + 4.0*yHe) * xm_p * xnH;

    auto funcL = [&](double Td) {
        return 4.9e-2 * detail::vol_gr(rho, Td)
             * std::sqrt(T_gas * 1.0e-3)
             * (1.0 - 0.8 * std::exp(-75.0 / T_gas))
             * (T_gas - Td);
    };
    auto func = [&](double Td) {
        double kp = detail::xkp_gr(rho, Td);
        return kp * std::pow(Td, 4.0) * esc_cnt
             - xnH * funcL(Td)
             - kp * (std::pow(T_gas, 4.0) * esc_fact
                     + std::pow(T_rad, 4.0) * std::exp(-tau_cnt));
    };

    double f0 = func(T_gr);
    if (f0 == 0.0) { T_gr = T_gas; return; }  // dust evaporated: kp≈0,vol_gr≈0 → T_gr=T_gas (matches Fortran)

    double x_l = T_gr, x_u = T_gr;
    if (f0 > 0.0) {
        for (int i = 0; i < 1000; ++i) {
            x_l *= 0.99;
            if (func(x_l) < 0.0) break;
        }
    } else {
        for (int i = 0; i < 1000; ++i) {
            x_u *= 1.01;
            if (func(x_u) > 0.0) break;
        }
    }

    double fmid = func(x_l);
    double f    = func(x_u);
    if (f == 0.0) { T_gr = std::min(x_u, T_gas); return; }
    if (f * fmid >= 0.0) {
        // bracketing failed — return initial guess
        T_gr = std::min(T_gr, T_gas);
        return;
    }

    double dx;
    if (f < 0.0) { T_gr = x_u; dx = x_l - x_u; }
    else         { T_gr = x_l; dx = x_u - x_l; }

    for (int j = 0; j < 100; ++j) {
        dx *= 0.5;
        double xmid = T_gr + dx;
        double fm   = func(xmid);
        if (fm <= 0.0) T_gr = xmid;
        if (T_gr != 0.0) {
            if (std::abs(dx / T_gr) < 1.0e-12 || fm == 0.0) break;
        } else {
            if (fm == 0.0) break;
        }
    }
    T_gr = std::min(T_gr, T_gas);
}

// ---------------------------------------------------------------------------
// phelectr — photoelectric heating / cooling [erg/g/s]
// Bakes & Tielens 1994 fitting formula
//
// D   = vol_gr(rho,T_K) / vol_gr(rho,0) * Z_metal  (grain volume fraction)
// x   = G_0 * sqrt(T_K) / (y_e * xnH)
// eps = epsilon efficiency
// ---------------------------------------------------------------------------
inline void phelectr(double xnH, double T_K, double y_e,
                     double Z_metal, double G_0, double A_v,
                     double& Gam_pe, double& xLmbd_pe)
{
    constexpr double yHe = abundance_ref::yHe;
    constexpr double xm_p = phys::legacy::xm_p;
    const double rho = (1.0 + 4.0*yHe) * xm_p * xnH;

    // Grain volume fraction (relative to T→0 reference)
    double vol0 = detail::vol_gr(rho, 0.0);
    double D    = (vol0 > 0.0) ? detail::vol_gr(rho, T_K) / vol0 * Z_metal
                               : 0.0;

    double x   = G_0 * std::sqrt(T_K) / (y_e * xnH + 1.0e-30);
    double eps  = 4.87e-2 / (1.0 + 4.0e-3 * std::pow(x, 0.73))
                + 3.65e-2 * std::pow(T_K * 1.0e-4, 0.7) / (1.0 + 2.0e-4 * x);

    double norm  = 1.0 / (1.0 + 4.0*yHe);
    Gam_pe   = 0.60 * eps * G_0 * std::exp(-1.8 * A_v) * D * norm;
    double beta  = 0.735 / std::pow(T_K, 0.068);
    xLmbd_pe = 2.78e-6 * std::pow(T_K, 0.944) * std::pow(x, beta)
             * (y_e * xnH) * D * norm;
}

} // namespace chemistry
