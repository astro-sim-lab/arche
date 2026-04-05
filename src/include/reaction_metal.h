// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// reaction_metal.h
//
// Template specialization of compute_base_rates<89, 1200> for the
// metal_grain chemistry network.
//
// Port of react_coef() in:
//   fortran/metal_grain/subs/reaction_metal_grain.f  (lines 122–2986)
//
// Includes:
//   - H/He/D reactions  xk[0..99]   (GA08, Glover&Abel 2008)
//   - Metal reactions   xk[100..542] (UMIST 2012)
//   - Additional        xk[600..644]
//   - Li reactions      xk[800..829] (Bovino+2011, Lepp+2002, Mizusawa+2005)
//   - K/Na/Mg reactions xk[700..729] (KIDA)
//   - CR reactions      xk[543..551], xk[655..681]
//   - Grain charge/surface rates mapped from react_grain_rates / grain_coef_rates
//   - Reverse rates xk[1200..2399] via detailed balance
//
// Index convention: Fortran xk(N) → C++ xk[N-1]  (0-based)
// ---------------------------------------------------------------------------
#include <array>
#include <cmath>
#include <algorithm>
#include "reaction.h"
#include "grain.h"
#include "reaction_metal_grain_charge.h"
#include "reaction_metal_grain_surface.h"

namespace chemistry {

template<>
inline void compute_base_rates<89, 1200>(
    double xnH, double T_K, double xmu,
    const ChemParams& params,
    const ReactionTable<89, 1200>& tbl,
    std::array<double, 2400>& xk)
{
    xk.fill(0.0);

    // Physical constants (CGS)
    constexpr double xk_B = phys::xk_B;
    constexpr double h_P  = phys::h_P;
    constexpr double pi   = phys::pi;
    constexpr double G    = phys::G;
    constexpr double xm_p = phys::xm_p;
    constexpr double yHe  = abundance_ref::yHe;

    // External parameters
    const double zeta    = params.zeta;
    const double T_gr_K  = params.T_gr_K;
    const double T_cr_desorp = params.T_cr_desorp;
    const double Z_metal = params.Z_metal;
    const double xJH2    = params.xJH2;
    const double xJH2O   = params.xJH2O;
    const double xJtot   = params.xJtot;

    // Temperature variables
    double T_eV    = phys::k_B_eV * T_K;
    double xlnT_eV = std::log(T_eV);
    double T300    = T_K / 300.0;
    double xlnT    = std::log(T_K);
    double xlgT    = std::log10(T_K);
    double xlgT2   = xlgT*xlgT;
    double xlgT3   = xlgT2*xlgT;
    double xlgT4   = xlgT3*xlgT;
    double xlgT5   = xlgT4*xlgT;
    double xlgT6   = xlgT5*xlgT;
    double xlT4    = std::log10(T_K / 1.0e4);

    // Critical densities for LTE interpolation (GA08)
    double xncr_H  = std::pow(10.0,  3.0 - 0.416*xlT4 - 0.327*xlT4*xlT4);
    double xncr_H2 = std::pow(10.0,  4.845 - 1.3*xlT4 + 1.62*xlT4*xlT4);
    double xncr_He = std::pow(10.0,  5.0792*(1.0 - 1.23e-5*(T_K - 2000.0)));
    double xncr    = 1.0 / (params.xH/xncr_H + params.xH2/xncr_H2
                                               + params.xHe/xncr_He);
    double xcr     = xnH / xncr;
    double xncr_HD = 1.0e2 * xncr;
    double xcr_HD  = xnH / xncr_HD;

    // Density, Jeans length, optical depth
    double rho     = xnH * (1.0 + 4.0*yHe) * xm_p;
    double xlmbd_J = std::sqrt(pi * xk_B * T_K / (G * xmu * xm_p * rho));
    double kap     = detail::eval_opacity(T_K, rho);
    double xk_gr   = detail::xkp_gr(rho, T_gr_K) * Z_metal;
    double tau_cnt = (kap + xk_gr) * rho * xlmbd_J;

    // Partition functions
    std::array<double, 102> pf;
    eval_partition_functions<89, 1200>(T_K, tbl, pf);

    // -----------------------------------------------------------------------
    //  H, He reactions  xk[0..45]  (Reactions 1..46, GA08)
    // -----------------------------------------------------------------------
    // 1) H + e -> H+ + 2e
    if (T_K < 1000.0)
        xk[0] = 0.0;
    else
        xk[0] = std::exp(-32.71396786
            + xlnT_eV*(13.536556
            + xlnT_eV*(-5.73932875
            + xlnT_eV*(1.56315498
            + xlnT_eV*(-0.2877056
            + xlnT_eV*(3.48255977e-2
            + xlnT_eV*(-2.63197617e-3
            + xlnT_eV*(1.11954395e-4
            + xlnT_eV*(-2.03914985e-6)))))))));

    // 2) H+ + e -> H + ph.  (case B)
    xk[1] = 2.753e-14 * std::pow(315614.0/T_K, 1.500)
            * std::pow(1.0 + std::pow(115188.0/T_K, 0.407), -2.242);

    // 3) He + e -> He+ + 2e
    if (T_K < 1000.0)
        xk[2] = 0.0;
    else
        xk[2] = std::exp(-44.09864886
            + xlnT_eV*(23.91596563
            + xlnT_eV*(-10.7532302
            + xlnT_eV*(3.05803875
            + xlnT_eV*(-0.56851189
            + xlnT_eV*(6.79539123e-2
            + xlnT_eV*(-5.00905610e-3
            + xlnT_eV*(2.06723616e-4
            + xlnT_eV*(-3.64916141e-6)))))))));

    // 4) He+ + e -> He + ph.  (optically thick)
    {
        double xkrrA = 1.0e-11*std::pow(T_K,-0.5)
                      *(12.72 - 1.615*xlgT - 0.3162*xlgT2 + 0.0493*xlgT3);
        double xkrrB = 1.0e-11*std::pow(T_K,-0.5)
                      *(11.19 - 1.676*xlgT - 0.2852*xlgT2 + 0.04433*xlgT3);
        double xkdi  = 1.9e-3*std::pow(T_K,-1.5)*std::exp(-473421.0/T_K)
                      *(1.0 + 0.3*std::exp(-94684.0/T_K));
        xk[3] = 0.68*xkrrA + 0.32*xkrrB + xkdi;
    }

    // 5) He+ + e -> He++ + 2e
    if (T_K < 1000.0)
        xk[4] = 0.0;
    else
        xk[4] = std::exp(-68.71040990
            + xlnT_eV*(43.93347633
            + xlnT_eV*(-18.4806699
            + xlnT_eV*(4.70162649
            + xlnT_eV*(-0.76924663
            + xlnT_eV*(8.113042e-2
            + xlnT_eV*(-5.32402063e-3
            + xlnT_eV*(1.97570531e-4
            + xlnT_eV*(-3.16558106e-6)))))))));

    // 6) He++ + e -> He+ + ph.  (case B)
    xk[5] = 5.506e-14 * std::pow(1262456.0/T_K, 1.500)
            * std::pow(1.0 + std::pow(460752.0/T_K, 0.407), -2.242);

    // 7) H + e -> H- + ph.
    if (T_K < 6000.0)
        xk[6] = std::pow(10.0, -17.845 + 0.762*xlgT
                              + 0.1523*xlgT2 - 0.03274*xlgT3);
    else
        xk[6] = std::pow(10.0, -16.4199 + 0.1998*xlgT2
                              - 5.447e-3*xlgT4 + 4.0415e-5*xlgT6);

    // 8) H- + H -> H2 + e  (Kreckel+2010)
    xk[7] = 1.35e-9 * (std::pow(T_K,9.8493e-2)
                       + 3.2852e-1*std::pow(T_K,5.561e-1)
                       + 2.771e-7*std::pow(T_K,2.1826))
           / (1.0 + 6.191e-3*std::pow(T_K,1.0461)
                  + 8.9712e-11*std::pow(T_K,3.0424)
                  + 3.2576e-14*std::pow(T_K,3.7741));

    // 9) H + H+ -> H2+ + ph.  (Coppola+2011, Glover+2015)
    xk[8] = std::pow(10.0, -18.2 - 3.194*xlgT + 1.786*xlgT2 - 0.2072*xlgT3);

    // 10) H2+ + H -> H2 + H+
    xk[9] = 6.4e-10;

    // 11) H2 + H+ -> H2+ + H  (GA08 7)
    if (T_K < 100.0)
        xk[10] = 0.0;
    else
        xk[10] = (-3.3232183e-7 + 3.3735382e-7*xlnT
                 - 1.4491368e-7*xlnT*xlnT
                 + 3.4172805e-8*std::pow(xlnT,3)
                 - 4.7813720e-9*std::pow(xlnT,4)
                 + 3.9731542e-10*std::pow(xlnT,5)
                 - 1.8171411e-11*std::pow(xlnT,6)
                 + 3.5311932e-13*std::pow(xlnT,7))
                * std::exp(-21237.15/T_K);

    // 12) H2 + e -> 2H + e  (GA08 8, LTE/low-density blend)
    double xk12v0  = 4.49e-9 * std::pow(T_K,0.11) * std::exp(-101858.0/T_K);
    double xk12LTE = 1.91e-9 * std::pow(T_K,0.136)* std::exp(-53407.1/T_K);
    if (xk12v0 == 0.0 || xk12LTE == 0.0)
        xk[11] = 0.0;
    else {
        double xlgk12 = (xcr/(1.0+xcr))*std::log10(xk12LTE)
                      + (1.0/(1.0+xcr))*std::log10(xk12v0);
        xk[11] = std::pow(10.0, xlgk12);
    }

    // 13) H2 + H -> 3H  (GA08 9, LTE/low-density blend)
    double xk13v0  = 6.67e-12 * std::sqrt(T_K) * std::exp(-(1.0 + 63593.0/T_K));
    double xk13LTE = 3.52e-9 * std::exp(-43900.0/T_K);
    {
        double xlgk13 = (xcr/(1.0+xcr))*std::log10(xk13LTE)
                      + (1.0/(1.0+xcr))*std::log10(xk13v0);
        xk[12] = std::pow(10.0, xlgk13);
    }

    // 14) H- + e -> H + 2e  (GA08 14)
    if (T_K <= 40.0)
        xk[13] = 0.0;
    else
        xk[13] = std::exp(-18.01849334
            + xlnT_eV*(2.3608522
            + xlnT_eV*(-0.28274430
            + xlnT_eV*(1.62331664e-2
            + xlnT_eV*(-3.36501203e-2
            + xlnT_eV*(1.17832978e-2
            + xlnT_eV*(-1.65619470e-3
            + xlnT_eV*(1.06827520e-4
            + xlnT_eV*(-2.63128581e-6)))))))));

    // 15) H- + H+ -> 2H  (GA08 5)
    xk[14] = 2.4e-6 * std::pow(T_K,-0.5) * (1.0 + T_K/20000.0);

    // 16) H- + H+ -> H2+ + e  (GA08 16)
    xk[15] = (T_K <= 8000.0) ? 6.9e-9*std::pow(T_K,-0.35)
                              : 9.6e-7*std::pow(T_K,-0.90);

    // 17) H2+ + e -> 2H  (GA08 6)
    xk[16] = (T_K < 617.0) ? 1.0e-8 : 1.32e-6*std::pow(T_K,-0.76);

    // 18) H2+ + H- -> H2 + H  (GA08 21)
    xk[17] = 1.4e-7 * std::pow(T300,-0.5);

    // 19) 3H -> H2 + H  (Forrey+2013, Glover+2015)
    double xk19v0 = 6.0e-32*std::pow(T_K,-0.25) + 2.0e-31*std::pow(T_K,-0.5);
    xk[18] = xk19v0;

    // 20) 2H + H2 -> 2H2  (GA08 31)
    xk[19] = xk[18] / 8.0;

    // 21) 2H2 -> 2H + H2  (GA08 10, LTE/low-density blend)
    double xk21v0  = (5.996e-30*std::pow(T_K,4.1881)
                      / std::pow(1.0+6.761e-6*T_K, 5.6881))
                    * std::exp(-54657.4/T_K);
    double xk21LTE = 1.3e-9 * std::exp(-53300.0/T_K);
    if (xk21v0 == 0.0 || xk21LTE == 0.0)
        xk[20] = 0.0;
    else {
        double xlgk21 = (xcr/(1.0+xcr))*std::log10(xk21LTE)
                      + (1.0/(1.0+xcr))*std::log10(xk21v0);
        xk[20] = std::pow(10.0, xlgk21);
    }

    // 22) H + H -> H + e + H+  (Lenzuni+1991)
    xk[21] = 1.2e-17 * std::pow(T_K,1.2)
           * std::exp(-std::abs(tbl.delE[21]) / (xk_B*T_K));

    // 23) 2H + grain -> H2  (grain surface; set by grain_coef_rates)
    xk[22] = 0.0;

    // 24) He+ + H2 -> H+ + H + He  (GA08 24)
    xk[23] = 3.70e-14 * std::exp(-35.0/T_K);

    // 25) H2+ + He -> HeH+ + H  (GS09 TR33)
    xk[24] = 3.0e-10 * std::exp(-6717.0/T_K);

    // 26) H2+ + H2 -> H3+ + H  (GS09 TR1)
    xk[25] = 2.24e-9 * std::pow(T300,0.042) * std::exp(-T_K/4.66e4);

    // 27) H3+ + H- -> 2H2
    xk[26] = 2.30e-7 * std::pow(T300,-0.50);

    // 28) He+ + H -> H+ + He  (GA08 26)
    xk[27] = 1.2e-15 * std::pow(T300,0.25);

    // 29) He+ + H- -> H + He  (GA08 28)
    xk[28] = 2.32e-7 * std::pow(T300,-0.52) * std::exp(T_K/22400.0);

    // 30) He+ + H2 -> H2+ + He  (GA08 25)
    xk[29] = 7.20e-15;

    // 31) HeH+ + H -> H2+ + He  (GS09 TR37)
    xk[30] = 1.04e-9 * std::pow(T300,0.13) * std::exp(-T_K/3.31e4);

    // 32) HeH+ + H2 -> H3+ + He  (GS09 TR39)
    xk[31] = 1.53e-9 * std::pow(T300,0.24) * std::exp(-T_K/1.48e4);

    // 33) H2+ + H- -> H3+ + e  (GS09 AD13)
    xk[32] = 2.7e-10 * std::pow(T300,-0.485) * std::exp(T_K/3.12e4);

    // 34) H3+ + e -> H2 + H
    xk[33] = 2.34e-8 * std::pow(T300,-0.52);

    // 35) H3+ + e -> 3H
    xk[34] = 4.36e-8 * std::pow(T300,-0.52);

    // 36) HeH+ + e -> H + He  (GS09 DR14)
    xk[35] = 3.0e-8 * std::pow(T300,-0.47);

    // 37) H2 + He -> 2H + He  (GA08 11, LTE/low-density blend)
    double xk37v0  = std::pow(10.0, -27.029+3.801*xlgT-29487.0/T_K);
    double xk37LTE = std::pow(10.0, -2.729-1.75*xlgT-23474.0/T_K);
    if (xk37v0 == 0.0 || xk37LTE == 0.0)
        xk[36] = 0.0;
    else {
        double xlgk37 = (xcr/(1.0+xcr))*std::log10(xk37LTE)
                      + (1.0/(1.0+xcr))*std::log10(xk37v0);
        xk[36] = std::pow(10.0, xlgk37);
    }

    // 38) H- + H -> 2H + e  (GA08 15, zero below 0.1 eV)
    if (T_eV < 0.1)
        xk[37] = 0.0;
    else
        xk[37] = std::exp(-20.372609
            + xlnT_eV*(1.13944933
            + xlnT_eV*(-1.4210135e-1
            + xlnT_eV*(8.4644554e-3
            + xlnT_eV*(-1.4327641e-3
            + xlnT_eV*(2.0122503e-4
            + xlnT_eV*(8.6639632e-5
            + xlnT_eV*(-2.5850097e-5
            + xlnT_eV*(2.4555012e-6
            + xlnT_eV*(-8.0683825e-8))))))))));

    // 39) H- + H2+ -> 3H  (GA08 22)
    xk[38] = 1.4e-7 * std::pow(T300,-0.5);

    // 40) H2 + e -> H- + H  (GA08 23)
    double xk40v0 = 2.7e-8 * std::pow(T_K,-1.27) * std::exp(-43000.0/T_K);
    xk[39] = xk40v0;   // stored; zeroed by zero-forcing below

    // 41) He + H+ -> He+ + H  (GA08 27)
    xk[40] = (T_K < 1.0e4)
           ? 1.26e-9*std::pow(T_K,-0.75)*std::exp(-127500.0/T_K)
           : 4.0e-37*std::pow(T_K,4.74);

    // 42) He + H- -> He + H + e  (GA08 29)
    xk[41] = 4.1e-17 * T_K*T_K * std::exp(-19870.0/T_K);

    // 43) 2H + He -> H2 + He  (GA08 32)
    xk[42] = 6.9e-32 * std::pow(T_K,-0.4);

    // 44) H + H2+ -> H3+ + ph.  (GS09 RA6)
    xk[43] = 1.5e-17 * std::pow(T300,1.8) * std::exp(20.0/T_K);

    // 45) H2 + H+ -> H3+ + ph.  (Stancil+98)
    xk[44] = 1.0e-20;

    // 46) He + H+ -> HeH+ + ph.  (GS09 RA25)
    xk[45] = 8.0e-20 * std::pow(T300,-0.24) * std::exp(-T_K/4.0e3);

    // -----------------------------------------------------------------------
    //  H/He zero-forcing (L513–524)
    // -----------------------------------------------------------------------
    xk[10] = 0.0;   // 11
    xk[12] = 0.0;   // 13
    xk[19] = 0.0;   // 20
    xk[24] = 0.0;   // 25
    xk[39] = 0.0;   // 40
    xk[40] = 0.0;   // 41
    xk[42] = 0.0;   // 43

    // -----------------------------------------------------------------------
    //  D reactions  xk[46..99]  (Reactions 47..100, GA08)
    // -----------------------------------------------------------------------
    // 47) HD + e -> H + D + e  (GA08 111, LTE/low blend with xcr_HD)
    {
        double xk47v0  = 5.09e-9*std::pow(T_K,0.128)*std::exp(-103258.0/T_K);
        double xk47LTE = 1.04e-9*std::pow(T_K,0.218)*std::exp(-53070.7/T_K);
        if (xk47v0 == 0.0 || xk47LTE == 0.0)
            xk[46] = 0.0;
        else {
            double xlgk47 = (xcr_HD/(1.0+xcr_HD))*std::log10(xk47LTE)
                          + (1.0/(1.0+xcr_HD))*std::log10(xk47v0);
            xk[46] = std::pow(10.0, xlgk47);
        }
    }

    // 48) HD + He -> H + D + He  (GA08 110, uses xk37 branches)
    if (xk37v0 == 0.0 || xk37LTE == 0.0)
        xk[47] = 0.0;
    else {
        double xlgk48 = (xcr_HD/(1.0+xcr_HD))*std::log10(xk37LTE)
                      + (1.0/(1.0+xcr_HD))*std::log10(xk37v0);
        xk[47] = std::pow(10.0, xlgk48);
    }

    // 49) HD + H2 -> H + D + H2  (GA08 109, uses xk21 branches)
    if (xk21v0 == 0.0 || xk21LTE == 0.0)
        xk[48] = 0.0;
    else {
        double xlgk49 = (xcr_HD/(1.0+xcr_HD))*std::log10(xk21LTE)
                      + (1.0/(1.0+xcr_HD))*std::log10(xk21v0);
        xk[48] = std::pow(10.0, xlgk49);
    }

    // 50) HD + H -> 2H + D  (GA08 108, uses xk13 branches)
    if (xk13v0 == 0.0 || xk13LTE == 0.0)
        xk[49] = 0.0;
    else {
        double xlgk50 = (xcr_HD/(1.0+xcr_HD))*std::log10(xk13LTE)
                      + (1.0/(1.0+xcr_HD))*std::log10(xk13v0);
        xk[49] = std::pow(10.0, xlgk50);
    }

    // 51) D+ + e -> D + ph.  = xk[1]
    xk[50] = xk[1];

    // 52) D + H+ -> D+ + H  (GA08 34)
    if (T_K < 2.0e5)
        xk[51] = 2.0e-10*std::pow(T_K,0.402)*std::exp(-37.1/T_K)
                 - 3.31e-17*std::pow(T_K,1.48);
    else
        xk[51] = 3.44e-10*std::pow(T_K,0.35);

    // 53) D+ + H -> D + H+  (GA08 35)
    xk[52] = 2.06e-10*std::pow(T_K,0.396)*std::exp(-33.0/T_K)
           + 2.03e-9*std::pow(T_K,-0.332);

    // 54) D + H -> HD + ph.  (GA08 36)
    if (T_K > 10.0 && T_K <= 200.0)
        xk[53] = 1.0e-25*(2.80202 - 6.63697*xlnT + 4.75619*xlnT*xlnT
                          - 1.39325*std::pow(xlnT,3) + 0.178259*std::pow(xlnT,4)
                          - 0.00817097*std::pow(xlnT,5));
    else if (T_K > 200.0)
        xk[53] = 1.0e-25*std::exp(507.207 - 370.889*xlnT + 104.854*xlnT*xlnT
                                  - 14.4192*std::pow(xlnT,3)
                                  + 0.971469*std::pow(xlnT,4)
                                  - 0.0258076*std::pow(xlnT,5));
    else
        xk[53] = 0.0;

    // 55) D + H2 -> H + HD  (GA08 37)
    if (T_K < 2000.0)
        xk[54] = std::pow(10.0, -56.4737 + 5.88886*xlgT + 7.19692*xlgT2
                               + 2.25069*xlgT3 - 2.16903*xlgT4
                               + 0.317887*xlgT5);
    else
        xk[54] = 3.17e-10 * std::exp(-5207.0/T_K);

    // 56) HD+ + H -> H+ + HD  = xk[9]
    xk[55] = xk[9];

    // 57) D+ + H2 -> H+ + HD  (GA08 39)
    xk[56] = (0.417 + 0.846*xlgT - 0.137*xlgT2) * 1.0e-9;

    // 58) HD + H -> H2 + D  (GA08 40)
    xk[57] = (T_K < 200.0)
           ? 5.25e-11*std::exp(-4430.0/T_K)
           : 5.25e-11*std::exp(-4430.0/T_K + 173900.0/(T_K*T_K));

    // 59) HD + H+ -> H2 + D+  (GA08 41)
    xk[58] = 1.1e-9 * std::exp(-488.0/T_K);

    // 60) D + H+ -> HD+ + ph.  (GA08 42)
    xk[59] = 3.9e-19 * std::pow(T300,1.8) * std::exp(20.0/T_K);

    // 61) D+ + H -> HD+ + ph.  (GA08 43)
    xk[60] = 3.9e-19 * std::pow(T300,1.8) * std::exp(20.0/T_K);

    // 62) HD+ + e -> H + D  (GA08 44)
    xk[61] = 7.2e-8 * std::pow(T_K,-0.5);

    // 63) D + e -> D- + ph.  = xk[6]
    xk[62] = xk[6];

    // 64) D+ + D- -> 2D  = xk[14]
    xk[63] = xk[14];

    // 65) H+ + D- -> D + H  = xk[14]
    xk[64] = xk[14];

    // 66) H- + D -> H + D-  (GA08 53)
    xk[65] = 6.4e-9 * std::pow(T300,0.41);

    // 67) D- + H -> D + H-  (GA08 52)
    xk[66] = 6.4e-9 * std::pow(T300,0.41);

    // 68) D- + H -> HD + e  (GA08 55)
    xk[67] = 0.5 * xk[7];

    // 69) D + e -> D+ + 2e  = xk[0]
    xk[68] = xk[0];

    // 70) He+ + D -> D+ + He  (GA08 46)
    xk[69] = 1.1e-15 * std::pow(T300,0.25);

    // 71) He + D+ -> D + He+  (GA08 47)
    xk[70] = (T_K < 10000.0)
           ? 1.85e-9*std::pow(T_K,-0.75)*std::exp(-127500.0/T_K)
           : 5.9e-37*std::pow(T_K,4.74);

    // 72) H2+ + D -> HD+ + H  (GA08 48)
    xk[71] = 1.07e-9 * std::pow(T300,0.062) * std::exp(-T_K/41400.0);

    // 73) HD+ + D -> HD + D+  = xk[9]
    xk[72] = xk[9];

    // 74) HD+ + H -> H2+ + D  (GA08 50)
    xk[73] = 1.0e-9 * std::exp(-154.0/T_K);

    // 75) HD + e -> H + D-  (GA08 57)
    double xk75v0 = 1.35e-9*std::pow(T_K,-1.27)*std::exp(-43000.0/T_K);
    xk[74] = xk75v0;   // stored; zeroed below

    // 76) HD + e -> D + H-  (GA08 58)
    xk[75] = 1.35e-9 * std::pow(T_K,-1.27) * std::exp(-43000.0/T_K);

    // 77) H+ + D- -> HD+ + e  (GA08 60)
    xk[76] = 1.1e-9 * std::pow(T300,-0.4);

    // 78) D+ + H- -> HD+ + e  (GA08 61)
    xk[77] = 1.1e-9 * std::pow(T300,-0.4);

    // 79) D- + e -> D + 2e  = xk[13]
    xk[78] = xk[13];

    // 80) D- + H -> D + H + e  = xk[37]
    xk[79] = xk[37];

    // 81) D- + He -> D + He + e  (GA08 65)
    xk[80] = 1.5e-17 * T_K*T_K * std::exp(-19870.0/T_K);

    // 82) D+ + H- -> D + H  = xk[14]
    xk[81] = xk[14];

    // 83) H2+ + D- -> H2 + D  (GA08 69)
    xk[82] = 1.7e-7 * std::pow(T300,-0.5);

    // 84) H2+ + D- -> 2H + D  (GA08 70)
    xk[83] = 1.7e-7 * std::pow(T300,-0.5);

    // 85) HD+ + H- -> HD + H  (GA08 71)
    xk[84] = 1.5e-7 * std::pow(T300,-0.5);

    // 86) HD+ + H- -> D + 2H  (GA08 72)
    xk[85] = 1.5e-7 * std::pow(T300,-0.5);

    // 87) HD+ + D- -> HD + D  (GA08 73)
    xk[86] = 1.9e-7 * std::pow(T300,-0.5);

    // 88) HD+ + D- -> 2D + H  (GA08 74)
    xk[87] = 1.9e-7 * std::pow(T300,-0.5);

    // 89) He+ + D- -> He + D  (GA08 79)
    xk[88] = 3.03e-7 * std::pow(T300,-0.52) * std::exp(T_K/22400.0);

    // 90) D + H2+ -> H2 + D+  = xk[9]
    xk[89] = xk[9];

    // 91) H2+ + D -> HD + H+  (GA08 82)
    xk[90] = 1.0e-9;

    // 92) HD+ + H -> H2 + D+  (GA08 83)
    xk[91] = 1.0e-9;

    // 93) H2 + D+ -> H2+ + D  = xk[10] (but xk[10]=0 from zero-forcing)
    xk[92] = xk[10];

    // 94) H2 + D+ -> HD+ + H  (GA08 91)
    xk[93] = (1.04e-9 + 9.52e-9*(T_K/10000.0)
              - 1.81e-9*std::pow(T_K/10000.0,2)) * std::exp(-21000.0/T_K);

    // 95) HD + H+ -> HD+ + H  = xk[10] (=0)
    xk[94] = xk[10];

    // 96) HD + H+ -> H2+ + D  (GA08 93)
    xk[95] = 1.0e-9 * std::exp(-21600.0/T_K);

    // 97) HD + D+ -> HD+ + D  = xk[10] (=0)
    xk[96] = xk[10];

    // 98) HD + He+ -> HD+ + He  = xk[29]
    xk[97] = xk[29];

    // 99) HD + He+ -> He + H+ + D  (GA08 102)
    xk[98] = 1.85e-14 * std::exp(-35.0/T_K);

    // 100) HD + He+ -> He + H + D+  (GA08 103)
    xk[99] = 1.85e-14 * std::exp(-35.0/T_K);

    // -----------------------------------------------------------------------
    //  D zero-forcing (L811–825)
    // -----------------------------------------------------------------------
    xk[52] = 0.0;   // 53
    xk[57] = 0.0;   // 58
    xk[58] = 0.0;   // 59
    xk[66] = 0.0;   // 67
    xk[70] = 0.0;   // 71
    xk[73] = 0.0;   // 74
    xk[74] = 0.0;   // 75
    xk[92] = 0.0;   // 93
    xk[93] = 0.0;   // 94
    xk[94] = 0.0;   // 95
    xk[95] = 0.0;   // 96
    xk[96] = 0.0;   // 97

    // -----------------------------------------------------------------------
    //  Metal reactions  xk[100..542]  (UMIST 2012)
    // -----------------------------------------------------------------------
    xk[100]=1.31e-10*std::exp(-80.0/T_K);
    xk[101]=6.00e-9*std::exp(-40200.0/T_K);
    xk[102]=2.2e-10;
    xk[103]=1.00e-10*std::exp(-7600.0/T_K);
    xk[104]=5.94e-13*std::pow(T300,3.0)*std::exp(-4045.0/T_K);
    xk[105]=6.99e-14*std::pow(T300,2.80)*std::exp(-1950.0/T_K);
    xk[106]=6.00e-9*std::exp(-50900.0/T_K);
    xk[107]=1.59e-11*std::pow(T300,1.2)*std::exp(-9610.0/T_K);
    xk[108]=5.80e-9*std::exp(-52900.0/T_K);
    xk[109]=4.67e-10*std::pow(T300,0.50)*std::exp(-30450.0/T_K);
    xk[110]=1.10e-10*std::pow(T300,0.50)*std::exp(-77700.0/T_K);
    xk[111]=4.85e-12*std::pow(T300,1.9)*std::exp(-1379.0/T_K);
    xk[112]=2.61e-10*std::exp(-8156.0/T_K);
    xk[113]=6.00e-9*std::exp(-52300.0/T_K);
    xk[114]=5.00e-11*std::exp(-866.0/T_K);
    xk[115]=2.06e-11*std::pow(T300,0.84)*std::exp(-277.0/T_K);
    xk[116]=1.66e-10*std::exp(-413.0/T_K);
    xk[117]=8.00e-11*std::exp(-4000.0/T_K);
    xk[118]=3.38e-10*std::exp(-13163.0/T_K);
    xk[119]=6.64e-10*std::exp(-11700.0/T_K);
    xk[120]=2.25e-11*std::pow(T300,0.50)*std::exp(-14800.0/T_K);
    xk[121]=2.94e-11*std::pow(T300,0.50)*std::exp(-58025.0/T_K);
    xk[122]=3.14e-13*std::pow(T300,2.7)*std::exp(-3150.0/T_K);
    xk[123]=2.52e-11*std::exp(-2381.0/T_K);
    xk[124]=4.98e-10*std::exp(-6000.0/T_K);
    xk[125]=5.01e-11;
    xk[126]=2.29e-12*std::pow(T300,2.2)*std::exp(-3820.0/T_K);
    xk[127]=1.85e-11*std::pow(T300,0.95)*std::exp(-8571.0/T_K);
    xk[128]=1.07e-11*std::pow(T300,1.17)*std::exp(-1242.0/T_K);
    xk[129]=8.54e-14*std::pow(T300,3.25)*std::exp(-1200.0/T_K);
    xk[130]=2.46e-11*std::exp(-26567.0/T_K);
    xk[131]=6.86e-10*std::pow(T300,0.26)*std::exp(-224.3/T_K);
    xk[132]=5.46e-10*std::exp(-1943.0/T_K);
    xk[133]=6.00e-9*std::exp(-40200.0/T_K);
    xk[134]=5.18e-11*std::pow(T300,0.17)*std::exp(-6400.0/T_K);
    xk[135]=6.86e-14*std::pow(T300,2.74)*std::exp(-4740.0/T_K);
    xk[136]=2.05e-12*std::pow(T300,1.52)*std::exp(-1736.0/T_K);
    xk[137]=6.00e-9*std::exp(-50900.0/T_K);
    xk[138]=5.80e-9*std::exp(-52900.0/T_K);
    xk[139]=2.40e-10*std::exp(-28500.0/T_K);
    xk[140]=3.16e-10*std::exp(-21890.0/T_K);
    xk[141]=6.00e-9*std::exp(-52300.0/T_K);
    xk[142]=4.38e-12*std::exp(-10751.0/T_K);
    xk[143]=0.0;   // 144: HD+grain->HD (grain surface)
    xk[144]=9.30e-10*std::exp(-100.0/T_K);
    xk[145]=1.00e-10*std::exp(-4640.0/T_K);
    xk[146]=0.0;   // 147: CH+CH4->CH2+CH3 (commented out)
    xk[147]=1.44e-11*std::pow(T300,0.50)*std::exp(-5000.0/T_K);
    xk[148]=2.87e-12*std::pow(T300,0.70)*std::exp(-500.0/T_K);
    xk[149]=9.21e-12*std::pow(T300,0.70)*std::exp(-2000.0/T_K);
    // 151) CH + O2 -> HCO + O  (piece-wise)
    if (T_K <= 10.0)
        xk[150] = 7.6e-12 * std::pow(10.0/300.0,-0.48);
    else if (T_K >= 300.0)
        xk[150] = 7.6e-12;
    else
        xk[150] = 7.6e-12 * std::pow(T300,-0.48);
    xk[151]=2.94e-13*std::pow(T300,0.50)*std::exp(-7550.0/T_K);
    xk[152]=1.44e-11*std::pow(T300,0.50)*std::exp(-3000.0/T_K);
    xk[153]=2.94e-13*std::pow(T300,0.50)*std::exp(-3000.0/T_K);
    xk[154]=4.00e-10*std::exp(-5000.0/T_K);
    xk[155]=7.13e-12*std::exp(-5050.0/T_K);
    xk[156]=3.00e-11;
    xk[157]=3.30e-13*std::exp(-3270.0/T_K);
    xk[158]=1.00e-9*std::exp(-7080.0/T_K);
    xk[159]=1.34e-15*std::pow(T300,5.05)*std::exp(-1636.0/T_K);
    xk[160]=7.00e-10*std::exp(-10560.0/T_K);
    xk[161]=1.44e-11*std::pow(T300,0.50)*std::exp(-3000.0/T_K);
    xk[162]=1.44e-11*std::pow(T300,0.50)*std::exp(-3000.0/T_K);
    xk[163]=3.00e-11;
    xk[164]=1.20e-10*std::exp(-1400.0/T_K);
    xk[165]=3.77e-13*std::pow(T300,2.42)*std::exp(-1162.0/T_K);
    xk[166]=1.65e-12*std::pow(T300,1.14)*std::exp(-50.0/T_K);
    xk[167]=2.81e-13*std::exp(-176.0/T_K);
    xk[168]=5.26e-12*std::exp(-307.0/T_K);
    xk[169]=2.30e-15*std::pow(T300,3.47)*std::exp(-6681.0/T_K);
    xk[170]=5.99e-12*std::exp(-24075.0/T_K);
    xk[171]=5.60e-10*std::exp(-12160.0/T_K);
    xk[172]=4.10e-11*std::exp(-750.0/T_K);
    xk[173]=3.65e-11*std::pow(T300,-3.3)*std::exp(-1443.0/T_K);
    xk[174]=5.30e-12*std::exp(-34975.0/T_K);
    xk[175]=5.64e-13*std::exp(-4500.0/T_K);
    xk[176]=6.70e-11*std::exp(-28640.0/T_K);
    xk[177]=4.64e-12*std::pow(T300,0.7)*std::exp(25.6/T_K);
    xk[178]=6.00e-12;
    xk[179]=4.65e-11*std::exp(-16500.0/T_K);
    xk[180]=5.00e-11;
    xk[181]=3.30e-12*std::exp(-5870.0/T_K);
    xk[182]=1.0e-12;
    xk[183]=1.50e-10;
    xk[184]=1.00e-17;
    xk[185]=4.36e-18*std::pow(T300,0.35)*std::exp(-161.3/T_K);
    xk[186]=4.69e-19*std::pow(T300,1.52)*std::exp(50.50/T_K);
    xk[187]=6.59e-11;
    xk[188]=1.00e-10;
    xk[189]=1.00e-10;
    xk[190]=5.56e-11*std::pow(T300,0.41)*std::exp(26.9/T_K);
    xk[191]=2.00e-10;
    xk[192]=9.90e-19*std::pow(T300,-0.38);
    xk[193]=4.90e-20*std::pow(T300,1.58);
    // 195) O + CH -> CO + H  (piece-wise 2000K)
    xk[194] = (T_K < 2000.0)
            ? 6.02e-11*std::pow(T300,0.10)*std::exp(4.50/T_K)
            : 1.02e-10*std::exp(-914.0/T_K);
    xk[195]=1.09e-11*std::pow(T300,-2.19)*std::exp(-165.1/T_K);
    xk[196]=1.33e-10;
    xk[197]=1.30e-10;
    // 199) O + OH -> O2 + H  (completely commented out → 0)
    // xk[198] = 0.0;  // already zero from fill
    xk[199]=2.00e-10*std::pow(T300,-0.12);
    xk[200]=5.00e-11;
    xk[201]=5.00e-11;
    // 203) O + O2H -> OH + O2  (piece-wise 200K)
    if (T_K > 200.0)
        xk[202] = 5.76e-11*std::pow(T300,-0.3)*std::exp(-7.5/T_K);
    else
        xk[202] = 5.76e-11*std::pow(200.0/300.0,-0.3)*std::exp(-7.5/200.0);
    xk[203]=1.70e-10;
    // 205) OH + H2CO -> H2O + HCO  (piece-wise 200K)
    if (T_K > 200.0)
        xk[204] = 7.76e-12*std::pow(T300,0.82)*std::exp(30.6/T_K);
    else
        xk[204] = 7.76e-12*std::pow(200.0/300.0,0.82)*std::exp(30.6/200.0);
    // 206) OH + O2H -> H2O + O2  (piece-wise 200K)
    if (T_K > 200.0)
        xk[205] = 8.58e-11*std::pow(T300,-0.56)*std::exp(-14.8/T_K);
    else
        xk[205] = 8.58e-11*std::pow(200.0/300.0,-0.56)*std::exp(-14.8/200.0);
    xk[206]=3.00e-11;
    xk[207]=1.90e-9*std::pow(T300,-0.5);
    xk[208]=1.40e-9;
    xk[209]=1.40e-9;
    xk[210]=3.40e-9;
    xk[211]=2.30e-9;
    xk[212]=1.50e-9;
    xk[213]=2.10e-9*std::pow(T300,-0.5);
    xk[214]=6.90e-9*std::pow(T300,-0.5);
    xk[215]=3.10e-9;
    xk[216]=9.40e-10*std::pow(T300,-0.5);
    xk[217]=9.40e-10*std::pow(T300,-0.5);
    xk[218]=9.40e-10*std::pow(T300,-0.5);
    xk[219]=2.96e-9*std::pow(T300,-0.5);
    xk[220]=3.57e-9*std::pow(T300,-0.5);
    xk[221]=2.00e-9;
    xk[222]=3.50e-9;
    xk[223]=1.00e-9;
    xk[224]=1.00e-9;
    xk[225]=1.00e-10;
    xk[226]=1.00e-9;
    xk[227]=1.00e-9;
    xk[228]=1.00e-10;
    xk[229]=2.00e-11;
    xk[230]=1.00e-9;
    xk[231]=2.40e-9;
    xk[232]=1.50e-9;
    xk[233]=7.10e-10*std::pow(T300,-0.5);
    xk[234]=7.10e-10*std::pow(T300,-0.5);
    xk[235]=1.00e-9;
    xk[236]=1.00e-9;
    xk[237]=1.40e-9;
    xk[238]=1.14e-10;
    xk[239]=2.30e-9;
    xk[240]=7.60e-10*std::pow(T300,-0.5);
    xk[241]=7.60e-10*std::pow(T300,-0.5);
    xk[242]=3.90e-9*std::pow(T300,-0.5);
    xk[243]=3.40e-9*std::pow(T300,-0.5);
    xk[244]=1.10e-9;
    xk[245]=6.44e-10;
    xk[246]=2.16e-9;
    xk[247]=1.00e-9*std::pow(T300,-0.5);
    xk[248]=1.00e-9*std::pow(T300,-0.5);
    xk[249]=1.40e-9*std::pow(T300,-0.5);
    xk[250]=1.40e-9*std::pow(T300,-0.5);
    xk[251]=1.90e-9;
    xk[252]=8.00e-10;
    xk[253]=2.35e-9;
    xk[254]=2.00e-9;
    // 256) H3+ + O -> OH+ + H2  (piece-wise 400K)
    if (T_K < 400.0)
        xk[255] = 7.98e-10*std::pow(T300,-0.16)*std::exp(-1.4/T_K);
    else
        xk[255] = 7.98e-10*std::pow(400.0/300.0,-0.16)*std::exp(-1.4/400.0);
    xk[256]=1.20e-9*std::pow(T300,-0.5);
    xk[257]=1.70e-9;
    xk[258]=2.10e-9;
    xk[259]=2.40e-9;
    xk[260]=1.30e-9*std::pow(T300,-0.5);
    xk[261]=5.90e-9*std::pow(T300,-0.5);
    // 263) H3+ + CO -> HCO+ + H2  (piece-wise 400K)
    if (T_K < 400.0)
        xk[262] = 1.36e-9*std::pow(T300,-0.14)*std::exp(3.4/T_K);
    else
        xk[262] = 1.36e-9*std::pow(400.0/300.0,-0.14)*std::exp(3.4/400.0);
    xk[263]=1.70e-9*std::pow(T300,-0.5);
    xk[264]=6.30e-9*std::pow(T300,-0.5);
    xk[265]=2.00e-9;
    xk[266]=1.10e-9*std::pow(T300,-0.5);
    xk[267]=5.00e-10*std::pow(T300,-0.5);
    xk[268]=7.50e-10;
    xk[269]=7.50e-10;
    xk[270]=1.80e-9;
    xk[271]=4.80e-10;
    xk[272]=2.40e-10;
    xk[273]=9.50e-10;
    xk[274]=8.50e-11;
    xk[275]=5.10e-11;
    xk[276]=1.10e-9*std::pow(T300,-0.5);
    xk[277]=2.04e-10*std::pow(T300,-0.5);
    xk[278]=2.86e-10*std::pow(T300,-0.5);
    xk[279]=6.05e-11*std::pow(T300,-0.5);
    xk[280]=1.60e-9;
    xk[281]=5.00e-10;
    xk[282]=1.60e-9;
    xk[283]=4.90e-10*std::pow(T300,-0.5);
    xk[284]=4.90e-10*std::pow(T300,-0.5);
    xk[285]=3.00e-10*std::pow(T300,-0.5);
    xk[286]=1.88e-9*std::pow(T300,-0.5);
    xk[287]=1.14e-9*std::pow(T300,-0.5);
    xk[288]=1.10e-9;
    xk[289]=3.30e-11;
    xk[290]=1.10e-11;
    xk[291]=1.00e-10;
    xk[292]=8.70e-10;
    xk[293]=4.00e-11;
    xk[294]=1.70e-17;
    xk[295]=3.14e-18*std::pow(T300,-0.15)*std::exp(-68.0/T_K);
    xk[296]=7.51e-8*std::pow(T300,-0.50);
    xk[297]=2.00e-16*std::pow(T300,-1.30)*std::exp(-23.0/T_K);
    xk[298]=3.80e-10*std::pow(T300,-0.50);
    xk[299]=3.80e-10*std::pow(T300,-0.50);
    xk[300]=5.20e-10;
    xk[301]=7.70e-10*std::pow(T300,-0.50);
    xk[302]=9.00e-10*std::pow(T300,-0.50);
    xk[303]=4.80e-10*std::pow(T300,-0.50);
    xk[304]=4.80e-10*std::pow(T300,-0.50);
    xk[305]=2.34e-9*std::pow(T300,-0.50);
    xk[306]=7.80e-10*std::pow(T300,-0.50);
    xk[307]=7.80e-10*std::pow(T300,-0.50);
    xk[308]=3.42e-10;
    xk[309]=4.54e-10;
    xk[310]=1.10e-9;
    // 312) CH+ + H -> C+ + H2  (piece-wise 1000K)
    if (T_K < 1000.0)
        xk[311] = 9.06e-10*std::pow(T300,-0.37)*std::exp(-29.1/T_K);
    else
        xk[311] = 9.06e-10*std::pow(1000.0/300.0,-0.37)*std::exp(-29.1/1000.0);
    xk[312]=1.20e-9;
    xk[313]=3.50e-10;
    xk[314]=1.20e-9;
    xk[315]=7.40e-10*std::pow(T300,-0.50);
    xk[316]=7.50e-10*std::pow(T300,-0.50);
    xk[317]=2.90e-9*std::pow(T300,-0.50);
    xk[318]=5.80e-10*std::pow(T300,-0.50);
    xk[319]=5.80e-10*std::pow(T300,-0.50);
    xk[320]=0.0;   // 321: none
    xk[321]=4.60e-10*std::pow(T300,-0.50);
    xk[322]=4.60e-10*std::pow(T300,-0.50);
    xk[323]=9.60e-10*std::pow(T300,-0.50);
    xk[324]=9.60e-10*std::pow(T300,-0.50);
    xk[325]=9.60e-10*std::pow(T300,-0.50);
    xk[326]=9.70e-10;
    xk[327]=1.00e-11;
    xk[328]=1.00e-11;
    xk[329]=1.60e-9;
    xk[330]=7.50e-10;
    xk[331]=1.60e-9;
    xk[332]=1.20e-9*std::pow(T300,-0.50);
    xk[333]=4.50e-10*std::pow(T300,-0.50);
    xk[334]=2.81e-9*std::pow(T300,-0.50);
    xk[335]=9.10e-10;
    xk[336]=1.60e-9;
    xk[337]=4.00e-10;
    xk[338]=4.00e-11;
    xk[339]=3.92e-16*std::pow(T300,-2.29)*std::exp(-21.3/T_K);
    xk[340]=7.20e-10*std::pow(T300,-0.50);
    xk[341]=4.40e-10*std::pow(T300,-0.50);
    xk[342]=4.40e-10*std::pow(T300,-0.50);
    xk[343]=1.60e-9*std::pow(T300,-0.50);
    xk[344]=5.00e-12;
    xk[345]=5.66e-10*std::pow(T300,0.36)*std::exp(8.6/T_K);
    xk[346]=7.51e-8*std::pow(T300,-0.50);
    xk[347]=1.70e-9;
    xk[348]=3.50e-10*std::pow(T300,-0.50);
    xk[349]=3.50e-10*std::pow(T300,-0.50);
    xk[350]=9.70e-10;
    xk[351]=1.10e-10;
    xk[352]=8.90e-10;
    xk[353]=3.60e-10*std::pow(T300,-0.50);
    xk[354]=3.60e-10*std::pow(T300,-0.50);
    xk[355]=3.20e-9*std::pow(T300,-0.50);
    xk[356]=4.80e-10;
    xk[357]=4.80e-10;
    xk[358]=4.30e-10*std::pow(T300,-0.50);
    xk[359]=4.30e-10*std::pow(T300,-0.50);
    xk[360]=1.40e-9*std::pow(T300,-0.50);
    xk[361]=2.10e-9*std::pow(T300,-0.50);
    xk[362]=1.90e-11;
    xk[363]=9.40e-10;
    xk[364]=1.00e-11;
    xk[365]=1.00e-9;
    xk[366]=4.89e-11*std::pow(T300,-0.14)*std::exp(36.1/T_K);
    xk[367]=1.50e-9;
    xk[368]=2.60e-9*std::pow(T300,-0.50);
    xk[369]=1.40e-9;
    xk[370]=1.62e-9*std::pow(T300,-0.50);
    xk[371]=1.98e-9*std::pow(T300,-0.50);
    xk[372]=3.90e-10;
    xk[373]=1.20e-9;
    xk[374]=1.20e-9;
    xk[375]=7.10e-10;
    xk[376]=1.01e-9;
    xk[377]=3.50e-10*std::pow(T300,-0.50);
    xk[378]=3.50e-10*std::pow(T300,-0.50);
    xk[379]=4.80e-10;
    xk[380]=4.80e-10;
    xk[381]=1.31e-9;
    xk[382]=1.95e-10;
    xk[383]=7.00e-10*std::pow(T300,-0.50);
    xk[384]=1.59e-9*std::pow(T300,-0.50);
    xk[385]=1.30e-9*std::pow(T300,-0.50);
    xk[386]=4.80e-10;
    xk[387]=1.05e-9;
    xk[388]=2.80e-10*std::pow(T300,-0.50);
    xk[389]=2.80e-10*std::pow(T300,-0.50);
    xk[390]=2.80e-10*std::pow(T300,-0.50);
    xk[391]=1.12e-9*std::pow(T300,-0.50);
    xk[392]=7.44e-10*std::pow(T300,-0.50);
    xk[393]=5.90e-10;
    xk[394]=1.44e-9;
    xk[395]=1.50e-10;
    xk[396]=1.20e-9;
    xk[397]=2.20e-10;
    xk[398]=4.40e-12;
    xk[399]=6.90e-10*std::pow(T300,-0.50);
    xk[400]=9.60e-10;
    xk[401]=7.00e-10*std::pow(T300,-0.50);
    xk[402]=3.70e-9*std::pow(T300,-0.50);
    xk[403]=1.00e-9;
    xk[404]=8.50e-10*std::pow(T300,-0.50);
    xk[405]=4.50e-9*std::pow(T300,-0.50);
    if (T_K <= 10.0)
        xk[406] = 0.0;
    else
        xk[406] = 3.20e-11;
    xk[407]=1.10e-9;
    xk[408]=4.00e-11;
    xk[409]=6.40e-10;
    xk[410]=3.40e-10*std::pow(T300,-0.50);
    xk[411]=3.40e-10*std::pow(T300,-0.50);
    xk[412]=4.70e-10;
    xk[413]=4.70e-10;
    xk[414]=1.40e-9;
    xk[415]=6.90e-10*std::pow(T300,-0.50);
    xk[416]=2.10e-9*std::pow(T300,-0.50);
    xk[417]=4.70e-10;
    xk[418]=5.00e-10;
    xk[419]=2.80e-10*std::pow(T300,-0.50);
    xk[420]=2.80e-10*std::pow(T300,-0.50);
    xk[421]=2.80e-10*std::pow(T300,-0.50);
    xk[422]=1.41e-9*std::pow(T300,-0.50);
    xk[423]=6.62e-10*std::pow(T300,-0.50);
    xk[424]=4.60e-10;
    xk[425]=1.00e-11;
    xk[426]=7.51e-8*std::pow(T300,-0.50);
    xk[427]=7.51e-8*std::pow(T300,-0.50);
    xk[428]=6.80e-10*std::pow(T300,-0.50);
    xk[429]=9.40e-10;
    xk[430]=3.40e-9*std::pow(T300,-0.50);
    xk[431]=1.10e-10;
    xk[432]=3.10e-10;
    xk[433]=3.20e-10*std::pow(T300,-0.50);
    xk[434]=4.50e-10;
    xk[435]=6.50e-10*std::pow(T300,-0.50);
    xk[436]=3.80e-10*std::pow(T300,-0.50);
    xk[437]=8.00e-10;
    xk[438]=7.50e-10;
    xk[439]=1.10e-10;
    xk[440]=1.40e-10;
    xk[441]=7.50e-10;
    xk[442]=3.20e-10*std::pow(T300,-0.50);
    xk[443]=3.20e-10*std::pow(T300,-0.50);
    xk[444]=4.30e-10;
    xk[445]=4.30e-10;
    xk[446]=7.93e-10;
    xk[447]=4.55e-10;
    xk[448]=3.10e-10*std::pow(T300,-0.50);
    xk[449]=3.10e-10*std::pow(T300,-0.50);
    xk[450]=1.72e-9*std::pow(T300,-0.50);
    xk[451]=8.84e-10*std::pow(T300,-0.50);
    xk[452]=8.40e-10;
    xk[453]=7.40e-10*std::pow(T300,-0.50);
    xk[454]=1.65e-9*std::pow(T300,-0.50);
    xk[455]=1.35e-9*std::pow(T300,-0.50);
    xk[456]=1.20e-10;
    xk[457]=1.10e-9;
    xk[458]=3.76e-8*std::pow(T300,-0.50);
    xk[459]=6.30e-10*std::pow(T300,-0.50);
    xk[460]=8.60e-10;
    xk[461]=6.20e-10*std::pow(T300,-0.50);
    xk[462]=1.00e-9*std::pow(T300,-0.50);
    xk[463]=2.50e-9*std::pow(T300,-0.50);
    xk[464]=7.30e-10*std::pow(T300,-0.50);
    xk[465]=3.30e-9*std::pow(T300,-0.50);
    xk[466]=3.10e-10*std::pow(T300,-0.50);
    xk[467]=3.10e-10*std::pow(T300,-0.50);
    xk[468]=4.30e-10;
    xk[469]=4.30e-10;
    xk[470]=9.35e-11;
    xk[471]=2.60e-9*std::pow(T300,-0.50);
    xk[472]=3.60e-10*std::pow(T300,-0.50);
    xk[473]=3.60e-10*std::pow(T300,-0.50);
    xk[474]=3.20e-9*std::pow(T300,-0.50);
    xk[475]=7.70e-11;
    xk[476]=6.20e-10*std::pow(T300,-0.50);
    xk[477]=2.30e-10*std::pow(T300,-0.50);
    xk[478]=5.20e-11;
    xk[479]=5.20e-11;
    xk[480]=3.10e-10*std::pow(T300,-0.50);
    xk[481]=3.10e-10*std::pow(T300,-0.50);
    xk[482]=4.30e-10;
    xk[483]=4.30e-10;
    xk[484]=4.10e-10;
    xk[485]=4.10e-10;
    xk[486]=3.60e-10*std::pow(T300,-0.50);
    xk[487]=3.60e-10*std::pow(T300,-0.50);
    xk[488]=2.30e-10*std::pow(T300,-0.50);
    xk[489]=2.07e-9*std::pow(T300,-0.50);
    xk[490]=1.00e-9;
    xk[491]=6.20e-10;
    xk[492]=6.40e-10;
    xk[493]=6.20e-10*std::pow(T300,-0.50);
    xk[494]=8.50e-10;
    xk[495]=6.10e-10*std::pow(T300,-0.50);
    xk[496]=8.20e-10*std::pow(T300,-0.50);
    xk[497]=8.40e-10;
    xk[498]=7.10e-10*std::pow(T300,-0.50);
    xk[499]=9.80e-10*std::pow(T300,-0.50);
    xk[500]=1.10e-9;
    xk[501]=1.00e-9;
    xk[502]=1.00e-9;
    xk[503]=7.80e-10;
    xk[504]=2.30e-9*std::pow(T300,-0.50);
    xk[505]=7.80e-10;
    xk[506]=2.36e-12*std::pow(T300,-0.29)*std::exp(17.6/T_K);
    xk[507]=1.50e-7*std::pow(T300,-0.42);
    xk[508]=7.68e-8*std::pow(T300,-0.60);
    xk[509]=1.60e-7*std::pow(T300,-0.60);
    xk[510]=7.75e-8*std::pow(T300,-0.50);
    xk[511]=1.95e-7*std::pow(T300,-0.50);
    xk[512]=2.00e-7*std::pow(T300,-0.40);
    xk[513]=1.10e-10*std::pow(T300,-0.50);
    xk[514]=3.24e-12*std::pow(T300,-0.66);
    xk[515]=1.75e-7*std::pow(T300,-0.50);
    xk[516]=1.75e-7*std::pow(T300,-0.50);
    xk[517]=3.75e-8*std::pow(T300,-0.50);
    xk[518]=1.40e-8*std::pow(T300,-0.52);
    xk[519]=1.40e-8*std::pow(T300,-0.52);
    xk[520]=8.60e-8*std::pow(T300,-0.50);
    xk[521]=3.90e-8*std::pow(T300,-0.50);
    xk[522]=7.09e-8*std::pow(T300,-0.50);
    xk[523]=3.05e-7*std::pow(T300,-0.50);
    xk[524]=3.00e-7*std::pow(T300,-0.50);
    xk[525]=2.00e-7*std::pow(T300,-0.48);
    xk[526]=2.40e-7*std::pow(T300,-0.69);
    xk[527]=1.60e-7*std::pow(T300,-0.70);
    xk[528]=2.50e-7*std::pow(T300,-0.70);
    xk[529]=1.10e-10*std::pow(T300,-0.70);
    xk[530]=2.10e-7*std::pow(T300,-0.78);
    xk[531]=2.17e-7*std::pow(T300,-0.78);
    xk[532]=2.17e-7*std::pow(T300,-0.78);
    xk[533]=1.95e-7*std::pow(T300,-0.70);
    xk[534]=3.00e-7*std::pow(T300,-0.50);
    xk[535]=6.00e-8*std::pow(T300,-0.64);
    xk[536]=3.20e-7*std::pow(T300,-0.64);
    xk[537]=1.00e-17;
    xk[538]=3.27e-14*std::pow(T300,2.20)*std::exp(-2240.0/T_K);
    xk[539]=1.06e-9*std::pow(T300,-0.5);
    xk[540]=6.30e-15*std::pow(T300,0.75);  // 541: He+ + C -> C+ + He
    xk[541]=9.69e-10*std::pow(T300,-0.5);
    xk[542]=1.71e-9*std::pow(T300,-0.5);

    // -----------------------------------------------------------------------
    //  Additional reactions  xk[600..644]  (601..645, various UMIST)
    // -----------------------------------------------------------------------
    xk[600]=6.61e-11*std::exp(-51598.0/T_K);
    xk[601]=1.70e-11*std::exp(-1800.0/T_K);
    xk[602]=2.69e-12*std::exp(-23550.0/T_K);
    xk[603]=7.60e-12;
    xk[604]=0.0;    // 605: error duplicate, explicitly zero
    xk[605]=3.65e-11*std::pow(T300,-3.3)*std::exp(-1443.0/T_K);
    xk[606]=2.92e-11*std::pow(T300,-3.3)*std::exp(-1443.0/T_K);
    xk[607]=2.48e-10*std::pow(T300,-3.3)*std::exp(-1443.0/T_K);
    xk[608]=7.13e-12*std::exp(-5052.0/T_K);
    xk[609]=3.60e-11*std::exp(-202.0/T_K);
    xk[610]=1.70e-12;
    xk[611]=1.66e-12;
    xk[612]=1.00e-10;
    xk[613]=1.50e-11*std::exp(-4300.0/T_K);
    xk[614]=3.63e-11;
    xk[615]=7.60e-13;
    xk[616]=3.69e-11*std::pow(T300,-0.27)*std::exp(-12.9/T_K);
    // 621) H3+ + O -> H2O+ + H  (piece-wise 400K)
    if (T_K <= 400.0)
        xk[620] = 3.42e-10*std::pow(T300,-0.16)*std::exp(-1.4/T_K);
    else
        xk[620] = 3.42e-10*std::pow(400.0/300.0,-0.16)*std::exp(-1.4/400.0);
    xk[625]=4.90e-12*std::pow(T300,0.5)*std::exp(-4580.0/T_K);
    xk[631]=4.03e-7*std::pow(T300,-0.6);
    xk[632]=1.96e-7*std::pow(T300,-0.52);
    xk[633]=4.76e-8*std::pow(T300,-0.52);
    xk[634]=8.40e-9*std::pow(T300,-0.52);
    xk[635]=3.05e-7*std::pow(T300,-0.5);
    xk[636]=5.60e-9*std::pow(T300,-0.5);
    xk[637]=5.37e-8*std::pow(T300,-0.5);
    xk[639]=8.10e-7*std::pow(T300,-0.64);
    // xk[640]=0.0;  // 641: commented out
    xk[641]=5.26e-18*std::pow(T300,-5.22)*std::exp(-90.0/T_K);
    xk[642]=5.09e-18*std::pow(T300,-0.71)*std::exp(-11.6/T_K);
    xk[643]=4.01e-18*std::pow(T300,0.17)*std::exp(-101.5/T_K);
    // 645) C + O+ -> CO+ + ph.  (piece-wise 2000K)
    if (T_K > 2000.0)
        xk[644] = 5.0e-10*std::pow(T300,-3.70)*std::exp(-800.0/T_K);
    else
        xk[644] = 5.0e-10*std::pow(2000.0/300.0,-3.70)*std::exp(-800.0/2000.0);

    // -----------------------------------------------------------------------
    //  Zero-forcing list  (L2332–2387)
    // -----------------------------------------------------------------------
    xk[119]=0.0;  xk[122]=0.0;  xk[120]=0.0;  xk[132]=0.0;
    xk[103]=0.0;  xk[104]=0.0;  xk[107]=0.0;  xk[139]=0.0;
    xk[142]=0.0;  xk[127]=0.0;  xk[118]=0.0;  xk[169]=0.0;
    xk[130]=0.0;  xk[176]=0.0;  xk[179]=0.0;  xk[109]=0.0;
    xk[110]=0.0;  xk[121]=0.0;  xk[145]=0.0;  xk[158]=0.0;
    xk[160]=0.0;  xk[131]=0.0;  xk[395]=0.0;  xk[435]=0.0;
    xk[320]=0.0;  xk[461]=0.0;  xk[477]=0.0;  xk[144]=0.0;
    xk[406]=0.0;  xk[126]=0.0;  xk[600]=0.0;  xk[608]=0.0;
    xk[112]=0.0;  xk[625]=0.0;

    // -----------------------------------------------------------------------
    //  Li reactions  xk[800..829]  (Reactions 801..830)
    // -----------------------------------------------------------------------
    // Bovino+2011: 801-820
    xk[800] = 1.036e-11 / (std::sqrt(T_K/107.7)
                           * std::pow(1.0 + std::sqrt(T_K/107.7), 0.6612)
                           * std::pow(1.0 + std::sqrt(T_K/1.177e7), 1.3388));
    xk[801] = 6.3e-9 * std::pow(T_K,-0.5) * (1.0 + T_K/14000.0);
    xk[802] = 2.3e-6 * std::pow(T_K,-0.5);
    xk[803] = 6.1e-17 * std::pow(T_K,0.58) * std::exp(-T_K/1.72e4);
    xk[804] = 2.5e-40 * std::pow(T_K,7.9) * std::exp(-T_K/1210.0);
    xk[805] = 1.7e-13 * std::pow(T_K,-0.051) * std::exp(-T_K/282000.0);
    xk[806] = 4.0e-10;
    xk[807] = 4.0e-10;
    xk[808] = 1.0e-11 * std::exp(-67900.0/T_K);
    xk[809] = 1.0e-9;
    xk[810] = 2.0e-12 * T_K * std::exp(-T_K/1200.0);
    xk[811] = 1.4e-20 * std::pow(T_K,-0.9) * std::exp(-T_K/7000.0);
    xk[812] = 5.3e-14 * std::pow(T_K,-0.49);
    xk[813] = 1.0e-9;
    xk[814] = 3.9e-6 * std::pow(T_K,-0.70) * std::exp(-T_K/1200.0);
    xk[815] = 9.0e-10 * std::exp(-66400.0/T_K);
    xk[816] = 8.7e-10 * std::pow(T_K,0.040) * std::exp(T_K/5.92e8);
    xk[817] = 4.0e-20 * std::exp(-T_K/4065.0
                                  + std::pow(T_K/13193.0, 3.0));
    if (T_K < 500.0)
        xk[818] = 6.3e-10 * std::exp(-2553.0/T_K);
    else
        xk[818] = 7.2e-14 * std::pow(T_K,1.18) * std::exp(-1470.0/T_K);
    xk[819] = 2.9e-10*std::pow(T_K,0.59) - 2.6e-10*std::pow(T_K,0.6)*std::exp(-400.0/T_K);
    // Lepp+2002: 821-827
    xk[820] = 5.34e-8 * std::pow(T300,-1.23) * std::exp(-T_K/9.23e5);
    xk[821] = 4.83e-11 * std::pow(T300,-0.621) * std::exp(-T_K/1.67e6);
    xk[822] = 3.71e-7 * std::pow(T300,-0.51) * std::exp(T_K/4.41e4);
    xk[823] = 2.28e-7 * std::pow(T300,-0.51) * std::exp(T_K/4.41e4);
    xk[824] = 3.11e-8 * std::pow(T300,0.163) * std::exp(-6.27e4/T_K);
    xk[825] = 5.67e-12 * std::pow(T300,0.715) * std::exp(-8.77e5/T_K);
    xk[826] = 1.70e-12 * std::pow(T300,0.709) * std::exp(-1.42e6/T_K);
    // Mizusawa+2005: 828-829
    xk[827] = 2.5e-29 / T_K;
    xk[828] = 4.1e-30 / T_K;
    // 830) H2 + Li -> H2 + Li+ + e  (posit)
    xk[829] = 9.9e-9 * std::sqrt(T_K) * std::exp(tbl.delE[634] / (xk_B*T_K));

    // Li zero-forcing
    xk[808] = 0.0;   // 809
    xk[819] = 0.0;   // 820

    // -----------------------------------------------------------------------
    //  K, Na, Mg reactions  xk[700..729]  (Reactions 701..730)
    // -----------------------------------------------------------------------
    xk[700] = 3.0e-11 / std::sqrt(T_K);
    xk[701] = 9.9e-9 * std::sqrt(T_K) * std::exp(tbl.delE[576] / (xk_B*T_K));
    xk[702] = 2.76e-12 * std::pow(T300,-0.68);
    xk[703] = 1.20e-9;
    xk[704] = 7.51e-8 * std::pow(T300,-0.50);
    xk[705] = 2.1e-9;
    xk[706] = 1.1e-9;
    xk[707] = 3.5e-10;
    xk[708] = 3.4e-9;
    xk[709] = 7.1e-10;
    xk[710] = 6.2e-9;
    xk[711] = 3.1e-9;
    xk[712] = 2.6e-9;
    xk[713] = 2.6e-9;
    xk[714] = 2.6e-9;
    xk[715] = 1.0e-11;
    xk[716] = 9.9e-9 * std::sqrt(T_K) * std::exp(tbl.delE[591] / (xk_B*T_K));
    xk[717] = 2.78e-12 * std::pow(T300,-0.68);
    xk[718] = 1.1e-9;
    xk[719] = 7.51e-8 * std::pow(T300,-0.50);
    xk[720] = 1.0e-9;
    xk[721] = 1.1e-9;
    xk[722] = 3.6e-10;
    xk[723] = 3.5e-9;
    xk[724] = 1.4e-9;
    xk[725] = 1.2e-9;
    xk[726] = 2.2e-9;
    xk[727] = 0.0;    // 728: = 0 from Fortran
    xk[728] = 2.9e-9;
    xk[729] = 2.9e-9;

    // -----------------------------------------------------------------------
    //  Reverse reactions  (xk[1199..2399] = xk[num-1 + 1200])
    //  Loop over reaction table entries ire = 0..634 (Fortran 1..635)
    //  Uses detailed balance: lnK_rev = lnK_for + lnKeqb
    // -----------------------------------------------------------------------
    // Save lnKeqb for post-loop corrections
    std::array<double, 1200> lnKeqb{};

    const double lnC1_base = std::log(2.0*pi*xk_B*T_K / (h_P*h_P));

    for (int ire = 0; ire < 635; ++ire) {
        int num = tbl.rean[ire];   // Fortran 1-based reaction number
        int r1  = tbl.rea1[ire];
        int r2  = tbl.rea2[ire];
        int r3  = tbl.rea3[ire];
        int p1  = tbl.pro1[ire];
        int p2  = tbl.pro2[ire];
        int p3  = tbl.pro3[ire];
        int nr  = tbl.nrea[ire];
        int np  = tbl.npro[ire];
        double Cm = tbl.Cmass[ire];
        double dE = tbl.delE[ire];

        double xlnC1  = 1.5 * (nr - np) * lnC1_base;
        double xlnCm  = std::log(Cm);
        double xlnCpf = std::log(pf[r1]) + std::log(pf[r2]) + std::log(pf[r3])
                      - std::log(pf[p1]) - std::log(pf[p2]) - std::log(pf[p3]);

        double lnKeqb_num = xlnC1 + xlnCm + xlnCpf - dE/(xk_B*T_K);
        lnKeqb[num-1] = lnKeqb_num;   // store (0-based)

        int n = num - 1;   // 0-based index into xk
        if (p3 == 101) {
            // photodissociation reverse
            if (tau_cnt <= 1.0e-10) {
                xk[n + 1200] = 0.0;
            } else {
                double esc_fact = (tau_cnt <= 1.0e-3)
                                ? tau_cnt - 0.5*tau_cnt*tau_cnt
                                : 1.0 - std::exp(-tau_cnt);
                lnKeqb_num += std::log(esc_fact);
                lnKeqb[n] = lnKeqb_num;
                xk[n + 1200] = std::exp(std::log(xk[n]) + lnKeqb_num);
            }
        } else {
            xk[n + 1200] = std::exp(std::log(xk[n]) + lnKeqb_num);
        }
        if (xk[n] == 0.0) xk[n + 1200] = 0.0;
    }

    // -----------------------------------------------------------------------
    //  Special case: reaction 273  He+ + CH4 -> CH+ + H2 + He + H  (4-body)
    //  Overrides the value set in the loop above.
    // -----------------------------------------------------------------------
    {
        constexpr double Cm273 = 1.01e71;
        constexpr double dE273 = 7.80e-12;
        double lnC1_273  = -3.0 * lnC1_base;
        double lnCpf_273 = std::log(pf[9])  + std::log(pf[22])
                         - std::log(pf[1])  - std::log(pf[2])
                         - std::log(pf[8])  - std::log(pf[25]);
        double lnK273 = lnC1_273 + std::log(Cm273) + lnCpf_273 - dE273/(xk_B*T_K);
        xk[272 + 1200] = std::exp(std::log(xk[272]) + lnK273);
    }

    // -----------------------------------------------------------------------
    //  Post-loop LTE/low-density corrections for specific reverse rates
    // -----------------------------------------------------------------------
    // xk(8+N_react): H- + H -> H2 + e  reverse  [reaction 40: H2+e->H-+H]
    {
        double xk40LTE = std::exp(std::log(xk[7]) + lnKeqb[7]);  // xk(8)+lnKeqb(8)
        if (xk40v0 == 0.0 || xk40LTE == 0.0) {
            xk[7 + 1200] = 0.0;
        } else {
            double xlgk40 = (xcr/(1.0+xcr))*std::log10(xk40LTE)
                          + (1.0/(1.0+xcr))*std::log10(xk40v0);
            xk[7 + 1200] = std::pow(10.0, xlgk40);
        }
    }

    // xk(12+N_react): H2+e->2H+e reverse
    {
        double xlnk12rev = std::log(xk12LTE) + lnKeqb[11];
        xk[11 + 1200] = std::exp(xlnk12rev);
    }

    // xk(19+N_react): H2+H->3H reverse  (use xk13LTE as LTE branch of 13)
    {
        double xk13LTE_rev = std::exp(std::log(xk[18]) + lnKeqb[18]);
        if (xk13v0 == 0.0 || xk13LTE_rev == 0.0) {
            xk[18 + 1200] = 0.0;
        } else {
            double xlgk13r = (xcr/(1.0+xcr))*std::log10(xk13LTE_rev)
                           + (1.0/(1.0+xcr))*std::log10(xk13v0);
            xk[18 + 1200] = std::pow(10.0, xlgk13r);
        }
    }

    // xk(21+N_react): 2H2->2H+H2 reverse
    {
        double xlnk21rev = std::log(xk21LTE) + lnKeqb[20];
        xk[20 + 1200] = std::exp(xlnk21rev);
    }

    // xk(37+N_react): H2+He->2H+He reverse
    {
        double xlnk37rev = std::log(xk37LTE) + lnKeqb[36];
        xk[36 + 1200] = std::exp(xlnk37rev);
    }

    // xk(47+N_react): HD+e->H+D+e reverse
    {
        double xk47LTE_val = std::exp(std::log(xk[46]) + lnKeqb[46]);
        double xk47v0_val  = 5.09e-9*std::pow(T_K,0.128)*std::exp(-103258.0/T_K);
        if (xk47v0_val == 0.0 || xk47LTE_val == 0.0)
            xk[46 + 1200] = 0.0;
        else {
            double xlgk47r = (xcr_HD/(1.0+xcr_HD))*std::log10(xk47LTE_val)
                           + (1.0/(1.0+xcr_HD))*std::log10(xk47v0_val);
            xk[46 + 1200] = std::pow(10.0, xlgk47r);
        }
    }

    // xk(48+N_react): HD+He->H+D+He  reverse  (uses xk37LTE)
    {
        double xlnk48rev = std::log(xk37LTE) + lnKeqb[47];
        xk[47 + 1200] = std::exp(xlnk48rev);
    }

    // xk(49+N_react): HD+H2->H+D+H2  reverse  (uses xk21LTE)
    {
        double xlnk49rev = std::log(xk21LTE) + lnKeqb[48];
        xk[48 + 1200] = std::exp(xlnk49rev);
    }

    // xk(50+N_react): HD+H->2H+D  reverse
    // Fortran reassigns xk13LTE at line 2741 to exp(log(xk(19))+xlnKeqb(19))
    // before using it here.  Replicate that: use detailed-balance LTE from xk[18].
    {
        double xk13LTE_db = std::exp(std::log(xk[18]) + lnKeqb[18]);
        double xlnk50rev  = std::log(xk13LTE_db) + lnKeqb[49];
        xk[49 + 1200] = std::exp(xlnk50rev);
    }

    // xk(68+N_react): D-+H->HD+e  reverse  (xcr_HD blend)
    {
        double xk75LTE_val = std::exp(std::log(xk[67]) + lnKeqb[67]);
        if (xk75v0 == 0.0 || xk75LTE_val == 0.0) {
            xk[67 + 1200] = 0.0;
        } else {
            double xlgk75 = (xcr_HD/(1.0+xcr_HD))*std::log10(xk75LTE_val)
                          + (1.0/(1.0+xcr_HD))*std::log10(xk75v0);
            xk[67 + 1200] = std::pow(10.0, xlgk75);
        }
    }

    // xcr < 1: zero reverse of reactions 1 and 3
    if (xcr < 1.0) {
        xk[0 + 1200] = 0.0;
        xk[2 + 1200] = 0.0;
    }

    // T <= 300K: zero specific reverse rates
    if (T_K <= 300.0) {
        xk[101 + 1200] = 0.0;
        xk[106 + 1200] = 0.0;
        xk[108 + 1200] = 0.0;
        xk[113 + 1200] = 0.0;
        xk[133 + 1200] = 0.0;
        xk[137 + 1200] = 0.0;
        xk[138 + 1200] = 0.0;
        xk[140 + 1200] = 0.0;
        xk[141 + 1200] = 0.0;
        xk[174 + 1200] = 0.0;
        xk[181 + 1200] = 0.0;
        xk[229 + 1200] = 0.0;
        xk[539 + 1200] = 0.0;
        xk[602 + 1200] = 0.0;
    }

    // -----------------------------------------------------------------------
    //  Dust grain reactions
    //  react_grain_rates() -> xk[840..989] (Reactions 841..990)
    //  grain_coef_rates()  -> xk[990..1140] (Reactions 991..1141)
    // -----------------------------------------------------------------------
    {
        // Vaporization temperature check
        double T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol;
        detail::vaptemp(rho, T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol);
        double T_evap = std::max({T_ir, T_pyr, T_ol});

        std::array<double,150> xk_charge;
        std::array<double,151> xkgr;
        react_grain_rates(T_K, T_gr_K, xk_charge);
        grain_coef_rates(xnH, T_K, T_gr_K, Z_metal,
                         xJH2, xJH2O, xJtot, zeta, T_cr_desorp, xkgr);

        if (T_gr_K >= 1.025 * T_evap) {
            xk_charge.fill(0.0);
            xkgr.fill(0.0);
        }

        // Map charge rates: rean(672..821) -> xk[840..989]
        for (int i = 0; i < 150; ++i)
            xk[840 + i] = xk_charge[i];

        // Electron sticking coefficient = 0.6
        xk[872] *= 0.6;   // Fortran xk(873)
        xk[946] *= 0.6;   // Fortran xk(947)

        // Map surface rates: xkgr[0..150] -> xk[990..1140]
        for (int i = 0; i < 151; ++i)
            xk[990 + i] = xkgr[i];
    }

    // -----------------------------------------------------------------------
    //  CR reactions  (set after grain mapping, as in Fortran)
    //  xk[543..551]  (Reactions 544..552)
    //  xk[655..681]  (Reactions 656..682)
    // -----------------------------------------------------------------------
    const double zeta_fac = zeta / 1.36e-17;
    xk[543]=5.98e-18*zeta_fac;
    xk[544]=6.50e-18*zeta_fac;
    xk[545]=2.30e-17*zeta_fac;
    xk[546]=3.40e-17*zeta_fac;
    xk[547]=2.86e-19*zeta_fac;
    xk[548]=1.20e-17*zeta_fac;
    xk[549]=1.30e-18*zeta_fac;
    xk[550]=3.90e-21*zeta_fac;
    xk[551]=3.90e-17*zeta_fac;

    // CR-induced photoreactions (omega = 0.6)
    constexpr double omega = model::cr_photo_albedo;
    const double cr_ph = 1.3e-17 * zeta_fac / (1.0 - omega);
    xk[655]=cr_ph*255.0;
    xk[656]=cr_ph*365.0;
    xk[657]=cr_ph*88.0;
    xk[658]=cr_ph*250.0;
    xk[659]=cr_ph*250.0;
    xk[660]=cr_ph*250.0;
    xk[661]=cr_ph*250.0;
    xk[662]=cr_ph*250.0;
    xk[663]=cr_ph*1169.5;
    xk[664]=cr_ph*254.5;
    xk[665]=cr_ph*485.5;
    xk[666]=cr_ph*119.5;
    // 668) CO + CR ph. -> O + C  (T-dependent)
    xk[667]=1.3e-17*zeta_fac*std::pow(T300,1.17)*105.0/(1.0-omega);
    xk[668]=cr_ph*210.5;
    xk[669]=cr_ph*584.5;
    xk[670]=cr_ph*1329.5;
    xk[671]=cr_ph*375.5;
    xk[672]=cr_ph*58.5;
    xk[673]=cr_ph*750.0;
    xk[674]=cr_ph*854.0;
    xk[675]=cr_ph*250.0;
    xk[676]=cr_ph*0.2;
    xk[677]=cr_ph*0.2;
    xk[678]=cr_ph*1.4;
    xk[679]=cr_ph*375.0;
    xk[680]=cr_ph*8.5;
    xk[681]=cr_ph*66.5;
}

// ---------------------------------------------------------------------------
// compute_rates<89,1200>  — template specialization for metal_grain network
//
// Mirrors react_rat() in reaction_metal_grain.f (lines 2989-3429).
// Four loop categories:
//   Loop 1  ire=0..634  (Fortran 1..635):  standard reactions, fwd+rev
//   Loop 2  ire=635..670 (Fortran 636..671): CR first-order reactions, no rev
//   Loop 3  ire=671..820 (Fortran 672..821): gas-grain charge transfer, fwd only
//   Loop 4  ire=0..n_grain-1:               grain surface reactions, fwd only
//
// Requires:
//   - tbl arrays (rean/rea1/..) stored sequentially by file position (load_reaction_table)
//   - tbl grain surface arrays (grea1/grea2/gpro1/gpro2/gnrea) loaded from react_grain_surface.dat
// ---------------------------------------------------------------------------
template<>
inline void compute_rates<89, 1200>(
    const std::array<double, 2400>& xk,
    double xnH,
    const std::array<double, 89>& y,
    const ReactionTable<89, 1200>& tbl,
    std::array<double, 89>& r_f,
    std::array<double, 89*89>& dr_fdy,
    std::array<double, 2400>& var)
{
    r_f.fill(0.0);
    dr_fdy.fill(0.0);
    var.fill(0.0);

    // Extended species array: 1-based species indices, 0 and 101 are sentinels
    std::array<double, 102> y_ext;
    y_ext.fill(0.0);
    y_ext[0]   = 1.0;   // vacant slot sentinel
    y_ext[101] = 1.0;   // photon sentinel
    for (int i = 0; i < 89; ++i) y_ext[i+1] = y[i];

    std::array<double, 102> r_f_dum{};
    r_f_dum.fill(0.0);
    std::array<double, 102*102> dr_fdy_dum{};
    dr_fdy_dum.fill(0.0);

    auto dJ = [&](int i, int j) -> double& {
        return dr_fdy_dum[i * 102 + j];
    };

    // ─── Loop 1: standard reactions ire=0..634 (Fortran do 102 ire=1,635) ─────
    for (int ire = 0; ire < 635; ++ire) {
        int num = tbl.rean[ire];
        if (num < 1 || num > 1200) continue;
        int r1 = tbl.rea1[ire], r2 = tbl.rea2[ire], r3 = tbl.rea3[ire];
        int p1 = tbl.pro1[ire], p2 = tbl.pro2[ire], p3 = tbl.pro3[ire];
        int nr = tbl.nrea[ire], np = tbl.npro[ire];

        double nHf = std::pow(xnH, double(nr - 1));

        // Forward rate: special cases for reactions 23 and 144
        double rate_fwd;
        if (num == 23) {
            rate_fwd = xk[22] * y_ext[1] * xnH;          // 2H+grain→H2
        } else if (num == 144) {
            rate_fwd = xk[143] * y_ext[12] * xnH;         // 2D+grain→HD
        } else {
            rate_fwd = xk[num-1] * y_ext[r1] * y_ext[r2] * y_ext[r3] * nHf;
        }

        // Reverse rate: special case for reaction 273 (4-body: extra y_H factor)
        double rate_rev, nHr;
        if (num == 273) {
            nHr = std::pow(xnH, double(np));
            rate_rev = xk[num-1+1200]
                     * y_ext[p1] * y_ext[p2] * y_ext[p3] * y_ext[1] * nHr;
        } else {
            nHr = std::pow(xnH, double(np - 1));
            rate_rev = xk[num-1+1200] * y_ext[p1] * y_ext[p2] * y_ext[p3] * nHr;
        }

        double rate = rate_fwd - rate_rev;

        // Species update (r_f_dum uses 1-based sentinel indices)
        r_f_dum[r1] -= rate;  r_f_dum[r2] -= rate;  r_f_dum[r3] -= rate;
        r_f_dum[p1] += rate;  r_f_dum[p2] += rate;  r_f_dum[p3] += rate;
        if (num == 273) r_f_dum[1] += rate;  // 4th product H (species 1)

        var[num-1]        = rate_fwd;
        var[num-1 + 1200] = rate_rev;

        // ── Jacobian: forward contributions ──────────────────────────────────
        if (num == 23) {
            // rate_fwd = xk[22]*y_ext[1]*xnH; only depends on y_ext[1] (H)
            double d1 = xk[22] * xnH;
            dJ(r1,1) -= d1;  dJ(r2,1) -= d1;  dJ(r3,1) -= d1;
            dJ(p1,1) += d1;  dJ(p2,1) += d1;  dJ(p3,1) += d1;
        } else if (num == 144) {
            // rate_fwd = xk[143]*y_ext[12]*xnH; only depends on y_ext[12] (D)
            double d12 = xk[143] * xnH;
            dJ(r1,12) -= d12;  dJ(r2,12) -= d12;  dJ(r3,12) -= d12;
            dJ(p1,12) += d12;  dJ(p2,12) += d12;  dJ(p3,12) += d12;
        } else {
            double fJ_r1 = xk[num-1] * y_ext[r2] * y_ext[r3] * nHf;
            dJ(r1,r1) -= fJ_r1;  dJ(r2,r1) -= fJ_r1;  dJ(r3,r1) -= fJ_r1;
            dJ(p1,r1) += fJ_r1;  dJ(p2,r1) += fJ_r1;  dJ(p3,r1) += fJ_r1;
            if (num == 273) dJ(1,r1) += fJ_r1;

            double fJ_r2 = xk[num-1] * y_ext[r1] * y_ext[r3] * nHf;
            dJ(r1,r2) -= fJ_r2;  dJ(r2,r2) -= fJ_r2;  dJ(r3,r2) -= fJ_r2;
            dJ(p1,r2) += fJ_r2;  dJ(p2,r2) += fJ_r2;  dJ(p3,r2) += fJ_r2;
            if (num == 273) dJ(1,r2) += fJ_r2;

            double fJ_r3 = xk[num-1] * y_ext[r1] * y_ext[r2] * nHf;
            dJ(r1,r3) -= fJ_r3;  dJ(r2,r3) -= fJ_r3;  dJ(r3,r3) -= fJ_r3;
            dJ(p1,r3) += fJ_r3;  dJ(p2,r3) += fJ_r3;  dJ(p3,r3) += fJ_r3;
            if (num == 273) dJ(1,r3) += fJ_r3;
        }

        // ── Jacobian: reverse contributions ──────────────────────────────────
        if (num == 273) {
            // rate_rev includes extra y_ext[1] (H) and uses nHr = xnH^np
            double rJ_p1 = xk[num-1+1200] * y_ext[p2] * y_ext[p3] * y_ext[1] * nHr;
            dJ(r1,p1) += rJ_p1;  dJ(r2,p1) += rJ_p1;  dJ(r3,p1) += rJ_p1;
            dJ(p1,p1) -= rJ_p1;  dJ(p2,p1) -= rJ_p1;  dJ(p3,p1) -= rJ_p1;
            dJ(1,p1)  -= rJ_p1;

            double rJ_p2 = xk[num-1+1200] * y_ext[p1] * y_ext[p3] * y_ext[1] * nHr;
            dJ(r1,p2) += rJ_p2;  dJ(r2,p2) += rJ_p2;  dJ(r3,p2) += rJ_p2;
            dJ(p1,p2) -= rJ_p2;  dJ(p2,p2) -= rJ_p2;  dJ(p3,p2) -= rJ_p2;
            dJ(1,p2)  -= rJ_p2;

            double rJ_p3 = xk[num-1+1200] * y_ext[p1] * y_ext[p2] * y_ext[1] * nHr;
            dJ(r1,p3) += rJ_p3;  dJ(r2,p3) += rJ_p3;  dJ(r3,p3) += rJ_p3;
            dJ(p1,p3) -= rJ_p3;  dJ(p2,p3) -= rJ_p3;  dJ(p3,p3) -= rJ_p3;
            dJ(1,p3)  -= rJ_p3;

            // d/d(y_ext[1]): rate_rev ∝ y_ext[1]
            double rJ_1 = xk[num-1+1200] * y_ext[p1] * y_ext[p2] * y_ext[p3] * nHr;
            dJ(r1,1) += rJ_1;  dJ(r2,1) += rJ_1;  dJ(r3,1) += rJ_1;
            dJ(p1,1) -= rJ_1;  dJ(p2,1) -= rJ_1;  dJ(p3,1) -= rJ_1;
            dJ(1,1)  -= rJ_1;
        } else {
            double rJ_p1 = xk[num-1+1200] * y_ext[p2] * y_ext[p3] * nHr;
            dJ(r1,p1) += rJ_p1;  dJ(r2,p1) += rJ_p1;  dJ(r3,p1) += rJ_p1;
            dJ(p1,p1) -= rJ_p1;  dJ(p2,p1) -= rJ_p1;  dJ(p3,p1) -= rJ_p1;

            double rJ_p2 = xk[num-1+1200] * y_ext[p1] * y_ext[p3] * nHr;
            dJ(r1,p2) += rJ_p2;  dJ(r2,p2) += rJ_p2;  dJ(r3,p2) += rJ_p2;
            dJ(p1,p2) -= rJ_p2;  dJ(p2,p2) -= rJ_p2;  dJ(p3,p2) -= rJ_p2;

            double rJ_p3 = xk[num-1+1200] * y_ext[p1] * y_ext[p2] * nHr;
            dJ(r1,p3) += rJ_p3;  dJ(r2,p3) += rJ_p3;  dJ(r3,p3) += rJ_p3;
            dJ(p1,p3) -= rJ_p3;  dJ(p2,p3) -= rJ_p3;  dJ(p3,p3) -= rJ_p3;
        }
    }

    // ─── Loop 2: CR first-order reactions ire=635..670 (Fortran do 103 ire=636,671) ─
    for (int ire = 635; ire < 671; ++ire) {
        int num = tbl.rean[ire];
        if (num < 1 || num > 1200) continue;
        int r2 = tbl.rea2[ire];   // r1 = 0 (CR particle, sentinel)
        int p1 = tbl.pro1[ire], p2 = tbl.pro2[ire], p3 = tbl.pro3[ire];

        double rate = xk[num-1] * y_ext[r2];

        r_f_dum[r2] -= rate;
        r_f_dum[p1] += rate;  r_f_dum[p2] += rate;  r_f_dum[p3] += rate;

        var[num-1] = rate;

        double d2 = xk[num-1];
        dJ(r2,r2) -= d2;
        dJ(p1,r2) += d2;  dJ(p2,r2) += d2;  dJ(p3,r2) += d2;
    }

    // ─── Loop 3: gas-grain charge transfer ire=671..820 (Fortran do 104 ire=672,821) ─
    for (int ire = 671; ire < 821; ++ire) {
        int num = tbl.rean[ire];
        if (num < 1 || num > 1200) continue;
        int r1 = tbl.rea1[ire], r2 = tbl.rea2[ire], r3 = tbl.rea3[ire];
        int p1 = tbl.pro1[ire], p2 = tbl.pro2[ire], p3 = tbl.pro3[ire];
        int nr = tbl.nrea[ire];

        double nHf = std::pow(xnH, double(nr - 1));
        double rate = xk[num-1] * y_ext[r1] * y_ext[r2] * y_ext[r3] * nHf;

        r_f_dum[r1] -= rate;  r_f_dum[r2] -= rate;  r_f_dum[r3] -= rate;
        r_f_dum[p1] += rate;  r_f_dum[p2] += rate;  r_f_dum[p3] += rate;

        var[num-1] = rate;

        double fJ_r1 = xk[num-1] * y_ext[r2] * y_ext[r3] * nHf;
        dJ(r1,r1) -= fJ_r1;  dJ(r2,r1) -= fJ_r1;  dJ(r3,r1) -= fJ_r1;
        dJ(p1,r1) += fJ_r1;  dJ(p2,r1) += fJ_r1;  dJ(p3,r1) += fJ_r1;

        double fJ_r2 = xk[num-1] * y_ext[r1] * y_ext[r3] * nHf;
        dJ(r1,r2) -= fJ_r2;  dJ(r2,r2) -= fJ_r2;  dJ(r3,r2) -= fJ_r2;
        dJ(p1,r2) += fJ_r2;  dJ(p2,r2) += fJ_r2;  dJ(p3,r2) += fJ_r2;

        double fJ_r3 = xk[num-1] * y_ext[r1] * y_ext[r2] * nHf;
        dJ(r1,r3) -= fJ_r3;  dJ(r2,r3) -= fJ_r3;  dJ(r3,r3) -= fJ_r3;
        dJ(p1,r3) += fJ_r3;  dJ(p2,r3) += fJ_r3;  dJ(p3,r3) += fJ_r3;
    }

    // ─── Loop 4: grain surface reactions (Fortran do 106 ire=1,151) ──────────
    // num increments from tbl.rean[820] (last charge reaction) + 1 = 991, 992, ..., 1141
    {
        const int num_start = (tbl.n_loaded >= 821 && tbl.rean[820] > 0)
                              ? tbl.rean[820]   // = 990 from data file
                              : 990;
        for (int ire = 0; ire < tbl.n_grain; ++ire) {
            int num = num_start + 1 + ire;      // 991, 992, ..., 1141
            if (num < 1 || num > 1200) continue;
            int r1 = tbl.grea1[ire], r2 = tbl.grea2[ire];
            int p1 = tbl.gpro1[ire], p2 = tbl.gpro2[ire];
            int nr = tbl.gnrea[ire];   // exponent directly: rate ∝ xnH^nr (not nr-1)

            double nHr = std::pow(xnH, double(nr));
            double rate = xk[num-1] * y_ext[r1] * y_ext[r2] * nHr;

            r_f_dum[r1] -= rate;  r_f_dum[r2] -= rate;
            r_f_dum[p1] += rate;  r_f_dum[p2] += rate;

            var[num-1] = rate;

            double fJ_r1 = xk[num-1] * y_ext[r2] * nHr;
            dJ(r1,r1) -= fJ_r1;  dJ(r2,r1) -= fJ_r1;
            dJ(p1,r1) += fJ_r1;  dJ(p2,r1) += fJ_r1;

            double fJ_r2 = xk[num-1] * y_ext[r1] * nHr;
            dJ(r1,r2) -= fJ_r2;  dJ(r2,r2) -= fJ_r2;
            dJ(p1,r2) += fJ_r2;  dJ(p2,r2) += fJ_r2;
        }
    }

    // ── Copy extended arrays [1..89] back to compact 0-based output ──────────
    for (int i = 0; i < 89; ++i)
        r_f[i] = r_f_dum[i+1];

    for (int i = 0; i < 89; ++i)
        for (int j = 0; j < 89; ++j)
            dr_fdy[i*89 + j] = dr_fdy_dum[(i+1)*102 + (j+1)];
}

} // namespace chemistry
