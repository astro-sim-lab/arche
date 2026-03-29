#pragma once
// ---------------------------------------------------------------------------
// reaction_metal_grain_surface.h
//
// Grain surface reaction rate coefficients for the metal_grain network.
//
// Port of:
//   fortran/metal_grain/subs/reaction_metal_grain.f
//     SUBROUTINE grain_coef  (L4088–5492)
//
// Output: grain_coef_rates() fills xkgr[0..150]  (151 entries, 0-based).
// These map to tbl.react[820..970] in the main reaction table
// (Fortran xkgr(1..151), giving global xk(991..1141)).
//
// NOTE: xkgr[129]  (Fortran xkgr(130)) is absent — never assigned, stays 0.
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
#include <cmath>
#include <array>
#include <algorithm>
#include "grain.h"   // vol_gr()

namespace chemistry {

// ---------------------------------------------------------------------------
// grain_coef_rates
//
//   Computes 151 grain surface reaction rate coefficients.
//
//   Fortran reaction index → C++ array index:  xkgr(n) → xkgr[n-1]
//   Note: index 130 (C++ [129]) is absent/zero.
//
//   Group summary (Fortran 1-based):
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
//     xnH    : hydrogen number density [cm^{-3}]
//     T_K    : gas temperature [K]
//     T_gr_K : grain temperature [K]
//     Z_metal: metallicity relative to solar
//     xJH2   : y[H2_p]  = H2 physisorbed abundance (per H)
//     xJH2O  : y[H2O_p] = H2O physisorbed abundance (per H)
//     xJtot  : y[H_c] + y[D_c]  = total chemisorbed H+D abundance (per H)
//     zeta   : cosmic-ray ionization rate [s^{-1}]
//     xkgr   : output array, length 151
// ---------------------------------------------------------------------------
inline void grain_coef_rates(double xnH,   double T_K,    double T_gr_K,
                              double Z_metal,
                              double xJH2,  double xJH2O,  double xJtot,
                              double zeta, double T_cr_eff,
                              std::array<double, 151>& xkgr)
{
    // ── Physical constants ────────────────────────────────────────────────
    constexpr double xk_B = phys::xk_B;
    constexpr double pi   = phys::pi;
    constexpr double xm_p = phys::xm_p;
    constexpr double h_b  = phys::hbar;   // ℏ [erg·s]
    constexpr double yHe = abundance_ref::yHe;

    // ── Species masses (Fortran mass(1..21), C++ 0-based) ─────────────────
    // 1:H(p), 2:H(c), 3:H2(p), 4:D(p), 5:D(c), 6:HD(p),
    // 7:O(p), 8:O2(p), 9:OH(p), 10:CO(p), 11:CO2(p), 12:H2O(p),
    // 13:HO2(p), 14:H2O2(p), 15:HCO(p), 16:H2CO(p),
    // 17:C(p), 18:CH(p), 19:CH2(p), 20:CH3(p), 21:CH4(p)
    static constexpr double mass[21] = {
        1.00783,  1.00783,  2.01565,  // H, H, H2
        2.01410,  2.01410,  3.02204,  // D, D, HD
        15.99491, 31.98983, 17.00274, 27.99491,  // O, O2, OH, CO
        43.98983, 18.01056, 32.99765, 34.00548,  // CO2, H2O, HO2, H2O2
        29.00274, 30.01056,                       // HCO, H2CO
        12.0,     13.00783, 14.01565, 15.02348, 16.0313  // C,CH,CH2,CH3,CH4
    };

    // ── Derived quantities ────────────────────────────────────────────────
    double rho       = xnH * (1.0 + 4.0 * yHe) * xm_p;
    double stick     = 1.0 / (1.0 + 4.0e-2 * std::sqrt(T_K + T_gr_K)
                               + 2.0e-3 * T_K + 8.0e-6 * T_K * T_K);
    double nu_0      = 1.0e12;
    double alpha_MRN = 1.0e-21 * Z_metal
                     * (detail::vol_gr(rho, T_gr_K) / detail::vol_gr(rho, 2.725));
    double vel_H     = std::sqrt(8.0 * xk_B * T_K / (pi * xm_p));
    double xnd_nsite = 4.0e-5 * xnH / 9.0;

    double f_ice  = std::min(1.0, xJH2O * xnH / xnd_nsite);
    double f_bare = std::max(0.0, 1.0 - f_ice);
    double f_h2i  = std::min(1.0, xJH2  * xnH / xnd_nsite);
    double f_h2b  = std::max(0.0, 1.0 - f_h2i);
    double f_chem = xJtot * xnH / xnd_nsite;

    // Helper: tunneling + activation barrier P_max for given Ea [K] and masses
    auto P_tunnel = [&](double m_red_amu, double Ea_K) -> double {
        double m_red = xm_p * m_red_amu;
        return std::max(
            std::exp(-(2.0e-8 / h_b) * std::sqrt(2.0 * m_red * xk_B * Ea_K)),
            std::exp(-Ea_K / T_gr_K)
        );
    };

    xkgr.fill(0.0);

    // ==========================================================================
    // (1) Adsorption  xkgr(1..21)  →  C++ [0..20]
    //     X(g) → X(p):  k = α_MRN · stick · vel_H / √mass
    //     xkgr(2): H(g) → H(c)  includes chemisorption probability
    //     xkgr(5): D(g) → D(c)  same
    // ==========================================================================

    // 1: H(g) → H(p)
    xkgr[0] = alpha_MRN * stick * vel_H / std::sqrt(mass[0]);

    // 2: H(g) → H(c)  Esplugues 2016 appendix B
    {
        double xm_red = xm_p * mass[0];
        double P = std::exp(-(2.0e-8 / h_b) * std::sqrt(2.0 * xm_red * xk_B * 1.0e3))
                 + std::exp(-1.0e3 / T_gr_K);
        xkgr[1] = alpha_MRN * stick * (vel_H / std::sqrt(mass[0])) * P * (1.0 - f_chem);
    }

    // 3: H2(g) → H2(p)
    xkgr[2] = alpha_MRN * stick * vel_H / std::sqrt(mass[2]);

    // 4: D(g) → D(p)
    xkgr[3] = alpha_MRN * stick * vel_H / std::sqrt(mass[3]);

    // 5: D(g) → D(c)
    {
        double xm_red = xm_p * mass[3];
        double P = std::exp(-(2.0e-8 / h_b) * std::sqrt(2.0 * xm_red * xk_B * 1.0e3))
                 + std::exp(-1.0e3 / T_gr_K);
        xkgr[4] = alpha_MRN * stick * (vel_H / std::sqrt(mass[3])) * P * (1.0 - f_chem);
    }

    // 6: HD(g) → HD(p)
    xkgr[5] = alpha_MRN * stick * vel_H / std::sqrt(mass[5]);

    // 7..16: O, O2, OH, CO, CO2, H2O, HO2, H2O2, HCO, H2CO
    for (int i = 6; i <= 15; ++i)
        xkgr[i] = alpha_MRN * stick * vel_H / std::sqrt(mass[i]);

    // 17..21: C, CH, CH2, CH3, CH4
    for (int i = 16; i <= 20; ++i)
        xkgr[i] = alpha_MRN * stick * vel_H / std::sqrt(mass[i]);

    // ==========================================================================
    // (2) Thermal desorption  xkgr(22..42)  →  C++ [21..41]
    // ==========================================================================

    // 22: H(p) → H(g)
    xkgr[21] = nu_0 * (f_bare * std::exp(-500.0   / T_gr_K)
                     + f_ice  * std::exp(-650.0   / T_gr_K));
    // 23: H(c) → H(g)
    xkgr[22] = nu_0 * std::exp(-1.0e4 / T_gr_K);
    // 24: H2(p) → H2(g)
    xkgr[23] = nu_0 * (f_h2b * (f_bare * std::exp(-300.0 / T_gr_K)
                              +  f_ice  * std::exp(-500.0 / T_gr_K))
                     + f_h2i * std::exp(-100.0 / T_gr_K));
    // 25: D(p) → D(g)   [H + 58K: Thi et al.]
    xkgr[24] = nu_0 * (f_bare * std::exp(-558.0   / T_gr_K)
                     + f_ice  * std::exp(-708.0   / T_gr_K));
    // 26: D(c) → D(g)
    xkgr[25] = nu_0 * std::exp(-1.0e4 / T_gr_K);
    // 27: HD(p) → HD(g)
    xkgr[26] = nu_0 * (f_h2b * (f_bare * std::exp(-358.0 / T_gr_K)
                              +  f_ice  * std::exp(-558.0 / T_gr_K))
                     + f_h2i * std::exp(-158.0 / T_gr_K));
    // 28: O(p) → O(g)
    xkgr[27] = nu_0 * (f_bare * std::exp(-1500.0  / T_gr_K)
                     + f_ice  * std::exp(-1420.0  / T_gr_K));
    // 29: O2(p) → O2(g)
    xkgr[28] = nu_0 * (f_bare * std::exp(-1250.0  / T_gr_K)
                     + f_ice  * std::exp(-1160.0  / T_gr_K));
    // 30: OH(p) → OH(g)
    xkgr[29] = nu_0 * std::exp(-4600.0 / T_gr_K);
    // 31: CO(p) → CO(g)
    xkgr[30] = nu_0 * (f_bare * std::exp(-1200.0  / T_gr_K)
                     + f_ice  * std::exp(-1300.0  / T_gr_K));
    // 32: CO2(p) → CO2(g)
    xkgr[31] = nu_0 * (f_bare * std::exp(-3000.0  / T_gr_K)
                     + f_ice  * std::exp(-2670.0  / T_gr_K));
    // 33: H2O(p) → H2O(g)
    xkgr[32] = nu_0 * (f_bare * std::exp(-4800.0  / T_gr_K)
                     + f_ice  * std::exp(-5700.0  / T_gr_K));
    // 34: HO2(p) → O(g) + OH(g)
    xkgr[33] = nu_0 * std::exp(-4000.0 / T_gr_K);
    // 35: H2O2(p) → H2O2(g)
    xkgr[34] = nu_0 * std::exp(-6000.0 / T_gr_K);
    // 36: HCO(p) → HCO(g)
    xkgr[35] = nu_0 * std::exp(-1600.0 / T_gr_K);
    // 37: H2CO(p) → H2CO(g)
    xkgr[36] = nu_0 * (f_bare * std::exp(-3700.0  / T_gr_K)
                     + f_ice  * std::exp(-3250.0  / T_gr_K));
    // 38: C(p) → C(g)
    xkgr[37] = nu_0 * std::exp(-800.0  / T_gr_K);
    // 39: CH(p) → CH(g)
    xkgr[38] = nu_0 * std::exp(-870.0  / T_gr_K);
    // 40: CH2(p) → CH2(g)
    xkgr[39] = nu_0 * std::exp(-945.0  / T_gr_K);
    // 41: CH3(p) → CH3(g)
    xkgr[40] = nu_0 * std::exp(-1017.0 / T_gr_K);
    // 42: CH4(p) → CH4(g)
    xkgr[41] = nu_0 * std::exp(-1090.0 / T_gr_K);

    // ==========================================================================
    // (3) Surface reactions, single product stays on grain  xkgr(43..58) → [42..57]
    //   P_bare = f_bare*(exp(-2E1/(3T)) + exp(-2E2/(3T)))
    //   P_ice  = f_ice *(exp(-2E1'/(3T))+ exp(-2E2'/(3T)))
    //   k = (nu_0/N_site)*(P_bare*(1-alpha_b)+P_ice*(1-alpha_i)) [/2 homonuclear]
    // ==========================================================================

    auto Pbr = [&](double E1b, double E2b) -> double {
        return f_bare * (std::exp(-2.0*E1b/(3.0*T_gr_K))
                       + std::exp(-2.0*E2b/(3.0*T_gr_K)));
    };
    auto Pic = [&](double E1i, double E2i) -> double {
        return f_ice  * (std::exp(-2.0*E1i/(3.0*T_gr_K))
                       + std::exp(-2.0*E2i/(3.0*T_gr_K)));
    };

    // 43: H(p)+H(p) → H2(p)
    xkgr[42] = (nu_0/xnd_nsite) *
               (Pbr(500.,500.)*(1.-9.630e-1) + Pic(650.,650.)*(1.-9.640e-2)) / 2.0;
    // 44: H(p)+D(p) → HD(p)
    xkgr[43] = (nu_0/xnd_nsite) *
               (Pbr(500.,500.)*(1.-9.630e-1) + Pic(650.,650.)*(1.-9.640e-2));
    // 45: H(p)+O(p) → OH(p)
    xkgr[44] = (nu_0/xnd_nsite) *
               (Pbr(500.,1500.)*(1.-3.875e-1) + Pic(650.,1420.)*(1.-3.880e-2));
    // 46: H(p)+OH(p) → H2O(p)
    xkgr[45] = (nu_0/xnd_nsite) *
               (Pbr(500.,4600.)*(1.-2.677e-1) + Pic(650.,4600.)*(1.-2.680e-2));
    // 47: H(p)+O2(p) → HO2(p)
    xkgr[46] = (nu_0/xnd_nsite) *
               (Pbr(500.,1250.)*(1.-1.380e-2) + Pic(650.,1160.)*(1.-1.400e-3));
    // 48: H(p)+CO(p) → HCO(p)   (with tunneling, Ea=2000K)
    {
        double Pb = Pbr(500.,1200.);
        double Pi = Pic(650.,1300.);
        double mu = mass[0]*mass[9]/(mass[0]+mass[9]);
        double Pm = P_tunnel(mu, 2000.0);
        double Pr = Pm / (Pm + Pb + Pi);
        xkgr[47] = (nu_0/xnd_nsite) * Pr *
                   (Pb*(1.-6.700e-3) + Pi*(1.-7.000e-4));
    }
    // 49: H(p)+HO2(p) → H2O2(p)
    xkgr[48] = (nu_0/xnd_nsite) *
               (Pbr(500.,4000.)*(1.-4.600e-3) + Pic(650.,4000.)*(1.-5.000e-4));
    // 50: H(p)+HCO(p) → H2CO(p)
    xkgr[49] = (nu_0/xnd_nsite) *
               (Pbr(500.,1600.)*(1.-6.610e-2) + Pic(650.,1600.)*(1.-6.700e-3));
    // 51: H(p)+C(p) → CH(p)
    xkgr[50] = (nu_0/xnd_nsite) *
               (Pbr(500.,800.)*(1.-8.212e-1) + Pic(650.,800.)*(1.-8.220e-2));
    // 52: H(p)+CH(p) → CH2(p)
    xkgr[51] = (nu_0/xnd_nsite) *
               (Pbr(500.,870.)*(1.-7.668e-1) + Pic(650.,870.)*(1.-7.670e-2));
    // 53: H(p)+CH2(p) → CH3(p)
    xkgr[52] = (nu_0/xnd_nsite) *
               (Pbr(500.,945.)*(1.-6.937e-1) + Pic(650.,945.)*(1.-6.940e-2));
    // 54: H(p)+CH3(p) → CH4(p)
    xkgr[53] = (nu_0/xnd_nsite) *
               (Pbr(500.,1017.)*(1.-5.886e-1) + Pic(650.,1017.)*(1.-5.890e-2));
    // 55: O(p)+O(p) → O2(p)
    xkgr[54] = (nu_0/xnd_nsite) *
               (Pbr(1500.,1500.)*(1.-6.884e-1) + Pic(1420.,1420.)*(1.-6.890e-2)) / 2.0;
    // 56: O(p)+C(p) → CO(p)
    xkgr[55] = (nu_0/xnd_nsite) *
               (Pbr(1500.,800.)*(1.-8.659e-1) + Pic(1420.,800.)*(1.-8.660e-2));
    // 57: O(p)+CO(p) → CO2(p)  (with tunneling, Ea=650K)
    // NOTE: Fortran line 4514 has typo: dexp(-6.5+2/T_gr_K) instead of dexp(-650/T_gr_K)
    {
        double Pb = Pbr(1500.,1200.);
        double Pi = Pic(1420.,1300.);
        double mu = mass[6]*mass[9]/(mass[6]+mass[9]);
        // Reproducing Fortran typo exactly: exp(-6.5 + 2.0/T_gr_K)
        double Pm = std::max(
            std::exp(-(2.0e-8/h_b) * std::sqrt(2.0 * xm_p*mu * xk_B * 650.0)),
            std::exp(-6.5 + 2.0/T_gr_K));   // Fortran typo: should be exp(-650/T_gr_K)
        double Pr = Pm / (Pm + Pb + Pi);
        xkgr[56] = (nu_0/xnd_nsite) * Pr *
                   (Pb*(1.-1.403e-1) + Pi*(1.-1.400e-2));
    }
    // 58: OH(p)+OH(p) → H2O2(p)
    xkgr[57] = (nu_0/xnd_nsite) *
               (Pbr(4600.,4600.)*(1.-2.000e-4) + Pic(4600.,4600.)*(1.-1.000e-4)) / 2.0;

    // ==========================================================================
    // (4) Surface reactions, single product desorbs  xkgr(59..74) → [58..73]
    //   k = (nu_0/N_site)*(P_bare*alpha_b + P_ice*alpha_i) [/2 homonuclear]
    // ==========================================================================

    // 59: H(p)+H(p) → H2(g)
    xkgr[58] = (nu_0/xnd_nsite) *
               (Pbr(500.,500.)*9.630e-1 + Pic(650.,650.)*9.640e-2) / 2.0;
    // 60: H(p)+D(p) → HD(g)
    xkgr[59] = (nu_0/xnd_nsite) *
               (Pbr(500.,500.)*9.630e-1 + Pic(650.,650.)*9.640e-2);
    // 61: H(p)+O(p) → OH(g)
    xkgr[60] = (nu_0/xnd_nsite) *
               (Pbr(500.,1500.)*3.875e-1 + Pic(650.,1420.)*3.880e-2);
    // 62: H(p)+OH(p) → H2O(g)
    xkgr[61] = (nu_0/xnd_nsite) *
               (Pbr(500.,4600.)*2.677e-1 + Pic(650.,4600.)*2.680e-2);
    // 63: H(p)+O2(p) → HO2(g)
    xkgr[62] = (nu_0/xnd_nsite) *
               (Pbr(500.,1250.)*1.380e-2 + Pic(650.,1160.)*1.400e-3);
    // 64: H(p)+CO(p) → HCO(g)  (with tunneling, Ea=2000K)
    {
        double Pb = Pbr(500.,1200.);
        double Pi = Pic(650.,1300.);
        double mu = mass[0]*mass[9]/(mass[0]+mass[9]);
        double Pm = P_tunnel(mu, 2000.0);
        double Pr = Pm / (Pm + Pb + Pi);
        xkgr[63] = (nu_0/xnd_nsite) * Pr *
                   (Pb*6.700e-3 + Pi*7.000e-4);
    }
    // 65: H(p)+HO2(p) → H2O2(g)
    xkgr[64] = (nu_0/xnd_nsite) *
               (Pbr(500.,4000.)*4.600e-3 + Pic(650.,4000.)*5.000e-4);
    // 66: H(p)+HCO(p) → H2CO(g)
    xkgr[65] = (nu_0/xnd_nsite) *
               (Pbr(500.,1600.)*6.610e-2 + Pic(650.,1600.)*6.700e-3);
    // 67: H(p)+C(p) → CH(g)
    xkgr[66] = (nu_0/xnd_nsite) *
               (Pbr(500.,800.)*8.212e-1 + Pic(650.,800.)*8.220e-2);
    // 68: H(p)+CH(p) → CH2(g)
    xkgr[67] = (nu_0/xnd_nsite) *
               (Pbr(500.,870.)*7.668e-1 + Pic(650.,870.)*7.670e-2);
    // 69: H(p)+CH2(p) → CH3(g)
    xkgr[68] = (nu_0/xnd_nsite) *
               (Pbr(500.,945.)*6.937e-1 + Pic(650.,945.)*6.940e-2);
    // 70: H(p)+CH3(p) → CH4(g)
    xkgr[69] = (nu_0/xnd_nsite) *
               (Pbr(500.,1017.)*5.886e-1 + Pic(650.,1017.)*5.890e-2);
    // 71: O(p)+O(p) → O2(g)
    xkgr[70] = (nu_0/xnd_nsite) *
               (Pbr(1500.,1500.)*6.884e-1 + Pic(1420.,1420.)*6.890e-2) / 2.0;
    // 72: O(p)+C(p) → CO(g)
    xkgr[71] = (nu_0/xnd_nsite) *
               (Pbr(1500.,800.)*8.659e-1 + Pic(1420.,800.)*8.660e-2);
    // 73: O(p)+CO(p) → CO2(g)  (with tunneling, Ea=650K)
    // NOTE: same Fortran typo as reaction 57
    {
        double Pb = Pbr(1500.,1200.);
        double Pi = Pic(1420.,1300.);
        double mu = mass[6]*mass[9]/(mass[6]+mass[9]);
        double Pm = std::max(
            std::exp(-(2.0e-8/h_b) * std::sqrt(2.0 * xm_p*mu * xk_B * 650.0)),
            std::exp(-6.5 + 2.0/T_gr_K));   // Fortran typo, faithful port
        double Pr = Pm / (Pm + Pb + Pi);
        xkgr[72] = (nu_0/xnd_nsite) * Pr *
                   (Pb*1.403e-1 + Pi*1.400e-2);
    }
    // 74: OH(p)+OH(p) → H2O2(g)
    xkgr[73] = (nu_0/xnd_nsite) *
               (Pbr(4600.,4600.)*2.000e-4 + Pic(4600.,4600.)*1.000e-4) / 2.0;

    // ==========================================================================
    // (5) Surface reactions, two products  xkgr(75..116) → [74..115]
    // ==========================================================================

    // Helper: tunneling P_react given (mu_amu, Ea_K, Pb, Pi)
    auto Pr_tun = [&](double mu_amu, double Ea_K, double Pb, double Pi) -> double {
        double Pm = P_tunnel(mu_amu, Ea_K);
        return Pm / (Pm + Pb + Pi);
    };

    // 75: H(p)+H2O(p) → H2(p)+OH(p)   Ea=9600K
    {
        double Pb = Pbr(500.,4800.);
        double Pi = Pic(650.,5700.);
        double mu = mass[0]*mass[11]/(mass[0]+mass[11]);
        double Pr = Pr_tun(mu, 9600.0, Pb, Pi);
        xkgr[74] = (nu_0/xnd_nsite) * Pr * (Pb + Pi);
    }
    // 76: H(p)+HO2(p) → OH(p)+OH(p)
    xkgr[75] = (nu_0/xnd_nsite) *
               (Pbr(500.,4000.)*(1.-3.400e-3) + Pic(650.,4000.)*(1.-4.000e-4));
    // 77: H(p)+HO2(p) → OH(g)+OH(g)
    xkgr[76] = (nu_0/xnd_nsite) *
               (Pbr(500.,4000.)*3.400e-3 + Pic(650.,4000.)*4.000e-4);
    // 78: H(p)+H2O2(p) → OH(p)+H2O(p)  Ea=1000K
    {
        double Pb = Pbr(500.,6000.);
        double Pi = Pic(650.,6000.);
        double mu = mass[0]*mass[13]/(mass[0]+mass[13]);
        double Pr = Pr_tun(mu, 1000.0, Pb, Pi);
        xkgr[77] = (nu_0/xnd_nsite) * Pr *
                   (Pb*9.716e-1 + Pi*9.971e-1);
    }
    // 79: H(p)+H2O2(p) → OH(g)+H2O(p)  Ea=1000K
    {
        double Pb = Pbr(500.,6000.);
        double Pi = Pic(650.,6000.);
        double mu = mass[0]*mass[13]/(mass[0]+mass[13]);
        double Pr = Pr_tun(mu, 1000.0, Pb, Pi);
        xkgr[78] = (nu_0/xnd_nsite) * Pr *
                   (Pb*7.200e-3 + Pi*8.000e-4);
    }
    // 80: H(p)+H2O2(p) → OH(g)+H2O(g)  Ea=1000K
    {
        double Pb = Pbr(500.,6000.);
        double Pi = Pic(650.,6000.);
        double mu = mass[0]*mass[13]/(mass[0]+mass[13]);
        double Pr = Pr_tun(mu, 1000.0, Pb, Pi);
        xkgr[79] = (nu_0/xnd_nsite) * Pr *
                   (Pb*2.120e-2 + Pi*2.100e-3);
    }
    // 81: H(p)+HCO(p) → H2(p)+CO(p)
    xkgr[80] = (nu_0/xnd_nsite) *
               (Pbr(500.,1600.)*8.260e-2 + Pic(650.,1600.)*9.082e-1);
    // 82: H(p)+HCO(p) → H2(g)+CO(p)
    xkgr[81] = (nu_0/xnd_nsite) *
               (Pbr(500.,1600.)*4.827e-1 + Pic(650.,1600.)*4.820e-2);
    // 83: H(p)+HCO(p) → H2(g)+CO(g)
    xkgr[82] = (nu_0/xnd_nsite) *
               (Pbr(500.,1600.)*4.347e-1 + Pic(650.,1600.)*4.360e-2);
    // 84: H(p)+H2CO(p) → H2(p)+HCO(p)  Ea=2200K
    {
        double Pb = Pbr(500.,3700.);
        double Pi = Pic(650.,3250.);
        double mu = mass[0]*mass[15]/(mass[0]+mass[15]);
        double Pr = Pr_tun(mu, 2200.0, Pb, Pi);
        xkgr[83] = (nu_0/xnd_nsite) * Pr *
                   (Pb*4.948e-1 + Pi*9.494e-1);
    }
    // 85: H(p)+H2CO(p) → H2(g)+HCO(p)  Ea=2200K
    {
        double Pb = Pbr(500.,3700.);
        double Pi = Pic(650.,3250.);
        double mu = mass[0]*mass[15]/(mass[0]+mass[15]);
        double Pr = Pr_tun(mu, 2200.0, Pb, Pi);
        xkgr[84] = (nu_0/xnd_nsite) * Pr *
                   (Pb*5.050e-1 + Pi*5.050e-2);
    }
    // 86: H(p)+H2CO(p) → H2(g)+HCO(g)  Ea=2200K
    {
        double Pb = Pbr(500.,3700.);
        double Pi = Pic(650.,3250.);
        double mu = mass[0]*mass[15]/(mass[0]+mass[15]);
        double Pr = Pr_tun(mu, 2200.0, Pb, Pi);
        xkgr[85] = (nu_0/xnd_nsite) * Pr *
                   (Pb*2.000e-4 + Pi*1.000e-4);
    }
    // 87: H(p)+CO2(p) → OH(p)+CO(p)  Ea=10000K
    {
        double Pb = Pbr(500.,3000.);
        double Pi = Pic(650.,2670.);
        double mu = mass[0]*mass[10]/(mass[0]+mass[10]);
        double Pr = Pr_tun(mu, 10000.0, Pb, Pi);
        xkgr[86] = (nu_0/xnd_nsite) * Pr * (Pb + Pi);
    }
    // 88: H(p)+CH(p) → H2(p)+C(p)
    xkgr[87] = (nu_0/xnd_nsite) *
               (Pbr(500.,870.)*2.226e-1 + Pic(650.,870.)*9.222e-1);
    // 89: H(p)+CH(p) → H2(g)+C(p)
    xkgr[88] = (nu_0/xnd_nsite) *
               (Pbr(500.,870.)*3.859e-1 + Pic(650.,870.)*3.865e-2);
    // 90: H(p)+CH(p) → H2(g)+C(g)
    xkgr[89] = (nu_0/xnd_nsite) *
               (Pbr(500.,870.)*3.915e-1 + Pic(650.,870.)*3.915e-2);
    // 91: H(p)+CH2(p) → H2(p)+CH(p)
    xkgr[90] = (nu_0/xnd_nsite) *
               (Pbr(500.,945.)*9.554e-1 + Pic(650.,945.)*9.955e-1);
    // 92: H(p)+CH2(p) → H2(g)+CH(p)
    xkgr[91] = (nu_0/xnd_nsite) *
               (Pbr(500.,945.)*4.452e-2 + Pic(650.,945.)*4.452e-3);
    // 93: H(p)+CH2(p) → H2(g)+CH(g)
    xkgr[92] = (nu_0/xnd_nsite) *
               (Pbr(500.,945.)*8.000e-5 + Pic(650.,945.)*4.800e-5);
    // 94: H(p)+CH3(p) → H2(p)+CH2(p)
    xkgr[93] = (nu_0/xnd_nsite) * (Pbr(500.,1017.) + Pic(650.,1017.));
    // 95: H(p)+CH4(p) → H2(p)+CH3(p)
    xkgr[94] = (nu_0/xnd_nsite) * (Pbr(500.,1090.) + Pic(650.,1090.));
    // 96: O(p)+OH(p) → H(p)+O2(p)
    xkgr[95] = (nu_0/xnd_nsite) *
               (Pbr(1500.,4600.)*4.547e-1 + Pic(1420.,4600.)*9.454e-1);
    // 97: O(p)+OH(p) → H(g)+O2(p)
    xkgr[96] = (nu_0/xnd_nsite) *
               (Pbr(1500.,4600.)*5.264e-1 + Pic(1420.,4600.)*5.260e-2);
    // 98: O(p)+OH(p) → H(g)+O2(g)
    xkgr[97] = (nu_0/xnd_nsite) *
               (Pbr(1500.,4600.)*1.890e-2 + Pic(1420.,4600.)*2.000e-3);
    // 99: O(p)+HO2(p) → O2(p)+OH(p)
    xkgr[98] = (nu_0/xnd_nsite) *
               (Pbr(1500.,4000.)*8.265e-1 + Pic(1420.,4000.)*9.826e-1);
    // 100: O(p)+HO2(p) → O2(g)+OH(p)
    xkgr[99] = (nu_0/xnd_nsite) *
               (Pbr(1500.,4000.)*1.516e-1 + Pic(1420.,4000.)*1.510e-2);
    // 101: O(p)+HO2(p) → O2(g)+OH(g)
    xkgr[100] = (nu_0/xnd_nsite) *
                (Pbr(1500.,4000.)*2.190e-2 + Pic(1420.,4000.)*2.300e-3);
    // 102: O(p)+HCO(p) → H(p)+CO2(p)
    xkgr[101] = (nu_0/xnd_nsite) *
                (Pbr(1500.,1600.)*1.141e-1 + Pic(1420.,1600.)*9.114e-1);
    // 103: O(p)+HCO(p) → H(g)+CO2(p)
    xkgr[102] = (nu_0/xnd_nsite) *
                (Pbr(1500.,1600.)*8.348e-1 + Pic(1420.,1600.)*8.340e-2);
    // 104: O(p)+HCO(p) → H(g)+CO2(g)
    xkgr[103] = (nu_0/xnd_nsite) *
                (Pbr(1500.,1600.)*5.110e-2 + Pic(1420.,1600.)*5.200e-3);
    // 105: O(p)+H2CO(p) → H2(p)+CO2(p)  Ea=335K
    {
        double Pb = Pbr(1500.,3700.);
        double Pi = Pic(1420.,3250.);
        double mu = mass[6]*mass[15]/(mass[6]+mass[15]);
        double Pr = Pr_tun(mu, 335.0, Pb, Pi);
        xkgr[104] = (nu_0/xnd_nsite) * Pr *
                    (Pb*7.310e-2 + Pi*9.073e-1);
    }
    // 106: O(p)+H2CO(p) → H2(g)+CO2(p)  Ea=335K
    {
        double Pb = Pbr(1500.,3700.);
        double Pi = Pic(1420.,3250.);
        double mu = mass[6]*mass[15]/(mass[6]+mass[15]);
        double Pr = Pr_tun(mu, 335.0, Pb, Pi);
        xkgr[105] = (nu_0/xnd_nsite) * Pr *
                    (Pb*8.901e-1 + Pi*8.900e-2);
    }
    // 107: O(p)+H2CO(p) → H2(g)+CO2(g)  Ea=335K
    {
        double Pb = Pbr(1500.,3700.);
        double Pi = Pic(1420.,3250.);
        double mu = mass[6]*mass[15]/(mass[6]+mass[15]);
        double Pr = Pr_tun(mu, 335.0, Pb, Pi);
        xkgr[106] = (nu_0/xnd_nsite) * Pr *
                    (Pb*3.680e-2 + Pi*3.700e-3);
    }
    // 108: H2(p)+OH(p) → H(p)+H2O(p)  Ea=2100K
    {
        double Pb = Pbr(300.,4600.);
        double Pi = Pic(500.,4600.);
        double mu = mass[2]*mass[8]/(mass[2]+mass[8]);
        double Pr = Pr_tun(mu, 2100.0, Pb, Pi);
        xkgr[107] = (nu_0/xnd_nsite) * Pr *
                    (Pb*5.948e-1 + Pi*9.594e-1);
    }
    // 109: H2(p)+OH(p) → H(g)+H2O(p)  Ea=2100K
    {
        double Pb = Pbr(300.,4600.);
        double Pi = Pic(500.,4600.);
        double mu = mass[2]*mass[8]/(mass[2]+mass[8]);
        double Pr = Pr_tun(mu, 2100.0, Pb, Pi);
        xkgr[108] = (nu_0/xnd_nsite) * Pr *
                    (Pb*4.052e-1 + Pi*4.060e-2);
    }
    // 110: OH(p)+CO(p) → H(p)+CO2(p)  Ea=400K
    {
        double Pb = Pbr(4600.,1200.);
        double Pi = Pic(4600.,1300.);
        double mu = mass[8]*mass[9]/(mass[8]+mass[9]);
        double Pr = Pr_tun(mu, 400.0, Pb, Pi);
        xkgr[109] = (nu_0/xnd_nsite) * Pr *
                    (Pb*4.205e-1 + Pi*9.420e-1);
    }
    // 111: OH(p)+CO(p) → H(g)+CO2(p)  Ea=400K
    {
        double Pb = Pbr(4600.,1200.);
        double Pi = Pic(4600.,1300.);
        double mu = mass[8]*mass[9]/(mass[8]+mass[9]);
        double Pr = Pr_tun(mu, 400.0, Pb, Pi);
        xkgr[110] = (nu_0/xnd_nsite) * Pr *
                    (Pb*5.794e-1 + Pi*5.794e-2);
    }
    // 112: OH(p)+CO(p) → H(g)+CO2(g)  Ea=400K
    {
        double Pb = Pbr(4600.,1200.);
        double Pi = Pic(4600.,1300.);
        double mu = mass[8]*mass[9]/(mass[8]+mass[9]);
        double Pr = Pr_tun(mu, 400.0, Pb, Pi);
        xkgr[111] = (nu_0/xnd_nsite) * Pr *
                    (Pb*1.000e-4 + Pi*6.000e-5);
    }
    // 113: OH(p)+HCO(p) → H2(p)+CO2(p)
    xkgr[112] = (nu_0/xnd_nsite) *
                (Pbr(4600.,1600.)*8.063e-2 + Pic(4600.,1600.)*9.080e-1);
    // 114: OH(p)+HCO(p) → H2(g)+CO2(p)
    xkgr[113] = (nu_0/xnd_nsite) *
                (Pbr(4600.,1600.)*8.936e-1 + Pic(4600.,1600.)*8.936e-2);
    // 115: OH(p)+HCO(p) → H2(g)+CO2(g)
    xkgr[114] = (nu_0/xnd_nsite) *
                (Pbr(4600.,1600.)*2.577e-2 + Pic(4600.,1600.)*2.640e-3);
    // 116: H2(p)+HO2(p) → H(p)+H2O2(p)  Ea=5000K
    {
        double Pb = Pbr(300.,4000.);
        double Pi = Pic(500.,4000.);
        double mu = mass[2]*mass[12]/(mass[2]+mass[12]);
        double Pr = Pr_tun(mu, 5000.0, Pb, Pi);
        xkgr[115] = (nu_0/xnd_nsite) * Pr * (Pb + Pi);
    }

    // ==========================================================================
    // (6) Chemisorption reactions  xkgr(117..129) → [116..128]
    // ==========================================================================

    // Cazaux & Tielens 2010 helper for physisorbed → chemisorbed transition
    // Parameters: E_Hc [K], E_Hp [K], E_s [K], mass [amu]
    // Returns: k_transition factor (not yet multiplied by alpha_MRN)
    auto CT_rate = [&](double E_Hc, double E_Hp, double E_s, double m_amu) -> double {
        double xm_red = xm_p * m_amu;
        double P = std::exp(-(6.0e-8/h_b) * std::sqrt(2.0 * xm_red * xk_B * (E_Hp - E_s)));
        return nu_0 * (8.0 * std::sqrt(pi * T_gr_K) * P * std::sqrt(E_Hc - E_Hp)
                           / (E_Hc - E_s)
                     + 4.0 * std::sqrt((E_Hp - E_s) / (E_Hc - E_s))
                           * std::exp(-(E_Hp - E_s) / T_gr_K));
    };

    // 117: H(p) → H(c)
    xkgr[116] = CT_rate(1.0e4, 500.0, 200.0, mass[0]) * (1.0 - f_chem);

    // 118: H(g) + H(c) → H2(g)
    {
        double xm_red = xm_p * mass[0];
        double P = std::exp(-(2.0e-8/h_b) * std::sqrt(2.0 * xm_red * xk_B * 1.0e3))
                 + std::exp(-1.0e3 / T_gr_K);
        xkgr[117] = alpha_MRN * stick * (vel_H / std::sqrt(mass[0])) * P / xnd_nsite;
    }

    // 119: H(p) + H(c) → H2(g)
    xkgr[118] = CT_rate(1.0e4, 500.0, 200.0, mass[0]) / xnd_nsite;

    // 120: D(p) → D(c)
    xkgr[119] = CT_rate(1.0e4, 558.0, 200.0, mass[3]) * (1.0 - f_chem);

    // 121: H(g) + D(c) → HD(g)  [same rate as 118]
    xkgr[120] = xkgr[117];

    // 122: D(g) + H(c) → HD(g)
    {
        double xm_red = xm_p * mass[3];
        double P = std::exp(-(2.0e-8/h_b) * std::sqrt(2.0 * xm_red * xk_B * 1.0e3))
                 + std::exp(-1.0e3 / T_gr_K);
        xkgr[121] = alpha_MRN * stick * (vel_H / std::sqrt(mass[3])) * P / xnd_nsite;
    }

    // 123: H(p) + D(c) → HD(g)  [same rate as 119]
    xkgr[122] = xkgr[118];

    // 124: D(p) + H(c) → HD(g)
    xkgr[123] = CT_rate(1.0e4, 558.0, 200.0, mass[3]) / xnd_nsite;

    // 125: H(g) + H(p) → H2(g)  = xkgr[0] / xnd_nsite
    xkgr[124] = xkgr[0] / xnd_nsite;

    // 126: H(g) + D(p) → HD(g)  = xkgr[0] / xnd_nsite
    xkgr[125] = xkgr[0] / xnd_nsite;

    // 127: D(g) + H(p) → HD(g)  = xkgr[3] / xnd_nsite
    xkgr[126] = xkgr[3] / xnd_nsite;

    // 128: H(g) + H(g) → H2(g)  (simple formula, Cazaux&Tielens)
    {
        double xm_red = xm_p * mass[0];
        double Ptun   = std::exp(-(6.0e-8/h_b) * std::sqrt(2.0 * xm_red * xk_B * 300.0));
        // Reuse CT_rate-like expression for effective chemisorption probability
        double P_eff  = 8.0 * std::sqrt(pi * T_gr_K) * Ptun
                            * std::sqrt(1.0e4 - 500.0) / (1.0e4 - 200.0)
                      + 4.0 * std::sqrt((500.0 - 200.0) / (1.0e4 - 200.0))
                            * std::exp(-(500.0 - 200.0) / T_gr_K);
        xkgr[127] = alpha_MRN * stick * vel_H / std::sqrt(mass[0]) / 2.0
                  / (1.0 + std::exp(-500.0 / T_gr_K) / P_eff);
    }

    // 129: H(g) + D(g) → HD(g)
    {
        double xm_red = xm_p * mass[3];
        double Ptun   = std::exp(-(6.0e-8/h_b) * std::sqrt(2.0 * xm_red * xk_B * 358.0));
        double P_eff  = 8.0 * std::sqrt(pi * T_gr_K) * Ptun
                            * std::sqrt(1.0e4 - 558.0) / (1.0e4 - 200.0)
                      + 4.0 * std::sqrt((558.0 - 200.0) / (1.0e4 - 200.0))
                            * std::exp(-(558.0 - 200.0) / T_gr_K);
        xkgr[128] = alpha_MRN * stick * vel_H / std::sqrt(mass[3])
                  / (1.0 + std::exp(-558.0 / T_gr_K) / P_eff);
    }

    // xkgr[129] (Fortran 130) intentionally absent — stays 0.

    // ==========================================================================
    // (7) CR + CRUV desorption  xkgr(131..151) → [130..150]
    //   Pattern:  k_base = nu_0 * f(T_eff=T_cr_eff)
    //   k_CR = (3.16e-19 * k_base + 2.07e-16) * (zeta / 1.36e-17)
    //   [Hasegawa & Herbst 1993 eq. 9 + CRUV term]
    // ==========================================================================

    auto CR = [&](double k_base) -> double {
        return (3.16e-19 * k_base + 2.07e-16) * (zeta / 1.36e-17);
    };

    // 131: H(p) → H(g)
    xkgr[130] = CR(nu_0*(f_bare*std::exp(-500.0/T_cr_eff) + f_ice*std::exp(-650.0/T_cr_eff)));
    // 132: H(c) → H(g)
    xkgr[131] = CR(nu_0 * std::exp(-1.0e4/T_cr_eff));
    // 133: H2(p) → H2(g)
    xkgr[132] = CR(nu_0*(f_h2b*(f_bare*std::exp(-300.0/T_cr_eff)
                               + f_ice*std::exp(-500.0/T_cr_eff))
                       + f_h2i*std::exp(-100.0/T_cr_eff)));
    // 134: D(p) → D(g)
    xkgr[133] = CR(nu_0*(f_bare*std::exp(-558.0/T_cr_eff) + f_ice*std::exp(-708.0/T_cr_eff)));
    // 135: D(c) → D(g)
    xkgr[134] = CR(nu_0 * std::exp(-1.0e4/T_cr_eff));
    // 136: HD(p) → HD(g)
    xkgr[135] = CR(nu_0*(f_h2b*(f_bare*std::exp(-358.0/T_cr_eff)
                               + f_ice*std::exp(-558.0/T_cr_eff))
                       + f_h2i*std::exp(-158.0/T_cr_eff)));
    // 137: O(p) → O(g)
    xkgr[136] = CR(nu_0*(f_bare*std::exp(-1500.0/T_cr_eff) + f_ice*std::exp(-1420.0/T_cr_eff)));
    // 138: O2(p) → O2(g)
    xkgr[137] = CR(nu_0*(f_bare*std::exp(-1250.0/T_cr_eff) + f_ice*std::exp(-1160.0/T_cr_eff)));
    // 139: OH(p) → OH(g)
    xkgr[138] = CR(nu_0 * std::exp(-4600.0/T_cr_eff));
    // 140: CO(p) → CO(g)
    xkgr[139] = CR(nu_0*(f_bare*std::exp(-1200.0/T_cr_eff) + f_ice*std::exp(-1300.0/T_cr_eff)));
    // 141: CO2(p) → CO2(g)
    xkgr[140] = CR(nu_0*(f_bare*std::exp(-3000.0/T_cr_eff) + f_ice*std::exp(-2670.0/T_cr_eff)));
    // 142: H2O(p) → H2O(g)
    xkgr[141] = CR(nu_0*(f_bare*std::exp(-4800.0/T_cr_eff) + f_ice*std::exp(-5700.0/T_cr_eff)));
    // 143: HO2(p) → O(g)+OH(g)
    xkgr[142] = CR(nu_0 * std::exp(-4000.0/T_cr_eff));
    // 144: H2O2(p) → H2O2(g)
    xkgr[143] = CR(nu_0 * std::exp(-6000.0/T_cr_eff));
    // 145: HCO(p) → HCO(g)
    xkgr[144] = CR(nu_0 * std::exp(-1600.0/T_cr_eff));
    // 146: H2CO(p) → H2CO(g)
    xkgr[145] = CR(nu_0*(f_bare*std::exp(-3700.0/T_cr_eff) + f_ice*std::exp(-3250.0/T_cr_eff)));
    // 147: C(p) → C(g)
    xkgr[146] = CR(nu_0 * std::exp(-800.0/T_cr_eff));
    // 148: CH(p) → CH(g)
    xkgr[147] = CR(nu_0 * std::exp(-870.0/T_cr_eff));
    // 149: CH2(p) → CH2(g)
    xkgr[148] = CR(nu_0 * std::exp(-945.0/T_cr_eff));
    // 150: CH3(p) → CH3(g)
    xkgr[149] = CR(nu_0 * std::exp(-1017.0/T_cr_eff));
    // 151: CH4(p) → CH4(g)
    xkgr[150] = CR(nu_0 * std::exp(-1090.0/T_cr_eff));
}

} // namespace chemistry
