// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
#include <array>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include "reaction_table.h"
#include "state.h"

// ---------------------------------------------------------------------------
// Partition function evaluation for zero_metal species (N_sp = 23).
//
// Mirrors Fortran part_fnc_prm.f, H2_new.f (H2_pf / E_H2_BFM).
// Sources referenced in the Fortran comments:
//   P09   = Pagano+09
//   BC16  = Barklem & Collet 2016 (diatomic, 32-point hardcoded table)
//   PJ16  = Popovas & Jorgensen 2016 (H2 polynomial fits)
//   BFM89 = Borysow, Frommhold & Moraldi 1989 (H2 energy levels)
//   HITRAN / ExoMol — tabulated data read from pf_*.dat files
//
// eval_partition_functions() fills pf[1..N_sp] (1-based to match Fortran
// species indices); pf[0] = pf[101] = 1.0 act as sentinels for "vacant"
// reactant/product slots and photons in the reaction table.
// ---------------------------------------------------------------------------

namespace chemistry {
namespace detail {

// Linear interpolation on a sorted (T, value) vector.
// Clamps to boundary values outside the tabulated range.
inline double pf_interp(const std::vector<std::pair<double,double>>& tab,
                        double T_K)
{
    if (tab.empty()) return 1.0;
    if (T_K <= tab.front().first) return tab.front().second;
    if (T_K >= tab.back().first)  return tab.back().second;
    // Find first entry with T >= T_K
    auto it = std::lower_bound(tab.begin(), tab.end(),
                               std::pair<double,double>{T_K, 0.0},
                               [](const auto& a, const auto& b){
                                   return a.first < b.first; });
    auto lo = std::prev(it);
    double t = (T_K - lo->first) / (it->first - lo->first);
    return (1.0 - t) * lo->second + t * it->second;
}

// BC16 32-point temperature grid (shared by H2+, HeH+, Li, LiH)
constexpr std::array<double,32> kTBC = {
     1.0,   1.3,   1.7,    2.0,    3.0,    5.0,    7.0,   10.0,
    15.0,  20.0,  30.0,   50.0,   70.0,  100.0,  130.0,  170.0,
   200.0, 250.0, 300.0,  500.0,  700.0, 1000.0, 1500.0, 2000.0,
  3000.0,4000.0,5000.0, 6000.0, 7000.0, 8000.0, 9000.0,10000.0};

// Linear interpolation on the 32-point BC16 table
inline double bc16_interp(const std::array<double,32>& ytab, double T_K)
{
    int ms = 31;
    for (int i = 0; i < 32; ++i) {
        if (T_K <= kTBC[i]) { ms = i; break; }
    }
    if (ms == 0) ms = 1;
    double t = (T_K - kTBC[ms-1]) / (kTBC[ms] - kTBC[ms-1]);
    return (1.0-t)*ytab[ms-1] + t*ytab[ms];
}

// H2 vibrational+rotational energy levels, BFM 1989 (ApJ 336,495) eq.(A11)
// Returns E[K] for vibrational quantum number iv and rotational quantum J.
inline double E_H2_BFM(int iv, int J)
{
    constexpr double Ev0= 0.38496,  Ev1=-0.04609,  Ev2= 0.00178,
                     Ev3=-7.0e-5,   Ev4= 2.9511e-6;
    constexpr double Bv0=54.438,    Bv1= 4.6063,   Bv2=-2.0050,
                     Bv3= 0.19260,  Bv4=-6.2953e-3;
    constexpr double Dv0=-0.72593,  Dv1= 5.9990,   Dv2=-1.5187,
                     Dv3= 0.12721,  Dv4=-3.3391e-3;
    constexpr double Fv0=-12.662,   Fv1=20.047,    Fv2=-4.7873,
                     Fv3= 0.36900,  Fv4=-9.003e-3;
    constexpr double Gv0=-24.006,   Gv1=33.989,    Gv2=-7.6841,
                     Gv3= 0.52413,  Gv4=-1.1297e-2;
    constexpr double Hv0=-22.384,   Hv1=31.150,    Hv2=-6.6139,
                     Hv3= 0.39226,  Hv4=-9.496e-3;
    constexpr double Ov0=-10.541,   Ov1=14.746,    Ov2=-2.9476,
                     Ov3= 0.16016,  Ov4=-6.5005e-3;
    constexpr double Pv0=-2.0021,   Pv1= 2.8400,   Pv2=-0.54654,
                     Pv3= 0.031636, Pv4=-2.4398e-3;

    double vv = iv + 0.5;
    double vv2=vv*vv, vv3=vv2*vv, vv4=vv3*vv;
    double Ev = Ev0+Ev1*vv+Ev2*vv2+Ev3*vv3+Ev4*vv4;
    double Bv = Bv0+Bv1*vv+Bv2*vv2+Bv3*vv3+Bv4*vv4;
    double Dv = Dv0+Dv1*vv+Dv2*vv2+Dv3*vv3+Dv4*vv4;
    double Fv = Fv0+Fv1*vv+Fv2*vv2+Fv3*vv3+Fv4*vv4;
    double Gv = Gv0+Gv1*vv+Gv2*vv2+Gv3*vv3+Gv4*vv4;
    double Hv = Hv0+Hv1*vv+Hv2*vv2+Hv3*vv3+Hv4*vv4;
    double Ov = Ov0+Ov1*vv+Ov2*vv2+Ov3*vv3+Ov4*vv4;
    double Pv = Pv0+Pv1*vv+Pv2*vv2+Pv3*vv3+Pv4*vv4;

    double rrot = double(J*(J+1));
    return 1.43879 * (-Ev*1e5 + Bv*rrot
                      - Dv*1e-2*rrot*rrot
                      + Fv*1e-5*rrot*rrot*rrot
                      - Gv*1e-8*rrot*rrot*rrot*rrot
                      + Hv*1e-11*rrot*rrot*rrot*rrot*rrot
                      - Ov*1e-14*rrot*rrot*rrot*rrot*rrot*rrot
                      + Pv*1e-17*rrot*rrot*rrot*rrot*rrot*rrot*rrot);
}

// H2 partition function (equilibrium ortho/para 1:3).
// From H2_new.f: H2_pf()
inline double H2_pf(double T_K)
{
    constexpr int iv_max = 5, j_max = 25;
    double ET[iv_max+1][j_max+1];
    for (int iv = 0; iv <= iv_max; ++iv)
        for (int j = 0; j <= j_max; ++j)
            ET[iv][j] = E_H2_BFM(iv, j);

    const double ET_diss = 5.196e4;
    double z_p = 0.0, z_o = 0.0;
    for (int iv = 0; iv <= iv_max; ++iv)
        for (int j = 0; j <= j_max; ++j)
            if (ET[iv][j] - ET[0][0] <= ET_diss) {
                if (j % 2 == 0)
                    z_p += (2*j+1) * std::exp(-(ET[iv][j]-ET[0][0])/T_K);
                else
                    z_o += (2*j+1) * std::exp(-(ET[iv][j]-ET[0][1])/T_K);
            }
    return 0.25*z_p + 0.75*z_o;
}

} // namespace detail

// ---------------------------------------------------------------------------
// eval_partition_functions()
//
// Fills pf[0..101]:
//   pf[0]   = 1.0  (vacant-slot sentinel)
//   pf[1..N_sp] = species partition functions (1-based)
//   pf[101] = 1.0  (photon sentinel)
//
// Templated but currently only provides complete data for N_sp = 23
// (zero_metal).  Metal-grain species 24..89 would need extension.
// ---------------------------------------------------------------------------
template<int N_sp, int N_react>
void eval_partition_functions(double T_K,
                              const ReactionTable<N_sp, N_react>& tbl,
                              std::array<double, 102>& pf)
{
    pf.fill(1.0);   // default: sentinel value for vacant/photon slots

    double T_K2=T_K*T_K, T_K3=T_K2*T_K, T_K4=T_K3*T_K, T_K5=T_K4*T_K;
    constexpr double xk_B = phys::xk_B;

    // 1: H (P09) — ground state 2, + excited levels for T > 1000 K
    pf[1] = 2.0;
    if (T_K > 1.0e3) {
        for (int i = 2; i <= 5; ++i) {
            double di = double(i);
            pf[1] += (2.0*di*di)
                   * std::exp(-13.598*phys::eV_to_erg*(1.0 - 1.0/(di*di))/xk_B/T_K);
        }
    }

    // 2: H2 (Popovas & Jorgensen 2016 equilibrium fit + BFM cap)
    double xq_H2;
    if (T_K <= 200.0) {
        xq_H2 = 2.673071615415136e-1 - 3.495586051757211e-3*T_K
               + 1.227901619954258e-4*T_K2 - 5.776440695273789e-7*T_K3
               + 9.251224490175610e-10*T_K4;
    } else if (T_K <= 1000.0) {
        xq_H2 = 1.410033600294133e-1 + 6.085738724141971e-3*T_K
               - 4.096994909866605e-7*T_K2 + 4.220221708082499e-10*T_K3
               - 8.790594164685680e-14*T_K4;
    } else {
        xq_H2 = -9.661842638994980e-1 + 7.302127874247883e-3*T_K
               - 6.760893004505151e-7*T_K2 + 3.128741080316710e-10*T_K3
               - 1.645206030945910e-14*T_K4 + 2.788597060472472e-19*T_K5;
        xq_H2 = std::min(xq_H2, detail::H2_pf(T_K));
    }
    pf[2] = xq_H2;

    // 3: e (P09)
    pf[3] = 2.0;

    // 4: H+ (P09)
    pf[4] = 1.0;

    // 5: H2+ (BC16, 32-point)
    static constexpr std::array<double,32> xqh2p = {
        5.00000e-1, 5.00000e-1, 5.00000e-1, 5.00000e-1, 5.00000e-1,
        5.00000e-1, 5.00024e-1, 5.00921e-1, 5.15632e-1, 5.64397e-1,
        7.65738e-1, 1.33880e+0, 1.91022e+0, 2.68536e+0, 3.40918e+0,
        4.35256e+0, 5.05717e+0, 6.23172e+0, 7.40726e+0, 1.21361e+1,
        1.70214e+1, 2.50338e+1, 4.11228e+1, 6.13523e+1, 1.16205e+2,
        1.93054e+2, 2.94603e+2, 4.22340e+2, 5.76489e+2, 7.56297e+2,
        9.60417e+2, 1.18728e+3};
    pf[5] = detail::bc16_interp(xqh2p, T_K);

    // 6: H3+ (ExoMol file, nuclear-spin deg removed: /8)
    {
        double raw = detail::pf_interp(tbl.pf_table[6], T_K);
        pf[6] = std::max(1.0, raw / 8.0);
    }

    // 7: H- (P09)
    pf[7] = 1.0;

    // 8: He (P09)
    pf[8] = 1.0;

    // 9: He+ (P09)
    pf[9] = 2.0;

    // 10: He++ (P09)
    pf[10] = 1.0;

    // 11: HeH+ (BC16, 32-point)
    static constexpr std::array<double,32> xqhehp = {
        1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0,
        1.00000e+0, 1.00000e+0, 1.00019e+0, 1.00482e+0, 1.02408e+0,
        1.12058e+0, 1.45094e+0, 1.83823e+0, 2.44356e+0, 3.05868e+0,
        3.88521e+0, 4.50784e+0, 5.54901e+0, 6.59355e+0, 1.08010e+1,
        1.50831e+1, 2.18432e+1, 3.48084e+1, 5.08408e+1, 9.58110e+1,
        1.63926e+2, 2.59055e+2, 3.80397e+2, 5.24517e+2, 6.87325e+2,
        8.65041e+2, 1.05451e+3};
    pf[11] = detail::bc16_interp(xqhehp, T_K);

    // 12: D (P09, same as H ground state)
    pf[12] = 2.0;

    // 13: HD (HITRAN file, nuclear-spin deg removed: /6)
    pf[13] = detail::pf_interp(tbl.pf_table[13], T_K) / 6.0;

    // 14: D+
    pf[14] = 1.0;

    // 15: HD+ (ExoMol file, nuclear-spin deg removed: /6)
    pf[15] = detail::pf_interp(tbl.pf_table[15], T_K) / 6.0;

    // 16: D- (same as H-)
    pf[16] = 1.0;

    // Li species (17-23) — only populated for zero_metal (N_sp == 23)
    if constexpr (N_sp == 23) {
        // 17: Li (BC16, 32-point)
        static constexpr std::array<double,32> xqLi = {
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00013e+0, 2.00473e+0,
            2.02869e+0, 2.08800e+0, 2.19942e+0, 2.39054e+0, 2.70188e+0,
            3.17983e+0, 3.86752e+0};
        pf[17] = detail::bc16_interp(xqLi, T_K);

        // 18: LiH (BC16, 32-point)
        static constexpr std::array<double,32> xqLiH = {
            1.00000e+0, 1.00000e+0, 1.00001e+0, 1.00007e+0, 1.00247e+0,
            1.04234e+0, 1.14353e+0, 1.36475e+0, 1.79698e+0, 2.25093e+0,
            3.17620e+0, 5.04593e+0, 6.92359e+0, 9.74658e+0, 1.25749e+1,
            1.63530e+1, 1.91923e+1, 2.39405e+1, 2.87248e+1, 4.88043e+1,
            7.16412e+1, 1.13107e+2, 2.04843e+2, 3.28787e+2, 6.95286e+2,
            1.26125e+3, 2.06412e+3, 3.13695e+3, 4.53795e+3, 6.36625e+3,
            8.75776e+3, 1.18704e+4};
        pf[18] = detail::bc16_interp(xqLiH, T_K);

        // 19: Li+
        pf[19] = 1.0;

        // 20: Li- (JANAF-NIST)
        pf[20] = 1.0;

        // 21: LiH+ (ExoMol, nuclear-spin deg removed: /2)
        pf[21] = detail::pf_interp(tbl.pf_table[21], T_K) / 2.0;

        // 22: Li++
        pf[22] = 2.0;

        // 23: Li+++
        pf[23] = 1.0;
    }

    // -----------------------------------------------------------------------
    // Metal-grain species 17..63  (N_sp == 89, metal_grain network)
    // Species 64..89 are grain/ice pseudo-species → default 1.0 from fill().
    // -----------------------------------------------------------------------
    if constexpr (N_sp >= 89) {
        constexpr double h_P = phys::h_P; // erg·s
        constexpr double pi_ = phys::pi;
        constexpr double MHz = 1.0e6;

        // xlg_th = log10(5040/T) for I88 polynomial fits
        double xlg_th  = std::log10(5040.0 / T_K);
        double xlg_th2 = xlg_th*xlg_th;
        double xlg_th3 = xlg_th2*xlg_th;
        double xlg_th4 = xlg_th3*xlg_th;
        double xlg_th5 = xlg_th4*xlg_th;
        double xlg_th6 = xlg_th5*xlg_th;

        // Classical asymmetric-top partition function
        //   pf = sqrt( (kT/h)^3 * pi / (A*B*C) ) * factor
        auto asym_top = [&](double A_MHz, double B_MHz, double C_MHz,
                            double factor = 1.0) {
            double A = A_MHz*MHz, B = B_MHz*MHz, C = C_MHz*MHz;
            double kT_h = xk_B * T_K / h_P;
            return std::sqrt(kT_h*kT_h*kT_h * pi_ / (A*B*C)) * factor;
        };

        // 17: C (BC16)
        static constexpr std::array<double,32> xqc = {
            1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00002e+0, 1.00115e+0,
            1.02680e+0, 1.10379e+0, 1.29317e+0, 1.70020e+0, 2.14254e+0,
            2.99030e+0, 4.30589e+0, 5.19094e+0, 6.04754e+0, 6.59516e+0,
            7.07437e+0, 7.32544e+0, 7.62484e+0, 7.83366e+0, 8.27478e+0,
            8.47391e+0, 8.62742e+0, 8.74962e+0, 8.81441e+0, 8.91124e+0,
            9.03328e+0, 9.19239e+0, 9.37771e+0, 9.57770e+0, 9.78474e+0,
            9.99517e+0, 1.02090e+1};
        pf[17] = detail::bc16_interp(xqc, T_K);

        // 18: C2 (BC16)
        static constexpr std::array<double,32> xqC2 = {
            1.00000e+0, 1.00003e+0, 1.00051e+0, 1.00201e+0, 1.02728e+0,
            1.21957e+0, 1.54110e+0, 2.09649e+0, 3.05116e+0, 4.00914e+0,
            5.92681e+0, 9.76432e+0, 1.36033e+1, 1.93817e+1, 2.53157e+1,
            3.40450e+1, 4.17374e+1, 5.77217e+1, 7.83868e+1, 2.08677e+2,
            4.04331e+2, 8.01464e+2, 1.72802e+3, 3.00284e+3, 6.75852e+3,
            1.24327e+4, 2.03901e+4, 3.09829e+4, 4.45725e+4, 6.15554e+4,
            8.23790e+4, 1.07544e+5};
        pf[18] = detail::bc16_interp(xqC2, T_K);

        // 19: CH (BC16)
        static constexpr std::array<double,32> xqch = {
            1.20000e+1, 1.20000e+1, 1.20000e+1, 1.20000e+1, 1.20000e+1,
            1.20000e+1, 1.20002e+1, 1.20057e+1, 1.20868e+1, 1.23391e+1,
            1.33487e+1, 1.64087e+1, 1.99589e+1, 2.55718e+1, 3.13181e+1,
            3.90674e+1, 4.49138e+1, 5.46956e+1, 6.45085e+1, 1.03985e+2,
            1.44136e+2, 2.07690e+2, 3.30136e+2, 4.82164e+2, 9.04412e+2,
            1.52375e+3, 2.38711e+3, 3.54843e+3, 5.06630e+3, 6.99662e+3,
            9.38682e+3, 1.22732e+4};
        pf[19] = detail::bc16_interp(xqch, T_K);

        // 20: CH2 (CDMS rotational for T<=600, I88 polynomial for T>600)
        if (T_K <= 600.0) {
            pf[20] = asym_top(2211494.0, 253618.1, 215102.2, 1.5);
        } else {
            pf[20] = std::pow(10.0, 5.05057221 - 3.7971851*xlg_th
                              + 0.928883*xlg_th2 + 0.44249*xlg_th3
                              + 0.14807*xlg_th4  - 0.37003*xlg_th5);
        }

        // 21: CH3 (HITRAN file / 8 for T<=500, I88 for T>500)
        if (T_K <= 500.0) {
            pf[21] = detail::pf_interp(tbl.pf_table[21], T_K) / 8.0;
        } else {
            pf[21] = std::pow(10.0, 6.2439699776 - 5.909173844*xlg_th
                              + 1.5818662*xlg_th2 + 0.767593*xlg_th3
                              + 0.449532*xlg_th4  - 1.04287*xlg_th5
                              + 0.264771*xlg_th6);
        }

        // 22: CH4 (HITRAN file / 16 for T<=3500, I88 for T>3500)
        if (T_K <= 3500.0) {
            pf[22] = detail::pf_interp(tbl.pf_table[22], T_K) / 16.0;
        } else {
            pf[22] = std::pow(10.0, 6.942294193 - 8.66259486*xlg_th
                              + 3.218163*xlg_th2 + 0.69805*xlg_th3
                              + 0.92444*xlg_th4  - 1.3878*xlg_th5
                              + 0.1889*xlg_th6);
        }

        // 23: C+ (BC16)
        static constexpr std::array<double,32> xqcp = {
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00001e+0, 2.00044e+0, 2.00913e+0, 2.04177e+0,
            2.19111e+0, 2.64504e+0, 3.08645e+0, 3.60629e+0, 3.98273e+0,
            4.33873e+0, 4.53479e+0, 4.77694e+0, 4.95108e+0, 5.33283e+0,
            5.51120e+0, 5.65121e+0, 5.76395e+0, 5.82163e+0, 5.88018e+0,
            5.90980e+0, 5.92772e+0, 5.94003e+0, 5.94994e+0, 5.95988e+0,
            5.97207e+0, 5.98845e+0};
        pf[23] = detail::bc16_interp(xqcp, T_K);

        // 24: C2+ (BC16)
        static constexpr std::array<double,32> xqc2p = {
            6.00071e+0, 6.00646e+0, 6.03640e+0, 6.08459e+0, 6.41978e+0,
            7.60333e+0, 9.05993e+0, 1.14110e+1, 1.54767e+1, 1.96058e+1,
            2.79266e+1, 4.46440e+1, 6.13923e+1, 8.65365e+1, 1.11695e+2,
            1.45258e+2, 1.70449e+2, 2.12517e+2, 2.54831e+2, 4.31572e+2,
            6.30643e+2, 9.86295e+2, 1.74817e+3, 2.72947e+3, 5.36363e+3,
            8.90825e+3, 1.33809e+4, 1.88008e+4, 2.51901e+4, 3.25762e+4,
            4.09927e+4, 5.04806e+4};
        pf[24] = detail::bc16_interp(xqc2p, T_K);

        // 25: CH+ (BC16)
        static constexpr std::array<double,32> xqchp = {
            1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0,
            1.00099e+0, 1.00977e+0, 1.05449e+0, 1.20886e+0, 1.41647e+0,
            1.88153e+0, 2.85679e+0, 3.84577e+0, 5.33613e+0, 6.82949e+0,
            8.82249e+0, 1.03179e+1, 1.28110e+1, 1.53046e+1, 2.52948e+1,
            3.54289e+1, 5.14755e+1, 8.22832e+1, 1.20289e+2, 2.30672e+2,
            4.15819e+2, 7.12556e+2, 1.15545e+3, 1.77578e+3, 2.60123e+3,
            3.65584e+3, 4.96011e+3};
        pf[25] = detail::bc16_interp(xqchp, T_K);

        // 26: CH2+ (asymmetric top)
        pf[26] = asym_top(2075377.0, 236712.0, 212478.0);

        // 27: CH3+ (symmetric top / 6)
        pf[27] = asym_top(279591.0, 279591.0, 139796.0, 1.0/6.0);

        // 28: CH4+ (asymmetric top / 2)
        pf[28] = asym_top(199364.0, 154791.0, 113699.0, 0.5);

        // 29: CH5+ (asymmetric top)
        pf[29] = asym_top(113321.5, 113921.1, 114820.5);

        // 30: O (BC16)
        static constexpr std::array<double,32> xqo = {
            5.00000e+0, 5.00000e+0, 5.00000e+0, 5.00000e+0, 5.00000e+0,
            5.00000e+0, 5.00000e+0, 5.00000e+0, 5.00000e+0, 5.00003e+0,
            5.00154e+0, 5.03305e+0, 5.12544e+0, 5.34602e+0, 5.60172e+0,
            5.93258e+0, 6.15640e+0, 6.47757e+0, 6.74123e+0, 7.42310e+0,
            7.79423e+0, 8.11056e+0, 8.38189e+0, 8.52662e+0, 8.68009e+0,
            8.77224e+0, 8.85532e+0, 8.94697e+0, 9.05112e+0, 9.16637e+0,
            9.28982e+0, 9.41864e+0};
        pf[30] = detail::bc16_interp(xqo, T_K);

        // 31: O2 (BC16)
        static constexpr std::array<double,32> xqo2 = {
            9.00000e+0, 9.00000e+0, 9.00011e+0, 9.00068e+0, 9.02128e+0,
            9.33579e+0, 1.01023e+1, 1.17557e+1, 1.50103e+1, 1.84654e+1,
            2.55565e+1, 3.99407e+1, 5.43989e+1, 7.61289e+1, 9.78808e+1,
            1.26901e+2, 1.48678e+2, 1.85004e+2, 2.21427e+2, 3.70966e+2,
            5.34408e+2, 8.19779e+2, 1.42263e+3, 2.19951e+3, 4.34427e+3,
            7.41870e+3, 1.16070e+4, 1.71056e+4, 2.41444e+4, 3.30098e+4,
            4.40540e+4, 5.76869e+4};
        pf[31] = detail::bc16_interp(xqo2, T_K);

        // 32: OH (BC16)
        static constexpr std::array<double,32> xqoh = {
            1.20000e+1, 1.20000e+1, 1.20000e+1, 1.20000e+1, 1.20000e+1,
            1.20000e+1, 1.20000e+1, 1.20005e+1, 1.20163e+1, 1.20966e+1,
            1.25753e+1, 1.45073e+1, 1.70195e+1, 2.11581e+1, 2.54704e+1,
            3.13320e+1, 3.57710e+1, 4.32133e+1, 5.06892e+1, 8.07652e+1,
            1.11070e+2, 1.57530e+2, 2.41294e+2, 3.37500e+2, 5.77700e+2,
            8.91168e+2, 1.28807e+3, 1.78078e+3, 2.38461e+3, 3.11647e+3,
            3.99262e+3, 5.02698e+3};
        pf[32] = detail::bc16_interp(xqoh, T_K);

        // 33: CO (BC16)
        static constexpr std::array<double,32> xqco = {
            1.01186e+0, 1.04254e+0, 1.11606e+0, 1.18989e+0, 1.49427e+0,
            2.18214e+0, 2.89218e+0, 3.96754e+0, 5.76822e+0, 7.57241e+0,
            1.11842e+1, 1.84123e+1, 2.56425e+1, 3.64899e+1, 4.73391e+1,
            6.18073e+1, 7.26603e+1, 9.07524e+1, 1.08852e+2, 1.81659e+2,
            2.56937e+2, 3.80242e+2, 6.25544e+2, 9.28195e+2, 1.71706e+3,
            2.75995e+3, 4.06593e+3, 5.64445e+3, 7.50751e+3, 9.67381e+3,
            1.21758e+4, 1.50689e+4};
        pf[33] = detail::bc16_interp(xqco, T_K);

        // 34: H2O (Vidler & Tennyson 2000 polynomial in log10(T))
        {
            double xlgT  = std::log10(T_K);
            double xlgT2 = xlgT*xlgT, xlgT3 = xlgT2*xlgT,
                   xlgT4 = xlgT3*xlgT, xlgT5 = xlgT4*xlgT,
                   xlgT6 = xlgT5*xlgT;
            double xlgq = -14.0874691574179 + 37.9243248539882*xlgT
                          - 42.6817978731789*xlgT2 + 25.3302448517916*xlgT3
                          - 8.10851262935532*xlgT4 +  1.33106871720535*xlgT5
                          - 0.0872981051095757*xlgT6;
            pf[34] = std::max(0.25, std::pow(10.0, xlgq));
        }

        // 35: HCO (JPL rotational formula for T<=500, I88 for T>500)
        if (T_K <= 500.0) {
            pf[35] = asym_top(670485.83, 44788.0, 41930.4, 2.0);
        } else {
            pf[35] = std::pow(10.0, 6.298781639 - 3.85672804*xlg_th
                              + 0.8551678*xlg_th2 + 0.321901*xlg_th3
                              + 0.020274*xlg_th4  + 0.15254*xlg_th5
                              - 0.25298*xlg_th6);
        }

        // 36: HO2 (HITRAN file / 2)
        pf[36] = detail::pf_interp(tbl.pf_table[36], T_K) / 2.0;

        // 37: CO2 (HITRAN file for T<=5000, I88 for T>5000)
        if (T_K <= 5000.0) {
            pf[37] = detail::pf_interp(tbl.pf_table[37], T_K);
        } else {
            pf[37] = std::pow(10.0, 6.01081285 - 4.438833*xlg_th
                              + 0.840462*xlg_th2 + 0.2945*xlg_th3
                              + 0.3694*xlg_th4   - 0.273*xlg_th5);
        }

        // 38: H2CO (HITRAN file / 4)
        pf[38] = detail::pf_interp(tbl.pf_table[38], T_K) / 4.0;

        // 39: H2O2 (HITRAN file / 4)
        pf[39] = detail::pf_interp(tbl.pf_table[39], T_K) / 4.0;

        // 40: O+ (BC16)
        static constexpr std::array<double,32> xqop = {
            4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0,
            4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0,
            4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0,
            4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0,
            4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00003e+0,
            4.00065e+0, 4.00451e+0, 4.01649e+0, 4.04186e+0, 4.08460e+0,
            4.14679e+0, 4.22885e+0};
        pf[40] = detail::bc16_interp(xqop, T_K);

        // 41: O2+ (BC16)
        static constexpr std::array<double,32> xqo2p = {
            3.00032e+0, 3.00294e+0, 3.01694e+0, 3.03979e+0, 3.20139e+0,
            3.78014e+0, 4.49647e+0, 5.65503e+0, 7.66052e+0, 9.69811e+0,
            1.38059e+1, 2.21326e+1, 3.08526e+1, 4.52466e+1, 6.13856e+1,
            8.52491e+1, 1.04528e+2, 1.38583e+2, 1.74381e+2, 3.28202e+2,
            4.96506e+2, 7.83125e+2, 1.37289e+3, 2.11772e+3, 4.10188e+3,
            6.77817e+3, 1.01888e+4, 1.44008e+4, 1.95417e+4, 2.58386e+4,
            3.36476e+4, 4.34606e+4};
        pf[41] = detail::bc16_interp(xqo2p, T_K);

        // 42: OH+ (BC16)
        static constexpr std::array<double,32> xqohp = {
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00047e+0, 2.00703e+0, 2.05326e+0, 2.25799e+0, 2.57365e+0,
            3.33233e+0, 4.97055e+0, 6.64415e+0, 9.17264e+0, 1.17099e+1,
            1.50998e+1, 1.76456e+1, 2.18936e+1, 2.61470e+1, 4.32173e+1,
            6.04874e+1, 8.74853e+1, 1.38201e+2, 1.98838e+2, 3.56403e+2,
            5.71245e+2, 8.57404e+2, 1.23389e+3, 1.72209e+3, 2.34214e+3,
            3.11051e+3, 4.03909e+3};
        pf[42] = detail::bc16_interp(xqohp, T_K);

        // 43: CO+ (BC16)
        static constexpr std::array<double,32> xqcop = {
            2.02087e+0, 2.07708e+0, 2.21516e+0, 2.35588e+0, 2.94397e+0,
            4.28449e+0, 5.67104e+0, 7.77209e+0, 1.12911e+1, 1.48173e+1,
            2.18769e+1, 3.60049e+1, 5.01372e+1, 7.13398e+1, 9.25460e+1,
            1.20826e+2, 1.42039e+2, 1.77403e+2, 2.12779e+2, 3.55006e+2,
            5.01729e+2, 7.41213e+2, 1.21592e+3, 1.80037e+3, 3.32256e+3,
            5.34428e+3, 7.92672e+3, 1.11813e+4, 1.52681e+4, 2.03854e+4,
            2.67626e+4, 3.46562e+4};
        pf[43] = detail::bc16_interp(xqcop, T_K);

        // 44: H2O+ (asymmetric top)
        pf[44] = asym_top(870581.0, 372365.0, 253880.0);

        // 45: HCO+ (linear rotor for T<=350, I88 for T>350)
        if (T_K <= 350.0) {
            double B = 44594.43 * MHz;
            pf[45] = (xk_B * T_K / h_P) / B;
        } else {
            pf[45] = std::pow(10.0, 5.453394934 - 4.14568927*xlg_th
                              + 0.8632023*xlg_th2 + 0.482875*xlg_th3
                              + 0.14863*xlg_th4   - 0.281734*xlg_th5);
        }

        // 46: O2H+ (asymmetric top × 3)
        pf[46] = asym_top(663778.0, 38723.0, 36588.0, 3.0);

        // 47: H3O+ (symmetric top / 3 for T<=500, I88 for T>500)
        if (T_K <= 500.0) {
            pf[47] = asym_top(334405.0, 334405.0, 184725.0, 1.0/3.0);
        } else {
            pf[47] = std::pow(10.0, 5.6188225644 - 5.599410492*xlg_th
                              + 1.8194835*xlg_th2 + 0.774831*xlg_th3
                              + 0.404714*xlg_th4  - 1.53688*xlg_th5
                              + 0.655507*xlg_th6);
        }

        // 48: H2CO+ (asymmetric top)
        pf[48] = asym_top(266035.83, 40232.15, 34416.17);

        // 49: HOCO+ (asymmetric top)
        pf[49] = asym_top(789947.8, 10773.733, 10609.431);

        // 50: H2COH+ (asymmetric top)
        pf[50] = asym_top(197581.56, 34350.55, 29172.65);

        // 51: Li (BC16)
        static constexpr std::array<double,32> xqLi_mg = {
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00013e+0, 2.00473e+0,
            2.02869e+0, 2.08800e+0, 2.19942e+0, 2.39054e+0, 2.70188e+0,
            3.17983e+0, 3.86752e+0};
        pf[51] = detail::bc16_interp(xqLi_mg, T_K);

        // 52: LiH (BC16)
        static constexpr std::array<double,32> xqLiH_mg = {
            1.00000e+0, 1.00000e+0, 1.00001e+0, 1.00007e+0, 1.00247e+0,
            1.04234e+0, 1.14353e+0, 1.36475e+0, 1.79698e+0, 2.25093e+0,
            3.17620e+0, 5.04593e+0, 6.92359e+0, 9.74658e+0, 1.25749e+1,
            1.63530e+1, 1.91923e+1, 2.39405e+1, 2.87248e+1, 4.88043e+1,
            7.16412e+1, 1.13107e+2, 2.04843e+2, 3.28787e+2, 6.95286e+2,
            1.26125e+3, 2.06412e+3, 3.13695e+3, 4.53795e+3, 6.36625e+3,
            8.75776e+3, 1.18704e+4};
        pf[52] = detail::bc16_interp(xqLiH_mg, T_K);

        // 53: Li+  = 1.0 (fill default)
        // 54: Li-  = 1.0 (fill default)

        // 55: LiH+ (ExoMol file, nuclear-spin deg removed: /2)
        pf[55] = detail::pf_interp(tbl.pf_table[55], T_K) / 2.0;

        // 56: Li++
        pf[56] = 2.0;

        // 57: Li+++ = 1.0 (fill default)

        // 58: K (BC16)
        static constexpr std::array<double,32> xqK = {
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00002e+0, 2.00052e+0, 2.01222e+0,
            2.06702e+0, 2.22592e+0, 2.61039e+0, 3.39777e+0, 4.77353e+0,
            6.88524e+0, 9.82105e+0};
        pf[58] = detail::bc16_interp(xqK, T_K);

        // 59: K+  = 1.0 (fill default)

        // 60: Na (BC16)
        static constexpr std::array<double,32> xqNa = {
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00003e+0, 2.00178e+0,
            2.01488e+0, 2.06447e+0, 2.21609e+0, 2.60177e+0, 3.40984e+0,
            4.84529e+0, 7.08960e+0};
        pf[60] = detail::bc16_interp(xqNa, T_K);

        // 61: Na+ = 1.0 (fill default)

        // 62: Mg (BC16)
        static constexpr std::array<double,32> xqMg = {
            1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0,
            1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0,
            1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0,
            1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0,
            1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00025e+0,
            1.00344e+0, 1.01678e+0, 1.04919e+0, 1.11004e+0, 1.21285e+0,
            1.37994e+0, 1.64434e+0};
        pf[62] = detail::bc16_interp(xqMg, T_K);

        // 63: Mg+ (BC16)
        static constexpr std::array<double,32> xqMgp = {
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
            2.00002e+0, 2.00021e+0, 2.00114e+0, 2.00389e+0, 2.00976e+0,
            2.02002e+0, 2.03571e+0};
        pf[63] = detail::bc16_interp(xqMgp, T_K);

        // 64..89: grain/ice pseudo-species = 1.0 (fill default)
    }
}

} // namespace chemistry
