// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// reaction_metal_grain_charge.h
//
// MRN size-averaged gas-grain and grain-grain charge transfer rates.
//
// Port of:
//   fortran/metal_grain/subs/reaction_metal_grain.f
//     SUBROUTINE react_grain    (L3476–3977)  — MRN quadrature driver
//     SUBROUTINE rate_gsneu     (L3981–3995)  — neutral grain collision
//     SUBROUTINE rate_gsatt     (L3998–4014)  — attractive ion-grain
//     SUBROUTINE rate_gsrep     (L4017–4036)  — repulsive ion-grain
//     SUBROUTINE rate_TE        (L4040–4055)  — thermionic emission
//     SUBROUTINE rate_gg        (L4059–4084)  — grain-grain collision
//
// Output: react_grain_rates() fills xk_charge[0..149]  (150 entries).
// These map to tbl.react[670..819] in the main reaction table
// (Fortran xk(841..990)).
//
// Grain size distribution: broken MRN power law
//   n(a) ∝ a^{-3.5}  for a ∈ [5e-7, 1e-4] cm   (Nint1 = 15 midpoints)
//   n(a) ∝ a^{-5.5}  for a ∈ [1e-4,  5e-4] cm   (Nint2 =  5 midpoints)
//
// All units CGS.
// ---------------------------------------------------------------------------
#include <cmath>
#include <array>
#include "species.h"
#include "state.h"

namespace chemistry {
namespace detail {

// ---------------------------------------------------------------------------
// rate_gsneu — ion/anion + neutral grain  (L3981–3995)
//   Draine & Sutin 1987, eq. for J_neu
//   τ = a·kB·T / (qs·e)²
//   J = 1 + √(π/2τ)
//   k = π a² v_th J
// ---------------------------------------------------------------------------
inline double rate_gsneu(double T_K, double xm_gs, double xq_gs, double a_gr)
{
    constexpr double xk_B = phys::xk_B;
    constexpr double pi   = phys::pi;
    constexpr double xm_p = phys::xm_p;
    constexpr double qe   = phys::qe;

    double tau_gr = a_gr * xk_B * T_K / (xq_gs * qe * (xq_gs * qe));
    double J_neu  = 1.0 + std::sqrt(pi / (2.0 * tau_gr));
    return pi * a_gr * a_gr
           * std::sqrt(8.0 * xk_B * T_K / (pi * xm_gs * xm_p))
           * J_neu;
}

// ---------------------------------------------------------------------------
// rate_gsatt — ion + oppositely charged grain  (L3998–4014)
//   J_att = (1 - ν/τ)(1 + √(2/(τ - 2ν)))    ν = q_gr/q_gs
// ---------------------------------------------------------------------------
inline double rate_gsatt(double T_K, double xm_gs, double xq_gs,
                          double a_gr, double xq_gr)
{
    constexpr double xk_B = phys::xk_B;
    constexpr double pi   = phys::pi;
    constexpr double xm_p = phys::xm_p;
    constexpr double qe   = phys::qe;

    double tau_gr = a_gr * xk_B * T_K / (xq_gs * qe * (xq_gs * qe));
    double nu_gr  = xq_gr / xq_gs;
    double J_att  = (1.0 - nu_gr / tau_gr)
                  * (1.0 + std::sqrt(2.0 / (tau_gr - 2.0 * nu_gr)));
    return pi * a_gr * a_gr
           * std::sqrt(8.0 * xk_B * T_K / (pi * xm_gs * xm_p))
           * J_att;
}

// ---------------------------------------------------------------------------
// rate_gsrep — ion + same-sign charged grain  (L4017–4036)
//   θ_ν = ν / (1 + 1/√ν)
//   J_rep = (1 + 1/√(4τ + 3ν))² · exp(−θ_ν/τ)
// ---------------------------------------------------------------------------
inline double rate_gsrep(double T_K, double xm_gs, double xq_gs,
                          double a_gr, double xq_gr)
{
    constexpr double xk_B = phys::xk_B;
    constexpr double pi   = phys::pi;
    constexpr double xm_p = phys::xm_p;
    constexpr double qe   = phys::qe;

    double tau_gr   = a_gr * xk_B * T_K / (xq_gs * qe * (xq_gs * qe));
    double nu_gr    = xq_gr / xq_gs;
    double theta_nu = nu_gr / (1.0 + 1.0 / std::sqrt(nu_gr));
    double tmp      = 1.0 + 1.0 / std::sqrt(4.0 * tau_gr + 3.0 * nu_gr);
    double J_rep    = tmp * tmp * std::exp(-theta_nu / tau_gr);
    return pi * a_gr * a_gr
           * std::sqrt(8.0 * xk_B * T_K / (pi * xm_gs * xm_p))
           * J_rep;
}

// ---------------------------------------------------------------------------
// rate_TE — thermionic electron emission from grain  (L4040–4055)
//   Richardson-Dushman with Richardson factor λ_R = 0.5, WF = 5 eV
//   j_TE = λ_R (4π m_e / h³) (kB T_gr)² exp(−(WF + q_gr e²/a) / kB T_gr)
//   k_TE = 4π a² j_TE   [s^{-1}]
// ---------------------------------------------------------------------------
inline double rate_TE(double T_gr_K, double a_gr, double xq_gr)
{
    constexpr double xk_B   = phys::xk_B;
    constexpr double pi     = phys::pi;
    constexpr double xm_e   = phys::xm_e;
    constexpr double qe   = phys::qe;
    constexpr double h_p    = phys::h_P;
    constexpr double lmbd_R = 0.5;
    constexpr double WF     = 5.0 * phys::eV_to_erg;   // 5 eV in erg

    double h3    = h_p * h_p * h_p;
    double kT    = xk_B * T_gr_K;
    double j_TE  = lmbd_R * (4.0 * pi * xm_e / h3)
                 * kT * kT
                 * std::exp(-(WF + xq_gr * qe * qe / a_gr) / kT);
    return 4.0 * pi * a_gr * a_gr * j_TE;
}

// ---------------------------------------------------------------------------
// rate_gg — grain-grain charge transfer  (L4059–4084)
//   Attractive geometry (q1·q2 < 0).
//   Symmetrizes J by averaging the two Coulomb focusing factors.
//   rho_gr = 3 g/cm³ (silicate bulk density)
// ---------------------------------------------------------------------------
inline double rate_gg(double T_K, double a1, double a2, double q1, double q2)
{
    constexpr double xk_B   = phys::xk_B;
    constexpr double pi     = phys::pi;
    constexpr double qe   = phys::qe;
    constexpr double rho_gr = 3.0;   // [g/cm³]

    double a_sum = a1 + a2;
    double m1    = (4.0 * pi / 3.0) * a1 * a1 * a1 * rho_gr;
    double m2    = (4.0 * pi / 3.0) * a2 * a2 * a2 * rho_gr;
    double mu    = m1 * m2 / (m1 + m2);   // reduced mass

    // J from q1 perspective: grain 2 attracts grain 1
    double tau1 = a_sum * xk_B * T_K / (q1 * qe * (q1 * qe));
    double nu1  = q2 / q1;
    double J1   = (1.0 - nu1 / tau1)
                * (1.0 + std::sqrt(2.0 / (tau1 - 2.0 * nu1)));

    // J from q2 perspective
    double tau2 = a_sum * xk_B * T_K / (q2 * qe * (q2 * qe));
    double nu2  = q1 / q2;
    double J2   = (1.0 - nu2 / tau2)
                * (1.0 + std::sqrt(2.0 / (tau2 - 2.0 * nu2)));

    return pi * a_sum * a_sum
           * std::sqrt(8.0 * xk_B * T_K / (pi * mu))
           * (J1 + J2) / 2.0;
}

} // namespace detail

// ---------------------------------------------------------------------------
// react_grain_rates
//
//   Computes 150 MRN size-averaged charge-transfer rate coefficients and
//   stores them in xk_charge[0..149].
//
//   Ordering (matches Fortran react_grain local num = 1..150):
//     [  0.. 29]  X¹⁺ + gr⁰   (30 species, rate_gsneu)
//     [ 30.. 31]  X²⁺ + gr⁰   ( 2 species, rate_gsneu)
//     [ 32.. 35]  X⁻  + gr⁰   ( 4 species, rate_gsneu)
//     [ 36.. 65]  X¹⁺ + gr⁺   (30 species, rate_gsrep)
//     [ 66.. 69]  X⁻  + gr⁺   ( 4 species, rate_gsatt)
//     [ 70.. 73]  X⁻  + gr²⁺  ( 4 species, rate_gsatt)
//     [ 74..103]  X¹⁺ + gr⁻   (30 species, rate_gsatt)
//     [104..105]  X²⁺ + gr⁻   ( 2 species, rate_gsatt)
//     [106..109]  X⁻  + gr⁻   ( 4 species, rate_gsrep)
//     [110..139]  X¹⁺ + gr²⁻  (30 species, rate_gsatt)
//     [140..141]  X²⁺ + gr²⁻  ( 2 species, rate_gsatt)
//     [142..145]  thermionic emission  (xq_gr = −1,0,+1,+2)
//     [146]       gr⁺  + gr⁻   (rate_gg, q1=+1, q2=−1)
//     [147]       gr⁺  + gr²⁻  (rate_gg, q1=+1, q2=−2)
//     [148]       gr²⁺ + gr²⁻  (rate_gg, q1=+2, q2=−2)
//     [149]       gr²⁺ + gr⁻   (rate_gg, q1=+2, q2=−1)
//
//   Parameters
//     T_K    : gas temperature [K]
//     T_gr_K : grain temperature [K]  (used only for thermionic emission)
//     xk_charge : output array, length 150
// ---------------------------------------------------------------------------
inline void react_grain_rates(double T_K, double T_gr_K,
                               std::array<double, 150>& xk_charge)
{
    using namespace detail;

    // ── Species mass table (Fortran mass(1..63), C++ 0-based) ──────────────
    // Covers the first 63 gas-phase species of the metal_grain network.
    static constexpr double mass[63] = {
        1.00783,  2.01565,  5.5e-4,   1.00728,  2.0151,   3.02293,  // 1-6
        1.00837,  4.0026,   4.00205,  4.0,      5.00988,  2.0141,   // 7-12
        3.02204,  2.01355,  3.02149,  2.01465,  12.0,     24.0,     // 13-18
        13.00783, 14.01565, 15.02348, 16.0313,  11.99945, 23.99945, // 19-24
        13.00728, 14.0151,  15.02293, 16.03075, 17.03858, 15.99491, // 25-30
        31.98983, 17.00274, 27.99491, 18.01056, 29.00274, 32.99765, // 31-36
        43.98983, 30.01056, 34.00548, 15.99437, 31.98928, 17.00219, // 37-42
        27.99437, 18.01002, 29.00219, 32.99711, 19.01784, 30.01002, // 43-48
        44.99711, 31.01784,  6.941,    7.949,    6.94,     6.942,   // 49-54
         7.949,    6.94,     6.94,    39.0983,  39.09775, 22.98977, // 55-60
        22.98922, 24.305,   24.30445                                 // 61-63
    };

    // ── Charged species index lists (C++ 0-based into mass[]) ─────────────
    // Fortran x1p(1..30) = {4,5,6,9,11,14,15,23,24,25,26,27,28,29,
    //                        40,41,42,43,44,45,46,47,48,49,50,53,55,59,61,63}
    static constexpr int x1p[30] = {
         3,  4,  5,  8, 10, 13, 14, 22, 23, 24,
        25, 26, 27, 28, 39, 40, 41, 42, 43, 44,
        45, 46, 47, 48, 49, 52, 54, 58, 60, 62
    };
    // Fortran x2p(1..2) = {10, 56}
    static constexpr int x2p[2] = { 9, 55 };
    // Fortran x1m(1..4) = {3, 7, 16, 54}
    static constexpr int x1m[4] = { 2, 6, 15, 53 };

    // ── MRN grain size quadrature parameters ──────────────────────────────
    constexpr double a_min = chemistry::metal_grain::mrn_a_min;
    constexpr double a_mid = chemistry::metal_grain::mrn_a_mid;
    constexpr double a_max = chemistry::metal_grain::mrn_a_max;
    constexpr double ind1  = chemistry::metal_grain::mrn_ind1;
    constexpr double ind2  = chemistry::metal_grain::mrn_ind2;
    constexpr int    Nint1 = chemistry::metal_grain::mrn_nint1;
    constexpr int    Nint2 = chemistry::metal_grain::mrn_nint2;

    const double x_min  = std::log(a_min);
    const double x_mid  = std::log(a_mid);   // NOTE: Fortran typo (line 3520 overwrites
    const double x_max  = std::log(a_max);   //   x_max twice); intended value is log(a_mid)
    const double d1_lna = (x_mid - x_min) / double(Nint1);
    const double d2_lna = (x_max - x_mid) / double(Nint2);

    // Analytical normalization ∫ n(a)/n_ref da  (power-law integrals)
    const double sum_ana =
          std::pow(a_mid, ind1)
        * (std::pow(a_mid, 1.0 - ind1) - std::pow(a_min, 1.0 - ind1))
        / (1.0 - ind1)
        + std::pow(a_mid, ind2)
        * (std::pow(a_max, 1.0 - ind2) - std::pow(a_mid, 1.0 - ind2))
        / (1.0 - ind2);

    // ── 1-D MRN average: ∫ rate(a) · n(a) da / ∫ n(a) da ─────────────────
    // Uses midpoint rule over the two power-law segments.
    auto mrn_avg1 = [&](auto rate_fn) -> double {
        double val = 0.0;
        for (int itn = 0; itn < Nint1; ++itn) {
            double a = std::exp((x_min + 0.5 * d1_lna) + itn * d1_lna);
            val += rate_fn(a) * std::pow(a, 1.0 - ind1) * d1_lna
                              * std::pow(a_mid, ind1);
        }
        for (int itn = 0; itn < Nint2; ++itn) {
            double a = std::exp((x_mid + 0.5 * d2_lna) + itn * d2_lna);
            val += rate_fn(a) * std::pow(a, 1.0 - ind2) * d2_lna
                              * std::pow(a_mid, ind2);
        }
        return val / sum_ana;
    };

    // ── 2-D MRN average for grain-grain reactions ─────────────────────────
    // ∫∫ rate(a1,a2) · n(a1) n(a2) da1 da2 / [∫ n da]²
    auto mrn_avg2 = [&](auto rate_fn2) -> double {
        double val = 0.0;
        // Build the four cross-segment pairs: (seg1,seg1), (seg1,seg2),
        // (seg2,seg1), (seg2,seg2)
        auto seg = [&](int seg_id, int idx) -> double {
            // Returns midpoint log-radius for segment seg_id (0=small, 1=large)
            return (seg_id == 0)
                ? (x_min + 0.5 * d1_lna) + idx * d1_lna
                : (x_mid + 0.5 * d2_lna) + idx * d2_lna;
        };
        auto wt = [&](int seg_id, int idx) -> double {
            double a   = std::exp(seg(seg_id, idx));
            double ind = (seg_id == 0) ? ind1 : ind2;
            double dl  = (seg_id == 0) ? d1_lna : d2_lna;
            return std::pow(a, 1.0 - ind) * dl * std::pow(a_mid, ind);
        };
        int N[2] = { Nint1, Nint2 };
        for (int s1 = 0; s1 < 2; ++s1)
            for (int s2 = 0; s2 < 2; ++s2)
                for (int i = 0; i < N[s1]; ++i)
                    for (int j = 0; j < N[s2]; ++j) {
                        double a1 = std::exp(seg(s1, i));
                        double a2 = std::exp(seg(s2, j));
                        val += rate_fn2(a1, a2) * wt(s1, i) * wt(s2, j);
                    }
        return val / (sum_ana * sum_ana);
    };

    // ── Fill xk_charge[0..149] ─────────────────────────────────────────────
    xk_charge.fill(0.0);
    int num = 0;

    // ── (1) Neutral grain: all ions/anions with gr⁰ ─────────────────────
    // X¹⁺ + gr⁰  (30)
    for (int i = 0; i < 30; ++i) {
        const double m = mass[x1p[i]], q = 1.0;
        xk_charge[num++] = mrn_avg1([&](double a){ return rate_gsneu(T_K,m,q,a); });
    }
    // X²⁺ + gr⁰  (2)
    for (int i = 0; i < 2; ++i) {
        const double m = mass[x2p[i]], q = 2.0;
        xk_charge[num++] = mrn_avg1([&](double a){ return rate_gsneu(T_K,m,q,a); });
    }
    // X⁻  + gr⁰  (4)
    for (int i = 0; i < 4; ++i) {
        const double m = mass[x1m[i]], q = -1.0;
        xk_charge[num++] = mrn_avg1([&](double a){ return rate_gsneu(T_K,m,q,a); });
    }

    // ── (2) Positive grain: gr⁺ ─────────────────────────────────────────
    // X¹⁺ + gr⁺  repulsion (30)
    for (int i = 0; i < 30; ++i) {
        const double m = mass[x1p[i]], q = 1.0, qg = 1.0;
        xk_charge[num++] = mrn_avg1([&](double a){ return rate_gsrep(T_K,m,q,a,qg); });
    }
    // X⁻  + gr⁺  attraction (4)
    for (int i = 0; i < 4; ++i) {
        const double m = mass[x1m[i]], q = -1.0, qg = 1.0;
        xk_charge[num++] = mrn_avg1([&](double a){ return rate_gsatt(T_K,m,q,a,qg); });
    }
    // X⁻  + gr²⁺ attraction (4)
    for (int i = 0; i < 4; ++i) {
        const double m = mass[x1m[i]], q = -1.0, qg = 2.0;
        xk_charge[num++] = mrn_avg1([&](double a){ return rate_gsatt(T_K,m,q,a,qg); });
    }

    // ── (3) Negative grain: gr⁻ ─────────────────────────────────────────
    // X¹⁺ + gr⁻  attraction (30)
    for (int i = 0; i < 30; ++i) {
        const double m = mass[x1p[i]], q = 1.0, qg = -1.0;
        xk_charge[num++] = mrn_avg1([&](double a){ return rate_gsatt(T_K,m,q,a,qg); });
    }
    // X²⁺ + gr⁻  attraction (2)
    for (int i = 0; i < 2; ++i) {
        const double m = mass[x2p[i]], q = 2.0, qg = -1.0;
        xk_charge[num++] = mrn_avg1([&](double a){ return rate_gsatt(T_K,m,q,a,qg); });
    }
    // X⁻  + gr⁻  repulsion (4)
    for (int i = 0; i < 4; ++i) {
        const double m = mass[x1m[i]], q = -1.0, qg = -1.0;
        xk_charge[num++] = mrn_avg1([&](double a){ return rate_gsrep(T_K,m,q,a,qg); });
    }

    // ── (4) Doubly negative grain: gr²⁻ ────────────────────────────────
    // X¹⁺ + gr²⁻ attraction (30)
    for (int i = 0; i < 30; ++i) {
        const double m = mass[x1p[i]], q = 1.0, qg = -2.0;
        xk_charge[num++] = mrn_avg1([&](double a){ return rate_gsatt(T_K,m,q,a,qg); });
    }
    // X²⁺ + gr²⁻ attraction (2)
    for (int i = 0; i < 2; ++i) {
        const double m = mass[x2p[i]], q = 2.0, qg = -2.0;
        xk_charge[num++] = mrn_avg1([&](double a){ return rate_gsatt(T_K,m,q,a,qg); });
    }

    // ── (5) Thermionic emission: gr^{q} → gr^{q+1} + e⁻  (4 entries) ─────
    // xq_gr cycles: −1, 0, +1, +2
    for (int iq = 0; iq < 4; ++iq) {
        const double xq_gr = -1.0 + iq;
        xk_charge[num++] = mrn_avg1([&](double a){ return rate_TE(T_gr_K, a, xq_gr); });
    }

    // ── (6) Grain-grain charge transfer  (4 entries) ──────────────────────
    // gr⁺  + gr⁻   → gr⁰  + gr⁰
    xk_charge[num++] = mrn_avg2([&](double a1, double a2){
        return rate_gg(T_K, a1, a2, 1.0, -1.0);
    });
    // gr⁺  + gr²⁻  → gr⁰  + gr⁻
    xk_charge[num++] = mrn_avg2([&](double a1, double a2){
        return rate_gg(T_K, a1, a2, 1.0, -2.0);
    });
    // gr²⁺ + gr²⁻  → gr⁺  + gr⁻
    xk_charge[num++] = mrn_avg2([&](double a1, double a2){
        return rate_gg(T_K, a1, a2, 2.0, -2.0);
    });
    // gr²⁺ + gr⁻   → gr⁺  + gr⁰   (last entry: num not incremented, matching Fortran)
    xk_charge[num] = mrn_avg2([&](double a1, double a2){
        return rate_gg(T_K, a1, a2, 2.0, -1.0);
    });
    // num == 149 here; xk_charge[0..149] fully populated (150 entries total)
}

} // namespace chemistry
