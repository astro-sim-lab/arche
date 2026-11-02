// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// reaction_metal_grain_charge.h
//
// MRN size-averaged gas-grain and grain-grain charge transfer rates.
//
// Computes, by MRN quadrature, the gas-grain and grain-grain charge-transfer
// rates: neutral grain collision, attractive and repulsive ion-grain
// collisions, thermionic emission, and grain-grain collision.
//
// Output: react_grain_rates() fills k_charge[0..149]  (150 entries).
// These map to tbl.react[670..819] in the main reaction table
// (1-based reaction ids 841..990).
//
// Grain size distribution: broken MRN power law
//   n(a) ∝ a^{-3.5}  for a ∈ [5e-7, 1e-4] cm   (Nint1 = 15 midpoints)
//   n(a) ∝ a^{-5.5}  for a ∈ [1e-4,  5e-4] cm   (Nint2 =  5 midpoints)
//
// All units CGS.
// ---------------------------------------------------------------------------
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "core/species_index.h"
#include "core/state.h"

namespace arche {
namespace detail {

// ---------------------------------------------------------------------------
// rate_gsneu — ion/anion + neutral grain
//   Draine & Sutin 1987, eq. for J_neu
//   τ = a·kB·T / (qs·e)²
//   J = 1 + √(π/2τ)
//   k = π a² v_th J
// ---------------------------------------------------------------------------
inline double rate_gsneu(double T_K, double m_gs, double q_gs, double a_gr) {
  constexpr double k_B = phys::k_B;
  constexpr double pi = phys::pi;
  constexpr double m_p = phys::m_p;
  constexpr double qe = phys::qe;

  double tau_gr = a_gr * k_B * T_K / (q_gs * qe * (q_gs * qe));
  double J_neu = 1.0 + std::sqrt(pi / (2.0 * tau_gr));
  return pi * a_gr * a_gr * std::sqrt(8.0 * k_B * T_K / (pi * m_gs * m_p)) *
         J_neu;
}

// ---------------------------------------------------------------------------
// rate_gsatt — ion + oppositely charged grain
//   J_att = (1 - ν/τ)(1 + √(2/(τ - 2ν)))    ν = q_gr/q_gs
// ---------------------------------------------------------------------------
inline double rate_gsatt(double T_K, double m_gs, double q_gs, double a_gr,
                         double q_gr) {
  constexpr double k_B = phys::k_B;
  constexpr double pi = phys::pi;
  constexpr double m_p = phys::m_p;
  constexpr double qe = phys::qe;

  double tau_gr = a_gr * k_B * T_K / (q_gs * qe * (q_gs * qe));
  double nu_gr = q_gr / q_gs;
  double J_att =
      (1.0 - nu_gr / tau_gr) * (1.0 + std::sqrt(2.0 / (tau_gr - 2.0 * nu_gr)));
  return pi * a_gr * a_gr * std::sqrt(8.0 * k_B * T_K / (pi * m_gs * m_p)) *
         J_att;
}

// ---------------------------------------------------------------------------
// rate_gsrep — ion + same-sign charged grain
//   θ_ν = ν / (1 + 1/√ν)
//   J_rep = (1 + 1/√(4τ + 3ν))² · exp(−θ_ν/τ)
// ---------------------------------------------------------------------------
inline double rate_gsrep(double T_K, double m_gs, double q_gs, double a_gr,
                         double q_gr) {
  constexpr double k_B = phys::k_B;
  constexpr double pi = phys::pi;
  constexpr double m_p = phys::m_p;
  constexpr double qe = phys::qe;

  double tau_gr = a_gr * k_B * T_K / (q_gs * qe * (q_gs * qe));
  double nu_gr = q_gr / q_gs;
  double theta_nu = nu_gr / (1.0 + 1.0 / std::sqrt(nu_gr));
  double tmp = 1.0 + 1.0 / std::sqrt(4.0 * tau_gr + 3.0 * nu_gr);
  double J_rep = tmp * tmp * std::exp(-theta_nu / tau_gr);
  return pi * a_gr * a_gr * std::sqrt(8.0 * k_B * T_K / (pi * m_gs * m_p)) *
         J_rep;
}

// ---------------------------------------------------------------------------
// rate_TE — thermionic electron emission from grain
//   Richardson-Dushman with Richardson factor λ_R = 0.5, WF = 5 eV
//   j_TE = λ_R (4π m_e / h³) (kB T_gr)² exp(−(WF + q_gr e²/a) / kB T_gr)
//   k_TE = 4π a² j_TE   [s^{-1}]
// ---------------------------------------------------------------------------
inline double rate_TE(double T_gr_K, double a_gr, double q_gr) {
  constexpr double k_B = phys::k_B;
  constexpr double pi = phys::pi;
  constexpr double m_e = phys::m_e;
  constexpr double qe = phys::qe;
  constexpr double h_p = phys::h_P;
  constexpr double lmbd_R = 0.5;
  constexpr double WF = 5.0 * phys::eV_to_erg;  // 5 eV in erg

  double h3 = h_p * h_p * h_p;
  double kT = k_B * T_gr_K;
  double j_TE = lmbd_R * (4.0 * pi * m_e / h3) * kT * kT *
                std::exp(-(WF + q_gr * qe * qe / a_gr) / kT);
  return 4.0 * pi * a_gr * a_gr * j_TE;
}

// ---------------------------------------------------------------------------
// rate_gg — grain-grain charge transfer
//   Attractive geometry (q1·q2 < 0).
//   Symmetrizes J by averaging the two Coulomb focusing factors.
//   rho_gr = 3 g/cm³ (silicate bulk density)
// ---------------------------------------------------------------------------
inline double rate_gg(double T_K, double a1, double a2, double q1, double q2) {
  constexpr double k_B = phys::k_B;
  constexpr double pi = phys::pi;
  constexpr double qe = phys::qe;
  constexpr double rho_gr = 3.0;  // [g/cm³]

  double a_sum = a1 + a2;
  double m1 = (4.0 * pi / 3.0) * a1 * a1 * a1 * rho_gr;
  double m2 = (4.0 * pi / 3.0) * a2 * a2 * a2 * rho_gr;
  double mu = m1 * m2 / (m1 + m2);  // reduced mass

  // J from q1 perspective: grain 2 attracts grain 1
  double tau1 = a_sum * k_B * T_K / (q1 * qe * (q1 * qe));
  double nu1 = q2 / q1;
  double J1 = (1.0 - nu1 / tau1) * (1.0 + std::sqrt(2.0 / (tau1 - 2.0 * nu1)));

  // J from q2 perspective
  double tau2 = a_sum * k_B * T_K / (q2 * qe * (q2 * qe));
  double nu2 = q1 / q2;
  double J2 = (1.0 - nu2 / tau2) * (1.0 + std::sqrt(2.0 / (tau2 - 2.0 * nu2)));

  return pi * a_sum * a_sum * std::sqrt(8.0 * k_B * T_K / (pi * mu)) *
         (J1 + J2) / 2.0;
}

// ---------------------------------------------------------------------------
// mrn_quadrature — midpoint quadrature of the MRN size distribution over its two
//   power-law segments (seg 0: a_min..a_mid, seg 1: a_mid..a_max), in ln a:
//     a[s][i]  = exp(x_s + (i + 1/2) dl[s])      (bin-centre radius)
//     pa[s][i] = a[s][i]^(1 - ind_s),  pm[s] = a_mid^ind_s
//     wt[s][i] = pa[s][i] * dl[s] * pm[s]        (n(a) da weight of the bin)
//   and sum_ana = ∫ n(a)/n_ref da.  State-independent: built once.
// ---------------------------------------------------------------------------
struct MrnQuadrature {
  static constexpr int N1 = arche::metal_grain::mrn_nint1;
  static constexpr int N2 = arche::metal_grain::mrn_nint2;
  static constexpr int Nmax = N1 > N2 ? N1 : N2;
  double a[2][Nmax] = {};
  double pa[2][Nmax] = {};
  double wt[2][Nmax] = {};
  double dl[2] = {};
  double pm[2] = {};
  double sum_ana = 0.0;
};

inline MrnQuadrature make_mrn_quadrature() {
  constexpr double a_min = arche::metal_grain::mrn_a_min;
  constexpr double a_mid = arche::metal_grain::mrn_a_mid;
  constexpr double a_max = arche::metal_grain::mrn_a_max;
  constexpr double ind1 = arche::metal_grain::mrn_ind1;
  constexpr double ind2 = arche::metal_grain::mrn_ind2;
  constexpr int Nint1 = MrnQuadrature::N1;
  constexpr int Nint2 = MrnQuadrature::N2;
  MrnQuadrature q;
  const double x_min = std::log(a_min);
  const double x_mid = std::log(a_mid);
  const double x_max = std::log(a_max);
  const double d1_lna = (x_mid - x_min) / static_cast<double>(Nint1);
  const double d2_lna = (x_max - x_mid) / static_cast<double>(Nint2);
  q.dl[0] = d1_lna;
  q.dl[1] = d2_lna;
  q.pm[0] = std::pow(a_mid, ind1);
  q.pm[1] = std::pow(a_mid, ind2);
  q.sum_ana = std::pow(a_mid, ind1) *
                  (std::pow(a_mid, 1.0 - ind1) - std::pow(a_min, 1.0 - ind1)) /
                  (1.0 - ind1) +
              std::pow(a_mid, ind2) *
                  (std::pow(a_max, 1.0 - ind2) - std::pow(a_mid, 1.0 - ind2)) /
                  (1.0 - ind2);
  for (int i = 0; i < Nint1; ++i) {
    q.a[0][i] = std::exp((x_min + 0.5 * d1_lna) + i * d1_lna);
    q.pa[0][i] = std::pow(q.a[0][i], 1.0 - ind1);
    q.wt[0][i] =
        std::pow(q.a[0][i], 1.0 - ind1) * d1_lna * std::pow(a_mid, ind1);
  }
  for (int i = 0; i < Nint2; ++i) {
    q.a[1][i] = std::exp((x_mid + 0.5 * d2_lna) + i * d2_lna);
    q.pa[1][i] = std::pow(q.a[1][i], 1.0 - ind2);
    q.wt[1][i] =
        std::pow(q.a[1][i], 1.0 - ind2) * d2_lna * std::pow(a_mid, ind2);
  }
  return q;
}

inline const MrnQuadrature& mrn_quadrature() {
  static const MrnQuadrature q = make_mrn_quadrature();
  return q;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// react_grain_rates
//
//   Computes 150 MRN size-averaged charge-transfer rate coefficients and
//   stores them in k_charge[0..149].
//
//   Ordering (local index num = 1..150):
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
//     [142..145]  thermionic emission  (q_gr = −1,0,+1,+2)
//     [146]       gr⁺  + gr⁻   (rate_gg, q1=+1, q2=−1)
//     [147]       gr⁺  + gr²⁻  (rate_gg, q1=+1, q2=−2)
//     [148]       gr²⁺ + gr²⁻  (rate_gg, q1=+2, q2=−2)
//     [149]       gr²⁺ + gr⁻   (rate_gg, q1=+2, q2=−1)
//
//   Parameters
//     T_K    : gas temperature [K]
//     T_gr_K : grain temperature [K]  (used only for thermionic emission)
//     k_charge : output array, length 150
// ---------------------------------------------------------------------------
//   want   : want(i) asks for entry i (the caller's kept slots); the others are
//            not evaluated and stay 0.
template <class Want>
inline void react_grain_rates(double T_K, double T_gr_K,
                              std::array<double, 150>& k_charge, Want want) {
  // ── Species mass table [g], derived from the canonical catalog ────────
  // Indexed by metal_grain Sp index.  This kernel only touches the gas-phase
  // species 0..62 via the x1p/x2p/x1m lists below; the remaining grain/surface
  // entries are present but unused here.
  static constexpr auto mass = make_mass_array<metal_grain::Species>();

  // ── Charged species index lists (C++ 0-based into mass[]) ─────────────
  // 1-based x1p(1..30) = {4,5,6,9,11,14,15,23,24,25,26,27,28,29,
  //                        40,41,42,43,44,45,46,47,48,49,50,53,55,59,61,63}
  static constexpr int x1p[30] = {3,  4,  5,  8,  10, 13, 14, 22, 23, 24,
                                  25, 26, 27, 28, 39, 40, 41, 42, 43, 44,
                                  45, 46, 47, 48, 49, 52, 54, 58, 60, 62};
  // 1-based x2p(1..2) = {10, 56}
  static constexpr int x2p[2] = {9, 55};
  // 1-based x1m(1..4) = {3, 7, 16, 54}
  static constexpr int x1m[4] = {2, 6, 15, 53};

  // ── MRN grain size quadrature ─────────────────────────────────────────
  // The radii, weights and normalisation depend only on the MRN parameters, so
  // they are computed once (detail::mrn_quadrature); the averages multiply them in
  // a fixed order, so the coefficients equal a per-call evaluation bit for bit.
  constexpr int Nint1 = arche::metal_grain::mrn_nint1;
  constexpr int Nint2 = arche::metal_grain::mrn_nint2;
  const detail::MrnQuadrature& q = detail::mrn_quadrature();

  // ── 1-D MRN average: ∫ rate(a) · n(a) da / ∫ n(a) da ─────────────────
  // Uses midpoint rule over the two power-law segments.
  auto mrn_avg1 = [&](auto rate_fn) -> double {
    double val = 0.0;
    for (int itn = 0; itn < Nint1; ++itn)
      val += rate_fn(q.a[0][itn]) * q.pa[0][itn] * q.dl[0] * q.pm[0];
    for (int itn = 0; itn < Nint2; ++itn)
      val += rate_fn(q.a[1][itn]) * q.pa[1][itn] * q.dl[1] * q.pm[1];
    return val / q.sum_ana;
  };

  // ── Fill k_charge[0..149] ─────────────────────────────────────────────
  // An ion-grain rate is pi a^2 v_th J(tau, nu) with v_th = sqrt(8 kB T/(pi m m_p)) and
  // J from rate_gsneu / rate_gsatt / rate_gsrep: the ion's mass enters through v_th
  // alone, so each charge combination's MRN average of pi a^2 J is taken once and
  // multiplied by each ion's v_th.
  k_charge.fill(0.0);
  int num = 0;
  constexpr double k_B = phys::k_B, pi = phys::pi, m_p = phys::m_p,
                   qe = phys::qe;
  enum class Focus { neutral, attractive, repulsive };
  auto avg_area_J = [&](Focus focus, double q_gs, double q_gr) {
    return mrn_avg1([&](double a) {
      const double tau = a * k_B * T_K / (q_gs * qe * (q_gs * qe));
      const double nu = q_gr / q_gs;
      double J;
      if (focus == Focus::neutral) {
        J = 1.0 + std::sqrt(pi / (2.0 * tau));
      } else if (focus == Focus::attractive) {
        J = (1.0 - nu / tau) * (1.0 + std::sqrt(2.0 / (tau - 2.0 * nu)));
      } else {
        const double theta_nu = nu / (1.0 + 1.0 / std::sqrt(nu));
        const double t = 1.0 + 1.0 / std::sqrt(4.0 * tau + 3.0 * nu);
        J = t * t * std::exp(-theta_nu / tau);
      }
      return pi * a * a * J;
    });
  };
  // n ions (species x[0..n-1]) of charge q_gs on grains of charge q_gr.
  auto ion_group = [&](const int* x, int n, Focus focus, double q_gs,
                       double q_gr) {
    bool any = false;
    for (int i = 0; i < n; ++i) any = any || want(num + i);
    if (any) {
      const double area_J = avg_area_J(focus, q_gs, q_gr);
      for (int i = 0; i < n; ++i) {
        if (want(num + i))
          k_charge[num + i] =
              area_J * std::sqrt(8.0 * k_B * T_K / (pi * mass[x[i]] * m_p));
      }
    }
    num += n;
  };

  // ── (1) Neutral grain: all ions/anions with gr⁰ ─────────────────────
  ion_group(x1p, 30, Focus::neutral, 1.0, 0.0);  // X¹⁺ + gr⁰
  ion_group(x2p, 2, Focus::neutral, 2.0, 0.0);   // X²⁺ + gr⁰
  ion_group(x1m, 4, Focus::neutral, -1.0, 0.0);  // X⁻  + gr⁰
  // ── (2) Positive grain: gr⁺ ─────────────────────────────────────────
  ion_group(x1p, 30, Focus::repulsive, 1.0, 1.0);   // X¹⁺ + gr⁺
  ion_group(x1m, 4, Focus::attractive, -1.0, 1.0);  // X⁻  + gr⁺
  ion_group(x1m, 4, Focus::attractive, -1.0, 2.0);  // X⁻  + gr²⁺
  // ── (3) Negative grain: gr⁻ ─────────────────────────────────────────
  ion_group(x1p, 30, Focus::attractive, 1.0, -1.0);  // X¹⁺ + gr⁻
  ion_group(x2p, 2, Focus::attractive, 2.0, -1.0);   // X²⁺ + gr⁻
  ion_group(x1m, 4, Focus::repulsive, -1.0, -1.0);   // X⁻  + gr⁻
  // ── (4) Doubly negative grain: gr²⁻ ────────────────────────────────
  ion_group(x1p, 30, Focus::attractive, 1.0, -2.0);  // X¹⁺ + gr²⁻
  ion_group(x2p, 2, Focus::attractive, 2.0, -2.0);   // X²⁺ + gr²⁻

  // ── (5) Thermionic emission: gr^{q} → gr^{q+1} + e⁻  (4 entries) ─────
  // q_gr cycles: −1, 0, +1, +2
  for (int iq = 0; iq < 4; ++iq) {
    const double q_gr = -1.0 + iq;
    if (want(num))
      k_charge[num] =
          mrn_avg1([&](double a) { return detail::rate_TE(T_gr_K, a, q_gr); });
    ++num;
  }

  // ── (6) Grain-grain charge transfer  (4 entries) ──────────────────────
  // The 2-D MRN average of rate_gg.  rate_gg is symmetric in (a1, a2), so each
  // unordered pair of quadrature points is summed once (twice its weight off the
  // diagonal), and its size-only factors pi (a1+a2)^2 sqrt(8 kB/(pi mu)) w1 w2 are
  // computed once; per temperature only the two focusing factors remain.
  struct GrainPair {
    double a_sum, factor;
  };
  static const std::vector<GrainPair> pairs = [&q] {
    constexpr double rho_gr = 3.0;               // [g/cm³], as in rate_gg
    std::vector<std::pair<double, double>> pts;  // (radius, weight)
    for (int itn = 0; itn < Nint1; ++itn)
      pts.push_back({q.a[0][itn], q.wt[0][itn]});
    for (int itn = 0; itn < Nint2; ++itn)
      pts.push_back({q.a[1][itn], q.wt[1][itn]});
    std::vector<GrainPair> v;
    for (std::size_t i = 0; i < pts.size(); ++i) {
      for (std::size_t j = i; j < pts.size(); ++j) {
        const double a1 = pts[i].first, a2 = pts[j].first;
        const double m1 = (4.0 * pi / 3.0) * a1 * a1 * a1 * rho_gr;
        const double m2 = (4.0 * pi / 3.0) * a2 * a2 * a2 * rho_gr;
        const double mu = m1 * m2 / (m1 + m2);
        const double w = pts[i].second * pts[j].second * (i == j ? 1.0 : 2.0);
        const double a_sum = a1 + a2;
        v.push_back(
            {a_sum, pi * a_sum * a_sum * std::sqrt(8.0 * k_B / (pi * mu)) * w});
      }
    }
    return v;
  }();
  auto grain_grain = [&](double q1, double q2) {
    double val = 0.0;
    for (const GrainPair& pr : pairs) {
      const double tau1 = pr.a_sum * k_B * T_K / (q1 * qe * (q1 * qe));
      const double nu1 = q2 / q1;
      const double J1 =
          (1.0 - nu1 / tau1) * (1.0 + std::sqrt(2.0 / (tau1 - 2.0 * nu1)));
      const double tau2 = pr.a_sum * k_B * T_K / (q2 * qe * (q2 * qe));
      const double nu2 = q1 / q2;
      const double J2 =
          (1.0 - nu2 / tau2) * (1.0 + std::sqrt(2.0 / (tau2 - 2.0 * nu2)));
      val += pr.factor * (J1 + J2);
    }
    return std::sqrt(T_K) * val / 2.0 / (q.sum_ana * q.sum_ana);
  };
  // gr⁺  + gr⁻   → gr⁰  + gr⁰
  if (want(num)) k_charge[num] = grain_grain(1.0, -1.0);
  ++num;
  // gr⁺  + gr²⁻  → gr⁰  + gr⁻
  if (want(num)) k_charge[num] = grain_grain(1.0, -2.0);
  ++num;
  // gr²⁺ + gr²⁻  → gr⁺  + gr⁻
  if (want(num)) k_charge[num] = grain_grain(2.0, -2.0);
  ++num;
  // gr²⁺ + gr⁻   → gr⁺  + gr⁰   (last entry: num intentionally not
  // incremented)
  if (want(num)) k_charge[num] = grain_grain(2.0, -1.0);
  // num == 149 here; k_charge[0..149] fully populated (150 entries total)
}

inline void react_grain_rates(double T_K, double T_gr_K,
                              std::array<double, 150>& k_charge) {
  react_grain_rates(T_K, T_gr_K, k_charge, [](int) { return true; });
}

}  // namespace arche
