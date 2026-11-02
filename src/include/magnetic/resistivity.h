// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// resistivity.h — Ohmic and ambipolar magnetic diffusivities from a cell's
// abundances
//
// eta_O and eta_A [cm^2 s^-1] of a weakly ionised gas whose charged particles
// collide with the neutrals H, H2 and He (Wardle & Ng 1999; Nakauchi, Omukai &
// Susa 2021, eqs. 37-39):
//   eta_O = c^2 / (4 pi sigma_O),
//   eta_A = c^2 / (4 pi) * [sigma_P / (sigma_P^2 + sigma_H^2) - 1 / sigma_O],
//   sigma_O = sum_nu g_nu,             g_nu = n_nu q_nu^2 tau_nu / m_nu,
//   sigma_P = sum_nu p_nu,             p_nu = g_nu / (1 + beta_nu^2),
//   sigma_H = (c/B) sum_nu q_nu n_nu / (1 + beta_nu^2),
//   beta_nu = q_nu B tau_nu / (m_nu c)   (signed Hall parameter),
//   1/tau_nu = sum_n rho_n <sigma v>_nu,n / (m_nu + m_n),   n in {H, H2, He}.
//
// With charge neutrality, sum_nu h_nu = 0 for h_nu = (c/B) q_nu n_nu, so for
// any constant w
//   sigma_H = sum_nu h_nu (r_nu - w),   r_nu = 1 / (1 + beta_nu^2).
// w = 0 is the expression above and w = 1 gives -sum_nu p_nu beta_nu.  The
// terms h_nu r_nu cancel when every |beta| << 1 and the terms h_nu (r_nu - 1)
// cancel when every |beta| >> 1.  w is the median of r_nu weighted by |h_nu|,
// which minimises sum_nu |h_nu| |r_nu - w|, the bound of the rounding error of
// the sum.  The ambipolar numerator sigma_O sigma_P - sigma_P^2 - sigma_H^2
// equals the Lagrange-identity sum (1/2) sum_ij p_i p_j (beta_i - beta_j)^2,
// which is P^2 Var_p(beta) with P = sigma_P, Var_p the p-weighted variance of
// beta and <beta>_p = -sigma_H / P its mean.  Hence
//   eta_A = eta_O * Var_p(beta) / (1 + <beta>_p^2),
// with Var_p summed about <beta>_p, every term non-negative.  The identities
// hold for a neutral composition; a charge residual of the input enters
// sigma_H as -w (c/B) sum_nu q_nu n_nu.
//
// Charged particles: every species of the model with a non-zero charge
// (core/species_composition.h), with its catalog mass (core/species_catalog.h,
// times m_p).  The grain charge states (composition::grain_count == 1) are
// spread over the MRN size distribution of the metal_grain network
// (core/species_index.h: n(a) ~ a^-3.5 on [a_min, a_mid], a^-5.5 on
// [a_mid, a_max]; 15 + 5 bins equally spaced in ln a), each bin represented by
// its centre radius a_k with mass (4/3) pi a_k^3 rho_gr, rho_gr = 3 g cm^-3.
// The fraction of grains in a bin is the exact integral of n(a) over the bin
// divided by the integral over [a_min, a_max] (the fractions sum to 1); the
// fraction is the same for every charge state.
//
// Momentum-transfer rate coefficients <sigma v>:
//   e-H, e-H2, e-He, H+-H, H+-H2, H+-He, H3+-H2, HCO+-H2, C+-H
//       tabulated (magnetic/momentum_transfer_rates.h);
//   other ion-neutral pairs
//       Langevin, 2.41 pi |Z| e sqrt(alpha_n / mu) (Osterbrock 1961), with
//       alpha_H 0.667, alpha_H2 0.804, alpha_He 0.207 A^3 and mu the reduced
//       mass; the polarisation potential of a charge Ze scales as Z^2, so the
//       rate scales as |Z|;
//   grain-neutral
//       (4/3) pi a^2 delta sqrt(8 k T / (pi m_n)), delta = 1.3 (Pinto & Galli
//       2008, A&A 484, 17, eq. 25).
//
// Abundances below zero contribute nothing, and neutrals absent from a model
// collide with nothing.
// ---------------------------------------------------------------------------
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "core/species_catalog.h"
#include "core/species_composition.h"
#include "core/species_index.h"
#include "core/state.h"
#include "magnetic/momentum_transfer_rates.h"

namespace arche {
namespace magnetic {

// Ohmic and ambipolar diffusivities [cm^2 s^-1].
struct Diffusivities {
  double eta_O = std::numeric_limits<double>::quiet_NaN();
  double eta_A = std::numeric_limits<double>::quiet_NaN();
};

// Rate coefficients of the tabulated pairs.  A policy with the same static
// member can replace the tables (e.g. a published fit for a comparison).
struct TabulatedRates {
  static double rate(MtPair pair, double T_K) noexcept {
    return mt_rate(pair, T_K);
  }
};

inline constexpr double kSpeedOfLight = 2.99792458e10;  // [cm s^-1] (exact)

namespace detail {

inline constexpr int kNNeutral = 3;
inline constexpr std::array<SpId, kNNeutral> kNeutralIds = {
    {SpId::H, SpId::H2, SpId::He}};
// Static dipole polarisabilities [cm^3].
inline constexpr std::array<double, kNNeutral> kPolarizability = {
    {0.667e-24, 0.804e-24, 0.207e-24}};

inline constexpr double kGrainMaterialDensity = 3.0;  // [g cm^-3]
inline constexpr double kGrainNeutralDelta = 1.3;
inline constexpr int kGrainBins =
    metal_grain::mrn_nint1 + metal_grain::mrn_nint2;
inline constexpr int kGrainChargeStates = 5;  // Gr, Gr+, Gr2+, Gr-, Gr2-

// Tabulated pair of a charged species and a neutral, or -1.
inline constexpr int mt_pair_index(SpId charged, SpId neutral) {
  if (charged == SpId::e) {
    if (neutral == SpId::H) return static_cast<int>(MtPair::e_H);
    if (neutral == SpId::H2) return static_cast<int>(MtPair::e_H2);
    if (neutral == SpId::He) return static_cast<int>(MtPair::e_He);
  } else if (charged == SpId::Hp) {
    if (neutral == SpId::H) return static_cast<int>(MtPair::Hp_H);
    if (neutral == SpId::H2) return static_cast<int>(MtPair::Hp_H2);
    if (neutral == SpId::He) return static_cast<int>(MtPair::Hp_He);
  } else if (charged == SpId::H3p) {
    if (neutral == SpId::H2) return static_cast<int>(MtPair::H3p_H2);
  } else if (charged == SpId::HCOp) {
    if (neutral == SpId::H2) return static_cast<int>(MtPair::HCOp_H2);
  } else if (charged == SpId::Cp) {
    if (neutral == SpId::H) return static_cast<int>(MtPair::Cp_H);
  }
  return -1;
}

// T-independent data of one model, built once.
template <class Model>
struct Carriers {
  using Species = typename Model::Species;

  // Neutral colliders: local index (-1 if absent) and mass [g].
  std::array<int, kNNeutral> n_idx{};
  std::array<double, kNNeutral> n_mass{};

  // Gas-phase charged species.
  struct Gas {
    int idx = 0;
    double q = 0.0;                            // charge [statC], signed
    double m = 0.0;                            // mass [g]
    std::array<int, kNNeutral> pair{};         // MtPair index or -1
    std::array<double, kNNeutral> langevin{};  // [cm^3 s^-1]
    std::array<double, kNNeutral> inv_msum{};  // 1 / (m + m_n) [g^-1]
  };
  std::array<Gas, Model::N_sp> gas{};
  int n_gas = 0;

  // Grain charge states and size bins.
  struct GrainState {
    int idx = 0;
    double q = 0.0;  // [statC], signed
  };
  std::array<GrainState, kGrainChargeStates> grain{};
  int n_grain = 0;
  struct Bin {
    double m = 0.0;     // [g]
    double frac = 0.0;  // fraction of the grains in this bin
    double geo = 0.0;   // (4/3) pi a^2 delta [cm^2]
    std::array<double, kNNeutral> inv_msum{};
  };
  std::array<Bin, kGrainBins> bin{};

  Carriers() {
    for (int k = 0; k < kNNeutral; ++k) {
      n_idx[k] = Species::local(kNeutralIds[k]);
      n_mass[k] = species_mass(kNeutralIds[k]) * phys::m_p;
    }
    for (int i = 0; i < Model::N_sp; ++i) {
      const SpId id = Species::canonical(i);
      const int z = composition::composition_of(id).q;
      if (z == 0) continue;
      const double q = z * phys::qe;
      if (composition::grain_count(id) == 1) {
        grain[n_grain++] = GrainState{i, q};
        continue;
      }
      Gas& g = gas[n_gas++];
      g.idx = i;
      g.q = q;
      g.m = species_mass(id) * phys::m_p;
      for (int k = 0; k < kNNeutral; ++k) {
        const double mu = g.m * n_mass[k] / (g.m + n_mass[k]);
        g.pair[k] = mt_pair_index(id, kNeutralIds[k]);
        g.langevin[k] = 2.41 * phys::pi * std::abs(z) * phys::qe *
                        std::sqrt(kPolarizability[k] / mu);
        g.inv_msum[k] = 1.0 / (g.m + n_mass[k]);
      }
    }
    if (n_grain > 0) build_bins();
  }

  void build_bins() {
    constexpr double a_min = metal_grain::mrn_a_min;
    constexpr double a_mid = metal_grain::mrn_a_mid;
    constexpr double a_max = metal_grain::mrn_a_max;
    constexpr double ind[2] = {metal_grain::mrn_ind1, metal_grain::mrn_ind2};
    constexpr int nint[2] = {metal_grain::mrn_nint1, metal_grain::mrn_nint2};
    const double x_lo[2] = {std::log(a_min), std::log(a_mid)};
    const double x_hi[2] = {std::log(a_mid), std::log(a_max)};
    // Integral of n(a) = (a_mid / a)^ind over [lo, hi].
    auto count = [&](int s, double lo, double hi) {
      return std::pow(a_mid, ind[s]) *
             (std::pow(hi, 1.0 - ind[s]) - std::pow(lo, 1.0 - ind[s])) /
             (1.0 - ind[s]);
    };
    const double total = count(0, a_min, a_mid) + count(1, a_mid, a_max);
    int b = 0;
    for (int s = 0; s < 2; ++s) {
      const double dx = (x_hi[s] - x_lo[s]) / nint[s];
      for (int i = 0; i < nint[s]; ++i, ++b) {
        const double lo = std::exp(x_lo[s] + i * dx);
        const double hi = (i + 1 == nint[s]) ? std::exp(x_hi[s])
                                             : std::exp(x_lo[s] + (i + 1) * dx);
        const double a = std::exp(x_lo[s] + (i + 0.5) * dx);
        Bin& bn = bin[b];
        bn.m = (4.0 * phys::pi / 3.0) * a * a * a * kGrainMaterialDensity;
        bn.frac = count(s, lo, hi) / total;
        bn.geo = (4.0 * phys::pi / 3.0) * a * a * kGrainNeutralDelta;
        for (int k = 0; k < kNNeutral; ++k)
          bn.inv_msum[k] = 1.0 / (bn.m + n_mass[k]);
      }
    }
  }
};

template <class Model>
const Carriers<Model>& carriers() {
  static const Carriers<Model> c;
  return c;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// diffusivities — eta_O and eta_A of one cell.
//   y    abundances per H nucleus, in Model's local species order
//   nH   H-nucleus number density [cm^-3]
//   T_K  gas temperature [K]
//   B_G  magnetic field strength [G]
// Quiet NaN for both unless nH and T_K are positive and finite and B_G is
// non-negative and finite, and also when a charged particle is present but
// none of H, H2 and He is (no collision time).  With no charged particle
// (sigma_O = 0) eta_O = inf, and eta_A = inf for B_G > 0.  eta_A = 0 for
// B_G = 0.
// ---------------------------------------------------------------------------
template <class Model, class PairRates = TabulatedRates>
Diffusivities diffusivities(const double* y, double nH, double T_K,
                            double B_G) noexcept {
  Diffusivities out;
  if (!(nH > 0.0) || !std::isfinite(nH) || !(T_K > 0.0) ||
      !std::isfinite(T_K) || !(B_G >= 0.0) || !std::isfinite(B_G))
    return out;
  using detail::kNNeutral;
  const detail::Carriers<Model>& cr = detail::carriers<Model>();

  std::array<double, kNNeutral> rho_n{};
  std::array<double, kNNeutral> vth_n{};
  for (int k = 0; k < kNNeutral; ++k) {
    if (cr.n_idx[k] < 0) continue;
    rho_n[k] = std::max(y[cr.n_idx[k]], 0.0) * nH * cr.n_mass[k];
    vth_n[k] = std::sqrt(8.0 * phys::k_B * T_K / (phys::pi * cr.n_mass[k]));
  }
  std::array<double, static_cast<int>(MtPair::Count)> tab{};
  for (int p = 0; p < static_cast<int>(MtPair::Count); ++p)
    tab[p] = PairRates::rate(static_cast<MtPair>(p), T_K);

  // g = n q^2 tau / m, beta = q B tau / (m c) and h = (c/B) q n of every
  // carrier.
  constexpr int kMax =
      Model::N_sp + detail::kGrainChargeStates * detail::kGrainBins;
  std::array<double, kMax> g{};
  std::array<double, kMax> beta{};
  std::array<double, kMax> h{};
  int n = 0;
  bool collisionless = false;
  auto add = [&](double num, double q, double m, double inv_tau) {
    if (!(num > 0.0)) return;
    if (!(inv_tau > 0.0)) {
      collisionless = true;
      return;
    }
    const double tau = 1.0 / inv_tau;
    g[n] = num * q * q * tau / m;
    beta[n] = q * B_G * tau / (m * kSpeedOfLight);
    h[n] = B_G > 0.0 ? kSpeedOfLight * q * num / B_G : 0.0;
    ++n;
  };
  for (int c = 0; c < cr.n_gas; ++c) {
    const auto& s = cr.gas[c];
    double inv_tau = 0.0;
    for (int k = 0; k < kNNeutral; ++k) {
      if (!(rho_n[k] > 0.0)) continue;
      const double r = s.pair[k] >= 0 ? tab[s.pair[k]] : s.langevin[k];
      inv_tau += rho_n[k] * r * s.inv_msum[k];
    }
    add(y[s.idx] * nH, s.q, s.m, inv_tau);
  }
  for (int b = 0; b < (cr.n_grain > 0 ? detail::kGrainBins : 0); ++b) {
    const auto& bn = cr.bin[b];
    double inv_tau = 0.0;
    for (int k = 0; k < kNNeutral; ++k)
      if (rho_n[k] > 0.0)
        inv_tau += rho_n[k] * bn.geo * vth_n[k] * bn.inv_msum[k];
    for (int c = 0; c < cr.n_grain; ++c)
      add(y[cr.grain[c].idx] * nH * bn.frac, cr.grain[c].q, bn.m, inv_tau);
  }
  if (collisionless) return out;

  constexpr double kC2Over4Pi =
      kSpeedOfLight * kSpeedOfLight / (4.0 * phys::pi);
  constexpr double kInf = std::numeric_limits<double>::infinity();
  double sigma_O = 0.0;
  for (int i = 0; i < n; ++i) sigma_O += g[i];
  if (!(sigma_O > 0.0)) {
    out.eta_O = kInf;
    out.eta_A = B_G > 0.0 ? kInf : 0.0;
    return out;
  }
  out.eta_O = kC2Over4Pi / sigma_O;
  if (!(B_G > 0.0)) {
    out.eta_A = 0.0;
    return out;
  }

  // r = 1 / (1 + beta^2) and 1 - r = beta^2 / (1 + beta^2), each to full
  // precision.
  std::array<double, kMax> p{};
  std::array<double, kMax> r{};
  std::array<double, kMax> s{};
  std::array<int, kMax> order{};
  double P = 0.0;
  double h_abs = 0.0;
  for (int i = 0; i < n; ++i) {
    const double b2 = beta[i] * beta[i];
    r[i] = 1.0 / (1.0 + b2);
    s[i] = b2 / (1.0 + b2);
    p[i] = g[i] * r[i];
    P += p[i];
    h_abs += std::abs(h[i]);
    order[i] = i;
  }
  if (!(P > 0.0)) {
    out.eta_A = kInf;
    return out;
  }
  // |h|-weighted median of r (carriers ordered by |beta|).
  std::sort(order.begin(), order.begin() + n, [&](int a, int b) {
    return std::abs(beta[a]) < std::abs(beta[b]);
  });
  int km = order[n - 1];
  double acc = 0.0;
  for (int j = 0; j < n; ++j) {
    acc += std::abs(h[order[j]]);
    if (acc >= 0.5 * h_abs) {
      km = order[j];
      break;
    }
  }
  // r_i - r_k, written as s_k - s_i when r_k > 1/2.
  double sigma_H = 0.0;
  if (r[km] <= 0.5) {
    for (int i = 0; i < n; ++i) sigma_H += h[i] * (r[i] - r[km]);
  } else {
    for (int i = 0; i < n; ++i) sigma_H += h[i] * (s[km] - s[i]);
  }
  const double mean = -sigma_H / P;
  double var = 0.0;
  for (int i = 0; i < n; ++i) {
    const double d = beta[i] - mean;
    var += p[i] * d * d;
  }
  var /= P;
  out.eta_A = out.eta_O * var / (1.0 + mean * mean);
  return out;
}

}  // namespace magnetic
}  // namespace arche
