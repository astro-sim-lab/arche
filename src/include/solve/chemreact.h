// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// chemreact.h — implicit chemistry ODE integrator
//
// Provides the implicit chemistry solver chemreact (Newton–Raphson), the
// generic Newton driver newton_solve, and the dense linear solve gaussj.  The
// high-density Saha-equilibrium initial guess is dispatched through
// Model::SahaPolicy, defined in the per-model equilibrium.h headers.
//
// Optional Eigen support: define ARCHE_USE_EIGEN before including, or
// add -DARCHE_USE_EIGEN to the compiler flags.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <utility>

#include "cooling/cooling.h"  // c_H2
#include "core/newton.h"  // gaussj_solve, newton_solve, NewtonOpts, LinSolve
#include "kinetics/partition_function.h"  // PfProvider, pf_eval, detail primitives
#include "kinetics/reaction_index.h"
#include "kinetics/topology.h"
#include "models/rate_kernel.h"  // compute_base_rates, compute_rates

namespace arche {

// The high-density Saha branch of chemreact() dispatches to Model::SahaPolicy
// (FullPrimSaha / MinimalSaha / MetalSaha), forward-declared in the per-model
// traits.h and defined in the per-model equilibrium.h.  The including
// translation unit (chemistry.h) pulls in those definitions so the template
// instantiations resolve.

// ─────────────────────────────────────────────────────────────────────────────
// Error-handling policy
//
// Library layer (solve/, kinetics/, models/, cooling/ headers):
//   - Invalid caller input (e.g. unset T_rad) → throw std::invalid_argument.
//   - A diverged NR step (non-finite / unphysically large y) is reported
//   through
//     ChemSolveResult::converged, never via std::exit.
//   - gaussj_solve returns false on a singular system (interface honesty) and
//     zeroes the update; a transient singular Jacobian is absorbed as a
//     recoverable step (like the catastrophic / in-loop guards) and is NOT
//     escalated to converged==false — doing so regressed CR0 metal collapses.
// Application layer (collapse_*.h, apps/*/main.cc):
//   - May call std::exit, but only after writing a diagnostic to stderr.
// ─────────────────────────────────────────────────────────────────────────────

// The model-independent linear/Newton solvers (gaussj_solve, NewtonOpts,
// CramerSolve2, GaussjLinSolve, newton_solve) live in core/newton.h so the
// per-model equilibrium solvers can share them without a solve↔models cycle.

// ─────────────────────────────────────────────────────────────────────────────
// compute_chemistry_cooling<Model> — chemistry heating/cooling rate [erg g^-1
// s^-1]
//   Handles both zero_metal (primordial) and metal_grain networks via
//   if constexpr on Model traits.
//
//   Plain `inline` (not force-inlined): benchmarking showed no wall-clock
//   benefit from __attribute__((always_inline)) on the rate/cooling kernels
//   (forced inlining was marginally slower, within noise), so inlining is left
//   to the compiler.
// ─────────────────────────────────────────────────────────────────────────────
template <class Model>
inline void compute_chemistry_cooling(
    const std::array<double, Model::N_sp>& y,
    const std::array<double, Model::N_sp>& dy,
    const std::array<double, 2 * Model::N_react>& k_rxn, double nH, double T_K,
    double dt, double& Lambda_chem) {
  using Sp = typename Model::Sp;  // model-specific species names: y[Sp::H] etc.
  constexpr int N_react = Model::N_react;
  constexpr double yHe = abundance_ref::yHe;
  constexpr double m_p = phys::m_p;

  double L_H2 = 0.0;
  double dyHp = dy[Sp::Hp];
  // He+/He++ ionization-energy terms are a composable contribution: models that
  // carry He+ provide dyHep (has_he_ion), models that carry He++ provide dyHepp
  // (has_he_pp); a model without them leaves the term zero and never names
  // Sp::Hep / Sp::Hepp.  The additive energy sum below is left intact so the
  // full models stay bit-for-bit (both flags true ⇒ same assignments, in
  // order).
  double dyHep = 0.0;
  double dyHepp = 0.0;
  if constexpr (Model::cooling::has_he_ion) {
    dyHep = dy[Sp::Hep];
  }
  if constexpr (Model::cooling::has_he_pp) {
    dyHepp = dy[Sp::Hepp];
  }

  if (nH < 1.0e13) {
    double denom_cr = 1.6 * y[Sp::H] * std::exp(-std::pow(400.0 / T_K, 2.0)) +
                      1.4 * y[Sp::H2] * std::exp(-12000.0 / (T_K + 1200.0));
    double n_cr =
        (denom_cr > 1.0e-30) ? 1.0e6 / std::sqrt(T_K) / denom_cr : 1.0e50;
    double crit = 1.0 / (1.0 + n_cr / nH);

    // k_rxn[] indices are addressed through the model's cooling rate-index map
    // so a compact model (different reaction numbering) reaches the same
    // reactions; for the full models these equal the shared rxn:: ids, so the
    // expressions stay bit-for-bit.  The H2–He dissociation channel is a
    // composable term (a model without it sets has_h2_he_cool = false); each
    // accumulation preserves the original left-to-right order.
    using cidx = typename Model::cooling::rate_idx;
    double rtHm = k_rxn[cidx::Hm_H_to_H2_e] * y[Sp::Hm] * y[Sp::H];
    double rtH2p = k_rxn[cidx::H2p_H_to_H2_Hp] * y[Sp::H2p] * y[Sp::H];
    double rt3b_sum =
        k_rxn[cidx::three_H] * y[Sp::H] * y[Sp::H] * y[Sp::H] +
        k_rxn[cidx::H2H2_dis + N_react] * y[Sp::H] * y[Sp::H] * y[Sp::H2];
    if constexpr (Model::cooling::has_h2_he_cool)
      rt3b_sum +=
          k_rxn[cidx::H2_He_dis + N_react] * y[Sp::H] * y[Sp::H] * y[Sp::He];
    double rt3b = rt3b_sum * nH;
    double rtdis = k_rxn[cidx::three_H + N_react] * y[Sp::H2] * y[Sp::H] +
                   k_rxn[cidx::H2H2_dis] * y[Sp::H2] * y[Sp::H2];
    if constexpr (Model::cooling::has_h2_he_cool)
      rtdis += k_rxn[cidx::H2_He_dis] * y[Sp::H2] * y[Sp::He];

    if constexpr (Model::has_grain) {
      // Grain H2 formation cooling.  The surface channels (physisorbed /
      // chemisorbed H) need H_p / H_c, which the compact network drops, so they
      // are composable (gh2::has_surface); the full expression is kept verbatim
      // when present so the full model stays bit-for-bit.
      using gh2 = typename Model::cooling::grain_H2;
      double rtgrain;
      if constexpr (gh2::has_surface)
        rtgrain = k_rxn[gh2::two_H] * y[Sp::H] +
                  k_rxn[gh2::surface_grgr] * y[Sp::H_p] * y[Sp::H_p] +
                  k_rxn[gh2::surface_grdust] * y[Sp::H_p] * y[Sp::H_c] +
                  k_rxn[gh2::H_dust_3body] * y[Sp::H] * y[Sp::H_c] * nH +
                  k_rxn[gh2::H_gr_3body] * y[Sp::H] * y[Sp::H_p] * nH;
      else
        rtgrain = k_rxn[gh2::two_H] * y[Sp::H];

      L_H2 = -(rtgrain * (0.2 + 4.2 * crit) + rtHm * 3.53 * crit +
               rtH2p * 1.83 * crit + (rt3b * crit - rtdis) * 4.48) *
             nH * 1.60219e-12;
    } else {
      L_H2 = -(rtHm * 3.53 * crit + rtH2p * 1.83 * crit +
               (rt3b * crit - rtdis) * 4.48) *
             nH * 1.60219e-12;
    }

    // H+ recombination minus the CR-producer losses; the CR term is composable
    // (a CR-free model omits it).  Each branch is the model's full expression
    // verbatim, so the additive energy sum stays bit-for-bit.
    if constexpr (Model::cooling::has_cr_loss) {
      using Hp = typename Model::cooling::Hp_producers;
      // The second H2 CR channel (CR-induced photon) is composable: a model
      // whose keep-set lacks it sets H2_CR_ch2 = -1.  Full models keep both
      // channels, so the sum (ch1 + ch2) stays byte-identical.
      double h2_cr_loss = k_rxn[Hp::H2_CR_ch1];
      if constexpr (Hp::H2_CR_ch2 >= 0) h2_cr_loss += k_rxn[Hp::H2_CR_ch2];
      dyHp =
          dy[Sp::Hp] + (k_rxn[cidx::Hp_rec_caseB] * y[Sp::Hp] * y[Sp::e] * nH -
                        (k_rxn[Hp::H_CR] + k_rxn[Hp::H_CRph]) * y[Sp::H] -
                        h2_cr_loss * y[Sp::H2]) *
                           dt;
    } else {
      dyHp = dy[Sp::Hp] +
             (k_rxn[cidx::Hp_rec_caseB] * y[Sp::Hp] * y[Sp::e] * nH) * dt;
    }
    // He+ recombination minus CR losses — only for models carrying He ions; the
    // CR-producer loss is composable (a He+-but-CR-free model omits it,
    // mirroring the dyHp split above).  Full models keep both terms,
    // byte-identical.
    if constexpr (Model::cooling::has_he_ion) {
      if constexpr (Model::cooling::has_cr_loss) {
        using Hep = typename Model::cooling::Hep_producers;
        dyHep = dy[Sp::Hep] +
                (k_rxn[cidx::Hep_rec] * y[Sp::Hep] * y[Sp::e] * nH -
                 (k_rxn[Hep::He_CR] + k_rxn[Hep::He_CRph]) * y[Sp::He]) *
                    dt;
      } else {
        dyHep = dy[Sp::Hep] +
                (k_rxn[cidx::Hep_rec] * y[Sp::Hep] * y[Sp::e] * nH) * dt;
      }
    }

  } else {
    double dyH2 = dy[Sp::H2];
    L_H2 = -7.18e-12 * dyH2 / dt;
    constexpr double L_H2_max = 1.0e20;
    L_H2 = std::max(-L_H2_max, std::min(L_H2, L_H2_max));
  }

  Lambda_chem =
      ((2.18e-11 * dyHp + 3.94e-11 * dyHep + 12.66e-11 * dyHepp) / dt + L_H2) /
      ((1.0 + 4.0 * yHe) * m_p);
  if (!std::isfinite(Lambda_chem)) Lambda_chem = 0.0;
}

// ─────────────────────────────────────────────────────────────────────────────
// eos_particle_sum<Model> — free-particle count for the equation of state.
//   Sums the neutral and ionic species the model carries (H [+H2] + e + H+ +
//   He), adding He+/He++ only for models that have them.  The summation order
//   matches the historical inline expressions so mu/gamma stay bit-for-bit for
//   the full models; a model without He ions simply omits the last two terms.
// ─────────────────────────────────────────────────────────────────────────────
template <class Model>
inline double eos_particle_sum(const std::array<double, Model::N_sp>& y,
                               bool include_H2) {
  using Sp = typename Model::Sp;
  double s = y[Sp::H];
  if (include_H2) s += y[Sp::H2];
  s += y[Sp::e];
  s += y[Sp::Hp];
  s += y[Sp::He];
  if constexpr (Model::cooling::has_he_ion) {
    s += y[Sp::Hep];
  }
  if constexpr (Model::cooling::has_he_pp) {
    s += y[Sp::Hepp];
  }
  return s;
}

// ─────────────────────────────────────────────────────────────────────────────
// ChemSolveResult — scalar results of chemreact / chemcool.
//   y (abundances) and var (rate cache) are updated in place via reference
//   parameters; these computed scalars are returned by value.
// ─────────────────────────────────────────────────────────────────────────────
struct ChemSolveResult {
  double mu = 0.0;
  double gamma = 0.0;
  double Lambda_chem = 0.0;  // chemistry cooling [erg g^-1 s^-1]
  bool converged = false;
};

// ─────────────────────────────────────────────────────────────────────────────
// chemreact<Model> — one timestep NR chemistry solver
//   Solves dy/dt = r_f(y) implicitly, optionally via Saha for nH > 1e18.
//   Updates y and var in place; returns {mu, gamma, Lambda_chem, converged}.
// ─────────────────────────────────────────────────────────────────────────────
template <class Model>
ChemSolveResult chemreact(double nH, double T_K,
                          std::array<double, Model::N_sp>& y, double dt,
                          std::array<double, 2 * Model::N_react>& var,
                          const ReactionTable<Model::N_sp, Model::N_react>& tbl,
                          const ChemParams& params) {
  constexpr int N_sp = Model::N_sp;
  constexpr int N_react = Model::N_react;
  using Sp = typename Model::Sp;  // model-specific species names: y[Sp::H] etc.
  constexpr double yHe = abundance_ref::yHe;
  constexpr double m_p = phys::m_p;
  constexpr double eps_y = numerics::eps_y;
  constexpr double nH_eq = numerics::nH_eq;

  // Scalar results returned by value; y and var are updated in place.
  double mu = 0.0, gamma = 0.0, Lambda_chem = 0.0;
  bool converged = false;

  std::array<double, N_sp> y_init = y;
  std::array<double, N_sp> dy{};
  dy.fill(0.0);
  std::array<double, N_sp> ddy{};
  std::array<double, N_sp> r_f{};
  std::array<double, N_sp * N_sp> dr_fdy{};
  std::array<double, N_sp * N_sp> A_mat{};
  std::array<double, 2 * N_react> k_rxn{};

  // Non-thermal ionization (CR + X-ray) maintains ye above pure Saha thermal
  // equilibrium at high density.  Skip the Saha branch when either CR or
  // X-ray ionization is active so the NR solver remains active.
  const bool use_nr =
      (nH <= nH_eq) || (params.zeta > 0.0) || (params.zeta_X > 0.0);
  if (use_nr) {
    // ─── Non-Equilibrium (low-density) NR solver ───────────────────────
    ChemParams p_loc = params;
    p_loc.H = y[Sp::H];
    p_loc.H2 = 2.0 * y[Sp::H2];
    p_loc.He = y[Sp::He];

    if constexpr (Model::has_uv_shield) {
      p_loc.J_H2 = y[Sp::H2_p];
      p_loc.J_H2O = y[Sp::H2O_p];
      p_loc.J_tot = y[Sp::H_c] + y[Sp::D_c];
    }

    // Compute reaction rate coefficients
    mu = (1.0 + 4.0 * yHe) / eos_particle_sum<Model>(y, /*include_H2=*/true);
    compute_base_rates<Model>(nH, T_K, mu, p_loc, tbl, k_rxn);

    // Catastrophic-detection guard
    const bool use_cat_detect = [&]() {
      if constexpr (Model::catastrophic_detect_always) {
        return true;
      } else {
        return (nH < 4.0e22);
      }
    }();
    double err_fnc_prev = 0.0;
    for (int itr = 0; itr < Model::nr_max_iter; ++itr) {
      r_f.fill(0.0);
      dr_fdy.fill(0.0);

      compute_rates<Model>(k_rxn, nH, y, tbl, r_f, dr_fdy, var);

      // Build A = I - dt * J  (J = dr_fdy, row-major)
      for (int isp = 0; isp < N_sp; ++isp)
        for (int jsp = 0; jsp < N_sp; ++jsp)
          A_mat[isp * N_sp + jsp] =
              (isp == jsp ? 1.0 : 0.0) - dt * dr_fdy[isp * N_sp + jsp];

      // rhs: ddy = r_f*dt - dy
      double err_fnc = 0.0;
      for (int isp = 0; isp < N_sp; ++isp) {
        ddy[isp] = r_f[isp] * dt - dy[isp];
        err_fnc += std::abs(ddy[isp]);
      }

      if (use_cat_detect && itr > 0 &&
          err_fnc > 1.0e4 * std::max(err_fnc_prev, 1.0e-12)) {
        y = y_init;
        dy.fill(0.0);
        break;
      }
      err_fnc_prev = err_fnc;

      // Solve A * ddy = ddy  (in-place).  A singular Jacobian is absorbed
      // by gaussj_solve (it zeroes the update and returns false); the step
      // is treated as a recoverable transient — consistent with the
      // catastrophic / in-loop guards — so the bool return is not escalated
      // to a convergence failure here.
      gaussj_solve<N_sp>(A_mat, ddy);

      // Update abundances
      double err_y = 0.0;
      for (int isp = 0; isp < N_sp; ++isp) {
        if (y[isp] + ddy[isp] < 0.0) {
          ddy[isp] = -0.1 * y[isp];
          err_y += 1.0;
        }
        dy[isp] += ddy[isp];
        y[isp] += ddy[isp];
        err_y += std::abs(ddy[isp]);
      }

      // In-loop divergence guard (metal only): non-finite or unphysically large
      if constexpr (Model::has_in_loop_divergence_guard) {
        bool nr_ok_guard = true;
        for (int isp = 0; isp < N_sp && nr_ok_guard; ++isp)
          if (!std::isfinite(y[isp]) || y[isp] > 1.0e50) nr_ok_guard = false;
        if (!nr_ok_guard) {
          y = y_init;
          dy.fill(0.0);
          break;
        }
      }

      if (err_y <= eps_y && err_fnc <= 1.0e-4) break;
    }
    // Convergence notification — reported regardless of build options so the
    // caller is never left unaware of a diverged step.
    converged = true;
    for (int isp = 0; isp < N_sp && converged; ++isp)
      if (!std::isfinite(y[isp]) || y[isp] > 1.0e50) converged = false;
#ifdef ARCHE_SUBCYCLE
    // Subcycle builds roll the diverged step back so the caller can retry at
    // a smaller dt; non-subcycle numerics are intentionally left untouched.
    if (!converged) y = y_init;
#endif
  } else {
    // ─── Saha equilibrium (high-density) ───────────────────────────────
    // Each model names its Saha-equilibrium solver via Model::SahaPolicy
    // (FullPrimSaha / MinimalSaha / MetalSaha), defined in the per-model
    // equilibrium.h.  equichem_dy_count bounds the species the Saha branch
    // writes (the full network for primordial, y[0..62] for metal_grain).
    Model::SahaPolicy::template solve<N_sp, N_react>(nH, T_K, params, y, tbl);
    for (int i = 0; i < Model::equichem_dy_count; ++i) dy[i] = y[i] - y_init[i];
    converged = true;
  }

  // Update thermodynamic variables
  mu = (1.0 + 4.0 * yHe) / eos_particle_sum<Model>(y, /*include_H2=*/true);

  gamma =
      1.0 + (1.0 + 4.0 * yHe) /
                (mu * (1.5 * eos_particle_sum<Model>(y, /*include_H2=*/false) +
                       c_H2(T_K) * y[Sp::H2]));

  // ── Chemistry cooling Lambda_chem ──────────────────────────────────────────
  compute_chemistry_cooling<Model>(y, dy, k_rxn, nH, T_K, dt, Lambda_chem);

  return ChemSolveResult{mu, gamma, Lambda_chem, converged};
}

// ─────────────────────────────────────────────────────────────────────────────
// chemcool<Model> — secant method wrapper for temperature self-consistent
// integration
//   Finds T* such that T_K - T* - (γ-1)*μ*mp*ΛΔt/kB = 0,
//   then updates y in place and returns {mu, gamma, Lambda_chem, converged} at
//   T*.
//
//   secant_skip_high_density controls the early-exit condition:
//     primordial:   skip when  nH  > 1e16  AND  T <= 1650 K
//     metal_grain:  skip when  nH <= 1e16  AND  T <= 1650 K
// ─────────────────────────────────────────────────────────────────────────────
template <class Model>
ChemSolveResult chemcool(double nH, double Tp,
                         std::array<double, Model::N_sp>& y, double dt,
                         std::array<double, 2 * Model::N_react>& var,
                         const ReactionTable<Model::N_sp, Model::N_react>& tbl,
                         const ChemParams& params) {
  constexpr int N_sp = Model::N_sp;
  constexpr int N_react = Model::N_react;
  constexpr double m_p = phys::m_p;
  constexpr double k_B = phys::k_B;
  constexpr int maxit = 100;

  // Residual using updated gamma/mu from chemreact output
  auto func = [&](double t, double g, double mu, double Lc) -> double {
    return Tp - t - (g - 1.0) * mu * m_p * Lc * dt / k_B;
  };

  // First evaluation at Tp1 = Tp
  double Tp1 = Tp;
  std::array<double, N_sp> ytmp = y;
  double mu1 = 0.0, gamma1 = 0.0, Lambdach1 = 0.0;
  // Optimistic; flipped to false if any inner chemreact fails to converge.
  bool converged = true;
  ChemSolveResult r = chemreact<Model>(nH, Tp1, ytmp, dt, var, tbl, params);
  mu1 = r.mu;
  gamma1 = r.gamma;
  Lambdach1 = r.Lambda_chem;
  if (!r.converged) {
    converged = false;
#ifdef ARCHE_SUBCYCLE
    return ChemSolveResult{mu1, gamma1, Lambdach1, false};  // subcycle: bail
#endif
  }
  double fL = func(Tp1, gamma1, mu1, Lambdach1);

  // Skip secant: condition depends on model
  {
    bool skip;
    if constexpr (Model::secant_skip_high_density) {
      skip = (nH > 1.0e16 && Tp <= 1.65e3);
    } else {
      skip = (nH <= 1.0e16 && Tp <= 1.65e3);
    }
    if (skip) {
      y = ytmp;  // converged already reflects the single chemreact above
      return ChemSolveResult{mu1, gamma1, Lambdach1, converged};
    }
  }

  // Second evaluation at Tp2 = 0.999*Tp1
  double Tp2 = 0.999 * Tp1;
  ytmp = y;
  double mu2 = 0.0, gamma2 = 0.0, Lambdach2 = 0.0;
  r = chemreact<Model>(nH, Tp2, ytmp, dt, var, tbl, params);
  mu2 = r.mu;
  gamma2 = r.gamma;
  Lambdach2 = r.Lambda_chem;
  if (!r.converged) {
    converged = false;
#ifdef ARCHE_SUBCYCLE
    return ChemSolveResult{mu2, gamma2, Lambdach2, false};
#endif
  }
  double f = func(Tp2, gamma2, mu2, Lambdach2);

  // Ensure |fL| >= |f|  (larger residual at TpL)
  double TpL, Tpsec;
  if (std::abs(fL) < std::abs(f)) {
    TpL = Tp2;
    Tpsec = Tp1;
    std::swap(fL, f);
  } else {
    TpL = Tp1;
    Tpsec = Tp2;
  }

  // Secant iterations
  for (int i = 0; i < maxit; ++i) {
    // Guard against f ≈ fL which causes division by zero → NaN
    double df = f - fL;
    if (std::abs(df) < 1.0e-30 * (std::abs(f) + std::abs(fL) + 1.0e-30)) break;
    double dTp = (TpL - Tpsec) * f / df;
    TpL = Tpsec;
    fL = f;
    Tpsec += dTp;
    // Physical temperature bounds: prevent secant divergence from
    // propagating NaN cooling rates into the main loop
    if (!std::isfinite(Tpsec) || Tpsec < 1.0) Tpsec = 1.0;
    if (Tpsec > 1.0e8) Tpsec = 1.0e8;

    ytmp = y;
    r = chemreact<Model>(nH, Tpsec, ytmp, dt, var, tbl, params);
    mu1 = r.mu;
    gamma1 = r.gamma;
    Lambdach1 = r.Lambda_chem;
    if (!r.converged) {
      converged = false;
#ifdef ARCHE_SUBCYCLE
      return ChemSolveResult{mu1, gamma1, Lambdach1, false};
#endif
    }
    f = func(Tpsec, gamma1, mu1, Lambdach1);

    if (std::abs(f / Tpsec) <= 1.0e-7) break;
  }

  y = ytmp;
  // converged is true here unless a non-subcycle build saw a diverged step.
  return ChemSolveResult{mu1, gamma1, Lambdach1, converged};
}

}  // namespace arche
