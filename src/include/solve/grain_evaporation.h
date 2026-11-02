// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// grain_evaporation.h — release of the metals locked in the dust as the grains
// evaporate (grain models), over one call of chem_integrate_adaptive.
//
// The dust amount follows the grain temperature: each constituent evaporates above
// its vaporization temperature, and the grains keep their size distribution, so the
// grain number in every charge state scales with the grain volume.  Above the
// silicate vaporization temperature the depleted metals return to the gas in
// proportion to the decrease of the grain volume [Nakauchi, Omukai & Susa 2021,
// MNRAS 502, 3394, sec. 2; Finocchi & Gail 1997].  The release is irreversible: it
// follows the smallest grain volume reached so far.
//
// State (species of the model):
//   y[Sp::X_s]  metal X locked in the dust, per H nucleus (X = C, O, Mg; the
//               full network also Na and K).  Counted as element X by the
//               conservation rows.
//   y[Sp::Gr_v] smallest grain volume reached so far, relative to the volume at
//               T -> 0 (vol_gr(rho, T_gr) / vol_gr(rho, 0)); 0 = not yet set.
//
// One call over dt, planned from its input state (the grain temperature T_gr at
// the input gas temperature, rho):
//   V      = vol_gr(rho, T_gr) / vol_gr(rho, 0)
//   f_all  = min(V / Gr_v, 1);  f_rel = f_all if T_gr >= 0.975 T_pyr(rho), else 1
//   vol_gr(rho, T_gr) < 1e-60 (complete evaporation): f_all = f_rel = 0
// and over the call
//   - each X_s releases the fraction 1 - f_rel to the neutral X at a constant
//     rate (X_s ends at f_rel X_s);
//   - the grain species (neutral and charged) decay at the first-order rate
//     -ln(f_all) / dt; their net charge goes to the electrons (full network,
//     while the metals are released: a net positive charge goes to K+, taken
//     from K);
//   - Gr_v becomes min(Gr_v, V) at a constant rate.
// The constant rates and the first-order sinks enter every backward-Euler step of
// the call (solve/integrator.h), so the released metals react while they are
// released.  A complete evaporation sets the remaining grain species to 0 at the
// end of the call (charge as above).  Gr_v = 0 with grains present is set to V
// at the start of the call, with no release in that call.
//
// The release is taken from the call's input state because the caller applies
// the compression heating between calls: within a call the temperature moves
// only through the chemistry and the radiation.
// ---------------------------------------------------------------------------
#include <algorithm>
#include <array>
#include <cmath>

#include "cooling/grain.h"       // detail::vaptemp, detail::vol_gr
#include "core/species_index.h"  // abundance_ref
#include "core/state.h"          // phys

namespace arche {

template <class Model>
struct GrainRelease {
  static constexpr int N = Model::N_sp;
  bool active = false;
  std::array<double, N> src{};  // constant dy/dt over the call [1/s]
  double k = 0.0;               // first-order sink of the grain species [1/s]
  bool charge_to_k = false;     // net positive grain charge -> K+ (from K)
  bool complete = false;        // complete evaporation in this call
};

namespace grain_evap {

// Relative volume below which the grains are gone (complete evaporation).
constexpr double kVolumeGone = 1.0e-60;
// Fraction of T_pyr above which the locked metals are released.
constexpr double kReleaseFraction = 0.975;
// Floor of a remaining fraction in the first-order rate -ln(f) / dt.
constexpr double kFractionFloor = 1.0e-300;

// Net charge of the grain species, in units of e.
template <class Model>
double grain_charge(const std::array<double, Model::N_sp>& y) {
  using Sp = typename Model::Sp;
  double q = y[Sp::Grp] - y[Sp::Grm] - 2.0 * y[Sp::Gr2m];
  if constexpr (Model::grain_evaporation::has_gr2p) q += 2.0 * y[Sp::Gr2p];
  return q;
}

template <class Model>
bool has_grains(const std::array<double, Model::N_sp>& y) {
  using Sp = typename Model::Sp;
  bool any = y[Sp::Gr] > 0.0 || y[Sp::Grp] > 0.0 || y[Sp::Grm] > 0.0 ||
             y[Sp::Gr2m] > 0.0;
  if constexpr (Model::grain_evaporation::has_gr2p)
    any = any || y[Sp::Gr2p] > 0.0;
  return any;
}

// Adds rate * (the grains' sink at y) to out: the grain species at -k y, their
// charge to the electrons or to K+ (from K).
template <class Model>
void add_sink(const GrainRelease<Model>& rel,
              const std::array<double, Model::N_sp>& y, double rate,
              std::array<double, Model::N_sp>& out) {
  using Sp = typename Model::Sp;
  const double k = rate * rel.k;
  out[Sp::Gr] -= k * y[Sp::Gr];
  out[Sp::Grp] -= k * y[Sp::Grp];
  out[Sp::Grm] -= k * y[Sp::Grm];
  out[Sp::Gr2m] -= k * y[Sp::Gr2m];
  if constexpr (Model::grain_evaporation::has_gr2p)
    out[Sp::Gr2p] -= k * y[Sp::Gr2p];
  const double dq = k * grain_charge<Model>(y);  // positive charge lost
  if constexpr (Model::grain_evaporation::has_na_k) {
    if (rel.charge_to_k) {
      out[Sp::Kp] += dq;
      out[Sp::K] -= dq;
      return;
    }
  }
  out[Sp::e] -= dq;
}

}  // namespace grain_evap

// The release over a call of width dt from the input state.  Sets Gr_v in y
// when it is not yet set.
template <class Model>
GrainRelease<Model> plan_grain_release(double nH, double T_gr, double dt,
                                       std::array<double, Model::N_sp>& y) {
  using Sp = typename Model::Sp;
  constexpr int N = Model::N_sp;
  GrainRelease<Model> rel;
  if (!grain_evap::has_grains<Model>(y)) return rel;
  const double rho = nH * ((1.0 + 4.0 * abundance_ref::yHe) * phys::m_p);
  const double vol = detail::vol_gr(rho, T_gr);
  const double V = vol / detail::vol_gr(rho, 0.0);
  if (!(y[Sp::Gr_v] > 0.0)) {
    y[Sp::Gr_v] = V;
    return rel;
  }
  double T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol;
  detail::vaptemp(rho, T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol);
  double f_all, f_rel;
  if (vol < grain_evap::kVolumeGone) {
    rel.complete = true;
    f_all = f_rel = 0.0;
  } else {
    f_all = std::min(V / y[Sp::Gr_v], 1.0);
    f_rel = (T_gr >= grain_evap::kReleaseFraction * T_pyr) ? f_all : 1.0;
  }
  auto rate_of = [dt](double f) {
    return f < 1.0 ? -std::log(std::max(f, grain_evap::kFractionFloor)) / dt
                   : 0.0;
  };
  rel.k = rate_of(f_all);
  auto release = [&](int locked, int gas) {
    const double r = y[locked] * (1.0 - f_rel) / dt;
    rel.src[gas] += r;
    rel.src[locked] -= r;
  };
  release(Sp::C_s, Sp::C);
  release(Sp::O_s, Sp::O);
  release(Sp::Mg_s, Sp::Mg);
  if constexpr (Model::grain_evaporation::has_na_k) {
    release(Sp::Na_s, Sp::Na);
    release(Sp::K_s, Sp::K);
    rel.charge_to_k = f_rel < 1.0 && grain_evap::grain_charge<Model>(y) > 0.0;
  }
  const double V_end = rel.complete ? 0.0 : std::min(y[Sp::Gr_v], V);
  rel.src[Sp::Gr_v] = (V_end - y[Sp::Gr_v]) / dt;
  bool any = rel.k > 0.0;
  for (int i = 0; i < N && !any; ++i) any = rel.src[i] != 0.0;
  rel.active = any;
  return rel;
}

namespace grain_evap {

// dr/dy of the sinks, column by column: the sink is linear in the grain species.
template <class Model>
void add_sink_jacobian(const GrainRelease<Model>& rel,
                       std::array<double, Model::N_sp * Model::N_sp>& J) {
  using Sp = typename Model::Sp;
  constexpr int N = Model::N_sp;
  auto column = [&](int j, double q) {
    J[j * N + j] -= rel.k;
    const double dq = rel.k * q;  // d(positive charge lost)/dy_j
    if constexpr (Model::grain_evaporation::has_na_k) {
      if (rel.charge_to_k) {
        J[Sp::Kp * N + j] += dq;
        J[Sp::K * N + j] -= dq;
        return;
      }
    }
    J[Sp::e * N + j] -= dq;
  };
  column(Sp::Gr, 0.0);
  column(Sp::Grp, 1.0);
  column(Sp::Grm, -1.0);
  column(Sp::Gr2m, -2.0);
  if constexpr (Model::grain_evaporation::has_gr2p) column(Sp::Gr2p, 2.0);
}

template <class Model>
void add_sinks_and_jacobian(const GrainRelease<Model>& rel,
                            const std::array<double, Model::N_sp>& y,
                            std::array<double, Model::N_sp>& r,
                            std::array<double, Model::N_sp * Model::N_sp>& J) {
  add_sink<Model>(rel, y, 1.0, r);
  add_sink_jacobian<Model>(rel, J);
}

template <class Model>
void release_start(const GrainRelease<Model>& rel, double h,
                   std::array<double, Model::N_sp>& y) {
  using Sp = typename Model::Sp;
  for (int i = 0; i < Model::N_sp; ++i) y[i] += h * rel.src[i];
  const double decay = 1.0 + h * rel.k;
  y[Sp::Gr] /= decay;
  y[Sp::Grp] /= decay;
  y[Sp::Grm] /= decay;
  y[Sp::Gr2m] /= decay;
  if constexpr (Model::grain_evaporation::has_gr2p) y[Sp::Gr2p] /= decay;
}

template <class Model>
void remove_grains(std::array<double, Model::N_sp>& y) {
  using Sp = typename Model::Sp;
  const double q = grain_charge<Model>(y);
  bool to_k = false;
  if constexpr (Model::grain_evaporation::has_na_k) {
    to_k = q > 0.0;
    if (to_k) {
      y[Sp::Kp] += q;
      y[Sp::K] -= q;
    }
  }
  if (!to_k) y[Sp::e] -= q;
  y[Sp::Gr] = y[Sp::Grp] = y[Sp::Grm] = y[Sp::Gr2m] = 0.0;
  if constexpr (Model::grain_evaporation::has_gr2p) y[Sp::Gr2p] = 0.0;
}

}  // namespace grain_evap

// The functions below are called for every model; they act only on grain models.

// The grains' first-order sinks at y, added to the rates r and their Jacobian J.
template <class Model>
void add_grain_sinks(const GrainRelease<Model>& rel,
                     const std::array<double, Model::N_sp>& y,
                     std::array<double, Model::N_sp>& r,
                     std::array<double, Model::N_sp * Model::N_sp>& J) {
  if constexpr (Model::has_grain)
    grain_evap::add_sinks_and_jacobian<Model>(rel, y, r, J);
}

// The same, split: the sinks added to the rates alone, and their Jacobian alone.
template <class Model>
void add_grain_sink_rates(const GrainRelease<Model>& rel,
                          const std::array<double, Model::N_sp>& y,
                          std::array<double, Model::N_sp>& r) {
  if constexpr (Model::has_grain) grain_evap::add_sink<Model>(rel, y, 1.0, r);
}

template <class Model>
void add_grain_sink_jacobian(const GrainRelease<Model>& rel,
                             std::array<double, Model::N_sp * Model::N_sp>& J) {
  if constexpr (Model::has_grain) grain_evap::add_sink_jacobian<Model>(rel, J);
}

// out += h * (the grains' sink at y): the change a backward-Euler step of width h
// ending at y makes through the sinks.
template <class Model>
void add_grain_sink_change(const GrainRelease<Model>& rel,
                           const std::array<double, Model::N_sp>& y, double h,
                           std::array<double, Model::N_sp>& out) {
  if constexpr (Model::has_grain) grain_evap::add_sink<Model>(rel, y, h, out);
}

// The start state of a Saha step of width h: the constant rates applied, the grain
// species decayed as a backward-Euler step does (the Saha state sets the
// electrons from the grains' charge).
template <class Model>
void apply_grain_release_to_start(const GrainRelease<Model>& rel, double h,
                                  std::array<double, Model::N_sp>& y) {
  if constexpr (Model::has_grain) grain_evap::release_start<Model>(rel, h, y);
}

// After a call with complete evaporation: the remaining grain species to 0, their
// net charge to the electrons (or to K+, from K).
template <class Model>
void finish_grain_release(const GrainRelease<Model>& rel,
                          std::array<double, Model::N_sp>& y) {
  if constexpr (Model::has_grain) {
    if (rel.complete) grain_evap::remove_grains<Model>(y);
  }
}

}  // namespace arche
