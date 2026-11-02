// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// chemistry.h — portable chemistry module umbrella header
//
// Provides a simplified interface for embedding the chemistry network
// into a fluid/hydrodynamics code. Key components:
//
//   ChemCell<Model>             — per-cell state (ChemState + var cache)
//   ZeroMetalCell               — alias for ChemCell<Nakauchi2019>
//   MetalGrainCell              — alias for ChemCell<Nakauchi2021>
//
//   chem_full_step()            — one chemistry interval with its cooling/heating
//
//   make_zero_metal_table()     — load ZeroMetalTable from data files
//   make_metal_grain_table()    — load MetalGrainTable from data files
//
// Minimal usage (zero-metal network):
// ─────────────────────────────────────────────────────────────────────────
//   #include "solve/chemistry.h"
//   using namespace arche;
//
//   // Load tables once at startup (shared across cells/threads)
//   ZeroMetalTable tbl = make_zero_metal_table("/path/to/data");
//
//   // Allocate per-cell state (heap for OpenMP/GPU parallel loops)
//   auto cell = std::make_unique<ZeroMetalCell>();
//   cell->state.nH = 1.0e2;  // H number density [cm^-3]
//   cell->state.T_K = 100.0;  // temperature [K]
//   // ... initialise cell->state.y[0..22] ...
//   cell->state.mu    = 1.22;   // mean molecular weight [m_p]
//   cell->state.gamma = 5.0/3.0;
//
//   ChemParams params;
//   params.T_rad = 2.725;     // CMB temperature [K]
//   ChemShielding shield;
//   shield.zeta = 1.0e-17;    // CR ionization rate [s^-1]
//   shield.cnt_column = rho * L;  // continuum shielding column [g cm^-2]
//
//   // Per time step:
//   ChemFullRates rates = chem_full_step(*cell, dt, params, shield, tbl);
//   // rates.Lambda_net — interval-mean net cooling [erg g^-1 s^-1];
//   // the caller applies e -= dt * rates.Lambda_net
// ─────────────────────────────────────────────────────────────────────────
//
// Note on stack vs heap:
//   MetalGrainCell contains var[2400] (~19 KB) + state (~0.8 KB).
//   Inside chem_full_step the integrator allocates its (89+1)^2 Newton
//   matrices (~65 KB each) on the stack and its point cache on the heap.
//   For OpenMP: declare MetalGrainCell with thread-private storage or heap.
//   For GPU:    allocate cells in device memory.
// ---------------------------------------------------------------------------

#include "cooling/cooling.h"        // line_cool, cnt_cool
#include "cooling/cooling_grain.h"  // cnt_cool_metal
#include "cooling/cooling_molecular.h"  // line_cool_metal, LineCoolRates, EscapeState
#include "core/state.h"  // ChemState, ChemParams, ChemFullRates, ChemShielding,
#include "kinetics/topology.h"   // ReactionTable + load_pf_tables_h5
#include "solve/chemreact.h"     // compute_chemistry_cooling
#include "solve/integrator.h"    // integ:: error-controlled backward Euler
#include "solve/conservation.h"  // conservation::project / fill_invariants
#include "models/metal_grain/equilibrium.h"  // equichem_metal; pulls in rates.h
#include "models/metal_grain/reactions.h"    // metal_grain::net::init_topology
#include "models/primordial/equilibrium.h"   // equichem, equichem_minimal
#include "models/primordial/reactions.h"     // zero_metal::net::init_topology
                                             // ChemFullRates, phys::
#include "kinetics/reaction_index.h"  // zero_metal::, metal_grain:: constants
#include "model_traits.h"             // Nakauchi2019, Nakauchi2021
#ifdef ARCHE_XRAY
#include "cooling/xray_secondary.h"  // X-ray secondary ionization / heating
#endif
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <limits>
#include <memory>  // std::make_unique (integ::PointCache)

namespace arche {

// ---------------------------------------------------------------------------
// Empty — placeholder for models without escape-probability state.
// ---------------------------------------------------------------------------
struct Empty {};

// ---------------------------------------------------------------------------
// ChemCell<Model>
//
// Bundles a per-cell chemical state with the inter-step reaction rate cache
// (var[2*N_react]).  For models with escape-probability cooling (has_escape),
// an EscapeState member is included; otherwise an empty placeholder.
//
// Thread safety: each thread/cell must own its own ChemCell.
// ReactionTable and ChemParams are read-only and can be shared.
// ---------------------------------------------------------------------------
template <class Model>
struct ChemCell {
  ChemState<Model::N_sp> state{};
  std::array<double, 2 * Model::N_react> var{};
  std::conditional_t<Model::has_escape, EscapeState, Empty> es{};
  // Sub-ulp remainder of the element-conservation projection, one entry per
  // invariant row (solve/conservation.h).  Cleared with var: it belongs to the
  // integration history of this cell, not to its physical state.
  // Because it is history and not state, restoring a saved ChemState does NOT
  // restore it, and a caller that reuses one ChemCell as scratch across many
  // fluid cells carries the previous cell's remainder into the next cell's
  // first step.  Call reset_var() between cells; the un-reset error is bounded
  // by half an ulp of the element total.
  std::array<double, conservation::kMaxRows> cons_carry{};
  // The step size the integrator's error control proposed at the end of the last call
  // and the conditions it was found under (nH, T_K, and the shielding inputs zeta,
  // J_LW21, cnt_column).  The next call starts with it where those conditions changed
  // by at most kStepHintMaxChange (chem_integrate_adaptive); h = 0: none.  History,
  // like cons_carry: reset_var() clears it.
  struct StepHint {
    double h = 0.0, nH = 0.0, T_K = 0.0, zeta = 0.0, J_LW21 = 0.0,
           cnt_column = 0.0;
  } step_hint{};

  void reset_var() noexcept {
    var.fill(0.0);
    cons_carry.fill(0.0);
    step_hint = StepHint{};
  }
};

// Largest relative change of nH, T and the shielding inputs between calls for which the
// last call's proposed step size starts the next call.
constexpr double kStepHintMaxChange = 0.1;

// Convenience aliases matching the two supported networks
using ZeroMetalCell = ChemCell<Nakauchi2019>;
using MetalGrainCell = ChemCell<Nakauchi2021>;

// The rate fields of ChemFullRates that chem_integrate_adaptive averages over the
// accepted steps (weighted by their widths).
inline constexpr double ChemFullRates::*thermal_components[] = {
    &ChemFullRates::Lambda_net, &ChemFullRates::Lambda_line,
    &ChemFullRates::Lambda_cnt, &ChemFullRates::Lambda_chem,
    &ChemFullRates::Gamma_CR,   &ChemFullRates::Gamma_X,
    &ChemFullRates::Lambda_H2,  &ChemFullRates::Lambda_HD,
    &ChemFullRates::Lambda_Lya, &ChemFullRates::Lambda_CO,
    &ChemFullRates::Lambda_OH,  &ChemFullRates::Lambda_H2O,
    &ChemFullRates::Lambda_CII, &ChemFullRates::Lambda_CI,
    &ChemFullRates::Lambda_OI,  &ChemFullRates::Lambda_gr,
    &ChemFullRates::Lambda_gas};

// The output of a completed interval: finite, non-negative abundances, a usable mu and
// gamma, finite rates, and non-negative opacities and grain temperature.
template <class Cell>
bool interval_output_valid(const Cell& cell, const ChemFullRates& rates) {
  if (rates.solver_failed || !rates.conservation_projected) return false;
  for (double y : cell.state.y)
    if (!std::isfinite(y) || y < 0) return false;
  if (!std::isfinite(cell.state.mu) || cell.state.mu <= 0 ||
      !std::isfinite(cell.state.gamma) || cell.state.gamma <= 1)
    return false;
  for (auto field : thermal_components)
    if (!std::isfinite(rates.*field)) return false;
  return std::isfinite(rates.k_gas) && rates.k_gas >= 0 &&
         std::isfinite(rates.k_gr) && rates.k_gr >= 0 &&
         std::isfinite(rates.T_gr_K) && rates.T_gr_K >= 0;
}

// ---------------------------------------------------------------------------
// apply_photo_split<Model>
//
// The operator-split photo channels of one chemistry interval of width dt, applied
// to cell.state.y after each accepted step and its projection: (ARCHE_XRAY) X-ray
// ionisation, and, when lw, LW photodissociation / photodetachment (after a Saha step;
// the network steps carry LW in their rates). Sets rates.Gamma_X.
// ---------------------------------------------------------------------------
template <class Model>
void apply_photo_split(ChemCell<Model>& cell, double dt,
                       const ChemShielding& shield, ChemFullRates& rates,
                       bool lw) {
  using Sp = typename Model::Sp;
  (void)rates;
  // ── Lyman-Werner photodissociation (operator-split) ───────────────────────
  if (lw && shield.J_LW21 > 0.0) {
    auto& y = cell.state.y;
    const integ::LwRates k =
        integ::lw_rates(cell.state.nH, cell.state.T_K, shield);
    const double k_H2_LW = k.H2;
    const double k_HD_LW = k.HD;
    const double k_Hm_LW = k.Hm;

    // H2 + hν(LW) → H + H
    const double dy_H2 = y[Sp::H2] * (1.0 - std::exp(-k_H2_LW * dt));
    y[Sp::H2] -= dy_H2;
    y[Sp::H] += 2.0 * dy_H2;

    // HD + hν(LW) → H + D
    const double dy_HD = y[Sp::HD] * (1.0 - std::exp(-k_HD_LW * dt));
    y[Sp::HD] -= dy_HD;
    y[Sp::H] += dy_HD;
    y[Sp::D] += dy_HD;

    // H- + hν → H + e-
    const double dy_Hm = y[Sp::Hm] * (1.0 - std::exp(-k_Hm_LW * dt));
    y[Sp::Hm] -= dy_Hm;
    y[Sp::H] += dy_Hm;
    y[Sp::e] += dy_Hm;
  }

  // ── X-ray photoionization (operator-split) ────────────────────────────────
#ifdef ARCHE_XRAY
  if (shield.zeta_X > 0.0) {
    auto& y = cell.state.y;
    const double dt_x = dt;
    const double E_X = shield.E_X_eV;
    constexpr double yHe_ab = abundance_ref::yHe;

    const double xe = y[Sp::e] / (1.0 + yHe_ab);
    const double e_cl = std::max(0.0, std::min(1.0, xe));

    const double sig_HI = xray::sigma_HI_Verner96(E_X);
    const double sig_HeI = xray::sigma_HeI_Verner96(E_X);
    const double zeta_X_He =
        (sig_HI > 0.0) ? shield.zeta_X * (sig_HeI / sig_HI) : 0.0;

    const double E_prim_HI = E_X - 13.6;
    const double E_prim_HeI = E_X - 24.6;

    const double Phi_HI_ch1 =
        (E_prim_HI >= 28.0) ? xray::secion_HI2nd(E_prim_HI, e_cl) : 0.0;
    const double Phi_HI_ch2 =
        (E_prim_HeI >= 28.0) ? xray::secion_HI2nd(E_prim_HeI, e_cl) : 0.0;
    const double Phi_HeI_ch3 =
        (E_prim_HI >= 28.0) ? xray::secion_HeI2nd(E_prim_HI, e_cl) : 0.0;
    const double Phi_HeI_ch4 =
        (E_prim_HeI >= 28.0) ? xray::secion_HeI2nd(E_prim_HeI, e_cl) : 0.0;

    const double zeta_eff_H =
        shield.zeta_X * (1.0 + Phi_HI_ch1) + zeta_X_He * Phi_HI_ch2;
    const double zeta_eff_He =
        zeta_X_He * (1.0 + Phi_HeI_ch4) + shield.zeta_X * Phi_HeI_ch3;

    const double dy_H = y[Sp::H] * (1.0 - std::exp(-zeta_eff_H * dt_x));
    const double dy_He = y[Sp::He] * (1.0 - std::exp(-zeta_eff_He * dt_x));

    y[Sp::H] -= dy_H;
    y[Sp::Hp] += dy_H;
    y[Sp::He] -= dy_He;
    y[Sp::Hep] += dy_He;
    y[Sp::e] += dy_H + dy_He;

    // He+ -> He++ X-ray ionization is a composable channel: it needs Sp::Hepp,
    // so models without He++ (has_he_pp == false, e.g. the compact minimal
    // network) omit it and keep the dominant HI/HeI channels above/below.
    [[maybe_unused]] double dy_Hep = 0.0;
    if constexpr (Model::cooling::has_he_pp) {
      if (E_X > 54.42 && y[Sp::Hep] > 0.0) {
        const double sig_HeII = xray::sigma_HeII_Verner96(E_X);
        const double zeta_X_HeII =
            (sig_HI > 0.0) ? shield.zeta_X * (sig_HeII / sig_HI) : 0.0;
        dy_Hep = y[Sp::Hep] * (1.0 - std::exp(-zeta_X_HeII * dt_x));
        y[Sp::Hep] -= dy_Hep;
        y[Sp::Hepp] += dy_Hep;
        y[Sp::e] += dy_Hep;
      }
    }

    const double Eh_HI =
        (E_prim_HI > 0.0) ? xray::secion_heat(E_prim_HI, e_cl) : 0.0;
    const double Eh_HeI =
        (E_prim_HeI > 0.0) ? xray::secion_heat(E_prim_HeI, e_cl) : 0.0;

    const double pre = phys::eV_to_erg / ((1.0 + 4.0 * yHe_ab) * phys::m_p);
    rates.Gamma_X =
        (shield.zeta_X * y[Sp::H] * Eh_HI + zeta_X_He * y[Sp::He] * Eh_HeI) *
        pre;

    if constexpr (Model::cooling::has_he_pp) {
      const double E_prim_HeII = E_X - 54.42;
      if (E_X > 54.42 && dy_Hep > 0.0) {
        const double sig_HeII = xray::sigma_HeII_Verner96(E_X);
        const double zeta_X_HeII =
            (sig_HI > 0.0) ? shield.zeta_X * (sig_HeII / sig_HI) : 0.0;
        rates.Gamma_X += zeta_X_HeII * y[Sp::Hep] * E_prim_HeII * pre;
      }
    }
  }
#endif  // ARCHE_XRAY
}

// ---------------------------------------------------------------------------
// chem_integrate_adaptive<Model>
//
// The complete interval dt by the error-controlled backward Euler of
// solve/integrator.h: species and internal energy in one Newton, step doubling for
// the local error, weights ChemParams::rtol / atol / rtol_T. Each accepted step is
// projected onto the element/charge invariants and then gets the operator-split
// photo channels; the returned rates are the h-weighted means over the accepted
// steps, so e0 - dt*Lambda_net is the end state's EOS energy (to the Newton
// tolerance) when no photo channel acts. The internal steps advance their own
// temperature, but cell.state.T_K is returned at its input value (the caller owns the
// energy update). On failure the cell is restored to its input.
// Grain models: every step solves the grain temperature at its iterates, seeded from
// the one at the step's start; the returned T_gr_K, k_gr, k_gas, tau_cnt and esc_cnt
// are those of the interval's end state, and the cell's escape state is the end
// state's solution. The continuum depth is ChemShielding::cnt_column times the opacity
// at each evaluated state (0 = optically thin). The grains evaporate over the interval
// as solve/grain_evaporation.h describes: the release is planned from the input state
// and enters every step, and the projection's targets include it.
// ---------------------------------------------------------------------------
template <class Model>
ChemFullRates chem_integrate_adaptive(
    ChemCell<Model>& cell, double dt, const ChemParams& params,
    const ChemShielding& shield,
    const ReactionTable<Model::N_sp, Model::N_react>& tbl) {
  constexpr int N_sp = Model::N_sp;
  if (!std::isfinite(params.T_rad) || params.T_rad <= 0.0) {
    throw std::invalid_argument(
        "chem_full_step: params.T_rad must be set to a finite positive value "
        "[K]");
  }
  if (params.use_saha && !(std::isfinite(params.saha_threshold_nH) &&
                           params.saha_threshold_nH > 0.0)) {
    throw std::invalid_argument(
        "chem_full_step: params.saha_threshold_nH must be finite and positive "
        "[cm^-3] when params.use_saha is true");
  }
  if (!(params.rtol > 0.0) || !(params.atol > 0.0) || !(params.rtol_T > 0.0) ||
      params.internal_max_steps < 1 || !(params.internal_max_step > 0)) {
    throw std::invalid_argument(
        "chem_full_step: rtol, atol, rtol_T and internal_max_step must be "
        "positive and "
        "internal_max_steps >= 1");
  }
  if (!std::isfinite(shield.cnt_column) || shield.cnt_column < 0.0) {
    throw std::invalid_argument(
        "chem_full_step: shield.cnt_column must be finite and non-negative [g "
        "cm^-2]");
  }
  if constexpr (Model::has_grain) {
    if (!std::isfinite(params.T_gr_K) || params.T_gr_K <= 0.0) {
      throw std::invalid_argument(
          "chem_full_step: params.T_gr_K must be set to a finite positive "
          "value [K] for "
          "metal_grain");
    }
  }
  // The integrator keeps every abundance non-negative and accepts a step only if it
  // ends that way, so a negative input could never be advanced.
  for (int i = 0; i < N_sp; ++i) {
    if (!std::isfinite(cell.state.y[i]) || cell.state.y[i] < 0.0) {
      throw std::invalid_argument(
          "chem_full_step: every abundance y must be finite and non-negative");
    }
  }
  ChemParams p = params;
  p.zeta = shield.zeta;
#ifdef ARCHE_XRAY
  p.zeta_X = shield.zeta_X;
#endif
  const integ::Tolerance tol{params.rtol, params.atol, params.rtol_T};
  const ChemCell<Model> cell_in = cell;
  const double nH = cell.state.nH;
  SolveStatus st;
  auto fail = [&](SolveFailure reason) {
    cell = cell_in;
    ChemFullRates r;
    r.solver_failed = true;
    st.reason = reason;
    st.accepted_substeps = 0;
    st.completed_time = 0.0;
    r.status = st;
    return r;
  };
  if (!std::isfinite(dt) || dt <= 0.0) return fail(SolveFailure::no_progress);
  // Quasi-neutral input.  Over a fluid cell of size L a plasma can hold a net charge of
  // only (lambda_D / L)^2 of its electrons, decades below double precision, so a
  // charged input is a defect of the caller's state (transport, floors, an initial
  // condition), not physics.  Every step projects its end state onto zero net charge
  // and is refused when that repair exceeds the error weights or the projection's
  // limit (each species shifted by at most conservation::kMaxRelShift of itself, so a
  // charge residual above kMaxRelShift / max|q| of SUM q^2 y -- a multiply charged
  // species counts even when it is absent); an input imbalance makes the repair
  // independent of the step width, and no step could be accepted.  Any charged input
  // is therefore tried against the same element-preserving projection here:
  //   the projection accepts it: an imbalance whose neutralisation exceeds the error
  //     weights is projected once here, and SolveStatus reports it; a smaller one is
  //     left to the steps' own projection, so a neutral input is integrated exactly as
  //     before;
  //   the projection declines and |net| <= atol: the imbalance is moved onto the
  //     electrons (y_e += net), which carry no element; the change is inside every
  //     error weight;
  //   otherwise the input is rejected.
  if (tbl.n_invariants > 0 && tbl.charge_invariant_row >= 0) {
    const double* q = tbl.invariants.data() +
                      static_cast<std::size_t>(tbl.charge_invariant_row) * N_sp;
    double net = 0.0, carried = 0.0;
    for (int i = 0; i < N_sp; ++i) {
      net += q[i] * cell.state.y[i];
      carried += std::abs(q[i]) * cell.state.y[i];
    }
    st.input_charge_imbalance = carried > 0.0 ? std::abs(net) / carried : 0.0;
    std::array<double, N_sp> yn = cell.state.y;
    if (net != 0.0) {
      std::array<double, conservation::kMaxRows> no_carry{};
      if (!conservation::project<N_sp>(tbl.invariants.data(), tbl.n_invariants,
                                       tbl.charge_invariant_row, cell.state.y,
                                       yn, no_carry)) {
        using Sp = typename Model::Sp;
        const double ye = cell.state.y[Sp::e] + net;
        if (std::abs(net) <= tol.atol && ye >= 0.0) {
          cell.state.y[Sp::e] = ye;
        } else {
          throw std::invalid_argument(
              "chem_full_step: the input's net charge (" +
              std::to_string(st.input_charge_imbalance) +
              " of the charge carried) cannot be neutralised within "
              "conservation::kMaxRelShift; the input state must be "
              "charge-neutral");
        }
      }
    }
    if (std::abs(net) > tol.atol) {
      for (int i = 0; i < N_sp && !st.input_neutralized; ++i) {
        const double w =
            tol.atol +
            tol.rtol * std::max(std::abs(yn[i]), std::abs(cell.state.y[i]));
        if (std::abs(yn[i] - cell.state.y[i]) > w) st.input_neutralized = true;
      }
      if (st.input_neutralized) cell.state.y = yn;
    }
  }

  ChemFullRates acc;     // sums of h * rate over the accepted steps
  integ::Radiative end;  // the last accepted step's end state
  end.T_gr = params.T_gr_K;
  integ::EscOf<Model> es{};
  if constexpr (Model::has_escape) es = cell.es;
  // Point evaluations shared by the three steps of one attempt (integrator.h).
  auto cache = std::make_unique<integ::PointCache<Model>>();
  // Grain evaporation over the interval (solve/grain_evaporation.h), from its input
  // state: the grain temperature at the input temperature, seeded from params.
  GrainRelease<Model> release;
  if constexpr (Model::has_grain) {
    integ::Radiative r_in;
    integ::continuum_at<Model>(nH, cell.state.T_K, p, shield, r_in);
    release = plan_grain_release<Model>(nH, r_in.T_gr, dt, cell.state.y);
  }
  const GrainRelease<Model>* rel = release.active ? &release : nullptr;
  double t = 0.0;
  double h = std::min(dt, params.internal_max_step);
  // Start with the step size the error control proposed at the end of this cell's last
  // call, where the conditions it was found under have not changed much.
  {
    const auto& sh = cell.step_hint;
    auto close = [](double a, double b) {
      return std::abs(a - b) <=
             kStepHintMaxChange * std::max(std::abs(a), std::abs(b));
    };
    if (sh.h > 0.0 && close(sh.nH, nH) && close(sh.T_K, cell.state.T_K) &&
        close(sh.zeta, shield.zeta) && close(sh.J_LW21, shield.J_LW21) &&
        close(sh.cnt_column, shield.cnt_column))
      h = std::min(h, sh.h);
  }
  // The error control's step size, before the interval end cuts it.
  double h_ctrl = h;
  // The latest accepted step, for the temperature seed of the next one.
  double h_prev = 0.0, T_prev = 0.0;
  // The Newton matrix of the first half step, which the second (of the same width)
  // starts from.
  auto lu_half = std::make_unique<integ::DenseLU<N_sp + 1>>();
  // The interval's resolution: the loop ends once less than this remains.  The error
  // control may take the step below it -- a stiff start can need its first steps many
  // decades below dt (e.g. ~1e-11 s at n_H = 1e19 cm^-3) -- but a step that a rejection
  // for any other cause (a failed solve, a refused projection, a non-finite error)
  // shortened below it ends the call.  A cause that does not shrink with the step width
  // (a defect of the state rather than of the step) therefore ends the call after at
  // most log4(1 / (8 eps)) ~ 25 shortenings.
  const double h_min = 8.0 * std::numeric_limits<double>::epsilon() * dt;
  // Whether the last attempt was refused for a cause other than a finite local error.
  bool failure_shortened = false;
  auto reject = [&](SolveFailure why) {
    ++st.rejected_substeps;
    ++st.rejections_by_reason[static_cast<int>(why)];
    if (st.last_failure_run > 0 && why == st.last_failure) {
      ++st.last_failure_run;
    } else {
      st.last_failure_run = 1;
      st.last_failure_run_h = h;
    }
    st.last_failure = why;
  };
  while (dt - t > h_min) {
    if (st.attempted_substeps >= params.internal_max_steps)
      return fail(SolveFailure::interval_budget);
    h_ctrl = h;
    h = std::min(h, dt - t);
    if (!(h > (failure_shortened ? h_min : 0.0)) || t + h == t)
      return fail(SolveFailure::no_progress);
    ++st.attempted_substeps;
    const std::array<double, N_sp> y0 = cell.state.y;
    const double T0 = cell.state.T_K;
    // The grain temperature at the attempt's start: the seed of the steps' solves, and
    // the value the Saha test's rate coefficients use.
    ChemParams pa = p;
    if constexpr (Model::has_grain) {
      integ::Radiative r0;
      integ::continuum_at<Model>(nH, T0, p, shield, r0);
      pa.T_gr_K = r0.T_gr;
    }
    // Saha or network for this attempt, decided at its start.
    double shift = -1.0;
    const bool saha =
        integ::saha_regime<Model>(nH, y0, T0, pa, tbl, cell.var, &shift);
    if (shift >= 0.0) st.saha_nonthermal_shift = shift;
    cache->clear();
    auto step = [&](const std::array<double, N_sp>& ys, double Ts, double hs,
                    const ChemParams& ps, const integ::EscOf<Model>& es0,
                    bool end_rates, const integ::StepGuess<Model>& g,
                    integ::DenseLU<N_sp + 1>* lu) {
      return saha
                 ? integ::saha_step<Model>(nH, ys, Ts, hs, ps, shield, es0, tbl,
                                           tol, cache.get(), rel, g.T_seed)
                 : integ::be_step<Model>(nH, ys, Ts, hs, ps, shield, es0, tbl,
                                         tol, cache.get(), end_rates, rel, &g,
                                         lu);
    };
    // Where the steps start their iterations (integ::StepGuess): the full step's
    // temperature root from the trend of the last accepted step, the first half's from
    // the midpoint of the start and the full step's end, and the second half -- which
    // ends where the full step ends -- from the full step's end state, with the first
    // half's Newton matrix (the same width).
    integ::StepGuess<Model> ga, gb1, gb2;
    if (h_prev > 0.0) ga.T_seed = T0 + (T0 - T_prev) * (h / h_prev);
    // The full step serves only the error estimate: its end state, not its rates.
    auto a = step(y0, T0, h, pa, es, false, ga, nullptr);
    if (a.ok) gb1.T_seed = 0.5 * (T0 + a.T);
    auto b1 = step(y0, T0, 0.5 * h, pa, es, true, gb1, lu_half.get());
    integ::StepResult<Model> b2;
    if (a.ok && b1.ok) {
      ChemParams pm = pa;  // the second half starts at the first half's end
      pm.T_gr_K = b1.rad_end.T_gr;
      gb2.y = &a.y;
      gb2.T = a.T;
      gb2.lu = b1.lu_valid;
      gb2.T_seed = a.T;
      b2 = step(b1.y, b1.T, 0.5 * h, pm, b1.es_end, true, gb2, lu_half.get());
    }
    const int newton = a.newton + b1.newton + b2.newton;
    st.newton_total += newton;
    st.newton_iterations =
        std::max({st.newton_iterations, a.newton, b1.newton, b2.newton});
    // Steps that needed integrator.h's fixed-temperature fallback.
    st.temperature_iterations +=
        int(a.fallback) + int(b1.fallback) + int(b2.fallback);
    if (!(a.ok && b1.ok && b2.ok)) {
      reject(!a.ok ? a.reason : (!b1.ok ? b1.reason : b2.reason));
      failure_shortened = true;
      h *= integ::kShrinkFail;
      continue;
    }
    const double err = integ::weighted_rms<N_sp>(a.y, a.T, b2.y, b2.T, tol);
    st.update_norm = err;
    if (!(err <= 1.0)) {
      reject(SolveFailure::output_rejected);
      failure_shortened = !std::isfinite(err);
      h *= std::isfinite(err)
               ? std::max(integ::kShrinkMin,
                          std::min(0.9, integ::kSafety / std::sqrt(err)))
               : integ::kShrinkFail;
      continue;
    }
    // Accept the pair's end state: project it onto the invariants of the step start.
    // The repair must itself be inside the tolerance: a step whose solve broke an element
    // or the charge by more than the weights allow is refused (and shortened), not
    // silently repaired -- the projection's own 1e-2 limit is far looser than rtol.
    std::array<double, N_sp> y1 = b2.y;
    // The grain release changes y outside the reactions: the targets are the start
    // state's totals plus what the release added over the two half steps.
    std::array<double, N_sp> y0_target = y0;
    if (rel) {
      for (int i = 0; i < N_sp; ++i) y0_target[i] += h * rel->src[i];
      add_grain_sink_change<Model>(*rel, b1.y, 0.5 * h, y0_target);
      add_grain_sink_change<Model>(*rel, b2.y, 0.5 * h, y0_target);
    }
    const auto carry_before = cell.cons_carry;
    bool projected = conservation::project<N_sp>(
        tbl.invariants.data(), tbl.n_invariants, tbl.charge_invariant_row,
        y0_target, y1, cell.cons_carry);
    if (projected) {
      for (int i = 0; i < N_sp && projected; ++i) {
        const double w =
            tol.atol + tol.rtol * std::max(std::abs(y1[i]), std::abs(b2.y[i]));
        if (std::abs(y1[i] - b2.y[i]) > w) projected = false;
      }
      if (!projected) cell.cons_carry = carry_before;
    }
    if (!projected) {
      reject(tbl.n_invariants == 0 ? SolveFailure::invariants_unavailable
                                   : SolveFailure::projection_rejected);
      failure_shortened = true;
      h *= integ::kShrinkFail;
      continue;
    }
    for (int i = 0; i < N_sp; ++i)
      st.projection_max_relative = std::max(
          st.projection_max_relative,
          std::abs(y1[i] - b2.y[i]) / std::max(std::abs(b2.y[i]), 1.e-300));
    cell.state.y = y1;
    cell.var = b2.var;
    ChemFullRates sr;  // this step's rates: mean of its two half steps
    integ::add_step_rates<Model>(b1, 0.5, sr);
    integ::add_step_rates<Model>(b2, 0.5, sr);
    // Operator-split photo channels at fixed specific internal energy (plus X-ray
    // heating); the end temperature follows from the EOS when they changed anything.
    // LW is split only after a Saha step: the network steps carry it in their rates.
    double T1 = b2.T;
    const double e_before = integ::eos_energy<Model>(cell.state.y, T1);
    const auto y_before_split = cell.state.y;
    apply_photo_split<Model>(cell, h, shield, sr, saha);
    if (cell.state.y != y_before_split || sr.Gamma_X != 0.0) {
      T1 = thermo::temperature_from_specific_internal_energy<Model>(
          cell.state.y, e_before + h * sr.Gamma_X);
      if (!std::isfinite(T1) || !(T1 > 0.0))
        return fail(SolveFailure::nonfinite);
    }
    cell.state.T_K = T1;
    for (auto field : thermal_components) acc.*field += h * (sr.*field);
    end = b2.rad_end;
    es = b2.es_end;
    if constexpr (Model::has_grain) p.T_gr_K = end.T_gr;
    t += h;
    ++st.accepted_substeps;
    failure_shortened = false;
    st.last_failure_run = 0;
    st.last_failure_run_h = 0.0;
    st.internal_step = h;
    h_prev = h;
    T_prev = T0;
    const bool cut = h < h_ctrl;  // the interval end shortened this step
    h *= std::min(integ::kGrowMax,
                  integ::kSafety / std::sqrt(std::max(err, 1.0e-10)));
    h = std::min(h, params.internal_max_step);
    if (cut) h = std::max(h, h_ctrl);
    cell.step_hint = {
        h, nH, cell.state.T_K, shield.zeta, shield.J_LW21, shield.cnt_column};
  }

  ChemFullRates r;
  for (auto field : thermal_components) r.*field = acc.*field / dt;
  r.Lambda_net = r.Lambda_line + r.Lambda_cnt + r.Lambda_chem - r.Gamma_CR
#ifdef ARCHE_XRAY
                 - r.Gamma_X
#endif
      ;
  r.k_gas = end.k_gas;
  r.k_gr = end.k_gr;
  r.tau_cnt = end.tau;
  r.esc_cnt = end.esc;
  r.T_gr_K = Model::has_grain ? end.T_gr : params.T_gr_K;
  if constexpr (Model::has_escape) cell.es = es;
  r.conservation_projected = true;
  // The API leaves the gas temperature to the caller (ChemFullRates: the input
  // temperature stays fixed; the caller applies -dt*Lambda_net to its energy), so the
  // integrator's end temperature is not written back.
  cell.state.T_K = cell_in.state.T_K;
  if (rel)
    finish_grain_release<Model>(*rel, cell.state.y);  // complete evaporation
  {
    const auto eos = thermo::make_eos_composition<Model>(cell.state.y);
    cell.state.mu = thermo::mean_molecular_weight(eos);
    cell.state.gamma = thermo::adiabatic_index(eos, cell.state.T_K);
  }
  st.reason = SolveFailure::none;
  st.completed_time = dt;
  r.status = st;
  if (!interval_output_valid(cell, r))
    return fail(SolveFailure::output_rejected);
  return r;
}

// ---------------------------------------------------------------------------
// chem_full_step<Model>
//
// Advances the chemistry of one cell over dt and returns the interval-mean
// cooling/heating breakdown (chem_integrate_adaptive). The step control, the
// retries and the rollback on failure all live here, never in the caller.
// ---------------------------------------------------------------------------
template <class Model>
ChemFullRates chem_full_step(
    ChemCell<Model>& cell, double dt, const ChemParams& params,
    const ChemShielding& shield,
    const ReactionTable<Model::N_sp, Model::N_react>& tbl) {
  return chem_integrate_adaptive<Model>(cell, dt, params, shield, tbl);
}

// ---------------------------------------------------------------------------
// make_zero_metal_table()
// ---------------------------------------------------------------------------
inline ZeroMetalTable make_zero_metal_table(
    const std::string& prim_chem_table) {
  ZeroMetalTable tbl;
  zero_metal::net::init_topology(tbl);      // reaction topology from C++ source
  load_pf_tables_h5(tbl, prim_chem_table);  // partition functions from HDF5
  // Element/charge rows for the conservation projection, verified against the
  // topology just loaded (solve/conservation.h).
  conservation::fill_invariants<zero_metal::Species>(tbl);
  return tbl;
}

// ---------------------------------------------------------------------------
// make_minimal_table() — compact 15-species / 33-reaction Nakauchi2019_Minimal
//   topology (built from the full keep-set via the catalog remap) with the
//   partition-function tables placed into the compact slots by canonical
//   species (the full PF tables are loaded once and copied across).
// ---------------------------------------------------------------------------
inline zero_metal_minimal::MinimalTable make_minimal_table(
    const std::string& prim_chem_table) {
  zero_metal_minimal::MinimalTable tbl;
  zero_metal_minimal::build_topology(tbl);  // compact reaction topology

  ZeroMetalTable full_tmp;
  load_pf_tables_h5(full_tmp, prim_chem_table);  // full PF tables (by species)
  for (int i = 0; i < zero_metal_minimal::N_sp; ++i) {
    int full_i =
        zero_metal::Species::local(zero_metal_minimal::Species::canonical(i));
    tbl.pf_table[i] = full_tmp.pf_table[full_i];
  }
  conservation::fill_invariants<zero_metal_minimal::Species>(tbl);
  return tbl;
}

// ---------------------------------------------------------------------------
// make_metal_grain_table()
// ---------------------------------------------------------------------------
inline MetalGrainTable make_metal_grain_table(
    const std::string& metal_chem_table) {
  MetalGrainTable tbl;
  metal_grain::net::init_topology(tbl);  // reaction topology from C++ source
  // Partition functions from HDF5 for the species 0..88; the dust reservoirs
  // after them have none.
  load_pf_tables_h5(tbl, metal_chem_table, metal_grain::N_sp_table);
  // Element/charge rows for the conservation projection, verified against the
  // topology just loaded (solve/conservation.h).  core/species_composition.h
  // registers all 95 catalog species, so the rows are derived and checked
  // rather than skipped.
  conservation::fill_invariants<metal_grain::Species>(tbl);
  return tbl;
}

// ---------------------------------------------------------------------------
// make_metal_minimal_table() — compact 44-species Nakauchi2021_Minimal topology
//   built from the full keep-set via the catalog remap, with the
//   partition-function tables placed into the compact slots by canonical
//   species.  The full PF-loaded table is retained (aux_full_metal): the compact
//   rate path evaluates the kept reactions' laws on it and the compact Saha the
//   compact species' partition functions (see compute_base_rates /
//   MinimalMetalSaha).
// ---------------------------------------------------------------------------
inline metal_grain_minimal::MinimalTable make_metal_minimal_table(
    const MetalGrainTable& full) {
  metal_grain_minimal::MinimalTable tbl;
  metal_grain_minimal::build_topology(tbl);  // compact reaction topology
  for (int i = 0; i < metal_grain_minimal::N_sp; ++i)
    tbl.pf_table[i] = full.pf_table[metal_grain::Species::local(
        metal_grain_minimal::Species::canonical(i))];
  tbl.aux_full_metal = std::make_shared<MetalGrainTable>(full);
  conservation::fill_invariants<metal_grain_minimal::Species>(tbl);
  return tbl;
}

}  // namespace arche
