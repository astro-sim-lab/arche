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
//   chem_step()                 — single time-step chemistry wrapper
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
//   params.zeta  = 1.0e-17;   // CR ionization rate [s^-1]
//   params.T_rad = 2.725;     // CMB temperature [K]
//
//   // Per time step:
//   ChemRates rates = chem_step(*cell, dt, params, tbl);
//   // rates.Lambda_chem — chemistry cooling [erg g^-1 s^-1]
//   // rates.Gamma_CR — CR heating        [erg g^-1 s^-1]
// ─────────────────────────────────────────────────────────────────────────
//
// Note on stack vs heap:
//   MetalGrainCell contains var[2400] (~19 KB) + state (~0.8 KB).
//   Inside chem_step, chemreact<Model> also allocates A_mat[89×89]
//   and dr_fdy[89×89] (~63 KB each) on the stack.
//   For OpenMP: declare MetalGrainCell with thread-private storage or heap.
//   For GPU:    allocate cells in device memory.
// ---------------------------------------------------------------------------

#include "cooling/cooling.h"        // line_cool, cnt_cool
#include "cooling/cooling_grain.h"  // cnt_cool_metal
#include "cooling/cooling_molecular.h"  // line_cool_metal, LineCoolRates, EscapeState
#include "core/state.h"  // ChemState, ChemParams, ChemRates, ChemShielding,
#include "kinetics/topology.h"  // ReactionTable + load_pf_tables_h5
#include "models/metal_grain/equilibrium.h"  // equichem_metal; pulls in rates.h
#include "models/metal_grain/reactions.h"    // metal_grain::net::init_topology
#include "models/primordial/equilibrium.h"   // equichem, equichem_minimal
#include "models/primordial/reactions.h"     // zero_metal::net::init_topology
#include "solve/chemreact.h"  // chemreact, chemcool, compute_chemistry_cooling
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

namespace arche {

#ifdef ARCHE_SUBCYCLE
constexpr int kSubcycleMaxLevel = 6;
#endif

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

  void reset_var() noexcept { var.fill(0.0); }
};

// Convenience aliases matching the two supported networks
using ZeroMetalCell = ChemCell<Nakauchi2019>;
using MetalGrainCell = ChemCell<Nakauchi2021>;

// ---------------------------------------------------------------------------
// chem_step<Model>
//
// Advances the chemistry network by one time step dt for a single cell.
// ---------------------------------------------------------------------------
template <class Model>
ChemRates chem_step(ChemCell<Model>& cell, double dt, const ChemParams& params,
                    const ReactionTable<Model::N_sp, Model::N_react>& tbl) {
  if (!std::isfinite(params.T_rad) || params.T_rad <= 0.0) {
    throw std::invalid_argument(
        "chem_step: params.T_rad must be set to a finite positive value [K]");
  }
  if constexpr (Model::has_grain) {
    if (!std::isfinite(params.T_gr_K) || params.T_gr_K <= 0.0) {
      throw std::invalid_argument(
          "chem_step: params.T_gr_K must be set to a finite positive value [K] "
          "for metal_grain");
    }
  }

  auto r = chemcool<Model>(cell.state.nH, cell.state.T_K, cell.state.y, dt,
                           cell.var, tbl, params);
  cell.state.mu = r.mu;
  cell.state.gamma = r.gamma;
  const double Lambda_chem = r.Lambda_chem;

  // CR ionization heating
  double rate_CR = 0.0;
  for (int i = Model::cr_heat_var_begin; i < Model::cr_heat_var_end; ++i)
    rate_CR += cell.var[i];

  constexpr double yHe = abundance_ref::yHe;
  double Gamma_CR =
      3.4 * phys::eV_to_erg * rate_CR / ((1.0 + 4.0 * yHe) * phys::m_p);

  return ChemRates{Lambda_chem, Gamma_CR};
}

// ---------------------------------------------------------------------------
// chem_full_step<Model>
//
// Advances the chemistry network by one time step AND computes all
// cooling/heating rate components for a single cell.
// ---------------------------------------------------------------------------
template <class Model>
ChemFullRates chem_full_step(
    ChemCell<Model>& cell, double dt, const ChemParams& params,
    const ChemShielding& shield,
    const ReactionTable<Model::N_sp, Model::N_react>& tbl) {
  constexpr int N_sp = Model::N_sp;
  constexpr int N_react = Model::N_react;
  using Sp = typename Model::Sp;  // model-specific species names: y[Sp::H] etc.

  if (!std::isfinite(params.T_rad) || params.T_rad <= 0.0) {
    throw std::invalid_argument(
        "chem_full_step: params.T_rad must be set to a finite positive value "
        "[K]");
  }
  if constexpr (Model::has_grain) {
    if (!std::isfinite(params.T_gr_K) || params.T_gr_K <= 0.0) {
      throw std::invalid_argument(
          "chem_full_step: params.T_gr_K must be set to a finite positive "
          "value [K] for metal_grain");
    }
  }

  // ChemParams is read-only for the kernel (state.h). The shield→params
  // transcription and the in-step grain-temperature update happen on this
  // private working copy, leaving the caller's `params` untouched.
  ChemParams p = params;
  p.zeta = shield.zeta;
#ifdef ARCHE_XRAY
  p.zeta_X = shield.zeta_X;
#endif

  ChemFullRates rates;

  // ── Line cooling ──────────────────────────────────────────────────────────
  if constexpr (Model::has_escape) {
    // Metal: H2/HD/Lya primitive part from zero-metal line_cool
    double Lambda_line_prim;
    LineCoolRates mcool{};
    line_cool<Sp, N_sp>(cell.state.y, shield.Nc_H2, shield.Nc_HD, cell.state.nH,
                        cell.state.T_K, p.T_rad, shield.tau_cnt,
                        Lambda_line_prim, mcool.H2, mcool.HD, mcool.Lya);
    const auto& y = cell.state.y;
    line_cool_metal(cell.state.nH, cell.state.T_K, p.T_rad, y[Sp::H], y[Sp::H2],
                    y[Sp::e], y[Sp::Hp], y[Sp::He], y[Sp::C], y[Sp::Cp],
                    y[Sp::O], y[Sp::OH], y[Sp::CO], y[Sp::H2O], shield.tau_cnt,
                    shield.Nc_CO, shield.Nc_OH, shield.Nc_H2O, shield.Nc_CII,
                    shield.Nc_CI, shield.Nc_OI, cell.es, mcool);
    rates.Lambda_H2 = mcool.H2;
    rates.Lambda_HD = mcool.HD;
    rates.Lambda_Lya = mcool.Lya;
    rates.Lambda_CO = mcool.CO;
    rates.Lambda_OH = mcool.OH;
    rates.Lambda_H2O = mcool.H2O;
    rates.Lambda_CII = mcool.CII;
    rates.Lambda_CI = mcool.CI;
    rates.Lambda_OI = mcool.OI;
    rates.Lambda_line = mcool.total();
  } else {
    // Zero-metal: H2 / HD / Ly-alpha only
    double Lambda_line;
    line_cool<Sp, N_sp>(cell.state.y, shield.Nc_H2, shield.Nc_HD, cell.state.nH,
                        cell.state.T_K, p.T_rad, shield.tau_cnt, Lambda_line,
                        rates.Lambda_H2, rates.Lambda_HD, rates.Lambda_Lya);
    rates.Lambda_line = Lambda_line;
  }

  // ── Continuum cooling ─────────────────────────────────────────────────────
  if constexpr (Model::has_grain) {
    // cnt_cool_metal warm-starts the grain-temperature solve from p.T_gr_K
    // and overwrites it; the updated value then feeds the grain reaction
    // rates inside chemcool below, so it must be the same `p`.
    cnt_cool_metal(cell.state.nH, cell.state.T_K, p.T_rad, shield.tau_cnt,
                   shield.esc_cnt, p.Z_metal, p.T_gr_K, rates.k_gr, rates.k_gas,
                   rates.Lambda_gr, rates.Lambda_gas, rates.Lambda_cnt);
  } else {
    cnt_cool(cell.state.nH, cell.state.T_K, p.T_rad, shield.esc_cnt,
             rates.k_gas, rates.Lambda_gas, rates.Lambda_cnt);
  }
  rates.T_gr_K = p.T_gr_K;  // export solved grain temperature to the caller

  // ── Chemistry ─────────────────────────────────────────────────────────────
  double Lambda_ch = 0.0;
#ifdef ARCHE_SUBCYCLE
  {
    const auto y_save = cell.state.y;
    const auto mu_save = cell.state.mu;
    const auto gam_save = cell.state.gamma;
    bool all_ok = false;
    for (int level = 0; level <= kSubcycleMaxLevel; ++level) {
      const int n_sub = 1 << level;
      const double dt_sub = dt / n_sub;
      cell.state.y = y_save;
      cell.state.mu = mu_save;
      cell.state.gamma = gam_save;
      double Lambda_ch_acc = 0.0;
      bool attempt_ok = true;
      for (int isub = 0; isub < n_sub; ++isub) {
        auto r = chemcool<Model>(cell.state.nH, cell.state.T_K, cell.state.y,
                                 dt_sub, cell.var, tbl, p);
        if (!r.converged) {
          attempt_ok = false;
          break;
        }
        cell.state.mu = r.mu;
        cell.state.gamma = r.gamma;
        Lambda_ch_acc += r.Lambda_chem;
      }
      if (attempt_ok) {
        Lambda_ch = Lambda_ch_acc / n_sub;
        all_ok = true;
        break;
      }
    }
    if (!all_ok) {
      cell.state.y = y_save;
      cell.state.mu = mu_save;
      cell.state.gamma = gam_save;
      rates.solver_failed = true;
      return rates;
    }
  }
#else
  auto r = chemcool<Model>(cell.state.nH, cell.state.T_K, cell.state.y, dt,
                           cell.var, tbl, p);
  // Report a surface solver failure even without subcycling so the caller can
  // react instead of silently consuming a rolled-back step.
  if (!r.converged) {
    rates.solver_failed = true;
    return rates;
  }
  cell.state.mu = r.mu;
  cell.state.gamma = r.gamma;
  Lambda_ch = r.Lambda_chem;
#endif
  rates.Lambda_chem = Lambda_ch;

  // ── Lyman-Werner photodissociation (operator-split) ───────────────────────
  if (shield.J_LW21 > 0.0) {
    auto& y = cell.state.y;
    const double T_K = cell.state.T_K;
    const double nH = cell.state.nH;

    const double b5_H2 = std::sqrt(phys::k_B * T_K / phys::m_p) / 1.0e5;
    const double b5_HD = std::sqrt(phys::k_B * T_K / (1.5 * phys::m_p)) / 1.0e5;

    const double log10_nH = std::log10(nH > 0.0 ? nH : 1.0);
    const double log10_TK = std::log10(T_K > 0.0 ? T_K : 1.0);
    const double A1_H2 = 0.8711 * log10_TK - 1.928;
    const double A2_H2 = -0.9639 * log10_TK + 3.892;
    const double alpha_H2 = A1_H2 * std::pow(10.0, -0.2856 * log10_nH) + A2_H2;
    const double x_H2 = shield.Nc_H2 / 5.0e14;
    const double fsh_H2 = 0.965 / std::pow(1.0 + x_H2 / b5_H2, alpha_H2) +
                          0.035 / std::sqrt(1.0 + x_H2) *
                              std::exp(-8.5e-4 * std::sqrt(1.0 + x_H2));
    const double k_H2_LW = 1.38e-12 * shield.J_LW21 * fsh_H2;

    const double x_HD = shield.Nc_HD / 5.0e14;
    const double fsh_HD = 0.965 / std::pow(1.0 + x_HD / b5_HD, 2.0) +
                          0.035 / std::sqrt(1.0 + x_HD) *
                              std::exp(-8.5e-4 * std::sqrt(1.0 + x_HD));
    const double k_HD_LW = 1.38e-12 * shield.J_LW21 * fsh_HD;

    const double k_Hm_LW = 1.1e-10 * shield.J_LW21;

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

  // ── CR heating ────────────────────────────────────────────────────────────
  double rate_CR = 0.0;
  for (int i = Model::cr_heat_var_begin; i < Model::cr_heat_var_end; ++i)
    rate_CR += cell.var[i];
  constexpr double yHe = abundance_ref::yHe;
  rates.Gamma_CR =
      3.4 * phys::eV_to_erg * rate_CR / ((1.0 + 4.0 * yHe) * phys::m_p);

  // ── Net cooling ───────────────────────────────────────────────────────────
  rates.Lambda_net = rates.Lambda_line + rates.Lambda_cnt + rates.Lambda_chem -
                     rates.Gamma_CR
#ifdef ARCHE_XRAY
                     - rates.Gamma_X
#endif
      ;

  return rates;
}

// ---------------------------------------------------------------------------
// make_zero_metal_table()
// ---------------------------------------------------------------------------
inline ZeroMetalTable make_zero_metal_table(
    const std::string& prim_chem_table) {
  ZeroMetalTable tbl;
  zero_metal::net::init_topology(tbl);      // reaction topology from C++ source
  load_pf_tables_h5(tbl, prim_chem_table);  // partition functions from HDF5
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
  return tbl;
}

// ---------------------------------------------------------------------------
// make_metal_grain_table()
// ---------------------------------------------------------------------------
inline MetalGrainTable make_metal_grain_table(
    const std::string& metal_chem_table) {
  MetalGrainTable tbl;
  metal_grain::net::init_topology(tbl);  // reaction topology from C++ source
  load_pf_tables_h5(tbl, metal_chem_table);  // partition functions from HDF5
  return tbl;
}

// ---------------------------------------------------------------------------
// make_metal_minimal_table() — compact 40-species Nakauchi2021_Minimal topology
//   built from the full keep-set via the catalog remap, with the
//   partition-function tables placed into the compact slots by canonical
//   species.  The full PF-loaded table is retained (aux_full_metal) so the
//   strategy-1 rate path and the compact Saha can gather full coefficients /
//   partition functions from it (see compute_base_rates / MinimalMetalSaha).
// ---------------------------------------------------------------------------
inline metal_grain_minimal::MinimalTable make_metal_minimal_table(
    const MetalGrainTable& full) {
  metal_grain_minimal::MinimalTable tbl;
  metal_grain_minimal::build_topology(tbl);  // compact reaction topology
  for (int i = 0; i < metal_grain_minimal::N_sp; ++i)
    tbl.pf_table[i] = full.pf_table[metal_grain::Species::local(
        metal_grain_minimal::Species::canonical(i))];
  tbl.aux_full_metal = std::make_shared<MetalGrainTable>(full);
  return tbl;
}

}  // namespace arche
