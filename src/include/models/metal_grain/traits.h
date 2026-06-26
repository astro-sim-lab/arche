// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// traits.h — compile-time model trait for the metal_grain model
//
// Bundles the network dimensions, solver behaviour flags, the Saha-equilibrium
// policy, and cooling-index accessors.  Kernel templates parameterised on
// `template<class Model>` use these via Model::N_sp, Model::SahaPolicy, etc.
//
//   Nakauchi2021          — metal + grain (metal_grain, 89 species, 1200
//                           reactions)
//   Nakauchi2021_Minimal  — compact metal_grain network (40 species, 103
//                           gas-phase + 10 grain-surface reactions); a
//                           first-class composed model, not a mask over
//                           Nakauchi2021
// ---------------------------------------------------------------------------
#include "kinetics/reaction_index.h"
#include "models/metal_grain/minimal.h"  // metal_grain_minimal (compact model)

namespace arche {

// Forward declarations of the metal_grain Saha-equilibrium policies, defined in
// models/metal_grain/equilibrium.h.  MinimalMetalSaha is introduced with the
// compact 4D Newton in a later phase; the alias below only needs the name.
struct MetalSaha;
struct MinimalMetalSaha;

// ---------------------------------------------------------------------------
// Nakauchi2021  (metal_grain)
// ---------------------------------------------------------------------------
struct Nakauchi2021 {
  using Sp = metal_grain::Sp;            // species-index enum (y[Sp::H] etc.)
  using Species = metal_grain::Species;  // catalog species set (names/masses)
  using SahaPolicy = MetalSaha;          // high-density Saha-equilibrium solver
  static constexpr int N_sp = metal_grain::N_sp;
  static constexpr int N_react = metal_grain::N_react;
  static constexpr RateForm rate_form = RateForm::MetalGrain;

  static constexpr bool has_grain = true;
  static constexpr bool is_compact_prim = false;
  static constexpr bool is_compact_metal = false;
  static constexpr bool has_escape = true;
  static constexpr bool has_uv_shield = true;

  static constexpr int nr_max_iter = 60;
  static constexpr bool catastrophic_detect_always = false;
  static constexpr bool has_in_loop_divergence_guard = true;
  static constexpr int equichem_dy_count = 63;
  static constexpr bool secant_skip_high_density = false;

  static constexpr int cr_heat_var_begin = metal_grain::cr_heat_var_begin;
  static constexpr int cr_heat_var_end = metal_grain::cr_heat_var_end;

  struct cooling {
    // Composable cooling terms carried by this model (see
    // compute_chemistry_cooling): He+/He++ ionization energy, the CR-producer
    // loss correction, and the H2-He dissociation cooling channel.
    static constexpr bool has_he_ion = true;
    static constexpr bool has_he_pp = true;
    static constexpr bool has_cr_loss = true;
    static constexpr bool has_h2_he_cool = true;
    // k_rxn[] indices for the cooling reactions, in this model's numbering
    // (here the full network's shared rxn:: ids).
    struct rate_idx {
      static constexpr int Hp_rec_caseB = rxn::Hp_rec_caseB;
      static constexpr int Hep_rec = rxn::Hep_rec;
      static constexpr int Hm_H_to_H2_e = rxn::Hm_H_to_H2_e;
      static constexpr int H2p_H_to_H2_Hp = rxn::H2p_H_to_H2_Hp;
      static constexpr int three_H = rxn::three_H;
      static constexpr int H2H2_dis = rxn::H2H2_dis;
      static constexpr int H2_He_dis = rxn::H2_He_dis;
    };
    struct Hp_producers {
      static constexpr int H_CR = metal_grain::cooling_idx::Hp_producers::H_CR;
      static constexpr int H_CRph =
          metal_grain::cooling_idx::Hp_producers::H_CRph;
      static constexpr int H2_CR_ch1 =
          metal_grain::cooling_idx::Hp_producers::H2_CR;
      static constexpr int H2_CR_ch2 =
          metal_grain::cooling_idx::Hp_producers::H2_CRph;
    };
    struct Hep_producers {
      static constexpr int He_CR =
          metal_grain::cooling_idx::Hep_producers::He_CR;
      static constexpr int He_CRph =
          metal_grain::cooling_idx::Hep_producers::He_CRph;
    };
    // Grain H2 formation cooling: the 2H+grain channel plus the surface
    // (physisorbed / chemisorbed H) channels.
    struct grain_H2 {
      static constexpr bool has_surface = true;
      static constexpr int two_H = metal_grain::cooling_idx::grain_H2::two_H;
      static constexpr int surface_grgr =
          metal_grain::cooling_idx::grain_H2::surface_grgr;
      static constexpr int surface_grdust =
          metal_grain::cooling_idx::grain_H2::surface_grdust;
      static constexpr int H_dust_3body =
          metal_grain::cooling_idx::grain_H2::H_dust_3body;
      static constexpr int H_gr_3body =
          metal_grain::cooling_idx::grain_H2::H_gr_3body;
    };
  };
};

// ---------------------------------------------------------------------------
// Nakauchi2021_Minimal  (compact metal_grain reduced network)
//
// A first-class composed model carrying only its 40 species and 103 gas-phase
// reactions (68 standard + 8 cosmic-ray + 27 ion-grain charge transfer), plus
// 10 grain-surface freeze-out reactions (metal_grain_minimal), so the linear
// solve is 40x40.  Its forward rates are gathered in compute_base_rates'
// is_compact_metal branch from the full metal_grain kernel run on a scratch
// layout; its Saha uses MinimalMetalSaha; its cooling carries the He+
// recombination term (has_he_ion) and the CR-producer loss / CR heating terms
// (has_cr_loss), but omits the He++ and H2-He terms (those species/reactions
// are not in the keep-set).
// ---------------------------------------------------------------------------
struct Nakauchi2021_Minimal {
  using Sp = metal_grain_minimal::Sp;            // compact 40-species enum
  using Species = metal_grain_minimal::Species;  // 40-species catalog selection
  using SahaPolicy = MinimalMetalSaha;           // compact 4D Newton (Saha)
  static constexpr int N_sp = metal_grain_minimal::N_sp;  // 40
  // 113 k_rxn stride: 103 gas-phase reactions + a 10-wide grain-surface band.
  static constexpr int N_react = metal_grain_minimal::N_react;  // 113
  static constexpr RateForm rate_form = RateForm::MetalGrain;

  static constexpr bool has_grain = true;
  static constexpr bool is_compact_prim = false;
  static constexpr bool is_compact_metal = true;
  static constexpr bool has_escape = true;
  // The grain-coverage feedback (params.J_H2 / J_tot from the physisorbed H2 and
  // chemisorbed H/D abundances) needs species the compact network drops
  // (H2_p, H_c, D_c), so it is disabled; the compact model carries only the
  // O/OH/CO/H2O/C ice-mantle freeze-out, not H2 physisorption shielding.
  static constexpr bool has_uv_shield = false;

  static constexpr int nr_max_iter = 60;
  static constexpr bool catastrophic_detect_always = false;
  static constexpr bool has_in_loop_divergence_guard = true;
  // The compact 4D Newton (MinimalMetalSaha) fills the 31 gas-phase species
  // (indices 0..30 = H..Mg+); the grain charge states and ice-mantle species
  // (31..39) are left unchanged, so the Saha dy-tracking spans 0..30.
  static constexpr int equichem_dy_count = 31;
  static constexpr bool secant_skip_high_density = false;

  // CR heating sums the four direct cosmic-ray ionization channels at the head
  // of the compact CR block (slots 68..71).
  static constexpr int cr_heat_var_begin =
      metal_grain_minimal::cr_heat_var_begin;  // 68
  static constexpr int cr_heat_var_end =
      metal_grain_minimal::cr_heat_var_end;  // 72

  struct cooling {
    // He+ recombination cooling (has_he_ion) and the CR-producer loss / CR
    // heating terms (has_cr_loss) are carried; He++ is absent (no He++ species)
    // and the H2-He dissociation channel is not in the keep-set.
    static constexpr bool has_he_ion = true;
    static constexpr bool has_he_pp = false;
    static constexpr bool has_cr_loss = true;
    static constexpr bool has_h2_he_cool = false;
    // Cooling reactions in the compact numbering (slot = keep-set position in
    // metal_grain_minimal::net::kMetalMinimalKeep).
    struct rate_idx {
      static constexpr int Hp_rec_caseB = 0;    // num 2
      static constexpr int Hep_rec = 1;         // num 4
      static constexpr int Hm_H_to_H2_e = 3;    // num 8
      static constexpr int H2p_H_to_H2_Hp = 5;  // num 10
      static constexpr int three_H = 7;         // num 19
      static constexpr int H2H2_dis = 8;        // num 21
      static constexpr int H2_He_dis = -1;      // absent (has_h2_he_cool=false)
    };
    // CR producers of H+ / He+ in the compact CR numbering (slots 68..75).
    struct Hp_producers {
      static constexpr int H_CR =
          metal_grain_minimal::cooling_idx::Hp_producers::H_CR;
      static constexpr int H_CRph =
          metal_grain_minimal::cooling_idx::Hp_producers::H_CRph;
      static constexpr int H2_CR_ch1 =
          metal_grain_minimal::cooling_idx::Hp_producers::H2_CR;
      // The H2 CR-induced-photon channel (full id 551) is not in the keep-set.
      // The cooling consumer currently reads both H2 channels unconditionally,
      // so the absent second channel must be guarded when the compact cooling
      // path is wired (cooling phase); -1 marks it absent.
      static constexpr int H2_CR_ch2 = -1;
    };
    struct Hep_producers {
      static constexpr int He_CR =
          metal_grain_minimal::cooling_idx::Hep_producers::He_CR;
      static constexpr int He_CRph =
          metal_grain_minimal::cooling_idx::Hep_producers::He_CRph;
    };
    // Grain H2 formation cooling.  Only the grain-catalysed 2H+grain channel
    // (compact slot 9) survives; the surface H2 formation channels need the
    // physisorbed / chemisorbed H species (H_p / H_c) the network drops.
    struct grain_H2 {
      static constexpr bool has_surface = false;
      static constexpr int two_H = metal_grain_minimal::special::rxn_2H_grain - 1;
      static constexpr int surface_grgr = -1;
      static constexpr int surface_grdust = -1;
      static constexpr int H_dust_3body = -1;
      static constexpr int H_gr_3body = -1;
    };
  };
};

}  // namespace arche
