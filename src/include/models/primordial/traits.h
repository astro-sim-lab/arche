// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// traits.h — compile-time model traits for the primordial (zero_metal) models
//
// Each struct bundles the network dimensions, solver behaviour flags, the
// Saha-equilibrium policy, and cooling-index accessors for one primordial
// model.  Kernel templates parameterised on `template<class Model>` use these
// via Model::N_sp, Model::SahaPolicy, etc.
//
//   Nakauchi2019          — primordial (zero_metal, 23 species, 140 reactions)
//   Nakauchi2019_Minimal  — compact primordial network (15 species, 33
//                           reactions); a first-class composed model, not a
//                           mask over Nakauchi2019
// ---------------------------------------------------------------------------
#include "kinetics/reaction_index.h"
#include "models/primordial/minimal.h"  // zero_metal_minimal (compact model)

namespace arche {

// Forward declarations of the primordial Saha-equilibrium policies, defined in
// models/primordial/equilibrium.h.
struct FullPrimSaha;
struct MinimalSaha;

// ---------------------------------------------------------------------------
// Nakauchi2019  (primordial / zero_metal)
// ---------------------------------------------------------------------------
struct Nakauchi2019 {
  using Sp = zero_metal::Sp;            // species-index enum (y[Sp::H] etc.)
  using Species = zero_metal::Species;  // catalog species set (names/masses)
  using SahaPolicy = FullPrimSaha;      // high-density Saha-equilibrium solver
  static constexpr int N_sp = zero_metal::N_sp;
  static constexpr int N_react = zero_metal::N_react;
  static constexpr RateForm rate_form = RateForm::Primordial;

  static constexpr bool has_grain = false;
  static constexpr bool is_compact_prim = false;
  static constexpr bool is_compact_metal = false;
  static constexpr bool has_escape = false;
  static constexpr bool has_uv_shield = false;

  static constexpr int nr_max_iter = 30;
  static constexpr bool catastrophic_detect_always = true;
  static constexpr bool has_in_loop_divergence_guard = false;
  static constexpr int equichem_dy_count = 23;
  static constexpr bool secant_skip_high_density = true;

  static constexpr int cr_heat_var_begin = zero_metal::cr_heat_var_begin;
  static constexpr int cr_heat_var_end = zero_metal::cr_heat_var_end;

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
      static constexpr int H_CR = zero_metal::cooling_idx::Hp_producers::H_CR;
      static constexpr int H_CRph =
          zero_metal::cooling_idx::Hp_producers::H_CRph;
      static constexpr int H2_CR_ch1 =
          zero_metal::cooling_idx::Hp_producers::H2_CR_HpH;
      static constexpr int H2_CR_ch2 =
          zero_metal::cooling_idx::Hp_producers::H2_CR_HpHm;
    };
    struct Hep_producers {
      static constexpr int He_CR =
          zero_metal::cooling_idx::Hep_producers::He_CR;
      static constexpr int He_CRph =
          zero_metal::cooling_idx::Hep_producers::He_CRph;
    };
  };
};

// ---------------------------------------------------------------------------
// Nakauchi2019_Minimal  (compact primordial reduced network)
//
// The true-compacted Nakauchi2019_Minimal: a first-class composed model
// carrying only its 15 species and 33 reactions (24 standard + 9 cosmic-ray;
// zero_metal_minimal), so the linear solve is 15x15 and the standard reaction
// loop runs 24 entries.  Its forward
// rates are assembled in compute_base_rates' is_compact_prim branch from the
// shared ratelaw::prim laws; its Saha uses MinimalSaha; its cooling carries the
// He+ recombination term (has_he_ion) and the CR-producer loss / CR heating
// terms (has_cr_loss), but omits the He++ and H2-He terms (those
// species/reactions are not present).
// ---------------------------------------------------------------------------
struct Nakauchi2019_Minimal {
  using Sp = zero_metal_minimal::Sp;  // compact 15-species enum
  using Species = zero_metal_minimal::Species;
  using SahaPolicy = MinimalSaha;
  static constexpr int N_sp = zero_metal_minimal::N_sp;        // 15
  static constexpr int N_react = zero_metal_minimal::N_react;  // 33
  static constexpr RateForm rate_form = RateForm::Primordial;

  static constexpr bool has_grain = false;
  static constexpr bool is_compact_prim = true;
  static constexpr bool is_compact_metal = false;
  static constexpr bool has_escape = false;
  static constexpr bool has_uv_shield = false;

  static constexpr int nr_max_iter = 30;
  static constexpr bool catastrophic_detect_always = true;
  static constexpr bool has_in_loop_divergence_guard = false;
  static constexpr int equichem_dy_count = zero_metal_minimal::N_sp;  // 15
  static constexpr bool secant_skip_high_density = true;

  // CR heating sums the 6 direct ionization channels (compact slots 24..29).
  static constexpr int cr_heat_var_begin =
      zero_metal_minimal::cr_heat_var_begin;
  static constexpr int cr_heat_var_end = zero_metal_minimal::cr_heat_var_end;

  struct cooling {
    // He+ recombination cooling (has_he_ion) and the CR-producer loss / CR
    // heating terms (has_cr_loss) are carried; He++ and the H2-He dissociation
    // channel are not.
    static constexpr bool has_he_ion = true;
    static constexpr bool has_he_pp = false;
    static constexpr bool has_cr_loss = true;
    static constexpr bool has_h2_he_cool = false;
    // Cooling reactions in the compact numbering (slot = keep-set position).
    struct rate_idx {
      static constexpr int Hp_rec_caseB = 0;    // num 2
      static constexpr int Hep_rec = 19;        // num 4 (appended)
      static constexpr int Hm_H_to_H2_e = 2;    // num 8
      static constexpr int H2p_H_to_H2_Hp = 4;  // num 10
      static constexpr int three_H = 5;         // num 19
      static constexpr int H2H2_dis = 6;        // num 21
      static constexpr int H2_He_dis = -1;      // absent
    };
    // CR producers of H+ / He+ in the compact CR numbering (slots 24..32).
    struct Hp_producers {
      static constexpr int H_CR =
          zero_metal_minimal::cooling_idx::Hp_producers::H_CR;
      static constexpr int H_CRph =
          zero_metal_minimal::cooling_idx::Hp_producers::H_CRph;
      static constexpr int H2_CR_ch1 =
          zero_metal_minimal::cooling_idx::Hp_producers::H2_CR_HpH;
      static constexpr int H2_CR_ch2 =
          zero_metal_minimal::cooling_idx::Hp_producers::H2_CR_HpHm;
    };
    struct Hep_producers {
      static constexpr int He_CR =
          zero_metal_minimal::cooling_idx::Hep_producers::He_CR;
      static constexpr int He_CRph =
          zero_metal_minimal::cooling_idx::Hep_producers::He_CRph;
    };
  };
};

}  // namespace arche
