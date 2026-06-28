// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// minimal.h — compact Nakauchi2019_Minimal network (15 species, 33 reactions).
//
// A first-class composed model: a 15-species selection from the canonical
// catalog, carrying the minimal reaction set plus the ion-processing reactions
// (ids 4, 15, 17, 24, 38) and the 9 cosmic-ray ionization channels (ids
// 131..139) that let it track the full network's H / D / Li chemistry both
// thermally and under cosmic-ray ionization.  He+ is the one species the
// compact set tracks only under cosmic rays: it is produced solely by the CR
// channels (ids 136/139) and removed by recombination (id 4) and the H2 + He+
// charge transfer (id 24); the collisional He ionization (full id 3) and the
// H/He charge-transfer reactions (full ids 28, 41) are not carried, so without
// cosmic rays He+ stays ~0 instead of forming collisionally at high T.  Its
// topology is the full primordial network's keep-set relabelled into the
// compact local index space, so the linear solve is 15x15 and the reaction loop
// runs 24 standard + 9 CR entries.
//
// The compact reaction / Saha tables are derived from the full network
// (models/primordial/reactions.h) by selecting the keep-set, remapping species
// through the canonical<->local maps, and renumbering 1..N.  Deriving rather
// than re-typing guarantees the kept reactions carry exactly the full network's
// stoichiometry, Cmass, and delE.
// ---------------------------------------------------------------------------
#include <array>

#include "core/species_catalog.h"
#include "kinetics/reaction_index.h"
#include "kinetics/topology.h"
#include "models/primordial/reactions.h"

namespace arche {
namespace zero_metal_minimal {

// 15-species selection (canonical order); He retained as neutral background,
// He+ added so the He recombination (id 4), the H2 + He+ charge transfer
// (id 24), and the cosmic-ray He ionization (ids 136/139) are represented.
// No collisional He ionization (full id 3) is carried, so without cosmic rays
// He+ stays ~0.  Removed vs the full set: He++, HeH+, HD+, D-, Li-, LiH+, Li++,
// Li+++.
using Species = SpeciesSet<SpId::H, SpId::H2, SpId::e, SpId::Hp, SpId::H2p,
                           SpId::H3p, SpId::Hm, SpId::He, SpId::Hep, SpId::D,
                           SpId::HD, SpId::Dp, SpId::Li, SpId::LiH, SpId::Lip>;

constexpr int N_sp = Species::N;  // 15
constexpr int N_react = 33;       // 24 standard + 9 cosmic-ray channels
constexpr int N_saha = 9;         // high-density Saha sub-table

// Compact species indices, derived from the catalog selection.
enum Sp : int {
  H = Species::local(SpId::H),
  H2 = Species::local(SpId::H2),
  e = Species::local(SpId::e),
  Hp = Species::local(SpId::Hp),
  H2p = Species::local(SpId::H2p),
  H3p = Species::local(SpId::H3p),
  Hm = Species::local(SpId::Hm),
  He = Species::local(SpId::He),
  Hep = Species::local(SpId::Hep),
  D = Species::local(SpId::D),
  HD = Species::local(SpId::HD),
  Dp = Species::local(SpId::Dp),
  Li = Species::local(SpId::Li),
  LiH = Species::local(SpId::LiH),
  Lip = Species::local(SpId::Lip),
};

static_assert(N_sp == 15, "minimal network must carry 15 species");
static_assert(H == 0 && He == 7 && Hep == 8 && Lip == 14,
              "minimal Sp indices must follow the catalog selection order");

// ── Loop bounds: 24 standard reactions (slots 0..23) + 9 CR reactions
//    (slots 24..32).  compute_rates runs the standard forward+reverse loop over
//    [0, n_std) and the CR first-order loop over [n_std, n_std+n_cr_react). ──
namespace loop {
constexpr int n_std = 24;
constexpr int n_cr_react = 9;
}  // namespace loop

// var[] range summed for CR heating in chemistry.h: the 6 direct ionization
// channels (compact slots 24..29 = full ids 131..136).
constexpr int cr_heat_var_begin = 24;
constexpr int cr_heat_var_end = 30;  // exclusive

// Hot-spot k_rxn[] indices for the cooling/heating block, in the compact
// numbering (CR keep-set order: 131->24, 132->25, ..., 139->32).
namespace cooling_idx {
namespace Hp_producers {
constexpr int H_CR = 24;        // id 131  H  + CR   -> e- + H+
constexpr int H_CRph = 30;      // id 137  H  + CRph -> e- + H+
constexpr int H2_CR_HpH = 25;   // id 132  H2 + CR   -> H + e- + H+
constexpr int H2_CR_HpHm = 28;  // id 135  H2 + CR   -> H+ + H-
}  // namespace Hp_producers
namespace Hep_producers {
constexpr int He_CR = 29;    // id 136  He + CR   -> e- + He+
constexpr int He_CRph = 32;  // id 139  He + CRph -> e- + He+
}  // namespace Hep_producers
}  // namespace cooling_idx

using MinimalTable = ReactionTable<N_sp, N_react>;

// High-density Saha sub-table (canonical reaction ids), 9 rows.  Re-keyed to
// compact ids 1..9 in selection order when the table is built; equichem reads
// the compact ids.  (D+ via id 51 and HD via id 54 are needed even though they
// are not in the kinetic keep-set.)
inline constexpr std::array<int, N_saha> kSahaKeep = {{
    2,    // e- + H+  -> H  + γ
    7,    // H  + e-  -> H- + γ
    8,    // H  + H-  -> H2 + e-
    9,    // H  + H+  -> H2+ + γ
    26,   // H2 + H2+ -> H  + H3+
    51,   // e- + D+  -> D  + γ
    54,   // H  + D   -> HD + γ
    101,  // e- + Li+ -> Li + γ
    118,  // H  + Li  -> LiH + γ
}};

// Map a full-network table index (real species 0..N_sp_full-1, or a sentinel
// VACANT/PHOTON/CR at N_sp_full + {0,1,2}) into the compact index space.  A
// real species that is not part of the minimal selection maps to -1 (a kept
// reaction must never reference one — the builder asserts this).
inline constexpr int remap_index(int full_idx) {
  constexpr int Nf = zero_metal::N_sp;
  if (full_idx < Nf) {
    SpId canon = zero_metal::Species::canonical(full_idx);
    return Species::local(canon);
  }
  // Sentinel: VACANT(Nf)->N_sp, PHOTON(Nf+1)->N_sp+1, CR(Nf+2)->N_sp+2.
  return (full_idx - Nf) + N_sp;
}

// Build the compact topology by relabelling the full network's keep-set.
// Returns true on success; sets *bad_species to the offending full index if a
// kept reaction references a species outside the minimal selection.
inline bool build_topology(MinimalTable& tbl, int* bad_species = nullptr) {
  auto remap_slots = [](const std::array<int, 3>& src, std::array<int, 3>& dst,
                        int* bad) -> bool {
    bool ok = true;
    for (int i = 0; i < 3; ++i) {
      int m = remap_index(src[i]);
      if (m < 0) {
        ok = false;
        if (bad) *bad = src[i];
      }
      dst[i] = m;
    }
    return ok;
  };

  // Build one compact reaction row from a full-network keep id (lookup +
  // remap).
  auto build_row = [&](int keep_num, int slot, Reaction& d) -> bool {
    const Reaction* src = nullptr;
    for (const auto& r : zero_metal::net::kReactions)
      if (r.num == keep_num) {
        src = &r;
        break;
      }
    if (!src) return false;  // keep-set id absent from the full network
    d.num = slot + 1;
    if (!remap_slots(src->reactants, d.reactants, bad_species)) return false;
    if (!remap_slots(src->products, d.products, bad_species)) return false;
    d.n_reactants = src->n_reactants;
    d.n_products = src->n_products;
    d.Cmass = src->Cmass;
    d.delE = src->delE;
    return true;
  };

  // ── standard reactions: keep-set, renumbered 1..n_std (slots 0..23) ──
  int slot = 0;
  for (int keep_num : zero_metal::net::kMinimalKeep) {
    Reaction d;
    if (!build_row(keep_num, slot, d)) return false;
    tbl.reactions[slot++] = d;
  }
  // ── CR reactions: keep-set, continuing the numbering (slots 24..32) ──
  for (int keep_num : zero_metal::net::kMinimalCRKeep) {
    Reaction d;
    if (!build_row(keep_num, slot, d)) return false;
    tbl.reactions[slot++] = d;
  }
  tbl.n_loaded = slot;

  // ── Saha sub-table: keep-set, renumbered 1..N_saha ──
  int ss = 0;
  for (int keep_num : kSahaKeep) {
    const SahaReaction* src = nullptr;
    for (const auto& s : zero_metal::net::kSaha)
      if (s.num == keep_num) {
        src = &s;
        break;
      }
    if (!src) return false;
    SahaReaction d;
    d.num = ss + 1;
    {
      std::array<int, 3> r3{src->reactants[0], src->reactants[1],
                            zero_metal::N_sp};
      std::array<int, 3> p3{src->products[0], src->products[1],
                            zero_metal::N_sp};
      std::array<int, 3> mr{}, mp{};
      if (!remap_slots(r3, mr, bad_species)) return false;
      if (!remap_slots(p3, mp, bad_species)) return false;
      d.reactants = {mr[0], mr[1]};
      d.products = {mp[0], mp[1]};
    }
    d.n_reactants = src->n_reactants;
    d.n_products = src->n_products;
    d.Cmass = src->Cmass;
    d.delE = src->delE;
    tbl.saha[ss++] = d;
  }
  tbl.n_saha = ss;

  tbl.n_grain = 0;
  return true;
}

}  // namespace zero_metal_minimal
}  // namespace arche
