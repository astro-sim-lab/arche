// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// minimal.h — compact metal_grain Minimal network (40 species, 104 gas-phase
// reactions + 10 grain-surface freeze-out reactions).
//
// A first-class composed model: a 40-species selection from the canonical
// catalog carrying the reduced chemical network of a low-metallicity
// star-forming cloud (H/He/D core, a compact C/O molecular network, Li and Mg,
// four grain charge states, and five ice-mantle surface species).  The
// gas-phase reactions split into three solver blocks that mirror the full
// network's partition: standard two-body reactions with detailed-balance
// reverse, the cosmic-ray channels, and the ion-grain charge-transfer reactions
// (both first-order).  Grain-surface freeze-out (adsorption and thermal
// desorption) is carried in the separate grain-reaction table.
//
// Its topology is the full metal_grain network's keep-set relabelled into the
// compact local index space, so the linear solve is 42x42.  The compact
// reaction / grain / Saha tables are derived from the full network
// (models/metal_grain/reactions.h) by selecting the keep-set, remapping species
// through the canonical<->local maps, and renumbering 1..N.  Deriving rather
// than re-typing guarantees the kept reactions carry exactly the full network's
// stoichiometry, Cmass, and delE.
//
// He and the carbon ion / neutral chemistry are both retained: the network is
// the union of the two reduced sub-networks (with and without an ionization
// source), so a single model spans the neutral-CH and the C+ regimes.
// ---------------------------------------------------------------------------
#include <array>

#include "core/species_catalog.h"
#include "kinetics/reaction_index.h"
#include "kinetics/topology.h"
#include "models/metal_grain/reactions.h"

namespace arche {
namespace metal_grain_minimal {

// 40-species selection (canonical order).  Core H/He/D (12), compact carbon
// (C, CH, CH2, C+), neutral oxygen (O, O2, OH, CO, H2O, HCO), O/C ions (O+,
// OH+, H2O+, HCO+, H3O+), Li/Mg, four grain charge states (no gr2+), and five
// ice-mantle species (O, OH, CO, H2O, C on grain).
using Species =
    SpeciesSet<SpId::H, SpId::H2, SpId::e, SpId::Hp, SpId::H2p, SpId::H3p,
               SpId::Hm, SpId::He, SpId::Hep, SpId::D, SpId::HD, SpId::Dp,
               SpId::C, SpId::CH, SpId::CH2, SpId::Cp, SpId::O, SpId::O2,
               SpId::OH, SpId::CO, SpId::H2O, SpId::HCO, SpId::Op, SpId::OHp,
               SpId::H2Op, SpId::HCOp, SpId::H3Op, SpId::Li, SpId::Lip,
               SpId::Mg, SpId::Mgp, SpId::Gr, SpId::Grp, SpId::Grm, SpId::Gr2m,
               SpId::O_p, SpId::OH_p, SpId::CO_p, SpId::H2O_p, SpId::C_p>;

constexpr int N_sp = Species::N;  // 40
constexpr int N_gas = 103;        // 68 standard + 8 CR + 27 grain-charge
constexpr int N_grain = 10;       // 5 adsorption + 5 thermal desorption
// k_rxn stride.  Like the full network (metal_grain::N_react = 1200 >> 821
// gas-phase reactions), the stride exceeds the gas-phase reaction count so the
// grain-surface forward rates occupy a band above the gas block: gas forward
// rates fill k_rxn[0, N_gas), the grain freeze-out forward rates fill
// k_rxn[N_gas, N_gas + N_grain), and the detailed-balance reverse rates start
// at k_rxn[N_react] = k_rxn[N_gas + N_grain].  Sizing the stride to N_gas alone
// would push the grain band into the reverse region and the loop guard
// (num > N_react) would silently drop every grain reaction.
constexpr int N_react = N_gas + N_grain;  // 113
constexpr int N_saha = 23;  // high-density Saha sub-table (provisional)

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
  C = Species::local(SpId::C),
  CH = Species::local(SpId::CH),
  CH2 = Species::local(SpId::CH2),
  Cp = Species::local(SpId::Cp),
  O = Species::local(SpId::O),
  O2 = Species::local(SpId::O2),
  OH = Species::local(SpId::OH),
  CO = Species::local(SpId::CO),
  H2O = Species::local(SpId::H2O),
  HCO = Species::local(SpId::HCO),
  Op = Species::local(SpId::Op),
  OHp = Species::local(SpId::OHp),
  H2Op = Species::local(SpId::H2Op),
  HCOp = Species::local(SpId::HCOp),
  H3Op = Species::local(SpId::H3Op),
  Li = Species::local(SpId::Li),
  Lip = Species::local(SpId::Lip),
  Mg = Species::local(SpId::Mg),
  Mgp = Species::local(SpId::Mgp),
  Gr = Species::local(SpId::Gr),
  Grp = Species::local(SpId::Grp),
  Grm = Species::local(SpId::Grm),
  Gr2m = Species::local(SpId::Gr2m),
  O_p = Species::local(SpId::O_p),
  OH_p = Species::local(SpId::OH_p),
  CO_p = Species::local(SpId::CO_p),
  H2O_p = Species::local(SpId::H2O_p),
  C_p = Species::local(SpId::C_p),
};

static_assert(N_sp == 40, "metal Minimal network must carry 40 species");
static_assert(
    H == 0 && He == 7 && C == 12 && O == 16 && Gr == 31 && Gr2m == 34 &&
        C_p == 39,
    "metal Minimal Sp indices must follow the catalog selection order");

// ── Loop bounds.  compute_rates runs the standard forward+reverse loop over
//    [0, n_std), the cosmic-ray first-order loop over [n_std,
//    n_std+n_cr_react), and the grain-charge first-order loop over the
//    remaining slots.  The grain freeze-out reactions live in the separate
//    grain-reaction table. ──
namespace loop {
constexpr int n_std = 68;
constexpr int n_cr_react = 8;
constexpr int n_charge_react = 27;
constexpr int n_grain_react = 10;
}  // namespace loop

static_assert(loop::n_std + loop::n_cr_react + loop::n_charge_react == N_gas,
              "compact gas-phase loop blocks must tile the gas reaction block");
static_assert(loop::n_grain_react == N_grain,
              "grain loop count must match the grain-surface band width");

// ── Special-form gas-phase reactions (compact renumbering).  The grain-
//    catalysed three-body reactions 2H+grain->H2 and 2D+grain->HD use a
//    non-product-form rate law (rate = k * y[H or D] * nH); compute_rates
//    selects these by compact num.  The full network's 4-body H reaction
//    (full num 273) is not in the keep-set, so its sentinel never matches. ──
namespace special {
constexpr int rxn_2H_grain = 10;  // full num 23 -> compact slot 9 (k_rxn[9])
constexpr int rxn_2D_grain = 26;  // full num 144 -> compact slot 25 (k_rxn[25])
constexpr int rxn_4body_H = -1;   // absent from the Minimal keep-set
}  // namespace special

// var[] range summed for CR heating in chemistry.h: the four direct cosmic-ray
// ionization channels at the head of the CR block (compact slots 68..71 = full
// ids 544, 545, 548, 549).
constexpr int cr_heat_var_begin = loop::n_std;    // 68
constexpr int cr_heat_var_end = loop::n_std + 4;  // 72 (exclusive)

// Hot-spot k_rxn[] indices for the cooling/heating block, in the compact
// numbering.  Channels absent from the keep-set (e.g. the H2 CR-induced-photon
// ionization) have no compact slot and are simply not summed.
namespace cooling_idx {
namespace Hp_producers {
constexpr int H_CR = 68;    // id 544  H  + CR   -> e + H+
constexpr int H2_CR = 70;   // id 548  H2 + CR   -> H + e + H+
constexpr int H_CRph = 74;  // id 677  H  + CRph -> e + H+
}  // namespace Hp_producers
namespace Hep_producers {
constexpr int He_CR = 69;    // id 545  He + CR   -> e + He+
constexpr int He_CRph = 75;  // id 678  He + CRph -> e + He+
}  // namespace Hep_producers
}  // namespace cooling_idx

using MinimalTable = ReactionTable<N_sp, N_react>;

// Map a full-network table index (real species 0..N_sp_full-1, or a sentinel
// VACANT/PHOTON/CR at N_sp_full + {0,1,2}) into the compact index space.  A
// real species that is not part of the Minimal selection maps to -1 (a kept
// reaction must never reference one — the builder asserts this).
inline constexpr int remap_index(int full_idx) {
  constexpr int Nf = metal_grain::N_sp;
  if (full_idx < Nf) {
    SpId canon = metal_grain::Species::canonical(full_idx);
    return Species::local(canon);
  }
  // Sentinel: VACANT(Nf)->N_sp, PHOTON(Nf+1)->N_sp+1, CR(Nf+2)->N_sp+2.
  return (full_idx - Nf) + N_sp;
}

// Build the compact topology by relabelling the full network's keep-set.
// Returns true on success; sets *bad_species to the offending full index if a
// kept reaction references a species outside the Minimal selection.
inline bool build_topology(MinimalTable& tbl, int* bad_species = nullptr) {
  auto remap_slots = [](const int* src, int n, int* dst, int* bad) -> bool {
    bool ok = true;
    for (int i = 0; i < n; ++i) {
      int m = remap_index(src[i]);
      if (m < 0) {
        ok = false;
        if (bad) *bad = src[i];
      }
      dst[i] = m;
    }
    return ok;
  };

  // Build one compact gas-phase reaction row from a full-network keep id.
  auto build_row = [&](int keep_num, int slot, Reaction& d) -> bool {
    const Reaction* src = nullptr;
    for (const auto& r : metal_grain::net::kReactions)
      if (r.num == keep_num) {
        src = &r;
        break;
      }
    if (!src) return false;  // keep-set id absent from the full network
    d.num = slot + 1;
    if (!remap_slots(src->reactants.data(), 3, d.reactants.data(), bad_species))
      return false;
    if (!remap_slots(src->products.data(), 3, d.products.data(), bad_species))
      return false;
    d.n_reactants = src->n_reactants;
    d.n_products = src->n_products;
    d.Cmass = src->Cmass;
    d.delE = src->delE;
    return true;
  };

  // ── gas-phase reactions, in solver-block order: standard (detailed-balance
  //    reverse), then cosmic-ray, then grain charge transfer.  Renumbered
  //    1..N_react in slot order. ──
  int slot = 0;
  for (int keep_num : metal_grain::net::kMetalMinimalKeep) {
    Reaction d;
    if (!build_row(keep_num, slot, d)) return false;
    tbl.reactions[slot++] = d;
  }
  for (int keep_num : metal_grain::net::kMetalMinimalCRKeep) {
    Reaction d;
    if (!build_row(keep_num, slot, d)) return false;
    tbl.reactions[slot++] = d;
  }
  for (int keep_num : metal_grain::net::kMetalMinimalChargeKeep) {
    Reaction d;
    if (!build_row(keep_num, slot, d)) return false;
    tbl.reactions[slot++] = d;
  }
  tbl.n_loaded = slot;

  // ── grain-surface freeze-out: adsorption + thermal desorption, renumbered
  //    1..N_grain. ──
  int gslot = 0;
  for (int keep_num : metal_grain::net::kMetalMinimalGrainKeep) {
    const GrainReaction* src = nullptr;
    for (const auto& g : metal_grain::net::kGrain)
      if (g.num == keep_num) {
        src = &g;
        break;
      }
    if (!src) return false;
    GrainReaction d;
    d.num = gslot + 1;
    if (!remap_slots(src->reactants.data(), 2, d.reactants.data(), bad_species))
      return false;
    if (!remap_slots(src->products.data(), 2, d.products.data(), bad_species))
      return false;
    d.n_reactants = src->n_reactants;
    tbl.grain_reactions[gslot++] = d;
  }
  tbl.n_grain = gslot;

  // ── Saha sub-table: keep-set, renumbered 1..N_saha. ──
  int ss = 0;
  for (int keep_num : metal_grain::net::kMetalMinimalSahaKeep) {
    const SahaReaction* src = nullptr;
    for (const auto& s : metal_grain::net::kSaha)
      if (s.num == keep_num) {
        src = &s;
        break;
      }
    if (!src) return false;
    SahaReaction d;
    d.num = ss + 1;
    if (!remap_slots(src->reactants.data(), 2, d.reactants.data(), bad_species))
      return false;
    if (!remap_slots(src->products.data(), 2, d.products.data(), bad_species))
      return false;
    d.n_reactants = src->n_reactants;
    d.n_products = src->n_products;
    d.Cmass = src->Cmass;
    d.delE = src->delE;
    tbl.saha[ss++] = d;
  }
  tbl.n_saha = ss;

  return true;
}

}  // namespace metal_grain_minimal
}  // namespace arche
