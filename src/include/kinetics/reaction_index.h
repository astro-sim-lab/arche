// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// reaction_index.h — internal reaction-rate bookkeeping (NOT a public header).
//
// These constants describe how the rate kernel lays out the forward-rate array
// k_rxn[] and which slots the cooling/heating expressions read.  They are
// tightly coupled to the reaction numbering stored in
// data/{primordial,metal_grain}.h5 (reactions/rean[]) and to the cooling/solver
// code, so they are deliberately kept out of the public species_index.h: a
// consumer's contract must not move when this internal layout changes.
//
//   * rxn::                   — cross-network shared reaction indices
//   * <model>::slot::         — forward-rate slot allocation in k_rxn[]
//   * <model>::loop::         — detailed-balance / first-order loop bounds
//   * <model>::cooling_idx::  — k_rxn[] hot-spots used by cooling/heating
//   * <model>::special::      — reactions with a non-trivial forward/Jacobian
//   * <model>::cr_heat_var_*  — var[] range summed for CR heating
//   * IDX_NONE / IDX_PHOTON   — reaction-table sentinel encoding
//   * RateForm                — expression-form selector in rate routines
//
// The trailing static_asserts encode the invariants between these constants and
// the HDF5 row order; a violation means the data file layout and the C++
// constants have drifted out of sync (update both together).
// ---------------------------------------------------------------------------
#include "core/species_index.h"

namespace arche {

// ---------------------------------------------------------------------------
// Cross-network shared reaction indices (rxn:: namespace)
//
// k_rxn[0..99] in both zero_metal::N_react=140 and metal_grain::N_react=1200
// stores the SAME 100 reactions (GA08 H/He/D network).  These named constants
// are the canonical source for any code (rate assignment in the model
// reaction headers, cooling expressions in solve/chemreact.h, internal alias
// lines like k_rxn[63]=k_rxn[14]) that refers to one of these shared slots.
//
// IMPORTANT: These names are coupled to the reaction numbering
// stored in data/{primordial,metal_grain}.h5 (reactions/rean[]).  Both files
// must list the same chemistry in the same row positions for indices to
// agree across networks.  Changes here require synchronised changes to the
// HDF5 data files.
// ---------------------------------------------------------------------------
namespace rxn {
// ── Ionisation / recombination ──
constexpr int H_e_to_Hp_2e = 0;  // k_rxn[0]: H + e -> H+ + 2e
constexpr int Hp_rec_caseB = 1;  // k_rxn[1]: H+ + e -> H + ph. (case B)
constexpr int Hep_rec = 3;  // k_rxn[3]: He+ + e -> He + ph. (cool/cool fit)

// ── H- formation / destruction ──
constexpr int H_e_to_Hm = 6;      // k_rxn[6]: H + e -> H- + ph.
constexpr int Hm_H_to_H2_e = 7;   // k_rxn[7]: H- + H -> H2 + e  (Kreckel+2010)
constexpr int Hm_e_to_H_2e = 13;  // k_rxn[13]: H- + e -> H + 2e
constexpr int Hm_Hp_to_2H = 14;   // k_rxn[14]: H- + H+ -> 2H
constexpr int Hm_H_to_2H_e = 37;  // k_rxn[37]: H- + H -> 2H + e

// ── H2+ ──
constexpr int H2p_H_to_H2_Hp = 9;  // k_rxn[9]: H2+ + H -> H2 + H+

// ── H2 formation (3-body) / dissociation (forward storage) ──
// Forward index is the dissociation direction unless noted.
// The formation contribution uses k_rxn[i + N_react] (detailed-balance
// reverse).
constexpr int three_H = 18;    // k_rxn[18]: 3H -> H2 + H (Forrey+2013)
constexpr int H2H2_dis = 20;   // k_rxn[20]: 2H2 -> 2H + H2 (LTE blend)
constexpr int H2_He_dis = 36;  // k_rxn[36]: H2 + He -> 2H + He

// ── He+ ──
constexpr int Hep_H2_to_H2p_He = 29;  // k_rxn[29]: He+ + H2 -> H2+ + He
}  // namespace rxn

// ---------------------------------------------------------------------------
// zero_metal rate-kernel bookkeeping
// ---------------------------------------------------------------------------
namespace zero_metal {

// var[] range summed for CR heating in chemistry.h (zero_metal CR block:
// reactions 131..136 -> var[130..135]).  Mirrors metal_grain::cr_heat_var_*.
constexpr int cr_heat_var_begin = 130;
constexpr int cr_heat_var_end = 136;  // exclusive

// ---------------------------------------------------------------------------
// Network layout constants (forward-rate slot allocation in
// k_rxn[0..N_react-1]). Index convention: 0-based k_rxn[] index N == 1-based
// reaction id N+1. The generic compute_base_rates<N_sp,N_react> /
// compute_rates<N_sp,N_react> templates in models/rate_kernel.h are
// effectively zero_metal-shaped; these constants name that shape so any new
// model can re-anchor or specialise.
// ---------------------------------------------------------------------------
namespace slot {
// Standard reactions (1..130, GA08 H/He/D/Li) -> k_rxn[0..129]
constexpr int std_begin = 0;
constexpr int std_end = 130;  // exclusive
// CR / CR-photo (131..139) -> k_rxn[130..138]
constexpr int cr_begin = 130;
constexpr int cr_end = 139;  // exclusive (k_rxn[139] unused; N_react=140)
}  // namespace slot

namespace loop {
constexpr int n_std = 130;     // forward+reverse loop bound
constexpr int n_cr_react = 9;  // CR first-order loop bound (139 - 130)
}  // namespace loop

// Hot-spot k_rxn[] indices referenced by the cooling/heating block in
// solve/chemreact.h, organised by use in the cooling expression:
//   dyHp  = recomb - Σ(Hp_producers)
//   dyHep = recomb - Σ(Hep_producers)
// Recombination indices live in arche::rxn:: (shared with metal_grain).
namespace cooling_idx {
// k_rxn[] channels producing H+ from neutral H/H2 (dyHp source sum).
namespace Hp_producers {
constexpr int H_CR = 130;        // k_rxn[130] H  + CR   -> H+ + e
constexpr int H_CRph = 136;      // k_rxn[136] H  + CRph -> H+ + e
constexpr int H2_CR_HpH = 131;   // k_rxn[131] H2 + CR   -> H+ + H + e
constexpr int H2_CR_HpHm = 134;  // k_rxn[134] H2 + CR   -> H+ + H-
}  // namespace Hp_producers

// k_rxn[] channels producing He+ (dyHep source sum).
namespace Hep_producers {
constexpr int He_CR = 135;    // k_rxn[135] He + CR   -> He+ + e
constexpr int He_CRph = 138;  // k_rxn[138] He + CRph -> He+ + e
}  // namespace Hep_producers
}  // namespace cooling_idx

}  // namespace zero_metal

// ---------------------------------------------------------------------------
// metal_grain rate-kernel bookkeeping
// ---------------------------------------------------------------------------
namespace metal_grain {

// ---------------------------------------------------------------------------
// Network layout constants (forward-rate slot allocation in
// k_rxn[0..N_react-1]). Index convention: 0-based k_rxn[] index N == 1-based
// reaction id N+1. Block boundaries reflect the metal_grain reaction-number
// partition and the row layout of data/metal_grain.h5 (reactions/*).  A new
// model can re-anchor the same names with different ranges without searching
// call sites.
// ---------------------------------------------------------------------------
namespace slot {
// H, He, D, HD network: reactions 1..100  -> k_rxn[0..99]
constexpr int hhe_d_begin = 0;
constexpr int hhe_d_end = 100;  // exclusive

// UMIST 2012 metal network: reactions 101..543 -> k_rxn[100..542]
constexpr int metal_begin = 100;
constexpr int metal_end = 543;

// CR ionisation: reactions 544..552 -> k_rxn[543..551]
constexpr int cr_begin = 543;
constexpr int cr_end = 552;

// Additional UMIST: reactions 601..645 -> k_rxn[600..644]
constexpr int extra_begin = 600;
constexpr int extra_end = 645;

// CR photo-induced: reactions 656..682 -> k_rxn[655..681]
constexpr int cr_photo_begin = 655;
constexpr int cr_photo_end = 682;

// K / Na / Mg (KIDA): reactions 701..730 -> k_rxn[700..729]
constexpr int knamg_begin = 700;
constexpr int knamg_end = 730;

// Li (Bovino+2011 / Lepp+2002 / Mizusawa+2005): 801..830
constexpr int li_begin = 800;
constexpr int li_end = 830;

// Grain charge transfer: reactions 841..990 -> k_rxn[840..989] (150 entries)
constexpr int grain_charge_begin = 840;
constexpr int grain_charge_count = 150;

// Grain surface: reactions 991..1141 -> k_rxn[990..1140] (151 entries)
constexpr int grain_surface_begin = 990;
constexpr int grain_surface_count = 151;
}  // namespace slot

// ---------------------------------------------------------------------------
// Reaction-table row-position boundaries used in compute_rates<89,1200>.
// These reflect the row layout in data/metal_grain.h5: standard reactions
// occupy rows [0, n_std), CR first-order rows [n_std, n_std + n_cr_react),
// gas-grain charge rows [n_std + n_cr_react, +n_charge_react).
// ---------------------------------------------------------------------------
namespace loop {
constexpr int n_std = 635;
constexpr int n_cr_react = 36;       // 671 - 635
constexpr int n_charge_react = 150;  // 821 - 671
}  // namespace loop

// ---------------------------------------------------------------------------
// Reactions whose forward / Jacobian use a non-trivial form (3-body on grain,
// 4-body H-producing).  Kept here so any new model can swap or remove.
// ---------------------------------------------------------------------------
namespace special {
constexpr int rxn_2H_grain = 23;   // 2H + grain -> H2
constexpr int rxn_2D_grain = 144;  // 2D + grain -> HD
constexpr int rxn_4body_H = 273;   // adds extra H to product side
}  // namespace special

// ---------------------------------------------------------------------------
// Hot-spot k_rxn[] indices referenced by cooling / heating routines
// (solve/chemreact.h, solve/chemistry.h), organised by use in cooling
// expression:
//   dyHp     = recomb - Σ(Hp_producers)
//   dyHep    = recomb - Σ(Hep_producers)
//   rtgrain  = Σ(grain_H2 channels)
// Recombination indices live in arche::rxn:: (shared with zero_metal).
// ---------------------------------------------------------------------------
namespace cooling_idx {
// k_rxn[] channels producing H+ from neutral H/H2 (dyHp source sum).
namespace Hp_producers {
constexpr int H_CR = 543;     // k_rxn[543] H  + CR   -> H+ + e
constexpr int H_CRph = 676;   // k_rxn[676] H  + CRph -> H+ + e
constexpr int H2_CR = 547;    // k_rxn[547] H2 + CR   -> H2+ + e (H+ via dissoc)
constexpr int H2_CRph = 550;  // k_rxn[550] H2 + CRph
}  // namespace Hp_producers

// k_rxn[] channels producing He+ (dyHep source sum).
namespace Hep_producers {
constexpr int He_CR = 544;    // k_rxn[544] He + CR   -> He+ + e
constexpr int He_CRph = 677;  // k_rxn[677] He + CRph -> He+ + e
}  // namespace Hep_producers

// k_rxn[] indices for grain-catalysed H2 formation (used in rtgrain sum).
namespace grain_H2 {
constexpr int two_H = 22;             // k_rxn[22]   2H + grain -> H2 (3-body)
constexpr int surface_grgr = 1048;    // k_rxn[1048] grain + grain
constexpr int surface_grdust = 1108;  // k_rxn[1108] grain + dust
constexpr int H_dust_3body = 1107;    // k_rxn[1107] H + dust + H2
constexpr int H_gr_3body = 1114;      // k_rxn[1114] H + grain + H2
}  // namespace grain_H2
}  // namespace cooling_idx

// var[] range summed for CR heating (chemistry.h chem_step / chem_full_step).
// Loop body: for (int i = cr_heat_var_begin; i < cr_heat_var_end; ++i)
constexpr int cr_heat_var_begin = 543;
constexpr int cr_heat_var_end = 552;  // exclusive

}  // namespace metal_grain

// ---------------------------------------------------------------------------
// Compile-time invariants: relationships between the namespace constants
// declared above must hold for the cooling/solver code to operate correctly.
// These are coupled to the HDF5 data file (reactions/rean[] row order); any
// violation indicates that the data file layout and the C++ constants have
// drifted out of sync.  Updating either side without the other is a bug.
// ---------------------------------------------------------------------------
// ── zero_metal ──
static_assert(zero_metal::slot::std_end - zero_metal::slot::std_begin ==
                  zero_metal::loop::n_std,
              "zero_metal::loop::n_std must equal the standard-block width");
static_assert(zero_metal::slot::cr_end - zero_metal::slot::cr_begin ==
                  zero_metal::loop::n_cr_react,
              "zero_metal::loop::n_cr_react must equal the CR-block width");
static_assert(zero_metal::slot::cr_begin == zero_metal::cr_heat_var_begin,
              "zero_metal CR-heating range must start at the CR-block head");
static_assert(zero_metal::cr_heat_var_end - zero_metal::cr_heat_var_begin == 6,
              "zero_metal CR-heating sums exactly 6 direct ionisation channels "
              "(Hp/Hep producers via CR); update if data layout changes.");
static_assert(zero_metal::slot::cr_end <= zero_metal::N_react,
              "CR block must fit inside k_rxn[0..N_react-1]");

// ── metal_grain ──
static_assert(
    metal_grain::slot::cr_end - metal_grain::slot::cr_begin == 9,
    "metal_grain direct-CR block must be 9 channels (reaction ids 544..552); "
    "update data/metal_grain.h5 row layout if this fails.");
static_assert(
    metal_grain::slot::cr_photo_end - metal_grain::slot::cr_photo_begin == 27,
    "metal_grain CR-photo block must be 27 channels (reaction ids 656..682).");
static_assert(metal_grain::slot::cr_begin == metal_grain::cr_heat_var_begin,
              "metal_grain CR-heating range must start at the CR-block head");
static_assert(
    metal_grain::slot::cr_end == metal_grain::cr_heat_var_end,
    "metal_grain CR-heating range must cover the entire direct-CR block");
static_assert(
    metal_grain::slot::grain_charge_begin ==
        metal_grain::slot::cr_photo_end + (840 - 682),
    "grain_charge block starts after a fixed offset from CR-photo end");
static_assert(metal_grain::slot::grain_surface_begin ==
                  metal_grain::slot::grain_charge_begin +
                      metal_grain::slot::grain_charge_count,
              "grain_surface block must immediately follow grain_charge block");
static_assert(metal_grain::slot::grain_surface_begin +
                      metal_grain::slot::grain_surface_count <=
                  metal_grain::N_react,
              "Grain blocks must fit inside k_rxn[0..N_react-1]");

// ── shared rxn:: invariants (consistent across networks) ──
// k_rxn[0..99] in both data files stores the same 100 GA08 reactions.
// rxn:: indices must fall inside this shared block.
static_assert(rxn::H_e_to_Hp_2e >= 0 && rxn::H_e_to_Hp_2e < 100);
static_assert(rxn::Hp_rec_caseB >= 0 && rxn::Hp_rec_caseB < 100);
static_assert(rxn::Hep_rec >= 0 && rxn::Hep_rec < 100);
static_assert(rxn::H_e_to_Hm >= 0 && rxn::H_e_to_Hm < 100);
static_assert(rxn::Hm_H_to_H2_e >= 0 && rxn::Hm_H_to_H2_e < 100);
static_assert(rxn::Hm_e_to_H_2e >= 0 && rxn::Hm_e_to_H_2e < 100);
static_assert(rxn::Hm_Hp_to_2H >= 0 && rxn::Hm_Hp_to_2H < 100);
static_assert(rxn::Hm_H_to_2H_e >= 0 && rxn::Hm_H_to_2H_e < 100);
static_assert(rxn::H2p_H_to_H2_Hp >= 0 && rxn::H2p_H_to_H2_Hp < 100);
static_assert(rxn::three_H >= 0 && rxn::three_H < 100);
static_assert(rxn::H2H2_dis >= 0 && rxn::H2H2_dis < 100);
static_assert(rxn::H2_He_dis >= 0 && rxn::H2_He_dis < 100);
static_assert(rxn::Hep_H2_to_H2p_He >= 0 && rxn::Hep_H2_to_H2p_He < 100);

// ---------------------------------------------------------------------------
// Photon/dummy index used in the reaction table.
// The table encodes "no reactant/product" as index 0 or 101 in 1-based arrays.
// In C++ we map: index 0 → NONE, index 101 → PHOTON (kept in y_dum with value
// 1)
// ---------------------------------------------------------------------------
constexpr int IDX_NONE = 0;      // vacant reactant/product slot
constexpr int IDX_PHOTON = 101;  // photon (y_dum[101] = 1)

// ---------------------------------------------------------------------------
// RateForm — selects expression-form differences in rate coefficient routines
// (e.g., compute_HHeD_rates, compute_reverse_loop).
// ---------------------------------------------------------------------------
enum class RateForm { Primordial, MetalGrain };

}  // namespace arche
