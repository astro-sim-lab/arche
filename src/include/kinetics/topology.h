// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
#include <hdf5.h>

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/hdf5_utils.h"
#include "core/species_composition.h"  // composition::kMaxInvariantRows
#include "core/species_index.h"

namespace arche {

// ---------------------------------------------------------------------------
// ReactionTable<N_sp, N_react>
//
// Read-only after initialization. Shared across threads/cells.
// Topology comes from the per-model C++ network (models/<model>/reactions.h);
// the partition-function tables are loaded once at startup from a per-network
// HDF5 file (data/primordial.h5 or data/metal_grain.h5).
//
// Index convention (0-based species, sentinels appended after real species):
//   real species i      → index 0..N_sp-1
//   vacant slot         → index N_sp
//   photon              → index N_sp + 1
//   CR particle         → index N_sp + 2
//
// Reaction::num is a physical reaction ID kept 1-based (used as an offset
// into k_rxn[] via k_rxn[num - 1]).
// ---------------------------------------------------------------------------

// Standard / CR / charge-transfer reaction (Array-of-Structs record).
//
// Replaces the former parallel arrays rean[], rea1..rea3[], pro1..pro3[],
// nrea[], npro[], Cmass[], delE[]. Species indices are 0-based with sentinels
// appended after the real species (vacant = N_sp, photon = N_sp+1, CR =
// N_sp+2).
struct Reaction {
  int num = 0;                     // physical reaction ID (1-based)
  std::array<int, 3> reactants{};  // reactant species indices (was rea1..rea3)
  std::array<int, 3> products{};   // product species indices  (was pro1..pro3)
  int n_reactants = 0;             // body count for nH scaling (was nrea)
  int n_products = 0;              // product body count        (was npro)
  // Detailed-balance mass factor.  For a reversible reaction it is the
  // translational-mass part of the equilibrium constant,
  //   Cmass = (Π reactant masses / Π product masses)^(3/2)   [masses in grams],
  // the (2π k_B T/h²)^(3/2·Δn) and partition-function parts being applied
  // separately when the reverse rate is assembled.  It also acts as a
  // reversibility control: 0.0 forces the reverse rate to zero (forward-only
  // reactions, e.g. grain charge transfer), and 1.0 is used as a unit sentinel
  // where no mass-weighted reverse is intended.  Stored as a literal rather
  // than derived at run time: the tabulated values are reduced-precision and
  // the 0.0 / 1.0 sentinels are not reproducible from the masses alone.
  double Cmass = 0.0;
  double delE = 0.0;  // reaction energy |ΔE| [erg]
};

// Grain-surface reaction (metal_grain only): A + B -> C + D on grain.
// Replaces grean[], grea1/grea2[], gpro1/gpro2[], gnrea[].
struct GrainReaction {
  int num = 0;                     // grain reaction ID (was grean)
  std::array<int, 2> reactants{};  // reactant species indices (was grea1/grea2)
  std::array<int, 2> products{};   // product species indices  (was gpro1/gpro2)
  int n_reactants = 0;             // body count for nH scaling (was gnrea)
};

// Saha-equilibrium reaction (high-density branch, nH > 1e18).
// Replaces saha_num[], saha_rs1/rs2[], saha_ps1/ps2[], saha_nsr/nsp[],
// saha_Cm[], saha_dE[].
struct SahaReaction {
  int num = 0;  // reaction ID (was saha_num)
  std::array<int, 2>
      reactants{};                // reactant species indices (was saha_rs1/rs2)
  std::array<int, 2> products{};  // product indices, photon = N_sp+1 (saha_ps*)
  int n_reactants = 0;            // was saha_nsr
  int n_products = 0;             // was saha_nsp
  double Cmass = 0.0;  // detailed-balance mass factor (see Reaction::Cmass)
  double delE = 0.0;   // was saha_dE
};

// ---------------------------------------------------------------------------
// Reaction-record builders (constexpr).
//
// Used by the hand-maintained models/<model>/reactions.h tables.  They let a
// reaction be written as just its chemistry — the participating species (and
// the PHOTON / CR sentinels) in slot order — and reconstruct the stored record:
//
//   * slots are RIGHT-justified into the fixed 3- (reactant/product) or
//     2-wide (saha/grain) arrays, padding the left with the vacant sentinel;
//   * n_reactants / n_products are the count of real species (index < n_sp);
//     vacant / photon / CR sentinels are not counted.
//
// n_sp is passed as a plain argument so a single non-template builder serves
// every network (the model wrapper binds it).  These mirror the historical
// HDF5 row layout exactly and are checked field-by-field at migration time.
// ---------------------------------------------------------------------------
namespace topo {

constexpr Reaction make_reaction(int n_sp, int num,
                                 std::initializer_list<int> reactants,
                                 std::initializer_list<int> products,
                                 double Cmass, double delE) {
  Reaction r{};
  r.num = num;
  for (int i = 0; i < 3; ++i) r.reactants[i] = n_sp;  // vacant
  for (int i = 0; i < 3; ++i) r.products[i] = n_sp;
  int slot = 3 - static_cast<int>(reactants.size());
  for (int s : reactants) {
    r.reactants[slot++] = s;
    if (s < n_sp) ++r.n_reactants;
  }
  slot = 3 - static_cast<int>(products.size());
  for (int s : products) {
    r.products[slot++] = s;
    if (s < n_sp) ++r.n_products;
  }
  r.Cmass = Cmass;
  r.delE = delE;
  return r;
}

constexpr SahaReaction make_saha(int n_sp, int num,
                                 std::initializer_list<int> reactants,
                                 std::initializer_list<int> products,
                                 double Cmass, double delE) {
  SahaReaction s{};
  s.num = num;
  for (int i = 0; i < 2; ++i) s.reactants[i] = n_sp;
  for (int i = 0; i < 2; ++i) s.products[i] = n_sp;
  int slot = 2 - static_cast<int>(reactants.size());
  for (int x : reactants) {
    s.reactants[slot++] = x;
    if (x < n_sp) ++s.n_reactants;
  }
  slot = 2 - static_cast<int>(products.size());
  for (int x : products) {
    s.products[slot++] = x;
    if (x < n_sp) ++s.n_products;
  }
  s.Cmass = Cmass;
  s.delE = delE;
  return s;
}

constexpr GrainReaction make_grain(int n_sp, int num,
                                   std::initializer_list<int> reactants,
                                   std::initializer_list<int> products,
                                   int n_reactants) {
  GrainReaction g{};
  g.num = num;
  for (int i = 0; i < 2; ++i) g.reactants[i] = n_sp;
  for (int i = 0; i < 2; ++i) g.products[i] = n_sp;
  int slot = 2 - static_cast<int>(reactants.size());
  for (int x : reactants) g.reactants[slot++] = x;
  slot = 2 - static_cast<int>(products.size());
  for (int x : products) g.products[slot++] = x;
  g.n_reactants = n_reactants;  // independent body count (not the slot count)
  return g;
}

}  // namespace topo

template <int N_sp_, int N_react_>
struct ReactionTable {
  static constexpr int N_sp = N_sp_;
  static constexpr int N_react = N_react_;

  // Sentinel positions inside the extended-species index space.
  //
  // Reactant/product slots that do not reference a real species point at one
  // of these three indices. In the extended abundance vector y_ext (built by
  // make_y_ext() in kinetics/rates.h) each sentinel slot holds 1.0, the neutral
  // element of multiplication, so an unused reactant slot contributes a factor
  // of 1.0 to the rate product y_ext[r1] * y_ext[r2] * y_ext[r3] (i.e. no
  // contribution). This lets every reaction use a fixed 3-reactant / 3-product
  // shape regardless of how many species actually participate.
  static constexpr int IDX_VACANT =
      N_sp;  // empty (unused) reactant/product slot
  static constexpr int IDX_PHOTON =
      N_sp + 1;  // photon (photodissociation/-ionization)
  static constexpr int IDX_CR = N_sp + 2;  // cosmic-ray particle
  static constexpr int N_ext = N_sp + 3;   // total entries in y_ext / pf

  // Standard / CR / charge reactions (Array-of-Structs).
  // 0-based species indices; sentinels at N_sp, N_sp+1, N_sp+2.
  std::array<Reaction, N_react> reactions{};

  // Partition function tables: pf_table[isp] = vector of (T_K, pf_value) pairs
  // isp is 0-based; only species whose PfProvider is a NumericTable kind have a
  // non-empty vector (constant / analytic / BC16-table species ignore it).
  std::array<std::vector<std::pair<double, double>>, N_sp> pf_table{};

  // Number of reactions actually loaded.
  int n_loaded = 0;

  // ── Linear invariants of the network (solve/conservation.h) ───────────
  // Per-species element counts and charge — a basis of the LEFT NULL SPACE of
  // the stoichiometric matrix, but deliberately NOT an orthonormal one.  An
  // orthonormalised basis of the same space was implemented and measured, and
  // rejected because it raised cond(C W C^T) from ~1 to 2.4e8 and made the Li
  // residual 1000x worse (see solve/conservation.h, "WHY THE ROWS ARE ELEMENT
  // COUNTS").  Stored row-major, n_invariants rows of N_sp.
  // Every row v satisfies v^T S = 0 and
  // therefore d(v·y)/dt = 0 exactly, so v·y is conserved by the kinetics
  // whatever the rate coefficients.  Derived from `reactions` above by
  // conservation::fill_invariants() when the table is built, so it can never
  // assert an invariant the loaded network does not have.  n_invariants = 0
  // means no projection is applied.
  // Sized kMaxInvariantRows x N_sp, which is the most build_invariants() can
  // ever write.  (It arrived as N_sp x N_sp, which over-allocates by 4.6x at
  // N_sp = 23 and 30x at N_sp = 150 for rows that are never filled.)
  std::array<double,
             static_cast<std::size_t>(composition::kMaxInvariantRows) * N_sp_>
      invariants{};
  int n_invariants = 0;
  // Index of the charge row inside `invariants` (always the last row when
  // present), or -1.  It is the one row whose target is absolute (neutrality)
  // rather than "the value this step was handed".
  int charge_invariant_row = -1;

  // ── Reactions whose product list does not fit the record ──────────────
  // A Reaction holds at most three products, but reaction 273 of the full
  // metal_grain network has four:
  //     He+ + CH4 -> H2 + He + CH+ + H
  // The rate kernel injects the fourth rather than widening every record by a
  // slot (models/rate_kernel.h, the `sp_4body_H` branch).
  // build_invariants() reads `reactions` only, so it MUST be told about the
  // injected product: without it the reaction reads as unbalanced in H (4 in,
  // 3 out) and the whole matrix fails closed with n_invariants = 0, silently
  // disabling the projection for the entire network.
  // `extra_product_rxn` is the 1-based Reaction::num, or -1 when the model has
  // no such reaction (every primordial model, and Nakauchi2021_Minimal whose
  // keep-set drops 273).  `extra_product_sp` is the species index it adds one
  // of.  Only a single-species, single-count injection is representable; a
  // model needing more must extend this rather than leave it unstated, or the
  // fail-closed guarantee is lost.
  int extra_product_rxn = -1;
  int extra_product_sp = -1;

  // ── Grain surface reaction table (metal_grain only) ───────────────────
  static constexpr int N_grain = 200;
  std::array<GrainReaction, N_grain> grain_reactions{};
  int n_grain = 0;

  // ── Saha equilibrium data (high-density branch, nH > 1e18) ───────────
  static constexpr int N_saha = 100;
  std::array<SahaReaction, N_saha> saha{};
  int n_saha = 0;

  // ── Optional full-network back-reference (compact metal Minimal model) ──
  // The compact metal network gathers its rate coefficients from a full
  // metal_grain coefficient run (strategy-1 rate path) and its compact-Saha
  // partition functions from a full PF evaluation; both read this PF-loaded
  // full table.  Held as shared_ptr so the table stays copyable; null for every
  // other model, including the full metal table itself.
  std::shared_ptr<const ReactionTable<metal_grain::N_sp, metal_grain::N_react>>
      aux_full_metal{};
};

// Convenience aliases
using ZeroMetalTable = ReactionTable<zero_metal::N_sp, zero_metal::N_react>;
using MetalGrainTable = ReactionTable<metal_grain::N_sp, metal_grain::N_react>;

// ---------------------------------------------------------------------------
// load_pf_tables_h5()
//
// Loads the partition-function numerical tables from an HDF5 file into an
// already-topology-initialised ReactionTable.  Reaction topology (reactions,
// saha, grain_surface, mass) is no longer read from HDF5 — it is the
// hand-maintained C++ network in models/<model>/reactions.h, applied by that
// model's init_topology().  Only the partition-function tables remain data,
// because they are genuine numerical spectra (hundreds of (T, Q(T)) points).
//
//   /partition_functions/sp{i}/{T,pf}        (i = 0-based species index)
//
// The file's n_species attribute is still checked against the compile-time
// N_sp as a guard against pairing a table file with the wrong network.
// ---------------------------------------------------------------------------
template <int N_sp, int N_react>
void load_pf_tables_h5(ReactionTable<N_sp, N_react>& tbl,
                       const std::string& h5_path) {
  using h5utils::H5FileGuard;
  using h5utils::H5GroupExists;
  using h5utils::H5GroupGuard;
  using h5utils::H5OpenRead;
  using h5utils::H5Read1d;
  using h5utils::H5ReadIntAttr;

  H5FileGuard fid{H5OpenRead(h5_path)};
  if (fid < 0)
    throw std::runtime_error("Cannot open HDF5 chem table: " + h5_path);

  // Sanity: n_species must match the compile-time N_sp.
  int n_sp_file = H5ReadIntAttr(fid, "n_species");
  if (n_sp_file != N_sp) {
    throw std::runtime_error(
        "HDF5 chem table: n_species=" + std::to_string(n_sp_file) + " in " +
        h5_path + " does not match compile-time N_sp=" + std::to_string(N_sp));
  }

  // ── partition_functions ───────────────────────────────────────────────
  if (H5GroupExists(fid, "partition_functions")) {
    H5GroupGuard pfg{H5Gopen2(fid, "partition_functions", H5P_DEFAULT)};
    // Iterate over child groups named "spK" with K in 0..N_sp-1.
    for (int isp = 0; isp < N_sp; ++isp) {
      std::string subname = "sp" + std::to_string(isp);
      if (!H5GroupExists(pfg, subname)) continue;
      H5GroupGuard sg{H5Gopen2(pfg, subname.c_str(), H5P_DEFAULT)};
      auto T = H5Read1d(sg, "T");
      auto pf = H5Read1d(sg, "pf");
      if (T.size() != pf.size()) {
        throw std::runtime_error("HDF5 chem table: pf sp" +
                                 std::to_string(isp) +
                                 " has T/pf length mismatch");
      }
      tbl.pf_table[isp].clear();
      tbl.pf_table[isp].reserve(T.size());
      for (size_t k = 0; k < T.size(); ++k)
        tbl.pf_table[isp].emplace_back(T[k], pf[k]);
    }
  }
}

}  // namespace arche
