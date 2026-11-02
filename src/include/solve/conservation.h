// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// conservation.h — restore the network's element and charge budgets after an
// implicit chemistry step.
//
// WHY THIS EXISTS
// ---------------
// The backward-Euler step solves  dy = r_f(y_init + dy)*dt  by Newton.  The
// increment it computes is conservative *by construction*: for any element (or
// charge) count vector c we have c·r_f = 0, hence c^T J = 0, hence
// c^T A = c^T for A = I - dt*J, hence c^T A^{-1} = c^T, and the Newton update
// ddy = A^{-1}(r_f*dt - dy) drives c·dy to zero at every iteration up to the
// rounding of the linear solve.
//
// The increment is nevertheless *applied* as  y[i] += ddy[i], which rounds at
// ulp(y[i]).  Species sharing an element sit decades apart in magnitude
// (ulp(y_D) = 3.4e-21 against ulp(y_HD) = 1.3e-23 and ulp(y_D+) = 1.3e-29 in a
// typical collapse), so those roundings cannot cancel between them, and the
// element total absorbs up to ~1 ulp of its dominant carrier per step.  Over the
// many steps of a collapse these roundings accumulate; the projection below
// bounds the accumulated drift whatever its statistics.
//
// WHAT THIS DOES
// --------------
// After the implicit solve, project y onto the affine set that restores every
// element total and neutralises the net charge:
//
//     find  Δ  minimising  Σ_i Δ_i² / w_i    subject to    C Δ = r
//
// with C the (m x N_sp) matrix whose rows are the per-species element counts,
// the number of grains and the per-species charge, w_i = max(y_i, 0), and the
// target r
//
//     element rows : r_e = Σ_i C[e][i] (y_in[i] − y[i])   (+ carry, below)
//     charge row   : r_q = − Σ_i q_i y[i]
//
// The solution is Δ = W C^T (C W C^T)^{-1} r: one m x m symmetric
// positive-definite solve, m <= 11 for the shipped networks.
//
// ELEMENTS ARE CONSERVED, CHARGE IS NEUTRALISED — a deliberate asymmetry.
// An element total has no reference the kernel knows (a fluid code advects
// composition in and out), so the only defensible target is "whatever this step
// was handed".  Net charge does have one: a plasma is quasi-neutral, arche's
// own initial conditions set y_e = Σ_ion q y, and the collapse app's
// charge_imbal diagnostic is written on the assumption that it is ~0.  Holding
// the input value instead would freeze whatever float64 net charge the early,
// electron-rich phase accumulated, while the normalisation Σ_{q>0} q y falls by
// orders of magnitude along the collapse.
// A charged input is neutralised once, before the first step, by this projection
// (chem_integrate_adaptive in solve/chemistry.h) and rejected when that needs more
// than kMaxRelShift.  The guard bounds the per-species RELATIVE shift, and for a pure charge residual
// shift_i = q_i * lambda_q with lambda_q = r_q / SUM_i y_i q_i^2, so the
// projection declines once
//     |SUM_i q_i y_i|  >  (kMaxRelShift / max_i|q_i|) * SUM_i y_i q_i^2
// = 3.33e-3 * SUM_i y_i q_i^2 for the 23-species network and 1.0e-2 * SUM_i
// y_i q_i^2 for the 15-species one.  The two differ only because the 23-species
// catalog carries Li+++ (q = +3) and the compact one stops at singly charged
// ions.
//
// THE SUB-ULP REMAINDER IS CARRIED (element rows only).  r_e is typically a
// fraction of ulp(y_dominant), so applying it rounds to zero as often as not, and
// the projection's own rounding would then accumulate like the drift it repairs.
// Feeding the unapplied part back into the next step's target — a first-order
// delta-sigma accumulator — bounds the drift at ~1 ulp of the element total for
// any number of steps.  The charge row needs no carry: its target is absolute,
// so each step already sees the full remaining error.
//
// THE WEIGHT IS A CHOICE, NOT A DERIVATION.  w_i = y_i makes the *relative*
// shift Δ_i/y_i identical for every species carrying a given element, so a
// trace species 15 decades below its element total moves by ~1e-16 of its own
// value instead of absorbing an absolute share of the correction — which is
// what an unweighted (w_i = 1) least-norm correction would do, and which would
// move such a species by many decades.  It also returns essentially the whole
// absolute correction to the dominant carrier, i.e. to the species whose
// rounding produced the residual.  It cannot flip a sign: Δ_i = y_i·s_i with
// |s_i| bounded by kMaxRelShift, and a species at exactly zero has w_i = 0 and
// is never moved off zero.
//
// ELEMENTS AND CHARGE ARE ENFORCED TOGETHER, NOT IN SEQUENCE.  They are
// independent constraints — the electron carries charge and no element, LiH
// carries two elements and no charge — and repairing them one after the other
// would undo the earlier repair.  Both are rows of the same C, so the single
// least-norm solution satisfies all of them simultaneously.
//
// WHY THE ROWS ARE ELEMENT COUNTS AND NOT SOME OTHER BASIS OF THE SAME SPACE.
// The conserved subspace is exactly the left null space of the stoichiometric
// matrix, and any basis of it defines the same constraint set in exact
// arithmetic.  It does not define the same computation: the element totals span
// 1 (H) to 4.65e-10 (Li), and a basis that mixes elements (an orthonormalised
// one, for example) adds rows of very different scale, which ill-conditions
// C W C^T and buries the small totals in the rounding of the large ones.  In the
// element/charge basis each row carries one element, so rows of different scale
// do not mix and the Jacobi scaling below removes the scale differences.
//
// SCOPE.  (1) Every shipped network.  The invariant matrix is built from
// core/species_composition.h, which registers all 95 catalog species over nine
// elements, plus one row for the number of grains (composition::grain_count: a
// grain changes its charge, never its number); a model whose species are not all
// registered gets n_invariants = 0, and then every step is rejected.
// n_invariants: 5 primordial (full and Minimal; no grains), 9
// Nakauchi2021_Minimal (the compact metal network carries no K or Na, so those
// two rows are pruned), 11 Nakauchi2021.  The metals locked in the dust (C_s,
// O_s, Mg_s, Na_s, K_s) count toward their elements, so the metal rows conserve
// the gas + dust totals.  (2) The implicit chemistry solve only: the projection
// runs before the operator-split channels of chem_full_step() (X-ray, and
// Lyman-Werner after a Saha step), whose stoichiometry is outside the reaction
// table this file validates against.  In the network steps the Lyman-Werner
// channels are part of the step's rates and keep every element and the charge.
//
// The residual floor is a property of the DIAGNOSTIC as well as of the solver:
// this file never forms an element total in float64 (it accumulates increments
// plus carry in long double), so the floor a checker reports depends on how it
// sums (compensated or naive left-to-right).  Name the summation convention when
// quoting one.
//
// KNOWN LIMITS.
//  * Above nH_eq = 1e18 cm^-3 with zeta == 0 and zeta_X == 0, the step takes
//    a Saha equilibrium branch that RE-ANCHORS the He, D and Li totals to
//    hard-coded reference abundances, while this projection's target is the
//    totals the step was handed.  When the totals agree with the references the
//    two coincide; when anything else has moved the totals, the projection keeps
//    the handed totals and the Saha branch's re-anchoring no longer holds.  The
//    branch drives the H total to 1.0 only to its Newton tolerance.
//  * An empty row (all carriers exactly zero) with nothing owed is dropped from
//    that step's solve (the ACTIVE SET below); otherwise it would make the
//    weighted matrix singular and decline the whole projection.  An empty row
//    that IS owed a correction, and a repair beyond kMaxRelShift, remain decline
//    paths.
//  * The projection changes abundances only, but those feed the ionisation terms
//    and hence the adaptive timestep, so the density schedule of a collapse can
//    move at the rounding level.
//  * The long double accumulations below are this codebase's only arithmetic
//    whose width depends on the ABI (80-bit extended on x86-64 SysV, but 64-bit
//    on some others).  Deliberately NOT guarded by a static_assert: failing the
//    build on such a target is worse than the degradation, which is confined to
//    the charge row (the element rows sum increments, not totals, and provably
//    do not need the extra width).
//  * carry lives in ChemCell, not in ChemState, so the documented copy-in /
//    copy-out cell API can neither see nor reset it.  A caller that reuses one
//    ChemCell as scratch across many fluid cells will replay a foreign
//    remainder into the next cell's first step; call ChemCell::reset_var()
//    between cells.  The error this bounds is <= half an ulp of the element
//    total (derived).
//  * build_invariants() balances all three stoichiometry tables: tbl.reactions,
//    tbl.saha (n_saha = 18 full / 9 minimal for primordial, 53 for
//    metal_grain) and tbl.grain_reactions (148 for metal_grain; primordial has
//    n_grain = 0).  saha drives the branch that rewrites y above nH_eq and
//    grain_reactions moves species between the gas and the mantle, so a typo
//    in either is as consequential as one in reactions.  Every shipped table
//    balances (n_invariants > 0 for all four models).
//
// ---------------------------------------------------------------------------
#include <array>
#include <cmath>
#include <cstddef>

#include "core/species_composition.h"
#include "kinetics/topology.h"  // Reaction, ReactionTable

namespace arche {
namespace conservation {

//: Largest relative per-species shift the projection will apply.  A violation
//: that would need a larger repair is a defect to be seen, not to be papered
//: over, so the projection declines and leaves both the state and the evidence
//: untouched.
inline constexpr double kMaxRelShift = 1.0e-2;

//: Upper bound on the number of invariant rows: one per tracked element, one for
//: the number of grains, plus charge.  Sizes the per-cell carry (ChemCell::cons_carry) and the scratch
//: normal matrix below.  Defined next to the element list itself so that this
//: and ReactionTable::invariants cannot drift apart.
inline constexpr int kMaxRows = composition::kMaxInvariantRows;

// ---------------------------------------------------------------------------
// build_invariants — element-count and charge rows for one model's species.
//
//   SpeciesSetT : the model's species selection (local index -> SpId).
//   tbl         : the loaded reaction table, used to VERIFY the rows.
//   out         : receives m rows of N_sp, row-major; needs N_sp*N_sp doubles.
//   charge_row  : receives the index of the charge row, or -1 if the model
//                 carries no charged species.  Element rows come first, so the
//                 charge row is always the last one when present.
//
// Returns m, or 0 if the rows cannot be trusted, in which case no projection
// must be attempted.  Zero is returned when
//   * any species of the model has no registered composition, or
//   * any loaded reaction fails to balance under those compositions.
// The second test is the important one: it re-derives, from the table actually
// loaded, that the rows really are conserved quantities of this network.  An
// element absent from the model (all counts zero) contributes no row.
// ---------------------------------------------------------------------------
template <class SpeciesSetT, class Table>
int build_invariants(const Table& tbl, double* out, int& charge_row) {
  charge_row = -1;
  constexpr int N_sp = Table::N_sp;
  // elements, the number of grains, charge (always last)
  constexpr int N_rows = composition::kNElements + 2;
  constexpr int kGrainRow = composition::kNElements;
  constexpr int kChargeRow = composition::kNElements + 1;

  // Per-species integer composition, in the model's local index order.
  std::array<std::array<int, N_rows>, N_sp> comp{};
  for (int i = 0; i < N_sp; ++i) {
    const composition::Composition c =
        composition::composition_of(SpeciesSetT::canonical(i));
    if (!c.known) return 0;
    for (int e = 0; e < composition::kNElements; ++e) comp[i][e] = c.n[e];
    comp[i][kGrainRow] = composition::grain_count(SpeciesSetT::canonical(i));
    comp[i][kChargeRow] = c.q;
  }

  // One record balances when every tracked row nets to zero.  Unused slots hold
  // the VACANT/PHOTON/CR sentinel at index >= N_sp and drop out; n_reactants /
  // n_products are body counts for the nH scaling, not slot counts, so they
  // must not bound these loops.
  // `extra_sp` is a product the record could not hold, injected by the rate
  // kernel instead (see ReactionTable::extra_product_rxn); it is >= 0 only for
  // the one declared exception.  Balancing that record on its stored slots
  // alone would read it as element-violating and fail the whole matrix closed,
  // which is why the table declares the injection rather than the kernel owning
  // it privately.
  // Slot counts come from the record types (.size()), never from a literal:
  // Reaction holds three slots per side, SahaReaction and GrainReaction two, and
  // a literal here would silently under-check if a record ever gained a slot.
  const auto balanced = [&comp](const int* rea, int n_rea, const int* pro,
                                int n_pro, int extra_sp) {
    std::array<int, N_rows> bal{};
    for (int k = 0; k < n_rea; ++k) {
      const int i = rea[k];
      if (i >= 0 && i < N_sp)
        for (int e = 0; e < N_rows; ++e) bal[e] -= comp[i][e];
    }
    for (int k = 0; k < n_pro; ++k) {
      const int i = pro[k];
      if (i >= 0 && i < N_sp)
        for (int e = 0; e < N_rows; ++e) bal[e] += comp[i][e];
    }
    if (extra_sp >= 0 && extra_sp < N_sp)
      for (int e = 0; e < N_rows; ++e) bal[e] += comp[extra_sp][e];
    for (int e = 0; e < N_rows; ++e)
      if (bal[e] != 0) return false;
    return true;
  };

  // Every loaded record of every stoichiometry table the ReactionTable carries
  // must balance.  Checking `reactions` alone leaves the guarantee incomplete:
  // `saha` drives the branch that rewrites y above nH_eq and `grain_reactions`
  // moves species between the gas and the mantle, so a typo in either is just
  // as consequential.
  for (int r = 0; r < tbl.n_loaded && r < Table::N_react; ++r) {
    const Reaction& rx = tbl.reactions[r];
    const int extra =
        (rx.num == tbl.extra_product_rxn) ? tbl.extra_product_sp : -1;
    if (!balanced(rx.reactants.data(), static_cast<int>(rx.reactants.size()),
                  rx.products.data(), static_cast<int>(rx.products.size()),
                  extra))
      return 0;
  }
  for (int r = 0; r < tbl.n_saha && r < Table::N_saha; ++r) {
    const SahaReaction& sr = tbl.saha[r];
    if (!balanced(sr.reactants.data(), static_cast<int>(sr.reactants.size()),
                  sr.products.data(), static_cast<int>(sr.products.size()), -1))
      return 0;
  }
  for (int r = 0; r < tbl.n_grain && r < Table::N_grain; ++r) {
    const GrainReaction& gr = tbl.grain_reactions[r];
    if (!balanced(gr.reactants.data(), static_cast<int>(gr.reactants.size()),
                  gr.products.data(), static_cast<int>(gr.products.size()), -1))
      return 0;
  }

  int m = 0;
  for (int e = 0; e < N_rows; ++e) {
    bool any = false;
    for (int i = 0; i < N_sp && !any; ++i) any = (comp[i][e] != 0);
    if (!any) continue;  // element (or charge) not carried by this model
    double* row = out + static_cast<std::size_t>(m) * N_sp;
    for (int i = 0; i < N_sp; ++i) row[i] = static_cast<double>(comp[i][e]);
    if (e == kChargeRow) charge_row = m;
    ++m;
  }
  return m;
}

// ---------------------------------------------------------------------------
// fill_invariants — populate tbl.invariants / tbl.n_invariants.  Call once,
// when the table is built: the result depends only on the topology.
// ---------------------------------------------------------------------------
template <class SpeciesSetT, class Table>
void fill_invariants(Table& tbl) {
  tbl.n_invariants = build_invariants<SpeciesSetT>(tbl, tbl.invariants.data(),
                                                   tbl.charge_invariant_row);
  if (tbl.n_invariants == 0) tbl.charge_invariant_row = -1;
}

// ---------------------------------------------------------------------------
// project — restore the element totals and neutralise the net charge.
//
//   inv        : m x N_sp element/charge rows, row-major (build_invariants)
//   charge_row : index of the neutrality row, or -1
//   y_in       : abundances entering the step
//   y          : abundances leaving the step, corrected in place
//   carry      : per-cell sub-ulp remainder, updated in place (kMaxRows)
//
// Returns true if the correction was computed and applied (including the case
// where there was nothing to correct), false if the projection declined: no
// invariant rows, a weighted normal matrix that is not positive definite, or a
// repair exceeding kMaxRelShift.  Declining leaves y and carry exactly as they
// were, i.e. the behaviour without this file, so it is never worse than not
// having called it.
// ---------------------------------------------------------------------------
template <int N_sp>
bool project_nonzero(const double* inv, int m, int charge_row,
                     const std::array<double, N_sp>& y_in,
                     std::array<double, N_sp>& y,
                     std::array<double, kMaxRows>& carry);

// An element whose carriers are all exactly zero entering the step, with nothing
// carried, has a total of exactly zero; its carriers are non-negative, so each of them
// leaves the step at zero.  Whatever the step left on them is rounding of the linear
// solve (a carrier's row couples to other carriers of the same element, all zero), and
// a correction of it would have to remove all of it -- a relative repair of -1, which
// project_nonzero declines.  They are set to zero first; the charge row then accounts
// for any charge they held.
template <int N_sp>
bool project(const double* inv, int m, int charge_row,
             const std::array<double, N_sp>& y_in, std::array<double, N_sp>& y,
             std::array<double, kMaxRows>& carry) {
  if (m <= 0 || m > N_sp || m > kMaxRows) return false;
  bool cleared = false;
  std::array<double, N_sp> yc = y;
  for (int k = 0; k < m; ++k) {
    if (k == charge_row || carry[k] != 0.0) continue;
    const double* c = inv + static_cast<std::size_t>(k) * N_sp;
    bool empty = true;
    for (int i = 0; i < N_sp && empty; ++i)
      if (c[i] != 0.0 && y_in[i] != 0.0) empty = false;
    if (!empty) continue;
    for (int i = 0; i < N_sp; ++i)
      if (c[i] != 0.0 && yc[i] != 0.0) {
        yc[i] = 0.0;
        cleared = true;
      }
  }
  if (!cleared)
    return project_nonzero<N_sp>(inv, m, charge_row, y_in, y, carry);
  // Declining leaves y as it was, so the cleared state is only kept on success.
  if (!project_nonzero<N_sp>(inv, m, charge_row, y_in, yc, carry)) return false;
  y = yc;
  return true;
}

template <int N_sp>
bool project_nonzero(const double* inv, int m, int charge_row,
                     const std::array<double, N_sp>& y_in,
                     std::array<double, N_sp>& y,
                     std::array<double, kMaxRows>& carry) {
  // Targets.  Element rows are formed from the per-step INCREMENTS and never as
  // the difference of two element totals: the increments are ~1e9 times smaller
  // than the totals here, so this expression carries the residual at full
  // precision while differencing two float64 totals would return their
  // summation noise instead.  The charge row is absolute (neutrality), so it
  // has no such shortcut and is accumulated in long double to keep the
  // cancellation between the electron and ion terms from setting the floor.
  std::array<double, kMaxRows> r{};
  std::array<bool, kMaxRows>
      moved_by_step{};  // element rows: this step's increment != 0
  bool all_zero = true;
  for (int k = 0; k < m; ++k) {
    const double* c = inv + static_cast<std::size_t>(k) * N_sp;
    long double s = 0.0L;
    if (k == charge_row) {
      for (int i = 0; i < N_sp; ++i)
        s -= static_cast<long double>(c[i]) * static_cast<long double>(y[i]);
    } else {
      for (int i = 0; i < N_sp; ++i)
        s += static_cast<long double>(c[i]) *
             (static_cast<long double>(y_in[i]) -
              static_cast<long double>(y[i]));
      moved_by_step[k] = (s != 0.0L);
      s += static_cast<long double>(carry[k]);
    }
    r[k] = static_cast<double>(s);
    if (r[k] != 0.0) all_zero = false;
  }
  if (all_zero) {
    // Nothing to correct, so nothing is owed.  The general path below closes
    // each element row with carry[k] = r[k] - done[k]; here r[k] == 0 and
    // done[k] == 0, so the same rule gives carry[k] = 0.  Returning without
    // clearing would leave the previous step's remainder in the accumulator and
    // request it a second time on the next step.
    carry.fill(0.0);
    return true;
  }

  // G = C W C^T, W = diag(max(y,0)); symmetric positive semi-definite.
  std::array<double, N_sp> wgt{};
  for (int i = 0; i < N_sp; ++i) wgt[i] = (y[i] > 0.0) ? y[i] : 0.0;
  // Leading dimension of the scratch normal matrix.  Only the m x m block is
  // ever touched (m <= kMaxRows), so this is sized kMaxRows^2 rather than the
  // N_sp^2 it arrived as -- 25 doubles instead of 22500 at N_sp = 150.  The
  // arithmetic is unchanged: same entries in the same order, only the stride
  // moves.
  constexpr int kLd = kMaxRows;
  std::array<double, kMaxRows * kMaxRows> g{};
  for (int k = 0; k < m; ++k) {
    const double* ck = inv + static_cast<std::size_t>(k) * N_sp;
    for (int l = 0; l <= k; ++l) {
      const double* cl = inv + static_cast<std::size_t>(l) * N_sp;
      double s = 0.0;
      for (int i = 0; i < N_sp; ++i) s += wgt[i] * ck[i] * cl[i];
      g[k * kLd + l] = s;
      g[l * kLd + k] = s;
    }
  }

  // ACTIVE SET.  A row whose weight is identically zero -- every carrier of that
  // element absent from y (a metal-free composition in a metal network, or the
  // grains after complete evaporation) -- puts a zero on G's diagonal.
  // Declining the whole projection there would discard the other rows too.
  //   gkk == 0 and r == 0 -> nothing to conserve and nothing owed; drop the row
  //                          and solve the rest.
  //   gkk == 0 and r != 0 -> the step changed an element that no surviving
  //                          species carries.  Unrepairable: decline.
  //   gkk == 0, only the carry owed -> the step moved none of it (increment
  //                          exactly 0), so the total is exactly zero on both
  //                          sides and the remainder belongs to a total that no
  //                          longer exists -- the caller removed every carrier
  //                          between steps (the library zeroes the grains at
  //                          complete evaporation, solve/grain_evaporation.h).
  //                          Drop the row and clear its carry.
  // A diagonal regulariser is NOT the fix.  It would manufacture a correction
  // and spread it over species that cannot carry the element at all.
  std::array<int, kMaxRows> act{};
  int n_act = 0;
  for (int k = 0; k < m; ++k) {
    if (g[k * kLd + k] > 0.0) {
      act[n_act++] = k;
    } else if (r[k] != 0.0) {
      if (k == charge_row || moved_by_step[k]) return false;
      r[k] = 0.0;  // the carry loop below then stores 0 for this row
    }
  }
  if (n_act == 0) {
    // Every row empty and nothing owed on any of them.
    carry.fill(0.0);
    return true;
  }
  // Compact the active block into the leading n_act x n_act corner.  act[] is
  // strictly increasing, so act[a] >= a and act[b] >= b: every cell read sits
  // at or after the cell written, and no live value is overwritten first.
  for (int a = 0; a < n_act; ++a)
    for (int b = 0; b < n_act; ++b) g[a * kLd + b] = g[act[a] * kLd + act[b]];
  // r itself must survive: the carry loop below is indexed by the ORIGINAL row.
  std::array<double, kMaxRows> ra{};
  for (int a = 0; a < n_act; ++a) ra[a] = r[act[a]];

  // Jacobi (symmetric diagonal) scaling.  The element totals span nine decades,
  // so G's raw condition number is ~1e9 while the scaled one is ~1; without
  // this the smallest element's residual is lost in the solve.
  std::array<double, kMaxRows> d{};
  for (int k = 0; k < n_act; ++k) {
    const double gkk = g[k * kLd + k];
    if (!(gkk > 0.0)) return false;
    d[k] = std::sqrt(gkk);
  }
  for (int k = 0; k < n_act; ++k)
    for (int l = 0; l < n_act; ++l) g[k * kLd + l] /= (d[k] * d[l]);

  // Cholesky of the scaled matrix, then solve for the scaled multipliers.
  for (int k = 0; k < n_act; ++k) {
    double v = g[k * kLd + k];
    for (int i = 0; i < k; ++i) v -= g[k * kLd + i] * g[k * kLd + i];
    if (!(v > 0.0)) return false;
    v = std::sqrt(v);
    g[k * kLd + k] = v;
    const double invv = 1.0 / v;
    for (int j = k + 1; j < n_act; ++j) {
      double s = g[j * kLd + k];
      for (int i = 0; i < k; ++i) s -= g[j * kLd + i] * g[k * kLd + i];
      g[j * kLd + k] = s * invv;
    }
  }
  std::array<double, kMaxRows> lam_act{};
  for (int k = 0; k < n_act; ++k) {
    double s = ra[k] / d[k];
    for (int i = 0; i < k; ++i) s -= g[k * kLd + i] * lam_act[i];
    lam_act[k] = s / g[k * kLd + k];
  }
  for (int k = n_act - 1; k >= 0; --k) {
    double s = lam_act[k];
    for (int i = k + 1; i < n_act; ++i) s -= g[i * kLd + k] * lam_act[i];
    lam_act[k] = s / g[k * kLd + k];
  }
  for (int k = 0; k < n_act; ++k) lam_act[k] /= d[k];  // undo the scaling
  // Scatter back to the original row numbering; a dropped row contributes
  // nothing to the shift, which is what "not conserved this step" means.
  std::array<double, kMaxRows> lam{};
  for (int k = 0; k < n_act; ++k) lam[act[k]] = lam_act[k];

  // shift[i] = (C^T lambda)_i is the RELATIVE change applied to species i.
  std::array<double, N_sp> shift{};
  for (int i = 0; i < N_sp; ++i) shift[i] = 0.0;
  for (int k = 0; k < m; ++k) {
    const double* c = inv + static_cast<std::size_t>(k) * N_sp;
    const double lk = lam[k];
    for (int i = 0; i < N_sp; ++i) shift[i] += c[i] * lk;
  }
  for (int i = 0; i < N_sp; ++i)
    if (!std::isfinite(shift[i]) || std::abs(shift[i]) > kMaxRelShift)
      return false;

  // Apply, and record what was actually stored.  y[i] and its corrected value
  // are within a factor of two of each other, so (y_new − y_old) is exact by
  // Sterbenz's lemma and `moved` is the change the file will show — not the
  // change that was requested.
  std::array<double, N_sp> moved{};
  for (int i = 0; i < N_sp; ++i) {
    const double before = y[i];
    y[i] += wgt[i] * shift[i];
    moved[i] = y[i] - before;
  }

  // Carry the part of each element target that rounding refused to store, so
  // the next step asks for it again.  Without this the leftovers random-walk
  // and the totals drift again, one decade lower.  The charge row's target is
  // absolute and is re-derived from y every step, so it carries nothing.
  for (int k = 0; k < m; ++k) {
    if (k == charge_row) {
      carry[k] = 0.0;
      continue;
    }
    const double* c = inv + static_cast<std::size_t>(k) * N_sp;
    long double done = 0.0L;
    for (int i = 0; i < N_sp; ++i)
      done +=
          static_cast<long double>(c[i]) * static_cast<long double>(moved[i]);
    carry[k] = static_cast<double>(static_cast<long double>(r[k]) - done);
  }
  return true;
}

}  // namespace conservation
}  // namespace arche
