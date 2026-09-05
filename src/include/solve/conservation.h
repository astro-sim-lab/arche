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
// chemreact() solves  dy = r_f(y_init + dy)*dt  by Newton-Raphson.  The
// increment it computes is conservative *by construction*: for any element (or
// charge) count vector c we have c·r_f = 0, hence c^T J = 0, hence
// c^T A = c^T for A = I - dt*J, hence c^T A^{-1} = c^T, and the Newton update
// ddy = A^{-1}(r_f*dt - dy) drives c·dy to zero at every iteration.  Measured
// on the one-zone primordial collapse, |c·dy| / (element total) stays between
// 1e-19 and 1e-17 at every step size tried — the linear algebra is not the
// problem.
//
// The increment is nevertheless *applied* as  y[i] += ddy[i], which rounds at
// ulp(y[i]).  Species sharing an element sit decades apart in magnitude
// (ulp(y_D) = 3.4e-21 against ulp(y_HD) = 1.3e-23 and ulp(y_D+) = 1.3e-29 in a
// typical collapse), so those roundings cannot cancel between them, and the
// element total absorbs ~1 ulp of its dominant carrier per step.  The per-step
// error is not white — its measured lag-1 autocorrelation reaches +0.48 — so
// the sum grows faster than sqrt(N_steps): at N = 2.3e6 the deuterium residual
// reaches 1.6e-11, about 3000x its value at N = 2.3e4.  Superlinear in step
// count, which is what the projection below exists to bound.
// Measured on 15-species primordial collapses at nH <= 1e8 cm^-3 with g++
// 13.3.0.  An independent measurement on the 89-species metal network at
// nH <= 1e23 with clang-22 reports a sqrt(N_steps) random walk instead; the two
// regimes are not reconciled.  The projection does not depend on which is right
// -- it re-imposes the constraint at the end of the step whatever produced the
// drift.
//
// WHAT THIS DOES
// --------------
// After the implicit solve, project y onto the affine set that restores every
// element total and neutralises the net charge:
//
//     find  Δ  minimising  Σ_i Δ_i² / w_i    subject to    C Δ = r
//
// with C the (m x N_sp) matrix whose rows are the per-species element counts
// and the per-species charge, w_i = max(y_i, 0), and the target r
//
//     element rows : r_e = Σ_i C[e][i] (y_in[i] − y[i])   (+ carry, below)
//     charge row   : r_q = − Σ_i q_i y[i]
//
// The solution is Δ = W C^T (C W C^T)^{-1} r: one m x m symmetric
// positive-definite solve, m <= 5 for the primordial networks.
//
// ELEMENTS ARE CONSERVED, CHARGE IS NEUTRALISED — a deliberate asymmetry.
// An element total has no reference the kernel knows (a fluid code advects
// composition in and out), so the only defensible target is "whatever this step
// was handed".  Net charge does have one: a plasma is quasi-neutral, arche's
// own initial conditions set y_e = Σ_ion q y, and the collapse app's
// charge_imbal diagnostic is written on the assumption that it is ~0.  Chasing
// the input value instead was implemented and measured: it freezes whatever
// float64 net charge the early, electron-rich phase happened to accumulate, and
// because the normalisation Σ_{q>0} q y then falls five decades along the
// collapse, the normalised charge residual came out 15x WORSE than leaving
// the charge row unprojected.
// A caller that really does hand in a charged cell is protected by
// kMaxRelShift.  The guard bounds the per-species RELATIVE shift, and for a
// pure charge residual shift_i = q_i * lambda_q with lambda_q = r_q / SUM_i
// y_i q_i^2, so the projection declines once
//     |SUM_i q_i y_i|  >  (kMaxRelShift / max_i|q_i|) * SUM_i y_i q_i^2
// = 3.33e-3 * SUM_i y_i q_i^2 for the 23-species network and 1.0e-2 * SUM_i
// y_i q_i^2 for the 15-species one.  The two differ only because the 23-species
// catalog carries Li+++ (q = +3) and the compact one stops at singly charged
// ions; the guard is therefore 3x tighter for the full network on account of a
// species whose abundance never exceeds ~2e-19.
//
// THE SUB-ULP REMAINDER IS CARRIED (element rows only).  r_e is typically a
// fraction of ulp(y_dominant), so applying it rounds to zero as often as not.
// Dropping that remainder reproduces the original failure one decade down: the
// projection's own rounding random-walks, and with the measured sign bias the
// deuterium residual still reached 1e-12 at 2.3e6 steps (16x better than
// no projection at all, but not a fix).  Feeding the unapplied part back into
// the next step's
// target — a first-order delta-sigma accumulator — bounds the drift at ~1 ulp
// of the element total for any number of steps.  The charge row needs no carry:
// its target is absolute, so each step already sees the full remaining error.
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
// arithmetic.  It does not define the same computation: with y taken from a
// real collapse, the Jacobi-scaled normal matrix C W C^T has condition number
// 1.02 in the element/charge basis and 2.4e8 in an orthonormalised basis of the
// same space, because the element totals span 1 (H) to 4.65e-10 (Li) and any
// mixing buries the small ones.  An orthonormalised basis was implemented,
// measured, and rejected: it made the Li residual 1000x *worse*.
//
// SCOPE.  (1) Every shipped network.  The invariant matrix is built from
// core/species_composition.h, which registers all 89 catalog species over nine
// elements; a model whose species are not all registered gets n_invariants = 0
// and no projection, i.e. exactly the behaviour without this file.  Measured
// n_invariants: 5 primordial (full and Minimal), 8 Nakauchi2021_Minimal (the
// compact metal network carries no K or Na, so those two rows are pruned), 10
// Nakauchi2021.  (2) The implicit chemistry solve only: the projection runs
// before the operator-split Lyman-Werner / X-ray channels of chem_full_step(),
// whose stoichiometry is outside the reaction table this file validates against.
//
// The projection is NOT result-neutral on metal_grain: measured at zeta = 0,
// Z = 1, 24 of 33 datasets differ from an unprojected run (attributes 24/24
// unchanged).  On the primordial networks it is result-neutral: 25/25
// datasets byte-equal.
//
// MEASURED RANGE.  PRECONDITIONS: J_LW21 == 0 and ARCHE_XRAY off.  Under those,
// element residuals stay at the float64 representation floor of the element
// total (0 .. 4.4e-16), and the net charge residual at (0 .. 6.0e-16), over
// nH = 0.1 .. 1e23 cm^-3 on both primordial networks, at zeta = 0 and
// zeta = 1e-17, under g++ 13.3.0, clang 22.1.6 and clang 22.1.8.  The
// same runs reach 7.6e-11 (H) and 1.06e-6 (net charge) without it.  No density
// was found
// at which the correction stops holding.
//   The number quoted is a property of the DIAGNOSTIC, not of the solver: this
// file never forms an element total in float64 (it accumulates increments plus
// carry in long double), so the floor depends on how the checker sums.  The same
// projected data gives 1.1e-16 (fsum) .. 4.4e-16 (naive left-to-right) for the
// H
// row.  4.4e-16 is the loosest convention and is therefore the value to gate on
// -- but name the convention when you do.
//
// WITH J_LW21 > 0 THE ELEMENT FLOOR IS 2.4 DECADES HIGHER.  Measured (clang
// 22.1.8, 4 runs): H reaches 1.5e-13 and D 5.5e-14, while He, Li and the
// net charge stay at the floor quoted above.  That split is the diagnosis: no
// Lyman-Werner channel touches He or Li, and the charge row's target is absolute
// and re-derived every step.  The cause is the placement recorded under SCOPE --
// the LW block of chem_full_step() runs AFTER this projection, so its own
// rounding is folded into the next step's y_in and is never corrected.  The
// effect is essentially independent of the field strength (factor 1.32 across 11
// decades of J_LW21, worst at 1e-8, not at the 1e3 of README Quick Start Case 3)
// and it SATURATES rather than accumulating: a sqrt(N_steps) walk, so reaching
// 1e-12 would take ~4e7 steps against the 5.6e4 of the longest run measured.
// Without the projection the same runs reach 4.6e-11 (H) and 6.3e-09 (D) -- it
// is still 3.6e2x (H) and 1.1e5x (D) better with the field on.
//
// The X-ray block is the second uncorrected operator-split channel and is
// UNMEASURED: no run anywhere has had J_LW21 > 0 with ARCHE_XRAY compiled in.
//
// KNOWN LIMITS.
//  * Above nH_eq = 1e18 cm^-3 with zeta == 0 and zeta_X == 0, chemreact() takes
//    a Saha equilibrium branch that RE-ANCHORS the He, D and Li totals to
//    hard-coded reference abundances, while this projection's target is the
//    totals the step was handed.  At J_LW21 == 0 the two agree to float64 at the
//    crossing and do not fight.  They DO fight as soon as anything else perturbs
//    the totals, and then the projection wins every step.  Measured (zeta = 0,
//    J_LW21 = 1e3, 764 output rows above 1e18): without the projection the
//    branch returns the D total to yD = 2.58e-5 every step (max deviation
//    2.0 ulp, 338/764 rows exactly 0); with it, D sits 4.0e-14 (307 ulp) below
//    the reference and never returns to it, frozen from nH = 2.5e19 to the
//    last row.  Over the whole run the projected D residual is still 7.2e4x
//    smaller, and 4e-14 relative on yD is ~12 decades below the BBN
//    uncertainty on D/H: a semantic change, not a physical error.  Note the
//    branch does not leave H alone either -- it drives the H total to 1.0 by 2D
//    Newton, but only to Newton tolerance (without the projection,
//    |H_tot - 1| = 6.6e-11 here).
//    A caller entering that branch with a composition away from the reference
//    values has not been tested.
//  * On metal_grain the singular-matrix branch fires on the MAJORITY of steps
//    unless the ACTIVE SET below is applied: the network starts with all of
//    its K and Na in
//    the application-side grain reservoir, so both rows were empty and the
//    joint solve declined for 3566 of 5537 output rows (64%, zeta = 1e-17,
//    Z = 1) and 4034 of 5612 (72%, zeta = 0, Z = 1e-3).  The ACTIVE SET below
//    is the fix: an empty row with nothing owed is dropped, and the remaining
//    rows are still corrected.  Both remaining decline paths -- kMaxRelShift
//    exceeded, and an empty row that IS owed a correction -- remain designed
//    behaviour rather than demonstrated behaviour.
//  * Conditioning after Jacobi scaling, measured on the 10-row metal system
//    (zeta = 1e-17, Z = 1, 1971 engaged output rows): cond(G) min 1.117,
//    median 6.37, max 6.725 at nH = 3.7e15.  Nine decades of element-total
//    spread do not threaten the Cholesky; rank deficiency does, which is what
//    the active set handles.
//  * The projection changes abundances only, but at zeta != 0 those feed the
//    ionisation terms and hence the adaptive timestep, so the density schedule
//    itself moves by up to 1.2e-14 relative.  At zeta == 0 the nH grid is
//    bit-identical.
//  * The long double accumulations below are this codebase's only arithmetic
//    whose width depends on the ABI (80-bit extended on x86-64 SysV, but 64-bit
//    on some others).  Deliberately NOT guarded by a static_assert: failing the
//    build on such a target is worse than the degradation, which is confined to
//    the charge row (the element rows sum increments, not totals, and provably
//    do not need the extra width).  If you port to an ABI where
//    LDBL_MANT_DIG == DBL_MANT_DIG, re-measure the charge residual before
//    quoting the range above.
//  * carry lives in ChemCell, not in ChemState, so the documented copy-in /
//    copy-out cell API can neither see nor reset it.  A caller that reuses one
//    ChemCell as scratch across many fluid cells will replay a foreign
//    remainder into the next cell's first step; call ChemCell::reset_var()
//    between cells.  The error this bounds is <= half an ulp of the element
//    total (derived, not measured).
//  * build_invariants() balances all three stoichiometry tables: tbl.reactions,
//    tbl.saha (n_saha = 18 full / 9 minimal for primordial, 53 for
//    metal_grain) and tbl.grain_reactions (148 for metal_grain; primordial has
//    n_grain = 0).  saha drives the branch that rewrites y above nH_eq and
//    grain_reactions moves species between the gas and the mantle, so a typo
//    in either is as consequential as one in reactions.  Measured: all 18
//    primordial saha records and all 53 + 148 metal records balance in nine
//    elements plus charge, so the check fails no shipped network closed.
//
// ---------------------------------------------------------------------------
#include <array>
#include <cmath>
#include <cstddef>

#include "core/species_composition.h"
#include "kinetics/topology.h"  // Reaction, ReactionTable

namespace arche {
namespace conservation {

//: Largest relative per-species shift the projection will apply.  Matched to
//: the 1e-2 element/charge residual bound the collapse app is validated
//: against: a violation that would need a larger repair is a defect to be seen,
//: not to be papered over, so the projection declines and leaves both the state
//: and the evidence untouched.
inline constexpr double kMaxRelShift = 1.0e-2;

//: Upper bound on the number of invariant rows: one per tracked element plus
//: charge.  Sizes the per-cell carry (ChemCell::cons_carry) and the scratch
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
  constexpr int N_rows = composition::kNElements + 1;  // elements + charge

  // Per-species integer composition, in the model's local index order.
  std::array<std::array<int, N_rows>, N_sp> comp{};
  for (int i = 0; i < N_sp; ++i) {
    const composition::Composition c =
        composition::composition_of(SpeciesSetT::canonical(i));
    if (!c.known) return 0;
    for (int e = 0; e < composition::kNElements; ++e) comp[i][e] = c.n[e];
    comp[i][composition::kNElements] = c.q;
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
    if (e == composition::kNElements) charge_row = m;
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
bool project(const double* inv, int m, int charge_row,
             const std::array<double, N_sp>& y_in, std::array<double, N_sp>& y,
             std::array<double, kMaxRows>& carry) {
  if (m <= 0 || m > N_sp || m > kMaxRows) return false;

  // Targets.  Element rows are formed from the per-step INCREMENTS and never as
  // the difference of two element totals: the increments are ~1e9 times smaller
  // than the totals here, so this expression carries the residual at full
  // precision while differencing two float64 totals would return their
  // summation noise instead.  The charge row is absolute (neutrality), so it
  // has no such shortcut and is accumulated in long double to keep the
  // cancellation between the electron and ion terms from setting the floor.
  std::array<double, kMaxRows> r{};
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
  // element absent from y -- puts a zero on G's diagonal.  Declining the whole
  // projection there is not a rare corner: the full metal network starts with
  // ALL of its K and Na in the application-side grain reservoir (yKs / yNas in
  // apps/collapse_metal_grain/main.cc), so y[K] = y[K+] = y[Na] = y[Na+] = 0
  // and both rows are empty until evaporation releases them.  Measured on
  // zeta = 1e-17, Z = 1: that is 3566 of 5537 output rows, 64% of the collapse,
  // during which the H / He / D / Li rows would be discarded too -- leaving the
  // projection configured (n_invariants = 10) and silently inert.
  //   gkk == 0 and r == 0 -> nothing to conserve and nothing owed; drop the row
  //                          and solve the rest.
  //   gkk == 0 and r != 0 -> the step changed an element that no surviving
  //                          species carries.  Unrepairable: decline.
  // A diagonal regulariser is NOT the fix.  It would manufacture a correction
  // and spread it over species that cannot carry the element at all.
  std::array<int, kMaxRows> act{};
  int n_act = 0;
  for (int k = 0; k < m; ++k) {
    if (g[k * kLd + k] > 0.0) {
      act[n_act++] = k;
    } else if (r[k] != 0.0) {
      return false;
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
