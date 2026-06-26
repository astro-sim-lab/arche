// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
//
// main.cc — minimal external C++ consumer of libarche.
//
// This is a STANDALONE project: it is built on its own (not as part of the
// arche tree) and finds the installed library with find_package(arche). It
// includes ONLY <api/arche_api.h> and links arche::arche — no Eigen, no HDF5 in
// this source.
//
// Usage:  cpp_minimal [primordial.h5]   (defaults to ./data/primordial.h5)

#include <cstdio>
#include <string>

#include "api/arche_api.h"

namespace zm = arche::zero_metal;

int main(int argc, char** argv) {
  const std::string table = (argc > 1) ? argv[1] : "data/primordial.h5";

  try {
    arche::PrimTablePtr tbl = arche::load_prim_table_owned(table);
    arche::PrimCellPtr cell = arche::prim_cell_create_owned();

    arche::ChemStateZM& s = arche::prim_cell_state(*cell);
    s.nH = 1.0e4;
    s.T_K = 200.0;
    s.mu = 1.22;
    s.gamma = 5.0 / 3.0;
    // Index y[] by species name (the order is model-specific; see
    // core/species_index.h).
    s.y.fill(0.0);
    s.y[zm::H] = 1.0 - 1.0e-4;               // atomic H (~all of nH)
    s.y[zm::e] = 1.0e-4;                     // electrons
    s.y[zm::Hp] = 1.0e-4;                    // H+ (balances e-)
    s.y[zm::He] = arche::abundance_ref::yHe;  // He, n(He)/nH

    arche::ChemParams params{};
    params.T_rad = 2.725;
    arche::ChemShielding shield{};
    shield.zeta = 1.0e-17;

    const arche::ChemFullRates r =
        arche::chem_full_step_prim(*cell, /*dt=*/1.0e8, params, shield, *tbl);

    std::printf(
        "[arche cpp_minimal] L_net=%.6e  Gam_CR=%.6e  solver_failed=%d\n",
        r.Lambda_net, r.Gamma_CR, static_cast<int>(r.solver_failed));
  } catch (const std::exception& e) {
    std::fprintf(stderr, "cpp_minimal: %s\n", e.what());
    return 1;
  }
  return 0;
}
