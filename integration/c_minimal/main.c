/* Copyright (C) 2026 Shingo Hirano and Sho Higashi
 * Licensed under the MIT found in the
 * https://github.com/astro-sim-lab/arche/blob/main/LICENSE
 *
 * main.c — minimal external C consumer of libarche.
 *
 * STANDALONE project compiled as C: includes ONLY <api/arche_capi.h> and links
 * arche::arche. Exceptions never cross the boundary; fallible calls return a
 * status code and arche_last_error() gives the message.
 *
 * Usage:  c_minimal [primordial.h5]   (defaults to ./data/primordial.h5)
 */

#include <stdio.h>

#include "api/arche_capi.h"

int main(int argc, char** argv) {
  const char* table = (argc > 1) ? argv[1] : "data/primordial.h5";

  ArchePrimTable* tbl = arche_prim_table_load(table);
  if (!tbl) {
    fprintf(stderr, "c_minimal: table load failed: %s\n", arche_last_error());
    return 1;
  }
  ArchePrimCell* cell = arche_prim_cell_create();
  if (!cell) {
    fprintf(stderr, "c_minimal: cell allocation failed\n");
    arche_prim_table_free(tbl);
    return 1;
  }

  ArchePrimState st;
  for (int i = 0; i < ARCHE_PRIM_NSP; ++i) st.y[i] = 0.0;
  st.y[0] = 1.0 - 1.0e-4; /* H  */
  st.y[2] = 1.0e-4;       /* e- */
  st.y[3] = 1.0e-4;       /* H+ */
  st.y[7] = 8.33e-2;      /* He (yHe) */
  st.nH = 1.0e4;
  st.T_K = 200.0;
  st.mu = 1.22;
  st.gamma = 5.0 / 3.0;
  arche_prim_cell_set_state(cell, &st);

  /* Initialise to the library defaults (NOT the same as = {0}: esc_cnt
   * defaults to 1.0 = optically thin), then set what you use. */
  ArcheChemParams params;
  arche_chem_params_init(&params);
  params.T_rad = 2.725;
  ArcheChemShielding shield;
  arche_chem_shielding_init(&shield);
  shield.zeta = 1.0e-17;

  ArcheChemFullRates r;
  int rc = arche_chem_full_step_prim(cell, 1.0e8, &params, &shield, tbl, &r);
  if (rc != ARCHE_OK) {
    fprintf(stderr, "c_minimal: step error %d: %s\n", rc, arche_last_error());
    arche_prim_cell_free(cell);
    arche_prim_table_free(tbl);
    return 1;
  }
  printf("[arche c_minimal] L_net=%.6e  Gam_CR=%.6e  solver_failed=%d\n",
         r.Lambda_net, r.Gamma_CR, r.solver_failed);

  arche_prim_cell_free(cell);
  arche_prim_table_free(tbl);
  return 0;
}
