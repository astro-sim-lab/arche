// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
// ---------------------------------------------------------------------------
// arche_capi.cc — implementation of the C ABI.
//
// Wraps the C++ facade (arche_api.h): marshals the plain-C mirror structs to
// and from the arche:: POD structs, reinterprets the opaque C handle pointers
// as the corresponding arche:: handle pointers, and catches every exception at
// the boundary, translating it to a status code plus a thread-local message.
// ---------------------------------------------------------------------------
#include "api/arche_capi.h"

#include <stdexcept>
#include <string>

#include "api/arche_api.h"

namespace {

// Per-thread last-error message (returned by arche_last_error()).
thread_local std::string g_last_error;
void set_last_error(const char* msg) { g_last_error = msg ? msg : ""; }

// ── Struct marshalling (field-by-field; order-independent, layout-agnostic) ──
arche::ChemParams to_cpp(const ArcheChemParams& c) {
  arche::ChemParams p;
  p.zeta = c.zeta;
  p.zeta_X = c.zeta_X;
  p.T_rad = c.T_rad;
  p.T_gr_K = c.T_gr_K;
  p.Z_metal = c.Z_metal;
  p.T_cr_desorp = c.T_cr_desorp;
  p.H = c.H;
  p.H2 = c.H2;
  p.He = c.He;
  p.J_H2 = c.J_H2;
  p.J_H2O = c.J_H2O;
  p.J_tot = c.J_tot;
  return p;
}

arche::ChemShielding to_cpp(const ArcheChemShielding& c) {
  arche::ChemShielding s;
  s.zeta = c.zeta;
  s.Nc_H = c.Nc_H;
  s.Nc_H2 = c.Nc_H2;
  s.Nc_HD = c.Nc_HD;
  s.tau_cnt = c.tau_cnt;
  s.esc_cnt = c.esc_cnt;
  s.J_LW21 = c.J_LW21;
  s.Nc_CO = c.Nc_CO;
  s.Nc_OH = c.Nc_OH;
  s.Nc_H2O = c.Nc_H2O;
  s.Nc_CII = c.Nc_CII;
  s.Nc_CI = c.Nc_CI;
  s.Nc_OI = c.Nc_OI;
  s.zeta_X = c.zeta_X;
  s.E_X_eV = c.E_X_eV;
  return s;
}

void from_cpp(const arche::ChemFullRates& r, ArcheChemFullRates& c) {
  c.Lambda_net = r.Lambda_net;
  c.Lambda_line = r.Lambda_line;
  c.Lambda_cnt = r.Lambda_cnt;
  c.Lambda_chem = r.Lambda_chem;
  c.Gamma_CR = r.Gamma_CR;
  c.Gamma_X = r.Gamma_X;
  c.Lambda_H2 = r.Lambda_H2;
  c.Lambda_HD = r.Lambda_HD;
  c.Lambda_Lya = r.Lambda_Lya;
  c.Lambda_CO = r.Lambda_CO;
  c.Lambda_OH = r.Lambda_OH;
  c.Lambda_H2O = r.Lambda_H2O;
  c.Lambda_CII = r.Lambda_CII;
  c.Lambda_CI = r.Lambda_CI;
  c.Lambda_OI = r.Lambda_OI;
  c.Lambda_gr = r.Lambda_gr;
  c.Lambda_gas = r.Lambda_gas;
  c.k_gas = r.k_gas;
  c.k_gr = r.k_gr;
  c.T_gr_K = r.T_gr_K;
  c.solver_failed = r.solver_failed ? 1 : 0;
}

// Map a caught exception to an ARCHE_ERR_* code (and record its message).
template <class Fn>
int guarded(Fn&& fn) {
  try {
    fn();
    return ARCHE_OK;
  } catch (const std::invalid_argument& e) {
    set_last_error(e.what());
    return ARCHE_ERR_INVALID_ARGUMENT;
  } catch (const std::runtime_error& e) {
    set_last_error(e.what());
    return ARCHE_ERR_RUNTIME;
  } catch (const std::exception& e) {
    set_last_error(e.what());
    return ARCHE_ERR_UNKNOWN;
  } catch (...) {
    set_last_error("unknown error");
    return ARCHE_ERR_UNKNOWN;
  }
}

}  // namespace

extern "C" {

// ── Struct initialisers ─────────────────────────────────────────────────────
// Copy the C++ defaults field-by-field so the C mirrors can never drift from
// core/state.h.
void arche_chem_params_init(ArcheChemParams* params) {
  if (!params) return;
  const arche::ChemParams d;
  params->zeta = d.zeta;
  params->zeta_X = d.zeta_X;
  params->T_rad = d.T_rad;
  params->T_gr_K = d.T_gr_K;
  params->Z_metal = d.Z_metal;
  params->T_cr_desorp = d.T_cr_desorp;
  params->H = d.H;
  params->H2 = d.H2;
  params->He = d.He;
  params->J_H2 = d.J_H2;
  params->J_H2O = d.J_H2O;
  params->J_tot = d.J_tot;
}
void arche_chem_shielding_init(ArcheChemShielding* shield) {
  if (!shield) return;
  const arche::ChemShielding d;
  shield->zeta = d.zeta;
  shield->Nc_H = d.Nc_H;
  shield->Nc_H2 = d.Nc_H2;
  shield->Nc_HD = d.Nc_HD;
  shield->tau_cnt = d.tau_cnt;
  shield->esc_cnt = d.esc_cnt;
  shield->J_LW21 = d.J_LW21;
  shield->Nc_CO = d.Nc_CO;
  shield->Nc_OH = d.Nc_OH;
  shield->Nc_H2O = d.Nc_H2O;
  shield->Nc_CII = d.Nc_CII;
  shield->Nc_CI = d.Nc_CI;
  shield->Nc_OI = d.Nc_OI;
  shield->zeta_X = d.zeta_X;
  shield->E_X_eV = d.E_X_eV;
}

// ── Cell lifecycle ─────────────────────────────────────────────────────────
ArchePrimCell* arche_prim_cell_create(void) {
  return reinterpret_cast<ArchePrimCell*>(arche::prim_cell_create());
}
void arche_prim_cell_free(ArchePrimCell* cell) {
  arche::prim_cell_destroy(reinterpret_cast<arche::PrimCell*>(cell));
}
ArcheMetalCell* arche_metal_cell_create(void) {
  return reinterpret_cast<ArcheMetalCell*>(arche::metal_cell_create());
}
void arche_metal_cell_free(ArcheMetalCell* cell) {
  arche::metal_cell_destroy(reinterpret_cast<arche::MetalCell*>(cell));
}

// ── Per-cell state ─────────────────────────────────────────────────────────
void arche_prim_cell_set_state(ArchePrimCell* cell, const ArchePrimState* st) {
  if (!cell || !st) return;
  arche::ChemStateZM& s =
      arche::prim_cell_state(*reinterpret_cast<arche::PrimCell*>(cell));
  for (int i = 0; i < ARCHE_PRIM_NSP; ++i) s.y[i] = st->y[i];
  s.nH = st->nH;
  s.T_K = st->T_K;
  s.mu = st->mu;
  s.gamma = st->gamma;
}
void arche_prim_cell_get_state(const ArchePrimCell* cell, ArchePrimState* st) {
  if (!cell || !st) return;
  const arche::ChemStateZM& s =
      arche::prim_cell_state(*reinterpret_cast<const arche::PrimCell*>(cell));
  for (int i = 0; i < ARCHE_PRIM_NSP; ++i) st->y[i] = s.y[i];
  st->nH = s.nH;
  st->T_K = s.T_K;
  st->mu = s.mu;
  st->gamma = s.gamma;
}
void arche_metal_cell_set_state(ArcheMetalCell* cell,
                                const ArcheMetalState* st) {
  if (!cell || !st) return;
  arche::ChemStateMG& s =
      arche::metal_cell_state(*reinterpret_cast<arche::MetalCell*>(cell));
  for (int i = 0; i < ARCHE_METAL_NSP; ++i) s.y[i] = st->y[i];
  s.nH = st->nH;
  s.T_K = st->T_K;
  s.mu = st->mu;
  s.gamma = st->gamma;
}
void arche_metal_cell_get_state(const ArcheMetalCell* cell,
                                ArcheMetalState* st) {
  if (!cell || !st) return;
  const arche::ChemStateMG& s =
      arche::metal_cell_state(*reinterpret_cast<const arche::MetalCell*>(cell));
  for (int i = 0; i < ARCHE_METAL_NSP; ++i) st->y[i] = s.y[i];
  st->nH = s.nH;
  st->T_K = s.T_K;
  st->mu = s.mu;
  st->gamma = s.gamma;
}

// ── Table loaders ──────────────────────────────────────────────────────────
ArchePrimTable* arche_prim_table_load(const char* h5_path) {
  if (!h5_path) {
    set_last_error("arche_prim_table_load: h5_path is NULL");
    return nullptr;
  }
  ArchePrimTable* out = nullptr;
  guarded([&] {
    out = reinterpret_cast<ArchePrimTable*>(arche::load_prim_table(h5_path));
  });
  return out;
}
void arche_prim_table_free(ArchePrimTable* tbl) {
  arche::prim_table_destroy(reinterpret_cast<arche::PrimTable*>(tbl));
}
ArcheMetalTable* arche_metal_table_load(const char* h5_path) {
  if (!h5_path) {
    set_last_error("arche_metal_table_load: h5_path is NULL");
    return nullptr;
  }
  ArcheMetalTable* out = nullptr;
  guarded([&] {
    out = reinterpret_cast<ArcheMetalTable*>(arche::load_metal_table(h5_path));
  });
  return out;
}
void arche_metal_table_free(ArcheMetalTable* tbl) {
  arche::metal_table_destroy(reinterpret_cast<arche::MetalTable*>(tbl));
}

// ── Stepping ───────────────────────────────────────────────────────────────
int arche_chem_full_step_prim(ArchePrimCell* cell, double dt,
                              const ArcheChemParams* params,
                              const ArcheChemShielding* shield,
                              const ArchePrimTable* tbl,
                              ArcheChemFullRates* out) {
  if (!cell || !params || !shield || !tbl) return ARCHE_ERR_NULL;
  return guarded([&] {
    arche::ChemFullRates r = arche::chem_full_step_prim(
        *reinterpret_cast<arche::PrimCell*>(cell), dt, to_cpp(*params),
        to_cpp(*shield), *reinterpret_cast<const arche::PrimTable*>(tbl));
    if (out) from_cpp(r, *out);
  });
}
int arche_chem_full_step_metal(ArcheMetalCell* cell, double dt,
                               const ArcheChemParams* params,
                               const ArcheChemShielding* shield,
                               const ArcheMetalTable* tbl,
                               ArcheChemFullRates* out) {
  if (!cell || !params || !shield || !tbl) return ARCHE_ERR_NULL;
  return guarded([&] {
    arche::ChemFullRates r = arche::chem_full_step_metal(
        *reinterpret_cast<arche::MetalCell*>(cell), dt, to_cpp(*params),
        to_cpp(*shield), *reinterpret_cast<const arche::MetalTable*>(tbl));
    if (out) from_cpp(r, *out);
  });
}

// ── Runtime model dispatch ──────────────────────────────────────────────────
// ArcheCell/ArcheTable map 1:1 onto arche::Cell/arche::Table (the C++ tagged
// dispatch handles); reinterpret the opaque pointers and delegate.
ArcheCell* arche_cell_create(ArcheModel model) {
  arche::Model m =
      (model == ARCHE_MODEL_METAL) ? arche::Model::Metal : arche::Model::Prim;
  return reinterpret_cast<ArcheCell*>(arche::cell_create(m));
}
void arche_cell_free(ArcheCell* cell) {
  arche::cell_destroy(reinterpret_cast<arche::Cell*>(cell));
}
ArcheTable* arche_table_load(ArcheModel model, const char* h5_path) {
  if (!h5_path) {
    set_last_error("arche_table_load: h5_path is NULL");
    return nullptr;
  }
  arche::Model m =
      (model == ARCHE_MODEL_METAL) ? arche::Model::Metal : arche::Model::Prim;
  ArcheTable* out = nullptr;
  guarded([&] {
    out = reinterpret_cast<ArcheTable*>(arche::load_table(m, h5_path));
  });
  return out;
}
void arche_table_free(ArcheTable* tbl) {
  arche::table_destroy(reinterpret_cast<arche::Table*>(tbl));
}

ArcheModel arche_cell_model(const ArcheCell* cell) {
  return arche::cell_model(*reinterpret_cast<const arche::Cell*>(cell)) ==
                 arche::Model::Metal
             ? ARCHE_MODEL_METAL
             : ARCHE_MODEL_PRIM;
}
int arche_cell_nsp(const ArcheCell* cell) {
  return arche::cell_n_sp(*reinterpret_cast<const arche::Cell*>(cell));
}
double* arche_cell_y(ArcheCell* cell) {
  return arche::cell_y(*reinterpret_cast<arche::Cell*>(cell));
}
void arche_cell_set_scalars(ArcheCell* cell, double nH, double T_K, double mu,
                            double gamma) {
  arche::Cell& c = *reinterpret_cast<arche::Cell*>(cell);
  arche::cell_nH(c) = nH;
  arche::cell_T_K(c) = T_K;
  arche::cell_mu(c) = mu;
  arche::cell_gamma(c) = gamma;
}
void arche_cell_get_scalars(const ArcheCell* cell, double* nH, double* T_K,
                            double* mu, double* gamma) {
  // cell_* accessors are non-const refs; cast away constness for read-only use.
  arche::Cell& c =
      *const_cast<arche::Cell*>(reinterpret_cast<const arche::Cell*>(cell));
  if (nH) *nH = arche::cell_nH(c);
  if (T_K) *T_K = arche::cell_T_K(c);
  if (mu) *mu = arche::cell_mu(c);
  if (gamma) *gamma = arche::cell_gamma(c);
}
int arche_chem_full_step(ArcheCell* cell, double dt,
                         const ArcheChemParams* params,
                         const ArcheChemShielding* shield,
                         const ArcheTable* tbl, ArcheChemFullRates* out) {
  if (!cell || !params || !shield || !tbl) return ARCHE_ERR_NULL;
  return guarded([&] {
    arche::ChemFullRates r = arche::chem_full_step(
        *reinterpret_cast<arche::Cell*>(cell), dt, to_cpp(*params),
        to_cpp(*shield), *reinterpret_cast<const arche::Table*>(tbl));
    if (out) from_cpp(r, *out);
  });
}

// ── Registry-backed, name-selected models ───────────────────────────────────
// ArcheModelRuntime/ArcheModelCell map 1:1 onto arche::ModelRuntime/ModelCell;
// reinterpret the opaque pointers and delegate to the C++ facade.
ArcheModelRuntime* arche_model_create(const char* name, const char* h5_path) {
  if (!name || !h5_path) {
    set_last_error("arche_model_create: name or h5_path is NULL");
    return nullptr;
  }
  ArcheModelRuntime* out = nullptr;
  guarded([&] {
    out = reinterpret_cast<ArcheModelRuntime*>(
        arche::model_create(name, h5_path));
  });
  return out;
}
void arche_model_destroy(ArcheModelRuntime* model) {
  arche::model_destroy(reinterpret_cast<arche::ModelRuntime*>(model));
}
int arche_model_n_species(const ArcheModelRuntime* model) {
  if (!model) return -1;
  return arche::model_n_species(
      *reinterpret_cast<const arche::ModelRuntime*>(model));
}
const char* const* arche_model_species(const ArcheModelRuntime* model) {
  if (!model) return nullptr;
  return arche::model_species(
      *reinterpret_cast<const arche::ModelRuntime*>(model));
}
const char* arche_model_name(const ArcheModelRuntime* model) {
  if (!model) return nullptr;
  return arche::model_name(
      *reinterpret_cast<const arche::ModelRuntime*>(model));
}
int arche_model_count(void) { return arche::model_count(); }
const char* arche_model_registry_name(int i) {
  return arche::model_registry_name(i);
}
ArcheModelCell* arche_model_cell_create(const ArcheModelRuntime* model) {
  if (!model) {
    set_last_error("arche_model_cell_create: model is NULL");
    return nullptr;
  }
  ArcheModelCell* out = nullptr;
  guarded([&] {
    out = reinterpret_cast<ArcheModelCell*>(arche::model_cell_create(
        *reinterpret_cast<const arche::ModelRuntime*>(model)));
  });
  return out;
}
void arche_model_cell_free(ArcheModelCell* cell) {
  arche::model_cell_destroy(reinterpret_cast<arche::ModelCell*>(cell));
}
double* arche_model_cell_y(ArcheModelCell* cell) {
  if (!cell) return nullptr;
  return arche::model_cell_y(*reinterpret_cast<arche::ModelCell*>(cell));
}
void arche_model_cell_set_scalars(ArcheModelCell* cell, double nH, double T_K,
                                  double mu, double gamma) {
  if (!cell) return;
  arche::ModelCell& c = *reinterpret_cast<arche::ModelCell*>(cell);
  arche::model_cell_nH(c) = nH;
  arche::model_cell_T_K(c) = T_K;
  arche::model_cell_mu(c) = mu;
  arche::model_cell_gamma(c) = gamma;
}
void arche_model_cell_get_scalars(const ArcheModelCell* cell, double* nH,
                                  double* T_K, double* mu, double* gamma) {
  if (!cell) return;
  arche::ModelCell& c = *const_cast<arche::ModelCell*>(
      reinterpret_cast<const arche::ModelCell*>(cell));
  if (nH) *nH = arche::model_cell_nH(c);
  if (T_K) *T_K = arche::model_cell_T_K(c);
  if (mu) *mu = arche::model_cell_mu(c);
  if (gamma) *gamma = arche::model_cell_gamma(c);
}
int arche_model_step(const ArcheModelRuntime* model, ArcheModelCell* cell,
                     double dt, const ArcheChemParams* params,
                     const ArcheChemShielding* shield,
                     ArcheChemFullRates* out) {
  if (!model || !cell || !params || !shield) return ARCHE_ERR_NULL;
  return guarded([&] {
    arche::ChemFullRates r =
        arche::model_step(*reinterpret_cast<const arche::ModelRuntime*>(model),
                          *reinterpret_cast<arche::ModelCell*>(cell), dt,
                          to_cpp(*params), to_cpp(*shield));
    if (out) from_cpp(r, *out);
  });
}

const char* arche_last_error(void) { return g_last_error.c_str(); }

}  // extern "C"
