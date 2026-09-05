// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#ifndef SRC_INCLUDE_API_ARCHE_CAPI_H_
#define SRC_INCLUDE_API_ARCHE_CAPI_H_
// ---------------------------------------------------------------------------
// arche_capi.h — C ABI for libarche
//
// A thin `extern "C"` skin over the C++ facade (arche_api.h) so the chemistry
// network can be driven from C, Fortran (ISO_C_BINDING), Python (ctypes/cffi),
// etc.  Pure C: this header parses under a C compiler and contains no C++
// types.  Compared to arche_api.h it differs in two ways suited to a stable C
// boundary:
//
//   * Exceptions do not cross the boundary.  Functions that can fail return an
//     int status (ARCHE_OK or an ARCHE_ERR_* code); the human-readable message
//     for the calling thread is available from arche_last_error().
//   * The C++ POD structs are mirrored as plain C structs (same fields) and
//     marshalled by value across the boundary; the sized kernel types stay
//     opaque (handle pointers), exactly as in arche_api.h.
//
// All numeric conventions (units, meaning of each field) match arche_api.h /
// docs/api_reference.md.  Output dataset schema is unrelated (see
// docs/output_schema.md).
// ---------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif

// Network sizes (length of the species vector for each model).
#define ARCHE_PRIM_NSP 23
#define ARCHE_METAL_NSP 89

// Status codes returned by fallible functions.
enum {
  ARCHE_OK = 0,
  ARCHE_ERR_INVALID_ARGUMENT = 1,  // std::invalid_argument at the boundary
  ARCHE_ERR_RUNTIME = 2,           // std::runtime_error (e.g. table load)
  ARCHE_ERR_UNKNOWN = 3,           // any other exception
  ARCHE_ERR_NULL = 4               // a required pointer argument was NULL
};

// ── Opaque handles ─────────────────────────────────────────────────────────
// Definitions live inside libarche; callers only hold/pass the pointers.
typedef struct ArchePrimCell ArchePrimCell;      // primordial (Nakauchi2019)
typedef struct ArcheMetalCell ArcheMetalCell;    // metal+grain (Nakauchi2021)
typedef struct ArchePrimTable ArchePrimTable;    // primordial reaction table
typedef struct ArcheMetalTable ArcheMetalTable;  // metal+grain reaction table

// ── POD mirror structs (mirror arche::* fields; see arche_api.h/state.h) ────
typedef struct {
  double zeta;         // cosmic-ray ionization rate [s^-1] (pre-attenuated)
  double zeta_X;       // X-ray secondary ionization rate [s^-1]
  double T_rad;        // CMB radiation temperature [K]  (must be finite > 0)
  double T_gr_K;       // grain temperature warm-start [K] (metal; must be >0)
  double Z_metal;      // metallicity [Z_sun] (metal)
  double T_cr_desorp;  // effective CR desorption spike temperature [K]
  double H, H2, He;    // LTE-interpolation abundances (filled by kernel)
  double J_H2, J_H2O, J_tot;  // UV radiation field strengths (metal)
} ArcheChemParams;

typedef struct {
  double zeta;                // effective CR ionization rate [s^-1]
  double Nc_H, Nc_H2, Nc_HD;  // column densities [cm^-2]
  double tau_cnt;             // continuum optical depth
  double esc_cnt;             // continuum escape fraction (1 = optically thin)
  double J_LW21;  // Lyman-Werner intensity [10^-21 erg/s/cm^2/Hz/sr]
  double Nc_CO, Nc_OH, Nc_H2O, Nc_CII, Nc_CI, Nc_OI;  // metal column densities
  double zeta_X;  // X-ray HI primary photoionization rate [s^-1]
  double E_X_eV;  // representative X-ray photon energy [eV]
} ArcheChemShielding;

typedef struct {
  double Lambda_net, Lambda_line, Lambda_cnt, Lambda_chem;  // aggregate cooling
  double Gamma_CR, Gamma_X;                                 // heating
  double Lambda_H2, Lambda_HD, Lambda_Lya;                  // per-line
  double Lambda_CO, Lambda_OH, Lambda_H2O;                  // per-line (metal)
  double Lambda_CII, Lambda_CI, Lambda_OI;                  // per-line (metal)
  double Lambda_gr, Lambda_gas;                             // continuum split
  double k_gas, k_gr;                                       // opacities
  double T_gr_K;                                            // solved grain temp
  int solver_failed;  // 0 = ok, 1 = failed
  // 0 = element/charge conservation was NOT enforced on this step, 1 = it was.
  // See ChemFullRates::conservation_projected in core/state.h for the four
  // reasons a 0 appears.  All four shipped networks register invariant rows
  // (5 / 5 / 8 / 10), so the "no rows registered" reason does not apply to
  // them: on a step that converged, a 0 means the projection itself declined.
  int conservation_projected;
} ArcheChemFullRates;

// ── Struct initialisers ─────────────────────────────────────────────────────
// Set every field to the library default (matching the C++ defaults declared
// in core/state.h).  Plain zero-initialisation (`= {0}`) is NOT equivalent:
// the defaults include esc_cnt = 1.0 (optically thin), E_X_eV = 300,
// T_cr_desorp = 70, and T_rad = NaN (a sentinel that forces the caller to set
// it).  Call the initialiser first, then overwrite the fields you use.
void arche_chem_params_init(ArcheChemParams* params);
void arche_chem_shielding_init(ArcheChemShielding* shield);

// Per-cell state the caller fills each step (mirrors arche::ChemState<N_sp>).
typedef struct {
  double y[ARCHE_PRIM_NSP];  // species abundances (relative to nH)
  double nH, T_K, mu, gamma;
} ArchePrimState;

typedef struct {
  double y[ARCHE_METAL_NSP];
  double nH, T_K, mu, gamma;
} ArcheMetalState;

// ── Cell lifecycle ─────────────────────────────────────────────────────────
// create returns NULL on allocation failure. Cells are NOT thread-safe (one
// per thread); tables are read-only and shareable.
ArchePrimCell* arche_prim_cell_create(void);
void arche_prim_cell_free(ArchePrimCell* cell);
ArcheMetalCell* arche_metal_cell_create(void);
void arche_metal_cell_free(ArcheMetalCell* cell);

// ── Per-cell state (copy in / copy out) ────────────────────────────────────
void arche_prim_cell_set_state(ArchePrimCell* cell, const ArchePrimState* st);
void arche_prim_cell_get_state(const ArchePrimCell* cell, ArchePrimState* st);
void arche_metal_cell_set_state(ArcheMetalCell* cell,
                                const ArcheMetalState* st);
void arche_metal_cell_get_state(const ArcheMetalCell* cell,
                                ArcheMetalState* st);

// ── Table loaders (return NULL on failure; see arche_last_error) ───────────
ArchePrimTable* arche_prim_table_load(const char* h5_path);
void arche_prim_table_free(ArchePrimTable* tbl);
ArcheMetalTable* arche_metal_table_load(const char* h5_path);
void arche_metal_table_free(ArcheMetalTable* tbl);

// ── Stepping (returns ARCHE_OK or an ARCHE_ERR_* code) ─────────────────────
// On success the full cooling/heating breakdown is written to *out (out may be
// NULL to discard it). A solver failure is reported via out->solver_failed = 1
// with a return of ARCHE_OK (it is a numerical outcome, not a boundary error).
int arche_chem_full_step_prim(ArchePrimCell* cell, double dt,
                              const ArcheChemParams* params,
                              const ArcheChemShielding* shield,
                              const ArchePrimTable* tbl,
                              ArcheChemFullRates* out);
int arche_chem_full_step_metal(ArcheMetalCell* cell, double dt,
                               const ArcheChemParams* params,
                               const ArcheChemShielding* shield,
                               const ArcheMetalTable* tbl,
                               ArcheChemFullRates* out);

// ── Runtime model dispatch ──────────────────────────────────────────────────
// C view of the C++ runtime-dispatch layer: choose the model at runtime, then
// drive it through one set of calls. The species-vector length differs by
// model (Prim = ARCHE_PRIM_NSP, Metal = ARCHE_METAL_NSP), so the unified cell
// exposes y as a pointer of length arche_cell_nsp().
typedef enum { ARCHE_MODEL_PRIM = 0, ARCHE_MODEL_METAL = 1 } ArcheModel;

typedef struct ArcheCell ArcheCell;    // tagged opaque cell
typedef struct ArcheTable ArcheTable;  // tagged opaque table

ArcheCell* arche_cell_create(ArcheModel model);
void arche_cell_free(ArcheCell* cell);
ArcheTable* arche_table_load(ArcheModel model,
                             const char* h5_path);  // NULL on err
void arche_table_free(ArcheTable* tbl);

ArcheModel arche_cell_model(const ArcheCell* cell);
int arche_cell_nsp(
    const ArcheCell* cell);  // species-vector length for this model

// Mutable pointer to the cell's species vector (length arche_cell_nsp()).
double* arche_cell_y(ArcheCell* cell);
void arche_cell_set_scalars(ArcheCell* cell, double nH, double T_K, double mu,
                            double gamma);
void arche_cell_get_scalars(const ArcheCell* cell, double* nH, double* T_K,
                            double* mu, double* gamma);

// Dispatches on the cell's model (cell and table models must match, else
// ARCHE_ERR_INVALID_ARGUMENT). Returns ARCHE_OK or an ARCHE_ERR_* code.
int arche_chem_full_step(ArcheCell* cell, double dt,
                         const ArcheChemParams* params,
                         const ArcheChemShielding* shield,
                         const ArcheTable* tbl, ArcheChemFullRates* out);

// ── Registry-backed, name-selected models (species metadata) ────────────────
// Choose a model by name at runtime and read its species count and host-facing
// species names, so a caller can map its own columns without knowing the
// compile-time species count.  Purely additive over the per-model and
// ArcheModel-enum entry points above.  A model handle owns its reaction table;
// cells are created from it.  The species-vector length is
// arche_model_n_species().
typedef struct ArcheModelRuntime ArcheModelRuntime;  // named model + its table
typedef struct ArcheModelCell ArcheModelCell;        // one cell of that model

// create returns NULL on unknown name or table-load failure (arche_last_error).
ArcheModelRuntime* arche_model_create(const char* name, const char* h5_path);
void arche_model_destroy(ArcheModelRuntime* model);

// Species metadata. n_species returns -1, the pointer getters NULL, if model
// is NULL. arche_model_species returns an array of n_species C strings.
int arche_model_n_species(const ArcheModelRuntime* model);
const char* const* arche_model_species(const ArcheModelRuntime* model);
const char* arche_model_name(const ArcheModelRuntime* model);

// Registry introspection: number of registered models and their names.
int arche_model_count(void);
const char* arche_model_registry_name(int i);  // NULL if i out of range

// Per-cell lifecycle + state. arche_model_cell_y returns a mutable pointer to
// the species vector (length arche_model_n_species()).
ArcheModelCell* arche_model_cell_create(const ArcheModelRuntime* model);
void arche_model_cell_free(ArcheModelCell* cell);
double* arche_model_cell_y(ArcheModelCell* cell);
void arche_model_cell_set_scalars(ArcheModelCell* cell, double nH, double T_K,
                                  double mu, double gamma);
void arche_model_cell_get_scalars(const ArcheModelCell* cell, double* nH,
                                  double* T_K, double* mu, double* gamma);

// Reset only integration history; visible state is unchanged. No-op on NULL.
void arche_model_cell_reset(ArcheModelCell* cell);

// Read-only thermodynamic maps over the current composition. mu is in proton-
// mass units, gamma is dimensionless, T_K is [K], and e_cgs is [erg g^-1].
// Each query returns a quiet NaN when cell is NULL. T_from_e also propagates a
// NaN e_cgs input. gamma_from_y returns quiet NaN unless T_K is positive and
// finite. e(T) is monotone over the search bracket, so T_from_e returns the
// unique root; see arche_api.h for the bracket and its measured accuracy.
double arche_model_mu_from_y(const ArcheModelCell* cell);
double arche_model_gamma_from_y(const ArcheModelCell* cell, double T_K);
double arche_model_T_from_e(const ArcheModelCell* cell, double e_cgs);

// Advance one cell by dt. Returns ARCHE_OK or an ARCHE_ERR_* code; on success
// the cooling/heating breakdown is written to *out (out may be NULL).  The
// cell must have been created from the same model (else
// ARCHE_ERR_INVALID_ARGUMENT).
int arche_model_step(const ArcheModelRuntime* model, ArcheModelCell* cell,
                     double dt, const ArcheChemParams* params,
                     const ArcheChemShielding* shield, ArcheChemFullRates* out);

// Human-readable message for the most recent failure on the calling thread
// (empty string if none).
const char* arche_last_error(void);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // SRC_INCLUDE_API_ARCHE_CAPI_H_
