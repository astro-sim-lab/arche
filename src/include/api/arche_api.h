// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// arche_api.h — stable, non-template public facade for libarche
//
// This header is the *only* include an external application needs in order to
// drive the chemistry network.  It is deliberately lightweight:
//
//   * NO Eigen          — the linear solver (core/newton.h) is hidden inside
//   libarche.
//   * NO HDF5           — the table loader (kinetics/topology.h) is hidden too.
//   * NO templates in the public signatures — every entry point is a concrete,
//                         per-model function, so callers do not instantiate the
//                         kernel and the network sizes (N_sp / N_react) never
//                         leak into the API/ABI.
//
// The implementation (explicit instantiation of chem_(full_)step<Model> and the
// opaque-handle bridges) lives in src/lib/arche_api.cc.  This header is
// declarations only.
//
// ── Design notes ──────────────────────────────────────────────────────────
//   POD structs (ChemParams / ChemShielding / ChemRates / ChemFullRates) are
//   already Eigen-free plain data and are exposed directly via "core/state.h",
//   which pulls in nothing heavier than <array>/<cmath>/<limits> and the public
//   species metadata (core/species_index.h: per-model Sp index enums and
//   abundance_ref::).  Including this header is therefore sufficient on its own
//   to index y[] by species name; a separate species include is optional.
//
//   The sized kernel types are the problem children:
//       ChemCell<Model>                 — holds var[2*N_react] + EscapeState
//       ReactionTable<N_sp, N_react>    — loaded from HDF5
//   Exposing them as concrete aliases (ZeroMetalCell, ZeroMetalTable, ...)
//   would drag Eigen (solve/chemistry.h) and HDF5 (kinetics/topology.h) back
//   into every
//   consumer TU and would bake the compile-time network sizes into the API.
//   They are therefore hidden behind *opaque handles* (PrimCell / MetalCell /
//   PrimTable / MetalTable).  Handles are heap-allocated by the create/load
//   functions and released by the matching destroy functions (or the RAII
//   typedefs below).  The opaque names intentionally differ from the internal
//   ZeroMetalCell / MetalGrainCell aliases so that arche_api.cc can include
//   both this header and chemistry.h without a redefinition clash.
//
//   The mutable per-cell state the caller must fill each step (y[], nH, T_K,
//   mu, gamma) is itself a plain POD (ChemState<N_sp>) with no heavy deps, so
//   it is *not* hidden: an accessor returns a reference to the concrete
//   per-model alias (ChemStateZM / ChemStateMG) embedded in the opaque cell,
//   letting callers read/write the species vector directly.
//
//   Errors are still reported as C++ exceptions (std::invalid_argument /
//   std::runtime_error), exactly as the kernel does today.  The C ABI in
//   arche_capi.h (extern "C") translates them to status codes for a
//   C/Fortran boundary.
//
// ── Minimal usage (primordial / zero-metal network) ─────────────────────────
//   #include "api/arche_api.h"
//   using namespace arche;
//
//   PrimTablePtr tbl = load_prim_table_owned("/path/to/primordial.h5");
//   PrimCellPtr  cell = prim_cell_create_owned();
//
//   ChemStateZM& s = prim_cell_state(*cell);
//   s.nH = 1.0e2; s.T_K = 100.0; s.mu = 1.22; s.gamma = 5.0/3.0;
//   // ... initialise s.y[0..22] ...
//
//   ChemParams    params;  params.zeta = 1.0e-17; params.T_rad = 2.725;
//   ChemShielding shield;  shield.zeta = params.zeta;
//
//   ChemFullRates r = chem_full_step_prim(*cell, dt, params, shield, *tbl);
//   // r.Lambda_net, r.Gamma_CR, ...
// ---------------------------------------------------------------------------

#include <memory>
#include <string>

#include "core/state.h"  // ChemParams, ChemShielding, ChemRates, ChemFullRates,
// ChemState<N_sp>, ChemStateZM / ChemStateMG, phys::...
// (no Eigen, no HDF5)

namespace arche {

// Compact-minimal per-cell state (15 species).  Declared as a fixed-size alias
// so this header stays free of the model headers; arche_api.cc static_asserts
// it equals zero_metal_minimal::N_sp.
using ChemStatePrimMinimal = ChemState<15>;

// Compact metal-minimal per-cell state (40 species).  Same convention as
// ChemStatePrimMinimal; arche_api.cc static_asserts it equals
// metal_grain_minimal::N_sp.
using ChemStateMetalMinimal = ChemState<40>;

// ---------------------------------------------------------------------------
// Opaque handles.
//
// Forward declarations only; the complete definitions
//   PrimCell         ≙ ChemCell<Nakauchi2019>
//   PrimMinimalCell  ≙ ChemCell<Nakauchi2019_Minimal>   (compact 15/33)
//   MetalCell        ≙ ChemCell<Nakauchi2021>
//   PrimTable        ≙ ReactionTable<zero_metal::N_sp,  zero_metal::N_react>
//   PrimMinimalTable ≙ ReactionTable<zero_metal_minimal::N_sp, 33>
//   MetalTable       ≙ ReactionTable<metal_grain::N_sp, metal_grain::N_react>
// are visible only inside libarche (arche_api.cc).  Consumers may
// only hold pointers/references to these types.
// ---------------------------------------------------------------------------
struct PrimCell;           // primordial         (Nakauchi2019, 23 species)
struct PrimMinimalCell;    // primordial minimal (Nakauchi2019_Minimal, compact
                           // 15 species / 33 reactions)
struct MetalCell;          // metal+grain        (Nakauchi2021, 89 species)
struct MetalMinimalCell;   // metal minimal      (Nakauchi2021_Minimal, compact
                           // 40 species / 103+10 reactions)
struct PrimTable;          // reaction table for the primordial network
struct PrimMinimalTable;   // reaction table for the compact minimal network
struct MetalTable;         // reaction table for the metal+grain network
struct MetalMinimalTable;  // reaction table for the compact metal minimal
                           // network

// ---------------------------------------------------------------------------
// Cell lifecycle.
//
// One cell owns the per-cell state plus the inter-step reaction-rate cache.
// Cells are NOT thread-safe: each thread must own its own cell.  Tables are
// read-only and may be shared across threads/cells.
// ---------------------------------------------------------------------------
PrimCell* prim_cell_create();
void prim_cell_destroy(PrimCell* cell) noexcept;

PrimMinimalCell* prim_minimal_cell_create();
void prim_minimal_cell_destroy(PrimMinimalCell* cell) noexcept;

MetalCell* metal_cell_create();
void metal_cell_destroy(MetalCell* cell) noexcept;

MetalMinimalCell* metal_minimal_cell_create();
void metal_minimal_cell_destroy(MetalMinimalCell* cell) noexcept;

// ---------------------------------------------------------------------------
// Direct access to the mutable per-cell state embedded in the opaque cell.
//
// Returns a reference to the concrete POD state (ChemStateZM / ChemStateMG) so
// callers can set y[], nH, T_K, mu, gamma without any accessor boilerplate.
// ---------------------------------------------------------------------------
ChemStateZM& prim_cell_state(PrimCell& cell) noexcept;
const ChemStateZM& prim_cell_state(const PrimCell& cell) noexcept;

// PrimMinimalCell carries the compact 15-species state.
ChemStatePrimMinimal& prim_minimal_cell_state(PrimMinimalCell& cell) noexcept;
const ChemStatePrimMinimal& prim_minimal_cell_state(
    const PrimMinimalCell& cell) noexcept;

ChemStateMG& metal_cell_state(MetalCell& cell) noexcept;
const ChemStateMG& metal_cell_state(const MetalCell& cell) noexcept;

// MetalMinimalCell carries the compact 40-species state.
ChemStateMetalMinimal& metal_minimal_cell_state(
    MetalMinimalCell& cell) noexcept;
const ChemStateMetalMinimal& metal_minimal_cell_state(
    const MetalMinimalCell& cell) noexcept;

// ---------------------------------------------------------------------------
// Table loaders.
//
// Build a model's reaction table: the topology comes from the C++ network
// source inside libarche, the partition-function tables from the model's HDF5
// data file (data/primordial.h5 / data/metal_grain.h5).  All HDF5 handling is
// confined to libarche; this header has no <hdf5.h> dependency.  Throws
// std::runtime_error on a missing file or a network-size mismatch.
// ---------------------------------------------------------------------------
PrimTable* load_prim_table(const std::string& h5_path);
void prim_table_destroy(PrimTable* tbl) noexcept;

MetalTable* load_metal_table(const std::string& h5_path);
void metal_table_destroy(MetalTable* tbl) noexcept;

// The compact minimal network has its own 15-species / 33-reaction table:
// a compact topology with partition functions remapped from the full
// primordial HDF5 file (data/primordial.h5).
PrimMinimalTable* load_minimal_table(const std::string& h5_path);
void prim_minimal_table_destroy(PrimMinimalTable* tbl) noexcept;

MetalMinimalTable* load_metal_minimal_table(const std::string& h5_path);
void metal_minimal_table_destroy(MetalMinimalTable* tbl) noexcept;

// ---------------------------------------------------------------------------
// Invariant-row count of a loaded table (see solve/conservation.h).
//
// Returns the number of element/charge rows the loaded table is CONFIGURED with:
// 0 means the network registers none, so chem_full_step_* leaves
// ChemFullRates::conservation_projected false on every step by design.
//
// ⚠ This is not the number of rows enforced on a given step.  The projection
// weights each species by its own abundance, so a row whose carrier species are
// all zero drops out of that step's solve — correctly, since it owes nothing —
// while the step still reports conservation_projected = true.
//
// Read through the table the caller actually holds, not from the make_*
// factory.  A factory copy that inlines the topology steps can drop the
// invariant rows, which would let the factory report the expected count while
// every app going through the facade ran with none.
// ---------------------------------------------------------------------------
int prim_table_n_invariants(const PrimTable& tbl) noexcept;
int prim_minimal_table_n_invariants(const PrimMinimalTable& tbl) noexcept;
int metal_table_n_invariants(const MetalTable& tbl) noexcept;
int metal_minimal_table_n_invariants(const MetalMinimalTable& tbl) noexcept;

// ---------------------------------------------------------------------------
// Stepping entry points (non-template; one pair per model).
//
//   chem_full_step_* — advance one dt AND return the full cooling/heating
//                      breakdown (line, continuum, chemistry, CR, X-ray, net,
//                      opacities, solved grain temperature).
//   chem_step_*      — advance one dt and return only the scalar chemistry
//                      cooling + CR heating (low-level; caller does line and
//                      continuum cooling itself).
//
// On a solver failure chem_full_step_* leaves the cell state rolled back and
// sets ChemFullRates::solver_failed = true.  T_rad (and, for metal, T_gr_K)
// must be finite and positive in `params` or std::invalid_argument is thrown.
// ---------------------------------------------------------------------------
ChemFullRates chem_full_step_prim(PrimCell& cell, double dt,
                                  const ChemParams& params,
                                  const ChemShielding& shield,
                                  const PrimTable& tbl);
ChemFullRates chem_full_step_prim_minimal(PrimMinimalCell& cell, double dt,
                                          const ChemParams& params,
                                          const ChemShielding& shield,
                                          const PrimMinimalTable& tbl);
ChemFullRates chem_full_step_metal_minimal(MetalMinimalCell& cell, double dt,
                                           const ChemParams& params,
                                           const ChemShielding& shield,
                                           const MetalMinimalTable& tbl);
ChemFullRates chem_full_step_metal(MetalCell& cell, double dt,
                                   const ChemParams& params,
                                   const ChemShielding& shield,
                                   const MetalTable& tbl);

ChemRates chem_step_prim(PrimCell& cell, double dt, const ChemParams& params,
                         const PrimTable& tbl);
ChemRates chem_step_prim_minimal(PrimMinimalCell& cell, double dt,
                                 const ChemParams& params,
                                 const PrimMinimalTable& tbl);
ChemRates chem_step_metal_minimal(MetalMinimalCell& cell, double dt,
                                  const ChemParams& params,
                                  const MetalMinimalTable& tbl);
ChemRates chem_step_metal(MetalCell& cell, double dt, const ChemParams& params,
                          const MetalTable& tbl);

// ---------------------------------------------------------------------------
// RAII convenience.
//
// Header-only smart-pointer wrappers around the create/load + destroy pairs.
// The custom deleters only forward a pointer to an incomplete type, so they
// compile without the handle definitions being visible here.
// ---------------------------------------------------------------------------
struct PrimCellDeleter {
  void operator()(PrimCell* p) const noexcept { prim_cell_destroy(p); }
};
struct PrimMinimalCellDeleter {
  void operator()(PrimMinimalCell* p) const noexcept {
    prim_minimal_cell_destroy(p);
  }
};
struct MetalCellDeleter {
  void operator()(MetalCell* p) const noexcept { metal_cell_destroy(p); }
};
struct MetalMinimalCellDeleter {
  void operator()(MetalMinimalCell* p) const noexcept {
    metal_minimal_cell_destroy(p);
  }
};
struct PrimTableDeleter {
  void operator()(PrimTable* p) const noexcept { prim_table_destroy(p); }
};
struct PrimMinimalTableDeleter {
  void operator()(PrimMinimalTable* p) const noexcept {
    prim_minimal_table_destroy(p);
  }
};
struct MetalTableDeleter {
  void operator()(MetalTable* p) const noexcept { metal_table_destroy(p); }
};
struct MetalMinimalTableDeleter {
  void operator()(MetalMinimalTable* p) const noexcept {
    metal_minimal_table_destroy(p);
  }
};

using PrimCellPtr = std::unique_ptr<PrimCell, PrimCellDeleter>;
using PrimMinimalCellPtr =
    std::unique_ptr<PrimMinimalCell, PrimMinimalCellDeleter>;
using MetalCellPtr = std::unique_ptr<MetalCell, MetalCellDeleter>;
using MetalMinimalCellPtr =
    std::unique_ptr<MetalMinimalCell, MetalMinimalCellDeleter>;
using PrimTablePtr = std::unique_ptr<PrimTable, PrimTableDeleter>;
using PrimMinimalTablePtr =
    std::unique_ptr<PrimMinimalTable, PrimMinimalTableDeleter>;
using MetalTablePtr = std::unique_ptr<MetalTable, MetalTableDeleter>;
using MetalMinimalTablePtr =
    std::unique_ptr<MetalMinimalTable, MetalMinimalTableDeleter>;

inline PrimCellPtr prim_cell_create_owned() {
  return PrimCellPtr{prim_cell_create()};
}
inline PrimMinimalCellPtr prim_minimal_cell_create_owned() {
  return PrimMinimalCellPtr{prim_minimal_cell_create()};
}
inline MetalCellPtr metal_cell_create_owned() {
  return MetalCellPtr{metal_cell_create()};
}
inline MetalMinimalCellPtr metal_minimal_cell_create_owned() {
  return MetalMinimalCellPtr{metal_minimal_cell_create()};
}
inline PrimTablePtr load_prim_table_owned(const std::string& h5_path) {
  return PrimTablePtr{load_prim_table(h5_path)};
}
inline PrimMinimalTablePtr load_minimal_table_owned(
    const std::string& h5_path) {
  return PrimMinimalTablePtr{load_minimal_table(h5_path)};
}
inline MetalTablePtr load_metal_table_owned(const std::string& h5_path) {
  return MetalTablePtr{load_metal_table(h5_path)};
}
inline MetalMinimalTablePtr load_metal_minimal_table_owned(
    const std::string& h5_path) {
  return MetalMinimalTablePtr{load_metal_minimal_table(h5_path)};
}

// ===========================================================================
// Runtime model dispatch
//
// The per-model entry points above (chem_full_step_prim / _metal, PrimCell /
// MetalCell, ...) fix the model at compile time. This layer adds a thin
// dispatch path for consumers that choose the model at *runtime* (e.g. from a
// config file):
// one Model enum, one tagged opaque Cell / Table handle, and one chem_full_step
// that delegates to the matching per-model entry. It is purely additive —
// numerically identical to calling chem_full_step_<model> directly — and the
// in-tree model-fixed apps need not use it.
//
// Because the species-vector length differs by model (Prim = 23, Metal = 89),
// the unified state is exposed as a pointer + length (cell_y / cell_n_sp)
// rather than a typed ChemState; scalars are reference accessors.
// ===========================================================================
enum class Model { Prim, Metal };

// Tagged opaque handles (definitions live in arche_api.cc).
struct Cell;
struct Table;

// ── lifecycle (model chosen at runtime) ────────────────────────────────────
Cell* cell_create(Model m);
void cell_destroy(Cell* cell) noexcept;
Table* load_table(Model m, const std::string& h5_path);
void table_destroy(Table* tbl) noexcept;

// ── introspection ──────────────────────────────────────────────────────────
Model cell_model(const Cell& cell) noexcept;
Model table_model(const Table& tbl) noexcept;
int cell_n_sp(const Cell& cell) noexcept;  // 23 (Prim) or 89 (Metal)

// ── per-cell state (length-erased: y is a pointer of length cell_n_sp) ─────
double* cell_y(Cell& cell) noexcept;
const double* cell_y(const Cell& cell) noexcept;
double& cell_nH(Cell& cell) noexcept;
double& cell_T_K(Cell& cell) noexcept;
double& cell_mu(Cell& cell) noexcept;
double& cell_gamma(Cell& cell) noexcept;

// ── stepping (dispatches on the cell's model; cell and table models must
//    match or std::invalid_argument is thrown) ─────────────────────────────
ChemFullRates chem_full_step(Cell& cell, double dt, const ChemParams& params,
                             const ChemShielding& shield, const Table& tbl);

// ── RAII convenience ───────────────────────────────────────────────────────
struct CellDeleter {
  void operator()(Cell* p) const noexcept { cell_destroy(p); }
};
struct TableDeleter {
  void operator()(Table* p) const noexcept { table_destroy(p); }
};
using CellPtr = std::unique_ptr<Cell, CellDeleter>;
using TablePtr = std::unique_ptr<Table, TableDeleter>;

inline CellPtr cell_create_owned(Model m) { return CellPtr{cell_create(m)}; }
inline TablePtr load_table_owned(Model m, const std::string& h5_path) {
  return TablePtr{load_table(m, h5_path)};
}

// ===========================================================================
// Registry-backed, name-selected models (type-erased; species metadata)
//
// A model is chosen by *name* at runtime (e.g. "Nakauchi2019_Minimal"); the
// handle owns its reaction table and exposes the species count and the
// host-facing species names, so a caller can map its own columns without
// knowing the compile-time N_sp.  Per-cell stepping reuses the same kernel path
// as the per-model entry points (numerically identical).  Adding a model is a
// single registry row (see arche_api.cc).  This layer is purely additive: the
// per-model symbols and the Model-enum dispatch above are unchanged.
//
//   ModelRuntime* m = model_create("Nakauchi2019_Minimal",
//   "data/primordial.h5"); int n = model_n_species(*m);                  // 15
//   const char* const* names = model_species(*m); // {"H","H2",...}
//   ModelCell* c = model_cell_create(*m);
//   double* y = model_cell_y(*c);                 // length n
//   model_cell_nH(*c) = 1.0e4; /* T_K, mu, gamma ... */
//   ChemFullRates r = model_step(*m, *c, dt, params, shield);
//   model_cell_destroy(c); model_destroy(m);
// ===========================================================================
struct ModelRuntime;  // opaque: a named model plus its loaded reaction table
struct ModelCell;     // opaque: one cell's state for a ModelRuntime

// Create a model by registered name; loads its reaction table from h5_path.
// Throws std::invalid_argument (unknown name) or std::runtime_error (load).
ModelRuntime* model_create(const std::string& name, const std::string& h5_path);
void model_destroy(ModelRuntime* m) noexcept;

// Species metadata (lets a host map columns without knowing N_sp).
int model_n_species(const ModelRuntime& m) noexcept;
const char* const* model_species(
    const ModelRuntime& m) noexcept;  // length N_sp
const char* model_name(const ModelRuntime& m) noexcept;

// Registry introspection: the set of registered model names.
int model_count() noexcept;
const char* model_registry_name(int i) noexcept;  // nullptr if out of range

// Per-cell lifecycle + state (y has length model_n_species).
ModelCell* model_cell_create(const ModelRuntime& m);
void model_cell_destroy(ModelCell* c) noexcept;
double* model_cell_y(ModelCell& c) noexcept;
const double* model_cell_y(const ModelCell& c) noexcept;
double& model_cell_nH(ModelCell& c) noexcept;
double& model_cell_T_K(ModelCell& c) noexcept;
double& model_cell_mu(ModelCell& c) noexcept;
double& model_cell_gamma(ModelCell& c) noexcept;

// Advance one cell by dt (same kernel path as chem_full_step_<model>).
// The cell must have been created from the same ModelRuntime (or one of the
// same model); otherwise std::invalid_argument is thrown.
ChemFullRates model_step(const ModelRuntime& m, ModelCell& c, double dt,
                         const ChemParams& params, const ChemShielding& shield);

struct ModelRuntimeDeleter {
  void operator()(ModelRuntime* p) const noexcept { model_destroy(p); }
};
struct ModelCellDeleter {
  void operator()(ModelCell* p) const noexcept { model_cell_destroy(p); }
};
using ModelRuntimePtr = std::unique_ptr<ModelRuntime, ModelRuntimeDeleter>;
using ModelCellPtr = std::unique_ptr<ModelCell, ModelCellDeleter>;

inline ModelRuntimePtr model_create_owned(const std::string& name,
                                          const std::string& h5_path) {
  return ModelRuntimePtr{model_create(name, h5_path)};
}
inline ModelCellPtr model_cell_create_owned(const ModelRuntime& m) {
  return ModelCellPtr{model_cell_create(m)};
}

}  // namespace arche
