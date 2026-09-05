// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
// ---------------------------------------------------------------------------
// arche_api.cc — implementation of the non-template public facade.
//
// This is the single translation unit that includes the full templated kernel
// (chemistry.h, which pulls in Eigen via chemreact.h and HDF5 via
// kinetics/topology.h) and instantiates it for the two supported models.
// Linking libarche therefore gives an external application the complete
// chemistry kernel without recompiling any template and without seeing
// Eigen/HDF5 in its own translation units (it includes only arche_api.h).
//
//   * The opaque handles declared in arche_api.h (PrimCell / MetalCell /
//     PrimTable / MetalTable) are defined here as thin wrappers over the
//     concrete kernel types (ChemCell<Model> / ReactionTable<...>).
//   * chem_(full_)step<Nakauchi2019|Nakauchi2021> are explicitly instantiated
//     below so the kernel object code lives here, in libarche, once.
// ---------------------------------------------------------------------------
#include "api/arche_api.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "solve/chemistry.h"  // ChemCell<Model>, chem_(full_)step<Model>,
                              // make_*_table, <model>::net::init_topology,
                              // load_pf_tables_h5, ReactionTable<...>

namespace arche {

// ---------------------------------------------------------------------------
// Opaque handle definitions.
//
// Composition (not inheritance) keeps the boundary explicit: the facade owns a
// concrete kernel object and exposes only what arche_api.h declares.  Because
// PrimCell::impl is a ChemCell<Nakauchi2019>, its `.state` is exactly the POD
// ChemStateZM that prim_cell_state() hands back by reference.
// ---------------------------------------------------------------------------
struct PrimCell {
  ZeroMetalCell impl;
};
struct PrimMinimalCell {
  // Compact 15-species / 33-reaction minimal network.
  ChemCell<Nakauchi2019_Minimal> impl;
};
struct MetalCell {
  MetalGrainCell impl;
};
struct MetalMinimalCell {
  // Compact 40-species / 103+10-reaction metal minimal network.
  ChemCell<Nakauchi2021_Minimal> impl;
};
struct PrimTable {
  ZeroMetalTable impl;
};
struct PrimMinimalTable {
  zero_metal_minimal::MinimalTable impl;
};
struct MetalTable {
  MetalGrainTable impl;
};
struct MetalMinimalTable {
  metal_grain_minimal::MinimalTable impl;
};
static_assert(zero_metal_minimal::N_sp == 15,
              "ChemStatePrimMinimal in arche_api.h assumes 15 compact species");
static_assert(
    metal_grain_minimal::N_sp == 40,
    "ChemStateMetalMinimal in arche_api.h assumes 40 compact species");

// ── Cell lifecycle ─────────────────────────────────────────────────────────
PrimCell* prim_cell_create() { return new PrimCell(); }
void prim_cell_destroy(PrimCell* cell) noexcept { delete cell; }

PrimMinimalCell* prim_minimal_cell_create() { return new PrimMinimalCell(); }
void prim_minimal_cell_destroy(PrimMinimalCell* cell) noexcept { delete cell; }

MetalCell* metal_cell_create() { return new MetalCell(); }
void metal_cell_destroy(MetalCell* cell) noexcept { delete cell; }

MetalMinimalCell* metal_minimal_cell_create() { return new MetalMinimalCell(); }
void metal_minimal_cell_destroy(MetalMinimalCell* cell) noexcept {
  delete cell;
}

// ── Per-cell state access ────────────────────────────────────────────────--
ChemStateZM& prim_cell_state(PrimCell& cell) noexcept {
  return cell.impl.state;
}
const ChemStateZM& prim_cell_state(const PrimCell& cell) noexcept {
  return cell.impl.state;
}
ChemStatePrimMinimal& prim_minimal_cell_state(PrimMinimalCell& cell) noexcept {
  return cell.impl.state;
}
const ChemStatePrimMinimal& prim_minimal_cell_state(
    const PrimMinimalCell& cell) noexcept {
  return cell.impl.state;
}
ChemStateMG& metal_cell_state(MetalCell& cell) noexcept {
  return cell.impl.state;
}
const ChemStateMG& metal_cell_state(const MetalCell& cell) noexcept {
  return cell.impl.state;
}
ChemStateMetalMinimal& metal_minimal_cell_state(
    MetalMinimalCell& cell) noexcept {
  return cell.impl.state;
}
const ChemStateMetalMinimal& metal_minimal_cell_state(
    const MetalMinimalCell& cell) noexcept {
  return cell.impl.state;
}

// ── Table loaders ────────────────────────────────────────────────────────--
PrimTable* load_prim_table(const std::string& h5_path) {
  auto* tbl = new PrimTable();
  try {
    // Go through make_zero_metal_table rather than repeating its body: it also
    // derives the element/charge invariant rows the conservation projection
    // needs (solve/conservation.h), and an inlined copy of the topology+PF
    // steps silently omitted them, leaving the full network unprojected while
    // the compact one (which does call its make_* function) was corrected.
    tbl->impl = make_zero_metal_table(h5_path);
  } catch (...) {
    delete tbl;
    throw;
  }
  return tbl;
}
void prim_table_destroy(PrimTable* tbl) noexcept { delete tbl; }

PrimMinimalTable* load_minimal_table(const std::string& h5_path) {
  auto* tbl = new PrimMinimalTable();
  try {
    tbl->impl = make_minimal_table(h5_path);  // compact topology + remapped PF
  } catch (...) {
    delete tbl;
    throw;
  }
  return tbl;
}
void prim_minimal_table_destroy(PrimMinimalTable* tbl) noexcept { delete tbl; }

MetalMinimalTable* load_metal_minimal_table(const std::string& h5_path) {
  auto* tbl = new MetalMinimalTable();
  try {
    // Build the full PF-loaded table once, then derive the compact table
    // (which retains the full table via aux_full_metal for the strategy-1 rate
    // path and the compact Saha).
    tbl->impl = make_metal_minimal_table(make_metal_grain_table(h5_path));
  } catch (...) {
    delete tbl;
    throw;
  }
  return tbl;
}
void metal_minimal_table_destroy(MetalMinimalTable* tbl) noexcept {
  delete tbl;
}

MetalTable* load_metal_table(const std::string& h5_path) {
  auto* tbl = new MetalTable();
  try {
    // Go through make_metal_grain_table for the same reason load_prim_table
    // goes through make_zero_metal_table: it also derives the element/charge
    // invariant rows.  Inlining the topology+PF steps here instead would omit
    // them silently, and every app reaching this network through the facade
    // would run unprojected.  Every load_* entry point calls its make_*
    // function for this reason.
    tbl->impl = make_metal_grain_table(h5_path);
  } catch (...) {
    delete tbl;
    throw;
  }
  return tbl;
}
void metal_table_destroy(MetalTable* tbl) noexcept { delete tbl; }

// ── Invariant-row count ────────────────────────────────────────────────────
// Reads the row count off the table the caller holds, so it reports the facade
// path rather than the make_* factory's return value.
int prim_table_n_invariants(const PrimTable& tbl) noexcept {
  return tbl.impl.n_invariants;
}
int prim_minimal_table_n_invariants(const PrimMinimalTable& tbl) noexcept {
  return tbl.impl.n_invariants;
}
int metal_table_n_invariants(const MetalTable& tbl) noexcept {
  return tbl.impl.n_invariants;
}
int metal_minimal_table_n_invariants(const MetalMinimalTable& tbl) noexcept {
  return tbl.impl.n_invariants;
}

// ── Stepping entry points ──────────────────────────────────────────────────
// Thin forwarders to the templated kernel.  No transformation of arguments or
// results, so the facade is numerically identical to calling the template.
ChemFullRates chem_full_step_prim(PrimCell& cell, double dt,
                                  const ChemParams& params,
                                  const ChemShielding& shield,
                                  const PrimTable& tbl) {
  return chem_full_step<Nakauchi2019>(cell.impl, dt, params, shield, tbl.impl);
}
ChemFullRates chem_full_step_prim_minimal(PrimMinimalCell& cell, double dt,
                                          const ChemParams& params,
                                          const ChemShielding& shield,
                                          const PrimMinimalTable& tbl) {
  return chem_full_step<Nakauchi2019_Minimal>(cell.impl, dt, params, shield,
                                              tbl.impl);
}
ChemFullRates chem_full_step_metal(MetalCell& cell, double dt,
                                   const ChemParams& params,
                                   const ChemShielding& shield,
                                   const MetalTable& tbl) {
  return chem_full_step<Nakauchi2021>(cell.impl, dt, params, shield, tbl.impl);
}
ChemFullRates chem_full_step_metal_minimal(MetalMinimalCell& cell, double dt,
                                           const ChemParams& params,
                                           const ChemShielding& shield,
                                           const MetalMinimalTable& tbl) {
  return chem_full_step<Nakauchi2021_Minimal>(cell.impl, dt, params, shield,
                                              tbl.impl);
}

ChemRates chem_step_prim(PrimCell& cell, double dt, const ChemParams& params,
                         const PrimTable& tbl) {
  return chem_step<Nakauchi2019>(cell.impl, dt, params, tbl.impl);
}
ChemRates chem_step_prim_minimal(PrimMinimalCell& cell, double dt,
                                 const ChemParams& params,
                                 const PrimMinimalTable& tbl) {
  return chem_step<Nakauchi2019_Minimal>(cell.impl, dt, params, tbl.impl);
}
ChemRates chem_step_metal(MetalCell& cell, double dt, const ChemParams& params,
                          const MetalTable& tbl) {
  return chem_step<Nakauchi2021>(cell.impl, dt, params, tbl.impl);
}
ChemRates chem_step_metal_minimal(MetalMinimalCell& cell, double dt,
                                  const ChemParams& params,
                                  const MetalMinimalTable& tbl) {
  return chem_step<Nakauchi2021_Minimal>(cell.impl, dt, params, tbl.impl);
}

// ---------------------------------------------------------------------------
// Runtime model dispatch.
//
// The unified Cell/Table own one of the per-model handles and carry a Model
// tag. Every operation switches on the tag and delegates to the per-model
// facade above, so the dispatch path is numerically identical to calling
// chem_full_step_<model> directly (bit-for-bit).
// ---------------------------------------------------------------------------
struct Cell {
  Model model;
  PrimCell* prim = nullptr;    // non-null iff model == Prim
  MetalCell* metal = nullptr;  // non-null iff model == Metal
};
struct Table {
  Model model;
  PrimTable* prim = nullptr;
  MetalTable* metal = nullptr;
};

Cell* cell_create(Model m) {
  auto* c = new Cell();
  c->model = m;
  if (m == Model::Prim)
    c->prim = prim_cell_create();
  else
    c->metal = metal_cell_create();
  return c;
}
void cell_destroy(Cell* c) noexcept {
  if (!c) return;
  prim_cell_destroy(c->prim);  // null-safe (delete nullptr is fine)
  metal_cell_destroy(c->metal);
  delete c;
}
Table* load_table(Model m, const std::string& h5_path) {
  auto* t = new Table();
  t->model = m;
  try {
    if (m == Model::Prim)
      t->prim = load_prim_table(h5_path);
    else
      t->metal = load_metal_table(h5_path);
  } catch (...) {
    delete t;  // inner loader threw before storing; nothing else to free
    throw;
  }
  return t;
}
void table_destroy(Table* t) noexcept {
  if (!t) return;
  prim_table_destroy(t->prim);
  metal_table_destroy(t->metal);
  delete t;
}

Model cell_model(const Cell& c) noexcept { return c.model; }
Model table_model(const Table& t) noexcept { return t.model; }
int cell_n_sp(const Cell& c) noexcept {
  return c.model == Model::Prim ? Nakauchi2019::N_sp : Nakauchi2021::N_sp;
}

double* cell_y(Cell& c) noexcept {
  return c.model == Model::Prim ? prim_cell_state(*c.prim).y.data()
                                : metal_cell_state(*c.metal).y.data();
}
const double* cell_y(const Cell& c) noexcept {
  return c.model == Model::Prim ? prim_cell_state(*c.prim).y.data()
                                : metal_cell_state(*c.metal).y.data();
}
double& cell_nH(Cell& c) noexcept {
  return c.model == Model::Prim ? prim_cell_state(*c.prim).nH
                                : metal_cell_state(*c.metal).nH;
}
double& cell_T_K(Cell& c) noexcept {
  return c.model == Model::Prim ? prim_cell_state(*c.prim).T_K
                                : metal_cell_state(*c.metal).T_K;
}
double& cell_mu(Cell& c) noexcept {
  return c.model == Model::Prim ? prim_cell_state(*c.prim).mu
                                : metal_cell_state(*c.metal).mu;
}
double& cell_gamma(Cell& c) noexcept {
  return c.model == Model::Prim ? prim_cell_state(*c.prim).gamma
                                : metal_cell_state(*c.metal).gamma;
}

ChemFullRates chem_full_step(Cell& cell, double dt, const ChemParams& params,
                             const ChemShielding& shield, const Table& tbl) {
  if (cell.model != tbl.model)
    throw std::invalid_argument(
        "arche::chem_full_step: cell and table belong to different models");
  if (cell.model == Model::Prim)
    return chem_full_step_prim(*cell.prim, dt, params, shield, *tbl.prim);
  return chem_full_step_metal(*cell.metal, dt, params, shield, *tbl.metal);
}

// ---------------------------------------------------------------------------
// Explicit template instantiation.
//
// Emits the full kernel object code for both models here, in libarche, with
// external linkage.  This is what lets consumer translation units avoid
// re-instantiating the kernel (they call the non-template facade above, or may
// use `extern template` against these symbols).
// ---------------------------------------------------------------------------
template ChemFullRates chem_full_step<Nakauchi2019>(ZeroMetalCell&, double,
                                                    const ChemParams&,
                                                    const ChemShielding&,
                                                    const ZeroMetalTable&);
template ChemFullRates chem_full_step<Nakauchi2019_Minimal>(
    ChemCell<Nakauchi2019_Minimal>&, double, const ChemParams&,
    const ChemShielding&, const zero_metal_minimal::MinimalTable&);
template ChemFullRates chem_full_step<Nakauchi2021>(MetalGrainCell&, double,
                                                    const ChemParams&,
                                                    const ChemShielding&,
                                                    const MetalGrainTable&);
template ChemRates chem_step<Nakauchi2019>(ZeroMetalCell&, double,
                                           const ChemParams&,
                                           const ZeroMetalTable&);
template ChemRates chem_step<Nakauchi2019_Minimal>(
    ChemCell<Nakauchi2019_Minimal>&, double, const ChemParams&,
    const zero_metal_minimal::MinimalTable&);
template ChemRates chem_step<Nakauchi2021>(MetalGrainCell&, double,
                                           const ChemParams&,
                                           const MetalGrainTable&);
template ChemFullRates chem_full_step<Nakauchi2021_Minimal>(
    ChemCell<Nakauchi2021_Minimal>&, double, const ChemParams&,
    const ChemShielding&, const metal_grain_minimal::MinimalTable&);
template ChemRates chem_step<Nakauchi2021_Minimal>(
    ChemCell<Nakauchi2021_Minimal>&, double, const ChemParams&,
    const metal_grain_minimal::MinimalTable&);

// ---------------------------------------------------------------------------
// Registry-backed type-erased models (name-selected; species metadata).
//
// IModelRuntime / IModelCell erase the compile-time model type behind a vtable;
// one concrete ModelRuntimeImpl<Model> per registered model holds its loaded
// reaction table and forwards stepping to the same chem_full_step<Model> path
// instantiated above (so the runtime-selected path is numerically identical to
// the per-model facade).  The registry is a single (name, factory) table —
// a new model is one row.
// ---------------------------------------------------------------------------
namespace {

struct IModelCell {
  virtual ~IModelCell() = default;
  virtual double* y() noexcept = 0;
  virtual double& nH() noexcept = 0;
  virtual double& T_K() noexcept = 0;
  virtual double& mu() noexcept = 0;
  virtual double& gamma() noexcept = 0;
};

struct IModelRuntime {
  virtual ~IModelRuntime() = default;
  virtual const char* name() const noexcept = 0;
  virtual int n_species() const noexcept = 0;
  virtual const char* const* species() const noexcept = 0;
  virtual std::unique_ptr<IModelCell> make_cell() const = 0;
  virtual ChemFullRates step(IModelCell& c, double dt, const ChemParams& p,
                             const ChemShielding& s) const = 0;
};

template <class Model>
struct ModelCellImpl final : IModelCell {
  ChemCell<Model> cell;
  double* y() noexcept override { return cell.state.y.data(); }
  double& nH() noexcept override { return cell.state.nH; }
  double& T_K() noexcept override { return cell.state.T_K; }
  double& mu() noexcept override { return cell.state.mu; }
  double& gamma() noexcept override { return cell.state.gamma; }
};

template <class Model>
struct ModelRuntimeImpl final : IModelRuntime {
  using Table = ReactionTable<Model::N_sp, Model::N_react>;
  std::string nm;
  Table tbl;
  std::array<const char*, Model::N_sp> names;

  ModelRuntimeImpl(std::string name, Table&& t)
      : nm(std::move(name)), tbl(std::move(t)) {
    // Host-facing species names in this model's local index order.
    for (int i = 0; i < Model::N_sp; ++i)
      names[i] = species_name(Model::Species::canonical(i));
  }
  const char* name() const noexcept override { return nm.c_str(); }
  int n_species() const noexcept override { return Model::N_sp; }
  const char* const* species() const noexcept override { return names.data(); }
  std::unique_ptr<IModelCell> make_cell() const override {
    return std::make_unique<ModelCellImpl<Model>>();
  }
  ChemFullRates step(IModelCell& c, double dt, const ChemParams& p,
                     const ChemShielding& s) const override {
    auto* mc = dynamic_cast<ModelCellImpl<Model>*>(&c);
    if (!mc)
      throw std::invalid_argument(
          "arche::model_step: cell was not created from this model");
    return chem_full_step<Model>(mc->cell, dt, p, s, tbl);
  }
};

using Factory = std::unique_ptr<IModelRuntime> (*)(const std::string& h5_path);

std::unique_ptr<IModelRuntime> make_prim(const std::string& h5) {
  return std::make_unique<ModelRuntimeImpl<Nakauchi2019>>(
      "Nakauchi2019", make_zero_metal_table(h5));
}
std::unique_ptr<IModelRuntime> make_prim_minimal(const std::string& h5) {
  return std::make_unique<ModelRuntimeImpl<Nakauchi2019_Minimal>>(
      "Nakauchi2019_Minimal", make_minimal_table(h5));
}
std::unique_ptr<IModelRuntime> make_metal(const std::string& h5) {
  return std::make_unique<ModelRuntimeImpl<Nakauchi2021>>(
      "Nakauchi2021", make_metal_grain_table(h5));
}
std::unique_ptr<IModelRuntime> make_metal_minimal(const std::string& h5) {
  // The compact table carries the full PF-loaded table (aux_full_metal) for the
  // strategy-1 rate path and the compact Saha.
  return std::make_unique<ModelRuntimeImpl<Nakauchi2021_Minimal>>(
      "Nakauchi2021_Minimal",
      make_metal_minimal_table(make_metal_grain_table(h5)));
}

// The model registry: name -> factory, in declaration order.  A new model is
// exposed to every host (C++ and the C ABI) by adding one row here.
const std::vector<std::pair<std::string, Factory>>& registry() {
  static const std::vector<std::pair<std::string, Factory>> r = {
      {"Nakauchi2019", &make_prim},
      {"Nakauchi2019_Minimal", &make_prim_minimal},
      {"Nakauchi2021", &make_metal},
      {"Nakauchi2021_Minimal", &make_metal_minimal},
  };
  return r;
}

}  // namespace

struct ModelRuntime {
  std::unique_ptr<IModelRuntime> impl;
};
struct ModelCell {
  std::unique_ptr<IModelCell> impl;
};

ModelRuntime* model_create(const std::string& name,
                           const std::string& h5_path) {
  for (const auto& [nm, factory] : registry()) {
    if (nm == name) {
      auto impl = factory(h5_path);  // may throw before anything is allocated
      auto* m = new ModelRuntime();
      m->impl = std::move(impl);
      return m;
    }
  }
  throw std::invalid_argument("arche::model_create: unknown model name '" +
                              name + "'");
}
void model_destroy(ModelRuntime* m) noexcept { delete m; }

int model_n_species(const ModelRuntime& m) noexcept {
  return m.impl->n_species();
}
const char* const* model_species(const ModelRuntime& m) noexcept {
  return m.impl->species();
}
const char* model_name(const ModelRuntime& m) noexcept {
  return m.impl->name();
}
int model_count() noexcept { return static_cast<int>(registry().size()); }
const char* model_registry_name(int i) noexcept {
  if (i < 0 || i >= static_cast<int>(registry().size())) return nullptr;
  return registry()[static_cast<std::size_t>(i)].first.c_str();
}

ModelCell* model_cell_create(const ModelRuntime& m) {
  auto* c = new ModelCell();
  c->impl = m.impl->make_cell();
  return c;
}
void model_cell_destroy(ModelCell* c) noexcept { delete c; }
double* model_cell_y(ModelCell& c) noexcept { return c.impl->y(); }
const double* model_cell_y(const ModelCell& c) noexcept { return c.impl->y(); }
double& model_cell_nH(ModelCell& c) noexcept { return c.impl->nH(); }
double& model_cell_T_K(ModelCell& c) noexcept { return c.impl->T_K(); }
double& model_cell_mu(ModelCell& c) noexcept { return c.impl->mu(); }
double& model_cell_gamma(ModelCell& c) noexcept { return c.impl->gamma(); }

ChemFullRates model_step(const ModelRuntime& m, ModelCell& c, double dt,
                         const ChemParams& params,
                         const ChemShielding& shield) {
  return m.impl->step(*c.impl, dt, params, shield);
}

}  // namespace arche
