// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// species_catalog.h — canonical species catalog shared by all models.
//
// A single place that names every species (primordial, metal, grain, and grain
// surface) and carries its per-species metadata (display name, mass).  Masses
// live here and nowhere else: a model selects a subset through
// SpeciesSet<...>, which yields a compact 0..N-1 local index together with the
// canonical<->local maps.  A model's own species enum, dimension N_sp, and mass
// array are then derived from its SpeciesSet rather than hand-numbered, so the
// local index space is a natural consequence of the catalog selection instead
// of a per-model special case.
//
//   SpId            — canonical species identifier (catalog key).
//   SpeciesInfo     — per-species metadata row.
//   species_info()  — canonical metadata lookup (name, mass).
//   SpeciesSet<...> — a model's species selection: N, local(), canonical().
//   make_mass_array — compact mass array in a SpeciesSet's local order.
// ---------------------------------------------------------------------------
#include <array>
#include <cstddef>

namespace arche {

// ---------------------------------------------------------------------------
// SpId — canonical identifiers for every species across all models.
//
// The primordial block (H..Lippp) leads in canonical reference order; a model
// that selects the whole primordial block in this order reproduces the
// historical zero_metal index layout.  The metal_grain block follows; each
// model's SpeciesSet selects and orders its own subset (see
// core/species_index.h).
// ---------------------------------------------------------------------------
enum class SpId : int {
  H,
  H2,
  e,
  Hp,
  H2p,
  H3p,
  Hm,
  He,
  Hep,
  Hepp,
  HeHp,
  D,
  HD,
  Dp,
  HDp,
  Dm,
  Li,
  LiH,
  Lip,
  Lim,
  LiHp,
  Lipp,
  Lippp,
  // ── metal_grain extension (declared in metal_grain Sp order) ──
  // Carbon family
  C,
  C2,
  CH,
  CH2,
  CH3,
  CH4,
  Cp,
  C2p,
  CHp,
  CH2p,
  CH3p,
  CH4p,
  CH5p,
  // Oxygen family
  O,
  O2,
  OH,
  CO,
  H2O,
  HCO,
  O2H,
  CO2,
  H2CO,
  H2O2,
  Op,
  O2p,
  OHp,
  COp,
  H2Op,
  HCOp,
  O2Hp,
  H3Op,
  H2COp,
  HOCOp,
  H2COHp,
  // K / Na / Mg
  K,
  Kp,
  Na,
  Nap,
  Mg,
  Mgp,
  // Dust grain charge states
  Gr,
  Grp,
  Gr2p,
  Grm,
  Gr2m,
  // Ice mantle / grain surface species (physisorbed/chemisorbed)
  H_p,
  H_c,
  H2_p,
  D_p,
  D_c,
  HD_p,
  O_p,
  O2_p,
  OH_p,
  CO_p,
  CO2_p,
  H2O_p,
  HO2_p,
  H2O2_p,
  HCO_p,
  H2CO_p,
  C_p,
  CH_p,
  CH2_p,
  CH3_p,
  CH4_p,
  Count
};

// Per-species metadata row.
struct SpeciesInfo {
  SpId id;
  const char* name;  // display name (e.g. host-facing species labels)
  double mass;       // [amu]
};

// ---------------------------------------------------------------------------
// Canonical catalog, indexed by SpId value (kept in SpId declaration order).
// Masses are the species masses in amu (consumers multiply by phys::m_p to
// obtain grams).
// ---------------------------------------------------------------------------
inline constexpr SpeciesInfo kSpeciesCatalog[] = {
    {SpId::H, "H", 1.00783},
    {SpId::H2, "H2", 2.01565},
    {SpId::e, "e", 0.00055},
    {SpId::Hp, "H+", 1.00728},
    {SpId::H2p, "H2+", 2.0151},
    {SpId::H3p, "H3+", 3.02293},
    {SpId::Hm, "H-", 1.00837},
    {SpId::He, "He", 4.0026},
    {SpId::Hep, "He+", 4.00205},
    {SpId::Hepp, "He++", 4.0},
    {SpId::HeHp, "HeH+", 5.00988},
    {SpId::D, "D", 2.0141},
    {SpId::HD, "HD", 3.02204},
    {SpId::Dp, "D+", 2.01355},
    {SpId::HDp, "HD+", 3.02149},
    {SpId::Dm, "D-", 2.01465},
    {SpId::Li, "Li", 6.941},
    {SpId::LiH, "LiH", 7.949},
    {SpId::Lip, "Li+", 6.94},
    {SpId::Lim, "Li-", 6.942},
    {SpId::LiHp, "LiH+", 7.949},
    {SpId::Lipp, "Li++", 6.94},
    {SpId::Lippp, "Li+++", 6.94},
    // ── metal_grain extension ──
    // Carbon family
    {SpId::C, "C", 12.0},
    {SpId::C2, "C2", 24.0},
    {SpId::CH, "CH", 13.00783},
    {SpId::CH2, "CH2", 14.01565},
    {SpId::CH3, "CH3", 15.02348},
    {SpId::CH4, "CH4", 16.0313},
    {SpId::Cp, "C+", 11.99945},
    {SpId::C2p, "C2+", 23.99945},
    {SpId::CHp, "CH+", 13.00728},
    {SpId::CH2p, "CH2+", 14.0151},
    {SpId::CH3p, "CH3+", 15.02293},
    {SpId::CH4p, "CH4+", 16.03075},
    {SpId::CH5p, "CH5+", 17.03858},
    // Oxygen family
    {SpId::O, "O", 15.99491},
    {SpId::O2, "O2", 31.98983},
    {SpId::OH, "OH", 17.00274},
    {SpId::CO, "CO", 27.99491},
    {SpId::H2O, "H2O", 18.01056},
    {SpId::HCO, "HCO", 29.00274},
    {SpId::O2H, "O2H", 32.99765},
    {SpId::CO2, "CO2", 43.98983},
    {SpId::H2CO, "H2CO", 30.01056},
    {SpId::H2O2, "H2O2", 34.00548},
    {SpId::Op, "O+", 15.99437},
    {SpId::O2p, "O2+", 31.98928},
    {SpId::OHp, "OH+", 17.00219},
    {SpId::COp, "CO+", 27.99437},
    {SpId::H2Op, "H2O+", 18.01002},
    {SpId::HCOp, "HCO+", 29.00219},
    {SpId::O2Hp, "O2H+", 32.99711},
    {SpId::H3Op, "H3O+", 19.01784},
    {SpId::H2COp, "H2CO+", 30.01002},
    {SpId::HOCOp, "HOCO+", 44.99711},
    {SpId::H2COHp, "H2COH+", 31.01784},
    // K / Na / Mg
    {SpId::K, "K", 39.0983},
    {SpId::Kp, "K+", 39.09775},
    {SpId::Na, "Na", 22.98977},
    {SpId::Nap, "Na+", 22.98922},
    {SpId::Mg, "Mg", 24.305},
    {SpId::Mgp, "Mg+", 24.30445},
    // Dust grain charge states (mass placeholder; not used by any live kernel)
    {SpId::Gr, "Gr", 1000000.0},
    {SpId::Grp, "Gr+", 1000000.0},
    {SpId::Gr2p, "Gr2+", 1000000.0},
    {SpId::Grm, "Gr-", 1000000.0},
    {SpId::Gr2m, "Gr2-", 1000000.0},
    // Ice mantle / grain surface species (mass mirrors the gas-phase parent)
    {SpId::H_p, "H_p", 1.00783},
    {SpId::H_c, "H_c", 1.00783},
    {SpId::H2_p, "H2_p", 2.01565},
    {SpId::D_p, "D_p", 2.0141},
    {SpId::D_c, "D_c", 2.0141},
    {SpId::HD_p, "HD_p", 3.02204},
    {SpId::O_p, "O_p", 15.99491},
    {SpId::O2_p, "O2_p", 31.98983},
    {SpId::OH_p, "OH_p", 17.00274},
    {SpId::CO_p, "CO_p", 27.99491},
    {SpId::CO2_p, "CO2_p", 43.98983},
    {SpId::H2O_p, "H2O_p", 18.01056},
    {SpId::HO2_p, "HO2_p", 32.99765},
    {SpId::H2O2_p, "H2O2_p", 34.00548},
    {SpId::HCO_p, "HCO_p", 29.00274},
    {SpId::H2CO_p, "H2CO_p", 30.01056},
    {SpId::C_p, "C_p", 12.0},
    {SpId::CH_p, "CH_p", 13.00783},
    {SpId::CH2_p, "CH2_p", 14.01565},
    {SpId::CH3_p, "CH3_p", 15.02348},
    {SpId::CH4_p, "CH4_p", 16.0313},
};

static_assert(std::size(kSpeciesCatalog) ==
                  static_cast<std::size_t>(SpId::Count),
              "kSpeciesCatalog must hold exactly one row per SpId");

// Canonical metadata lookup.  The catalog is stored in SpId order, so this is a
// direct index; a debug build can additionally assert id consistency.
inline constexpr const SpeciesInfo& species_info(SpId id) {
  return kSpeciesCatalog[static_cast<int>(id)];
}
inline constexpr const char* species_name(SpId id) {
  return species_info(id).name;
}
inline constexpr double species_mass(SpId id) { return species_info(id).mass; }

// ---------------------------------------------------------------------------
// SpeciesSet — a model's species selection.
//
// Sel... lists the canonical species the model carries, in the order they
// occupy local indices 0..N-1.  local()/canonical() are the compile-time maps
// between canonical SpId and the model's compact index space.
// ---------------------------------------------------------------------------
template <SpId... Sel>
struct SpeciesSet {
  static constexpr int N = sizeof...(Sel);
  static constexpr std::array<SpId, N> ids{{Sel...}};

  // canonical -> local (0..N-1), or -1 if the species is not in the set.
  static constexpr int local(SpId id) {
    for (int i = 0; i < N; ++i)
      if (ids[i] == id) return i;
    return -1;
  }
  // local -> canonical.
  static constexpr SpId canonical(int i) { return ids[i]; }

  static constexpr const char* name(int i) { return species_name(ids[i]); }
  static constexpr double mass(int i) { return species_mass(ids[i]); }
};

// Compact mass array [amu] in a SpeciesSet's local order (0..N-1).  Every mass
// consumer (e.g. the metal_grain grain-charge and grain-surface kernels)
// derives its array from this, so masses live only in the catalog.
template <class Set>
inline constexpr std::array<double, Set::N> make_mass_array() {
  std::array<double, Set::N> m{};
  for (int i = 0; i < Set::N; ++i) m[i] = species_mass(Set::canonical(i));
  return m;
}

}  // namespace arche
