// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// species_index.h — public species metadata for a network model.
//
// This is the consumer-facing slice of the species definitions: it answers
// "which species exist, in what y[] order, and what are their reference
// abundances".  It carries ONLY what an external application needs to fill the
// per-cell state vector by name:
//
//   * abundance_ref::                — cosmic/solar reference number ratios
//   * <model>::Sp                    — y[] index enums (e.g. zero_metal::H == 0)
//   * <model>::N_sp / N_react        — network sizes
//   * <model>:: abundance constants  — per-model cosmic/solar ratios
//   * abundance::                    — runtime abundance-preset selection
//
// The reaction-rate bookkeeping that the kernel uses internally (forward-rate
// slot allocation, cooling/heating hot-spot indices, shared reaction indices,
// detailed-balance loop bounds and the HDF5-row-layout static_asserts) is NOT
// here — it lives in the internal header kinetics/reaction_index.h, which is
// not part of the installed/public surface.  Keeping it out means a consumer's
// public contract does not move when the rate kernel's internal layout changes.
// ---------------------------------------------------------------------------
#include <stdexcept>
#include <string>
#include <string_view>

#include "core/species_catalog.h"

namespace arche {

namespace abundance_ref {
constexpr double yHe = 8.33e-2;   // He/H number ratio
constexpr double yD = 2.58e-5;    // D/H
constexpr double yLi = 4.65e-10;  // Li/H
constexpr double yC = 2.69e-4;    // C/H  (solar)
constexpr double yO = 4.90e-4;    // O/H  (solar)
constexpr double yNa = 1.74e-6;   // Na/H (solar)
constexpr double yMg = 3.98e-5;   // Mg/H (solar)
constexpr double yK = 1.07e-7;    // K/H  (solar)
}  // namespace abundance_ref

// ---------------------------------------------------------------------------
// zero_metal model  (N_sp=23, N_react=140)
// ---------------------------------------------------------------------------
namespace zero_metal {

constexpr int N_react = 140;

// Species selection from the canonical catalog.  Selecting the whole catalog in
// canonical order reproduces the historical 0..22 index layout; N_sp, the Sp
// enum, and the mass array are all derived from this set rather than
// hand-numbered (see core/species_catalog.h).
using Species =
    SpeciesSet<SpId::H, SpId::H2, SpId::e, SpId::Hp, SpId::H2p, SpId::H3p,
               SpId::Hm, SpId::He, SpId::Hep, SpId::Hepp, SpId::HeHp, SpId::D,
               SpId::HD, SpId::Dp, SpId::HDp, SpId::Dm, SpId::Li, SpId::LiH,
               SpId::Lip, SpId::Lim, SpId::LiHp, SpId::Lipp, SpId::Lippp>;

constexpr int N_sp = Species::N;

// Species indices (0-based), derived from the catalog selection above.
enum Sp : int {
  H = Species::local(SpId::H),
  H2 = Species::local(SpId::H2),
  e = Species::local(SpId::e),      // electron
  Hp = Species::local(SpId::Hp),    // H+
  H2p = Species::local(SpId::H2p),  // H2+
  H3p = Species::local(SpId::H3p),  // H3+
  Hm = Species::local(SpId::Hm),    // H-
  He = Species::local(SpId::He),
  Hep = Species::local(SpId::Hep),    // He+
  Hepp = Species::local(SpId::Hepp),  // He++
  HeHp = Species::local(SpId::HeHp),  // HeH+
  D = Species::local(SpId::D),
  HD = Species::local(SpId::HD),
  Dp = Species::local(SpId::Dp),    // D+
  HDp = Species::local(SpId::HDp),  // HD+
  Dm = Species::local(SpId::Dm),    // D-
  Li = Species::local(SpId::Li),
  LiH = Species::local(SpId::LiH),
  Lip = Species::local(SpId::Lip),      // Li+
  Lim = Species::local(SpId::Lim),      // Li-
  LiHp = Species::local(SpId::LiHp),    // LiH+
  Lipp = Species::local(SpId::Lipp),    // Li++
  Lippp = Species::local(SpId::Lippp),  // Li+++
};

// Layout lock: the full primordial model uses the catalog's canonical order, so
// the derived indices must reproduce the historical 23-species layout.
static_assert(N_sp == 23, "zero_metal must carry 23 species");
static_assert(H == 0 && H2 == 1 && e == 2 && He == 7 && D == 11 && Li == 16 &&
                  Lippp == 22,
              "zero_metal Sp indices must match the canonical catalog order");

// Cosmic abundances (used for initial conditions and conservation checks)
constexpr double yHe = abundance_ref::yHe;  // He/H number ratio
constexpr double yD = abundance_ref::yD;    // D/H
constexpr double yLi = abundance_ref::yLi;  // Li/H

}  // namespace zero_metal

// ---------------------------------------------------------------------------
// metal_grain model  (N_sp=89, N_react=1200)
// ---------------------------------------------------------------------------
namespace metal_grain {

constexpr int N_react = 1200;

// Species selection from the canonical catalog, in the historical metal_grain
// order (local indices 0..88).  N_sp, the Sp enum, and the grain-kernel mass
// arrays are all derived from this set rather than hand-numbered (see
// core/species_catalog.h).  Species 0-15 match zero_metal; 16+ add the
// metal/grain/surface species (Li is relocated to 50-56).
// clang-format off
using Species = SpeciesSet<
    SpId::H, SpId::H2, SpId::e, SpId::Hp, SpId::H2p, SpId::H3p, SpId::Hm,
    SpId::He, SpId::Hep, SpId::Hepp, SpId::HeHp, SpId::D, SpId::HD, SpId::Dp,
    SpId::HDp, SpId::Dm, SpId::C, SpId::C2, SpId::CH, SpId::CH2, SpId::CH3,
    SpId::CH4, SpId::Cp, SpId::C2p, SpId::CHp, SpId::CH2p, SpId::CH3p,
    SpId::CH4p, SpId::CH5p, SpId::O, SpId::O2, SpId::OH, SpId::CO, SpId::H2O,
    SpId::HCO, SpId::O2H, SpId::CO2, SpId::H2CO, SpId::H2O2, SpId::Op,
    SpId::O2p, SpId::OHp, SpId::COp, SpId::H2Op, SpId::HCOp, SpId::O2Hp,
    SpId::H3Op, SpId::H2COp, SpId::HOCOp, SpId::H2COHp, SpId::Li, SpId::LiH,
    SpId::Lip, SpId::Lim, SpId::LiHp, SpId::Lipp, SpId::Lippp, SpId::K,
    SpId::Kp, SpId::Na, SpId::Nap, SpId::Mg, SpId::Mgp, SpId::Gr, SpId::Grp,
    SpId::Gr2p, SpId::Grm, SpId::Gr2m, SpId::H_p, SpId::H_c, SpId::H2_p,
    SpId::D_p, SpId::D_c, SpId::HD_p, SpId::O_p, SpId::O2_p, SpId::OH_p,
    SpId::CO_p, SpId::CO2_p, SpId::H2O_p, SpId::HO2_p, SpId::H2O2_p,
    SpId::HCO_p, SpId::H2CO_p, SpId::C_p, SpId::CH_p, SpId::CH2_p,
    SpId::CH3_p, SpId::CH4_p>;

// clang-format on
constexpr int N_sp = Species::N;

// Species indices (0-based), derived from the catalog selection above.
// clang-format off
enum Sp : int {
  // H/He/D network (same layout as zero_metal 0-15)
  H = Species::local(SpId::H),
  H2 = Species::local(SpId::H2),
  e = Species::local(SpId::e),      // electron
  Hp = Species::local(SpId::Hp),    // H+
  H2p = Species::local(SpId::H2p),  // H2+
  H3p = Species::local(SpId::H3p),  // H3+
  Hm = Species::local(SpId::Hm),    // H-
  He = Species::local(SpId::He),
  Hep = Species::local(SpId::Hep),    // He+
  Hepp = Species::local(SpId::Hepp),  // He++
  HeHp = Species::local(SpId::HeHp),  // HeH+
  D = Species::local(SpId::D),
  HD = Species::local(SpId::HD),
  Dp = Species::local(SpId::Dp),    // D+
  HDp = Species::local(SpId::HDp),  // HD+
  Dm = Species::local(SpId::Dm),    // D-
  // Carbon species (C/CH/C2 family)
  C = Species::local(SpId::C),  // CI
  C2 = Species::local(SpId::C2),
  CH = Species::local(SpId::CH),
  CH2 = Species::local(SpId::CH2),
  CH3 = Species::local(SpId::CH3),
  CH4 = Species::local(SpId::CH4),
  Cp = Species::local(SpId::Cp),      // C+  (CII)
  C2p = Species::local(SpId::C2p),    // C2+
  CHp = Species::local(SpId::CHp),    // CH+
  CH2p = Species::local(SpId::CH2p),  // CH2+
  CH3p = Species::local(SpId::CH3p),  // CH3+
  CH4p = Species::local(SpId::CH4p),  // CH4+
  CH5p = Species::local(SpId::CH5p),  // CH5+
  // Oxygen species (O/OH/CO/H2O family)
  O = Species::local(SpId::O),  // OI
  O2 = Species::local(SpId::O2),
  OH = Species::local(SpId::OH),
  CO = Species::local(SpId::CO),
  H2O = Species::local(SpId::H2O),
  HCO = Species::local(SpId::HCO),
  O2H = Species::local(SpId::O2H),  // HO2
  CO2 = Species::local(SpId::CO2),
  H2CO = Species::local(SpId::H2CO),
  H2O2 = Species::local(SpId::H2O2),
  Op = Species::local(SpId::Op),        // O+
  O2p = Species::local(SpId::O2p),      // O2+
  OHp = Species::local(SpId::OHp),      // OH+
  COp = Species::local(SpId::COp),      // CO+
  H2Op = Species::local(SpId::H2Op),    // H2O+
  HCOp = Species::local(SpId::HCOp),    // HCO+
  O2Hp = Species::local(SpId::O2Hp),    // O2H+
  H3Op = Species::local(SpId::H3Op),    // H3O+
  H2COp = Species::local(SpId::H2COp),  // H2CO+
  HOCOp = Species::local(SpId::HOCOp),  // HOCO+
  H2COHp = Species::local(SpId::H2COHp),  // H2COH+
  // Li species
  Li = Species::local(SpId::Li),
  LiH = Species::local(SpId::LiH),
  Lip = Species::local(SpId::Lip),      // Li+
  Lim = Species::local(SpId::Lim),      // Li-
  LiHp = Species::local(SpId::LiHp),    // LiH+
  Lipp = Species::local(SpId::Lipp),    // Li++
  Lippp = Species::local(SpId::Lippp),  // Li+++
  // K species
  K = Species::local(SpId::K),
  Kp = Species::local(SpId::Kp),  // K+
  // Na species
  Na = Species::local(SpId::Na),
  Nap = Species::local(SpId::Nap),  // Na+
  // Mg species
  Mg = Species::local(SpId::Mg),
  Mgp = Species::local(SpId::Mgp),  // Mg+
  // Dust grain charge states
  Gr = Species::local(SpId::Gr),      // neutral grain
  Grp = Species::local(SpId::Grp),    // grain+
  Gr2p = Species::local(SpId::Gr2p),  // grain2+
  Grm = Species::local(SpId::Grm),    // grain-
  Gr2m = Species::local(SpId::Gr2m),  // grain2-
  // Ice mantle / grain surface species (physisorbed/chemisorbed)
  H_p = Species::local(SpId::H_p),    // H  physisorbed on grain
  H_c = Species::local(SpId::H_c),    // H  chemisorbed on grain
  H2_p = Species::local(SpId::H2_p),  // H2 physisorbed on grain
  D_p = Species::local(SpId::D_p),    // D  physisorbed
  D_c = Species::local(SpId::D_c),    // D  chemisorbed
  HD_p = Species::local(SpId::HD_p),  // HD physisorbed
  O_p = Species::local(SpId::O_p),    // O  on grain
  O2_p = Species::local(SpId::O2_p),  // O2 on grain
  OH_p = Species::local(SpId::OH_p),  // OH on grain
  CO_p = Species::local(SpId::CO_p),  // CO on grain
  CO2_p = Species::local(SpId::CO2_p),  // CO2 on grain
  H2O_p = Species::local(SpId::H2O_p),  // H2O on grain  ← JH2O source
  HO2_p = Species::local(SpId::HO2_p),    // HO2 on grain
  H2O2_p = Species::local(SpId::H2O2_p),  // H2O2 on grain
  HCO_p = Species::local(SpId::HCO_p),    // HCO on grain
  H2CO_p = Species::local(SpId::H2CO_p),  // H2CO on grain
  C_p = Species::local(SpId::C_p),    // C  on grain
  CH_p = Species::local(SpId::CH_p),  // CH on grain
  CH2_p = Species::local(SpId::CH2_p),  // CH2 on grain
  CH3_p = Species::local(SpId::CH3_p),  // CH3 on grain
  CH4_p = Species::local(SpId::CH4_p),  // CH4 on grain
};
// clang-format on

// Layout lock: the derived indices must reproduce the historical 89-species
// metal_grain layout used by the reaction tables and grain kernels.
static_assert(N_sp == 89, "metal_grain must carry 89 species");
static_assert(
    H == 0 && H2 == 1 && e == 2 && Dm == 15 && C == 16 && O == 29 &&
        H2COHp == 49 && Li == 50 && Lippp == 56 && K == 57 && Mgp == 62 &&
        Gr == 63 && Gr2m == 67 && H_p == 68 && H_c == 69 && H2_p == 70 &&
        D_c == 72 && H2O_p == 79 && CH4_p == 88,
    "metal_grain Sp indices must match the historical catalog layout");
// crit mapping (C++ 0-based):
//   JH2  = y[H2_p  = 70]
//   JH2O = y[H2O_p = 79]
//   Jtot = y[H_c   = 69] + y[D_c = 72]

// Cosmic and solar abundances
constexpr double yHe = abundance_ref::yHe;  // He/H
constexpr double yD = abundance_ref::yD;    // D/H
constexpr double yLi = abundance_ref::yLi;  // Li/H
constexpr double yC = abundance_ref::yC;    // C/H  (solar)
constexpr double yO = abundance_ref::yO;    // O/H  (solar)
constexpr double yNa = abundance_ref::yNa;  // Na/H (solar)
constexpr double yMg = abundance_ref::yMg;  // Mg/H (solar)
constexpr double yK = abundance_ref::yK;    // K/H  (solar)

// Default broken-MRN grain size distribution (used by metal-grain model)
constexpr double mrn_a_min = 5.0e-7;  // [cm]
constexpr double mrn_a_mid = 1.0e-4;  // [cm]
constexpr double mrn_a_max = 5.0e-4;  // [cm]
constexpr double mrn_ind1 = 3.5;      // n(a) ∝ a^-ind1 for a<a_mid
constexpr double mrn_ind2 = 5.5;      // n(a) ∝ a^-ind2 for a>a_mid
constexpr int mrn_nint1 = 15;         // quadrature bins (small-grain segment)
constexpr int mrn_nint2 = 5;          // quadrature bins (large-grain segment)
constexpr double mrn_norm_zsun = 5.333e-3;  // dust mass fraction per Z_sun

// Default initial gas-phase fractions (metal run)
constexpr double c_gas_frac_default = 0.28;
constexpr double o_gas_frac_default = 0.54;
constexpr double mg_gas_frac_default = 0.02;

}  // namespace metal_grain

// ---------------------------------------------------------------------------
// Abundance-set presets (runtime selection in app layer)
//
// Current policy:
// - Keep a single physically vetted baseline (`solar`).
// - Expose a preset selector now so additional sets can be added later
//   without changing app interfaces.
// ---------------------------------------------------------------------------
namespace abundance {

inline std::string to_lower_ascii(std::string_view s) {
  std::string out;
  out.reserve(s.size());
  for (const char c : s) {
    if (c >= 'A' && c <= 'Z')
      out.push_back(static_cast<char>(c - 'A' + 'a'));
    else
      out.push_back(c);
  }
  return out;
}

struct PrimordialSet {
  const char* name;
  double yHe;
  double yD;
  double yLi;
};

struct MetalSet {
  const char* name;
  double yHe;
  double yD;
  double yLi;
  double yC;
  double yO;
  double yNa;
  double yMg;
  double yK;
};

inline PrimordialSet get_primordial_set(std::string_view preset) {
  const std::string p = to_lower_ascii(preset);
  if (p.empty() || p == "solar" || p == "default" || p == "alpha-enhanced") {
    return PrimordialSet{
        (p == "default") ? "default"
                         : (p == "alpha-enhanced" ? "alpha-enhanced" : "solar"),
        zero_metal::yHe, zero_metal::yD, zero_metal::yLi};
  }
  throw std::invalid_argument("Unsupported primordial abundance preset: '" +
                              std::string(preset) +
                              "' (allowed: solar, default, alpha-enhanced)");
}

inline MetalSet get_metal_set(std::string_view preset) {
  const std::string p = to_lower_ascii(preset);
  if (p.empty() || p == "solar" || p == "default") {
    return MetalSet{(p == "default") ? "default" : "solar",
                    metal_grain::yHe,
                    metal_grain::yD,
                    metal_grain::yLi,
                    metal_grain::yC,
                    metal_grain::yO,
                    metal_grain::yNa,
                    metal_grain::yMg,
                    metal_grain::yK};
  }
  if (p == "alpha-enhanced") {
    // Sample preset for sensitivity tests (not a uniquely recommended
    // physical standard): alpha elements (O, Mg) enhanced by +0.4 dex
    // relative to solar, others unchanged.
    constexpr double f_alpha = 2.51188643150958;  // 10^0.4
    return MetalSet{"alpha-enhanced", metal_grain::yHe,
                    metal_grain::yD,  metal_grain::yLi,
                    metal_grain::yC,  metal_grain::yO * f_alpha,
                    metal_grain::yNa, metal_grain::yMg * f_alpha,
                    metal_grain::yK};
  }
  throw std::invalid_argument("Unsupported metal abundance preset: '" +
                              std::string(preset) +
                              "' (allowed: solar, default, alpha-enhanced)");
}

}  // namespace abundance

}  // namespace arche
