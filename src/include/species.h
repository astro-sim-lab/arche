// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
#include <stdexcept>
#include <string>
#include <string_view>

namespace chemistry {

namespace abundance_ref {
constexpr double yHe = 8.33e-2;   // He/H number ratio
constexpr double yD  = 2.58e-5;   // D/H
constexpr double yLi = 4.65e-10;  // Li/H
constexpr double yC  = 2.69e-4;   // C/H  (solar)
constexpr double yO  = 4.90e-4;   // O/H  (solar)
constexpr double yNa = 1.74e-6;   // Na/H (solar)
constexpr double yMg = 3.98e-5;   // Mg/H (solar)
constexpr double yK  = 1.07e-7;   // K/H  (solar)
} // namespace abundance_ref

// ---------------------------------------------------------------------------
// zero_metal model  (N_sp=23, N_react=140)
// ---------------------------------------------------------------------------
namespace zero_metal {

constexpr int N_sp    = 23;
constexpr int N_react = 140;

// Species indices (1-based in Fortran → 0-based here)
enum Sp : int {
    H    =  0,   // y(1)
    H2   =  1,   // y(2)
    e    =  2,   // y(3)   electron
    Hp   =  3,   // y(4)   H+
    H2p  =  4,   // y(5)   H2+
    H3p  =  5,   // y(6)   H3+
    Hm   =  6,   // y(7)   H-
    He   =  7,   // y(8)
    Hep  =  8,   // y(9)   He+
    Hepp =  9,   // y(10)  He++
    HeHp = 10,   // y(11)  HeH+
    D    = 11,   // y(12)
    HD   = 12,   // y(13)
    Dp   = 13,   // y(14)  D+
    HDp  = 14,   // y(15)  HD+
    Dm   = 15,   // y(16)  D-
    Li   = 16,   // y(17)
    LiH  = 17,   // y(18)
    Lip  = 18,   // y(19)  Li+
    Lim  = 19,   // y(20)  Li-
    LiHp = 20,   // y(21)  LiH+
    Lipp = 21,   // y(22)  Li++
    Lippp= 22,   // y(23)  Li+++
};

// Cosmic abundances (used for initial conditions and conservation checks)
constexpr double yHe  = abundance_ref::yHe;  // He/H number ratio
constexpr double yD   = abundance_ref::yD;   // D/H
constexpr double yLi  = abundance_ref::yLi;  // Li/H

} // namespace zero_metal

// ---------------------------------------------------------------------------
// metal_grain model  (N_sp=89, N_react=1200)
// ---------------------------------------------------------------------------
namespace metal_grain {

constexpr int N_sp    = 89;
constexpr int N_react = 1200;

// Species indices (Fortran 1-based → C++ 0-based).
// Species 0–15 are identical to zero_metal.
// Species 16+ diverge: metal/grain species replace Li.
enum Sp : int {
    // H/He/D network (same layout as zero_metal 0-15)
    H     =  0,   // y(1)
    H2    =  1,   // y(2)
    e     =  2,   // y(3)  electron
    Hp    =  3,   // y(4)  H+
    H2p   =  4,   // y(5)  H2+
    H3p   =  5,   // y(6)  H3+
    Hm    =  6,   // y(7)  H-
    He    =  7,   // y(8)
    Hep   =  8,   // y(9)  He+
    Hepp  =  9,   // y(10) He++
    HeHp  = 10,   // y(11) HeH+
    D     = 11,   // y(12)
    HD    = 12,   // y(13)
    Dp    = 13,   // y(14) D+
    HDp   = 14,   // y(15) HD+
    Dm    = 15,   // y(16) D-
    // Carbon species (C/CH/C2 family)
    C     = 16,   // y(17) CI
    C2    = 17,   // y(18)
    CH    = 18,   // y(19)
    CH2   = 19,   // y(20)
    CH3   = 20,   // y(21)
    CH4   = 21,   // y(22)
    Cp    = 22,   // y(23) C+  (CII)
    C2p   = 23,   // y(24) C2+
    CHp   = 24,   // y(25) CH+
    CH2p  = 25,   // y(26) CH2+
    CH3p  = 26,   // y(27) CH3+
    CH4p  = 27,   // y(28) CH4+
    CH5p  = 28,   // y(29) CH5+
    // Oxygen species (O/OH/CO/H2O family)
    O     = 29,   // y(30) OI
    O2    = 30,   // y(31)
    OH    = 31,   // y(32)
    CO    = 32,   // y(33)
    H2O   = 33,   // y(34)
    HCO   = 34,   // y(35)
    O2H   = 35,   // y(36) HO2
    CO2   = 36,   // y(37)
    H2CO  = 37,   // y(38)
    H2O2  = 38,   // y(39)
    Op    = 39,   // y(40) O+
    O2p   = 40,   // y(41) O2+
    OHp   = 41,   // y(42) OH+
    COp   = 42,   // y(43) CO+
    H2Op  = 43,   // y(44) H2O+
    HCOp  = 44,   // y(45) HCO+
    O2Hp  = 45,   // y(46) O2H+
    H3Op  = 46,   // y(47) H3O+
    H2COp = 47,   // y(48) H2CO+
    HOCOp = 48,   // y(49) HOCO+
    H2COHp= 49,   // y(50) H2COH+
    // Li species
    Li    = 50,   // y(51)
    LiH   = 51,   // y(52)
    Lip   = 52,   // y(53) Li+
    Lim   = 53,   // y(54) Li-
    LiHp  = 54,   // y(55) LiH+
    Lipp  = 55,   // y(56) Li++
    Lippp = 56,   // y(57) Li+++
    // K species
    K     = 57,   // y(58)
    Kp    = 58,   // y(59) K+
    // Na species
    Na    = 59,   // y(60)
    Nap   = 60,   // y(61) Na+
    // Mg species
    Mg    = 61,   // y(62)
    Mgp   = 62,   // y(63) Mg+
    // Dust grain charge states
    Gr    = 63,   // y(64) neutral grain
    Grp   = 64,   // y(65) grain+
    Gr2p  = 65,   // y(66) grain2+
    Grm   = 66,   // y(67) grain-
    Gr2m  = 67,   // y(68) grain2-
    // Ice mantle / grain surface species (physisorbed/chemisorbed)
    H_p   = 68,   // y(69) H  physisorbed on grain
    H_c   = 69,   // y(70) H  chemisorbed on grain
    H2_p  = 70,   // y(71) H2 physisorbed on grain
    D_p   = 71,   // y(72) D  physisorbed
    D_c   = 72,   // y(73) D  chemisorbed
    HD_p  = 73,   // y(74) HD physisorbed
    O_p   = 74,   // y(75) O  on grain
    O2_p  = 75,   // y(76) O2 on grain
    OH_p  = 76,   // y(77) OH on grain
    CO_p  = 77,   // y(78) CO on grain
    CO2_p = 78,   // y(79) CO2 on grain
    H2O_p = 79,   // y(80) H2O on grain  ← xJH2O source
    HO2_p = 80,   // y(81) HO2 on grain
    H2O2_p= 81,   // y(82) H2O2 on grain
    HCO_p = 82,   // y(83) HCO on grain
    H2CO_p= 83,   // y(84) H2CO on grain
    C_p   = 84,   // y(85) C  on grain
    CH_p  = 85,   // y(86) CH on grain
    CH2_p = 86,   // y(87) CH2 on grain
    CH3_p = 87,   // y(88) CH3 on grain
    CH4_p = 88,   // y(89) CH4 on grain
};
// xcrit mapping (C++ 0-based):
//   xJH2  = y[H2_p  = 70]   (Fortran y(71))
//   xJH2O = y[H2O_p = 79]   (Fortran y(80))
//   xJtot = y[H_c   = 69] + y[D_c = 72]  (Fortran y(70)+y(73))

// Cosmic and solar abundances
constexpr double yHe  = abundance_ref::yHe;  // He/H
constexpr double yD   = abundance_ref::yD;   // D/H
constexpr double yLi  = abundance_ref::yLi;  // Li/H
constexpr double yC   = abundance_ref::yC;   // C/H  (solar)
constexpr double yO   = abundance_ref::yO;   // O/H  (solar)
constexpr double yNa  = abundance_ref::yNa;  // Na/H (solar)
constexpr double yMg  = abundance_ref::yMg;  // Mg/H (solar)
constexpr double yK   = abundance_ref::yK;   // K/H  (solar)

// Default broken-MRN grain size distribution (used by metal-grain model)
constexpr double mrn_a_min = 5.0e-7;   // [cm]
constexpr double mrn_a_mid = 1.0e-4;   // [cm]
constexpr double mrn_a_max = 5.0e-4;   // [cm]
constexpr double mrn_ind1  = 3.5;      // n(a) ∝ a^-ind1 for a<a_mid
constexpr double mrn_ind2  = 5.5;      // n(a) ∝ a^-ind2 for a>a_mid
constexpr int    mrn_nint1 = 15;       // quadrature bins (small-grain segment)
constexpr int    mrn_nint2 = 5;        // quadrature bins (large-grain segment)
constexpr double mrn_norm_zsun = 5.333e-3; // dust mass fraction per Z_sun

// Default initial gas-phase fractions (metal run)
constexpr double c_gas_frac_default  = 0.28;
constexpr double o_gas_frac_default  = 0.54;
constexpr double mg_gas_frac_default = 0.02;

} // namespace metal_grain

// ---------------------------------------------------------------------------
// Photon/dummy index used in the Fortran reaction table
// Fortran encodes "no reactant/product" as index 0 or 101 in 1-based arrays.
// In C++ we map: index 0 → NONE, index 101 → PHOTON (kept in y_dum with value 1)
// ---------------------------------------------------------------------------
constexpr int IDX_NONE   = 0;    // vacant reactant/product slot
constexpr int IDX_PHOTON = 101;  // photon (y_dum[101] = 1 in Fortran)

// ---------------------------------------------------------------------------
// Abundance-set presets (runtime selection in app layer)
//
// Current policy:
// - Keep a single physically vetted baseline (`solar`).
// - Expose a preset selector now so additional sets can be added later
//   without changing app interfaces.
// ---------------------------------------------------------------------------
namespace abundance {

inline std::string to_lower_ascii(std::string_view s)
{
    std::string out;
    out.reserve(s.size());
    for (const char c : s) {
        if (c >= 'A' && c <= 'Z') out.push_back(static_cast<char>(c - 'A' + 'a'));
        else out.push_back(c);
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

inline PrimordialSet get_primordial_set(std::string_view preset)
{
    const std::string p = to_lower_ascii(preset);
    if (p.empty() || p == "solar" || p == "default" || p == "alpha-enhanced") {
        return PrimordialSet{
            (p == "default") ? "default" : (p == "alpha-enhanced" ? "alpha-enhanced" : "solar"),
            zero_metal::yHe,
            zero_metal::yD,
            zero_metal::yLi
        };
    }
    throw std::invalid_argument(
        "Unsupported primordial abundance preset: '" + std::string(preset) +
        "' (allowed: solar, default, alpha-enhanced)");
}

inline MetalSet get_metal_set(std::string_view preset)
{
    const std::string p = to_lower_ascii(preset);
    if (p.empty() || p == "solar" || p == "default") {
        return MetalSet{
            (p == "default") ? "default" : "solar",
            metal_grain::yHe,
            metal_grain::yD,
            metal_grain::yLi,
            metal_grain::yC,
            metal_grain::yO,
            metal_grain::yNa,
            metal_grain::yMg,
            metal_grain::yK
        };
    }
    if (p == "alpha-enhanced") {
        // Sample preset for sensitivity tests (not a uniquely recommended
        // physical standard): alpha elements (O, Mg) enhanced by +0.4 dex
        // relative to solar, others unchanged.
        constexpr double f_alpha = 2.51188643150958; // 10^0.4
        return MetalSet{
            "alpha-enhanced",
            metal_grain::yHe,
            metal_grain::yD,
            metal_grain::yLi,
            metal_grain::yC,
            metal_grain::yO  * f_alpha,
            metal_grain::yNa,
            metal_grain::yMg * f_alpha,
            metal_grain::yK
        };
    }
    throw std::invalid_argument(
        "Unsupported metal abundance preset: '" + std::string(preset) +
        "' (allowed: solar, default, alpha-enhanced)");
}

} // namespace abundance

} // namespace chemistry
