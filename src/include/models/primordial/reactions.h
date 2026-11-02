// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
#include <array>
#include <cstddef>

#include "core/species_index.h"
#include "kinetics/topology.h"

namespace arche {
namespace zero_metal {

// ---------------------------------------------------------------------------
// zero_metal reaction network — hand-maintained source of truth.
//
// Bootstrapped once from the legacy HDF5 topology and now edited directly
// here.  Species are
// referenced through the Sp:: enum so an index slip is a compile error.
// PHOTON / CR are sentinels; unused slots are auto-padded with VACANT and
// n_reactants / n_products are inferred (count of real species) by the
// rxn()/saha_rxn() helpers.  Cmass/delE are stored verbatim.
//
// Detailed balance.  Each reversible network row's delE and Cmass are the ones
// implied by the Saha formation reactions below (kSaha): delE is the difference
// of the species energies those rows define, and Cmass the matching
// (Π m_reactants / Π m_products)^{3/2}, to the Saha table's 6 significant
// digits.  The kinetic network's equilibrium then coincides with the Saha
// branch's.  With 3 digits (Cmass of ~1 rounded to 1.0) the network's equilibrium
// sits up to ~10% (H+) away from Saha at a fixed temperature.  Rows with Cmass
// 0.0 (forward-only), CR rows and delE = 0 rows keep their stored values.
// ---------------------------------------------------------------------------

namespace net {

using Sp = zero_metal::Sp;
inline constexpr int N_sp = zero_metal::N_sp;
inline constexpr int PHOTON = N_sp + 1;  // photon sentinel
inline constexpr int CR = N_sp + 2;      // cosmic-ray sentinel

// Thin wrappers binding the generic builders to this model's N_sp.
constexpr Reaction rxn(int num, std::initializer_list<int> reactants,
                       std::initializer_list<int> products, double Cmass,
                       double delE) {
  return topo::make_reaction(N_sp, num, reactants, products, Cmass, delE);
}
constexpr SahaReaction saha_rxn(int num, std::initializer_list<int> reactants,
                                std::initializer_list<int> products,
                                double Cmass, double delE) {
  return topo::make_saha(N_sp, num, reactants, products, Cmass, delE);
}

// Standard / CR / charge-transfer reactions (139 rows).
inline constexpr std::array<Reaction, 139> kReactions = {{
    rxn(1, {Sp::H, Sp::e}, {Sp::e, Sp::e, Sp::Hp}, 3.58686e+40,
        -2.17871e-11),  // H + e- -> e- + e- + H+
    rxn(2, {Sp::e, Sp::Hp}, {Sp::H, PHOTON}, 2.78795e-41,
        2.17871e-11),  // e- + H+ -> H + γ
    rxn(3, {Sp::e, Sp::He}, {Sp::e, Sp::e, Sp::Hep}, 3.58468e+40,
        -3.93933e-11),  // e- + He -> e- + e- + He+
    rxn(4, {Sp::e, Sp::Hep}, {Sp::He, PHOTON}, 2.78965e-41,
        3.93933e-11),  // e- + He+ -> He + γ
    rxn(5, {Sp::e, Sp::Hep}, {Sp::e, Sp::e, Sp::Hepp}, 3.58668e+40,
        -8.71859e-11),  // e- + He+ -> e- + e- + He++
    rxn(6, {Sp::e, Sp::Hepp}, {Sp::Hep, PHOTON}, 2.78809e-41,
        8.71859e-11),  // e- + He++ -> He+ + γ
    rxn(7, {Sp::H, Sp::e}, {Sp::Hm, PHOTON}, 2.78799e-41,
        1.20867e-12),  // H + e- -> H- + γ
    rxn(8, {Sp::H, Sp::Hm}, {Sp::H2, Sp::e}, 27755.2,
        5.96598e-12),  // H + H- -> H2 + e-
    rxn(9, {Sp::H, Sp::Hp}, {Sp::H2p, PHOTON}, 7.73495e-37,
        4.24688e-12),  // H + H+ -> H2+ + γ
    rxn(10, {Sp::H, Sp::H2p}, {Sp::H2, Sp::Hp}, 1.00041,
        2.92777e-12),  // H + H2+ -> H2 + H+
    rxn(11, {Sp::H2, Sp::Hp}, {Sp::H, Sp::H2p}, 0.99959,
        -2.92777e-12),  // H2 + H+ -> H + H2+
    rxn(12, {Sp::H2, Sp::e}, {Sp::H, Sp::H, Sp::e}, 1.2923e+36,
        -7.17465e-12),  // H2 + e- -> H + H + e-
    rxn(13, {Sp::H, Sp::H2}, {Sp::H, Sp::H, Sp::H}, 1.2923e+36,
        -7.17465e-12),  // H + H2 -> H + H + H
    rxn(14, {Sp::e, Sp::Hm}, {Sp::H, Sp::e, Sp::e}, 3.58681e+40,
        -1.20867e-12),  // e- + H- -> H + e- + e-
    rxn(15, {Sp::Hp, Sp::Hm}, {Sp::H, Sp::H}, 0.999986,
        2.05784e-11),  // H+ + H- -> H + H
    rxn(16, {Sp::Hp, Sp::Hm}, {Sp::e, Sp::H2p}, 27743.8,
        3.03821e-12),  // H+ + H- -> e- + H2+
    rxn(17, {Sp::e, Sp::H2p}, {Sp::H, Sp::H}, 3.60435e-05,
        1.75402e-11),  // e- + H2+ -> H + H
    rxn(18, {Sp::H2p, Sp::Hm}, {Sp::H, Sp::H2}, 1.0004,
        2.35062e-11),  // H2+ + H- -> H + H2
    rxn(19, {Sp::H, Sp::H, Sp::H}, {Sp::H, Sp::H2}, 7.73812e-37,
        7.17465e-12),  // H + H + H -> H + H2
    rxn(20, {Sp::H, Sp::H, Sp::H2}, {Sp::H2, Sp::H2}, 7.73812e-37,
        7.17465e-12),  // H + H + H2 -> H2 + H2
    rxn(21, {Sp::H2, Sp::H2}, {Sp::H, Sp::H, Sp::H2}, 1.2923e+36,
        -7.17465e-12),  // H2 + H2 -> H + H + H2
    rxn(22, {Sp::H, Sp::H}, {Sp::H, Sp::e, Sp::Hp}, 3.58686e+40,
        -2.17871e-11),          // H + H -> H + e- + H+
    rxn(23, {}, {}, 1.0, 0.0),  //  ->
    rxn(24, {Sp::H2, Sp::Hep}, {Sp::H, Sp::Hp, Sp::He}, 1.29309e+36,
        1.04315e-11),  // H2 + He+ -> H + H+ + He
    rxn(25, {Sp::H2p, Sp::He}, {Sp::H, Sp::HeHp}, 2.019,
        -1.34603e-12),  // H2+ + He -> H + HeH+
    rxn(26, {Sp::H2, Sp::H2p}, {Sp::H, Sp::H3p}, 1.53938,
        2.69052e-12),  // H2 + H2+ -> H + H3+
    rxn(27, {Sp::H3p, Sp::Hm}, {Sp::H2, Sp::H2}, 0.649869,
        2.08157e-11),  // H3+ + H- -> H2 + H2
    rxn(28, {Sp::H, Sp::Hep}, {Sp::Hp, Sp::He}, 1.00061,
        1.76062e-11),  // H + He+ -> H+ + He
    rxn(29, {Sp::Hm, Sp::Hep}, {Sp::H, Sp::He}, 1.0006,
        3.81846e-11),  // H- + He+ -> H + He
    rxn(30, {Sp::H2, Sp::Hep}, {Sp::H2p, Sp::He}, 1.0002,
        1.46784e-11),  // H2 + He+ -> H2+ + He
    rxn(31, {Sp::H, Sp::HeHp}, {Sp::H2p, Sp::He}, 0.495294,
        1.34603e-12),  // H + HeH+ -> H2+ + He
    rxn(32, {Sp::H2, Sp::HeHp}, {Sp::H3p, Sp::He}, 0.762445,
        4.03655e-12),  // H2 + HeH+ -> H3+ + He
    rxn(33, {Sp::H2p, Sp::Hm}, {Sp::e, Sp::H3p}, 42725.8,
        8.6565e-12),  // H2+ + H- -> e- + H3+
    rxn(34, {Sp::e, Sp::H3p}, {Sp::H, Sp::H2}, 2.34143e-05,
        1.48497e-11),  // e- + H3+ -> H + H2
    rxn(35, {Sp::e, Sp::H3p}, {Sp::H, Sp::H, Sp::H}, 3.02584e+31,
        7.67505e-12),  // e- + H3+ -> H + H + H
    rxn(36, {Sp::e, Sp::HeHp}, {Sp::H, Sp::He}, 1.78521e-05,
        1.88862e-11),  // e- + HeH+ -> H + He
    rxn(37, {Sp::H2, Sp::He}, {Sp::H, Sp::H, Sp::He}, 1.2923e+36,
        -7.17465e-12),  // H2 + He -> H + H + He
    rxn(38, {Sp::H, Sp::Hm}, {Sp::H, Sp::H, Sp::e}, 3.58681e+40,
        -1.20867e-12),  // H + H- -> H + H + e-
    rxn(39, {Sp::H2p, Sp::Hm}, {Sp::H, Sp::H, Sp::H}, 1.29281e+36,
        1.63315e-11),  // H2+ + H- -> H + H + H
    rxn(40, {Sp::H2, Sp::e}, {Sp::H, Sp::Hm}, 3.60293e-05,
        -5.96598e-12),  // H2 + e- -> H + H-
    rxn(41, {Sp::Hp, Sp::He}, {Sp::H, Sp::Hep}, 0.999391,
        -1.76062e-11),  // H+ + He -> H + He+
    rxn(42, {Sp::Hm, Sp::He}, {Sp::H, Sp::e, Sp::He}, 3.58681e+40,
        -1.20867e-12),  // H- + He -> H + e- + He
    rxn(43, {Sp::H, Sp::H, Sp::He}, {Sp::H2, Sp::He}, 7.73812e-37,
        7.17465e-12),  // H + H + He -> H2 + He
    rxn(44, {Sp::H, Sp::H2p}, {Sp::H3p, PHOTON}, 1.19119e-36,
        9.86517e-12),  // H + H2+ -> H3+ + γ
    rxn(45, {Sp::H2, Sp::Hp}, {Sp::H3p, PHOTON}, 1.1907e-36,
        6.9374e-12),  // H2 + H+ -> H3+ + γ
    rxn(46, {Sp::Hp, Sp::He}, {Sp::HeHp, PHOTON}, 1.56169e-36,
        2.90085e-12),  // H+ + He -> HeH+ + γ
    rxn(47, {Sp::e, Sp::HD}, {Sp::H, Sp::e, Sp::D}, 8.3975e+35,
        -7.23181e-12),  // e- + HD -> H + e- + D
    rxn(48, {Sp::He, Sp::HD}, {Sp::H, Sp::He, Sp::D}, 8.3975e+35,
        -7.23181e-12),  // He + HD -> H + He + D
    rxn(49, {Sp::H2, Sp::HD}, {Sp::H, Sp::H2, Sp::D}, 8.3975e+35,
        -7.23181e-12),  // H2 + HD -> H + H2 + D
    rxn(50, {Sp::H, Sp::HD}, {Sp::H, Sp::H, Sp::D}, 8.3975e+35,
        -7.23181e-12),  // H + HD -> H + H + D
    rxn(51, {Sp::e, Sp::Dp}, {Sp::D, PHOTON}, 2.78909e-41,
        2.1793e-11),  // e- + D+ -> D + γ
    rxn(52, {Sp::Hp, Sp::D}, {Sp::H, Sp::Dp}, 0.999591,
        -5.9e-15),  // H+ + D -> H + D+
    rxn(53, {Sp::H, Sp::Dp}, {Sp::Hp, Sp::D}, 1.00041,
        5.9e-15),  // H + D+ -> H+ + D
    rxn(54, {Sp::H, Sp::D}, {Sp::HD, PHOTON}, 1.19083e-36,
        7.23181e-12),  // H + D -> HD + γ
    rxn(55, {Sp::H2, Sp::D}, {Sp::H, Sp::HD}, 1.53891,
        5.716e-14),  // H2 + D -> H + HD
    rxn(56, {Sp::H, Sp::HDp}, {Sp::Hp, Sp::HD}, 1.00055,
        2.95775e-12),  // H + HD+ -> H+ + HD
    rxn(57, {Sp::H2, Sp::Dp}, {Sp::Hp, Sp::HD}, 1.53954,
        6.306e-14),  // H2 + D+ -> H+ + HD
    rxn(58, {Sp::H, Sp::HD}, {Sp::H2, Sp::D}, 0.649809,
        -5.716e-14),  // H + HD -> H2 + D
    rxn(59, {Sp::Hp, Sp::HD}, {Sp::H2, Sp::Dp}, 0.649544,
        -6.306e-14),  // H+ + HD -> H2 + D+
    rxn(60, {Sp::Hp, Sp::D}, {Sp::HDp, PHOTON}, 1.19018e-36,
        4.27406e-12),  // H+ + D -> HD+ + γ
    rxn(61, {Sp::H, Sp::Dp}, {Sp::HDp, PHOTON}, 1.19067e-36,
        4.27996e-12),  // H + D+ -> HD+ + γ
    rxn(62, {Sp::e, Sp::HDp}, {Sp::H, Sp::D}, 2.34246e-05,
        1.7513e-11),  // e- + HD+ -> H + D
    rxn(63, {Sp::e, Sp::D}, {Sp::Dm, PHOTON}, 2.78909e-41,
        1.20909e-12),  // e- + D -> D- + γ
    rxn(64, {Sp::Dp, Sp::Dm}, {Sp::D, Sp::D}, 1.0,
        2.05839e-11),  // D+ + D- -> D + D
    rxn(65, {Sp::Hp, Sp::Dm}, {Sp::H, Sp::D}, 0.999591,
        2.0578e-11),  // H+ + D- -> H + D
    rxn(66, {Sp::Hm, Sp::D}, {Sp::H, Sp::Dm}, 1.00039,
        4.2e-16),  // H- + D -> H + D-
    rxn(67, {Sp::H, Sp::Dm}, {Sp::Hm, Sp::D}, 0.999606,
        -4.2e-16),  // H + D- -> H- + D
    rxn(68, {Sp::H, Sp::Dm}, {Sp::e, Sp::HD}, 42696.0,
        6.02272e-12),  // H + D- -> e- + HD
    rxn(69, {Sp::e, Sp::D}, {Sp::e, Sp::e, Sp::Dp}, 3.5854e+40,
        -2.1793e-11),  // e- + D -> e- + e- + D+
    rxn(70, {Sp::Hep, Sp::D}, {Sp::He, Sp::Dp}, 1.0002,
        1.76003e-11),  // He+ + D -> He + D+
    rxn(71, {Sp::He, Sp::Dp}, {Sp::Hep, Sp::D}, 0.999799,
        -1.76003e-11),  // He + D+ -> He+ + D
    rxn(72, {Sp::H2p, Sp::D}, {Sp::H, Sp::HDp}, 1.5387,
        2.718e-14),  // H2+ + D -> H + HD+
    rxn(73, {Sp::D, Sp::HDp}, {Sp::HD, Sp::Dp}, 1.00014,
        2.95185e-12),  // D + HD+ -> HD + D+
    rxn(74, {Sp::H, Sp::HDp}, {Sp::H2p, Sp::D}, 0.649897,
        -2.718e-14),  // H + HD+ -> H2+ + D
    rxn(75, {Sp::e, Sp::HD}, {Sp::H, Sp::Dm}, 2.34214e-05,
        -6.02272e-12),  // e- + HD -> H + D-
    rxn(76, {Sp::e, Sp::HD}, {Sp::Hm, Sp::D}, 2.34122e-05,
        -6.02314e-12),  // e- + HD -> H- + D
    rxn(77, {Sp::Hp, Sp::Dm}, {Sp::e, Sp::HDp}, 42672.7,
        3.06497e-12),  // H+ + D- -> e- + HD+
    rxn(78, {Sp::Hm, Sp::Dp}, {Sp::e, Sp::HDp}, 42707.0,
        3.07129e-12),  // H- + D+ -> e- + HD+
    rxn(79, {Sp::e, Sp::Dm}, {Sp::e, Sp::e, Sp::D}, 3.5854e+40,
        -1.20909e-12),  // e- + D- -> e- + e- + D
    rxn(80, {Sp::H, Sp::Dm}, {Sp::H, Sp::e, Sp::D}, 3.5854e+40,
        -1.20909e-12),  // H + D- -> H + e- + D
    rxn(81, {Sp::He, Sp::Dm}, {Sp::e, Sp::He, Sp::D}, 3.5854e+40,
        -1.20909e-12),  // He + D- -> e- + He + D
    rxn(82, {Sp::Hm, Sp::Dp}, {Sp::H, Sp::D}, 1.00039,
        2.05843e-11),  // H- + D+ -> H + D
    rxn(83, {Sp::H2p, Sp::Dm}, {Sp::H2, Sp::D}, 1.0,
        2.35058e-11),  // H2+ + D- -> H2 + D
    rxn(84, {Sp::H2p, Sp::Dm}, {Sp::H, Sp::H, Sp::D}, 1.2923e+36,
        1.63311e-11),  // H2+ + D- -> H + H + D
    rxn(85, {Sp::Hm, Sp::HDp}, {Sp::H, Sp::HD}, 1.00053,
        2.35362e-11),  // H- + HD+ -> H + HD
    rxn(86, {Sp::Hm, Sp::HDp}, {Sp::H, Sp::H, Sp::D}, 8.40197e+35,
        1.63044e-11),  // H- + HD+ -> H + H + D
    rxn(87, {Sp::HDp, Sp::Dm}, {Sp::D, Sp::HD}, 1.00014,
        2.35358e-11),  // HD+ + D- -> D + HD
    rxn(88, {Sp::HDp, Sp::Dm}, {Sp::H, Sp::D, Sp::D}, 8.39866e+35,
        1.63039e-11),  // HD+ + D- -> H + D + D
    rxn(89, {Sp::Hep, Sp::Dm}, {Sp::He, Sp::D}, 1.0002,
        3.81842e-11),  // He+ + D- -> He + D
    rxn(90, {Sp::H2p, Sp::D}, {Sp::H2, Sp::Dp}, 1.0,
        2.92187e-12),  // H2+ + D -> H2 + D+
    rxn(91, {Sp::H2p, Sp::D}, {Sp::Hp, Sp::HD}, 1.53954,
        2.98493e-12),  // H2+ + D -> H+ + HD
    rxn(92, {Sp::H, Sp::HDp}, {Sp::H2, Sp::Dp}, 0.649898,
        2.89469e-12),  // H + HD+ -> H2 + D+
    rxn(93, {Sp::H2, Sp::Dp}, {Sp::H2p, Sp::D}, 0.999999,
        -2.92187e-12),  // H2 + D+ -> H2+ + D
    rxn(94, {Sp::H2, Sp::Dp}, {Sp::H, Sp::HDp}, 1.5387,
        -2.89469e-12),  // H2 + D+ -> H + HD+
    rxn(95, {Sp::Hp, Sp::HD}, {Sp::H, Sp::HDp}, 0.999454,
        -2.95775e-12),  // H+ + HD -> H + HD+
    rxn(96, {Sp::Hp, Sp::HD}, {Sp::H2p, Sp::D}, 0.649543,
        -2.98493e-12),  // H+ + HD -> H2+ + D
    rxn(97, {Sp::HD, Sp::Dp}, {Sp::D, Sp::HDp}, 0.999863,
        -2.95185e-12),  // HD + D+ -> D + HD+
    rxn(98, {Sp::Hep, Sp::HD}, {Sp::He, Sp::HDp}, 1.00006,
        1.46484e-11),  // He+ + HD -> He + HD+
    rxn(99, {Sp::Hep, Sp::HD}, {Sp::Hp, Sp::He, Sp::D}, 8.40262e+35,
        1.03744e-11),  // He+ + HD -> H+ + He + D
    rxn(100, {Sp::Hep, Sp::HD}, {Sp::H, Sp::He, Sp::Dp}, 8.39919e+35,
        1.03685e-11),  // He+ + HD -> H + He + D+
    rxn(101, {Sp::e, Sp::Lip}, {Sp::Li, PHOTON}, 2.78963e-41,
        8.63849e-12),  // e- + Li+ -> Li + γ
    rxn(102, {Sp::Hm, Sp::Lip}, {Sp::H, Sp::Li}, 1.00059,
        7.42982e-12),  // H- + Li+ -> H + Li
    rxn(103, {Sp::Hp, Sp::Lim}, {Sp::H, Sp::Li}, 0.999398,
        2.0797e-11),  // H+ + Li- -> H + Li
    rxn(104, {Sp::e, Sp::Li}, {Sp::Lim, PHOTON}, 2.78963e-41,
        9.90146e-13),  // e- + Li -> Li- + γ
    rxn(105, {Sp::Hp, Sp::Li}, {Sp::H, Sp::Lip}, 0.999398,
        1.31486e-11),  // H+ + Li -> H + Li+
    rxn(106, {Sp::Hp, Sp::Li}, {Sp::H, Sp::Lip, PHOTON}, 0.999398,
        1.31486e-11),  // H+ + Li -> H + Li+ + γ
    rxn(107, {Sp::Hm, Sp::Li}, {Sp::e, Sp::LiH}, 64054.4,
        2.6598e-12),  // H- + Li -> e- + LiH
    rxn(108, {Sp::H, Sp::Lim}, {Sp::e, Sp::LiH}, 64016.7,
        2.87832e-12),  // H + Li- -> e- + LiH
    rxn(109, {Sp::H, Sp::LiHp}, {Sp::Hp, Sp::LiH}, 1.00082,
        -9.51233e-12),  // H + LiH+ -> H+ + LiH
    rxn(110, {Sp::Hp, Sp::LiH}, {Sp::H, Sp::LiHp}, 0.999182,
        9.51233e-12),  // H+ + LiH -> H + LiH+
    rxn(111, {Sp::H, Sp::LiH}, {Sp::H2, Sp::Li}, 0.433307,
        3.30618e-12),  // H + LiH -> H2 + Li
    rxn(112, {Sp::H, Sp::Lip}, {Sp::LiHp, PHOTON}, 1.78545e-36,
        2.3219e-13),  // H + Li+ -> LiH+ + γ
    rxn(113, {Sp::Hp, Sp::Li}, {Sp::LiHp, PHOTON}, 1.78437e-36,
        1.33808e-11),  // H+ + Li -> LiH+ + γ
    rxn(114, {Sp::Hp, Sp::LiH}, {Sp::H2, Sp::Lip}, 0.433046,
        1.64548e-11),  // H+ + LiH -> H2 + Li+
    rxn(115, {Sp::e, Sp::LiHp}, {Sp::H, Sp::Li}, 1.56243e-05,
        8.4063e-12),  // e- + LiH+ -> H + Li
    rxn(116, {Sp::H, Sp::LiHp}, {Sp::H2p, Sp::Li}, 0.433484,
        -9.13392e-12),  // H + LiH+ -> H2+ + Li
    rxn(117, {Sp::H, Sp::LiHp}, {Sp::H2, Sp::Lip}, 0.4334,
        6.94246e-12),  // H + LiH+ -> H2 + Li+
    rxn(118, {Sp::H, Sp::Li}, {Sp::LiH, PHOTON}, 1.78583e-36,
        3.86847e-12),  // H + Li -> LiH + γ
    rxn(119, {Sp::H2p, Sp::Li}, {Sp::Hp, Sp::LiH}, 2.30878,
        -3.7841e-13),  // H2+ + Li -> H+ + LiH
    rxn(120, {Sp::Hp, Sp::LiH}, {Sp::H2p, Sp::Li}, 0.433129,
        3.7841e-13),  // H+ + LiH -> H2+ + Li
    rxn(121, {Sp::e, Sp::Lipp}, {Sp::Lip, PHOTON}, 2.79023e-41,
        1.21191e-10),  // e- + Li++ -> Li+ + γ
    rxn(122, {Sp::e, Sp::Lippp}, {Sp::Lipp, PHOTON}, 2.79023e-41,
        1.96196e-10),  // e- + Li+++ -> Li++ + γ
    rxn(123, {Sp::Dm, Sp::Lip}, {Sp::D, Sp::Li}, 1.00019,
        7.4294e-12),  // D- + Li+ -> D + Li
    rxn(124, {Sp::Dp, Sp::Lim}, {Sp::D, Sp::Li}, 0.999806,
        2.08029e-11),  // D+ + Li- -> D + Li
    rxn(125, {Sp::e, Sp::Li}, {Sp::e, Sp::e, Sp::Lip}, 3.5847e+40,
        -8.63849e-12),  // e- + Li -> e- + e- + Li+
    rxn(126, {Sp::e, Sp::Lip}, {Sp::e, Sp::e, Sp::Lipp}, 3.58393e+40,
        -1.21191e-10),  // e- + Li+ -> e- + e- + Li++
    rxn(127, {Sp::e, Sp::Lipp}, {Sp::e, Sp::e, Sp::Lippp}, 3.58393e+40,
        -1.96196e-10),  // e- + Li++ -> e- + e- + Li+++
    rxn(128, {Sp::H, Sp::H, Sp::Li}, {Sp::H, Sp::LiH}, 1.78583e-36,
        3.86847e-12),  // H + H + Li -> H + LiH
    rxn(129, {Sp::H, Sp::H2, Sp::Li}, {Sp::H2, Sp::LiH}, 1.78583e-36,
        3.86847e-12),  // H + H2 + Li -> H2 + LiH
    rxn(130, {Sp::H2, Sp::Li}, {Sp::H2, Sp::e, Sp::Lip}, 3.5847e+40,
        -8.63849e-12),  // H2 + Li -> H2 + e- + Li+
    rxn(131, {Sp::H, CR}, {Sp::e, Sp::Hp}, 0.0, 0.0),  // H + CR -> e- + H+
    rxn(132, {Sp::H2, CR}, {Sp::H, Sp::e, Sp::Hp}, 0.0,
        0.0),  // H2 + CR -> H + e- + H+
    rxn(133, {Sp::H2, CR}, {Sp::e, Sp::H2p}, 0.0, 0.0),  // H2 + CR -> e- + H2+
    rxn(134, {Sp::H2, CR}, {Sp::H, Sp::H}, 0.0, 0.0),    // H2 + CR -> H + H
    rxn(135, {Sp::H2, CR}, {Sp::Hp, Sp::Hm}, 0.0, 0.0),  // H2 + CR -> H+ + H-
    rxn(136, {Sp::He, CR}, {Sp::e, Sp::Hep}, 0.0, 0.0),  // He + CR -> e- + He+
    rxn(137, {Sp::H, CR}, {Sp::e, Sp::Hp}, 0.0, 0.0),    // H + CR -> e- + H+
    rxn(138, {Sp::Hm, CR}, {Sp::H, Sp::e}, 0.0, 0.0),    // H- + CR -> H + e-
    rxn(139, {Sp::He, CR}, {Sp::e, Sp::Hep}, 0.0, 0.0),  // He + CR -> e- + He+
}};

// ── Minimal-network keep-set (defines the compact Nakauchi2019_Minimal) ─────
// The compact Nakauchi2019_Minimal network is built from these reaction ids:
// models/primordial/minimal.h remaps the listed full-network reactions into the
// 15-species / 24-reaction compact table.  Retained set = the cosmic-ray-free
// minimal network (Nakauchi et al. 2019, Table 2's 18 reaction pairs) plus the
// collisional Li ionisation (id 130, which sources Li+ in the nH ~ 1e13–1e16
// band) and five ion-processing reactions (ids 4, 15, 17, 24, 38) that let the
// compact network better track the full network's ionization chemistry.  Of
// these, id 4 (He+ radiative recombination) and id 24 (H2 + He+ charge
// transfer) are He+ sinks; He+ itself is produced only by the cosmic-ray
// channels (ids 136/139, see kMinimalCRKeep).  The compact set omits the full
// network's collisional He ionization (id 3) and H/He charge transfer (ids 28,
// 41), so without cosmic rays He+ stays ~0 (the full network forms it
// collisionally at high T).  All are standard reversible reactions; the five
// appended ids keep slots 0..18 (and the cooling rate_idx map) unchanged.
inline constexpr std::array<int, 24> kMinimalKeep = {{
    2,    // e- + H+ -> H + γ
    7,    // H + e- -> H- + γ
    8,    // H + H- -> H2 + e-
    9,    // H + H+ -> H2+ + γ
    10,   // H + H2+ -> H2 + H+
    19,   // H + H + H -> H + H2          (3-body H2 formation)
    21,   // H2 + H2 -> H + H + H2
    26,   // H2 + H2+ -> H + H3+
    34,   // e- + H3+ -> H + H2
    35,   // e- + H3+ -> H + H + H
    52,   // H+ + D -> H + D+
    55,   // H2 + D -> H + HD
    57,   // H2 + D+ -> H+ + HD
    101,  // e- + Li+ -> Li + γ
    102,  // H- + Li+ -> H + Li
    106,  // H+ + Li -> H + Li+ + γ
    111,  // H + LiH -> H2 + Li
    118,  // H + Li -> LiH + γ
    130,  // H2 + Li -> H2 + e- + Li+     (collisional Li ionisation)
    4,    // e- + He+ -> He + γ           (He+ radiative recombination)
    15,   // H+ + H- -> H + H             (mutual neutralisation)
    17,   // e- + H2+ -> H + H            (dissociative recombination)
    24,   // H2 + He+ -> H + H+ + He      (He+ charge transfer)
    38,   // H + H- -> H + H + e-         (associative detachment)
}};

// ── Minimal-network CR keep-set ─────────────────────────────────────────────
// The 9 cosmic-ray (direct + CR-induced secondary photo) ionization channels
// of the full network, appended after the 24 standard reactions (compact slots
// 24..32) so the compact model carries the same CR source as the full network.
// Direct: 131..136 (H/H2/He + CR), CR-photo: 137..139.
inline constexpr std::array<int, 9> kMinimalCRKeep = {{
    131,  // H  + CR   -> e- + H+
    132,  // H2 + CR   -> H + e- + H+
    133,  // H2 + CR   -> e- + H2+
    134,  // H2 + CR   -> H + H
    135,  // H2 + CR   -> H+ + H-
    136,  // He + CR   -> e- + He+
    137,  // H  + CRph -> e- + H+
    138,  // H- + CRph -> H + e-
    139,  // He + CRph -> e- + He+
}};

// High-density Saha-equilibrium reactions (18 rows).
inline constexpr std::array<SahaReaction, 18> kSaha = {{
    saha_rxn(2, {Sp::e, Sp::Hp}, {Sp::H, PHOTON}, 2.78795e-41,
             2.17871e-11),  // e- + H+ -> H + γ
    saha_rxn(4, {Sp::e, Sp::Hep}, {Sp::He, PHOTON}, 2.78965e-41,
             3.93933e-11),  // e- + He+ -> He + γ
    saha_rxn(6, {Sp::e, Sp::Hepp}, {Sp::Hep, PHOTON}, 2.78809e-41,
             8.71859e-11),  // e- + He++ -> He+ + γ
    saha_rxn(7, {Sp::H, Sp::e}, {Sp::Hm, PHOTON}, 2.78799e-41,
             1.20867e-12),  // H + e- -> H- + γ
    saha_rxn(9, {Sp::H, Sp::Hp}, {Sp::H2p, PHOTON}, 7.73495e-37,
             4.24688e-12),  // H + H+ -> H2+ + γ
    saha_rxn(44, {Sp::Hp, Sp::He}, {Sp::HeHp, PHOTON}, 1.56169e-36,
             2.90085e-12),  // H+ + He -> HeH+ + γ
    saha_rxn(51, {Sp::e, Sp::Dp}, {Sp::D, PHOTON}, 2.78909e-41,
             2.1793e-11),  // e- + D+ -> D + γ
    saha_rxn(54, {Sp::H, Sp::D}, {Sp::HD, PHOTON}, 1.19083e-36,
             7.23181e-12),  // H + D -> HD + γ
    saha_rxn(60, {Sp::Hp, Sp::D}, {Sp::HDp, PHOTON}, 1.19018e-36,
             4.27406e-12),  // H+ + D -> HD+ + γ
    saha_rxn(63, {Sp::e, Sp::D}, {Sp::Dm, PHOTON}, 2.78909e-41,
             1.20909e-12),  // e- + D -> D- + γ
    saha_rxn(101, {Sp::e, Sp::Lip}, {Sp::Li, PHOTON}, 2.78963e-41,
             8.63849e-12),  // e- + Li+ -> Li + γ
    saha_rxn(104, {Sp::e, Sp::Li}, {Sp::Lim, PHOTON}, 2.78963e-41,
             9.90146e-13),  // e- + Li -> Li- + γ
    saha_rxn(113, {Sp::Hp, Sp::Li}, {Sp::LiHp, PHOTON}, 1.78437e-36,
             1.33808e-11),  // H+ + Li -> LiH+ + γ
    saha_rxn(118, {Sp::H, Sp::Li}, {Sp::LiH, PHOTON}, 1.78583e-36,
             3.86847e-12),  // H + Li -> LiH + γ
    saha_rxn(121, {Sp::e, Sp::Lipp}, {Sp::Lip, PHOTON}, 2.79023e-41,
             1.21191e-10),  // e- + Li++ -> Li+ + γ
    saha_rxn(122, {Sp::e, Sp::Lippp}, {Sp::Lipp, PHOTON}, 2.79023e-41,
             1.96196e-10),  // e- + Li+++ -> Li++ + γ
    saha_rxn(8, {Sp::H, Sp::Hm}, {Sp::H2, Sp::e}, 27755.2,
             5.96598e-12),  // H + H- -> H2 + e-
    saha_rxn(26, {Sp::H2, Sp::H2p}, {Sp::H, Sp::H3p}, 1.53938,
             2.69052e-12),  // H2 + H2+ -> H + H3+
}};

// Fill the topology fields of a ReactionTable from the constexpr network
// above.  Partition-function tables are loaded separately from HDF5.
// Species masses live only in core/species_catalog.h and are read directly by
// the kernels that need them (see make_mass_array).
inline void init_topology(ZeroMetalTable& tbl) {
  for (std::size_t i = 0; i < kReactions.size(); ++i)
    tbl.reactions[i] = kReactions[i];
  tbl.n_loaded = static_cast<int>(kReactions.size());
  for (std::size_t i = 0; i < kSaha.size(); ++i) tbl.saha[i] = kSaha[i];
  tbl.n_saha = static_cast<int>(kSaha.size());
  tbl.n_grain = 0;
}

}  // namespace net
}  // namespace zero_metal
}  // namespace arche
