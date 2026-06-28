// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
#include <array>
#include <cstddef>

#include "core/species_index.h"
#include "kinetics/topology.h"

namespace arche {
namespace metal_grain {

// ---------------------------------------------------------------------------
// metal_grain reaction network — hand-maintained source of truth.
//
// Bootstrapped once from the legacy HDF5 topology and now edited directly
// here.  Species are
// referenced through the Sp:: enum so an index slip is a compile error.
// PHOTON / CR are sentinels; unused slots are auto-padded with VACANT and
// n_reactants / n_products are inferred (count of real species) by the
// rxn()/saha_rxn() helpers.  Cmass/delE are stored verbatim.
// ---------------------------------------------------------------------------

namespace net {

using Sp = metal_grain::Sp;
inline constexpr int N_sp = metal_grain::N_sp;
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
constexpr GrainReaction grain_rxn(int num, std::initializer_list<int> reactants,
                                  std::initializer_list<int> products,
                                  int n_reactants) {
  return topo::make_grain(N_sp, num, reactants, products, n_reactants);
}

// Standard / CR / charge-transfer reactions (821 rows).
inline constexpr std::array<Reaction, 821> kReactions = {{
    rxn(1, {Sp::H, Sp::e}, {Sp::e, Sp::e, Sp::Hp}, 3.58687e+40,
        -2.17871e-11),  // H + e- -> e- + e- + H+
    rxn(2, {Sp::e, Sp::Hp}, {Sp::H, PHOTON}, 2.78795e-41,
        2.17871e-11),  // e- + H+ -> H + γ
    rxn(3, {Sp::e, Sp::He}, {Sp::e, Sp::e, Sp::Hep}, 3.58467e+40,
        -3.93933e-11),  // e- + He -> e- + e- + He+
    rxn(4, {Sp::e, Sp::Hep}, {Sp::He, PHOTON}, 2.78965e-41,
        3.93933e-11),  // e- + He+ -> He + γ
    rxn(5, {Sp::e, Sp::Hep}, {Sp::e, Sp::e, Sp::Hepp}, 3.58669e+40,
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
        2.92778e-12),  // H + H2+ -> H2 + H+
    rxn(11, {Sp::H2, Sp::Hp}, {Sp::H, Sp::H2p}, 0.999591,
        -2.92778e-12),  // H2 + H+ -> H + H2+
    rxn(12, {Sp::H2, Sp::e}, {Sp::H, Sp::H, Sp::e}, 1.2923e+36,
        -7.17466e-12),  // H2 + e- -> H + H + e-
    rxn(13, {Sp::H, Sp::H2}, {Sp::H, Sp::H, Sp::H}, 1.2923e+36,
        -7.17466e-12),  // H + H2 -> H + H + H
    rxn(14, {Sp::e, Sp::Hm}, {Sp::H, Sp::e, Sp::e}, 3.58682e+40,
        -1.20867e-12),  // e- + H- -> H + e- + e-
    rxn(15, {Sp::Hp, Sp::Hm}, {Sp::H, Sp::H}, 0.999985,
        2.05784e-11),  // H+ + H- -> H + H
    rxn(16, {Sp::Hp, Sp::Hm}, {Sp::e, Sp::H2p}, 27743.8,
        3.03821e-12),  // H+ + H- -> e- + H2+
    rxn(17, {Sp::e, Sp::H2p}, {Sp::H, Sp::H}, 3.60435e-05,
        1.75402e-11),  // e- + H2+ -> H + H
    rxn(18, {Sp::H2p, Sp::Hm}, {Sp::H, Sp::H2}, 1.00039,
        2.35062e-11),  // H2+ + H- -> H + H2
    rxn(19, {Sp::H, Sp::H, Sp::H}, {Sp::H, Sp::H2}, 7.73811e-37,
        7.17466e-12),  // H + H + H -> H + H2
    rxn(20, {Sp::H, Sp::H, Sp::H2}, {Sp::H2, Sp::H2}, 7.73811e-37,
        7.17466e-12),  // H + H + H2 -> H2 + H2
    rxn(21, {Sp::H2, Sp::H2}, {Sp::H, Sp::H, Sp::H2}, 1.2923e+36,
        -7.17466e-12),  // H2 + H2 -> H + H + H2
    rxn(22, {Sp::H, Sp::H}, {Sp::H, Sp::e, Sp::Hp}, 3.58687e+40,
        -2.17871e-11),                            // H + H -> H + e- + H+
    rxn(23, {Sp::H, Sp::H}, {Sp::H2}, 1.0, 0.0),  // H + H -> H2
    rxn(24, {Sp::H2, Sp::Hep}, {Sp::H, Sp::Hp, Sp::He}, 1.2931e+36,
        1.04316e-11),  // H2 + He+ -> H + H+ + He
    rxn(25, {Sp::H2p, Sp::He}, {Sp::H, Sp::HeHp}, 2.019,
        -1.34603e-12),  // H2+ + He -> H + HeH+
    rxn(26, {Sp::H2, Sp::H2p}, {Sp::H, Sp::H3p}, 1.53938,
        2.69052e-12),  // H2 + H2+ -> H + H3+
    rxn(27, {Sp::H3p, Sp::Hm}, {Sp::H2, Sp::H2}, 0.64987,
        2.08157e-11),  // H3+ + H- -> H2 + H2
    rxn(28, {Sp::H, Sp::Hep}, {Sp::Hp, Sp::He}, 1.00061,
        1.76062e-11),  // H + He+ -> H+ + He
    rxn(29, {Sp::Hm, Sp::Hep}, {Sp::H, Sp::He}, 1.0006,
        3.81847e-11),  // H- + He+ -> H + He
    rxn(30, {Sp::H2, Sp::Hep}, {Sp::H2p, Sp::He}, 1.0002,
        1.46785e-11),  // H2 + He+ -> H2+ + He
    rxn(31, {Sp::H, Sp::HeHp}, {Sp::H2p, Sp::He}, 0.495295,
        1.34603e-12),  // H + HeH+ -> H2+ + He
    rxn(32, {Sp::H2, Sp::HeHp}, {Sp::H3p, Sp::He}, 0.762445,
        4.03655e-12),  // H2 + HeH+ -> H3+ + He
    rxn(33, {Sp::H2p, Sp::Hm}, {Sp::e, Sp::H3p}, 42725.6,
        8.65651e-12),  // H2+ + H- -> e- + H3+
    rxn(34, {Sp::e, Sp::H3p}, {Sp::H, Sp::H2}, 2.34144e-05,
        1.48497e-11),  // e- + H3+ -> H + H2
    rxn(35, {Sp::e, Sp::H3p}, {Sp::H, Sp::H, Sp::H}, 3.02585e+31,
        7.67504e-12),  // e- + H3+ -> H + H + H
    rxn(36, {Sp::e, Sp::HeHp}, {Sp::H, Sp::He}, 1.78522e-05,
        1.88863e-11),  // e- + HeH+ -> H + He
    rxn(37, {Sp::H2, Sp::He}, {Sp::H, Sp::H, Sp::He}, 1.2923e+36,
        -7.17466e-12),  // H2 + He -> H + H + He
    rxn(38, {Sp::H, Sp::Hm}, {Sp::H, Sp::H, Sp::e}, 3.58682e+40,
        -1.20867e-12),  // H + H- -> H + H + e-
    rxn(39, {Sp::H2p, Sp::Hm}, {Sp::H, Sp::H, Sp::H}, 1.29281e+36,
        1.63316e-11),  // H2+ + H- -> H + H + H
    rxn(40, {Sp::H2, Sp::e}, {Sp::H, Sp::Hm}, 3.60293e-05,
        -5.96598e-12),  // H2 + e- -> H + H-
    rxn(41, {Sp::Hp, Sp::He}, {Sp::H, Sp::Hep}, 0.999388,
        -1.76062e-11),  // H+ + He -> H + He+
    rxn(42, {Sp::Hm, Sp::He}, {Sp::H, Sp::e, Sp::He}, 3.58682e+40,
        -1.20867e-12),  // H- + He -> H + e- + He
    rxn(43, {Sp::H, Sp::H, Sp::He}, {Sp::H2, Sp::He}, 7.73811e-37,
        7.17466e-12),  // H + H + He -> H2 + He
    rxn(44, {Sp::H, Sp::H2p}, {Sp::H3p, PHOTON}, 1.19119e-36,
        9.86518e-12),  // H + H2+ -> H3+ + γ
    rxn(45, {Sp::H2, Sp::Hp}, {Sp::H3p, PHOTON}, 1.1907e-36,
        6.9374e-12),  // H2 + H+ -> H3+ + γ
    rxn(46, {Sp::Hp, Sp::He}, {Sp::HeHp, PHOTON}, 1.56169e-36,
        2.90085e-12),  // H+ + He -> HeH+ + γ
    rxn(47, {Sp::e, Sp::HD}, {Sp::H, Sp::e, Sp::D}, 8.39754e+35,
        -7.23181e-12),  // e- + HD -> H + e- + D
    rxn(48, {Sp::He, Sp::HD}, {Sp::H, Sp::He, Sp::D}, 8.39754e+35,
        -7.23181e-12),  // He + HD -> H + He + D
    rxn(49, {Sp::H2, Sp::HD}, {Sp::H, Sp::H2, Sp::D}, 8.39754e+35,
        -7.23181e-12),  // H2 + HD -> H + H2 + D
    rxn(50, {Sp::H, Sp::HD}, {Sp::H, Sp::H, Sp::D}, 8.39754e+35,
        -7.23181e-12),  // H + HD -> H + H + D
    rxn(51, {Sp::e, Sp::Dp}, {Sp::D, PHOTON}, 2.78909e-41,
        2.1793e-11),  // e- + D+ -> D + γ
    rxn(52, {Sp::Hp, Sp::D}, {Sp::H, Sp::Dp}, 0.999591,
        -5.92812e-15),  // H+ + D -> H + D+
    rxn(53, {Sp::H, Sp::Dp}, {Sp::Hp, Sp::D}, 1.00041,
        5.92812e-15),  // H + D+ -> H+ + D
    rxn(54, {Sp::H, Sp::D}, {Sp::HD, PHOTON}, 1.19083e-36,
        7.23181e-12),  // H + D -> HD + γ
    rxn(55, {Sp::H2, Sp::D}, {Sp::H, Sp::HD}, 1.53891,
        5.71558e-14),  // H2 + D -> H + HD
    rxn(56, {Sp::H, Sp::HDp}, {Sp::Hp, Sp::HD}, 1.00055,
        2.95775e-12),  // H + HD+ -> H+ + HD
    rxn(57, {Sp::H2, Sp::Dp}, {Sp::Hp, Sp::HD}, 1.53954,
        6.30839e-14),  // H2 + D+ -> H+ + HD
    rxn(58, {Sp::H, Sp::HD}, {Sp::H2, Sp::D}, 0.649811,
        -5.71558e-14),  // H + HD -> H2 + D
    rxn(59, {Sp::Hp, Sp::HD}, {Sp::H2, Sp::Dp}, 0.649545,
        -6.30839e-14),  // H+ + HD -> H2 + D+
    rxn(60, {Sp::Hp, Sp::D}, {Sp::HDp, PHOTON}, 1.19018e-36,
        4.27406e-12),  // H+ + D -> HD+ + γ
    rxn(61, {Sp::H, Sp::Dp}, {Sp::HDp, PHOTON}, 1.19066e-36,
        4.27999e-12),  // H + D+ -> HD+ + γ
    rxn(62, {Sp::e, Sp::HDp}, {Sp::H, Sp::D}, 2.34247e-05,
        1.7513e-11),  // e- + HD+ -> H + D
    rxn(63, {Sp::e, Sp::D}, {Sp::Dm, PHOTON}, 2.78909e-41,
        1.20909e-12),  // e- + D -> D- + γ
    rxn(64, {Sp::Dp, Sp::Dm}, {Sp::D, Sp::D}, 1.0,
        2.05839e-11),  // D+ + D- -> D + D
    rxn(65, {Sp::Hp, Sp::Dm}, {Sp::H, Sp::D}, 0.999591,
        2.0578e-11),  // H+ + D- -> H + D
    rxn(66, {Sp::Hm, Sp::D}, {Sp::H, Sp::Dm}, 1.00039,
        4.15135e-16),  // H- + D -> H + D-
    rxn(67, {Sp::H, Sp::Dm}, {Sp::Hm, Sp::D}, 0.999606,
        -4.15135e-16),  // H + D- -> H- + D
    rxn(68, {Sp::H, Sp::Dm}, {Sp::e, Sp::HD}, 42695.9,
        6.02273e-12),  // H + D- -> e- + HD
    rxn(69, {Sp::e, Sp::D}, {Sp::e, Sp::e, Sp::Dp}, 3.5854e+40,
        -2.1793e-11),  // e- + D -> e- + e- + D+
    rxn(70, {Sp::Hep, Sp::D}, {Sp::He, Sp::Dp}, 1.0002,
        1.76003e-11),  // He+ + D -> He + D+
    rxn(71, {Sp::He, Sp::Dp}, {Sp::Hep, Sp::D}, 0.999796,
        -1.76003e-11),  // He + D+ -> He+ + D
    rxn(72, {Sp::H2p, Sp::D}, {Sp::H, Sp::HDp}, 1.5387,
        2.7183e-14),  // H2+ + D -> H + HD+
    rxn(73, {Sp::D, Sp::HDp}, {Sp::HD, Sp::Dp}, 1.00014,
        2.95182e-12),  // D + HD+ -> HD + D+
    rxn(74, {Sp::H, Sp::HDp}, {Sp::H2p, Sp::D}, 0.6499,
        -2.7183e-14),  // H + HD+ -> H2+ + D
    rxn(75, {Sp::e, Sp::HD}, {Sp::H, Sp::Dm}, 2.34215e-05,
        -6.02273e-12),  // e- + HD -> H + D-
    rxn(76, {Sp::e, Sp::HD}, {Sp::Hm, Sp::D}, 2.34122e-05,
        -6.02314e-12),  // e- + HD -> H- + D
    rxn(77, {Sp::Hp, Sp::Dm}, {Sp::e, Sp::HDp}, 42672.6,
        3.06497e-12),  // H+ + D- -> e- + HD+
    rxn(78, {Sp::Hm, Sp::Dp}, {Sp::e, Sp::HDp}, 42706.9,
        3.07132e-12),  // H- + D+ -> e- + HD+
    rxn(79, {Sp::e, Sp::Dm}, {Sp::e, Sp::e, Sp::D}, 3.5854e+40,
        -1.20909e-12),  // e- + D- -> e- + e- + D
    rxn(80, {Sp::H, Sp::Dm}, {Sp::H, Sp::e, Sp::D}, 3.5854e+40,
        -1.20909e-12),  // H + D- -> H + e- + D
    rxn(81, {Sp::He, Sp::Dm}, {Sp::e, Sp::He, Sp::D}, 3.5854e+40,
        -1.20909e-12),  // He + D- -> e- + He + D
    rxn(82, {Sp::Hm, Sp::Dp}, {Sp::H, Sp::D}, 1.00039,
        2.05844e-11),  // H- + D+ -> H + D
    rxn(83, {Sp::H2p, Sp::Dm}, {Sp::H2, Sp::D}, 1.0,
        2.35058e-11),  // H2+ + D- -> H2 + D
    rxn(84, {Sp::H2p, Sp::Dm}, {Sp::H, Sp::H, Sp::D}, 1.2923e+36,
        1.63311e-11),  // H2+ + D- -> H + H + D
    rxn(85, {Sp::Hm, Sp::HDp}, {Sp::H, Sp::HD}, 1.00053,
        2.35362e-11),  // H- + HD+ -> H + HD
    rxn(86, {Sp::Hm, Sp::HDp}, {Sp::H, Sp::H, Sp::D}, 8.40199e+35,
        1.63044e-11),  // H- + HD+ -> H + H + D
    rxn(87, {Sp::HDp, Sp::Dm}, {Sp::D, Sp::HD}, 1.00014,
        2.35358e-11),  // HD+ + D- -> D + HD
    rxn(88, {Sp::HDp, Sp::Dm}, {Sp::H, Sp::D, Sp::D}, 8.39868e+35,
        1.6304e-11),  // HD+ + D- -> H + D + D
    rxn(89, {Sp::Hep, Sp::Dm}, {Sp::He, Sp::D}, 1.0002,
        3.81842e-11),  // He+ + D- -> He + D
    rxn(90, {Sp::H2p, Sp::D}, {Sp::H2, Sp::Dp}, 1.0,
        2.92185e-12),  // H2+ + D -> H2 + D+
    rxn(91, {Sp::H2p, Sp::D}, {Sp::Hp, Sp::HD}, 1.53954,
        2.98494e-12),  // H2+ + D -> H+ + HD
    rxn(92, {Sp::H, Sp::HDp}, {Sp::H2, Sp::Dp}, 0.6499,
        2.89467e-12),  // H + HD+ -> H2 + D+
    rxn(93, {Sp::H2, Sp::Dp}, {Sp::H2p, Sp::D}, 1.0,
        -2.92185e-12),  // H2 + D+ -> H2+ + D
    rxn(94, {Sp::H2, Sp::Dp}, {Sp::H, Sp::HDp}, 1.5387,
        -2.89467e-12),  // H2 + D+ -> H + HD+
    rxn(95, {Sp::Hp, Sp::HD}, {Sp::H, Sp::HDp}, 0.999454,
        -2.95775e-12),  // H+ + HD -> H + HD+
    rxn(96, {Sp::Hp, Sp::HD}, {Sp::H2p, Sp::D}, 0.649545,
        -2.98494e-12),  // H+ + HD -> H2+ + D
    rxn(97, {Sp::HD, Sp::Dp}, {Sp::D, Sp::HDp}, 0.999863,
        -2.95182e-12),  // HD + D+ -> D + HD+
    rxn(98, {Sp::Hep, Sp::HD}, {Sp::He, Sp::HDp}, 1.00007,
        1.46485e-11),  // He+ + HD -> He + HD+
    rxn(99, {Sp::Hep, Sp::HD}, {Sp::Hp, Sp::He, Sp::D}, 8.40268e+35,
        1.03744e-11),  // He+ + HD -> H+ + He + D
    rxn(100, {Sp::Hep, Sp::HD}, {Sp::H, Sp::He, Sp::Dp}, 8.39925e+35,
        1.03685e-11),  // He+ + HD -> H + He + D+
    rxn(101, {Sp::H, Sp::CH}, {Sp::H2, Sp::C}, 0.399019,
        1.6163e-12),  // H + CH -> H2 + C
    rxn(102, {Sp::H, Sp::CH}, {Sp::H, Sp::H, Sp::C}, 5.15654e+35,
        -5.55836e-12),  // H + CH -> H + H + C
    rxn(103, {Sp::H, Sp::CH2}, {Sp::H2, Sp::CH}, 0.395431,
        2.36361e-13),  // H + CH2 -> H2 + CH
    rxn(104, {Sp::H, Sp::CH3}, {Sp::H2, Sp::CH2}, 0.392369,
        -4.09589e-13),  // H + CH3 -> H2 + CH2
    rxn(105, {Sp::H, Sp::CH4}, {Sp::H2, Sp::CH3}, 0.389723,
        -9.94663e-15),  // H + CH4 -> H2 + CH3
    rxn(106, {Sp::H, Sp::OH}, {Sp::H2, Sp::O}, 0.387493,
        1.05345e-13),  // H + OH -> H2 + O
    rxn(107, {Sp::H, Sp::OH}, {Sp::H, Sp::H, Sp::O}, 5.00759e+35,
        -7.06931e-12),  // H + OH -> H + H + O
    rxn(108, {Sp::H, Sp::H2O}, {Sp::H2, Sp::OH}, 0.385452,
        -9.96124e-13),  // H + H2O -> H2 + OH
    rxn(109, {Sp::H, Sp::H2O}, {Sp::H, Sp::H, Sp::OH}, 4.98122e+35,
        -8.17078e-12),  // H + H2O -> H + H + OH
    rxn(110, {Sp::H, Sp::C2}, {Sp::C, Sp::CH}, 0.0609984,
        -4.43952e-12),  // H + C2 -> C + CH
    rxn(111, {Sp::H, Sp::CO}, {Sp::C, Sp::OH}, 0.051422,
        -1.07284e-11),  // H + CO -> C + OH
    rxn(112, {Sp::H, Sp::H2CO}, {Sp::H2, Sp::HCO}, 0.372144,
        1.14219e-12),  // H + H2CO -> H2 + HCO
    rxn(113, {Sp::H, Sp::O2}, {Sp::O, Sp::OH}, 0.0408176,
        -1.12678e-12),  // H + O2 -> O + OH
    rxn(114, {Sp::H, Sp::O2}, {Sp::H, Sp::O, Sp::O}, 2.04398e+34,
        -8.19609e-12),  // H + O2 -> H + O + O
    rxn(115, {Sp::H, Sp::O2H}, {Sp::O, Sp::H2O}, 0.0392231,
        3.70921e-12),  // H + O2H -> O + H2O
    rxn(116, {Sp::H, Sp::O2H}, {Sp::H2, Sp::O2}, 0.370395,
        3.83986e-12),  // H + O2H -> H2 + O2
    rxn(117, {Sp::H, Sp::O2H}, {Sp::OH, Sp::OH}, 0.0390165,
        2.60774e-12),  // H + O2H -> OH + OH
    rxn(118, {Sp::H, Sp::H2O2}, {Sp::H2, Sp::O2H}, 0.369877,
        1.17792e-12),  // H + H2O2 -> H2 + O2H
    rxn(119, {Sp::H, Sp::CO2}, {Sp::OH, Sp::CO}, 0.0284258,
        -1.6671e-12),  // H + CO2 -> OH + CO
    rxn(120, {Sp::H2, Sp::C}, {Sp::H, Sp::CH}, 2.50615,
        -1.6163e-12),  // H2 + C -> H + CH
    rxn(121, {Sp::C, Sp::OH}, {Sp::CH, Sp::O}, 0.971114,
        -1.51096e-12),  // C + OH -> CH + O
    rxn(122, {Sp::C, Sp::CO}, {Sp::C2, Sp::O}, 0.818655,
        -7.7998e-12),  // C + CO -> C2 + O
    rxn(123, {Sp::H2, Sp::O}, {Sp::H, Sp::OH}, 2.58069,
        -1.05345e-13),  // H2 + O -> H + OH
    rxn(124, {Sp::CH, Sp::O}, {Sp::C, Sp::OH}, 1.02975,
        1.51096e-12),  // CH + O -> C + OH
    rxn(125, {Sp::CH2, Sp::O}, {Sp::CH, Sp::OH}, 1.02049,
        1.31017e-13),  // CH2 + O -> CH + OH
    rxn(126, {Sp::CH2, Sp::O}, {Sp::H, Sp::HCO}, 21.2399,
        6.30221e-12),  // CH2 + O -> H + HCO
    rxn(127, {Sp::CH4, Sp::O}, {Sp::CH3, Sp::OH}, 1.00575,
        -1.15291e-13),  // CH4 + O -> CH3 + OH
    rxn(128, {Sp::O, Sp::H2O}, {Sp::OH, Sp::OH}, 0.994734,
        -1.10147e-12),  // O + H2O -> OH + OH
    rxn(129, {Sp::O, Sp::H2CO}, {Sp::OH, Sp::HCO}, 0.960388,
        1.03684e-12),  // O + H2CO -> OH + HCO
    rxn(130, {Sp::O, Sp::H2O2}, {Sp::OH, Sp::O2H}, 0.954538,
        1.07258e-12),  // O + H2O2 -> OH + O2H
    rxn(131, {Sp::O, Sp::CO2}, {Sp::O2, Sp::CO}, 0.696409,
        -5.40323e-13),  // O + CO2 -> O2 + CO
    rxn(132, {Sp::Hp, Sp::O}, {Sp::H, Sp::Op}, 0.999232,
        -3.17495e-14),  // H+ + O -> H + O+
    rxn(133, {Sp::H2, Sp::CH}, {Sp::H, Sp::CH2}, 2.52889,
        -2.36361e-13),  // H2 + CH -> H + CH2
    rxn(134, {Sp::H2, Sp::CH}, {Sp::H, Sp::H2, Sp::C}, 5.15654e+35,
        -5.55836e-12),  // H2 + CH -> H + H2 + C
    rxn(135, {Sp::H2, Sp::CH2}, {Sp::H, Sp::CH3}, 2.54862,
        4.09589e-13),  // H2 + CH2 -> H + CH3
    rxn(136, {Sp::H2, Sp::CH3}, {Sp::H, Sp::CH4}, 2.56593,
        9.94663e-15),  // H2 + CH3 -> H + CH4
    rxn(137, {Sp::H2, Sp::OH}, {Sp::H, Sp::H2O}, 2.59435,
        9.96124e-13),  // H2 + OH -> H + H2O
    rxn(138, {Sp::H2, Sp::OH}, {Sp::H, Sp::H2, Sp::O}, 5.00759e+35,
        -7.06931e-12),  // H2 + OH -> H + H2 + O
    rxn(139, {Sp::H2, Sp::H2O}, {Sp::H, Sp::H2, Sp::OH}, 4.98122e+35,
        -8.17078e-12),  // H2 + H2O -> H + H2 + OH
    rxn(140, {Sp::H2, Sp::O2}, {Sp::H, Sp::O2H}, 2.69982,
        -3.83986e-12),  // H2 + O2 -> H + O2H
    rxn(141, {Sp::H2, Sp::O2}, {Sp::OH, Sp::OH}, 0.105338,
        -1.23212e-12),  // H2 + O2 -> OH + OH
    rxn(142, {Sp::H2, Sp::O2}, {Sp::H2, Sp::O, Sp::O}, 2.04398e+34,
        -8.19609e-12),  // H2 + O2 -> H2 + O + O
    rxn(143, {Sp::H2, Sp::O2H}, {Sp::H, Sp::H2O2}, 2.7036,
        -1.17792e-12),                             // H2 + O2H -> H + H2O2
    rxn(144, {Sp::H, Sp::D}, {Sp::HD}, 1.0, 0.0),  // H + D -> HD
    rxn(145, {Sp::H3p, Sp::O2}, {Sp::H2, Sp::O2Hp}, 1.75317,
        -4.33401e-15),  // H3+ + O2 -> H2 + O2H+
    rxn(146, {Sp::H2, Sp::Cp}, {Sp::H, Sp::CHp}, 2.50613,
        -5.33863e-13),  // H2 + C+ -> H + CH+
    rxn(147, {Sp::CH, Sp::CH4}, {Sp::CH2, Sp::CH3}, 0.985564,
        -2.46308e-13),  // CH + CH4 -> CH2 + CH3
    rxn(148, {Sp::CH, Sp::OH}, {Sp::H, Sp::HCO}, 20.8136,
        6.17119e-12),  // CH + OH -> H + HCO
    rxn(149, {Sp::CH, Sp::HCO}, {Sp::CH2, Sp::CO}, 0.942817,
        5.93711e-12),  // CH + HCO -> CH2 + CO
    rxn(150, {Sp::CH, Sp::H2CO}, {Sp::CH2, Sp::HCO}, 0.941109,
        9.05824e-13),  // CH + H2CO -> CH2 + HCO
    rxn(151, {Sp::CH, Sp::O2}, {Sp::O, Sp::HCO}, 0.84956,
        5.04442e-12),  // CH + O2 -> O + HCO
    rxn(152, {Sp::CH, Sp::O2H}, {Sp::CH2, Sp::O2}, 0.936686,
        3.6035e-12),  // CH + O2H -> CH2 + O2
    rxn(153, {Sp::CH, Sp::O2H}, {Sp::OH, Sp::HCO}, 0.812073,
        8.77894e-12),  // CH + O2H -> OH + HCO
    rxn(154, {Sp::CH, Sp::CO2}, {Sp::CO, Sp::HCO}, 0.591641,
        4.5041e-12),  // CH + CO2 -> CO + HCO
    rxn(155, {Sp::CH2, Sp::CH2}, {Sp::CH, Sp::CH3}, 1.00781,
        6.4595e-13),  // CH2 + CH2 -> CH + CH3
    rxn(156, {Sp::CH2, Sp::CH4}, {Sp::CH3, Sp::CH3}, 0.993256,
        3.99642e-13),  // CH2 + CH4 -> CH3 + CH3
    rxn(157, {Sp::CH2, Sp::HCO}, {Sp::CH3, Sp::CO}, 0.950176,
        6.58306e-12),  // CH2 + HCO -> CH3 + CO
    rxn(158, {Sp::CH2, Sp::H2CO}, {Sp::CH3, Sp::HCO}, 0.948455,
        1.55177e-12),  // CH2 + H2CO -> CH3 + HCO
    rxn(159, {Sp::H, Sp::CH2p}, {Sp::H2, Sp::CHp}, 0.395433,
        -1.82693e-13),  // H + CH2+ -> H2 + CH+
    rxn(160, {Sp::CH3, Sp::H2CO}, {Sp::CH4, Sp::HCO}, 0.954894,
        1.15213e-12),  // CH3 + H2CO -> CH4 + HCO
    rxn(161, {Sp::H, Sp::CH3p}, {Sp::H2, Sp::CH2p}, 0.39237,
        -1.2904e-12),  // H + CH3+ -> H2 + CH2+
    rxn(162, {Sp::CH2, Sp::OH}, {Sp::CH3, Sp::O}, 0.987574,
        5.14933e-13),  // CH2 + OH -> CH3 + O
    rxn(163, {Sp::CH2, Sp::OH}, {Sp::CH, Sp::H2O}, 1.02589,
        1.23249e-12),  // CH2 + OH -> CH + H2O
    rxn(164, {Sp::CH2, Sp::OH}, {Sp::H, Sp::H2CO}, 22.116,
        5.26537e-12),  // CH2 + OH -> H + H2CO
    rxn(165, {Sp::CH3, Sp::OH}, {Sp::CH2, Sp::H2O}, 1.01794,
        5.86536e-13),  // CH3 + OH -> CH2 + H2O
    rxn(166, {Sp::CH4, Sp::OH}, {Sp::CH3, Sp::H2O}, 1.01108,
        9.86178e-13),  // CH4 + OH -> CH3 + H2O
    rxn(167, {Sp::OH, Sp::OH}, {Sp::O, Sp::H2O}, 1.00529,
        1.10147e-12),  // OH + OH -> O + H2O
    rxn(168, {Sp::OH, Sp::CO}, {Sp::H, Sp::CO2}, 35.1793,
        1.6671e-12),  // OH + CO -> H + CO2
    rxn(169, {Sp::OH, Sp::H2O2}, {Sp::H2O, Sp::O2H}, 0.959592,
        2.17404e-12),  // OH + H2O2 -> H2O + O2H
    rxn(170, {Sp::CH3, Sp::H2O}, {Sp::CH4, Sp::OH}, 0.989043,
        -9.86178e-13),  // CH3 + H2O -> CH4 + OH
    rxn(171, {Sp::O2, Sp::CO}, {Sp::O, Sp::CO2}, 1.43594,
        5.40323e-13),  // O2 + CO -> O + CO2
    rxn(172, {Sp::CO, Sp::O2H}, {Sp::OH, Sp::CO2}, 1.37258,
        4.27484e-12),  // CO + O2H -> OH + CO2
    rxn(173, {Sp::CH2, Sp::O2}, {Sp::OH, Sp::HCO}, 0.866964,
        5.17544e-12),  // CH2 + O2 -> OH + HCO
    rxn(174, {Sp::CH2, Sp::O2}, {Sp::O, Sp::H2CO}, 0.902722,
        4.13859e-12),  // CH2 + O2 -> O + H2CO
    rxn(175, {Sp::CH3, Sp::O2}, {Sp::CH2, Sp::O2H}, 1.05933,
        -4.24945e-12),  // CH3 + O2 -> CH2 + O2H
    rxn(176, {Sp::CH3, Sp::O2}, {Sp::OH, Sp::H2CO}, 0.91408,
        3.62366e-12),  // CH3 + O2 -> OH + H2CO
    rxn(177, {Sp::CH4, Sp::O2}, {Sp::CH3, Sp::O2H}, 1.05218,
        -3.84981e-12),  // CH4 + O2 -> CH3 + O2H
    rxn(178, {Sp::O2, Sp::HCO}, {Sp::CO, Sp::O2H}, 1.00655,
        2.33361e-12),  // O2 + HCO -> CO + O2H
    rxn(179, {Sp::CH3, Sp::O2H}, {Sp::CH4, Sp::O2}, 0.950406,
        3.84981e-12),  // CH3 + O2H -> CH4 + O2
    rxn(180, {Sp::H2O, Sp::O2H}, {Sp::OH, Sp::H2O2}, 1.04211,
        -2.17404e-12),  // H2O + O2H -> OH + H2O2
    rxn(181, {Sp::HCO, Sp::O2H}, {Sp::O2, Sp::H2CO}, 0.9953,
        2.69768e-12),  // HCO + O2H -> O2 + H2CO
    rxn(182, {Sp::O2H, Sp::H2CO}, {Sp::HCO, Sp::H2O2}, 1.00613,
        -3.57348e-14),  // O2H + H2CO -> HCO + H2O2
    rxn(183, {Sp::O2H, Sp::O2H}, {Sp::O2, Sp::H2O2}, 1.0014,
        2.66194e-12),  // O2H + O2H -> O2 + H2O2
    rxn(184, {Sp::H, Sp::HCO}, {Sp::H2, Sp::CO}, 0.372819,
        6.17347e-12),  // H + HCO -> H2 + CO
    rxn(185, {Sp::H, Sp::C}, {Sp::CH, PHOTON}, 1.93928e-36,
        5.55836e-12),  // H + C -> CH + γ
    rxn(186, {Sp::C, Sp::C}, {Sp::C2, PHOTON}, 3.17924e-35,
        9.99787e-12),  // C + C -> C2 + γ
    rxn(187, {Sp::C, Sp::O}, {Sp::CO, PHOTON}, 3.88349e-35,
        1.77977e-11),  // C + O -> CO + γ
    rxn(188, {Sp::C, Sp::CH}, {Sp::H, Sp::C2}, 16.3939,
        4.43952e-12),  // C + CH -> H + C2
    rxn(189, {Sp::C, Sp::OH}, {Sp::H, Sp::CO}, 19.4469,
        1.07284e-11),  // C + OH -> H + CO
    rxn(190, {Sp::C, Sp::HCO}, {Sp::CH, Sp::CO}, 0.934339,
        4.55717e-12),  // C + HCO -> CH + CO
    rxn(191, {Sp::C, Sp::O2}, {Sp::O, Sp::CO}, 0.793777,
        9.60159e-12),  // C + O2 -> O + CO
    rxn(192, {Sp::CH3, Sp::HCO}, {Sp::CH4, Sp::CO}, 0.956627,
        6.18342e-12),  // CH3 + HCO -> CH4 + CO
    rxn(193, {Sp::H, Sp::O}, {Sp::OH, PHOTON}, 1.99697e-36,
        7.06931e-12),  // H + O -> OH + γ
    rxn(194, {Sp::O, Sp::O}, {Sp::O2, PHOTON}, 4.89242e-35,
        8.19609e-12),  // O + O -> O2 + γ
    rxn(195, {Sp::CH, Sp::O}, {Sp::H, Sp::CO}, 20.0254,
        1.22393e-11),  // CH + O -> H + CO
    rxn(196, {Sp::CH, Sp::O}, {Sp::e, Sp::HCOp}, 1489670.0,
        1.97438e-13),  // CH + O -> e- + HCO+
    rxn(197, {Sp::CH2, Sp::O}, {Sp::H, Sp::H, Sp::CO}, 1.02333e+37,
        5.30102e-12),  // CH2 + O -> H + H + CO
    rxn(198, {Sp::CH3, Sp::O}, {Sp::H, Sp::H2CO}, 22.3942,
        4.75044e-12),  // CH3 + O -> H + H2CO
    rxn(199, {Sp::O, Sp::OH}, {Sp::H, Sp::O2}, 24.4992,
        1.12678e-12),  // O + OH -> H + O2
    rxn(200, {Sp::C2, Sp::O}, {Sp::C, Sp::CO}, 1.22152,
        7.7998e-12),  // C2 + O -> C + CO
    rxn(201, {Sp::O, Sp::HCO}, {Sp::H, Sp::CO2}, 33.8471,
        7.73522e-12),  // O + HCO -> H + CO2
    rxn(202, {Sp::O, Sp::HCO}, {Sp::OH, Sp::CO}, 0.962131,
        6.06812e-12),  // O + HCO -> OH + CO
    rxn(203, {Sp::O, Sp::O2H}, {Sp::O2, Sp::OH}, 0.955875,
        3.73452e-12),  // O + O2H -> O2 + OH
    rxn(204, {Sp::OH, Sp::HCO}, {Sp::CO, Sp::H2O}, 0.967225,
        7.16959e-12),  // OH + HCO -> CO + H2O
    rxn(205, {Sp::OH, Sp::H2CO}, {Sp::H2O, Sp::HCO}, 0.965473,
        2.13831e-12),  // OH + H2CO -> H2O + HCO
    rxn(206, {Sp::OH, Sp::O2H}, {Sp::O2, Sp::H2O}, 0.960935,
        4.83599e-12),  // OH + O2H -> O2 + H2O
    rxn(207, {Sp::HCO, Sp::HCO}, {Sp::CO, Sp::H2CO}, 1.00181,
        5.03128e-12),  // HCO + HCO -> CO + H2CO
    rxn(208, {Sp::Hp, Sp::CH}, {Sp::H, Sp::CHp}, 0.999245,
        4.72257e-12),  // H+ + CH -> H + CH+
    rxn(209, {Sp::Hp, Sp::CH2}, {Sp::H2, Sp::CHp}, 0.395133,
        4.95893e-12),  // H+ + CH2 -> H2 + CH+
    rxn(210, {Sp::Hp, Sp::CH2}, {Sp::H, Sp::CH2p}, 0.99924,
        5.14163e-12),  // H+ + CH2 -> H + CH2+
    rxn(211, {Sp::Hp, Sp::CH3}, {Sp::H, Sp::CH3p}, 0.999236,
        6.02244e-12),  // H+ + CH3 -> H + CH3+
    rxn(212, {Sp::Hp, Sp::CH4}, {Sp::H2, Sp::CH3p}, 0.389425,
        6.0125e-12),  // H+ + CH4 -> H2 + CH3+
    rxn(213, {Sp::Hp, Sp::CH4}, {Sp::H, Sp::CH4p}, 0.999233,
        1.57016e-12),  // H+ + CH4 -> H + CH4+
    rxn(214, {Sp::Hp, Sp::OH}, {Sp::H, Sp::OHp}, 0.99923,
        9.32609e-13),  // H+ + OH -> H + OH+
    rxn(215, {Sp::Hp, Sp::H2O}, {Sp::H, Sp::H2Op}, 0.999226,
        1.57148e-12),  // H+ + H2O -> H + H2O+
    rxn(216, {Sp::Hp, Sp::C2}, {Sp::H, Sp::C2p}, 0.999216,
        2.72428e-12),  // H+ + C2 -> H + C2+
    rxn(217, {Sp::Hp, Sp::HCO}, {Sp::H2, Sp::COp}, 0.372525,
        5.50759e-12),  // H+ + HCO -> H2 + CO+
    rxn(218, {Sp::Hp, Sp::HCO}, {Sp::H2p, Sp::CO}, 0.372666,
        3.24569e-12),  // H+ + HCO -> H2+ + CO
    rxn(219, {Sp::Hp, Sp::HCO}, {Sp::H, Sp::HCOp}, 0.99921,
        8.74403e-12),  // H+ + HCO -> H + HCO+
    rxn(220, {Sp::Hp, Sp::H2CO}, {Sp::H, Sp::H2COp}, 0.999208,
        4.33258e-12),  // H+ + H2CO -> H + H2CO+
    rxn(221, {Sp::Hp, Sp::H2CO}, {Sp::H2, Sp::HCOp}, 0.37185,
        9.88622e-12),  // H+ + H2CO -> H2 + HCO+
    rxn(222, {Sp::Hp, Sp::O2}, {Sp::H, Sp::O2p}, 0.999207,
        2.44182e-12),  // H+ + O2 -> H + O2+
    rxn(223, {Sp::Hp, Sp::CO2}, {Sp::O, Sp::HCOp}, 0.0295213,
        1.00881e-12),  // H+ + CO2 -> O + HCO+
    rxn(224, {Sp::Hm, Sp::C}, {Sp::e, Sp::CH}, 69558.5,
        4.34968e-12),  // H- + C -> e- + CH
    rxn(225, {Sp::Hm, Sp::O}, {Sp::e, Sp::OH}, 71627.6,
        5.86064e-12),  // H- + O -> e- + OH
    rxn(226, {Sp::Hm, Sp::CH}, {Sp::e, Sp::CH2}, 70189.7,
        5.72962e-12),  // H- + CH -> e- + CH2
    rxn(227, {Sp::Hm, Sp::CH2}, {Sp::e, Sp::CH3}, 70737.5,
        6.37557e-12),  // H- + CH2 -> e- + CH3
    rxn(228, {Sp::Hm, Sp::CH3}, {Sp::e, Sp::CH4}, 71217.8,
        5.97593e-12),  // H- + CH3 -> e- + CH4
    rxn(229, {Sp::Hm, Sp::OH}, {Sp::e, Sp::H2O}, 72006.8,
        6.96211e-12),  // H- + OH -> e- + H2O
    rxn(230, {Sp::Hm, Sp::CO}, {Sp::e, Sp::HCO}, 74446.8,
        -2.07484e-13),  // H- + CO -> e- + HCO
    rxn(231, {Sp::Hm, Sp::HCO}, {Sp::e, Sp::H2CO}, 74581.9,
        4.8238e-12),  // H- + HCO -> e- + H2CO
    rxn(232, {Sp::H2p, Sp::C}, {Sp::H, Sp::CHp}, 2.50528,
        6.03405e-12),  // H2+ + C -> H + CH+
    rxn(233, {Sp::H2p, Sp::O}, {Sp::H, Sp::OHp}, 2.57976,
        3.75504e-12),  // H2+ + O -> H + OH+
    rxn(234, {Sp::H2p, Sp::CH}, {Sp::H2, Sp::CHp}, 0.999654,
        7.65035e-12),  // H2+ + CH -> H2 + CH+
    rxn(235, {Sp::H2p, Sp::CH}, {Sp::H, Sp::CH2p}, 2.528,
        7.83304e-12),  // H2+ + CH -> H + CH2+
    rxn(236, {Sp::H2p, Sp::CH2}, {Sp::H, Sp::CH3p}, 2.54772,
        9.35981e-12),  // H2+ + CH2 -> H + CH3+
    rxn(237, {Sp::H2p, Sp::CH2}, {Sp::H2, Sp::CH2p}, 0.99965,
        8.06941e-12),  // H2+ + CH2 -> H2 + CH2+
    rxn(238, {Sp::H2p, Sp::CH4}, {Sp::H2, Sp::CH4p}, 0.999642,
        4.49794e-12),  // H2+ + CH4 -> H2 + CH4+
    rxn(239, {Sp::H2p, Sp::CH4}, {Sp::H, Sp::CH5p}, 2.58028,
        4.70291e-12),  // H2+ + CH4 -> H + CH5+
    rxn(240, {Sp::H2p, Sp::CH4}, {Sp::H, Sp::H2, Sp::CH3p}, 5.03462e+35,
        1.76562e-12),  // H2+ + CH4 -> H + H2 + CH3+
    rxn(241, {Sp::H2p, Sp::OH}, {Sp::H2, Sp::OHp}, 0.999639,
        3.86039e-12),  // H2+ + OH -> H2 + OH+
    rxn(242, {Sp::H2p, Sp::OH}, {Sp::H, Sp::H2Op}, 2.59341,
        5.49539e-12),  // H2+ + OH -> H + H2O+
    rxn(243, {Sp::H2p, Sp::H2O}, {Sp::H2, Sp::H2Op}, 0.999636,
        4.49926e-12),  // H2+ + H2O -> H2 + H2O+
    rxn(244, {Sp::H2p, Sp::H2O}, {Sp::H, Sp::H3Op}, 2.60563,
        7.12693e-12),  // H2+ + H2O -> H + H3O+
    rxn(245, {Sp::H2p, Sp::C2}, {Sp::H2, Sp::C2p}, 0.999625,
        5.65206e-12),  // H2+ + C2 -> H2 + C2+
    rxn(246, {Sp::H2p, Sp::CO}, {Sp::H2, Sp::COp}, 0.99962,
        2.2619e-12),  // H2+ + CO -> H2 + CO+
    rxn(247, {Sp::H2p, Sp::CO}, {Sp::H, Sp::HCOp}, 2.68124,
        5.49834e-12),  // H2+ + CO -> H + HCO+
    rxn(248, {Sp::H2p, Sp::HCO}, {Sp::H3p, Sp::CO}, 0.573909,
        8.86399e-12),  // H2+ + HCO -> H3+ + CO
    rxn(249, {Sp::H2p, Sp::HCO}, {Sp::H2, Sp::HCOp}, 0.999619,
        1.16718e-11),  // H2+ + HCO -> H2 + HCO+
    rxn(250, {Sp::H2p, Sp::H2CO}, {Sp::H2, Sp::H2COp}, 0.999618,
        7.26036e-12),  // H2+ + H2CO -> H2 + H2CO+
    rxn(251, {Sp::H2p, Sp::H2CO}, {Sp::H, Sp::H2, Sp::HCOp}, 4.8074e+35,
        5.63934e-12),  // H2+ + H2CO -> H + H2 + HCO+
    rxn(252, {Sp::H2p, Sp::O2}, {Sp::H, Sp::O2Hp}, 2.69878,
        2.68619e-12),  // H2+ + O2 -> H + O2H+
    rxn(253, {Sp::H2p, Sp::O2}, {Sp::H2, Sp::O2p}, 0.999617,
        5.3696e-12),  // H2+ + O2 -> H2 + O2+
    rxn(254, {Sp::H2p, Sp::CO2}, {Sp::H, Sp::HOCOp}, 2.73285,
        4.62101e-12),  // H2+ + CO2 -> H + HOCO+
    rxn(255, {Sp::H3p, Sp::C}, {Sp::H2, Sp::CHp}, 1.62746,
        3.34353e-12),  // H3+ + C -> H2 + CH+
    rxn(256, {Sp::H3p, Sp::O}, {Sp::H2, Sp::OHp}, 1.67585,
        1.06452e-12),  // H3+ + O -> H2 + OH+
    rxn(257, {Sp::H3p, Sp::CH}, {Sp::H2, Sp::CH2p}, 1.64222,
        5.14252e-12),  // H3+ + CH -> H2 + CH2+
    rxn(258, {Sp::H3p, Sp::CH2}, {Sp::H2, Sp::CH3p}, 1.65504,
        6.66929e-12),  // H3+ + CH2 -> H2 + CH3+
    rxn(259, {Sp::H3p, Sp::CH3}, {Sp::H2, Sp::CH4p}, 1.66627,
        1.81736e-12),  // H3+ + CH3 -> H2 + CH4+
    rxn(260, {Sp::H3p, Sp::CH4}, {Sp::H2, Sp::CH5p}, 1.67619,
        2.01239e-12),  // H3+ + CH4 -> H2 + CH5+
    rxn(261, {Sp::H3p, Sp::OH}, {Sp::H2, Sp::H2Op}, 1.68472,
        2.80487e-12),  // H3+ + OH -> H2 + H2O+
    rxn(262, {Sp::H3p, Sp::H2O}, {Sp::H2, Sp::H3Op}, 1.69266,
        4.43641e-12),  // H3+ + H2O -> H2 + H3O+
    rxn(263, {Sp::H3p, Sp::CO}, {Sp::H2, Sp::HCOp}, 1.74177,
        2.80782e-12),  // H3+ + CO -> H2 + HCO+
    rxn(264, {Sp::H3p, Sp::HCO}, {Sp::H2, Sp::H2COp}, 1.74493,
        3.42765e-12),  // H3+ + HCO -> H2 + H2CO+
    rxn(265, {Sp::H3p, Sp::H2CO}, {Sp::H2, Sp::H2COHp}, 1.74789,
        4.76545e-12),  // H3+ + H2CO -> H2 + H2COH+
    rxn(266, {Sp::H3p, Sp::CO2}, {Sp::H2, Sp::HOCOp}, 1.7753,
        1.93049e-12),  // H3+ + CO2 -> H2 + HOCO+
    rxn(267, {Sp::Hep, Sp::CH}, {Sp::H, Sp::He, Sp::Cp}, 5.15583e+35,
        1.5688e-11),  // He+ + CH -> H + He + C+
    rxn(268, {Sp::Hep, Sp::CH}, {Sp::He, Sp::CHp}, 0.999857,
        2.23288e-11),  // He+ + CH -> He + CH+
    rxn(269, {Sp::Hep, Sp::CH2}, {Sp::H, Sp::He, Sp::CHp}, 5.10945e+35,
        1.53905e-11),  // He+ + CH2 -> H + He + CH+
    rxn(270, {Sp::Hep, Sp::CH2}, {Sp::H2, Sp::He, Sp::Cp}, 2.03878e+35,
        1.59244e-11),  // He+ + CH2 -> H2 + He + C+
    rxn(271, {Sp::Hep, Sp::CH3}, {Sp::H2, Sp::He, Sp::CHp}, 2.00479e+35,
        1.49809e-11),  // He+ + CH3 -> H2 + He + CH+
    rxn(272, {Sp::Hep, Sp::CH4}, {Sp::Hp, Sp::He, Sp::CH3}, 5.03949e+35,
        1.04216e-11),  // He+ + CH4 -> H+ + He + CH3
    rxn(273, {Sp::Hep, Sp::CH4}, {Sp::H2, Sp::He, Sp::CHp}, 1.0,
        0.0),  // He+ + CH4 -> H2 + He + CH+
    rxn(274, {Sp::Hep, Sp::CH4}, {Sp::H2, Sp::He, Sp::CH2p}, 1.97584e+35,
        1.51537e-11),  // He+ + CH4 -> H2 + He + CH2+
    rxn(275, {Sp::Hep, Sp::CH4}, {Sp::H, Sp::He, Sp::CH3p}, 5.03564e+35,
        1.64441e-11),  // He+ + CH4 -> H + He + CH3+
    rxn(276, {Sp::Hep, Sp::CH4}, {Sp::He, Sp::CH4p}, 0.999845,
        1.91764e-11),  // He+ + CH4 -> He + CH4+
    rxn(277, {Sp::Hep, Sp::OH}, {Sp::H, Sp::He, Sp::Op}, 5.00681e+35,
        1.05052e-11),  // He+ + OH -> H + He + O+
    rxn(278, {Sp::Hep, Sp::H2O}, {Sp::Hp, Sp::He, Sp::OH}, 4.98427e+35,
        9.43545e-12),  // He+ + H2O -> H+ + He + OH
    rxn(279, {Sp::Hep, Sp::H2O}, {Sp::H, Sp::He, Sp::OHp}, 4.98043e+35,
        1.03681e-11),  // He+ + H2O -> H + He + OH+
    rxn(280, {Sp::Hep, Sp::H2O}, {Sp::He, Sp::H2Op}, 0.999839,
        1.91777e-11),  // He+ + H2O -> He + H2O+
    rxn(281, {Sp::Hep, Sp::C2}, {Sp::He, Sp::C, Sp::Cp}, 3.14498e+34,
        1.12485e-11),  // He+ + C2 -> He + C + C+
    rxn(282, {Sp::Hep, Sp::C2}, {Sp::He, Sp::C2p}, 0.999828,
        2.03305e-11),  // He+ + C2 -> He + C2+
    rxn(283, {Sp::Hep, Sp::CO}, {Sp::He, Sp::Cp, Sp::O}, 2.57465e+34,
        3.44869e-12),  // He+ + CO -> He + C+ + O
    rxn(284, {Sp::Hep, Sp::HCO}, {Sp::H, Sp::He, Sp::COp}, 4.81711e+35,
        1.59392e-11),  // He+ + HCO -> H + He + CO+
    rxn(285, {Sp::Hep, Sp::HCO}, {Sp::He, Sp::CHp, Sp::O}, 2.40558e+34,
        9.0883e-12),  // He+ + HCO -> He + CH+ + O
    rxn(286, {Sp::Hep, Sp::HCO}, {Sp::HeHp, Sp::CO}, 0.752875,
        1.95059e-11),  // He+ + HCO -> HeH+ + CO
    rxn(287, {Sp::Hep, Sp::H2CO}, {Sp::H2, Sp::He, Sp::COp}, 1.79266e+35,
        1.70814e-11),  // He+ + H2CO -> H2 + He + CO+
    rxn(288, {Sp::Hep, Sp::H2CO}, {Sp::H, Sp::He, Sp::HCOp}, 4.80838e+35,
        2.03178e-11),  // He+ + H2CO -> H + He + HCO+
    rxn(289, {Sp::Hep, Sp::O2}, {Sp::He, Sp::O, Sp::Op}, 2.04366e+34,
        9.37839e-12),  // He+ + O2 -> He + O + O+
    rxn(290, {Sp::Hep, Sp::O2}, {Sp::He, Sp::O2p}, 0.99982,
        2.00481e-11),  // He+ + O2 -> He + O2+
    rxn(291, {Sp::Hep, Sp::CO2}, {Sp::He, Sp::C, Sp::O2p}, 1.79293e+34,
        1.71006e-12),  // He+ + CO2 -> He + C + O2+
    rxn(292, {Sp::Hep, Sp::CO2}, {Sp::He, Sp::CO, Sp::Op}, 1.42323e+34,
        8.83807e-12),  // He+ + CO2 -> He + CO + O+
    rxn(293, {Sp::Hep, Sp::CO2}, {Sp::He, Sp::O, Sp::COp}, 1.42319e+34,
        8.20394e-12),  // He+ + CO2 -> He + O + CO+
    rxn(294, {Sp::Hep, Sp::CO2}, {Sp::He, Sp::Cp, Sp::O2}, 1.79301e+34,
        2.90837e-12),  // He+ + CO2 -> He + C+ + O2
    rxn(295, {Sp::H, Sp::Cp}, {Sp::CHp, PHOTON}, 1.93927e-36,
        6.64079e-12),  // H + C+ -> CH+ + γ
    rxn(296, {Sp::Cp, Sp::O}, {Sp::COp, PHOTON}, 3.88333e-35,
        1.34917e-11),  // C+ + O -> CO+ + γ
    rxn(297, {Sp::Hm, Sp::Cp}, {Sp::H, Sp::C}, 1.00074,
        1.69383e-11),  // H- + C+ -> H + C
    rxn(298, {Sp::H2, Sp::Cp}, {Sp::CH2p, PHOTON}, 4.90418e-36,
        6.82349e-12),  // H2 + C+ -> CH2+ + γ
    rxn(299, {Sp::CH, Sp::Cp}, {Sp::H, Sp::C2p}, 16.3933,
        3.52366e-12),  // CH + C+ -> H + C2+
    rxn(300, {Sp::CH, Sp::Cp}, {Sp::C, Sp::CHp}, 0.999995,
        1.08244e-12),  // CH + C+ -> C + CH+
    rxn(301, {Sp::CH2, Sp::Cp}, {Sp::C, Sp::CH2p}, 0.99999,
        1.50149e-12),  // CH2 + C+ -> C + CH2+
    rxn(302, {Sp::Cp, Sp::OH}, {Sp::H, Sp::COp}, 19.4461,
        6.42235e-12),  // C+ + OH -> H + CO+
    rxn(303, {Sp::Cp, Sp::H2O}, {Sp::H, Sp::HCOp}, 20.1051,
        8.66267e-12),  // C+ + H2O -> H + HCO+
    rxn(304, {Sp::Cp, Sp::HCO}, {Sp::C, Sp::HCOp}, 0.99996,
        5.1039e-12),  // C+ + HCO -> C + HCO+
    rxn(305, {Sp::Cp, Sp::HCO}, {Sp::CHp, Sp::CO}, 0.934334,
        5.63961e-12),  // C+ + HCO -> CH+ + CO
    rxn(306, {Sp::Cp, Sp::H2CO}, {Sp::CH2p, Sp::CO}, 0.879306,
        6.96448e-12),  // C+ + H2CO -> CH2+ + CO
    rxn(307, {Sp::Cp, Sp::H2CO}, {Sp::CH, Sp::HCOp}, 0.932609,
        4.62978e-12),  // C+ + H2CO -> CH + HCO+
    rxn(308, {Sp::Cp, Sp::H2CO}, {Sp::C, Sp::H2COp}, 0.999958,
        6.92445e-13),  // C+ + H2CO -> C + H2CO+
    rxn(309, {Sp::Cp, Sp::O2}, {Sp::O, Sp::COp}, 0.793745,
        5.29558e-12),  // C+ + O2 -> O + CO+
    rxn(310, {Sp::Cp, Sp::O2}, {Sp::CO, Sp::Op}, 0.793763,
        5.9297e-12),  // C+ + O2 -> CO + O+
    rxn(311, {Sp::Cp, Sp::CO2}, {Sp::CO, Sp::COp}, 0.552772,
        4.75525e-12),  // C+ + CO2 -> CO + CO+
    rxn(312, {Sp::H, Sp::CHp}, {Sp::H2, Sp::Cp}, 0.399021,
        5.33863e-13),  // H + CH+ -> H2 + C+
    rxn(313, {Sp::C, Sp::CHp}, {Sp::H, Sp::C2p}, 16.3934,
        2.44122e-12),  // C + CH+ -> H + C2+
    rxn(314, {Sp::CHp, Sp::O}, {Sp::H, Sp::COp}, 20.0247,
        6.85087e-12),  // CH+ + O -> H + CO+
    rxn(315, {Sp::H2, Sp::CHp}, {Sp::H, Sp::CH2p}, 2.52887,
        1.82693e-13),  // H2 + CH+ -> H + CH2+
    rxn(316, {Sp::CH, Sp::CHp}, {Sp::H2, Sp::C2p}, 6.54127,
        4.05753e-12),  // CH + CH+ -> H2 + C2+
    rxn(317, {Sp::CHp, Sp::OH}, {Sp::H2, Sp::COp}, 7.75942,
        6.95621e-12),  // CH+ + OH -> H2 + CO+
    rxn(318, {Sp::CHp, Sp::H2O}, {Sp::H2, Sp::HCOp}, 8.02235,
        9.19653e-12),  // CH+ + H2O -> H2 + HCO+
    rxn(319, {Sp::CHp, Sp::H2O}, {Sp::C, Sp::H3Op}, 1.04006,
        1.09288e-12),  // CH+ + H2O -> C + H3O+
    rxn(320, {Sp::CHp, Sp::H2O}, {Sp::H, Sp::H2COp}, 21.5571,
        3.64289e-12),  // CH+ + H2O -> H + H2CO+
    rxn(321, {Sp::CHp, Sp::CO}, {Sp::C, Sp::HCOp}, 1.07024,
        -5.35706e-13),  // CH+ + CO -> C + HCO+
    rxn(322, {Sp::CHp, Sp::HCO}, {Sp::CH2p, Sp::CO}, 0.942813,
        6.35616e-12),  // CH+ + HCO -> CH2+ + CO
    rxn(323, {Sp::CHp, Sp::HCO}, {Sp::CH, Sp::HCOp}, 0.999965,
        4.02146e-12),  // CH+ + HCO -> CH + HCO+
    rxn(324, {Sp::CHp, Sp::H2CO}, {Sp::CH3p, Sp::CO}, 0.894212,
        8.78875e-12),  // CH+ + H2CO -> CH3+ + CO
    rxn(325, {Sp::CHp, Sp::H2CO}, {Sp::CH2, Sp::HCOp}, 0.941076,
        4.92728e-12),  // CH+ + H2CO -> CH2 + HCO+
    rxn(326, {Sp::CHp, Sp::H2CO}, {Sp::C, Sp::H2COHp}, 1.07399,
        1.42192e-12),  // CH+ + H2CO -> C + H2COH+
    rxn(327, {Sp::CHp, Sp::O2}, {Sp::O, Sp::HCOp}, 0.84953,
        9.06588e-12),  // CH+ + O2 -> O + HCO+
    rxn(328, {Sp::CHp, Sp::O2}, {Sp::HCO, Sp::Op}, 0.849549,
        2.90096e-13),  // CH+ + O2 -> HCO + O+
    rxn(329, {Sp::CHp, Sp::O2}, {Sp::OH, Sp::COp}, 0.81736,
        5.72409e-12),  // CH+ + O2 -> OH + CO+
    rxn(330, {Sp::CHp, Sp::CO2}, {Sp::CO, Sp::HCOp}, 0.591621,
        8.52556e-12),  // CH+ + CO2 -> CO + HCO+
    rxn(331, {Sp::CH2p, Sp::O}, {Sp::H, Sp::HCOp}, 21.2393,
        9.90462e-12),  // CH2+ + O -> H + HCO+
    rxn(332, {Sp::H2, Sp::CH2p}, {Sp::H, Sp::CH3p}, 2.54861,
        1.2904e-12),  // H2 + CH2+ -> H + CH3+
    rxn(333, {Sp::CH2p, Sp::H2O}, {Sp::H, Sp::H2COHp}, 22.9449,
        3.65581e-12),  // CH2+ + H2O -> H + H2COH+
    rxn(334, {Sp::CH2p, Sp::HCO}, {Sp::CH3p, Sp::CO}, 0.950172,
        7.46387e-12),  // CH2+ + HCO -> CH3+ + CO
    rxn(335, {Sp::CH2p, Sp::H2CO}, {Sp::CH3, Sp::HCOp}, 0.948426,
        5.15418e-12),  // CH2+ + H2CO -> CH3 + HCO+
    rxn(336, {Sp::CH2p, Sp::O2}, {Sp::OH, Sp::HCOp}, 0.866937,
        8.77784e-12),  // CH2+ + O2 -> OH + HCO+
    rxn(337, {Sp::CH2p, Sp::CO2}, {Sp::CO, Sp::H2COp}, 0.628644,
        2.78922e-12),  // CH2+ + CO2 -> CO + H2CO+
    rxn(338, {Sp::CH3p, Sp::O}, {Sp::H2, Sp::HCOp}, 8.33366,
        8.61421e-12),  // CH3+ + O -> H2 + HCO+
    rxn(339, {Sp::CH3p, Sp::O}, {Sp::H, Sp::H2COp}, 22.3936,
        3.06057e-12),  // CH3+ + O -> H + H2CO+
    rxn(340, {Sp::H2, Sp::CH3p}, {Sp::CH5p, PHOTON}, 5.12508e-36,
        2.93729e-12),  // H2 + CH3+ -> CH5+ + γ
    rxn(341, {Sp::CH3p, Sp::OH}, {Sp::H2, Sp::H2COp}, 8.67737,
        3.16592e-12),  // CH3+ + OH -> H2 + H2CO+
    rxn(342, {Sp::CH3p, Sp::HCO}, {Sp::CH4p, Sp::CO}, 0.956624,
        1.73113e-12),  // CH3+ + HCO -> CH4+ + CO
    rxn(343, {Sp::CH3p, Sp::HCO}, {Sp::CH3, Sp::HCOp}, 0.999974,
        2.72159e-12),  // CH3+ + HCO -> CH3 + HCO+
    rxn(344, {Sp::CH3p, Sp::H2CO}, {Sp::CH4, Sp::HCOp}, 0.954869,
        3.87372e-12),  // CH3+ + H2CO -> CH4 + HCO+
    rxn(345, {Sp::CH3p, Sp::O2}, {Sp::O, Sp::H2COHp}, 0.953366,
        2.23475e-12),  // CH3+ + O2 -> O + H2COH+
    rxn(346, {Sp::H, Sp::Op}, {Sp::Hp, Sp::O}, 1.00077,
        3.17495e-14),  // H + O+ -> H+ + O
    rxn(347, {Sp::Hm, Sp::Op}, {Sp::H, Sp::O}, 1.00075,
        2.06102e-11),  // H- + O+ -> H + O
    rxn(348, {Sp::H2, Sp::Op}, {Sp::H, Sp::OHp}, 2.58069,
        8.59013e-13),  // H2 + O+ -> H + OH+
    rxn(349, {Sp::CH, Sp::Op}, {Sp::CHp, Sp::O}, 1.00001,
        4.75432e-12),  // CH + O+ -> CH+ + O
    rxn(350, {Sp::CH, Sp::Op}, {Sp::H, Sp::COp}, 20.0249,
        1.16052e-11),  // CH + O+ -> H + CO+
    rxn(351, {Sp::CH2, Sp::Op}, {Sp::CH2p, Sp::O}, 1.00001,
        5.17338e-12),  // CH2 + O+ -> CH2+ + O
    rxn(352, {Sp::CH4, Sp::Op}, {Sp::CH3p, Sp::OH}, 1.00576,
        5.9389e-12),  // CH4 + O+ -> CH3+ + OH
    rxn(353, {Sp::CH4, Sp::Op}, {Sp::CH4p, Sp::O}, 1.0,
        1.60191e-12),  // CH4 + O+ -> CH4+ + O
    rxn(354, {Sp::OH, Sp::Op}, {Sp::O, Sp::OHp}, 0.999998,
        9.64358e-13),  // OH + O+ -> O + OH+
    rxn(355, {Sp::OH, Sp::Op}, {Sp::H, Sp::O2p}, 24.4986,
        3.60035e-12),  // OH + O+ -> H + O2+
    rxn(356, {Sp::H2O, Sp::Op}, {Sp::O, Sp::H2Op}, 0.999994,
        1.60323e-12),  // H2O + O+ -> O + H2O+
    rxn(357, {Sp::C2, Sp::Op}, {Sp::C2p, Sp::O}, 0.999984,
        2.75603e-12),  // C2 + O+ -> C2+ + O
    rxn(358, {Sp::C2, Sp::Op}, {Sp::C, Sp::COp}, 1.22149,
        7.16567e-12),  // C2 + O+ -> C + CO+
    rxn(359, {Sp::HCO, Sp::Op}, {Sp::CO, Sp::OHp}, 0.962129,
        7.03248e-12),  // HCO + O+ -> CO + OH+
    rxn(360, {Sp::HCO, Sp::Op}, {Sp::O, Sp::HCOp}, 0.999978,
        8.77578e-12),  // HCO + O+ -> O + HCO+
    rxn(361, {Sp::H2CO, Sp::Op}, {Sp::OH, Sp::HCOp}, 0.960367,
        9.81262e-12),  // H2CO + O+ -> OH + HCO+
    rxn(362, {Sp::H2CO, Sp::Op}, {Sp::O, Sp::H2COp}, 0.999976,
        4.36433e-12),  // H2CO + O+ -> O + H2CO+
    rxn(363, {Sp::O2, Sp::Op}, {Sp::O, Sp::O2p}, 0.999975,
        2.47357e-12),  // O2 + O+ -> O + O2+
    rxn(364, {Sp::CO2, Sp::Op}, {Sp::CO, Sp::O2p}, 0.696392,
        1.93325e-12),  // CO2 + O+ -> CO + O2+
    rxn(365, {Sp::H, Sp::CH4p}, {Sp::H2, Sp::CH3p}, 0.389724,
        4.44234e-12),  // H + CH4+ -> H2 + CH3+
    rxn(366, {Sp::CH4p, Sp::O}, {Sp::CH3p, Sp::OH}, 1.00576,
        4.337e-12),  // CH4+ + O -> CH3+ + OH
    rxn(367, {Sp::H2, Sp::CH4p}, {Sp::H, Sp::CH5p}, 2.5812,
        2.04977e-13),  // H2 + CH4+ -> H + CH5+
    rxn(368, {Sp::CH4, Sp::CH4p}, {Sp::CH3, Sp::CH5p}, 1.00595,
        1.9503e-13),  // CH4 + CH4+ -> CH3 + CH5+
    rxn(369, {Sp::CH4p, Sp::H2O}, {Sp::CH3, Sp::H3Op}, 1.01584,
        2.61905e-12),  // CH4+ + H2O -> CH3 + H3O+
    rxn(370, {Sp::CH4p, Sp::CO}, {Sp::CH3, Sp::HCOp}, 1.04532,
        9.90462e-13),  // CH4+ + CO -> CH3 + HCO+
    rxn(371, {Sp::CH4p, Sp::H2CO}, {Sp::CH4, Sp::H2COp}, 0.999976,
        2.76242e-12),  // CH4+ + H2CO -> CH4 + H2CO+
    rxn(372, {Sp::CH4p, Sp::H2CO}, {Sp::CH3, Sp::H2COHp}, 1.04898,
        2.94809e-12),  // CH4+ + H2CO -> CH3 + H2COH+
    rxn(373, {Sp::CH4p, Sp::O2}, {Sp::CH4, Sp::O2p}, 0.999974,
        8.71667e-13),  // CH4+ + O2 -> CH4 + O2+
    rxn(374, {Sp::CH4p, Sp::CO2}, {Sp::CH3, Sp::HOCOp}, 1.06543,
        1.13133e-13),  // CH4+ + CO2 -> CH3 + HOCO+
    rxn(375, {Sp::C, Sp::OHp}, {Sp::CHp, Sp::O}, 0.971129,
        2.27901e-12),  // C + OH+ -> CH+ + O
    rxn(376, {Sp::O, Sp::OHp}, {Sp::H, Sp::O2p}, 24.4987,
        2.63599e-12),  // O + OH+ -> H + O2+
    rxn(377, {Sp::H2, Sp::OHp}, {Sp::H, Sp::H2Op}, 2.59435,
        1.635e-12),  // H2 + OH+ -> H + H2O+
    rxn(378, {Sp::CH, Sp::OHp}, {Sp::CHp, Sp::OH}, 1.00001,
        3.78996e-12),  // CH + OH+ -> CH+ + OH
    rxn(379, {Sp::CH, Sp::OHp}, {Sp::CH2p, Sp::O}, 0.979936,
        4.078e-12),  // CH + OH+ -> CH2+ + O
    rxn(380, {Sp::CH2, Sp::OHp}, {Sp::CH2p, Sp::OH}, 1.00001,
        4.20902e-12),  // CH2 + OH+ -> CH2+ + OH
    rxn(381, {Sp::CH2, Sp::OHp}, {Sp::CH3p, Sp::O}, 0.98758,
        5.60477e-12),  // CH2 + OH+ -> CH3+ + O
    rxn(382, {Sp::CH4, Sp::OHp}, {Sp::CH2, Sp::H3Op}, 1.03407,
        3.84313e-12),  // CH4 + OH+ -> CH2 + H3O+
    rxn(383, {Sp::CH4, Sp::OHp}, {Sp::CH5p, Sp::O}, 1.0002,
        9.47869e-13),  // CH4 + OH+ -> CH5+ + O
    rxn(384, {Sp::OH, Sp::OHp}, {Sp::O, Sp::H2Op}, 1.00529,
        1.74034e-12),  // OH + OH+ -> O + H2O+
    rxn(385, {Sp::H2O, Sp::OHp}, {Sp::OH, Sp::H2Op}, 0.999996,
        6.38876e-13),  // H2O + OH+ -> OH + H2O+
    rxn(386, {Sp::H2O, Sp::OHp}, {Sp::O, Sp::H3Op}, 1.01003,
        3.37189e-12),  // H2O + OH+ -> O + H3O+
    rxn(387, {Sp::C2, Sp::OHp}, {Sp::C2p, Sp::OH}, 0.999986,
        1.79167e-12),  // C2 + OH+ -> C2+ + OH
    rxn(388, {Sp::CO, Sp::OHp}, {Sp::O, Sp::HCOp}, 1.03934,
        1.7433e-12),  // CO + OH+ -> O + HCO+
    rxn(389, {Sp::HCO, Sp::OHp}, {Sp::CO, Sp::H2Op}, 0.967221,
        7.80847e-12),  // HCO + OH+ -> CO + H2O+
    rxn(390, {Sp::HCO, Sp::OHp}, {Sp::OH, Sp::HCOp}, 0.99998,
        7.81142e-12),  // HCO + OH+ -> OH + HCO+
    rxn(391, {Sp::HCO, Sp::OHp}, {Sp::O, Sp::H2COp}, 1.04122,
        2.36313e-12),  // HCO + OH+ -> O + H2CO+
    rxn(392, {Sp::H2CO, Sp::OHp}, {Sp::O, Sp::H2COHp}, 1.04299,
        3.70093e-12),  // H2CO + OH+ -> O + H2COH+
    rxn(393, {Sp::H2CO, Sp::OHp}, {Sp::OH, Sp::H2COp}, 0.999978,
        3.39997e-12),  // H2CO + OH+ -> OH + H2CO+
    rxn(394, {Sp::O2, Sp::OHp}, {Sp::OH, Sp::O2p}, 0.999977,
        1.50921e-12),  // O2 + OH+ -> OH + O2+
    rxn(395, {Sp::CO2, Sp::OHp}, {Sp::O, Sp::HOCOp}, 1.05934,
        8.65971e-13),  // CO2 + OH+ -> O + HOCO+
    rxn(396, {Sp::H, Sp::CH5p}, {Sp::H2, Sp::CH4p}, 0.387416,
        -2.04977e-13),  // H + CH5+ -> H2 + CH4+
    rxn(397, {Sp::C, Sp::CH5p}, {Sp::CH4, Sp::CHp}, 0.970933,
        1.33114e-12),  // C + CH5+ -> CH4 + CH+
    rxn(398, {Sp::CH5p, Sp::O}, {Sp::CH2, Sp::H3Op}, 1.03386,
        2.89527e-12),  // CH5+ + O -> CH2 + H3O+
    rxn(399, {Sp::CH5p, Sp::O}, {Sp::H2, Sp::H2COHp}, 9.10087,
        7.49355e-12),  // CH5+ + O -> H2 + H2COH+
    rxn(400, {Sp::CH, Sp::CH5p}, {Sp::CH4, Sp::CH2p}, 0.979738,
        3.13013e-12),  // CH + CH5+ -> CH4 + CH2+
    rxn(401, {Sp::CH2, Sp::CH5p}, {Sp::CH4, Sp::CH3p}, 0.987381,
        4.6569e-12),  // CH2 + CH5+ -> CH4 + CH3+
    rxn(402, {Sp::CH5p, Sp::OH}, {Sp::CH4, Sp::H2Op}, 1.00509,
        7.92476e-13),  // CH5+ + OH -> CH4 + H2O+
    rxn(403, {Sp::CH5p, Sp::H2O}, {Sp::CH4, Sp::H3Op}, 1.00982,
        2.42402e-12),  // CH5+ + H2O -> CH4 + H3O+
    rxn(404, {Sp::CH5p, Sp::CO}, {Sp::CH4, Sp::HCOp}, 1.03913,
        7.95431e-13),  // CH5+ + CO -> CH4 + HCO+
    rxn(405, {Sp::CH5p, Sp::HCO}, {Sp::CH4, Sp::H2COp}, 1.04101,
        1.41526e-12),  // CH5+ + HCO -> CH4 + H2CO+
    rxn(406, {Sp::CH5p, Sp::H2CO}, {Sp::CH4, Sp::H2COHp}, 1.04278,
        2.75306e-12),  // CH5+ + H2CO -> CH4 + H2COH+
    rxn(407, {Sp::CH5p, Sp::CO2}, {Sp::CH4, Sp::HOCOp}, 1.05913,
        -8.18978e-14),  // CH5+ + CO2 -> CH4 + HOCO+
    rxn(408, {Sp::C, Sp::H2Op}, {Sp::CHp, Sp::OH}, 0.966018,
        5.38662e-13),  // C + H2O+ -> CH+ + OH
    rxn(409, {Sp::O, Sp::H2Op}, {Sp::H2, Sp::O2p}, 9.4431,
        1.00099e-12),  // O + H2O+ -> H2 + O2+
    rxn(410, {Sp::H2, Sp::H2Op}, {Sp::H, Sp::H3Op}, 2.60658,
        2.62767e-12),  // H2 + H2O+ -> H + H3O+
    rxn(411, {Sp::CH, Sp::H2Op}, {Sp::CHp, Sp::H2O}, 1.00002,
        3.15109e-12),  // CH + H2O+ -> CH+ + H2O
    rxn(412, {Sp::CH, Sp::H2Op}, {Sp::CH2p, Sp::OH}, 0.974778,
        2.33766e-12),  // CH + H2O+ -> CH2+ + OH
    rxn(413, {Sp::CH2, Sp::H2Op}, {Sp::CH3p, Sp::OH}, 0.982383,
        3.86442e-12),  // CH2 + H2O+ -> CH3+ + OH
    rxn(414, {Sp::CH2, Sp::H2Op}, {Sp::CH2p, Sp::H2O}, 1.00001,
        3.57014e-12),  // CH2 + H2O+ -> CH2+ + H2O
    rxn(415, {Sp::CH4, Sp::H2Op}, {Sp::CH3, Sp::H3Op}, 1.01584,
        2.61772e-12),  // CH4 + H2O+ -> CH3 + H3O+
    rxn(416, {Sp::OH, Sp::H2Op}, {Sp::O, Sp::H3Op}, 1.01003,
        2.73301e-12),  // OH + H2O+ -> O + H3O+
    rxn(417, {Sp::H2O, Sp::H2Op}, {Sp::OH, Sp::H3Op}, 1.00471,
        1.63155e-12),  // H2O + H2O+ -> OH + H3O+
    rxn(418, {Sp::C2, Sp::H2Op}, {Sp::C2p, Sp::H2O}, 0.999989,
        1.1528e-12),  // C2 + H2O+ -> C2+ + H2O
    rxn(419, {Sp::CO, Sp::H2Op}, {Sp::OH, Sp::HCOp}, 1.03387,
        2.95576e-15),  // CO + H2O+ -> OH + HCO+
    rxn(420, {Sp::HCO, Sp::H2Op}, {Sp::CO, Sp::H3Op}, 0.971783,
        8.80114e-12),  // HCO + H2O+ -> CO + H3O+
    rxn(421, {Sp::HCO, Sp::H2Op}, {Sp::H2O, Sp::HCOp}, 0.999983,
        7.17255e-12),  // HCO + H2O+ -> H2O + HCO+
    rxn(422, {Sp::HCO, Sp::H2Op}, {Sp::OH, Sp::H2COp}, 1.03574,
        6.22785e-13),  // HCO + H2O+ -> OH + H2CO+
    rxn(423, {Sp::H2CO, Sp::H2Op}, {Sp::H2O, Sp::H2COp}, 0.999982,
        2.76109e-12),  // H2CO + H2O+ -> H2O + H2CO+
    rxn(424, {Sp::H2CO, Sp::H2Op}, {Sp::OH, Sp::H2COHp}, 1.0375,
        1.96058e-12),  // H2CO + H2O+ -> OH + H2COH+
    rxn(425, {Sp::O2, Sp::H2Op}, {Sp::H2O, Sp::O2p}, 0.999981,
        8.70338e-13),  // O2 + H2O+ -> H2O + O2+
    rxn(426, {Sp::C, Sp::H3Op}, {Sp::H2, Sp::HCOp}, 7.71338,
        8.10365e-12),  // C + H3O+ -> H2 + HCO+
    rxn(427, {Sp::Hm, Sp::H3Op}, {Sp::H, Sp::H2, Sp::OH}, 1.91247e+35,
        8.20849e-12),  // H- + H3O+ -> H + H2 + OH
    rxn(428, {Sp::Hm, Sp::H3Op}, {Sp::H2, Sp::H2O}, 0.383935,
        1.63793e-11),  // H- + H3O+ -> H2 + H2O
    rxn(429, {Sp::CH, Sp::H3Op}, {Sp::CH2p, Sp::H2O}, 0.970206,
        7.06111e-13),  // CH + H3O+ -> CH2+ + H2O
    rxn(430, {Sp::CH2, Sp::H3Op}, {Sp::CH3p, Sp::H2O}, 0.977775,
        2.23288e-12),  // CH2 + H3O+ -> CH3+ + H2O
    rxn(431, {Sp::H2CO, Sp::H3Op}, {Sp::H2O, Sp::H2COHp}, 1.03263,
        3.29036e-13),  // H2CO + H3O+ -> H2O + H2COH+
    rxn(432, {Sp::C, Sp::C2p}, {Sp::C2, Sp::Cp}, 1.00003,
        9.15854e-13),  // C + C2+ -> C2 + C+
    rxn(433, {Sp::C2p, Sp::O}, {Sp::C, Sp::COp}, 1.22151,
        4.40964e-12),  // C2+ + O -> C + CO+
    rxn(434, {Sp::CH, Sp::C2p}, {Sp::C2, Sp::CHp}, 1.00003,
        1.99829e-12),  // CH + C2+ -> C2 + CH+
    rxn(435, {Sp::CH2, Sp::C2p}, {Sp::C2, Sp::CH2p}, 1.00002,
        2.41735e-12),  // CH2 + C2+ -> C2 + CH2+
    rxn(436, {Sp::C2p, Sp::OH}, {Sp::C2, Sp::OHp}, 1.00001,
        -1.79167e-12),  // C2+ + OH -> C2 + OH+
    rxn(437, {Sp::C2p, Sp::HCO}, {Sp::C2, Sp::HCOp}, 0.999994,
        6.01975e-12),  // C2+ + HCO -> C2 + HCO+
    rxn(438, {Sp::C2p, Sp::O2}, {Sp::CO, Sp::COp}, 0.969606,
        1.40112e-11),  // C2+ + O2 -> CO + CO+
    rxn(439, {Sp::H, Sp::COp}, {Sp::Hp, Sp::CO}, 1.00079,
        6.65876e-13),  // H + CO+ -> H+ + CO
    rxn(440, {Sp::C, Sp::COp}, {Sp::Cp, Sp::CO}, 1.00004,
        4.30601e-12),  // C + CO+ -> C+ + CO
    rxn(441, {Sp::O, Sp::COp}, {Sp::CO, Sp::Op}, 1.00002,
        6.34127e-13),  // O + CO+ -> CO + O+
    rxn(442, {Sp::H2, Sp::COp}, {Sp::H, Sp::HCOp}, 2.68226,
        3.23644e-12),  // H2 + CO+ -> H + HCO+
    rxn(443, {Sp::CH, Sp::COp}, {Sp::CHp, Sp::CO}, 1.00003,
        5.38845e-12),  // CH + CO+ -> CH+ + CO
    rxn(444, {Sp::CH, Sp::COp}, {Sp::C, Sp::HCOp}, 1.07027,
        4.85274e-12),  // CH + CO+ -> C + HCO+
    rxn(445, {Sp::CH2, Sp::COp}, {Sp::CH2p, Sp::CO}, 1.00003,
        5.8075e-12),  // CH2 + CO+ -> CH2+ + CO
    rxn(446, {Sp::CH2, Sp::COp}, {Sp::CH, Sp::HCOp}, 1.06065,
        3.4728e-12),  // CH2 + CO+ -> CH + HCO+
    rxn(447, {Sp::CH4, Sp::COp}, {Sp::CH4p, Sp::CO}, 1.00002,
        2.23603e-12),  // CH4 + CO+ -> CH4+ + CO
    rxn(448, {Sp::CH4, Sp::COp}, {Sp::CH3, Sp::HCOp}, 1.04534,
        3.22649e-12),  // CH4 + CO+ -> CH3 + HCO+
    rxn(449, {Sp::OH, Sp::COp}, {Sp::CO, Sp::OHp}, 1.00002,
        1.59848e-12),  // OH + CO+ -> CO + OH+
    rxn(450, {Sp::OH, Sp::COp}, {Sp::O, Sp::HCOp}, 1.03936,
        3.34178e-12),  // OH + CO+ -> O + HCO+
    rxn(451, {Sp::H2O, Sp::COp}, {Sp::CO, Sp::H2Op}, 1.00002,
        2.23736e-12),  // H2O + CO+ -> CO + H2O+
    rxn(452, {Sp::H2O, Sp::COp}, {Sp::OH, Sp::HCOp}, 1.03389,
        2.24032e-12),  // H2O + CO+ -> OH + HCO+
    rxn(453, {Sp::C2, Sp::COp}, {Sp::C2p, Sp::CO}, 1.00001,
        3.39016e-12),  // C2 + CO+ -> C2+ + CO
    rxn(454, {Sp::HCO, Sp::COp}, {Sp::CO, Sp::HCOp}, 1.0,
        9.40991e-12),  // HCO + CO+ -> CO + HCO+
    rxn(455, {Sp::H2CO, Sp::COp}, {Sp::HCO, Sp::HCOp}, 0.998188,
        4.37863e-12),  // H2CO + CO+ -> HCO + HCO+
    rxn(456, {Sp::H2CO, Sp::COp}, {Sp::CO, Sp::H2COp}, 0.999998,
        4.99845e-12),  // H2CO + CO+ -> CO + H2CO+
    rxn(457, {Sp::O2, Sp::COp}, {Sp::CO, Sp::O2p}, 0.999997,
        3.1077e-12),  // O2 + CO+ -> CO + O2+
    rxn(458, {Sp::C, Sp::HCOp}, {Sp::CHp, Sp::CO}, 0.934372,
        5.35706e-13),  // C + HCO+ -> CH+ + CO
    rxn(459, {Sp::Hm, Sp::HCOp}, {Sp::H2, Sp::CO}, 0.373108,
        1.80079e-11),  // H- + HCO+ -> H2 + CO
    rxn(460, {Sp::CH, Sp::HCOp}, {Sp::CH2p, Sp::CO}, 0.942846,
        2.3347e-12),  // CH + HCO+ -> CH2+ + CO
    rxn(461, {Sp::CH2, Sp::HCOp}, {Sp::CH3p, Sp::CO}, 0.950201,
        3.86147e-12),  // CH2 + HCO+ -> CH3+ + CO
    rxn(462, {Sp::OH, Sp::HCOp}, {Sp::CO, Sp::H2Op}, 0.967241,
        -2.95576e-15),  // OH + HCO+ -> CO + H2O+
    rxn(463, {Sp::OH, Sp::HCOp}, {Sp::H, Sp::HOCOp}, 35.8564,
        7.89769e-13),  // OH + HCO+ -> H + HOCO+
    rxn(464, {Sp::H2O, Sp::HCOp}, {Sp::CO, Sp::H3Op}, 0.971799,
        1.62859e-12),  // H2O + HCO+ -> CO + H3O+
    rxn(465, {Sp::HCO, Sp::HCOp}, {Sp::CO, Sp::H2COp}, 1.00181,
        6.19829e-13),  // HCO + HCO+ -> CO + H2CO+
    rxn(466, {Sp::H2CO, Sp::HCOp}, {Sp::CO, Sp::H2COHp}, 1.00351,
        1.95763e-12),  // H2CO + HCO+ -> CO + H2COH+
    rxn(467, {Sp::CH, Sp::H2COp}, {Sp::CHp, Sp::H2CO}, 1.00004,
        3.89994e-13),  // CH + H2CO+ -> CH+ + H2CO
    rxn(468, {Sp::CH, Sp::H2COp}, {Sp::CH2p, Sp::HCO}, 0.941139,
        1.71487e-12),  // CH + H2CO+ -> CH2+ + HCO
    rxn(469, {Sp::CH2, Sp::H2COp}, {Sp::CH3p, Sp::HCO}, 0.948481,
        3.24164e-12),  // CH2 + H2CO+ -> CH3+ + HCO
    rxn(470, {Sp::CH2, Sp::H2COp}, {Sp::CH2p, Sp::H2CO}, 1.00003,
        8.09048e-13),  // CH2 + H2CO+ -> CH2+ + H2CO
    rxn(471, {Sp::CH4, Sp::H2COp}, {Sp::CH3, Sp::H2COHp}, 1.04901,
        1.85665e-13),  // CH4 + H2CO+ -> CH3 + H2COH+
    rxn(472, {Sp::H2O, Sp::H2COp}, {Sp::HCO, Sp::H3Op}, 0.970041,
        1.00876e-12),  // H2O + H2CO+ -> HCO + H3O+
    rxn(473, {Sp::HCO, Sp::H2COp}, {Sp::H2CO, Sp::HCOp}, 1.0,
        4.41145e-12),  // HCO + H2CO+ -> H2CO + HCO+
    rxn(474, {Sp::HCO, Sp::H2COp}, {Sp::CO, Sp::H2COHp}, 1.00351,
        6.36908e-12),  // HCO + H2CO+ -> CO + H2COH+
    rxn(475, {Sp::H2CO, Sp::H2COp}, {Sp::HCO, Sp::H2COHp}, 1.00169,
        1.3378e-12),  // H2CO + H2CO+ -> HCO + H2COH+
    rxn(476, {Sp::O2, Sp::H2COp}, {Sp::O2H, Sp::HCOp}, 1.00472,
        1.71378e-12),  // O2 + H2CO+ -> O2H + HCO+
    rxn(477, {Sp::CH, Sp::H2COHp}, {Sp::CH2p, Sp::H2CO}, 0.939548,
        3.77075e-13),  // CH + H2COH+ -> CH2+ + H2CO
    rxn(478, {Sp::H2O, Sp::H2COHp}, {Sp::H2CO, Sp::H3Op}, 0.968401,
        -3.29036e-13),  // H2O + H2COH+ -> H2CO + H3O+
    rxn(479, {Sp::C, Sp::O2p}, {Sp::O, Sp::COp}, 0.793779,
        6.49389e-12),  // C + O2+ -> O + CO+
    rxn(480, {Sp::C, Sp::O2p}, {Sp::Cp, Sp::O2}, 1.00004,
        1.19831e-12),  // C + O2+ -> C+ + O2
    rxn(481, {Sp::CH, Sp::O2p}, {Sp::CHp, Sp::O2}, 1.00004,
        2.28075e-12),  // CH + O2+ -> CH+ + O2
    rxn(482, {Sp::CH, Sp::O2p}, {Sp::O, Sp::HCOp}, 0.849562,
        1.13466e-11),  // CH + O2+ -> O + HCO+
    rxn(483, {Sp::CH2, Sp::O2p}, {Sp::O, Sp::H2COp}, 0.902723,
        6.02935e-12),  // CH2 + O2+ -> O + H2CO+
    rxn(484, {Sp::CH2, Sp::O2p}, {Sp::CH2p, Sp::O2}, 1.00003,
        2.6998e-12),  // CH2 + O2+ -> CH2+ + O2
    rxn(485, {Sp::C2, Sp::O2p}, {Sp::CO, Sp::COp}, 0.969614,
        1.42937e-11),  // C2 + O2+ -> CO + CO+
    rxn(486, {Sp::C2, Sp::O2p}, {Sp::C2p, Sp::O2}, 1.00001,
        2.82458e-13),  // C2 + O2+ -> C2+ + O2
    rxn(487, {Sp::HCO, Sp::O2p}, {Sp::CO, Sp::O2Hp}, 1.00654,
        3.49005e-12),  // HCO + O2+ -> CO + O2H+
    rxn(488, {Sp::HCO, Sp::O2p}, {Sp::O2, Sp::HCOp}, 1.0,
        6.30221e-12),  // HCO + O2+ -> O2 + HCO+
    rxn(489, {Sp::H2CO, Sp::O2p}, {Sp::H, Sp::O2, Sp::HCOp}, 4.80925e+35,
        2.69738e-13),  // H2CO + O2+ -> H + O2 + HCO+
    rxn(490, {Sp::H2CO, Sp::O2p}, {Sp::O2, Sp::H2COp}, 1.0,
        1.89076e-12),  // H2CO + O2+ -> O2 + H2CO+
    rxn(491, {Sp::C, Sp::O2Hp}, {Sp::CHp, Sp::O2}, 0.928299,
        3.34786e-12),  // C + O2H+ -> CH+ + O2
    rxn(492, {Sp::O, Sp::O2Hp}, {Sp::O2, Sp::OHp}, 0.955898,
        1.06886e-12),  // O + O2H+ -> O2 + OH+
    rxn(493, {Sp::H2, Sp::O2Hp}, {Sp::H3p, Sp::O2}, 0.570396,
        4.33401e-15),  // H2 + O2H+ -> H3+ + O2
    rxn(494, {Sp::CH, Sp::O2Hp}, {Sp::CH2p, Sp::O2}, 0.936718,
        5.14686e-12),  // CH + O2H+ -> CH2+ + O2
    rxn(495, {Sp::CH2, Sp::O2Hp}, {Sp::CH3p, Sp::O2}, 0.944026,
        6.67362e-12),  // CH2 + O2H+ -> CH3+ + O2
    rxn(496, {Sp::OH, Sp::O2Hp}, {Sp::O2, Sp::H2Op}, 0.960955,
        2.8092e-12),  // OH + O2H+ -> O2 + H2O+
    rxn(497, {Sp::H2O, Sp::O2Hp}, {Sp::O2, Sp::H3Op}, 0.965484,
        4.44075e-12),  // H2O + O2H+ -> O2 + H3O+
    rxn(498, {Sp::CO, Sp::O2Hp}, {Sp::O2, Sp::HCOp}, 0.993501,
        2.81216e-12),  // CO + O2H+ -> O2 + HCO+
    rxn(499, {Sp::HCO, Sp::O2Hp}, {Sp::O2, Sp::H2COp}, 0.995302,
        3.43199e-12),  // HCO + O2H+ -> O2 + H2CO+
    rxn(500, {Sp::H2CO, Sp::O2Hp}, {Sp::O2, Sp::H2COHp}, 0.996988,
        4.76978e-12),  // H2CO + O2H+ -> O2 + H2COH+
    rxn(501, {Sp::CO2, Sp::O2Hp}, {Sp::O2, Sp::HOCOp}, 1.01262,
        1.93483e-12),  // CO2 + O2H+ -> O2 + HOCO+
    rxn(502, {Sp::C, Sp::HOCOp}, {Sp::CHp, Sp::CO2}, 0.916728,
        1.41304e-12),  // C + HOCO+ -> CH+ + CO2
    rxn(503, {Sp::O, Sp::HOCOp}, {Sp::O2, Sp::HCOp}, 0.683259,
        3.37006e-13),  // O + HOCO+ -> O2 + HCO+
    rxn(504, {Sp::CH4, Sp::HOCOp}, {Sp::CH5p, Sp::CO2}, 0.944173,
        8.18978e-14),  // CH4 + HOCO+ -> CH5+ + CO2
    rxn(505, {Sp::H2O, Sp::HOCOp}, {Sp::CO2, Sp::H3Op}, 0.953449,
        2.50592e-12),  // H2O + HOCO+ -> CO2 + H3O+
    rxn(506, {Sp::CO, Sp::HOCOp}, {Sp::CO2, Sp::HCOp}, 0.981117,
        8.77329e-13),  // CO + HOCO+ -> CO2 + HCO+
    rxn(507, {Sp::e, Sp::Cp}, {Sp::C, PHOTON}, 2.79004e-41,
        1.8147e-11),  // e- + C+ -> C + γ
    rxn(508, {Sp::e, Sp::CHp}, {Sp::H, Sp::C}, 1.4387e-05,
        1.15062e-11),  // e- + CH+ -> H + C
    rxn(509, {Sp::e, Sp::CH2p}, {Sp::H2, Sp::C}, 5.6891e-06,
        1.13235e-11),  // e- + CH2+ -> H2 + C
    rxn(510, {Sp::e, Sp::CH2p}, {Sp::H, Sp::CH}, 1.42577e-05,
        9.70718e-12),  // e- + CH2+ -> H + CH
    rxn(511, {Sp::e, Sp::CH3p}, {Sp::H, Sp::CH2}, 1.41474e-05,
        8.18041e-12),  // e- + CH3+ -> H + CH2
    rxn(512, {Sp::e, Sp::CH3p}, {Sp::H2, Sp::CH}, 5.5943e-06,
        8.41677e-12),  // e- + CH3+ -> H2 + CH
    rxn(513, {Sp::e, Sp::CH3p}, {Sp::H, Sp::H, Sp::CH}, 7.22954e+30,
        1.24212e-12),  // e- + CH3+ -> H + H + CH
    rxn(514, {Sp::e, Sp::CH3p}, {Sp::CH3, PHOTON}, 2.79008e-41,
        1.57647e-11),  // e- + CH3+ -> CH3 + γ
    rxn(515, {Sp::e, Sp::Op}, {Sp::O, PHOTON}, 2.79009e-41,
        2.18189e-11),  // e- + O+ -> O + γ
    rxn(516, {Sp::e, Sp::CH4p}, {Sp::H, Sp::CH3}, 1.4052e-05,
        1.30323e-11),  // e- + CH4+ -> H + CH3
    rxn(517, {Sp::e, Sp::CH4p}, {Sp::H, Sp::H, Sp::CH2}, 7.1252e+30,
        5.4481e-12),  // e- + CH4+ -> H + H + CH2
    rxn(518, {Sp::e, Sp::OHp}, {Sp::H, Sp::O}, 1.39717e-05,
        1.37852e-11),  // e- + OH+ -> H + O
    rxn(519, {Sp::e, Sp::CH5p}, {Sp::H, Sp::CH4}, 1.39688e-05,
        1.28373e-11),  // e- + CH5+ -> H + CH4
    rxn(520, {Sp::e, Sp::CH5p}, {Sp::H2, Sp::CH3}, 5.44397e-06,
        1.28274e-11),  // e- + CH5+ -> H2 + CH3
    rxn(521, {Sp::e, Sp::H2Op}, {Sp::H, Sp::OH}, 1.38981e-05,
        1.20448e-11),  // e- + H2O+ -> H + OH
    rxn(522, {Sp::e, Sp::H2Op}, {Sp::H2, Sp::O}, 5.38542e-06,
        1.21502e-11),  // e- + H2O+ -> H2 + O
    rxn(523, {Sp::e, Sp::H3Op}, {Sp::H, Sp::H2O}, 1.38329e-05,
        1.04133e-11),  // e- + H3O+ -> H + H2O
    rxn(524, {Sp::e, Sp::H3Op}, {Sp::H, Sp::H, Sp::OH}, 6.89048e+30,
        2.24251e-12),  // e- + H3O+ -> H + H + OH
    rxn(525, {Sp::e, Sp::C2p}, {Sp::C, Sp::C}, 8.77611e-07,
        9.06495e-12),  // e- + C2+ -> C + C
    rxn(526, {Sp::e, Sp::COp}, {Sp::C, Sp::O}, 7.18465e-07,
        4.6553e-12),  // e- + CO+ -> C + O
    rxn(527, {Sp::e, Sp::HCOp}, {Sp::H, Sp::CO}, 1.34428e-05,
        1.20419e-11),  // e- + HCO+ -> H + CO
    rxn(528, {Sp::e, Sp::H2COp}, {Sp::H, Sp::HCO}, 1.34185e-05,
        1.14221e-11),  // e- + H2CO+ -> H + HCO
    rxn(529, {Sp::e, Sp::H2COp}, {Sp::H, Sp::H, Sp::CO}, 6.46498e+30,
        1.04209e-11),  // e- + H2CO+ -> H + H + CO
    rxn(530, {Sp::e, Sp::H2COp}, {Sp::H2CO, PHOTON}, 2.79015e-41,
        1.74545e-11),  // e- + H2CO+ -> H2CO + γ
    rxn(531, {Sp::e, Sp::H2COHp}, {Sp::H, Sp::H2, Sp::CO}, 2.40184e+30,
        1.02253e-11),  // e- + H2COH+ -> H + H2 + CO
    rxn(532, {Sp::e, Sp::H2COHp}, {Sp::H, Sp::H, Sp::HCO}, 6.44236e+30,
        4.05178e-12),  // e- + H2COH+ -> H + H + HCO
    rxn(533, {Sp::e, Sp::H2COHp}, {Sp::H, Sp::H2CO}, 1.33958e-05,
        1.00843e-11),  // e- + H2COH+ -> H + H2CO
    rxn(534, {Sp::e, Sp::O2p}, {Sp::O, Sp::O}, 5.70303e-07,
        1.11492e-11),  // e- + O2+ -> O + O
    rxn(535, {Sp::e, Sp::O2Hp}, {Sp::H, Sp::O2}, 1.33555e-05,
        1.4854e-11),  // e- + O2H+ -> H + O2
    rxn(536, {Sp::e, Sp::HOCOp}, {Sp::H, Sp::CO2}, 1.3189e-05,
        1.29192e-11),  // e- + HOCO+ -> H + CO2
    rxn(537, {Sp::e, Sp::HOCOp}, {Sp::OH, Sp::CO}, 3.74907e-07,
        1.12521e-11),  // e- + HOCO+ -> OH + CO
    rxn(538, {Sp::H2, Sp::C}, {Sp::CH2, PHOTON}, 4.90423e-36,
        5.32199e-12),  // H2 + C -> CH2 + γ
    rxn(539, {Sp::CH3, Sp::OH}, {Sp::CH4, Sp::O}, 0.994279,
        1.15291e-13),  // CH3 + OH -> CH4 + O
    rxn(540, {Sp::Hp, Sp::H2CO}, {Sp::H, Sp::H2, Sp::COp}, 1.79156e+35,
        -5.2488e-13),  // H+ + H2CO -> H + H2 + CO+
    rxn(541, {Sp::Hep, Sp::C}, {Sp::He, Sp::Cp}, 0.999863,
        2.12464e-11),  // He+ + C -> He + C+
    rxn(542, {Sp::Hep, Sp::H2CO}, {Sp::He, Sp::H2COp}, 0.999821,
        2.19388e-11),  // He+ + H2CO -> He + H2CO+
    rxn(543, {Sp::Hep, Sp::H2CO}, {Sp::He, Sp::CH2p, Sp::O}, 2.26391e+34,
        1.04132e-11),  // He+ + H2CO -> He + CH2+ + O
    rxn(601, {Sp::H, Sp::HCO}, {Sp::CH2, Sp::O}, 0.0470811,
        -6.30221e-12),  // H + HCO -> CH2 + O
    rxn(602, {Sp::H, Sp::H2O2}, {Sp::OH, Sp::H2O}, 0.0374399,
        4.78179e-12),  // H + H2O2 -> OH + H2O
    rxn(603, {Sp::C, Sp::CH2}, {Sp::CH, Sp::CH}, 0.991008,
        -1.37994e-12),  // C + CH2 -> CH + CH
    rxn(604, {Sp::CH, Sp::O2}, {Sp::OH, Sp::CO}, 0.817388,
        1.11125e-11),  // CH + O2 -> OH + CO
    rxn(605, {Sp::H, Sp::HCO}, {Sp::CH2, Sp::O}, 0.0470811,
        -6.30221e-12),  // H + HCO -> CH2 + O
    rxn(606, {Sp::CH2, Sp::O2}, {Sp::H, Sp::H, Sp::CO2}, 1.46944e+37,
        5.84134e-12),  // CH2 + O2 -> H + H + CO2
    rxn(607, {Sp::CH2, Sp::O2}, {Sp::H2, Sp::CO2}, 11.3707,
        1.3016e-11),  // CH2 + O2 -> H2 + CO2
    rxn(608, {Sp::CH2, Sp::O2}, {Sp::CO, Sp::H2O}, 0.838549,
        1.2345e-11),  // CH2 + O2 -> CO + H2O
    rxn(609, {Sp::CH3, Sp::CH3}, {Sp::CH2, Sp::CH4}, 1.00679,
        -3.99642e-13),  // CH3 + CH3 -> CH2 + CH4
    rxn(610, {Sp::CH3, Sp::O}, {Sp::H, Sp::H2, Sp::CO}, 4.01523e+36,
        4.89143e-12),  // CH3 + O -> H + H2 + CO
    rxn(611, {Sp::CH3, Sp::OH}, {Sp::H2, Sp::H2CO}, 8.67761,
        4.85578e-12),  // CH3 + OH -> H2 + H2CO
    rxn(612, {Sp::CH3, Sp::O2}, {Sp::H2O, Sp::HCO}, 0.88252,
        5.76197e-12),  // CH3 + O2 -> H2O + HCO
    rxn(613, {Sp::O, Sp::H2CO}, {Sp::H, Sp::OH, Sp::CO}, 4.62711e+35,
        3.56518e-14),  // O + H2CO -> H + OH + CO
    rxn(614, {Sp::C2, Sp::O2}, {Sp::CO, Sp::CO}, 0.969611,
        1.74014e-11),  // C2 + O2 -> CO + CO
    rxn(615, {Sp::HCO, Sp::HCO}, {Sp::H2, Sp::CO, Sp::CO}, 1.79623e+35,
        5.17228e-12),  // HCO + HCO -> H2 + CO + CO
    rxn(616, {Sp::O2, Sp::HCO}, {Sp::OH, Sp::CO2}, 1.38156,
        6.60845e-12),  // O2 + HCO -> OH + CO2
    rxn(617, {Sp::O, Sp::OH}, {Sp::H, Sp::O2}, 24.4992,
        1.12678e-12),  // O + OH -> H + O2
    rxn(621, {Sp::H3p, Sp::O}, {Sp::H, Sp::H2Op}, 4.34773,
        2.69952e-12),  // H3+ + O -> H + H2O+
    rxn(626, {Sp::CO, Sp::Op}, {Sp::O, Sp::COp}, 0.999978,
        -6.34127e-13),  // CO + O+ -> O + CO+
    rxn(632, {Sp::e, Sp::CH2p}, {Sp::H, Sp::H, Sp::C}, 7.35205e+30,
        4.14882e-12),  // e- + CH2+ -> H + H + C
    rxn(633, {Sp::e, Sp::CH5p}, {Sp::H, Sp::H, Sp::CH3}, 7.03527e+30,
        5.65271e-12),  // e- + CH5+ -> H + H + CH3
    rxn(634, {Sp::e, Sp::CH5p}, {Sp::H, Sp::H2, Sp::CH2}, 2.76042e+30,
        5.24312e-12),  // e- + CH5+ -> H + H2 + CH2
    rxn(635, {Sp::e, Sp::CH5p}, {Sp::H2, Sp::H2, Sp::CH}, 1.09155e+30,
        5.47948e-12),  // e- + CH5+ -> H2 + H2 + CH
    rxn(636, {Sp::e, Sp::H2Op}, {Sp::H, Sp::H, Sp::O}, 6.95961e+30,
        4.97552e-12),  // e- + H2O+ -> H + H + O
    rxn(637, {Sp::e, Sp::H3Op}, {Sp::H, Sp::H2, Sp::O}, 2.67001e+30,
        2.34785e-12),  // e- + H3O+ -> H + H2 + O
    rxn(638, {Sp::e, Sp::H3Op}, {Sp::H2, Sp::OH}, 5.33193e-06,
        9.41717e-12),  // e- + H3O+ -> H2 + OH
    rxn(640, {Sp::e, Sp::HOCOp}, {Sp::H, Sp::O, Sp::CO}, 1.87738e+29,
        4.1828e-12),  // e- + HOCO+ -> H + O + CO
    rxn(641, {Sp::Hp, Sp::He}, {Sp::HeHp, PHOTON}, 1.56169e-36,
        2.90085e-12),  // H+ + He -> HeH+ + γ
    rxn(642, {Sp::H, Sp::OH}, {Sp::H2O, PHOTON}, 2.00754e-36,
        8.17078e-12),  // H + OH -> H2O + γ
    rxn(643, {Sp::H2, Sp::CH}, {Sp::CH3, PHOTON}, 4.98735e-36,
        7.34789e-12),  // H2 + CH -> CH3 + γ
    rxn(644, {Sp::C, Sp::Cp}, {Sp::C2p, PHOTON}, 3.17913e-35,
        9.08202e-12),  // C + C+ -> C2+ + γ
    rxn(645, {Sp::C, Sp::Op}, {Sp::COp, PHOTON}, 3.8834e-35,
        1.71635e-11),  // C + O+ -> CO+ + γ
    rxn(701, {Sp::e, Sp::Kp}, {Sp::K, PHOTON}, 2.79017e-41,
        6.95434e-12),  // e- + K+ -> K + γ
    rxn(702, {Sp::H2, Sp::K}, {Sp::H2, Sp::e, Sp::Kp}, 3.58401e+40,
        -6.95434e-12),  // H2 + K -> H2 + e- + K+
    rxn(703, {Sp::e, Sp::Nap}, {Sp::Na, PHOTON}, 2.79013e-41,
        8.23461e-12),  // e- + Na+ -> Na + γ
    rxn(704, {Sp::Hp, Sp::Na}, {Sp::H, Sp::Nap}, 0.999217,
        1.35525e-11),  // H+ + Na -> H + Na+
    rxn(705, {Sp::Hm, Sp::Nap}, {Sp::H, Sp::Na}, 1.00077,
        7.02594e-12),  // H- + Na+ -> H + Na
    rxn(706, {Sp::H3p, Sp::Na}, {Sp::H, Sp::H2, Sp::Nap}, 8.39186e+35,
        6.61509e-12),  // H3+ + Na -> H + H2 + Na+
    rxn(707, {Sp::Cp, Sp::Na}, {Sp::C, Sp::Nap}, 0.999967,
        9.91236e-12),  // C+ + Na -> C + Na+
    rxn(708, {Sp::CHp, Sp::Na}, {Sp::CH, Sp::Nap}, 0.999972,
        8.82992e-12),  // CH+ + Na -> CH + Na+
    rxn(709, {Sp::CH3p, Sp::Na}, {Sp::CH3, Sp::Nap}, 0.999981,
        7.53005e-12),  // CH3+ + Na -> CH3 + Na+
    rxn(710, {Sp::O2p, Sp::Na}, {Sp::O2, Sp::Nap}, 1.00001,
        1.11107e-11),  // O2+ + Na -> O2 + Na+
    rxn(711, {Sp::H2Op, Sp::Na}, {Sp::H2O, Sp::Nap}, 0.999991,
        1.1981e-11),  // H2O+ + Na -> H2O + Na+
    rxn(712, {Sp::H3Op, Sp::Na}, {Sp::H, Sp::H2O, Sp::Nap}, 4.95781e+35,
        2.17868e-12),  // H3O+ + Na -> H + H2O + Na+
    rxn(713, {Sp::HCOp, Sp::Na}, {Sp::HCO, Sp::Nap}, 1.00001,
        4.80846e-12),  // HCO+ + Na -> HCO + Na+
    rxn(714, {Sp::H2COp, Sp::Na}, {Sp::H2CO, Sp::Nap}, 1.00001,
        9.21991e-12),  // H2CO+ + Na -> H2CO + Na+
    rxn(715, {Sp::H2COHp, Sp::Na}, {Sp::H, Sp::H2CO, Sp::Nap}, 4.80115e+35,
        1.84964e-12),  // H2COH+ + Na -> H + H2CO + Na+
    rxn(716, {Sp::Na, Sp::Mgp}, {Sp::Nap, Sp::Mg}, 1.0,
        4.01601e-12),  // Na + Mg+ -> Na+ + Mg
    rxn(717, {Sp::H2, Sp::Na}, {Sp::H2, Sp::e, Sp::Nap}, 3.58406e+40,
        -8.23461e-12),  // H2 + Na -> H2 + e- + Na+
    rxn(718, {Sp::e, Sp::Mgp}, {Sp::Mg, PHOTON}, 2.79014e-41,
        1.22506e-11),  // e- + Mg+ -> Mg + γ
    rxn(719, {Sp::Hp, Sp::Mg}, {Sp::H, Sp::Mgp}, 0.999215,
        9.53648e-12),  // H+ + Mg -> H + Mg+
    rxn(720, {Sp::Hm, Sp::Mgp}, {Sp::H, Sp::Mg}, 1.00077,
        1.1042e-11),  // H- + Mg+ -> H + Mg
    rxn(721, {Sp::H3p, Sp::Mg}, {Sp::H, Sp::H2, Sp::Mgp}, 8.39184e+35,
        2.59908e-12),  // H3+ + Mg -> H + H2 + Mg+
    rxn(722, {Sp::Cp, Sp::Mg}, {Sp::C, Sp::Mgp}, 0.999965,
        5.89634e-12),  // C+ + Mg -> C + Mg+
    rxn(723, {Sp::CHp, Sp::Mg}, {Sp::CH, Sp::Mgp}, 0.999971,
        4.8139e-12),  // CH+ + Mg -> CH + Mg+
    rxn(724, {Sp::CH3p, Sp::Mg}, {Sp::CH3, Sp::Mgp}, 0.999979,
        3.51403e-12),  // CH3+ + Mg -> CH3 + Mg+
    rxn(725, {Sp::CH5p, Sp::Mg}, {Sp::H, Sp::CH4, Sp::Mgp}, 5.00651e+35,
        5.86685e-13),  // CH5+ + Mg -> H + CH4 + Mg+
    rxn(726, {Sp::O2p, Sp::Mg}, {Sp::O2, Sp::Mgp}, 1.00001,
        7.09465e-12),  // O2+ + Mg -> O2 + Mg+
    rxn(727, {Sp::H2Op, Sp::Mg}, {Sp::H2O, Sp::Mgp}, 0.999989,
        7.96499e-12),  // H2O+ + Mg -> H2O + Mg+
    rxn(728, {Sp::H3Op, Sp::Mg}, {Sp::H, Sp::H2O, Sp::Mgp}, 4.9578e+35,
        -1.83734e-12),  // H3O+ + Mg -> H + H2O + Mg+
    rxn(729, {Sp::HCOp, Sp::Mg}, {Sp::HCO, Sp::Mgp}, 1.00001,
        7.92442e-13),  // HCO+ + Mg -> HCO + Mg+
    rxn(730, {Sp::H2COp, Sp::Mg}, {Sp::H2CO, Sp::Mgp}, 1.00001,
        5.2039e-12),  // H2CO+ + Mg -> H2CO + Mg+
    rxn(801, {Sp::e, Sp::Lip}, {Sp::Li, PHOTON}, 2.78963e-41,
        8.63849e-12),  // e- + Li+ -> Li + γ
    rxn(802, {Sp::Hm, Sp::Lip}, {Sp::H, Sp::Li}, 1.00059,
        7.42982e-12),  // H- + Li+ -> H + Li
    rxn(803, {Sp::Hp, Sp::Lim}, {Sp::H, Sp::Li}, 0.999397,
        2.0797e-11),  // H+ + Li- -> H + Li
    rxn(804, {Sp::e, Sp::Li}, {Sp::Lim, PHOTON}, 2.78963e-41,
        9.90146e-13),  // e- + Li -> Li- + γ
    rxn(805, {Sp::Hp, Sp::Li}, {Sp::H, Sp::Lip}, 0.999397,
        1.31486e-11),  // H+ + Li -> H + Li+
    rxn(806, {Sp::Hp, Sp::Li}, {Sp::H, Sp::Lip, PHOTON}, 0.999397,
        1.31486e-11),  // H+ + Li -> H + Li+ + γ
    rxn(807, {Sp::Hm, Sp::Li}, {Sp::e, Sp::LiH}, 64054.6,
        2.6598e-12),  // H- + Li -> e- + LiH
    rxn(808, {Sp::H, Sp::Lim}, {Sp::e, Sp::LiH}, 64017.0,
        2.87833e-12),  // H + Li- -> e- + LiH
    rxn(809, {Sp::H, Sp::LiHp}, {Sp::Hp, Sp::LiH}, 1.00082,
        -9.5123e-12),  // H + LiH+ -> H+ + LiH
    rxn(810, {Sp::Hp, Sp::LiH}, {Sp::H, Sp::LiHp}, 0.999182,
        9.5123e-12),  // H+ + LiH -> H + LiH+
    rxn(811, {Sp::H, Sp::LiH}, {Sp::H2, Sp::Li}, 0.433305,
        3.30618e-12),  // H + LiH -> H2 + Li
    rxn(812, {Sp::H, Sp::Lip}, {Sp::LiHp, PHOTON}, 1.78545e-36,
        2.3216e-13),  // H + Li+ -> LiH+ + γ
    rxn(813, {Sp::Hp, Sp::Li}, {Sp::LiHp, PHOTON}, 1.78437e-36,
        1.33808e-11),  // H+ + Li -> LiH+ + γ
    rxn(814, {Sp::Hp, Sp::LiH}, {Sp::H2, Sp::Lip}, 0.433044,
        1.64548e-11),  // H+ + LiH -> H2 + Li+
    rxn(815, {Sp::e, Sp::LiHp}, {Sp::H, Sp::Li}, 1.56242e-05,
        8.40633e-12),  // e- + LiH+ -> H + Li
    rxn(816, {Sp::H, Sp::LiHp}, {Sp::H2p, Sp::Li}, 0.433483,
        -9.13389e-12),  // H + LiH+ -> H2+ + Li
    rxn(817, {Sp::H, Sp::LiHp}, {Sp::H2, Sp::Lip}, 0.433399,
        6.9425e-12),  // H + LiH+ -> H2 + Li+
    rxn(818, {Sp::H, Sp::Li}, {Sp::LiH, PHOTON}, 1.78583e-36,
        3.86847e-12),  // H + Li -> LiH + γ
    rxn(819, {Sp::H2p, Sp::Li}, {Sp::Hp, Sp::LiH}, 2.30879,
        -3.78404e-13),  // H2+ + Li -> H+ + LiH
    rxn(820, {Sp::Hp, Sp::LiH}, {Sp::H2p, Sp::Li}, 0.433128,
        3.78404e-13),  // H+ + LiH -> H2+ + Li
    rxn(821, {Sp::e, Sp::Lipp}, {Sp::Lip, PHOTON}, 2.79023e-41,
        1.21191e-10),  // e- + Li++ -> Li+ + γ
    rxn(822, {Sp::e, Sp::Lippp}, {Sp::Lipp, PHOTON}, 2.79023e-41,
        1.96196e-10),  // e- + Li+++ -> Li++ + γ
    rxn(823, {Sp::Dm, Sp::Lip}, {Sp::D, Sp::Li}, 1.00019,
        7.4294e-12),  // D- + Li+ -> D + Li
    rxn(824, {Sp::Dp, Sp::Lim}, {Sp::D, Sp::Li}, 0.999806,
        2.08029e-11),  // D+ + Li- -> D + Li
    rxn(825, {Sp::e, Sp::Li}, {Sp::e, Sp::e, Sp::Lip}, 3.58471e+40,
        -8.63849e-12),  // e- + Li -> e- + e- + Li+
    rxn(826, {Sp::e, Sp::Lip}, {Sp::e, Sp::e, Sp::Lipp}, 3.58393e+40,
        -1.21191e-10),  // e- + Li+ -> e- + e- + Li++
    rxn(827, {Sp::e, Sp::Lipp}, {Sp::e, Sp::e, Sp::Lippp}, 3.58393e+40,
        -1.96196e-10),  // e- + Li++ -> e- + e- + Li+++
    rxn(828, {Sp::H, Sp::H, Sp::Li}, {Sp::H, Sp::LiH}, 1.78583e-36,
        3.86847e-12),  // H + H + Li -> H + LiH
    rxn(829, {Sp::H, Sp::H2, Sp::Li}, {Sp::H2, Sp::LiH}, 1.78583e-36,
        3.86847e-12),  // H + H2 + Li -> H2 + LiH
    rxn(830, {Sp::H2, Sp::Li}, {Sp::H2, Sp::e, Sp::Lip}, 3.58471e+40,
        -8.63849e-12),  // H2 + Li -> H2 + e- + Li+
    rxn(544, {Sp::H, CR}, {Sp::e, Sp::Hp}, 0.0, 0.0),    // H + CR -> e- + H+
    rxn(545, {Sp::He, CR}, {Sp::e, Sp::Hep}, 0.0, 0.0),  // He + CR -> e- + He+
    rxn(546, {Sp::C, CR}, {Sp::e, Sp::Cp}, 0.0, 0.0),    // C + CR -> e- + C+
    rxn(547, {Sp::O, CR}, {Sp::e, Sp::Op}, 0.0, 0.0),    // O + CR -> e- + O+
    rxn(548, {Sp::H2, CR}, {Sp::H, Sp::e, Sp::Hp}, 0.0,
        0.0),  // H2 + CR -> H + e- + H+
    rxn(549, {Sp::H2, CR}, {Sp::e, Sp::H2p}, 0.0, 0.0),  // H2 + CR -> e- + H2+
    rxn(550, {Sp::H2, CR}, {Sp::H, Sp::H}, 0.0, 0.0),    // H2 + CR -> H + H
    rxn(551, {Sp::H2, CR}, {Sp::Hp, Sp::Hm}, 0.0, 0.0),  // H2 + CR -> H+ + H-
    rxn(552, {Sp::CO, CR}, {Sp::e, Sp::COp}, 0.0, 0.0),  // CO + CR -> e- + CO+
    rxn(656, {Sp::C, CR}, {Sp::e, Sp::Cp}, 0.0, 0.0),    // C + CR -> e- + C+
    rxn(657, {Sp::CH, CR}, {Sp::H, Sp::C}, 0.0, 0.0),    // CH + CR -> H + C
    rxn(658, {Sp::CHp, CR}, {Sp::H, Sp::Cp}, 0.0, 0.0),  // CH+ + CR -> H + C+
    rxn(659, {Sp::CH2, CR}, {Sp::e, Sp::CH2p}, 0.0,
        0.0),  // CH2 + CR -> e- + CH2+
    rxn(660, {Sp::CH2, CR}, {Sp::H, Sp::CH}, 0.0, 0.0),  // CH2 + CR -> H + CH
    rxn(661, {Sp::CH3, CR}, {Sp::e, Sp::CH3p}, 0.0,
        0.0),  // CH3 + CR -> e- + CH3+
    rxn(662, {Sp::CH3, CR}, {Sp::H, Sp::CH2}, 0.0, 0.0),  // CH3 + CR -> H + CH2
    rxn(663, {Sp::CH3, CR}, {Sp::H2, Sp::CH}, 0.0, 0.0),  // CH3 + CR -> H2 + CH
    rxn(664, {Sp::CH4, CR}, {Sp::H2, Sp::CH2}, 0.0,
        0.0),                                            // CH4 + CR -> H2 + CH2
    rxn(665, {Sp::OH, CR}, {Sp::H, Sp::O}, 0.0, 0.0),    // OH + CR -> H + O
    rxn(666, {Sp::H2O, CR}, {Sp::H, Sp::OH}, 0.0, 0.0),  // H2O + CR -> H + OH
    rxn(667, {Sp::C2, CR}, {Sp::C, Sp::C}, 0.0, 0.0),    // C2 + CR -> C + C
    rxn(668, {Sp::CO, CR}, {Sp::C, Sp::O}, 0.0, 0.0),    // CO + CR -> C + O
    rxn(669, {Sp::HCO, CR}, {Sp::H, Sp::CO}, 0.0, 0.0),  // HCO + CR -> H + CO
    rxn(670, {Sp::HCO, CR}, {Sp::e, Sp::HCOp}, 0.0,
        0.0),  // HCO + CR -> e- + HCO+
    rxn(671, {Sp::H2CO, CR}, {Sp::H2, Sp::CO}, 0.0,
        0.0),                                            // H2CO + CR -> H2 + CO
    rxn(672, {Sp::O2, CR}, {Sp::O, Sp::O}, 0.0, 0.0),    // O2 + CR -> O + O
    rxn(673, {Sp::O2, CR}, {Sp::e, Sp::O2p}, 0.0, 0.0),  // O2 + CR -> e- + O2+
    rxn(674, {Sp::H2O2, CR}, {Sp::OH, Sp::OH}, 0.0,
        0.0),                                            // H2O2 + CR -> OH + OH
    rxn(675, {Sp::CO2, CR}, {Sp::O, Sp::CO}, 0.0, 0.0),  // CO2 + CR -> O + CO
    rxn(676, {Sp::Hm, CR}, {Sp::H, Sp::e}, 0.0, 0.0),    // H- + CR -> H + e-
    rxn(677, {Sp::H, CR}, {Sp::e, Sp::Hp}, 0.0, 0.0),    // H + CR -> e- + H+
    rxn(678, {Sp::He, CR}, {Sp::e, Sp::Hep}, 0.0, 0.0),  // He + CR -> e- + He+
    rxn(679, {Sp::O, CR}, {Sp::e, Sp::Op}, 0.0, 0.0),    // O + CR -> e- + O+
    rxn(680, {Sp::O2H, CR}, {Sp::H, Sp::O2}, 0.0, 0.0),  // O2H + CR -> H + O2
    rxn(681, {Sp::Na, CR}, {Sp::e, Sp::Nap}, 0.0, 0.0),  // Na + CR -> e- + Na+
    rxn(682, {Sp::Mg, CR}, {Sp::e, Sp::Mgp}, 0.0, 0.0),  // Mg + CR -> e- + Mg+
    rxn(841, {Sp::Hp, Sp::Gr}, {Sp::H, Sp::Grp}, 0.0,
        0.0),  // H+ + Gr -> H + Gr+
    rxn(842, {Sp::H2p, Sp::Gr}, {Sp::H, Sp::H, Sp::Grp}, 0.0,
        0.0),  // H2+ + Gr -> H + H + Gr+
    rxn(843, {Sp::H3p, Sp::Gr}, {Sp::H, Sp::H2, Sp::Grp}, 0.0,
        0.0),  // H3+ + Gr -> H + H2 + Gr+
    rxn(844, {Sp::Hep, Sp::Gr}, {Sp::He, Sp::Grp}, 0.0,
        0.0),  // He+ + Gr -> He + Gr+
    rxn(845, {Sp::HeHp, Sp::Gr}, {Sp::H, Sp::He, Sp::Grp}, 0.0,
        0.0),  // HeH+ + Gr -> H + He + Gr+
    rxn(846, {Sp::Dp, Sp::Gr}, {Sp::D, Sp::Grp}, 0.0,
        0.0),  // D+ + Gr -> D + Gr+
    rxn(847, {Sp::HDp, Sp::Gr}, {Sp::H, Sp::D, Sp::Grp}, 0.0,
        0.0),  // HD+ + Gr -> H + D + Gr+
    rxn(848, {Sp::Cp, Sp::Gr}, {Sp::C, Sp::Grp}, 0.0,
        0.0),  // C+ + Gr -> C + Gr+
    rxn(849, {Sp::C2p, Sp::Gr}, {Sp::C, Sp::C, Sp::Grp}, 0.0,
        0.0),  // C2+ + Gr -> C + C + Gr+
    rxn(850, {Sp::CHp, Sp::Gr}, {Sp::H, Sp::C, Sp::Grp}, 0.0,
        0.0),  // CH+ + Gr -> H + C + Gr+
    rxn(851, {Sp::CH2p, Sp::Gr}, {Sp::H2, Sp::C, Sp::Grp}, 0.0,
        0.0),  // CH2+ + Gr -> H2 + C + Gr+
    rxn(852, {Sp::CH3p, Sp::Gr}, {Sp::H2, Sp::CH, Sp::Grp}, 0.0,
        0.0),  // CH3+ + Gr -> H2 + CH + Gr+
    rxn(853, {Sp::CH4p, Sp::Gr}, {Sp::H, Sp::CH3, Sp::Grp}, 0.0,
        0.0),  // CH4+ + Gr -> H + CH3 + Gr+
    rxn(854, {Sp::CH5p, Sp::Gr}, {Sp::H2, Sp::CH3, Sp::Grp}, 0.0,
        0.0),  // CH5+ + Gr -> H2 + CH3 + Gr+
    rxn(855, {Sp::Op, Sp::Gr}, {Sp::O, Sp::Grp}, 0.0,
        0.0),  // O+ + Gr -> O + Gr+
    rxn(856, {Sp::O2p, Sp::Gr}, {Sp::O, Sp::O, Sp::Grp}, 0.0,
        0.0),  // O2+ + Gr -> O + O + Gr+
    rxn(857, {Sp::OHp, Sp::Gr}, {Sp::H, Sp::O, Sp::Grp}, 0.0,
        0.0),  // OH+ + Gr -> H + O + Gr+
    rxn(858, {Sp::COp, Sp::Gr}, {Sp::C, Sp::O, Sp::Grp}, 0.0,
        0.0),  // CO+ + Gr -> C + O + Gr+
    rxn(859, {Sp::H2Op, Sp::Gr}, {Sp::H2, Sp::O, Sp::Grp}, 0.0,
        0.0),  // H2O+ + Gr -> H2 + O + Gr+
    rxn(860, {Sp::HCOp, Sp::Gr}, {Sp::H, Sp::CO, Sp::Grp}, 0.0,
        0.0),  // HCO+ + Gr -> H + CO + Gr+
    rxn(861, {Sp::O2Hp, Sp::Gr}, {Sp::H, Sp::O2, Sp::Grp}, 0.0,
        0.0),  // O2H+ + Gr -> H + O2 + Gr+
    rxn(862, {Sp::H3Op, Sp::Gr}, {Sp::H2, Sp::OH, Sp::Grp}, 0.0,
        0.0),  // H3O+ + Gr -> H2 + OH + Gr+
    rxn(863, {Sp::H2COp, Sp::Gr}, {Sp::H2CO, Sp::Grp}, 0.0,
        0.0),  // H2CO+ + Gr -> H2CO + Gr+
    rxn(864, {Sp::HOCOp, Sp::Gr}, {Sp::OH, Sp::CO, Sp::Grp}, 0.0,
        0.0),  // HOCO+ + Gr -> OH + CO + Gr+
    rxn(865, {Sp::H2COHp, Sp::Gr}, {Sp::H, Sp::H2CO, Sp::Grp}, 0.0,
        0.0),  // H2COH+ + Gr -> H + H2CO + Gr+
    rxn(866, {Sp::Lip, Sp::Gr}, {Sp::Li, Sp::Grp}, 0.0,
        0.0),  // Li+ + Gr -> Li + Gr+
    rxn(867, {Sp::LiHp, Sp::Gr}, {Sp::H, Sp::Li, Sp::Grp}, 0.0,
        0.0),  // LiH+ + Gr -> H + Li + Gr+
    rxn(868, {Sp::Kp, Sp::Gr}, {Sp::K, Sp::Grp}, 0.0,
        0.0),  // K+ + Gr -> K + Gr+
    rxn(869, {Sp::Nap, Sp::Gr}, {Sp::Na, Sp::Grp}, 0.0,
        0.0),  // Na+ + Gr -> Na + Gr+
    rxn(870, {Sp::Mgp, Sp::Gr}, {Sp::Mg, Sp::Grp}, 0.0,
        0.0),  // Mg+ + Gr -> Mg + Gr+
    rxn(871, {Sp::Hepp, Sp::Gr}, {Sp::He, Sp::Gr2p}, 0.0,
        0.0),  // He++ + Gr -> He + Gr2+
    rxn(872, {Sp::Lipp, Sp::Gr}, {Sp::Li, Sp::Gr2p}, 0.0,
        0.0),                                        // Li++ + Gr -> Li + Gr2+
    rxn(873, {Sp::e, Sp::Gr}, {Sp::Grm}, 0.0, 0.0),  // e- + Gr -> Gr-
    rxn(874, {Sp::Hm, Sp::Gr}, {Sp::H, Sp::Grm}, 0.0,
        0.0),  // H- + Gr -> H + Gr-
    rxn(875, {Sp::Dm, Sp::Gr}, {Sp::D, Sp::Grm}, 0.0,
        0.0),  // D- + Gr -> D + Gr-
    rxn(876, {Sp::Lim, Sp::Gr}, {Sp::Li, Sp::Grm}, 0.0,
        0.0),  // Li- + Gr -> Li + Gr-
    rxn(877, {Sp::Hp, Sp::Grp}, {Sp::H, Sp::Gr2p}, 0.0,
        0.0),  // H+ + Gr+ -> H + Gr2+
    rxn(878, {Sp::H2p, Sp::Grp}, {Sp::H, Sp::H, Sp::Gr2p}, 0.0,
        0.0),  // H2+ + Gr+ -> H + H + Gr2+
    rxn(879, {Sp::H3p, Sp::Grp}, {Sp::H, Sp::H2, Sp::Gr2p}, 0.0,
        0.0),  // H3+ + Gr+ -> H + H2 + Gr2+
    rxn(880, {Sp::Hep, Sp::Grp}, {Sp::He, Sp::Gr2p}, 0.0,
        0.0),  // He+ + Gr+ -> He + Gr2+
    rxn(881, {Sp::HeHp, Sp::Grp}, {Sp::H, Sp::He, Sp::Gr2p}, 0.0,
        0.0),  // HeH+ + Gr+ -> H + He + Gr2+
    rxn(882, {Sp::Dp, Sp::Grp}, {Sp::D, Sp::Gr2p}, 0.0,
        0.0),  // D+ + Gr+ -> D + Gr2+
    rxn(883, {Sp::HDp, Sp::Grp}, {Sp::H, Sp::D, Sp::Gr2p}, 0.0,
        0.0),  // HD+ + Gr+ -> H + D + Gr2+
    rxn(884, {Sp::Cp, Sp::Grp}, {Sp::C, Sp::Gr2p}, 0.0,
        0.0),  // C+ + Gr+ -> C + Gr2+
    rxn(885, {Sp::C2p, Sp::Grp}, {Sp::C, Sp::C, Sp::Gr2p}, 0.0,
        0.0),  // C2+ + Gr+ -> C + C + Gr2+
    rxn(886, {Sp::CHp, Sp::Grp}, {Sp::H, Sp::C, Sp::Gr2p}, 0.0,
        0.0),  // CH+ + Gr+ -> H + C + Gr2+
    rxn(887, {Sp::CH2p, Sp::Grp}, {Sp::H2, Sp::C, Sp::Gr2p}, 0.0,
        0.0),  // CH2+ + Gr+ -> H2 + C + Gr2+
    rxn(888, {Sp::CH3p, Sp::Grp}, {Sp::H2, Sp::CH, Sp::Gr2p}, 0.0,
        0.0),  // CH3+ + Gr+ -> H2 + CH + Gr2+
    rxn(889, {Sp::CH4p, Sp::Grp}, {Sp::H, Sp::CH3, Sp::Gr2p}, 0.0,
        0.0),  // CH4+ + Gr+ -> H + CH3 + Gr2+
    rxn(890, {Sp::CH5p, Sp::Grp}, {Sp::H2, Sp::CH3, Sp::Gr2p}, 0.0,
        0.0),  // CH5+ + Gr+ -> H2 + CH3 + Gr2+
    rxn(891, {Sp::Op, Sp::Grp}, {Sp::O, Sp::Gr2p}, 0.0,
        0.0),  // O+ + Gr+ -> O + Gr2+
    rxn(892, {Sp::O2p, Sp::Grp}, {Sp::O, Sp::O, Sp::Gr2p}, 0.0,
        0.0),  // O2+ + Gr+ -> O + O + Gr2+
    rxn(893, {Sp::OHp, Sp::Grp}, {Sp::H, Sp::O, Sp::Gr2p}, 0.0,
        0.0),  // OH+ + Gr+ -> H + O + Gr2+
    rxn(894, {Sp::COp, Sp::Grp}, {Sp::C, Sp::O, Sp::Gr2p}, 0.0,
        0.0),  // CO+ + Gr+ -> C + O + Gr2+
    rxn(895, {Sp::H2Op, Sp::Grp}, {Sp::H2, Sp::O, Sp::Gr2p}, 0.0,
        0.0),  // H2O+ + Gr+ -> H2 + O + Gr2+
    rxn(896, {Sp::HCOp, Sp::Grp}, {Sp::H, Sp::CO, Sp::Gr2p}, 0.0,
        0.0),  // HCO+ + Gr+ -> H + CO + Gr2+
    rxn(897, {Sp::O2Hp, Sp::Grp}, {Sp::H, Sp::O2, Sp::Gr2p}, 0.0,
        0.0),  // O2H+ + Gr+ -> H + O2 + Gr2+
    rxn(898, {Sp::H3Op, Sp::Grp}, {Sp::H2, Sp::OH, Sp::Gr2p}, 0.0,
        0.0),  // H3O+ + Gr+ -> H2 + OH + Gr2+
    rxn(899, {Sp::H2COp, Sp::Grp}, {Sp::H2CO, Sp::Gr2p}, 0.0,
        0.0),  // H2CO+ + Gr+ -> H2CO + Gr2+
    rxn(900, {Sp::HOCOp, Sp::Grp}, {Sp::OH, Sp::CO, Sp::Gr2p}, 0.0,
        0.0),  // HOCO+ + Gr+ -> OH + CO + Gr2+
    rxn(901, {Sp::H2COHp, Sp::Grp}, {Sp::H, Sp::H2CO, Sp::Gr2p}, 0.0,
        0.0),  // H2COH+ + Gr+ -> H + H2CO + Gr2+
    rxn(902, {Sp::Lip, Sp::Grp}, {Sp::Li, Sp::Gr2p}, 0.0,
        0.0),  // Li+ + Gr+ -> Li + Gr2+
    rxn(903, {Sp::LiHp, Sp::Grp}, {Sp::H, Sp::Li, Sp::Gr2p}, 0.0,
        0.0),  // LiH+ + Gr+ -> H + Li + Gr2+
    rxn(904, {Sp::Kp, Sp::Grp}, {Sp::K, Sp::Gr2p}, 0.0,
        0.0),  // K+ + Gr+ -> K + Gr2+
    rxn(905, {Sp::Nap, Sp::Grp}, {Sp::Na, Sp::Gr2p}, 0.0,
        0.0),  // Na+ + Gr+ -> Na + Gr2+
    rxn(906, {Sp::Mgp, Sp::Grp}, {Sp::Mg, Sp::Gr2p}, 0.0,
        0.0),                                        // Mg+ + Gr+ -> Mg + Gr2+
    rxn(907, {Sp::e, Sp::Grp}, {Sp::Gr}, 0.0, 0.0),  // e- + Gr+ -> Gr
    rxn(908, {Sp::Hm, Sp::Grp}, {Sp::H, Sp::Gr}, 0.0,
        0.0),  // H- + Gr+ -> H + Gr
    rxn(909, {Sp::Dm, Sp::Grp}, {Sp::D, Sp::Gr}, 0.0,
        0.0),  // D- + Gr+ -> D + Gr
    rxn(910, {Sp::Lim, Sp::Grp}, {Sp::Li, Sp::Gr}, 0.0,
        0.0),                                          // Li- + Gr+ -> Li + Gr
    rxn(911, {Sp::e, Sp::Gr2p}, {Sp::Grp}, 0.0, 0.0),  // e- + Gr2+ -> Gr+
    rxn(912, {Sp::Hm, Sp::Gr2p}, {Sp::H, Sp::Grp}, 0.0,
        0.0),  // H- + Gr2+ -> H + Gr+
    rxn(913, {Sp::Dm, Sp::Gr2p}, {Sp::D, Sp::Grp}, 0.0,
        0.0),  // D- + Gr2+ -> D + Gr+
    rxn(914, {Sp::Lim, Sp::Gr2p}, {Sp::Li, Sp::Grp}, 0.0,
        0.0),  // Li- + Gr2+ -> Li + Gr+
    rxn(915, {Sp::Hp, Sp::Grm}, {Sp::H, Sp::Gr}, 0.0,
        0.0),  // H+ + Gr- -> H + Gr
    rxn(916, {Sp::H2p, Sp::Grm}, {Sp::H, Sp::H, Sp::Gr}, 0.0,
        0.0),  // H2+ + Gr- -> H + H + Gr
    rxn(917, {Sp::H3p, Sp::Grm}, {Sp::H, Sp::H2, Sp::Gr}, 0.0,
        0.0),  // H3+ + Gr- -> H + H2 + Gr
    rxn(918, {Sp::Hep, Sp::Grm}, {Sp::He, Sp::Gr}, 0.0,
        0.0),  // He+ + Gr- -> He + Gr
    rxn(919, {Sp::HeHp, Sp::Grm}, {Sp::H, Sp::He, Sp::Gr}, 0.0,
        0.0),  // HeH+ + Gr- -> H + He + Gr
    rxn(920, {Sp::Dp, Sp::Grm}, {Sp::D, Sp::Gr}, 0.0,
        0.0),  // D+ + Gr- -> D + Gr
    rxn(921, {Sp::HDp, Sp::Grm}, {Sp::H, Sp::D, Sp::Gr}, 0.0,
        0.0),  // HD+ + Gr- -> H + D + Gr
    rxn(922, {Sp::Cp, Sp::Grm}, {Sp::C, Sp::Gr}, 0.0,
        0.0),  // C+ + Gr- -> C + Gr
    rxn(923, {Sp::C2p, Sp::Grm}, {Sp::C, Sp::C, Sp::Gr}, 0.0,
        0.0),  // C2+ + Gr- -> C + C + Gr
    rxn(924, {Sp::CHp, Sp::Grm}, {Sp::H, Sp::C, Sp::Gr}, 0.0,
        0.0),  // CH+ + Gr- -> H + C + Gr
    rxn(925, {Sp::CH2p, Sp::Grm}, {Sp::H2, Sp::C, Sp::Gr}, 0.0,
        0.0),  // CH2+ + Gr- -> H2 + C + Gr
    rxn(926, {Sp::CH3p, Sp::Grm}, {Sp::H2, Sp::CH, Sp::Gr}, 0.0,
        0.0),  // CH3+ + Gr- -> H2 + CH + Gr
    rxn(927, {Sp::CH4p, Sp::Grm}, {Sp::H, Sp::CH3, Sp::Gr}, 0.0,
        0.0),  // CH4+ + Gr- -> H + CH3 + Gr
    rxn(928, {Sp::CH5p, Sp::Grm}, {Sp::H2, Sp::CH3, Sp::Gr}, 0.0,
        0.0),  // CH5+ + Gr- -> H2 + CH3 + Gr
    rxn(929, {Sp::Op, Sp::Grm}, {Sp::O, Sp::Gr}, 0.0,
        0.0),  // O+ + Gr- -> O + Gr
    rxn(930, {Sp::O2p, Sp::Grm}, {Sp::O, Sp::O, Sp::Gr}, 0.0,
        0.0),  // O2+ + Gr- -> O + O + Gr
    rxn(931, {Sp::OHp, Sp::Grm}, {Sp::H, Sp::O, Sp::Gr}, 0.0,
        0.0),  // OH+ + Gr- -> H + O + Gr
    rxn(932, {Sp::COp, Sp::Grm}, {Sp::C, Sp::O, Sp::Gr}, 0.0,
        0.0),  // CO+ + Gr- -> C + O + Gr
    rxn(933, {Sp::H2Op, Sp::Grm}, {Sp::H2, Sp::O, Sp::Gr}, 0.0,
        0.0),  // H2O+ + Gr- -> H2 + O + Gr
    rxn(934, {Sp::HCOp, Sp::Grm}, {Sp::H, Sp::CO, Sp::Gr}, 0.0,
        0.0),  // HCO+ + Gr- -> H + CO + Gr
    rxn(935, {Sp::O2Hp, Sp::Grm}, {Sp::H, Sp::O2, Sp::Gr}, 0.0,
        0.0),  // O2H+ + Gr- -> H + O2 + Gr
    rxn(936, {Sp::H3Op, Sp::Grm}, {Sp::H2, Sp::OH, Sp::Gr}, 0.0,
        0.0),  // H3O+ + Gr- -> H2 + OH + Gr
    rxn(937, {Sp::H2COp, Sp::Grm}, {Sp::H2CO, Sp::Gr}, 0.0,
        0.0),  // H2CO+ + Gr- -> H2CO + Gr
    rxn(938, {Sp::HOCOp, Sp::Grm}, {Sp::OH, Sp::CO, Sp::Gr}, 0.0,
        0.0),  // HOCO+ + Gr- -> OH + CO + Gr
    rxn(939, {Sp::H2COHp, Sp::Grm}, {Sp::H, Sp::H2CO, Sp::Gr}, 0.0,
        0.0),  // H2COH+ + Gr- -> H + H2CO + Gr
    rxn(940, {Sp::Lip, Sp::Grm}, {Sp::Li, Sp::Gr}, 0.0,
        0.0),  // Li+ + Gr- -> Li + Gr
    rxn(941, {Sp::LiHp, Sp::Grm}, {Sp::H, Sp::Li, Sp::Gr}, 0.0,
        0.0),  // LiH+ + Gr- -> H + Li + Gr
    rxn(942, {Sp::Kp, Sp::Grm}, {Sp::K, Sp::Gr}, 0.0,
        0.0),  // K+ + Gr- -> K + Gr
    rxn(943, {Sp::Nap, Sp::Grm}, {Sp::Na, Sp::Gr}, 0.0,
        0.0),  // Na+ + Gr- -> Na + Gr
    rxn(944, {Sp::Mgp, Sp::Grm}, {Sp::Mg, Sp::Gr}, 0.0,
        0.0),  // Mg+ + Gr- -> Mg + Gr
    rxn(945, {Sp::Hepp, Sp::Grm}, {Sp::He, Sp::Grp}, 0.0,
        0.0),  // He++ + Gr- -> He + Gr+
    rxn(946, {Sp::Lipp, Sp::Grm}, {Sp::Li, Sp::Grp}, 0.0,
        0.0),                                          // Li++ + Gr- -> Li + Gr+
    rxn(947, {Sp::e, Sp::Grm}, {Sp::Gr2m}, 0.0, 0.0),  // e- + Gr- -> Gr2-
    rxn(948, {Sp::Hm, Sp::Grm}, {Sp::H, Sp::Gr2m}, 0.0,
        0.0),  // H- + Gr- -> H + Gr2-
    rxn(949, {Sp::Dm, Sp::Grm}, {Sp::D, Sp::Gr2m}, 0.0,
        0.0),  // D- + Gr- -> D + Gr2-
    rxn(950, {Sp::Lim, Sp::Grm}, {Sp::Li, Sp::Gr2m}, 0.0,
        0.0),  // Li- + Gr- -> Li + Gr2-
    rxn(951, {Sp::Hp, Sp::Gr2m}, {Sp::H, Sp::Grm}, 0.0,
        0.0),  // H+ + Gr2- -> H + Gr-
    rxn(952, {Sp::H2p, Sp::Gr2m}, {Sp::H, Sp::H, Sp::Grm}, 0.0,
        0.0),  // H2+ + Gr2- -> H + H + Gr-
    rxn(953, {Sp::H3p, Sp::Gr2m}, {Sp::H, Sp::H2, Sp::Grm}, 0.0,
        0.0),  // H3+ + Gr2- -> H + H2 + Gr-
    rxn(954, {Sp::Hep, Sp::Gr2m}, {Sp::He, Sp::Grm}, 0.0,
        0.0),  // He+ + Gr2- -> He + Gr-
    rxn(955, {Sp::HeHp, Sp::Gr2m}, {Sp::H, Sp::He, Sp::Grm}, 0.0,
        0.0),  // HeH+ + Gr2- -> H + He + Gr-
    rxn(956, {Sp::Dp, Sp::Gr2m}, {Sp::D, Sp::Grm}, 0.0,
        0.0),  // D+ + Gr2- -> D + Gr-
    rxn(957, {Sp::HDp, Sp::Gr2m}, {Sp::H, Sp::D, Sp::Grm}, 0.0,
        0.0),  // HD+ + Gr2- -> H + D + Gr-
    rxn(958, {Sp::Cp, Sp::Gr2m}, {Sp::C, Sp::Grm}, 0.0,
        0.0),  // C+ + Gr2- -> C + Gr-
    rxn(959, {Sp::C2p, Sp::Gr2m}, {Sp::C, Sp::C, Sp::Grm}, 0.0,
        0.0),  // C2+ + Gr2- -> C + C + Gr-
    rxn(960, {Sp::CHp, Sp::Gr2m}, {Sp::H, Sp::C, Sp::Grm}, 0.0,
        0.0),  // CH+ + Gr2- -> H + C + Gr-
    rxn(961, {Sp::CH2p, Sp::Gr2m}, {Sp::H2, Sp::C, Sp::Grm}, 0.0,
        0.0),  // CH2+ + Gr2- -> H2 + C + Gr-
    rxn(962, {Sp::CH3p, Sp::Gr2m}, {Sp::H2, Sp::CH, Sp::Grm}, 0.0,
        0.0),  // CH3+ + Gr2- -> H2 + CH + Gr-
    rxn(963, {Sp::CH4p, Sp::Gr2m}, {Sp::H, Sp::CH3, Sp::Grm}, 0.0,
        0.0),  // CH4+ + Gr2- -> H + CH3 + Gr-
    rxn(964, {Sp::CH5p, Sp::Gr2m}, {Sp::H2, Sp::CH3, Sp::Grm}, 0.0,
        0.0),  // CH5+ + Gr2- -> H2 + CH3 + Gr-
    rxn(965, {Sp::Op, Sp::Gr2m}, {Sp::O, Sp::Grm}, 0.0,
        0.0),  // O+ + Gr2- -> O + Gr-
    rxn(966, {Sp::O2p, Sp::Gr2m}, {Sp::O, Sp::O, Sp::Grm}, 0.0,
        0.0),  // O2+ + Gr2- -> O + O + Gr-
    rxn(967, {Sp::OHp, Sp::Gr2m}, {Sp::H, Sp::O, Sp::Grm}, 0.0,
        0.0),  // OH+ + Gr2- -> H + O + Gr-
    rxn(968, {Sp::COp, Sp::Gr2m}, {Sp::C, Sp::O, Sp::Grm}, 0.0,
        0.0),  // CO+ + Gr2- -> C + O + Gr-
    rxn(969, {Sp::H2Op, Sp::Gr2m}, {Sp::H2, Sp::O, Sp::Grm}, 0.0,
        0.0),  // H2O+ + Gr2- -> H2 + O + Gr-
    rxn(970, {Sp::HCOp, Sp::Gr2m}, {Sp::H, Sp::CO, Sp::Grm}, 0.0,
        0.0),  // HCO+ + Gr2- -> H + CO + Gr-
    rxn(971, {Sp::O2Hp, Sp::Gr2m}, {Sp::H, Sp::O2, Sp::Grm}, 0.0,
        0.0),  // O2H+ + Gr2- -> H + O2 + Gr-
    rxn(972, {Sp::H3Op, Sp::Gr2m}, {Sp::H2, Sp::OH, Sp::Grm}, 0.0,
        0.0),  // H3O+ + Gr2- -> H2 + OH + Gr-
    rxn(973, {Sp::H2COp, Sp::Gr2m}, {Sp::H2CO, Sp::Grm}, 0.0,
        0.0),  // H2CO+ + Gr2- -> H2CO + Gr-
    rxn(974, {Sp::HOCOp, Sp::Gr2m}, {Sp::OH, Sp::CO, Sp::Grm}, 0.0,
        0.0),  // HOCO+ + Gr2- -> OH + CO + Gr-
    rxn(975, {Sp::H2COHp, Sp::Gr2m}, {Sp::H, Sp::H2CO, Sp::Grm}, 0.0,
        0.0),  // H2COH+ + Gr2- -> H + H2CO + Gr-
    rxn(976, {Sp::Lip, Sp::Gr2m}, {Sp::Li, Sp::Grm}, 0.0,
        0.0),  // Li+ + Gr2- -> Li + Gr-
    rxn(977, {Sp::LiHp, Sp::Gr2m}, {Sp::H, Sp::Li, Sp::Grm}, 0.0,
        0.0),  // LiH+ + Gr2- -> H + Li + Gr-
    rxn(978, {Sp::Kp, Sp::Gr2m}, {Sp::K, Sp::Grm}, 0.0,
        0.0),  // K+ + Gr2- -> K + Gr-
    rxn(979, {Sp::Nap, Sp::Gr2m}, {Sp::Na, Sp::Grm}, 0.0,
        0.0),  // Na+ + Gr2- -> Na + Gr-
    rxn(980, {Sp::Mgp, Sp::Gr2m}, {Sp::Mg, Sp::Grm}, 0.0,
        0.0),  // Mg+ + Gr2- -> Mg + Gr-
    rxn(981, {Sp::Hepp, Sp::Gr2m}, {Sp::He, Sp::Gr}, 0.0,
        0.0),  // He++ + Gr2- -> He + Gr
    rxn(982, {Sp::Lipp, Sp::Gr2m}, {Sp::Li, Sp::Gr}, 0.0,
        0.0),                                          // Li++ + Gr2- -> Li + Gr
    rxn(983, {Sp::Gr2m}, {Sp::e, Sp::Grm}, 0.0, 0.0),  // Gr2- -> e- + Gr-
    rxn(984, {Sp::Grm}, {Sp::e, Sp::Gr}, 0.0, 0.0),    // Gr- -> e- + Gr
    rxn(985, {Sp::Gr}, {Sp::e, Sp::Grp}, 0.0, 0.0),    // Gr -> e- + Gr+
    rxn(986, {Sp::Grp}, {Sp::e, Sp::Gr2p}, 0.0, 0.0),  // Gr+ -> e- + Gr2+
    rxn(987, {Sp::Grp, Sp::Grm}, {Sp::Gr, Sp::Gr}, 0.0,
        0.0),  // Gr+ + Gr- -> Gr + Gr
    rxn(988, {Sp::Grp, Sp::Gr2m}, {Sp::Gr, Sp::Grm}, 0.0,
        0.0),  // Gr+ + Gr2- -> Gr + Gr-
    rxn(989, {Sp::Gr2p, Sp::Gr2m}, {Sp::Gr, Sp::Gr}, 0.0,
        0.0),  // Gr2+ + Gr2- -> Gr + Gr
    rxn(990, {Sp::Gr2p, Sp::Grm}, {Sp::Gr, Sp::Grp}, 0.0,
        0.0),  // Gr2+ + Gr- -> Gr + Gr+
}};

// High-density Saha-equilibrium reactions (53 rows).
inline constexpr std::array<SahaReaction, 53> kSaha = {{
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
    saha_rxn(51, {Sp::e, Sp::Dp}, {Sp::D, PHOTON}, 2.78909e-41,
             2.1793e-11),  // e- + D+ -> D + γ
    saha_rxn(54, {Sp::H, Sp::D}, {Sp::HD, PHOTON}, 1.19083e-36,
             7.23181e-12),  // H + D -> HD + γ
    saha_rxn(60, {Sp::Hp, Sp::D}, {Sp::HDp, PHOTON}, 1.19018e-36,
             4.27406e-12),  // H+ + D -> HD+ + γ
    saha_rxn(63, {Sp::e, Sp::D}, {Sp::Dm, PHOTON}, 2.78909e-41,
             1.20909e-12),  // e- + D -> D- + γ
    saha_rxn(185, {Sp::H, Sp::C}, {Sp::CH, PHOTON}, 1.93928e-36,
             5.55836e-12),  // H + C -> CH + γ
    saha_rxn(186, {Sp::C, Sp::C}, {Sp::C2, PHOTON}, 3.17924e-35,
             9.99787e-12),  // C + C -> C2 + γ
    saha_rxn(187, {Sp::C, Sp::O}, {Sp::CO, PHOTON}, 3.88349e-35,
             1.77977e-11),  // C + O -> CO + γ
    saha_rxn(193, {Sp::H, Sp::O}, {Sp::OH, PHOTON}, 1.99697e-36,
             7.06931e-12),  // H + O -> OH + γ
    saha_rxn(194, {Sp::O, Sp::O}, {Sp::O2, PHOTON}, 4.89242e-35,
             8.19609e-12),  // O + O -> O2 + γ
    saha_rxn(295, {Sp::H, Sp::Cp}, {Sp::CHp, PHOTON}, 1.93927e-36,
             6.64079e-12),  // H + C+ -> CH+ + γ
    saha_rxn(298, {Sp::H2, Sp::Cp}, {Sp::CH2p, PHOTON}, 4.90418e-36,
             6.82349e-12),  // H2 + C+ -> CH2+ + γ
    saha_rxn(340, {Sp::H2, Sp::CH3p}, {Sp::CH5p, PHOTON}, 5.12508e-36,
             2.93729e-12),  // H2 + CH3+ -> CH5+ + γ
    saha_rxn(507, {Sp::e, Sp::Cp}, {Sp::C, PHOTON}, 2.79004e-41,
             1.8147e-11),  // e- + C+ -> C + γ
    saha_rxn(514, {Sp::e, Sp::CH3p}, {Sp::CH3, PHOTON}, 2.79008e-41,
             1.57647e-11),  // e- + CH3+ -> CH3 + γ
    saha_rxn(515, {Sp::e, Sp::Op}, {Sp::O, PHOTON}, 2.79009e-41,
             2.18189e-11),  // e- + O+ -> O + γ
    saha_rxn(530, {Sp::e, Sp::H2COp}, {Sp::H2CO, PHOTON}, 2.79015e-41,
             1.74545e-11),  // e- + H2CO+ -> H2CO + γ
    saha_rxn(538, {Sp::H2, Sp::C}, {Sp::CH2, PHOTON}, 4.90423e-36,
             5.32199e-12),  // H2 + C -> CH2 + γ
    saha_rxn(641, {Sp::Hp, Sp::He}, {Sp::HeHp, PHOTON}, 1.56169e-36,
             2.90085e-12),  // H+ + He -> HeH+ + γ
    saha_rxn(642, {Sp::H, Sp::OH}, {Sp::H2O, PHOTON}, 2.00754e-36,
             8.17078e-12),  // H + OH -> H2O + γ
    saha_rxn(643, {Sp::H2, Sp::CH}, {Sp::CH3, PHOTON}, 4.98735e-36,
             7.34789e-12),  // H2 + CH -> CH3 + γ
    saha_rxn(644, {Sp::C, Sp::Cp}, {Sp::C2p, PHOTON}, 3.17913e-35,
             9.08202e-12),  // C + C+ -> C2+ + γ
    saha_rxn(645, {Sp::C, Sp::Op}, {Sp::COp, PHOTON}, 3.8834e-35,
             1.71635e-11),  // C + O+ -> CO+ + γ
    saha_rxn(701, {Sp::e, Sp::Kp}, {Sp::K, PHOTON}, 2.79017e-41,
             6.95434e-12),  // e- + K+ -> K + γ
    saha_rxn(703, {Sp::e, Sp::Nap}, {Sp::Na, PHOTON}, 2.79013e-41,
             8.23461e-12),  // e- + Na+ -> Na + γ
    saha_rxn(718, {Sp::e, Sp::Mgp}, {Sp::Mg, PHOTON}, 2.79014e-41,
             1.22506e-11),  // e- + Mg+ -> Mg + γ
    saha_rxn(801, {Sp::e, Sp::Lip}, {Sp::Li, PHOTON}, 2.78963e-41,
             8.63849e-12),  // e- + Li+ -> Li + γ
    saha_rxn(804, {Sp::e, Sp::Li}, {Sp::Lim, PHOTON}, 2.78963e-41,
             9.90146e-13),  // e- + Li -> Li- + γ
    saha_rxn(813, {Sp::Hp, Sp::Li}, {Sp::LiHp, PHOTON}, 1.78437e-36,
             1.33808e-11),  // H+ + Li -> LiH+ + γ
    saha_rxn(818, {Sp::H, Sp::Li}, {Sp::LiH, PHOTON}, 1.78583e-36,
             3.86847e-12),  // H + Li -> LiH + γ
    saha_rxn(821, {Sp::e, Sp::Lipp}, {Sp::Lip, PHOTON}, 2.79023e-41,
             1.21191e-10),  // e- + Li++ -> Li+ + γ
    saha_rxn(822, {Sp::e, Sp::Lippp}, {Sp::Lipp, PHOTON}, 2.79023e-41,
             1.96196e-10),  // e- + Li+++ -> Li++ + γ
    saha_rxn(8, {Sp::H, Sp::Hm}, {Sp::H2, Sp::e}, 27755.2,
             5.96598e-12),  // H + H- -> H2 + e-
    saha_rxn(26, {Sp::H2, Sp::H2p}, {Sp::H, Sp::H3p}, 1.53938,
             2.69052e-12),  // H2 + H2+ -> H + H3+
    saha_rxn(112, {Sp::H, Sp::H2CO}, {Sp::H2, Sp::HCO}, 0.372144,
             1.14219e-12),  // H + H2CO -> H2 + HCO
    saha_rxn(115, {Sp::H, Sp::O2H}, {Sp::O, Sp::H2O}, 0.0392231,
             3.70921e-12),  // H + O2H -> O + H2O
    saha_rxn(136, {Sp::H2, Sp::CH3}, {Sp::H, Sp::CH4}, 2.56593,
             9.94663e-15),  // H2 + CH3 -> H + CH4
    saha_rxn(184, {Sp::H, Sp::HCO}, {Sp::H2, Sp::CO}, 0.372819,
             6.17347e-12),  // H + HCO -> H2 + CO
    saha_rxn(201, {Sp::O, Sp::HCO}, {Sp::H, Sp::CO2}, 33.8471,
             7.73522e-12),  // O + HCO -> H + CO2
    saha_rxn(213, {Sp::Hp, Sp::CH4}, {Sp::H, Sp::CH4p}, 0.999233,
             1.57016e-12),  // H+ + CH4 -> H + CH4+
    saha_rxn(518, {Sp::e, Sp::OHp}, {Sp::H, Sp::O}, 1.39717e-05,
             1.37852e-11),  // e- + OH+ -> H + O
    saha_rxn(521, {Sp::e, Sp::H2Op}, {Sp::H, Sp::OH}, 1.38981e-05,
             1.20448e-11),  // e- + H2O+ -> H + OH
    saha_rxn(523, {Sp::e, Sp::H3Op}, {Sp::H, Sp::H2O}, 1.38329e-05,
             1.04133e-11),  // e- + H3O+ -> H + H2O
    saha_rxn(527, {Sp::e, Sp::HCOp}, {Sp::H, Sp::CO}, 1.34428e-05,
             1.20419e-11),  // e- + HCO+ -> H + CO
    saha_rxn(533, {Sp::e, Sp::H2COHp}, {Sp::H, Sp::H2CO}, 1.33958e-05,
             1.00843e-11),  // e- + H2COH+ -> H + H2CO
    saha_rxn(534, {Sp::e, Sp::O2p}, {Sp::O, Sp::O}, 5.70303e-07,
             1.11492e-11),  // e- + O2+ -> O + O
    saha_rxn(535, {Sp::e, Sp::O2Hp}, {Sp::H, Sp::O2}, 1.33555e-05,
             1.4854e-11),  // e- + O2H+ -> H + O2
    saha_rxn(536, {Sp::e, Sp::HOCOp}, {Sp::H, Sp::CO2}, 1.3189e-05,
             1.29192e-11),  // e- + HOCO+ -> H + CO2
    saha_rxn(602, {Sp::H, Sp::H2O2}, {Sp::OH, Sp::H2O}, 0.0374399,
             4.78179e-12),  // H + H2O2 -> OH + H2O
}};

// Grain-surface reactions (148 rows).
inline constexpr std::array<GrainReaction, 148> kGrain = {{
    grain_rxn(1, {Sp::H}, {Sp::H_p}, 1),                 // H -> H_p
    grain_rxn(2, {Sp::H}, {Sp::H_c}, 1),                 // H -> H_c
    grain_rxn(3, {Sp::H2}, {Sp::H2_p}, 1),               // H2 -> H2_p
    grain_rxn(4, {Sp::D}, {Sp::D_p}, 1),                 // D -> D_p
    grain_rxn(5, {Sp::D}, {Sp::D_c}, 1),                 // D -> D_c
    grain_rxn(6, {Sp::HD}, {Sp::HD_p}, 1),               // HD -> HD_p
    grain_rxn(7, {Sp::O}, {Sp::O_p}, 1),                 // O -> O_p
    grain_rxn(8, {Sp::O2}, {Sp::O2_p}, 1),               // O2 -> O2_p
    grain_rxn(9, {Sp::OH}, {Sp::OH_p}, 1),               // OH -> OH_p
    grain_rxn(10, {Sp::CO}, {Sp::CO_p}, 1),              // CO -> CO_p
    grain_rxn(11, {Sp::CO2}, {Sp::CO2_p}, 1),            // CO2 -> CO2_p
    grain_rxn(12, {Sp::H2O}, {Sp::H2O_p}, 1),            // H2O -> H2O_p
    grain_rxn(13, {Sp::O2H}, {Sp::HO2_p}, 1),            // O2H -> HO2_p
    grain_rxn(14, {Sp::H2O2}, {Sp::H2O2_p}, 1),          // H2O2 -> H2O2_p
    grain_rxn(15, {Sp::HCO}, {Sp::HCO_p}, 1),            // HCO -> HCO_p
    grain_rxn(16, {Sp::H2CO}, {Sp::H2CO_p}, 1),          // H2CO -> H2CO_p
    grain_rxn(17, {Sp::C}, {Sp::C_p}, 1),                // C -> C_p
    grain_rxn(18, {Sp::CH}, {Sp::CH_p}, 1),              // CH -> CH_p
    grain_rxn(19, {Sp::CH2}, {Sp::CH2_p}, 1),            // CH2 -> CH2_p
    grain_rxn(20, {Sp::CH3}, {Sp::CH3_p}, 1),            // CH3 -> CH3_p
    grain_rxn(21, {Sp::CH4}, {Sp::CH4_p}, 1),            // CH4 -> CH4_p
    grain_rxn(22, {Sp::H_p}, {Sp::H}, 0),                // H_p -> H
    grain_rxn(23, {Sp::H_c}, {Sp::H}, 0),                // H_c -> H
    grain_rxn(24, {Sp::H2_p}, {Sp::H2}, 0),              // H2_p -> H2
    grain_rxn(25, {Sp::D_p}, {Sp::D}, 0),                // D_p -> D
    grain_rxn(26, {Sp::D_c}, {Sp::D}, 0),                // D_c -> D
    grain_rxn(27, {Sp::HD_p}, {Sp::HD}, 0),              // HD_p -> HD
    grain_rxn(28, {Sp::O_p}, {Sp::O}, 0),                // O_p -> O
    grain_rxn(29, {Sp::O2_p}, {Sp::O2}, 0),              // O2_p -> O2
    grain_rxn(30, {Sp::OH_p}, {Sp::OH}, 0),              // OH_p -> OH
    grain_rxn(31, {Sp::CO_p}, {Sp::CO}, 0),              // CO_p -> CO
    grain_rxn(32, {Sp::CO2_p}, {Sp::CO2}, 0),            // CO2_p -> CO2
    grain_rxn(33, {Sp::H2O_p}, {Sp::H2O}, 0),            // H2O_p -> H2O
    grain_rxn(34, {Sp::HO2_p}, {Sp::O, Sp::OH}, 0),      // HO2_p -> O + OH
    grain_rxn(35, {Sp::H2O2_p}, {Sp::H2O2}, 0),          // H2O2_p -> H2O2
    grain_rxn(36, {Sp::HCO_p}, {Sp::HCO}, 0),            // HCO_p -> HCO
    grain_rxn(37, {Sp::H2CO_p}, {Sp::H2CO}, 0),          // H2CO_p -> H2CO
    grain_rxn(38, {Sp::C_p}, {Sp::C}, 0),                // C_p -> C
    grain_rxn(39, {Sp::CH_p}, {Sp::CH}, 0),              // CH_p -> CH
    grain_rxn(40, {Sp::CH2_p}, {Sp::CH2}, 0),            // CH2_p -> CH2
    grain_rxn(41, {Sp::CH3_p}, {Sp::CH3}, 0),            // CH3_p -> CH3
    grain_rxn(42, {Sp::CH4_p}, {Sp::CH4}, 0),            // CH4_p -> CH4
    grain_rxn(43, {Sp::H_p, Sp::H_p}, {Sp::H2_p}, 1),    // H_p + H_p -> H2_p
    grain_rxn(44, {Sp::H_p, Sp::D_p}, {Sp::HD_p}, 1),    // H_p + D_p -> HD_p
    grain_rxn(45, {Sp::H_p, Sp::O_p}, {Sp::OH_p}, 1),    // H_p + O_p -> OH_p
    grain_rxn(46, {Sp::H_p, Sp::OH_p}, {Sp::H2O_p}, 1),  // H_p + OH_p -> H2O_p
    grain_rxn(47, {Sp::H_p, Sp::O2_p}, {Sp::HO2_p}, 1),  // H_p + O2_p -> HO2_p
    grain_rxn(48, {Sp::H_p, Sp::CO_p}, {Sp::HCO_p}, 1),  // H_p + CO_p -> HCO_p
    grain_rxn(49, {Sp::H_p, Sp::HO2_p}, {Sp::H2O2_p},
              1),  // H_p + HO2_p -> H2O2_p
    grain_rxn(50, {Sp::H_p, Sp::HCO_p}, {Sp::H2CO_p},
              1),                                      // H_p + HCO_p -> H2CO_p
    grain_rxn(51, {Sp::H_p, Sp::C_p}, {Sp::CH_p}, 1),  // H_p + C_p -> CH_p
    grain_rxn(52, {Sp::H_p, Sp::CH_p}, {Sp::CH2_p}, 1),  // H_p + CH_p -> CH2_p
    grain_rxn(53, {Sp::H_p, Sp::CH2_p}, {Sp::CH3_p},
              1),  // H_p + CH2_p -> CH3_p
    grain_rxn(54, {Sp::H_p, Sp::CH3_p}, {Sp::CH4_p},
              1),                                        // H_p + CH3_p -> CH4_p
    grain_rxn(55, {Sp::O_p, Sp::O_p}, {Sp::O2_p}, 1),    // O_p + O_p -> O2_p
    grain_rxn(56, {Sp::O_p, Sp::C_p}, {Sp::CO_p}, 1),    // O_p + C_p -> CO_p
    grain_rxn(57, {Sp::O_p, Sp::CO_p}, {Sp::CO2_p}, 1),  // O_p + CO_p -> CO2_p
    grain_rxn(58, {Sp::OH_p, Sp::OH_p}, {Sp::H2O2_p},
              1),                                      // OH_p + OH_p -> H2O2_p
    grain_rxn(59, {Sp::H_p, Sp::H_p}, {Sp::H2}, 1),    // H_p + H_p -> H2
    grain_rxn(60, {Sp::H_p, Sp::D_p}, {Sp::HD}, 1),    // H_p + D_p -> HD
    grain_rxn(61, {Sp::H_p, Sp::O_p}, {Sp::OH}, 1),    // H_p + O_p -> OH
    grain_rxn(62, {Sp::H_p, Sp::OH_p}, {Sp::H2O}, 1),  // H_p + OH_p -> H2O
    grain_rxn(63, {Sp::H_p, Sp::O2_p}, {Sp::O2H}, 1),  // H_p + O2_p -> O2H
    grain_rxn(64, {Sp::H_p, Sp::CO_p}, {Sp::HCO}, 1),  // H_p + CO_p -> HCO
    grain_rxn(65, {Sp::H_p, Sp::HO2_p}, {Sp::H2O2}, 1),  // H_p + HO2_p -> H2O2
    grain_rxn(66, {Sp::H_p, Sp::HCO_p}, {Sp::H2CO}, 1),  // H_p + HCO_p -> H2CO
    grain_rxn(67, {Sp::H_p, Sp::C_p}, {Sp::CH}, 1),      // H_p + C_p -> CH
    grain_rxn(68, {Sp::H_p, Sp::CH_p}, {Sp::CH2}, 1),    // H_p + CH_p -> CH2
    grain_rxn(69, {Sp::H_p, Sp::CH2_p}, {Sp::CH3}, 1),   // H_p + CH2_p -> CH3
    grain_rxn(70, {Sp::H_p, Sp::CH3_p}, {Sp::CH4}, 1),   // H_p + CH3_p -> CH4
    grain_rxn(71, {Sp::O_p, Sp::O_p}, {Sp::O2}, 1),      // O_p + O_p -> O2
    grain_rxn(72, {Sp::O_p, Sp::C_p}, {Sp::CO}, 1),      // O_p + C_p -> CO
    grain_rxn(73, {Sp::O_p, Sp::CO_p}, {Sp::CO2}, 1),    // O_p + CO_p -> CO2
    grain_rxn(74, {Sp::OH_p, Sp::OH_p}, {Sp::H2O2}, 1),  // OH_p + OH_p -> H2O2
    grain_rxn(75, {Sp::H_p, Sp::H2O_p}, {Sp::H2_p, Sp::OH_p},
              1),  // H_p + H2O_p -> H2_p + OH_p
    grain_rxn(76, {Sp::H_p, Sp::HO2_p}, {Sp::OH_p, Sp::OH_p},
              1),  // H_p + HO2_p -> OH_p + OH_p
    grain_rxn(77, {Sp::H_p, Sp::HO2_p}, {Sp::OH, Sp::OH},
              1),  // H_p + HO2_p -> OH + OH
    grain_rxn(78, {Sp::H_p, Sp::H2O2_p}, {Sp::OH_p, Sp::H2O_p},
              1),  // H_p + H2O2_p -> OH_p + H2O_p
    grain_rxn(79, {Sp::H_p, Sp::H2O2_p}, {Sp::OH, Sp::H2O_p},
              1),  // H_p + H2O2_p -> OH + H2O_p
    grain_rxn(80, {Sp::H_p, Sp::H2O2_p}, {Sp::OH, Sp::H2O},
              1),  // H_p + H2O2_p -> OH + H2O
    grain_rxn(81, {Sp::H_p, Sp::HCO_p}, {Sp::H2_p, Sp::CO_p},
              1),  // H_p + HCO_p -> H2_p + CO_p
    grain_rxn(82, {Sp::H_p, Sp::HCO_p}, {Sp::H2, Sp::CO_p},
              1),  // H_p + HCO_p -> H2 + CO_p
    grain_rxn(83, {Sp::H_p, Sp::HCO_p}, {Sp::H2, Sp::CO},
              1),  // H_p + HCO_p -> H2 + CO
    grain_rxn(84, {Sp::H_p, Sp::H2CO_p}, {Sp::H2_p, Sp::HCO_p},
              1),  // H_p + H2CO_p -> H2_p + HCO_p
    grain_rxn(85, {Sp::H_p, Sp::H2CO_p}, {Sp::H2, Sp::HCO_p},
              1),  // H_p + H2CO_p -> H2 + HCO_p
    grain_rxn(86, {Sp::H_p, Sp::H2CO_p}, {Sp::H2, Sp::HCO},
              1),  // H_p + H2CO_p -> H2 + HCO
    grain_rxn(87, {Sp::H_p, Sp::CO2_p}, {Sp::OH_p, Sp::CO_p},
              1),  // H_p + CO2_p -> OH_p + CO_p
    grain_rxn(88, {Sp::H_p, Sp::CH_p}, {Sp::H2_p, Sp::C_p},
              1),  // H_p + CH_p -> H2_p + C_p
    grain_rxn(89, {Sp::H_p, Sp::CH_p}, {Sp::H2, Sp::C_p},
              1),  // H_p + CH_p -> H2 + C_p
    grain_rxn(90, {Sp::H_p, Sp::CH_p}, {Sp::H2, Sp::C},
              1),  // H_p + CH_p -> H2 + C
    grain_rxn(91, {Sp::H_p, Sp::CH2_p}, {Sp::H2_p, Sp::CH_p},
              1),  // H_p + CH2_p -> H2_p + CH_p
    grain_rxn(92, {Sp::H_p, Sp::CH2_p}, {Sp::H2, Sp::CH_p},
              1),  // H_p + CH2_p -> H2 + CH_p
    grain_rxn(93, {Sp::H_p, Sp::CH2_p}, {Sp::H2, Sp::CH},
              1),  // H_p + CH2_p -> H2 + CH
    grain_rxn(94, {Sp::H_p, Sp::CH3_p}, {Sp::H2_p, Sp::CH2_p},
              1),  // H_p + CH3_p -> H2_p + CH2_p
    grain_rxn(95, {Sp::H_p, Sp::CH4_p}, {Sp::H2_p, Sp::CH3_p},
              1),  // H_p + CH4_p -> H2_p + CH3_p
    grain_rxn(96, {Sp::O_p, Sp::OH_p}, {Sp::H_p, Sp::O2_p},
              1),  // O_p + OH_p -> H_p + O2_p
    grain_rxn(97, {Sp::O_p, Sp::OH_p}, {Sp::H, Sp::O2_p},
              1),  // O_p + OH_p -> H + O2_p
    grain_rxn(98, {Sp::O_p, Sp::OH_p}, {Sp::H, Sp::O2},
              1),  // O_p + OH_p -> H + O2
    grain_rxn(99, {Sp::O_p, Sp::HO2_p}, {Sp::O2_p, Sp::OH_p},
              1),  // O_p + HO2_p -> O2_p + OH_p
    grain_rxn(100, {Sp::O_p, Sp::HO2_p}, {Sp::O2, Sp::OH_p},
              1),  // O_p + HO2_p -> O2 + OH_p
    grain_rxn(101, {Sp::O_p, Sp::HO2_p}, {Sp::O2, Sp::OH},
              1),  // O_p + HO2_p -> O2 + OH
    grain_rxn(102, {Sp::O_p, Sp::HCO_p}, {Sp::H_p, Sp::CO2_p},
              1),  // O_p + HCO_p -> H_p + CO2_p
    grain_rxn(103, {Sp::O_p, Sp::HCO_p}, {Sp::H, Sp::CO2_p},
              1),  // O_p + HCO_p -> H + CO2_p
    grain_rxn(104, {Sp::O_p, Sp::HCO_p}, {Sp::H, Sp::CO2},
              1),  // O_p + HCO_p -> H + CO2
    grain_rxn(105, {Sp::O_p, Sp::H2CO_p}, {Sp::H2_p, Sp::CO2_p},
              1),  // O_p + H2CO_p -> H2_p + CO2_p
    grain_rxn(106, {Sp::O_p, Sp::H2CO_p}, {Sp::H2, Sp::CO2_p},
              1),  // O_p + H2CO_p -> H2 + CO2_p
    grain_rxn(107, {Sp::O_p, Sp::H2CO_p}, {Sp::H2, Sp::CO2},
              1),  // O_p + H2CO_p -> H2 + CO2
    grain_rxn(108, {Sp::H2_p, Sp::OH_p}, {Sp::H_p, Sp::H2O_p},
              1),  // H2_p + OH_p -> H_p + H2O_p
    grain_rxn(109, {Sp::H2_p, Sp::OH_p}, {Sp::H, Sp::H2O_p},
              1),  // H2_p + OH_p -> H + H2O_p
    grain_rxn(110, {Sp::OH_p, Sp::CO_p}, {Sp::H_p, Sp::CO2_p},
              1),  // OH_p + CO_p -> H_p + CO2_p
    grain_rxn(111, {Sp::OH_p, Sp::CO_p}, {Sp::H, Sp::CO2_p},
              1),  // OH_p + CO_p -> H + CO2_p
    grain_rxn(112, {Sp::OH_p, Sp::CO_p}, {Sp::H, Sp::CO2},
              1),  // OH_p + CO_p -> H + CO2
    grain_rxn(113, {Sp::OH_p, Sp::HCO_p}, {Sp::H2_p, Sp::CO2_p},
              1),  // OH_p + HCO_p -> H2_p + CO2_p
    grain_rxn(114, {Sp::OH_p, Sp::HCO_p}, {Sp::H2, Sp::CO2_p},
              1),  // OH_p + HCO_p -> H2 + CO2_p
    grain_rxn(115, {Sp::OH_p, Sp::HCO_p}, {Sp::H2, Sp::CO2},
              1),  // OH_p + HCO_p -> H2 + CO2
    grain_rxn(116, {Sp::H2_p, Sp::HO2_p}, {Sp::H_p, Sp::H2O2_p},
              1),                             // H2_p + HO2_p -> H_p + H2O2_p
    grain_rxn(117, {Sp::H_p}, {Sp::H_c}, 0),  // H_p -> H_c
    grain_rxn(118, {Sp::H, Sp::H_c}, {Sp::H2}, 2),    // H + H_c -> H2
    grain_rxn(119, {Sp::H_p, Sp::H_c}, {Sp::H2}, 1),  // H_p + H_c -> H2
    grain_rxn(120, {Sp::D_p}, {Sp::D_c}, 0),          // D_p -> D_c
    grain_rxn(121, {Sp::H, Sp::D_c}, {Sp::HD}, 2),    // H + D_c -> HD
    grain_rxn(122, {Sp::D, Sp::H_c}, {Sp::HD}, 2),    // D + H_c -> HD
    grain_rxn(123, {Sp::H_p, Sp::D_c}, {Sp::HD}, 1),  // H_p + D_c -> HD
    grain_rxn(124, {Sp::H_c, Sp::D_p}, {Sp::HD}, 1),  // H_c + D_p -> HD
    grain_rxn(125, {Sp::H, Sp::H_p}, {Sp::H2}, 2),    // H + H_p -> H2
    grain_rxn(126, {Sp::H, Sp::D_p}, {Sp::HD}, 2),    // H + D_p -> HD
    grain_rxn(127, {Sp::D, Sp::H_p}, {Sp::HD}, 2),    // D + H_p -> HD
    grain_rxn(131, {Sp::H_p}, {Sp::H}, 0),            // H_p -> H
    grain_rxn(132, {Sp::H_c}, {Sp::H}, 0),            // H_c -> H
    grain_rxn(133, {Sp::H2_p}, {Sp::H2}, 0),          // H2_p -> H2
    grain_rxn(134, {Sp::D_p}, {Sp::D}, 0),            // D_p -> D
    grain_rxn(135, {Sp::D_c}, {Sp::D}, 0),            // D_c -> D
    grain_rxn(136, {Sp::HD_p}, {Sp::HD}, 0),          // HD_p -> HD
    grain_rxn(137, {Sp::O_p}, {Sp::O}, 0),            // O_p -> O
    grain_rxn(138, {Sp::O2_p}, {Sp::O2}, 0),          // O2_p -> O2
    grain_rxn(139, {Sp::OH_p}, {Sp::OH}, 0),          // OH_p -> OH
    grain_rxn(140, {Sp::CO_p}, {Sp::CO}, 0),          // CO_p -> CO
    grain_rxn(141, {Sp::CO2_p}, {Sp::CO2}, 0),        // CO2_p -> CO2
    grain_rxn(142, {Sp::H2O_p}, {Sp::H2O}, 0),        // H2O_p -> H2O
    grain_rxn(143, {Sp::HO2_p}, {Sp::O, Sp::OH}, 0),  // HO2_p -> O + OH
    grain_rxn(144, {Sp::H2O2_p}, {Sp::H2O2}, 0),      // H2O2_p -> H2O2
    grain_rxn(145, {Sp::HCO_p}, {Sp::HCO}, 0),        // HCO_p -> HCO
    grain_rxn(146, {Sp::H2CO_p}, {Sp::H2CO}, 0),      // H2CO_p -> H2CO
    grain_rxn(147, {Sp::C_p}, {Sp::C}, 0),            // C_p -> C
    grain_rxn(148, {Sp::CH_p}, {Sp::CH}, 0),          // CH_p -> CH
    grain_rxn(149, {Sp::CH2_p}, {Sp::CH2}, 0),        // CH2_p -> CH2
    grain_rxn(150, {Sp::CH3_p}, {Sp::CH3}, 0),        // CH3_p -> CH3
    grain_rxn(151, {Sp::CH4_p}, {Sp::CH4}, 0),        // CH4_p -> CH4
}};

// ───────────────────────────────────────────────────────────────────────────
// Minimal-network keep-sets.
//
// Reaction-id selections that carve the compact metal_grain Minimal network
// (40 species) out of the full network above.  Each array lists the full
// reaction `num`s the Minimal model retains; models/metal_grain/minimal.h
// looks each id up here, remaps its species through the canonical<->local
// maps, and renumbers 1..N to build the compact topology.  Grouped by the
// compact solver block they land in so the build order matches the compact
// loop partition (standard reactions carry detailed-balance reverse; CR,
// grain-charge and grain-surface entries are first-order).
//
// C + O2 -> O + CO (Nakauchi R54, neutral) is carried by id 191; id 310 is the
// distinct C+ + O2 -> CO + O+ ion-neutral reaction (not in Nakauchi's reduced
// network) and is not kept.  O + OH -> H + O2 (R55) is carried by id 617 (the
// active rate); the duplicate id 199 holds an identically-zero coefficient and
// is omitted.  Ids 17 (e + H2+ -> 2H), 44 (H + H2+ -> H3+ + γ) and 233
// (H2+ + O -> H + OH+) are H2+ channels not in Nakauchi's table but kept for
// fidelity to the full network (H2+ is a live species here).  2H + grain -> H2
// (id 23) and H + D + grain -> HD (id 144) are grain-catalysed three-body
// reactions in this table, not surface reactions, so they sit in the standard
// block.  Grain charge transfer has no detailed-balance reverse, so its two
// reversible reactions keep both stored directions (ids 873/984, 907/985,
// 947/983).
// ───────────────────────────────────────────────────────────────────────────

// Standard reactions (detailed-balance reverse), 68.  Reaction 196
// (CH + O -> HCO+ + e-) is a chemi-ionization channel present in the arche full
// network but absent from Nakauchi's reduced network (neither the woIon nor the
// wiIon Table-1 keep-set carries it).  In the compact ζ=0 deep-recombination
// regime the reduced network retains an abundant CH/O pool (carbon is not
// drained to CO2 ice as in the full network), so reaction 196 acts as a
// spurious electron
// + HCO+ source: it props the gas charge, which the grains balance by symmetric
// charging (Gr+ ≈ Gr-), producing the gr- spike in fig10 (d1).  It is dropped
// to match Nakauchi's reduced network; the full network keeps it (CH ≈ 0 there,
// so it is inert).  Reaction 195 (CH + O -> H + CO) is its neutral sibling:
// also absent from Nakauchi's reduced network, and dropped for the same
// fidelity even though it carries no charge.  The slot it vacated holds
// reaction 54 (H + D -> HD + γ, radiative HD formation) — an off-table arche
// channel retained for full-network fidelity alongside 17, 44 and 233.  Keeping
// 54 in 195's slot leaves the grain-catalysed special-reaction positions
// (2H+grain->H2 at compact num 10, 2D+grain->HD at num 26) and the cosmic-ray
// block undisturbed.
inline constexpr std::array<int, 68> kMetalMinimalKeep = {{
    2,   4,   7,   8,   9,   10,  17,  19,  21,  23,  24,  26,  28,  30,
    34,  35,  44,  50,  52,  55,  57,  101, 103, 106, 137, 144, 184, 185,
    189, 191, 193, 54,  214, 215, 233, 256, 261, 262, 263, 279, 283, 346,
    348, 377, 410, 426, 464, 507, 518, 521, 522, 523, 524, 527, 538, 617,
    621, 636, 638, 642, 718, 719, 721, 722, 729, 801, 806, 830,
}};

// Cosmic-ray channels (first-order), 8.  The four direct-CR ionizations lead
// (they feed the CR-heating sum) and the four CR-induced-photon channels
// follow.
inline constexpr std::array<int, 8> kMetalMinimalCRKeep = {{
    544,
    545,
    548,
    549,  // direct CR: H, He, H2->H+e+H+, H2->e+H2+
    656,
    672,
    677,
    678,  // CRph: C, O2, H, He
}};

// Ion-grain charge transfer (first-order; reversible R67-R69 keep both
// rows), 27.
inline constexpr std::array<int, 27> kMetalMinimalChargeKeep = {{
    841, 843, 848, 862, 866, 870, 873, 907, 915, 917, 918, 922, 936, 940,
    944, 947, 951, 953, 958, 972, 976, 980, 983, 984, 985, 987, 988,
}};

// Grain-surface freeze-out: adsorption (5) + thermal desorption (5), 10.
inline constexpr std::array<int, 10> kMetalMinimalGrainKeep = {{
    7,
    9,
    10,
    12,
    17,  // adsorption  O, OH, CO, H2O, C
    28,
    30,
    31,
    33,
    38,  // thermal desorption  O_p, OH_p, CO_p, H2O_p, C_p
}};

// High-density Saha sub-table: equilibrium constants among the retained
// species, 23.  Provisional — the compact 4D Newton (equilibrium.h) selects
// the subset it references.
inline constexpr std::array<int, 23> kMetalMinimalSahaKeep = {{
    2,   4,   7,   8,   9,   26,  51,  54,  184, 185, 187, 193,
    194, 507, 515, 518, 521, 523, 527, 538, 642, 718, 801,
}};

// Fill the topology fields of a ReactionTable from the constexpr network
// above.  Partition-function tables are loaded separately from HDF5.
inline void init_topology(MetalGrainTable& tbl) {
  for (std::size_t i = 0; i < kReactions.size(); ++i)
    tbl.reactions[i] = kReactions[i];
  tbl.n_loaded = static_cast<int>(kReactions.size());
  for (std::size_t i = 0; i < kSaha.size(); ++i) tbl.saha[i] = kSaha[i];
  tbl.n_saha = static_cast<int>(kSaha.size());
  for (std::size_t i = 0; i < kGrain.size(); ++i)
    tbl.grain_reactions[i] = kGrain[i];
  tbl.n_grain = static_cast<int>(kGrain.size());
}

}  // namespace net
}  // namespace metal_grain
}  // namespace arche
