// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// reaction_metal.h — metal_grain-specific reaction rate block functions.
//
// Provides:
//   compute_metal_rates(...)    — UMIST metal reactions k_rxn[100..644].
//   compute_Li_rates(...)       — Li reactions k_rxn[800..829].
//   compute_KNaMg_rates(...)    — K/Na/Mg reactions k_rxn[700..729].
//   compute_CR_rates_metal(...) — CR + CR-photo reactions for metal_grain.
//   compute_grain_rates(...)    — grain charge + surface reactions.
//
// Shared block functions (H/He/D, reverse loop) live in reaction_blocks.h.
// ---------------------------------------------------------------------------
#include <algorithm>

#include "kinetics/reaction_index.h"
#include "cooling/grain.h"
#include "models/metal_grain/reaction_metal_grain_charge.h"
#include "models/metal_grain/reaction_metal_grain_surface.h"
#include "models/primordial/rate_laws.h"

namespace arche {
// ---------------------------------------------------------------------------
// compute_metal_rates — UMIST metal reactions k_rxn[100..542], k_rxn[600..644],
//                       and the associated zero-forcing list.
// Metal_grain only. Parameters: T_K [K], T300 = T_K/300.
//
// ---------------------------------------------------------------------------
template <class K>
inline void compute_metal_rates(K& k_rxn, double T_K, double T300) {
  if (k_rxn.want(100)) k_rxn[100] = 1.31e-10 * std::exp(-80.0 / T_K);
  if (k_rxn.want(101)) k_rxn[101] = 6.00e-9 * std::exp(-40200.0 / T_K);
  if (k_rxn.want(102)) k_rxn[102] = 2.2e-10;
  if (k_rxn.want(103)) k_rxn[103] = 0.0;
  if (k_rxn.want(104)) k_rxn[104] = 0.0;
  if (k_rxn.want(105))
    k_rxn[105] = 6.99e-14 * std::pow(T300, 2.80) * std::exp(-1950.0 / T_K);
  if (k_rxn.want(106)) k_rxn[106] = 6.00e-9 * std::exp(-50900.0 / T_K);
  if (k_rxn.want(107)) k_rxn[107] = 0.0;
  if (k_rxn.want(108)) k_rxn[108] = 5.80e-9 * std::exp(-52900.0 / T_K);
  if (k_rxn.want(109)) k_rxn[109] = 0.0;
  if (k_rxn.want(110)) k_rxn[110] = 0.0;
  if (k_rxn.want(111))
    k_rxn[111] = 4.85e-12 * std::pow(T300, 1.9) * std::exp(-1379.0 / T_K);
  if (k_rxn.want(112)) k_rxn[112] = 0.0;
  if (k_rxn.want(113)) k_rxn[113] = 6.00e-9 * std::exp(-52300.0 / T_K);
  if (k_rxn.want(114)) k_rxn[114] = 5.00e-11 * std::exp(-866.0 / T_K);
  if (k_rxn.want(115))
    k_rxn[115] = 2.06e-11 * std::pow(T300, 0.84) * std::exp(-277.0 / T_K);
  if (k_rxn.want(116)) k_rxn[116] = 1.66e-10 * std::exp(-413.0 / T_K);
  if (k_rxn.want(117)) k_rxn[117] = 8.00e-11 * std::exp(-4000.0 / T_K);
  if (k_rxn.want(118)) k_rxn[118] = 0.0;
  if (k_rxn.want(119)) k_rxn[119] = 0.0;
  if (k_rxn.want(120)) k_rxn[120] = 0.0;
  if (k_rxn.want(121)) k_rxn[121] = 0.0;
  if (k_rxn.want(122)) k_rxn[122] = 0.0;
  if (k_rxn.want(123)) k_rxn[123] = 2.52e-11 * std::exp(-2381.0 / T_K);
  if (k_rxn.want(124)) k_rxn[124] = 4.98e-10 * std::exp(-6000.0 / T_K);
  if (k_rxn.want(125)) k_rxn[125] = 5.01e-11;
  if (k_rxn.want(126)) k_rxn[126] = 0.0;
  if (k_rxn.want(127)) k_rxn[127] = 0.0;
  if (k_rxn.want(128))
    k_rxn[128] = 1.07e-11 * std::pow(T300, 1.17) * std::exp(-1242.0 / T_K);
  if (k_rxn.want(129))
    k_rxn[129] = 8.54e-14 * std::pow(T300, 3.25) * std::exp(-1200.0 / T_K);
  if (k_rxn.want(130)) k_rxn[130] = 0.0;
  if (k_rxn.want(131)) k_rxn[131] = 0.0;
  if (k_rxn.want(132)) k_rxn[132] = 0.0;
  if (k_rxn.want(133)) k_rxn[133] = 6.00e-9 * std::exp(-40200.0 / T_K);
  if (k_rxn.want(134))
    k_rxn[134] = 5.18e-11 * std::pow(T300, 0.17) * std::exp(-6400.0 / T_K);
  if (k_rxn.want(135))
    k_rxn[135] = 6.86e-14 * std::pow(T300, 2.74) * std::exp(-4740.0 / T_K);
  if (k_rxn.want(136))
    k_rxn[136] = 2.05e-12 * std::pow(T300, 1.52) * std::exp(-1736.0 / T_K);
  if (k_rxn.want(137)) k_rxn[137] = 6.00e-9 * std::exp(-50900.0 / T_K);
  if (k_rxn.want(138)) k_rxn[138] = 5.80e-9 * std::exp(-52900.0 / T_K);
  if (k_rxn.want(139)) k_rxn[139] = 0.0;
  if (k_rxn.want(140)) k_rxn[140] = 3.16e-10 * std::exp(-21890.0 / T_K);
  if (k_rxn.want(141)) k_rxn[141] = 6.00e-9 * std::exp(-52300.0 / T_K);
  if (k_rxn.want(142)) k_rxn[142] = 0.0;
  if (k_rxn.want(143)) k_rxn[143] = 0.0;  // HD+grain->HD (grain surface)
  if (k_rxn.want(144)) k_rxn[144] = 0.0;
  if (k_rxn.want(145)) k_rxn[145] = 0.0;
  if (k_rxn.want(146)) k_rxn[146] = 0.0;  // CH+CH4->CH2+CH3 (commented out)
  if (k_rxn.want(147))
    k_rxn[147] = 1.44e-11 * std::pow(T300, 0.50) * std::exp(-5000.0 / T_K);
  if (k_rxn.want(148))
    k_rxn[148] = 2.87e-12 * std::pow(T300, 0.70) * std::exp(-500.0 / T_K);
  if (k_rxn.want(149))
    k_rxn[149] = 9.21e-12 * std::pow(T300, 0.70) * std::exp(-2000.0 / T_K);
  // k_rxn[150] CH + O2 -> HCO + O  (piece-wise)
  if (k_rxn.want(150)) {
    if (T_K <= 10.0)
      k_rxn[150] = 7.6e-12 * std::pow(10.0 / 300.0, -0.48);
    else if (T_K >= 300.0)
      k_rxn[150] = 7.6e-12;
    else
      k_rxn[150] = 7.6e-12 * std::pow(T300, -0.48);
  }
  if (k_rxn.want(151))
    k_rxn[151] = 2.94e-13 * std::pow(T300, 0.50) * std::exp(-7550.0 / T_K);
  if (k_rxn.want(152))
    k_rxn[152] = 1.44e-11 * std::pow(T300, 0.50) * std::exp(-3000.0 / T_K);
  if (k_rxn.want(153))
    k_rxn[153] = 2.94e-13 * std::pow(T300, 0.50) * std::exp(-3000.0 / T_K);
  if (k_rxn.want(154)) k_rxn[154] = 4.00e-10 * std::exp(-5000.0 / T_K);
  if (k_rxn.want(155)) k_rxn[155] = 7.13e-12 * std::exp(-5050.0 / T_K);
  if (k_rxn.want(156)) k_rxn[156] = 3.00e-11;
  if (k_rxn.want(157)) k_rxn[157] = 3.30e-13 * std::exp(-3270.0 / T_K);
  if (k_rxn.want(158)) k_rxn[158] = 0.0;
  if (k_rxn.want(159))
    k_rxn[159] = 1.34e-15 * std::pow(T300, 5.05) * std::exp(-1636.0 / T_K);
  if (k_rxn.want(160)) k_rxn[160] = 0.0;
  if (k_rxn.want(161))
    k_rxn[161] = 1.44e-11 * std::pow(T300, 0.50) * std::exp(-3000.0 / T_K);
  if (k_rxn.want(162))
    k_rxn[162] = 1.44e-11 * std::pow(T300, 0.50) * std::exp(-3000.0 / T_K);
  if (k_rxn.want(163)) k_rxn[163] = 3.00e-11;
  if (k_rxn.want(164)) k_rxn[164] = 1.20e-10 * std::exp(-1400.0 / T_K);
  if (k_rxn.want(165))
    k_rxn[165] = 3.77e-13 * std::pow(T300, 2.42) * std::exp(-1162.0 / T_K);
  if (k_rxn.want(166))
    k_rxn[166] = 1.65e-12 * std::pow(T300, 1.14) * std::exp(-50.0 / T_K);
  if (k_rxn.want(167)) k_rxn[167] = 2.81e-13 * std::exp(-176.0 / T_K);
  if (k_rxn.want(168)) k_rxn[168] = 5.26e-12 * std::exp(-307.0 / T_K);
  if (k_rxn.want(169)) k_rxn[169] = 0.0;
  if (k_rxn.want(170)) k_rxn[170] = 5.99e-12 * std::exp(-24075.0 / T_K);
  if (k_rxn.want(171)) k_rxn[171] = 5.60e-10 * std::exp(-12160.0 / T_K);
  if (k_rxn.want(172)) k_rxn[172] = 4.10e-11 * std::exp(-750.0 / T_K);
  if (k_rxn.want(173))
    k_rxn[173] = 3.65e-11 * std::pow(T300, -3.3) * std::exp(-1443.0 / T_K);
  if (k_rxn.want(174)) k_rxn[174] = 5.30e-12 * std::exp(-34975.0 / T_K);
  if (k_rxn.want(175)) k_rxn[175] = 5.64e-13 * std::exp(-4500.0 / T_K);
  if (k_rxn.want(176)) k_rxn[176] = 0.0;
  if (k_rxn.want(177))
    k_rxn[177] = 4.64e-12 * std::pow(T300, 0.7) * std::exp(25.6 / T_K);
  if (k_rxn.want(178)) k_rxn[178] = 6.00e-12;
  if (k_rxn.want(179)) k_rxn[179] = 0.0;
  if (k_rxn.want(180)) k_rxn[180] = 5.00e-11;
  if (k_rxn.want(181)) k_rxn[181] = 3.30e-12 * std::exp(-5870.0 / T_K);
  if (k_rxn.want(182)) k_rxn[182] = 1.0e-12;
  if (k_rxn.want(183)) k_rxn[183] = 1.50e-10;
  if (k_rxn.want(184)) k_rxn[184] = 1.00e-17;
  if (k_rxn.want(185))
    k_rxn[185] = 4.36e-18 * std::pow(T300, 0.35) * std::exp(-161.3 / T_K);
  if (k_rxn.want(186))
    k_rxn[186] = 4.69e-19 * std::pow(T300, 1.52) * std::exp(50.50 / T_K);
  if (k_rxn.want(187)) k_rxn[187] = 6.59e-11;
  if (k_rxn.want(188)) k_rxn[188] = 1.00e-10;
  if (k_rxn.want(189)) k_rxn[189] = 1.00e-10;
  if (k_rxn.want(190))
    k_rxn[190] = 5.56e-11 * std::pow(T300, 0.41) * std::exp(26.9 / T_K);
  if (k_rxn.want(191)) k_rxn[191] = 2.00e-10;
  if (k_rxn.want(192)) k_rxn[192] = 9.90e-19 * std::pow(T300, -0.38);
  if (k_rxn.want(193)) k_rxn[193] = 4.90e-20 * std::pow(T300, 1.58);
  // k_rxn[194] O + CH -> CO + H  (piece-wise 2000K)
  if (k_rxn.want(194))
    k_rxn[194] = (T_K < 2000.0)
                     ? 6.02e-11 * std::pow(T300, 0.10) * std::exp(4.50 / T_K)
                     : 1.02e-10 * std::exp(-914.0 / T_K);
  if (k_rxn.want(195))
    k_rxn[195] = 1.09e-11 * std::pow(T300, -2.19) * std::exp(-165.1 / T_K);
  if (k_rxn.want(196)) k_rxn[196] = 1.33e-10;
  if (k_rxn.want(197)) k_rxn[197] = 1.30e-10;
  // k_rxn[198] O + OH -> O2 + H: absent from the source network, so the slot
  // keeps the zero from fill.
  if (k_rxn.want(199)) k_rxn[199] = 2.00e-10 * std::pow(T300, -0.12);
  if (k_rxn.want(200)) k_rxn[200] = 5.00e-11;
  if (k_rxn.want(201)) k_rxn[201] = 5.00e-11;
  // k_rxn[202] O + O2H -> OH + O2  (piece-wise 200K)
  if (k_rxn.want(202)) {
    if (T_K > 200.0)
      k_rxn[202] = 5.76e-11 * std::pow(T300, -0.3) * std::exp(-7.5 / T_K);
    else
      k_rxn[202] =
          5.76e-11 * std::pow(200.0 / 300.0, -0.3) * std::exp(-7.5 / 200.0);
  }
  if (k_rxn.want(203)) k_rxn[203] = 1.70e-10;
  // k_rxn[204] OH + H2CO -> H2O + HCO  (piece-wise 200K)
  if (k_rxn.want(204)) {
    if (T_K > 200.0)
      k_rxn[204] = 7.76e-12 * std::pow(T300, 0.82) * std::exp(30.6 / T_K);
    else
      k_rxn[204] =
          7.76e-12 * std::pow(200.0 / 300.0, 0.82) * std::exp(30.6 / 200.0);
  }
  // k_rxn[205] OH + O2H -> H2O + O2  (piece-wise 200K)
  if (k_rxn.want(205)) {
    if (T_K > 200.0)
      k_rxn[205] = 8.58e-11 * std::pow(T300, -0.56) * std::exp(-14.8 / T_K);
    else
      k_rxn[205] =
          8.58e-11 * std::pow(200.0 / 300.0, -0.56) * std::exp(-14.8 / 200.0);
  }
  if (k_rxn.want(206)) k_rxn[206] = 3.00e-11;
  if (k_rxn.want(207)) k_rxn[207] = 1.90e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(208)) k_rxn[208] = 1.40e-9;
  if (k_rxn.want(209)) k_rxn[209] = 1.40e-9;
  if (k_rxn.want(210)) k_rxn[210] = 3.40e-9;
  if (k_rxn.want(211)) k_rxn[211] = 2.30e-9;
  if (k_rxn.want(212)) k_rxn[212] = 1.50e-9;
  if (k_rxn.want(213)) k_rxn[213] = 2.10e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(214)) k_rxn[214] = 6.90e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(215)) k_rxn[215] = 3.10e-9;
  if (k_rxn.want(216)) k_rxn[216] = 9.40e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(217)) k_rxn[217] = 9.40e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(218)) k_rxn[218] = 9.40e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(219)) k_rxn[219] = 2.96e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(220)) k_rxn[220] = 3.57e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(221)) k_rxn[221] = 2.00e-9;
  if (k_rxn.want(222)) k_rxn[222] = 3.50e-9;
  if (k_rxn.want(223)) k_rxn[223] = 1.00e-9;
  if (k_rxn.want(224)) k_rxn[224] = 1.00e-9;
  if (k_rxn.want(225)) k_rxn[225] = 1.00e-10;
  if (k_rxn.want(226)) k_rxn[226] = 1.00e-9;
  if (k_rxn.want(227)) k_rxn[227] = 1.00e-9;
  if (k_rxn.want(228)) k_rxn[228] = 1.00e-10;
  if (k_rxn.want(229)) k_rxn[229] = 2.00e-11;
  if (k_rxn.want(230)) k_rxn[230] = 1.00e-9;
  if (k_rxn.want(231)) k_rxn[231] = 2.40e-9;
  if (k_rxn.want(232)) k_rxn[232] = 1.50e-9;
  if (k_rxn.want(233)) k_rxn[233] = 7.10e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(234)) k_rxn[234] = 7.10e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(235)) k_rxn[235] = 1.00e-9;
  if (k_rxn.want(236)) k_rxn[236] = 1.00e-9;
  if (k_rxn.want(237)) k_rxn[237] = 1.40e-9;
  if (k_rxn.want(238)) k_rxn[238] = 1.14e-10;
  if (k_rxn.want(239)) k_rxn[239] = 2.30e-9;
  if (k_rxn.want(240)) k_rxn[240] = 7.60e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(241)) k_rxn[241] = 7.60e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(242)) k_rxn[242] = 3.90e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(243)) k_rxn[243] = 3.40e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(244)) k_rxn[244] = 1.10e-9;
  if (k_rxn.want(245)) k_rxn[245] = 6.44e-10;
  if (k_rxn.want(246)) k_rxn[246] = 2.16e-9;
  if (k_rxn.want(247)) k_rxn[247] = 1.00e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(248)) k_rxn[248] = 1.00e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(249)) k_rxn[249] = 1.40e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(250)) k_rxn[250] = 1.40e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(251)) k_rxn[251] = 1.90e-9;
  if (k_rxn.want(252)) k_rxn[252] = 8.00e-10;
  if (k_rxn.want(253)) k_rxn[253] = 2.35e-9;
  if (k_rxn.want(254)) k_rxn[254] = 2.00e-9;
  // k_rxn[255] H3+ + O -> OH+ + H2  (piece-wise 400K)
  if (k_rxn.want(255)) {
    if (T_K < 400.0)
      k_rxn[255] = 7.98e-10 * std::pow(T300, -0.16) * std::exp(-1.4 / T_K);
    else
      k_rxn[255] =
          7.98e-10 * std::pow(400.0 / 300.0, -0.16) * std::exp(-1.4 / 400.0);
  }
  if (k_rxn.want(256)) k_rxn[256] = 1.20e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(257)) k_rxn[257] = 1.70e-9;
  if (k_rxn.want(258)) k_rxn[258] = 2.10e-9;
  if (k_rxn.want(259)) k_rxn[259] = 2.40e-9;
  if (k_rxn.want(260)) k_rxn[260] = 1.30e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(261)) k_rxn[261] = 5.90e-9 * std::pow(T300, -0.5);
  // k_rxn[262] H3+ + CO -> HCO+ + H2  (piece-wise 400K)
  if (k_rxn.want(262)) {
    if (T_K < 400.0)
      k_rxn[262] = 1.36e-9 * std::pow(T300, -0.14) * std::exp(3.4 / T_K);
    else
      k_rxn[262] =
          1.36e-9 * std::pow(400.0 / 300.0, -0.14) * std::exp(3.4 / 400.0);
  }
  if (k_rxn.want(263)) k_rxn[263] = 1.70e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(264)) k_rxn[264] = 6.30e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(265)) k_rxn[265] = 2.00e-9;
  if (k_rxn.want(266)) k_rxn[266] = 1.10e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(267)) k_rxn[267] = 5.00e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(268)) k_rxn[268] = 7.50e-10;
  if (k_rxn.want(269)) k_rxn[269] = 7.50e-10;
  if (k_rxn.want(270)) k_rxn[270] = 1.80e-9;
  if (k_rxn.want(271)) k_rxn[271] = 4.80e-10;
  if (k_rxn.want(272)) k_rxn[272] = 2.40e-10;
  if (k_rxn.want(273)) k_rxn[273] = 9.50e-10;
  if (k_rxn.want(274)) k_rxn[274] = 8.50e-11;
  if (k_rxn.want(275)) k_rxn[275] = 5.10e-11;
  if (k_rxn.want(276)) k_rxn[276] = 1.10e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(277)) k_rxn[277] = 2.04e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(278)) k_rxn[278] = 2.86e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(279)) k_rxn[279] = 6.05e-11 * std::pow(T300, -0.5);
  if (k_rxn.want(280)) k_rxn[280] = 1.60e-9;
  if (k_rxn.want(281)) k_rxn[281] = 5.00e-10;
  if (k_rxn.want(282)) k_rxn[282] = 1.60e-9;
  if (k_rxn.want(283)) k_rxn[283] = 4.90e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(284)) k_rxn[284] = 4.90e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(285)) k_rxn[285] = 3.00e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(286)) k_rxn[286] = 1.88e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(287)) k_rxn[287] = 1.14e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(288)) k_rxn[288] = 1.10e-9;
  if (k_rxn.want(289)) k_rxn[289] = 3.30e-11;
  if (k_rxn.want(290)) k_rxn[290] = 1.10e-11;
  if (k_rxn.want(291)) k_rxn[291] = 1.00e-10;
  if (k_rxn.want(292)) k_rxn[292] = 8.70e-10;
  if (k_rxn.want(293)) k_rxn[293] = 4.00e-11;
  if (k_rxn.want(294)) k_rxn[294] = 1.70e-17;
  if (k_rxn.want(295))
    k_rxn[295] = 3.14e-18 * std::pow(T300, -0.15) * std::exp(-68.0 / T_K);
  if (k_rxn.want(296)) k_rxn[296] = 7.51e-8 * std::pow(T300, -0.50);
  if (k_rxn.want(297))
    k_rxn[297] = 2.00e-16 * std::pow(T300, -1.30) * std::exp(-23.0 / T_K);
  if (k_rxn.want(298)) k_rxn[298] = 3.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(299)) k_rxn[299] = 3.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(300)) k_rxn[300] = 5.20e-10;
  if (k_rxn.want(301)) k_rxn[301] = 7.70e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(302)) k_rxn[302] = 9.00e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(303)) k_rxn[303] = 4.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(304)) k_rxn[304] = 4.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(305)) k_rxn[305] = 2.34e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(306)) k_rxn[306] = 7.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(307)) k_rxn[307] = 7.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(308)) k_rxn[308] = 3.42e-10;
  if (k_rxn.want(309)) k_rxn[309] = 4.54e-10;
  if (k_rxn.want(310)) k_rxn[310] = 1.10e-9;
  // k_rxn[311] CH+ + H -> C+ + H2  (piece-wise 1000K)
  if (k_rxn.want(311)) {
    if (T_K < 1000.0)
      k_rxn[311] = 9.06e-10 * std::pow(T300, -0.37) * std::exp(-29.1 / T_K);
    else
      k_rxn[311] =
          9.06e-10 * std::pow(1000.0 / 300.0, -0.37) * std::exp(-29.1 / 1000.0);
  }
  if (k_rxn.want(312)) k_rxn[312] = 1.20e-9;
  if (k_rxn.want(313)) k_rxn[313] = 3.50e-10;
  if (k_rxn.want(314)) k_rxn[314] = 1.20e-9;
  if (k_rxn.want(315)) k_rxn[315] = 7.40e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(316)) k_rxn[316] = 7.50e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(317)) k_rxn[317] = 2.90e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(318)) k_rxn[318] = 5.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(319)) k_rxn[319] = 5.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(320)) k_rxn[320] = 0.0;  // none
  if (k_rxn.want(321)) k_rxn[321] = 4.60e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(322)) k_rxn[322] = 4.60e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(323)) k_rxn[323] = 9.60e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(324)) k_rxn[324] = 9.60e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(325)) k_rxn[325] = 9.60e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(326)) k_rxn[326] = 9.70e-10;
  if (k_rxn.want(327)) k_rxn[327] = 1.00e-11;
  if (k_rxn.want(328)) k_rxn[328] = 1.00e-11;
  if (k_rxn.want(329)) k_rxn[329] = 1.60e-9;
  if (k_rxn.want(330)) k_rxn[330] = 7.50e-10;
  if (k_rxn.want(331)) k_rxn[331] = 1.60e-9;
  if (k_rxn.want(332)) k_rxn[332] = 1.20e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(333)) k_rxn[333] = 4.50e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(334)) k_rxn[334] = 2.81e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(335)) k_rxn[335] = 9.10e-10;
  if (k_rxn.want(336)) k_rxn[336] = 1.60e-9;
  if (k_rxn.want(337)) k_rxn[337] = 4.00e-10;
  if (k_rxn.want(338)) k_rxn[338] = 4.00e-11;
  if (k_rxn.want(339))
    k_rxn[339] = 3.92e-16 * std::pow(T300, -2.29) * std::exp(-21.3 / T_K);
  if (k_rxn.want(340)) k_rxn[340] = 7.20e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(341)) k_rxn[341] = 4.40e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(342)) k_rxn[342] = 4.40e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(343)) k_rxn[343] = 1.60e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(344)) k_rxn[344] = 5.00e-12;
  if (k_rxn.want(345))
    k_rxn[345] = 5.66e-10 * std::pow(T300, 0.36) * std::exp(8.6 / T_K);
  if (k_rxn.want(346)) k_rxn[346] = 7.51e-8 * std::pow(T300, -0.50);
  if (k_rxn.want(347)) k_rxn[347] = 1.70e-9;
  if (k_rxn.want(348)) k_rxn[348] = 3.50e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(349)) k_rxn[349] = 3.50e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(350)) k_rxn[350] = 9.70e-10;
  if (k_rxn.want(351)) k_rxn[351] = 1.10e-10;
  if (k_rxn.want(352)) k_rxn[352] = 8.90e-10;
  if (k_rxn.want(353)) k_rxn[353] = 3.60e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(354)) k_rxn[354] = 3.60e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(355)) k_rxn[355] = 3.20e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(356)) k_rxn[356] = 4.80e-10;
  if (k_rxn.want(357)) k_rxn[357] = 4.80e-10;
  if (k_rxn.want(358)) k_rxn[358] = 4.30e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(359)) k_rxn[359] = 4.30e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(360)) k_rxn[360] = 1.40e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(361)) k_rxn[361] = 2.10e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(362)) k_rxn[362] = 1.90e-11;
  if (k_rxn.want(363)) k_rxn[363] = 9.40e-10;
  if (k_rxn.want(364)) k_rxn[364] = 1.00e-11;
  if (k_rxn.want(365)) k_rxn[365] = 1.00e-9;
  if (k_rxn.want(366))
    k_rxn[366] = 4.89e-11 * std::pow(T300, -0.14) * std::exp(36.1 / T_K);
  if (k_rxn.want(367)) k_rxn[367] = 1.50e-9;
  if (k_rxn.want(368)) k_rxn[368] = 2.60e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(369)) k_rxn[369] = 1.40e-9;
  if (k_rxn.want(370)) k_rxn[370] = 1.62e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(371)) k_rxn[371] = 1.98e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(372)) k_rxn[372] = 3.90e-10;
  if (k_rxn.want(373)) k_rxn[373] = 1.20e-9;
  if (k_rxn.want(374)) k_rxn[374] = 1.20e-9;
  if (k_rxn.want(375)) k_rxn[375] = 7.10e-10;
  if (k_rxn.want(376)) k_rxn[376] = 1.01e-9;
  if (k_rxn.want(377)) k_rxn[377] = 3.50e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(378)) k_rxn[378] = 3.50e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(379)) k_rxn[379] = 4.80e-10;
  if (k_rxn.want(380)) k_rxn[380] = 4.80e-10;
  if (k_rxn.want(381)) k_rxn[381] = 1.31e-9;
  if (k_rxn.want(382)) k_rxn[382] = 1.95e-10;
  if (k_rxn.want(383)) k_rxn[383] = 7.00e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(384)) k_rxn[384] = 1.59e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(385)) k_rxn[385] = 1.30e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(386)) k_rxn[386] = 4.80e-10;
  if (k_rxn.want(387)) k_rxn[387] = 1.05e-9;
  if (k_rxn.want(388)) k_rxn[388] = 2.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(389)) k_rxn[389] = 2.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(390)) k_rxn[390] = 2.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(391)) k_rxn[391] = 1.12e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(392)) k_rxn[392] = 7.44e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(393)) k_rxn[393] = 5.90e-10;
  if (k_rxn.want(394)) k_rxn[394] = 1.44e-9;
  if (k_rxn.want(395)) k_rxn[395] = 0.0;
  if (k_rxn.want(396)) k_rxn[396] = 1.20e-9;
  if (k_rxn.want(397)) k_rxn[397] = 2.20e-10;
  if (k_rxn.want(398)) k_rxn[398] = 4.40e-12;
  if (k_rxn.want(399)) k_rxn[399] = 6.90e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(400)) k_rxn[400] = 9.60e-10;
  if (k_rxn.want(401)) k_rxn[401] = 7.00e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(402)) k_rxn[402] = 3.70e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(403)) k_rxn[403] = 1.00e-9;
  if (k_rxn.want(404)) k_rxn[404] = 8.50e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(405)) k_rxn[405] = 4.50e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(406)) k_rxn[406] = 0.0;
  if (k_rxn.want(407)) k_rxn[407] = 1.10e-9;
  if (k_rxn.want(408)) k_rxn[408] = 4.00e-11;
  if (k_rxn.want(409)) k_rxn[409] = 6.40e-10;
  if (k_rxn.want(410)) k_rxn[410] = 3.40e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(411)) k_rxn[411] = 3.40e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(412)) k_rxn[412] = 4.70e-10;
  if (k_rxn.want(413)) k_rxn[413] = 4.70e-10;
  if (k_rxn.want(414)) k_rxn[414] = 1.40e-9;
  if (k_rxn.want(415)) k_rxn[415] = 6.90e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(416)) k_rxn[416] = 2.10e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(417)) k_rxn[417] = 4.70e-10;
  if (k_rxn.want(418)) k_rxn[418] = 5.00e-10;
  if (k_rxn.want(419)) k_rxn[419] = 2.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(420)) k_rxn[420] = 2.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(421)) k_rxn[421] = 2.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(422)) k_rxn[422] = 1.41e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(423)) k_rxn[423] = 6.62e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(424)) k_rxn[424] = 4.60e-10;
  if (k_rxn.want(425)) k_rxn[425] = 1.00e-11;
  if (k_rxn.want(426)) k_rxn[426] = 7.51e-8 * std::pow(T300, -0.50);
  if (k_rxn.want(427)) k_rxn[427] = 7.51e-8 * std::pow(T300, -0.50);
  if (k_rxn.want(428)) k_rxn[428] = 6.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(429)) k_rxn[429] = 9.40e-10;
  if (k_rxn.want(430)) k_rxn[430] = 3.40e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(431)) k_rxn[431] = 1.10e-10;
  if (k_rxn.want(432)) k_rxn[432] = 3.10e-10;
  if (k_rxn.want(433)) k_rxn[433] = 3.20e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(434)) k_rxn[434] = 4.50e-10;
  if (k_rxn.want(435)) k_rxn[435] = 0.0;
  if (k_rxn.want(436)) k_rxn[436] = 3.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(437)) k_rxn[437] = 8.00e-10;
  if (k_rxn.want(438)) k_rxn[438] = 7.50e-10;
  if (k_rxn.want(439)) k_rxn[439] = 1.10e-10;
  if (k_rxn.want(440)) k_rxn[440] = 1.40e-10;
  if (k_rxn.want(441)) k_rxn[441] = 7.50e-10;
  if (k_rxn.want(442)) k_rxn[442] = 3.20e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(443)) k_rxn[443] = 3.20e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(444)) k_rxn[444] = 4.30e-10;
  if (k_rxn.want(445)) k_rxn[445] = 4.30e-10;
  if (k_rxn.want(446)) k_rxn[446] = 7.93e-10;
  if (k_rxn.want(447)) k_rxn[447] = 4.55e-10;
  if (k_rxn.want(448)) k_rxn[448] = 3.10e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(449)) k_rxn[449] = 3.10e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(450)) k_rxn[450] = 1.72e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(451)) k_rxn[451] = 8.84e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(452)) k_rxn[452] = 8.40e-10;
  if (k_rxn.want(453)) k_rxn[453] = 7.40e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(454)) k_rxn[454] = 1.65e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(455)) k_rxn[455] = 1.35e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(456)) k_rxn[456] = 1.20e-10;
  if (k_rxn.want(457)) k_rxn[457] = 1.10e-9;
  if (k_rxn.want(458)) k_rxn[458] = 3.76e-8 * std::pow(T300, -0.50);
  if (k_rxn.want(459)) k_rxn[459] = 6.30e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(460)) k_rxn[460] = 8.60e-10;
  if (k_rxn.want(461)) k_rxn[461] = 0.0;
  if (k_rxn.want(462)) k_rxn[462] = 1.00e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(463)) k_rxn[463] = 2.50e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(464)) k_rxn[464] = 7.30e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(465)) k_rxn[465] = 3.30e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(466)) k_rxn[466] = 3.10e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(467)) k_rxn[467] = 3.10e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(468)) k_rxn[468] = 4.30e-10;
  if (k_rxn.want(469)) k_rxn[469] = 4.30e-10;
  if (k_rxn.want(470)) k_rxn[470] = 9.35e-11;
  if (k_rxn.want(471)) k_rxn[471] = 2.60e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(472)) k_rxn[472] = 3.60e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(473)) k_rxn[473] = 3.60e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(474)) k_rxn[474] = 3.20e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(475)) k_rxn[475] = 7.70e-11;
  if (k_rxn.want(476)) k_rxn[476] = 6.20e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(477)) k_rxn[477] = 0.0;
  if (k_rxn.want(478)) k_rxn[478] = 5.20e-11;
  if (k_rxn.want(479)) k_rxn[479] = 5.20e-11;
  if (k_rxn.want(480)) k_rxn[480] = 3.10e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(481)) k_rxn[481] = 3.10e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(482)) k_rxn[482] = 4.30e-10;
  if (k_rxn.want(483)) k_rxn[483] = 4.30e-10;
  if (k_rxn.want(484)) k_rxn[484] = 4.10e-10;
  if (k_rxn.want(485)) k_rxn[485] = 4.10e-10;
  if (k_rxn.want(486)) k_rxn[486] = 3.60e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(487)) k_rxn[487] = 3.60e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(488)) k_rxn[488] = 2.30e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(489)) k_rxn[489] = 2.07e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(490)) k_rxn[490] = 1.00e-9;
  if (k_rxn.want(491)) k_rxn[491] = 6.20e-10;
  if (k_rxn.want(492)) k_rxn[492] = 6.40e-10;
  if (k_rxn.want(493)) k_rxn[493] = 6.20e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(494)) k_rxn[494] = 8.50e-10;
  if (k_rxn.want(495)) k_rxn[495] = 6.10e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(496)) k_rxn[496] = 8.20e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(497)) k_rxn[497] = 8.40e-10;
  if (k_rxn.want(498)) k_rxn[498] = 7.10e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(499)) k_rxn[499] = 9.80e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(500)) k_rxn[500] = 1.10e-9;
  if (k_rxn.want(501)) k_rxn[501] = 1.00e-9;
  if (k_rxn.want(502)) k_rxn[502] = 1.00e-9;
  if (k_rxn.want(503)) k_rxn[503] = 7.80e-10;
  if (k_rxn.want(504)) k_rxn[504] = 2.30e-9 * std::pow(T300, -0.50);
  if (k_rxn.want(505)) k_rxn[505] = 7.80e-10;
  if (k_rxn.want(506))
    k_rxn[506] = 2.36e-12 * std::pow(T300, -0.29) * std::exp(17.6 / T_K);
  if (k_rxn.want(507)) k_rxn[507] = 1.50e-7 * std::pow(T300, -0.42);
  if (k_rxn.want(508)) k_rxn[508] = 7.68e-8 * std::pow(T300, -0.60);
  if (k_rxn.want(509)) k_rxn[509] = 1.60e-7 * std::pow(T300, -0.60);
  if (k_rxn.want(510)) k_rxn[510] = 7.75e-8 * std::pow(T300, -0.50);
  if (k_rxn.want(511)) k_rxn[511] = 1.95e-7 * std::pow(T300, -0.50);
  if (k_rxn.want(512)) k_rxn[512] = 2.00e-7 * std::pow(T300, -0.40);
  if (k_rxn.want(513)) k_rxn[513] = 1.10e-10 * std::pow(T300, -0.50);
  if (k_rxn.want(514)) k_rxn[514] = 3.24e-12 * std::pow(T300, -0.66);
  if (k_rxn.want(515)) k_rxn[515] = 1.75e-7 * std::pow(T300, -0.50);
  if (k_rxn.want(516)) k_rxn[516] = 1.75e-7 * std::pow(T300, -0.50);
  if (k_rxn.want(517)) k_rxn[517] = 3.75e-8 * std::pow(T300, -0.50);
  if (k_rxn.want(518)) k_rxn[518] = 1.40e-8 * std::pow(T300, -0.52);
  if (k_rxn.want(519)) k_rxn[519] = 1.40e-8 * std::pow(T300, -0.52);
  if (k_rxn.want(520)) k_rxn[520] = 8.60e-8 * std::pow(T300, -0.50);
  if (k_rxn.want(521)) k_rxn[521] = 3.90e-8 * std::pow(T300, -0.50);
  if (k_rxn.want(522)) k_rxn[522] = 7.09e-8 * std::pow(T300, -0.50);
  if (k_rxn.want(523)) k_rxn[523] = 3.05e-7 * std::pow(T300, -0.50);
  if (k_rxn.want(524)) k_rxn[524] = 3.00e-7 * std::pow(T300, -0.50);
  if (k_rxn.want(525)) k_rxn[525] = 2.00e-7 * std::pow(T300, -0.48);
  if (k_rxn.want(526)) k_rxn[526] = 2.40e-7 * std::pow(T300, -0.69);
  if (k_rxn.want(527)) k_rxn[527] = 1.60e-7 * std::pow(T300, -0.70);
  if (k_rxn.want(528)) k_rxn[528] = 2.50e-7 * std::pow(T300, -0.70);
  if (k_rxn.want(529)) k_rxn[529] = 1.10e-10 * std::pow(T300, -0.70);
  if (k_rxn.want(530)) k_rxn[530] = 2.10e-7 * std::pow(T300, -0.78);
  if (k_rxn.want(531)) k_rxn[531] = 2.17e-7 * std::pow(T300, -0.78);
  if (k_rxn.want(532)) k_rxn[532] = 2.17e-7 * std::pow(T300, -0.78);
  if (k_rxn.want(533)) k_rxn[533] = 1.95e-7 * std::pow(T300, -0.70);
  if (k_rxn.want(534)) k_rxn[534] = 3.00e-7 * std::pow(T300, -0.50);
  if (k_rxn.want(535)) k_rxn[535] = 6.00e-8 * std::pow(T300, -0.64);
  if (k_rxn.want(536)) k_rxn[536] = 3.20e-7 * std::pow(T300, -0.64);
  if (k_rxn.want(537)) k_rxn[537] = 1.00e-17;
  if (k_rxn.want(538))
    k_rxn[538] = 3.27e-14 * std::pow(T300, 2.20) * std::exp(-2240.0 / T_K);
  if (k_rxn.want(539)) k_rxn[539] = 1.06e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(540))
    k_rxn[540] = 6.30e-15 * std::pow(T300, 0.75);  // He+ + C -> C+ + He
  if (k_rxn.want(541)) k_rxn[541] = 9.69e-10 * std::pow(T300, -0.5);
  if (k_rxn.want(542)) k_rxn[542] = 1.71e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(600)) k_rxn[600] = 0.0;
  if (k_rxn.want(601)) k_rxn[601] = 1.70e-11 * std::exp(-1800.0 / T_K);
  if (k_rxn.want(602)) k_rxn[602] = 2.69e-12 * std::exp(-23550.0 / T_K);
  if (k_rxn.want(603)) k_rxn[603] = 7.60e-12;
  if (k_rxn.want(604)) k_rxn[604] = 0.0;  // error duplicate, explicitly zero
  if (k_rxn.want(605))
    k_rxn[605] = 3.65e-11 * std::pow(T300, -3.3) * std::exp(-1443.0 / T_K);
  if (k_rxn.want(606))
    k_rxn[606] = 2.92e-11 * std::pow(T300, -3.3) * std::exp(-1443.0 / T_K);
  if (k_rxn.want(607))
    k_rxn[607] = 2.48e-10 * std::pow(T300, -3.3) * std::exp(-1443.0 / T_K);
  if (k_rxn.want(608)) k_rxn[608] = 0.0;
  if (k_rxn.want(609)) k_rxn[609] = 3.60e-11 * std::exp(-202.0 / T_K);
  if (k_rxn.want(610)) k_rxn[610] = 1.70e-12;
  if (k_rxn.want(611)) k_rxn[611] = 1.66e-12;
  if (k_rxn.want(612)) k_rxn[612] = 1.00e-10;
  if (k_rxn.want(613)) k_rxn[613] = 1.50e-11 * std::exp(-4300.0 / T_K);
  if (k_rxn.want(614)) k_rxn[614] = 3.63e-11;
  if (k_rxn.want(615)) k_rxn[615] = 7.60e-13;
  if (k_rxn.want(616))
    k_rxn[616] = 3.69e-11 * std::pow(T300, -0.27) * std::exp(-12.9 / T_K);
  // k_rxn[620] H3+ + O -> H2O+ + H  (piece-wise 400K)
  if (k_rxn.want(620)) {
    if (T_K <= 400.0)
      k_rxn[620] = 3.42e-10 * std::pow(T300, -0.16) * std::exp(-1.4 / T_K);
    else
      k_rxn[620] =
          3.42e-10 * std::pow(400.0 / 300.0, -0.16) * std::exp(-1.4 / 400.0);
  }
  if (k_rxn.want(625)) k_rxn[625] = 0.0;
  if (k_rxn.want(631)) k_rxn[631] = 4.03e-7 * std::pow(T300, -0.6);
  if (k_rxn.want(632)) k_rxn[632] = 1.96e-7 * std::pow(T300, -0.52);
  if (k_rxn.want(633)) k_rxn[633] = 4.76e-8 * std::pow(T300, -0.52);
  if (k_rxn.want(634)) k_rxn[634] = 8.40e-9 * std::pow(T300, -0.52);
  if (k_rxn.want(635)) k_rxn[635] = 3.05e-7 * std::pow(T300, -0.5);
  if (k_rxn.want(636)) k_rxn[636] = 5.60e-9 * std::pow(T300, -0.5);
  if (k_rxn.want(637)) k_rxn[637] = 5.37e-8 * std::pow(T300, -0.5);
  if (k_rxn.want(639)) k_rxn[639] = 8.10e-7 * std::pow(T300, -0.64);
  // k_rxn[640]=0.0;  // commented out
  if (k_rxn.want(641))
    k_rxn[641] = 5.26e-18 * std::pow(T300, -5.22) * std::exp(-90.0 / T_K);
  if (k_rxn.want(642))
    k_rxn[642] = 5.09e-18 * std::pow(T300, -0.71) * std::exp(-11.6 / T_K);
  if (k_rxn.want(643))
    k_rxn[643] = 4.01e-18 * std::pow(T300, 0.17) * std::exp(-101.5 / T_K);
  // k_rxn[644] C + O+ -> CO+ + ph.  (piece-wise 2000K)
  if (k_rxn.want(644)) {
    if (T_K > 2000.0)
      k_rxn[644] = 5.0e-10 * std::pow(T300, -3.70) * std::exp(-800.0 / T_K);
    else
      k_rxn[644] =
          5.0e-10 * std::pow(2000.0 / 300.0, -3.70) * std::exp(-800.0 / 2000.0);
  }
}

// ---------------------------------------------------------------------------
// compute_Li_rates — Li reactions k_rxn[800..829] and Li zero-forcing.
// Metal_grain only.
// Parameters: T_K [K], T300 = T_K/300, delE_634 = tbl.reactions[634].delE
// [erg].
// ---------------------------------------------------------------------------
template <class K>
inline void compute_Li_rates(K& k_rxn, double T_K, double T300,
                             double delE_634) {
  constexpr double k_B = phys::k_B;
  if (k_rxn.want(800))
    k_rxn[800] = 1.036e-11 / (std::sqrt(T_K / 107.7) *
                              std::pow(1.0 + std::sqrt(T_K / 107.7), 0.6612) *
                              std::pow(1.0 + std::sqrt(T_K / 1.177e7), 1.3388));
  if (k_rxn.want(801))
    k_rxn[801] = 6.3e-9 * std::pow(T_K, -0.5) * (1.0 + T_K / 14000.0);
  if (k_rxn.want(802)) k_rxn[802] = 2.3e-6 * std::pow(T_K, -0.5);
  if (k_rxn.want(803))
    k_rxn[803] = 6.1e-17 * std::pow(T_K, 0.58) * std::exp(-T_K / 1.72e4);
  if (k_rxn.want(804))
    k_rxn[804] = 2.5e-40 * std::pow(T_K, 7.9) * std::exp(-T_K / 1210.0);
  if (k_rxn.want(805))
    k_rxn[805] = 1.7e-13 * std::pow(T_K, -0.051) * std::exp(-T_K / 282000.0);
  if (k_rxn.want(806)) k_rxn[806] = 4.0e-10;
  if (k_rxn.want(807)) k_rxn[807] = 4.0e-10;
  if (k_rxn.want(808)) k_rxn[808] = 1.0e-11 * std::exp(-67900.0 / T_K);
  if (k_rxn.want(809)) k_rxn[809] = 1.0e-9;
  if (k_rxn.want(810)) k_rxn[810] = 2.0e-12 * T_K * std::exp(-T_K / 1200.0);
  if (k_rxn.want(811))
    k_rxn[811] = 1.4e-20 * std::pow(T_K, -0.9) * std::exp(-T_K / 7000.0);
  if (k_rxn.want(812)) k_rxn[812] = 5.3e-14 * std::pow(T_K, -0.49);
  if (k_rxn.want(813)) k_rxn[813] = 1.0e-9;
  if (k_rxn.want(814))
    k_rxn[814] = 3.9e-6 * std::pow(T_K, -0.70) * std::exp(-T_K / 1200.0);
  if (k_rxn.want(815)) k_rxn[815] = 9.0e-10 * std::exp(-66400.0 / T_K);
  if (k_rxn.want(816))
    k_rxn[816] = 8.7e-10 * std::pow(T_K, 0.040) * std::exp(T_K / 5.92e8);
  if (k_rxn.want(817))
    k_rxn[817] = ratelaw::prim::k_H_Li_to_LiH(T_K);  // shared, held above T*
  if (k_rxn.want(818)) {
    if (T_K < 500.0)
      k_rxn[818] = 6.3e-10 * std::exp(-2553.0 / T_K);
    else
      k_rxn[818] = 7.2e-14 * std::pow(T_K, 1.18) * std::exp(-1470.0 / T_K);
  }
  if (k_rxn.want(819))
    k_rxn[819] = 2.9e-10 * std::pow(T_K, 0.59) -
                 2.6e-10 * std::pow(T_K, 0.6) * std::exp(-400.0 / T_K);
  // Lepp+2002: k_rxn[820..826]
  if (k_rxn.want(820))
    k_rxn[820] = 5.34e-8 * std::pow(T300, -1.23) * std::exp(-T_K / 9.23e5);
  if (k_rxn.want(821))
    k_rxn[821] = 4.83e-11 * std::pow(T300, -0.621) * std::exp(-T_K / 1.67e6);
  if (k_rxn.want(822))
    k_rxn[822] = 3.71e-7 * std::pow(T300, -0.51) * std::exp(T_K / 4.41e4);
  if (k_rxn.want(823))
    k_rxn[823] = 2.28e-7 * std::pow(T300, -0.51) * std::exp(T_K / 4.41e4);
  if (k_rxn.want(824))
    k_rxn[824] = 3.11e-8 * std::pow(T300, 0.163) * std::exp(-6.27e4 / T_K);
  if (k_rxn.want(825))
    k_rxn[825] = 5.67e-12 * std::pow(T300, 0.715) * std::exp(-8.77e5 / T_K);
  if (k_rxn.want(826))
    k_rxn[826] = 1.70e-12 * std::pow(T300, 0.709) * std::exp(-1.42e6 / T_K);
  // Mizusawa+2005: k_rxn[827..828]
  if (k_rxn.want(827)) k_rxn[827] = 2.5e-29 / T_K;
  if (k_rxn.want(828)) k_rxn[828] = 4.1e-30 / T_K;
  // k_rxn[829] H2 + Li -> H2 + Li+ + e  (posit)
  if (k_rxn.want(829))
    k_rxn[829] = 9.9e-9 * std::sqrt(T_K) * std::exp(delE_634 / (k_B * T_K));

  // Li zero-forcing
  if (k_rxn.want(808)) k_rxn[808] = 0.0;
  if (k_rxn.want(819)) k_rxn[819] = 0.0;
}

// ---------------------------------------------------------------------------
// compute_KNaMg_rates — K/Na/Mg reactions k_rxn[700..729].
// Metal_grain only.
// Parameters: T_K [K], T300 = T_K/300,
//             delE_576 = tbl.reactions[576].delE, delE_591 =
//             tbl.reactions[591].delE [erg].
// ---------------------------------------------------------------------------
template <class K>
inline void compute_KNaMg_rates(K& k_rxn, double T_K, double T300,
                                double delE_576, double delE_591) {
  constexpr double k_B = phys::k_B;
  if (k_rxn.want(700)) k_rxn[700] = 3.0e-11 / std::sqrt(T_K);
  if (k_rxn.want(701))
    k_rxn[701] = 9.9e-9 * std::sqrt(T_K) * std::exp(delE_576 / (k_B * T_K));
  if (k_rxn.want(702)) k_rxn[702] = 2.76e-12 * std::pow(T300, -0.68);
  if (k_rxn.want(703)) k_rxn[703] = 1.20e-9;
  if (k_rxn.want(704)) k_rxn[704] = 7.51e-8 * std::pow(T300, -0.50);
  if (k_rxn.want(705)) k_rxn[705] = 2.1e-9;
  if (k_rxn.want(706)) k_rxn[706] = 1.1e-9;
  if (k_rxn.want(707)) k_rxn[707] = 3.5e-10;
  if (k_rxn.want(708)) k_rxn[708] = 3.4e-9;
  if (k_rxn.want(709)) k_rxn[709] = 7.1e-10;
  if (k_rxn.want(710)) k_rxn[710] = 6.2e-9;
  if (k_rxn.want(711)) k_rxn[711] = 3.1e-9;
  if (k_rxn.want(712)) k_rxn[712] = 2.6e-9;
  if (k_rxn.want(713)) k_rxn[713] = 2.6e-9;
  if (k_rxn.want(714)) k_rxn[714] = 2.6e-9;
  if (k_rxn.want(715)) k_rxn[715] = 1.0e-11;
  if (k_rxn.want(716))
    k_rxn[716] = 9.9e-9 * std::sqrt(T_K) * std::exp(delE_591 / (k_B * T_K));
  if (k_rxn.want(717)) k_rxn[717] = 2.78e-12 * std::pow(T300, -0.68);
  if (k_rxn.want(718)) k_rxn[718] = 1.1e-9;
  if (k_rxn.want(719)) k_rxn[719] = 7.51e-8 * std::pow(T300, -0.50);
  if (k_rxn.want(720)) k_rxn[720] = 1.0e-9;
  if (k_rxn.want(721)) k_rxn[721] = 1.1e-9;
  if (k_rxn.want(722)) k_rxn[722] = 3.6e-10;
  if (k_rxn.want(723)) k_rxn[723] = 3.5e-9;
  if (k_rxn.want(724)) k_rxn[724] = 1.4e-9;
  if (k_rxn.want(725)) k_rxn[725] = 1.2e-9;
  if (k_rxn.want(726)) k_rxn[726] = 2.2e-9;
  if (k_rxn.want(727)) k_rxn[727] = 0.0;  // reaction 728 has no rate (= 0)
  if (k_rxn.want(728)) k_rxn[728] = 2.9e-9;
  if (k_rxn.want(729)) k_rxn[729] = 2.9e-9;
}

// ---------------------------------------------------------------------------
// compute_CR_rates_metal — CR + CR-photo reactions for metal_grain.
// k_rxn[543..551] (direct CR), k_rxn[655..681] (CR-photoreactions).
// ---------------------------------------------------------------------------
template <class K>
inline void compute_CR_rates_metal(K& k_rxn, double zeta, double T300) {
  const double zeta_fac = zeta / 1.36e-17;
  if (k_rxn.want(metal_grain::slot::cr_begin))
    k_rxn[metal_grain::slot::cr_begin] = 5.98e-18 * zeta_fac;
  if (k_rxn.want(544)) k_rxn[544] = 6.50e-18 * zeta_fac;
  if (k_rxn.want(545)) k_rxn[545] = 2.30e-17 * zeta_fac;
  if (k_rxn.want(546)) k_rxn[546] = 3.40e-17 * zeta_fac;
  if (k_rxn.want(547)) k_rxn[547] = 2.86e-19 * zeta_fac;
  if (k_rxn.want(548)) k_rxn[548] = 1.20e-17 * zeta_fac;
  if (k_rxn.want(549)) k_rxn[549] = 1.30e-18 * zeta_fac;
  if (k_rxn.want(550)) k_rxn[550] = 3.90e-21 * zeta_fac;
  if (k_rxn.want(551)) k_rxn[551] = 3.90e-17 * zeta_fac;

  constexpr double omega = model::cr_photo_albedo;
  const double cr_ph = 1.3e-17 * zeta_fac / (1.0 - omega);
  if (k_rxn.want(metal_grain::slot::cr_photo_begin))
    k_rxn[metal_grain::slot::cr_photo_begin] = cr_ph * 255.0;
  if (k_rxn.want(656)) k_rxn[656] = cr_ph * 365.0;
  if (k_rxn.want(657)) k_rxn[657] = cr_ph * 88.0;
  if (k_rxn.want(658)) k_rxn[658] = cr_ph * 250.0;
  if (k_rxn.want(659)) k_rxn[659] = cr_ph * 250.0;
  if (k_rxn.want(660)) k_rxn[660] = cr_ph * 250.0;
  if (k_rxn.want(661)) k_rxn[661] = cr_ph * 250.0;
  if (k_rxn.want(662)) k_rxn[662] = cr_ph * 250.0;
  if (k_rxn.want(663)) k_rxn[663] = cr_ph * 1169.5;
  if (k_rxn.want(664)) k_rxn[664] = cr_ph * 254.5;
  if (k_rxn.want(665)) k_rxn[665] = cr_ph * 485.5;
  if (k_rxn.want(666)) k_rxn[666] = cr_ph * 119.5;
  if (k_rxn.want(667))
    k_rxn[667] =
        1.3e-17 * zeta_fac * std::pow(T300, 1.17) * 105.0 / (1.0 - omega);
  if (k_rxn.want(668)) k_rxn[668] = cr_ph * 210.5;
  if (k_rxn.want(669)) k_rxn[669] = cr_ph * 584.5;
  if (k_rxn.want(670)) k_rxn[670] = cr_ph * 1329.5;
  if (k_rxn.want(671)) k_rxn[671] = cr_ph * 375.5;
  if (k_rxn.want(672)) k_rxn[672] = cr_ph * 58.5;
  if (k_rxn.want(673)) k_rxn[673] = cr_ph * 750.0;
  if (k_rxn.want(674)) k_rxn[674] = cr_ph * 854.0;
  if (k_rxn.want(675)) k_rxn[675] = cr_ph * 250.0;
  if (k_rxn.want(676)) k_rxn[676] = cr_ph * 0.2;
  if (k_rxn.want(677)) k_rxn[677] = cr_ph * 0.2;
  if (k_rxn.want(678)) k_rxn[678] = cr_ph * 1.4;
  if (k_rxn.want(679)) k_rxn[679] = cr_ph * 375.0;
  if (k_rxn.want(680)) k_rxn[680] = cr_ph * 8.5;
  if (k_rxn.want(681)) k_rxn[681] = cr_ph * 66.5;
}

// ---------------------------------------------------------------------------
// compute_grain_rates — grain charge transfer + surface reactions.
// Metal_grain only. Maps react_grain_rates/grain_coef_rates output into
// k_rxn[].
// ---------------------------------------------------------------------------
template <class K>
inline void compute_grain_rates(K& k_rxn, double T_K, double T_gr_K, double rho,
                                double nH, double Z_metal, double JH2,
                                double JH2O, double Jtot, double zeta,
                                double T_cr_desorp) {
  constexpr int gcb = metal_grain::slot::grain_charge_begin;
  constexpr int gsb = metal_grain::slot::grain_surface_begin;
  // Grain volume, relative to the cold grains, below which the grain reactions
  // are switched off.  Above it the rates already scale with the remaining grains
  // (surface: alpha_MRN ∝ vol_gr; charge: the grain species decay with vol_gr),
  // so the switch is a negligible step.  Results do not depend on it for values
  // <= 1e-8; keeping the reactions to much smaller volumes only slows the solver.
  constexpr double kGrainVolumeFloor = 1.0e-10;

  // Only the entries the sink keeps are evaluated; kgr[128] also feeds k_rxn[143].
  std::array<double, metal_grain::slot::grain_charge_count> k_charge;
  std::array<double, metal_grain::slot::grain_surface_count> kgr;
  react_grain_rates(T_K, T_gr_K, k_charge,
                    [&](int i) { return k_rxn.want(gcb + i); });
  grain_coef_rates(nH, T_K, T_gr_K, Z_metal, JH2, JH2O, Jtot, zeta, T_cr_desorp,
                   kgr, [&](int i) {
                     return k_rxn.want(gsb + i) ||
                            (i == 128 && k_rxn.want(143));
                   });

  if (detail::vol_gr(rho, T_gr_K) <
      kGrainVolumeFloor * detail::vol_gr(rho, 2.725)) {
    k_charge.fill(0.0);
    kgr.fill(0.0);
  }

  for (int i = 0; i < metal_grain::slot::grain_charge_count; ++i)
    if (k_rxn.want(gcb + i)) k_rxn[gcb + i] = k_charge[i];

  if (k_rxn.want(872)) k_rxn[872] *= 0.6;
  if (k_rxn.want(946)) k_rxn[946] *= 0.6;

  for (int i = 0; i < metal_grain::slot::grain_surface_count; ++i)
    if (k_rxn.want(gsb + i)) k_rxn[gsb + i] = kgr[i];

  // Grain-catalysed HD formation (H + D -> HD, reaction 144) uses the
  // closed-form surface-formation rate (kgr[128]) directly from the gas-phase
  // D abundance, mirroring the H2 grain channel.  Because deuterated species
  // are not depleted onto grains (see grain_coef_rates), HD is formed without
  // routing through physisorbed intermediates, keeping HD in the gas phase.
  if (k_rxn.want(143)) k_rxn[143] = kgr[128];
}

}  // namespace arche
