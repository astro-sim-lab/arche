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
template <int N_react>
inline void compute_metal_rates(std::array<double, 2 * N_react>& k_rxn,
                                double T_K, double T300) {
  k_rxn[100] = 1.31e-10 * std::exp(-80.0 / T_K);
  k_rxn[101] = 6.00e-9 * std::exp(-40200.0 / T_K);
  k_rxn[102] = 2.2e-10;
  k_rxn[103] = 0.0;
  k_rxn[104] = 0.0;
  k_rxn[105] = 6.99e-14 * std::pow(T300, 2.80) * std::exp(-1950.0 / T_K);
  k_rxn[106] = 6.00e-9 * std::exp(-50900.0 / T_K);
  k_rxn[107] = 0.0;
  k_rxn[108] = 5.80e-9 * std::exp(-52900.0 / T_K);
  k_rxn[109] = 0.0;
  k_rxn[110] = 0.0;
  k_rxn[111] = 4.85e-12 * std::pow(T300, 1.9) * std::exp(-1379.0 / T_K);
  k_rxn[112] = 0.0;
  k_rxn[113] = 6.00e-9 * std::exp(-52300.0 / T_K);
  k_rxn[114] = 5.00e-11 * std::exp(-866.0 / T_K);
  k_rxn[115] = 2.06e-11 * std::pow(T300, 0.84) * std::exp(-277.0 / T_K);
  k_rxn[116] = 1.66e-10 * std::exp(-413.0 / T_K);
  k_rxn[117] = 8.00e-11 * std::exp(-4000.0 / T_K);
  k_rxn[118] = 0.0;
  k_rxn[119] = 0.0;
  k_rxn[120] = 0.0;
  k_rxn[121] = 0.0;
  k_rxn[122] = 0.0;
  k_rxn[123] = 2.52e-11 * std::exp(-2381.0 / T_K);
  k_rxn[124] = 4.98e-10 * std::exp(-6000.0 / T_K);
  k_rxn[125] = 5.01e-11;
  k_rxn[126] = 0.0;
  k_rxn[127] = 0.0;
  k_rxn[128] = 1.07e-11 * std::pow(T300, 1.17) * std::exp(-1242.0 / T_K);
  k_rxn[129] = 8.54e-14 * std::pow(T300, 3.25) * std::exp(-1200.0 / T_K);
  k_rxn[130] = 0.0;
  k_rxn[131] = 0.0;
  k_rxn[132] = 0.0;
  k_rxn[133] = 6.00e-9 * std::exp(-40200.0 / T_K);
  k_rxn[134] = 5.18e-11 * std::pow(T300, 0.17) * std::exp(-6400.0 / T_K);
  k_rxn[135] = 6.86e-14 * std::pow(T300, 2.74) * std::exp(-4740.0 / T_K);
  k_rxn[136] = 2.05e-12 * std::pow(T300, 1.52) * std::exp(-1736.0 / T_K);
  k_rxn[137] = 6.00e-9 * std::exp(-50900.0 / T_K);
  k_rxn[138] = 5.80e-9 * std::exp(-52900.0 / T_K);
  k_rxn[139] = 0.0;
  k_rxn[140] = 3.16e-10 * std::exp(-21890.0 / T_K);
  k_rxn[141] = 6.00e-9 * std::exp(-52300.0 / T_K);
  k_rxn[142] = 0.0;
  k_rxn[143] = 0.0;  // HD+grain->HD (grain surface)
  k_rxn[144] = 0.0;
  k_rxn[145] = 0.0;
  k_rxn[146] = 0.0;  // CH+CH4->CH2+CH3 (commented out)
  k_rxn[147] = 1.44e-11 * std::pow(T300, 0.50) * std::exp(-5000.0 / T_K);
  k_rxn[148] = 2.87e-12 * std::pow(T300, 0.70) * std::exp(-500.0 / T_K);
  k_rxn[149] = 9.21e-12 * std::pow(T300, 0.70) * std::exp(-2000.0 / T_K);
  // k_rxn[150] CH + O2 -> HCO + O  (piece-wise)
  if (T_K <= 10.0)
    k_rxn[150] = 7.6e-12 * std::pow(10.0 / 300.0, -0.48);
  else if (T_K >= 300.0)
    k_rxn[150] = 7.6e-12;
  else
    k_rxn[150] = 7.6e-12 * std::pow(T300, -0.48);
  k_rxn[151] = 2.94e-13 * std::pow(T300, 0.50) * std::exp(-7550.0 / T_K);
  k_rxn[152] = 1.44e-11 * std::pow(T300, 0.50) * std::exp(-3000.0 / T_K);
  k_rxn[153] = 2.94e-13 * std::pow(T300, 0.50) * std::exp(-3000.0 / T_K);
  k_rxn[154] = 4.00e-10 * std::exp(-5000.0 / T_K);
  k_rxn[155] = 7.13e-12 * std::exp(-5050.0 / T_K);
  k_rxn[156] = 3.00e-11;
  k_rxn[157] = 3.30e-13 * std::exp(-3270.0 / T_K);
  k_rxn[158] = 0.0;
  k_rxn[159] = 1.34e-15 * std::pow(T300, 5.05) * std::exp(-1636.0 / T_K);
  k_rxn[160] = 0.0;
  k_rxn[161] = 1.44e-11 * std::pow(T300, 0.50) * std::exp(-3000.0 / T_K);
  k_rxn[162] = 1.44e-11 * std::pow(T300, 0.50) * std::exp(-3000.0 / T_K);
  k_rxn[163] = 3.00e-11;
  k_rxn[164] = 1.20e-10 * std::exp(-1400.0 / T_K);
  k_rxn[165] = 3.77e-13 * std::pow(T300, 2.42) * std::exp(-1162.0 / T_K);
  k_rxn[166] = 1.65e-12 * std::pow(T300, 1.14) * std::exp(-50.0 / T_K);
  k_rxn[167] = 2.81e-13 * std::exp(-176.0 / T_K);
  k_rxn[168] = 5.26e-12 * std::exp(-307.0 / T_K);
  k_rxn[169] = 0.0;
  k_rxn[170] = 5.99e-12 * std::exp(-24075.0 / T_K);
  k_rxn[171] = 5.60e-10 * std::exp(-12160.0 / T_K);
  k_rxn[172] = 4.10e-11 * std::exp(-750.0 / T_K);
  k_rxn[173] = 3.65e-11 * std::pow(T300, -3.3) * std::exp(-1443.0 / T_K);
  k_rxn[174] = 5.30e-12 * std::exp(-34975.0 / T_K);
  k_rxn[175] = 5.64e-13 * std::exp(-4500.0 / T_K);
  k_rxn[176] = 0.0;
  k_rxn[177] = 4.64e-12 * std::pow(T300, 0.7) * std::exp(25.6 / T_K);
  k_rxn[178] = 6.00e-12;
  k_rxn[179] = 0.0;
  k_rxn[180] = 5.00e-11;
  k_rxn[181] = 3.30e-12 * std::exp(-5870.0 / T_K);
  k_rxn[182] = 1.0e-12;
  k_rxn[183] = 1.50e-10;
  k_rxn[184] = 1.00e-17;
  k_rxn[185] = 4.36e-18 * std::pow(T300, 0.35) * std::exp(-161.3 / T_K);
  k_rxn[186] = 4.69e-19 * std::pow(T300, 1.52) * std::exp(50.50 / T_K);
  k_rxn[187] = 6.59e-11;
  k_rxn[188] = 1.00e-10;
  k_rxn[189] = 1.00e-10;
  k_rxn[190] = 5.56e-11 * std::pow(T300, 0.41) * std::exp(26.9 / T_K);
  k_rxn[191] = 2.00e-10;
  k_rxn[192] = 9.90e-19 * std::pow(T300, -0.38);
  k_rxn[193] = 4.90e-20 * std::pow(T300, 1.58);
  // k_rxn[194] O + CH -> CO + H  (piece-wise 2000K)
  k_rxn[194] = (T_K < 2000.0)
                   ? 6.02e-11 * std::pow(T300, 0.10) * std::exp(4.50 / T_K)
                   : 1.02e-10 * std::exp(-914.0 / T_K);
  k_rxn[195] = 1.09e-11 * std::pow(T300, -2.19) * std::exp(-165.1 / T_K);
  k_rxn[196] = 1.33e-10;
  k_rxn[197] = 1.30e-10;
  // k_rxn[198] O + OH -> O2 + H: absent from the source network, so the slot
  // keeps the zero from fill.
  k_rxn[199] = 2.00e-10 * std::pow(T300, -0.12);
  k_rxn[200] = 5.00e-11;
  k_rxn[201] = 5.00e-11;
  // k_rxn[202] O + O2H -> OH + O2  (piece-wise 200K)
  if (T_K > 200.0)
    k_rxn[202] = 5.76e-11 * std::pow(T300, -0.3) * std::exp(-7.5 / T_K);
  else
    k_rxn[202] =
        5.76e-11 * std::pow(200.0 / 300.0, -0.3) * std::exp(-7.5 / 200.0);
  k_rxn[203] = 1.70e-10;
  // k_rxn[204] OH + H2CO -> H2O + HCO  (piece-wise 200K)
  if (T_K > 200.0)
    k_rxn[204] = 7.76e-12 * std::pow(T300, 0.82) * std::exp(30.6 / T_K);
  else
    k_rxn[204] =
        7.76e-12 * std::pow(200.0 / 300.0, 0.82) * std::exp(30.6 / 200.0);
  // k_rxn[205] OH + O2H -> H2O + O2  (piece-wise 200K)
  if (T_K > 200.0)
    k_rxn[205] = 8.58e-11 * std::pow(T300, -0.56) * std::exp(-14.8 / T_K);
  else
    k_rxn[205] =
        8.58e-11 * std::pow(200.0 / 300.0, -0.56) * std::exp(-14.8 / 200.0);
  k_rxn[206] = 3.00e-11;
  k_rxn[207] = 1.90e-9 * std::pow(T300, -0.5);
  k_rxn[208] = 1.40e-9;
  k_rxn[209] = 1.40e-9;
  k_rxn[210] = 3.40e-9;
  k_rxn[211] = 2.30e-9;
  k_rxn[212] = 1.50e-9;
  k_rxn[213] = 2.10e-9 * std::pow(T300, -0.5);
  k_rxn[214] = 6.90e-9 * std::pow(T300, -0.5);
  k_rxn[215] = 3.10e-9;
  k_rxn[216] = 9.40e-10 * std::pow(T300, -0.5);
  k_rxn[217] = 9.40e-10 * std::pow(T300, -0.5);
  k_rxn[218] = 9.40e-10 * std::pow(T300, -0.5);
  k_rxn[219] = 2.96e-9 * std::pow(T300, -0.5);
  k_rxn[220] = 3.57e-9 * std::pow(T300, -0.5);
  k_rxn[221] = 2.00e-9;
  k_rxn[222] = 3.50e-9;
  k_rxn[223] = 1.00e-9;
  k_rxn[224] = 1.00e-9;
  k_rxn[225] = 1.00e-10;
  k_rxn[226] = 1.00e-9;
  k_rxn[227] = 1.00e-9;
  k_rxn[228] = 1.00e-10;
  k_rxn[229] = 2.00e-11;
  k_rxn[230] = 1.00e-9;
  k_rxn[231] = 2.40e-9;
  k_rxn[232] = 1.50e-9;
  k_rxn[233] = 7.10e-10 * std::pow(T300, -0.5);
  k_rxn[234] = 7.10e-10 * std::pow(T300, -0.5);
  k_rxn[235] = 1.00e-9;
  k_rxn[236] = 1.00e-9;
  k_rxn[237] = 1.40e-9;
  k_rxn[238] = 1.14e-10;
  k_rxn[239] = 2.30e-9;
  k_rxn[240] = 7.60e-10 * std::pow(T300, -0.5);
  k_rxn[241] = 7.60e-10 * std::pow(T300, -0.5);
  k_rxn[242] = 3.90e-9 * std::pow(T300, -0.5);
  k_rxn[243] = 3.40e-9 * std::pow(T300, -0.5);
  k_rxn[244] = 1.10e-9;
  k_rxn[245] = 6.44e-10;
  k_rxn[246] = 2.16e-9;
  k_rxn[247] = 1.00e-9 * std::pow(T300, -0.5);
  k_rxn[248] = 1.00e-9 * std::pow(T300, -0.5);
  k_rxn[249] = 1.40e-9 * std::pow(T300, -0.5);
  k_rxn[250] = 1.40e-9 * std::pow(T300, -0.5);
  k_rxn[251] = 1.90e-9;
  k_rxn[252] = 8.00e-10;
  k_rxn[253] = 2.35e-9;
  k_rxn[254] = 2.00e-9;
  // k_rxn[255] H3+ + O -> OH+ + H2  (piece-wise 400K)
  if (T_K < 400.0)
    k_rxn[255] = 7.98e-10 * std::pow(T300, -0.16) * std::exp(-1.4 / T_K);
  else
    k_rxn[255] =
        7.98e-10 * std::pow(400.0 / 300.0, -0.16) * std::exp(-1.4 / 400.0);
  k_rxn[256] = 1.20e-9 * std::pow(T300, -0.5);
  k_rxn[257] = 1.70e-9;
  k_rxn[258] = 2.10e-9;
  k_rxn[259] = 2.40e-9;
  k_rxn[260] = 1.30e-9 * std::pow(T300, -0.5);
  k_rxn[261] = 5.90e-9 * std::pow(T300, -0.5);
  // k_rxn[262] H3+ + CO -> HCO+ + H2  (piece-wise 400K)
  if (T_K < 400.0)
    k_rxn[262] = 1.36e-9 * std::pow(T300, -0.14) * std::exp(3.4 / T_K);
  else
    k_rxn[262] =
        1.36e-9 * std::pow(400.0 / 300.0, -0.14) * std::exp(3.4 / 400.0);
  k_rxn[263] = 1.70e-9 * std::pow(T300, -0.5);
  k_rxn[264] = 6.30e-9 * std::pow(T300, -0.5);
  k_rxn[265] = 2.00e-9;
  k_rxn[266] = 1.10e-9 * std::pow(T300, -0.5);
  k_rxn[267] = 5.00e-10 * std::pow(T300, -0.5);
  k_rxn[268] = 7.50e-10;
  k_rxn[269] = 7.50e-10;
  k_rxn[270] = 1.80e-9;
  k_rxn[271] = 4.80e-10;
  k_rxn[272] = 2.40e-10;
  k_rxn[273] = 9.50e-10;
  k_rxn[274] = 8.50e-11;
  k_rxn[275] = 5.10e-11;
  k_rxn[276] = 1.10e-9 * std::pow(T300, -0.5);
  k_rxn[277] = 2.04e-10 * std::pow(T300, -0.5);
  k_rxn[278] = 2.86e-10 * std::pow(T300, -0.5);
  k_rxn[279] = 6.05e-11 * std::pow(T300, -0.5);
  k_rxn[280] = 1.60e-9;
  k_rxn[281] = 5.00e-10;
  k_rxn[282] = 1.60e-9;
  k_rxn[283] = 4.90e-10 * std::pow(T300, -0.5);
  k_rxn[284] = 4.90e-10 * std::pow(T300, -0.5);
  k_rxn[285] = 3.00e-10 * std::pow(T300, -0.5);
  k_rxn[286] = 1.88e-9 * std::pow(T300, -0.5);
  k_rxn[287] = 1.14e-9 * std::pow(T300, -0.5);
  k_rxn[288] = 1.10e-9;
  k_rxn[289] = 3.30e-11;
  k_rxn[290] = 1.10e-11;
  k_rxn[291] = 1.00e-10;
  k_rxn[292] = 8.70e-10;
  k_rxn[293] = 4.00e-11;
  k_rxn[294] = 1.70e-17;
  k_rxn[295] = 3.14e-18 * std::pow(T300, -0.15) * std::exp(-68.0 / T_K);
  k_rxn[296] = 7.51e-8 * std::pow(T300, -0.50);
  k_rxn[297] = 2.00e-16 * std::pow(T300, -1.30) * std::exp(-23.0 / T_K);
  k_rxn[298] = 3.80e-10 * std::pow(T300, -0.50);
  k_rxn[299] = 3.80e-10 * std::pow(T300, -0.50);
  k_rxn[300] = 5.20e-10;
  k_rxn[301] = 7.70e-10 * std::pow(T300, -0.50);
  k_rxn[302] = 9.00e-10 * std::pow(T300, -0.50);
  k_rxn[303] = 4.80e-10 * std::pow(T300, -0.50);
  k_rxn[304] = 4.80e-10 * std::pow(T300, -0.50);
  k_rxn[305] = 2.34e-9 * std::pow(T300, -0.50);
  k_rxn[306] = 7.80e-10 * std::pow(T300, -0.50);
  k_rxn[307] = 7.80e-10 * std::pow(T300, -0.50);
  k_rxn[308] = 3.42e-10;
  k_rxn[309] = 4.54e-10;
  k_rxn[310] = 1.10e-9;
  // k_rxn[311] CH+ + H -> C+ + H2  (piece-wise 1000K)
  if (T_K < 1000.0)
    k_rxn[311] = 9.06e-10 * std::pow(T300, -0.37) * std::exp(-29.1 / T_K);
  else
    k_rxn[311] =
        9.06e-10 * std::pow(1000.0 / 300.0, -0.37) * std::exp(-29.1 / 1000.0);
  k_rxn[312] = 1.20e-9;
  k_rxn[313] = 3.50e-10;
  k_rxn[314] = 1.20e-9;
  k_rxn[315] = 7.40e-10 * std::pow(T300, -0.50);
  k_rxn[316] = 7.50e-10 * std::pow(T300, -0.50);
  k_rxn[317] = 2.90e-9 * std::pow(T300, -0.50);
  k_rxn[318] = 5.80e-10 * std::pow(T300, -0.50);
  k_rxn[319] = 5.80e-10 * std::pow(T300, -0.50);
  k_rxn[320] = 0.0;  // none
  k_rxn[321] = 4.60e-10 * std::pow(T300, -0.50);
  k_rxn[322] = 4.60e-10 * std::pow(T300, -0.50);
  k_rxn[323] = 9.60e-10 * std::pow(T300, -0.50);
  k_rxn[324] = 9.60e-10 * std::pow(T300, -0.50);
  k_rxn[325] = 9.60e-10 * std::pow(T300, -0.50);
  k_rxn[326] = 9.70e-10;
  k_rxn[327] = 1.00e-11;
  k_rxn[328] = 1.00e-11;
  k_rxn[329] = 1.60e-9;
  k_rxn[330] = 7.50e-10;
  k_rxn[331] = 1.60e-9;
  k_rxn[332] = 1.20e-9 * std::pow(T300, -0.50);
  k_rxn[333] = 4.50e-10 * std::pow(T300, -0.50);
  k_rxn[334] = 2.81e-9 * std::pow(T300, -0.50);
  k_rxn[335] = 9.10e-10;
  k_rxn[336] = 1.60e-9;
  k_rxn[337] = 4.00e-10;
  k_rxn[338] = 4.00e-11;
  k_rxn[339] = 3.92e-16 * std::pow(T300, -2.29) * std::exp(-21.3 / T_K);
  k_rxn[340] = 7.20e-10 * std::pow(T300, -0.50);
  k_rxn[341] = 4.40e-10 * std::pow(T300, -0.50);
  k_rxn[342] = 4.40e-10 * std::pow(T300, -0.50);
  k_rxn[343] = 1.60e-9 * std::pow(T300, -0.50);
  k_rxn[344] = 5.00e-12;
  k_rxn[345] = 5.66e-10 * std::pow(T300, 0.36) * std::exp(8.6 / T_K);
  k_rxn[346] = 7.51e-8 * std::pow(T300, -0.50);
  k_rxn[347] = 1.70e-9;
  k_rxn[348] = 3.50e-10 * std::pow(T300, -0.50);
  k_rxn[349] = 3.50e-10 * std::pow(T300, -0.50);
  k_rxn[350] = 9.70e-10;
  k_rxn[351] = 1.10e-10;
  k_rxn[352] = 8.90e-10;
  k_rxn[353] = 3.60e-10 * std::pow(T300, -0.50);
  k_rxn[354] = 3.60e-10 * std::pow(T300, -0.50);
  k_rxn[355] = 3.20e-9 * std::pow(T300, -0.50);
  k_rxn[356] = 4.80e-10;
  k_rxn[357] = 4.80e-10;
  k_rxn[358] = 4.30e-10 * std::pow(T300, -0.50);
  k_rxn[359] = 4.30e-10 * std::pow(T300, -0.50);
  k_rxn[360] = 1.40e-9 * std::pow(T300, -0.50);
  k_rxn[361] = 2.10e-9 * std::pow(T300, -0.50);
  k_rxn[362] = 1.90e-11;
  k_rxn[363] = 9.40e-10;
  k_rxn[364] = 1.00e-11;
  k_rxn[365] = 1.00e-9;
  k_rxn[366] = 4.89e-11 * std::pow(T300, -0.14) * std::exp(36.1 / T_K);
  k_rxn[367] = 1.50e-9;
  k_rxn[368] = 2.60e-9 * std::pow(T300, -0.50);
  k_rxn[369] = 1.40e-9;
  k_rxn[370] = 1.62e-9 * std::pow(T300, -0.50);
  k_rxn[371] = 1.98e-9 * std::pow(T300, -0.50);
  k_rxn[372] = 3.90e-10;
  k_rxn[373] = 1.20e-9;
  k_rxn[374] = 1.20e-9;
  k_rxn[375] = 7.10e-10;
  k_rxn[376] = 1.01e-9;
  k_rxn[377] = 3.50e-10 * std::pow(T300, -0.50);
  k_rxn[378] = 3.50e-10 * std::pow(T300, -0.50);
  k_rxn[379] = 4.80e-10;
  k_rxn[380] = 4.80e-10;
  k_rxn[381] = 1.31e-9;
  k_rxn[382] = 1.95e-10;
  k_rxn[383] = 7.00e-10 * std::pow(T300, -0.50);
  k_rxn[384] = 1.59e-9 * std::pow(T300, -0.50);
  k_rxn[385] = 1.30e-9 * std::pow(T300, -0.50);
  k_rxn[386] = 4.80e-10;
  k_rxn[387] = 1.05e-9;
  k_rxn[388] = 2.80e-10 * std::pow(T300, -0.50);
  k_rxn[389] = 2.80e-10 * std::pow(T300, -0.50);
  k_rxn[390] = 2.80e-10 * std::pow(T300, -0.50);
  k_rxn[391] = 1.12e-9 * std::pow(T300, -0.50);
  k_rxn[392] = 7.44e-10 * std::pow(T300, -0.50);
  k_rxn[393] = 5.90e-10;
  k_rxn[394] = 1.44e-9;
  k_rxn[395] = 0.0;
  k_rxn[396] = 1.20e-9;
  k_rxn[397] = 2.20e-10;
  k_rxn[398] = 4.40e-12;
  k_rxn[399] = 6.90e-10 * std::pow(T300, -0.50);
  k_rxn[400] = 9.60e-10;
  k_rxn[401] = 7.00e-10 * std::pow(T300, -0.50);
  k_rxn[402] = 3.70e-9 * std::pow(T300, -0.50);
  k_rxn[403] = 1.00e-9;
  k_rxn[404] = 8.50e-10 * std::pow(T300, -0.50);
  k_rxn[405] = 4.50e-9 * std::pow(T300, -0.50);
  k_rxn[406] = 0.0;
  k_rxn[407] = 1.10e-9;
  k_rxn[408] = 4.00e-11;
  k_rxn[409] = 6.40e-10;
  k_rxn[410] = 3.40e-10 * std::pow(T300, -0.50);
  k_rxn[411] = 3.40e-10 * std::pow(T300, -0.50);
  k_rxn[412] = 4.70e-10;
  k_rxn[413] = 4.70e-10;
  k_rxn[414] = 1.40e-9;
  k_rxn[415] = 6.90e-10 * std::pow(T300, -0.50);
  k_rxn[416] = 2.10e-9 * std::pow(T300, -0.50);
  k_rxn[417] = 4.70e-10;
  k_rxn[418] = 5.00e-10;
  k_rxn[419] = 2.80e-10 * std::pow(T300, -0.50);
  k_rxn[420] = 2.80e-10 * std::pow(T300, -0.50);
  k_rxn[421] = 2.80e-10 * std::pow(T300, -0.50);
  k_rxn[422] = 1.41e-9 * std::pow(T300, -0.50);
  k_rxn[423] = 6.62e-10 * std::pow(T300, -0.50);
  k_rxn[424] = 4.60e-10;
  k_rxn[425] = 1.00e-11;
  k_rxn[426] = 7.51e-8 * std::pow(T300, -0.50);
  k_rxn[427] = 7.51e-8 * std::pow(T300, -0.50);
  k_rxn[428] = 6.80e-10 * std::pow(T300, -0.50);
  k_rxn[429] = 9.40e-10;
  k_rxn[430] = 3.40e-9 * std::pow(T300, -0.50);
  k_rxn[431] = 1.10e-10;
  k_rxn[432] = 3.10e-10;
  k_rxn[433] = 3.20e-10 * std::pow(T300, -0.50);
  k_rxn[434] = 4.50e-10;
  k_rxn[435] = 0.0;
  k_rxn[436] = 3.80e-10 * std::pow(T300, -0.50);
  k_rxn[437] = 8.00e-10;
  k_rxn[438] = 7.50e-10;
  k_rxn[439] = 1.10e-10;
  k_rxn[440] = 1.40e-10;
  k_rxn[441] = 7.50e-10;
  k_rxn[442] = 3.20e-10 * std::pow(T300, -0.50);
  k_rxn[443] = 3.20e-10 * std::pow(T300, -0.50);
  k_rxn[444] = 4.30e-10;
  k_rxn[445] = 4.30e-10;
  k_rxn[446] = 7.93e-10;
  k_rxn[447] = 4.55e-10;
  k_rxn[448] = 3.10e-10 * std::pow(T300, -0.50);
  k_rxn[449] = 3.10e-10 * std::pow(T300, -0.50);
  k_rxn[450] = 1.72e-9 * std::pow(T300, -0.50);
  k_rxn[451] = 8.84e-10 * std::pow(T300, -0.50);
  k_rxn[452] = 8.40e-10;
  k_rxn[453] = 7.40e-10 * std::pow(T300, -0.50);
  k_rxn[454] = 1.65e-9 * std::pow(T300, -0.50);
  k_rxn[455] = 1.35e-9 * std::pow(T300, -0.50);
  k_rxn[456] = 1.20e-10;
  k_rxn[457] = 1.10e-9;
  k_rxn[458] = 3.76e-8 * std::pow(T300, -0.50);
  k_rxn[459] = 6.30e-10 * std::pow(T300, -0.50);
  k_rxn[460] = 8.60e-10;
  k_rxn[461] = 0.0;
  k_rxn[462] = 1.00e-9 * std::pow(T300, -0.50);
  k_rxn[463] = 2.50e-9 * std::pow(T300, -0.50);
  k_rxn[464] = 7.30e-10 * std::pow(T300, -0.50);
  k_rxn[465] = 3.30e-9 * std::pow(T300, -0.50);
  k_rxn[466] = 3.10e-10 * std::pow(T300, -0.50);
  k_rxn[467] = 3.10e-10 * std::pow(T300, -0.50);
  k_rxn[468] = 4.30e-10;
  k_rxn[469] = 4.30e-10;
  k_rxn[470] = 9.35e-11;
  k_rxn[471] = 2.60e-9 * std::pow(T300, -0.50);
  k_rxn[472] = 3.60e-10 * std::pow(T300, -0.50);
  k_rxn[473] = 3.60e-10 * std::pow(T300, -0.50);
  k_rxn[474] = 3.20e-9 * std::pow(T300, -0.50);
  k_rxn[475] = 7.70e-11;
  k_rxn[476] = 6.20e-10 * std::pow(T300, -0.50);
  k_rxn[477] = 0.0;
  k_rxn[478] = 5.20e-11;
  k_rxn[479] = 5.20e-11;
  k_rxn[480] = 3.10e-10 * std::pow(T300, -0.50);
  k_rxn[481] = 3.10e-10 * std::pow(T300, -0.50);
  k_rxn[482] = 4.30e-10;
  k_rxn[483] = 4.30e-10;
  k_rxn[484] = 4.10e-10;
  k_rxn[485] = 4.10e-10;
  k_rxn[486] = 3.60e-10 * std::pow(T300, -0.50);
  k_rxn[487] = 3.60e-10 * std::pow(T300, -0.50);
  k_rxn[488] = 2.30e-10 * std::pow(T300, -0.50);
  k_rxn[489] = 2.07e-9 * std::pow(T300, -0.50);
  k_rxn[490] = 1.00e-9;
  k_rxn[491] = 6.20e-10;
  k_rxn[492] = 6.40e-10;
  k_rxn[493] = 6.20e-10 * std::pow(T300, -0.50);
  k_rxn[494] = 8.50e-10;
  k_rxn[495] = 6.10e-10 * std::pow(T300, -0.50);
  k_rxn[496] = 8.20e-10 * std::pow(T300, -0.50);
  k_rxn[497] = 8.40e-10;
  k_rxn[498] = 7.10e-10 * std::pow(T300, -0.50);
  k_rxn[499] = 9.80e-10 * std::pow(T300, -0.50);
  k_rxn[500] = 1.10e-9;
  k_rxn[501] = 1.00e-9;
  k_rxn[502] = 1.00e-9;
  k_rxn[503] = 7.80e-10;
  k_rxn[504] = 2.30e-9 * std::pow(T300, -0.50);
  k_rxn[505] = 7.80e-10;
  k_rxn[506] = 2.36e-12 * std::pow(T300, -0.29) * std::exp(17.6 / T_K);
  k_rxn[507] = 1.50e-7 * std::pow(T300, -0.42);
  k_rxn[508] = 7.68e-8 * std::pow(T300, -0.60);
  k_rxn[509] = 1.60e-7 * std::pow(T300, -0.60);
  k_rxn[510] = 7.75e-8 * std::pow(T300, -0.50);
  k_rxn[511] = 1.95e-7 * std::pow(T300, -0.50);
  k_rxn[512] = 2.00e-7 * std::pow(T300, -0.40);
  k_rxn[513] = 1.10e-10 * std::pow(T300, -0.50);
  k_rxn[514] = 3.24e-12 * std::pow(T300, -0.66);
  k_rxn[515] = 1.75e-7 * std::pow(T300, -0.50);
  k_rxn[516] = 1.75e-7 * std::pow(T300, -0.50);
  k_rxn[517] = 3.75e-8 * std::pow(T300, -0.50);
  k_rxn[518] = 1.40e-8 * std::pow(T300, -0.52);
  k_rxn[519] = 1.40e-8 * std::pow(T300, -0.52);
  k_rxn[520] = 8.60e-8 * std::pow(T300, -0.50);
  k_rxn[521] = 3.90e-8 * std::pow(T300, -0.50);
  k_rxn[522] = 7.09e-8 * std::pow(T300, -0.50);
  k_rxn[523] = 3.05e-7 * std::pow(T300, -0.50);
  k_rxn[524] = 3.00e-7 * std::pow(T300, -0.50);
  k_rxn[525] = 2.00e-7 * std::pow(T300, -0.48);
  k_rxn[526] = 2.40e-7 * std::pow(T300, -0.69);
  k_rxn[527] = 1.60e-7 * std::pow(T300, -0.70);
  k_rxn[528] = 2.50e-7 * std::pow(T300, -0.70);
  k_rxn[529] = 1.10e-10 * std::pow(T300, -0.70);
  k_rxn[530] = 2.10e-7 * std::pow(T300, -0.78);
  k_rxn[531] = 2.17e-7 * std::pow(T300, -0.78);
  k_rxn[532] = 2.17e-7 * std::pow(T300, -0.78);
  k_rxn[533] = 1.95e-7 * std::pow(T300, -0.70);
  k_rxn[534] = 3.00e-7 * std::pow(T300, -0.50);
  k_rxn[535] = 6.00e-8 * std::pow(T300, -0.64);
  k_rxn[536] = 3.20e-7 * std::pow(T300, -0.64);
  k_rxn[537] = 1.00e-17;
  k_rxn[538] = 3.27e-14 * std::pow(T300, 2.20) * std::exp(-2240.0 / T_K);
  k_rxn[539] = 1.06e-9 * std::pow(T300, -0.5);
  k_rxn[540] = 6.30e-15 * std::pow(T300, 0.75);  // He+ + C -> C+ + He
  k_rxn[541] = 9.69e-10 * std::pow(T300, -0.5);
  k_rxn[542] = 1.71e-9 * std::pow(T300, -0.5);
  k_rxn[600] = 0.0;
  k_rxn[601] = 1.70e-11 * std::exp(-1800.0 / T_K);
  k_rxn[602] = 2.69e-12 * std::exp(-23550.0 / T_K);
  k_rxn[603] = 7.60e-12;
  k_rxn[604] = 0.0;  // error duplicate, explicitly zero
  k_rxn[605] = 3.65e-11 * std::pow(T300, -3.3) * std::exp(-1443.0 / T_K);
  k_rxn[606] = 2.92e-11 * std::pow(T300, -3.3) * std::exp(-1443.0 / T_K);
  k_rxn[607] = 2.48e-10 * std::pow(T300, -3.3) * std::exp(-1443.0 / T_K);
  k_rxn[608] = 0.0;
  k_rxn[609] = 3.60e-11 * std::exp(-202.0 / T_K);
  k_rxn[610] = 1.70e-12;
  k_rxn[611] = 1.66e-12;
  k_rxn[612] = 1.00e-10;
  k_rxn[613] = 1.50e-11 * std::exp(-4300.0 / T_K);
  k_rxn[614] = 3.63e-11;
  k_rxn[615] = 7.60e-13;
  k_rxn[616] = 3.69e-11 * std::pow(T300, -0.27) * std::exp(-12.9 / T_K);
  // k_rxn[620] H3+ + O -> H2O+ + H  (piece-wise 400K)
  if (T_K <= 400.0)
    k_rxn[620] = 3.42e-10 * std::pow(T300, -0.16) * std::exp(-1.4 / T_K);
  else
    k_rxn[620] =
        3.42e-10 * std::pow(400.0 / 300.0, -0.16) * std::exp(-1.4 / 400.0);
  k_rxn[625] = 0.0;
  k_rxn[631] = 4.03e-7 * std::pow(T300, -0.6);
  k_rxn[632] = 1.96e-7 * std::pow(T300, -0.52);
  k_rxn[633] = 4.76e-8 * std::pow(T300, -0.52);
  k_rxn[634] = 8.40e-9 * std::pow(T300, -0.52);
  k_rxn[635] = 3.05e-7 * std::pow(T300, -0.5);
  k_rxn[636] = 5.60e-9 * std::pow(T300, -0.5);
  k_rxn[637] = 5.37e-8 * std::pow(T300, -0.5);
  k_rxn[639] = 8.10e-7 * std::pow(T300, -0.64);
  // k_rxn[640]=0.0;  // commented out
  k_rxn[641] = 5.26e-18 * std::pow(T300, -5.22) * std::exp(-90.0 / T_K);
  k_rxn[642] = 5.09e-18 * std::pow(T300, -0.71) * std::exp(-11.6 / T_K);
  k_rxn[643] = 4.01e-18 * std::pow(T300, 0.17) * std::exp(-101.5 / T_K);
  // k_rxn[644] C + O+ -> CO+ + ph.  (piece-wise 2000K)
  if (T_K > 2000.0)
    k_rxn[644] = 5.0e-10 * std::pow(T300, -3.70) * std::exp(-800.0 / T_K);
  else
    k_rxn[644] =
        5.0e-10 * std::pow(2000.0 / 300.0, -3.70) * std::exp(-800.0 / 2000.0);
}

// ---------------------------------------------------------------------------
// compute_Li_rates — Li reactions k_rxn[800..829] and Li zero-forcing.
// Metal_grain only.
// Parameters: T_K [K], T300 = T_K/300, delE_634 = tbl.reactions[634].delE
// [erg].
// ---------------------------------------------------------------------------
template <int N_react>
inline void compute_Li_rates(std::array<double, 2 * N_react>& k_rxn, double T_K,
                             double T300, double delE_634) {
  constexpr double k_B = phys::k_B;
  k_rxn[800] = 1.036e-11 / (std::sqrt(T_K / 107.7) *
                            std::pow(1.0 + std::sqrt(T_K / 107.7), 0.6612) *
                            std::pow(1.0 + std::sqrt(T_K / 1.177e7), 1.3388));
  k_rxn[801] = 6.3e-9 * std::pow(T_K, -0.5) * (1.0 + T_K / 14000.0);
  k_rxn[802] = 2.3e-6 * std::pow(T_K, -0.5);
  k_rxn[803] = 6.1e-17 * std::pow(T_K, 0.58) * std::exp(-T_K / 1.72e4);
  k_rxn[804] = 2.5e-40 * std::pow(T_K, 7.9) * std::exp(-T_K / 1210.0);
  k_rxn[805] = 1.7e-13 * std::pow(T_K, -0.051) * std::exp(-T_K / 282000.0);
  k_rxn[806] = 4.0e-10;
  k_rxn[807] = 4.0e-10;
  k_rxn[808] = 1.0e-11 * std::exp(-67900.0 / T_K);
  k_rxn[809] = 1.0e-9;
  k_rxn[810] = 2.0e-12 * T_K * std::exp(-T_K / 1200.0);
  k_rxn[811] = 1.4e-20 * std::pow(T_K, -0.9) * std::exp(-T_K / 7000.0);
  k_rxn[812] = 5.3e-14 * std::pow(T_K, -0.49);
  k_rxn[813] = 1.0e-9;
  k_rxn[814] = 3.9e-6 * std::pow(T_K, -0.70) * std::exp(-T_K / 1200.0);
  k_rxn[815] = 9.0e-10 * std::exp(-66400.0 / T_K);
  k_rxn[816] = 8.7e-10 * std::pow(T_K, 0.040) * std::exp(T_K / 5.92e8);
  k_rxn[817] = 4.0e-20 * std::exp(-T_K / 4065.0 + std::pow(T_K / 13193.0, 3.0));
  if (T_K < 500.0)
    k_rxn[818] = 6.3e-10 * std::exp(-2553.0 / T_K);
  else
    k_rxn[818] = 7.2e-14 * std::pow(T_K, 1.18) * std::exp(-1470.0 / T_K);
  k_rxn[819] = 2.9e-10 * std::pow(T_K, 0.59) -
               2.6e-10 * std::pow(T_K, 0.6) * std::exp(-400.0 / T_K);
  // Lepp+2002: k_rxn[820..826]
  k_rxn[820] = 5.34e-8 * std::pow(T300, -1.23) * std::exp(-T_K / 9.23e5);
  k_rxn[821] = 4.83e-11 * std::pow(T300, -0.621) * std::exp(-T_K / 1.67e6);
  k_rxn[822] = 3.71e-7 * std::pow(T300, -0.51) * std::exp(T_K / 4.41e4);
  k_rxn[823] = 2.28e-7 * std::pow(T300, -0.51) * std::exp(T_K / 4.41e4);
  k_rxn[824] = 3.11e-8 * std::pow(T300, 0.163) * std::exp(-6.27e4 / T_K);
  k_rxn[825] = 5.67e-12 * std::pow(T300, 0.715) * std::exp(-8.77e5 / T_K);
  k_rxn[826] = 1.70e-12 * std::pow(T300, 0.709) * std::exp(-1.42e6 / T_K);
  // Mizusawa+2005: k_rxn[827..828]
  k_rxn[827] = 2.5e-29 / T_K;
  k_rxn[828] = 4.1e-30 / T_K;
  // k_rxn[829] H2 + Li -> H2 + Li+ + e  (posit)
  k_rxn[829] = 9.9e-9 * std::sqrt(T_K) * std::exp(delE_634 / (k_B * T_K));

  // Li zero-forcing
  k_rxn[808] = 0.0;
  k_rxn[819] = 0.0;
}

// ---------------------------------------------------------------------------
// compute_KNaMg_rates — K/Na/Mg reactions k_rxn[700..729].
// Metal_grain only.
// Parameters: T_K [K], T300 = T_K/300,
//             delE_576 = tbl.reactions[576].delE, delE_591 =
//             tbl.reactions[591].delE [erg].
// ---------------------------------------------------------------------------
template <int N_react>
inline void compute_KNaMg_rates(std::array<double, 2 * N_react>& k_rxn,
                                double T_K, double T300, double delE_576,
                                double delE_591) {
  constexpr double k_B = phys::k_B;
  k_rxn[700] = 3.0e-11 / std::sqrt(T_K);
  k_rxn[701] = 9.9e-9 * std::sqrt(T_K) * std::exp(delE_576 / (k_B * T_K));
  k_rxn[702] = 2.76e-12 * std::pow(T300, -0.68);
  k_rxn[703] = 1.20e-9;
  k_rxn[704] = 7.51e-8 * std::pow(T300, -0.50);
  k_rxn[705] = 2.1e-9;
  k_rxn[706] = 1.1e-9;
  k_rxn[707] = 3.5e-10;
  k_rxn[708] = 3.4e-9;
  k_rxn[709] = 7.1e-10;
  k_rxn[710] = 6.2e-9;
  k_rxn[711] = 3.1e-9;
  k_rxn[712] = 2.6e-9;
  k_rxn[713] = 2.6e-9;
  k_rxn[714] = 2.6e-9;
  k_rxn[715] = 1.0e-11;
  k_rxn[716] = 9.9e-9 * std::sqrt(T_K) * std::exp(delE_591 / (k_B * T_K));
  k_rxn[717] = 2.78e-12 * std::pow(T300, -0.68);
  k_rxn[718] = 1.1e-9;
  k_rxn[719] = 7.51e-8 * std::pow(T300, -0.50);
  k_rxn[720] = 1.0e-9;
  k_rxn[721] = 1.1e-9;
  k_rxn[722] = 3.6e-10;
  k_rxn[723] = 3.5e-9;
  k_rxn[724] = 1.4e-9;
  k_rxn[725] = 1.2e-9;
  k_rxn[726] = 2.2e-9;
  k_rxn[727] = 0.0;  // reaction 728 has no rate (= 0)
  k_rxn[728] = 2.9e-9;
  k_rxn[729] = 2.9e-9;
}

// ---------------------------------------------------------------------------
// compute_CR_rates_metal — CR + CR-photo reactions for metal_grain.
// k_rxn[543..551] (direct CR), k_rxn[655..681] (CR-photoreactions).
// ---------------------------------------------------------------------------
template <int N_react>
inline void compute_CR_rates_metal(std::array<double, 2 * N_react>& k_rxn,
                                   double zeta, double T300) {
  const double zeta_fac = zeta / 1.36e-17;
  k_rxn[metal_grain::slot::cr_begin] = 5.98e-18 * zeta_fac;
  k_rxn[544] = 6.50e-18 * zeta_fac;
  k_rxn[545] = 2.30e-17 * zeta_fac;
  k_rxn[546] = 3.40e-17 * zeta_fac;
  k_rxn[547] = 2.86e-19 * zeta_fac;
  k_rxn[548] = 1.20e-17 * zeta_fac;
  k_rxn[549] = 1.30e-18 * zeta_fac;
  k_rxn[550] = 3.90e-21 * zeta_fac;
  k_rxn[551] = 3.90e-17 * zeta_fac;

  constexpr double omega = model::cr_photo_albedo;
  const double cr_ph = 1.3e-17 * zeta_fac / (1.0 - omega);
  k_rxn[metal_grain::slot::cr_photo_begin] = cr_ph * 255.0;
  k_rxn[656] = cr_ph * 365.0;
  k_rxn[657] = cr_ph * 88.0;
  k_rxn[658] = cr_ph * 250.0;
  k_rxn[659] = cr_ph * 250.0;
  k_rxn[660] = cr_ph * 250.0;
  k_rxn[661] = cr_ph * 250.0;
  k_rxn[662] = cr_ph * 250.0;
  k_rxn[663] = cr_ph * 1169.5;
  k_rxn[664] = cr_ph * 254.5;
  k_rxn[665] = cr_ph * 485.5;
  k_rxn[666] = cr_ph * 119.5;
  k_rxn[667] =
      1.3e-17 * zeta_fac * std::pow(T300, 1.17) * 105.0 / (1.0 - omega);
  k_rxn[668] = cr_ph * 210.5;
  k_rxn[669] = cr_ph * 584.5;
  k_rxn[670] = cr_ph * 1329.5;
  k_rxn[671] = cr_ph * 375.5;
  k_rxn[672] = cr_ph * 58.5;
  k_rxn[673] = cr_ph * 750.0;
  k_rxn[674] = cr_ph * 854.0;
  k_rxn[675] = cr_ph * 250.0;
  k_rxn[676] = cr_ph * 0.2;
  k_rxn[677] = cr_ph * 0.2;
  k_rxn[678] = cr_ph * 1.4;
  k_rxn[679] = cr_ph * 375.0;
  k_rxn[680] = cr_ph * 8.5;
  k_rxn[681] = cr_ph * 66.5;
}

// ---------------------------------------------------------------------------
// compute_grain_rates — grain charge transfer + surface reactions.
// Metal_grain only. Maps react_grain_rates/grain_coef_rates output into
// k_rxn[].
// ---------------------------------------------------------------------------
template <int N_react>
inline void compute_grain_rates(std::array<double, 2 * N_react>& k_rxn,
                                double T_K, double T_gr_K, double rho,
                                double nH, double Z_metal, double JH2,
                                double JH2O, double Jtot, double zeta,
                                double T_cr_desorp) {
  double T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol;
  detail::vaptemp(rho, T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol);
  double T_evap = std::max({T_ir, T_pyr, T_ol});

  std::array<double, metal_grain::slot::grain_charge_count> k_charge;
  std::array<double, metal_grain::slot::grain_surface_count> kgr;
  react_grain_rates(T_K, T_gr_K, k_charge);
  grain_coef_rates(nH, T_K, T_gr_K, Z_metal, JH2, JH2O, Jtot, zeta, T_cr_desorp,
                   kgr);

  if (T_gr_K >= 1.025 * T_evap) {
    k_charge.fill(0.0);
    kgr.fill(0.0);
  }

  for (int i = 0; i < metal_grain::slot::grain_charge_count; ++i)
    k_rxn[metal_grain::slot::grain_charge_begin + i] = k_charge[i];

  k_rxn[872] *= 0.6;
  k_rxn[946] *= 0.6;

  for (int i = 0; i < metal_grain::slot::grain_surface_count; ++i)
    k_rxn[metal_grain::slot::grain_surface_begin + i] = kgr[i];

  // Grain-catalysed HD formation (H + D -> HD, reaction 144) uses the
  // closed-form surface-formation rate (kgr[128]) directly from the gas-phase
  // D abundance, mirroring the H2 grain channel.  Because deuterated species
  // are not depleted onto grains (see grain_coef_rates), HD is formed without
  // routing through physisorbed intermediates, keeping HD in the gas phase.
  k_rxn[143] = kgr[128];
}

}  // namespace arche
