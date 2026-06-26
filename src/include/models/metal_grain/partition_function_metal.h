// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// partition_function_metal.h — partition-function evaluator for the metal_grain
// network (89 species).
//
// The first 16 species share the canonical primordial block (H..D-) and are
// delegated to the shared primordial provider (pf_prim::prim_provider), so that
// physics lives in one place.  Species 16..62 are the metal-specific gas-phase
// species (C/O/K/Na/Mg families and their ions); species 63..88 are grain and
// grain-surface pseudo-species whose partition function is the unit sentinel.
//
// eval_metal_partition_functions() fills pf[0..N_sp-1] and leaves the three
// trailing sentinel slots at 1.0 (vacant / photon / CR reactants act as
// multiplicative no-ops in detailed-balance terms).
// ---------------------------------------------------------------------------
#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

#include "core/species_index.h"
#include "core/species_catalog.h"
#include "core/state.h"
#include "kinetics/partition_function.h"  // PfProvider, pf_eval, detail::*
#include "kinetics/topology.h"            // ReactionTable
#include "models/primordial/partition_function_prim.h"  // pf_prim::prim_provider

namespace arche {

// The shared primordial block occupies metal_grain local indices 0..15 in
// canonical order, so delegating those slots to the primordial provider is a
// straight 1:1 mapping.
static_assert(metal_grain::Species::canonical(0) == SpId::H &&
                  metal_grain::Species::canonical(15) == SpId::Dm,
              "metal_grain species 0..15 must be the canonical primordial "
              "block H..D-");

template <int N_sp, int N_react>
void eval_metal_partition_functions(double T_K,
                                    const ReactionTable<N_sp, N_react>& tbl,
                                    std::array<double, N_sp + 3>& pf) {
  static_assert(N_sp == metal_grain::N_sp,
                "eval_metal_partition_functions is for the metal_grain network");
  constexpr double k_B = phys::k_B;

  pf.fill(1.0);  // default: sentinel value for vacant/photon and grain slots

  // Canonical primordial block (species 0..15 = H..D-) via the shared provider.
  for (int i = 0; i < 16; ++i)
    pf[i] = pf_eval(pf_prim::prim_provider(metal_grain::Species::canonical(i)),
                    T_K, tbl.pf_table[i]);

  // -----------------------------------------------------------------------
  // Metal-grain gas-phase species 16..62.  Species 63..88 are grain/ice
  // pseudo-species → default 1.0 from fill().
  // -----------------------------------------------------------------------
  constexpr double h_P = phys::h_P;  // erg·s
  constexpr double pi_ = phys::pi;
  constexpr double MHz = 1.0e6;

  // lg_th = log10(5040/T) for I88 polynomial fits
  double lg_th = std::log10(5040.0 / T_K);
  double lg_th2 = lg_th * lg_th;
  double lg_th3 = lg_th2 * lg_th;
  double lg_th4 = lg_th3 * lg_th;
  double lg_th5 = lg_th4 * lg_th;
  double lg_th6 = lg_th5 * lg_th;

  // Classical asymmetric-top partition function
  //   pf = sqrt( (kT/h)^3 * pi / (A*B*C) ) * factor
  auto asym_top = [&](double A_MHz, double B_MHz, double C_MHz,
                      double factor = 1.0) {
    double A = A_MHz * MHz, B = B_MHz * MHz, C = C_MHz * MHz;
    double kT_h = k_B * T_K / h_P;
    return std::sqrt(kT_h * kT_h * kT_h * pi_ / (A * B * C)) * factor;
  };

  // pf[16]: C (BC16)
  static constexpr std::array<double, 32> qc = {
      1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00002e+0, 1.00115e+0, 1.02680e+0,
      1.10379e+0, 1.29317e+0, 1.70020e+0, 2.14254e+0, 2.99030e+0, 4.30589e+0,
      5.19094e+0, 6.04754e+0, 6.59516e+0, 7.07437e+0, 7.32544e+0, 7.62484e+0,
      7.83366e+0, 8.27478e+0, 8.47391e+0, 8.62742e+0, 8.74962e+0, 8.81441e+0,
      8.91124e+0, 9.03328e+0, 9.19239e+0, 9.37771e+0, 9.57770e+0, 9.78474e+0,
      9.99517e+0, 1.02090e+1};
  pf[16] = detail::bc16_interp(qc, T_K);

  // pf[17]: C2 (BC16)
  static constexpr std::array<double, 32> qC2 = {
      1.00000e+0, 1.00003e+0, 1.00051e+0, 1.00201e+0, 1.02728e+0, 1.21957e+0,
      1.54110e+0, 2.09649e+0, 3.05116e+0, 4.00914e+0, 5.92681e+0, 9.76432e+0,
      1.36033e+1, 1.93817e+1, 2.53157e+1, 3.40450e+1, 4.17374e+1, 5.77217e+1,
      7.83868e+1, 2.08677e+2, 4.04331e+2, 8.01464e+2, 1.72802e+3, 3.00284e+3,
      6.75852e+3, 1.24327e+4, 2.03901e+4, 3.09829e+4, 4.45725e+4, 6.15554e+4,
      8.23790e+4, 1.07544e+5};
  pf[17] = detail::bc16_interp(qC2, T_K);

  // pf[18]: CH (BC16)
  static constexpr std::array<double, 32> qch = {
      1.20000e+1, 1.20000e+1, 1.20000e+1, 1.20000e+1, 1.20000e+1, 1.20000e+1,
      1.20002e+1, 1.20057e+1, 1.20868e+1, 1.23391e+1, 1.33487e+1, 1.64087e+1,
      1.99589e+1, 2.55718e+1, 3.13181e+1, 3.90674e+1, 4.49138e+1, 5.46956e+1,
      6.45085e+1, 1.03985e+2, 1.44136e+2, 2.07690e+2, 3.30136e+2, 4.82164e+2,
      9.04412e+2, 1.52375e+3, 2.38711e+3, 3.54843e+3, 5.06630e+3, 6.99662e+3,
      9.38682e+3, 1.22732e+4};
  pf[18] = detail::bc16_interp(qch, T_K);

  // pf[19]: CH2 (CDMS rotational for T<=600, I88 polynomial for T>600)
  if (T_K <= 600.0) {
    pf[19] = asym_top(2211494.0, 253618.1, 215102.2, 1.5);
  } else {
    pf[19] = std::pow(10.0, 5.05057221 - 3.7971851 * lg_th +
                                0.928883 * lg_th2 + 0.44249 * lg_th3 +
                                0.14807 * lg_th4 - 0.37003 * lg_th5);
  }

  // pf[20]: CH3 (HITRAN file / 8 for T<=500, I88 for T>500)
  if (T_K <= 500.0) {
    pf[20] = detail::pf_interp(tbl.pf_table[20], T_K) / 8.0;
  } else {
    pf[20] = std::pow(10.0, 6.2439699776 - 5.909173844 * lg_th +
                                1.5818662 * lg_th2 + 0.767593 * lg_th3 +
                                0.449532 * lg_th4 - 1.04287 * lg_th5 +
                                0.264771 * lg_th6);
  }

  // pf[21]: CH4 (HITRAN file / 16 for T<=3500, I88 for T>3500)
  if (T_K <= 3500.0) {
    pf[21] = detail::pf_interp(tbl.pf_table[21], T_K) / 16.0;
  } else {
    pf[21] =
        std::pow(10.0, 6.942294193 - 8.66259486 * lg_th + 3.218163 * lg_th2 +
                           0.69805 * lg_th3 + 0.92444 * lg_th4 -
                           1.3878 * lg_th5 + 0.1889 * lg_th6);
  }

  // pf[22]: C+ (BC16)
  static constexpr std::array<double, 32> qcp = {
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00001e+0, 2.00044e+0, 2.00913e+0, 2.04177e+0, 2.19111e+0, 2.64504e+0,
      3.08645e+0, 3.60629e+0, 3.98273e+0, 4.33873e+0, 4.53479e+0, 4.77694e+0,
      4.95108e+0, 5.33283e+0, 5.51120e+0, 5.65121e+0, 5.76395e+0, 5.82163e+0,
      5.88018e+0, 5.90980e+0, 5.92772e+0, 5.94003e+0, 5.94994e+0, 5.95988e+0,
      5.97207e+0, 5.98845e+0};
  pf[22] = detail::bc16_interp(qcp, T_K);

  // pf[23]: C2+ (BC16)
  static constexpr std::array<double, 32> qc2p = {
      6.00071e+0, 6.00646e+0, 6.03640e+0, 6.08459e+0, 6.41978e+0, 7.60333e+0,
      9.05993e+0, 1.14110e+1, 1.54767e+1, 1.96058e+1, 2.79266e+1, 4.46440e+1,
      6.13923e+1, 8.65365e+1, 1.11695e+2, 1.45258e+2, 1.70449e+2, 2.12517e+2,
      2.54831e+2, 4.31572e+2, 6.30643e+2, 9.86295e+2, 1.74817e+3, 2.72947e+3,
      5.36363e+3, 8.90825e+3, 1.33809e+4, 1.88008e+4, 2.51901e+4, 3.25762e+4,
      4.09927e+4, 5.04806e+4};
  pf[23] = detail::bc16_interp(qc2p, T_K);

  // pf[24]: CH+ (BC16)
  static constexpr std::array<double, 32> qchp = {
      1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00099e+0,
      1.00977e+0, 1.05449e+0, 1.20886e+0, 1.41647e+0, 1.88153e+0, 2.85679e+0,
      3.84577e+0, 5.33613e+0, 6.82949e+0, 8.82249e+0, 1.03179e+1, 1.28110e+1,
      1.53046e+1, 2.52948e+1, 3.54289e+1, 5.14755e+1, 8.22832e+1, 1.20289e+2,
      2.30672e+2, 4.15819e+2, 7.12556e+2, 1.15545e+3, 1.77578e+3, 2.60123e+3,
      3.65584e+3, 4.96011e+3};
  pf[24] = detail::bc16_interp(qchp, T_K);

  // pf[25]: CH2+ (asymmetric top)
  pf[25] = asym_top(2075377.0, 236712.0, 212478.0);

  // pf[26]: CH3+ (symmetric top / 6)
  pf[26] = asym_top(279591.0, 279591.0, 139796.0, 1.0 / 6.0);

  // pf[27]: CH4+ (asymmetric top / 2)
  pf[27] = asym_top(199364.0, 154791.0, 113699.0, 0.5);

  // pf[28]: CH5+ (asymmetric top)
  pf[28] = asym_top(113321.5, 113921.1, 114820.5);

  // pf[29]: O (BC16)
  static constexpr std::array<double, 32> qo = {
      5.00000e+0, 5.00000e+0, 5.00000e+0, 5.00000e+0, 5.00000e+0, 5.00000e+0,
      5.00000e+0, 5.00000e+0, 5.00000e+0, 5.00003e+0, 5.00154e+0, 5.03305e+0,
      5.12544e+0, 5.34602e+0, 5.60172e+0, 5.93258e+0, 6.15640e+0, 6.47757e+0,
      6.74123e+0, 7.42310e+0, 7.79423e+0, 8.11056e+0, 8.38189e+0, 8.52662e+0,
      8.68009e+0, 8.77224e+0, 8.85532e+0, 8.94697e+0, 9.05112e+0, 9.16637e+0,
      9.28982e+0, 9.41864e+0};
  pf[29] = detail::bc16_interp(qo, T_K);

  // pf[30]: O2 (BC16)
  static constexpr std::array<double, 32> qo2 = {
      9.00000e+0, 9.00000e+0, 9.00011e+0, 9.00068e+0, 9.02128e+0, 9.33579e+0,
      1.01023e+1, 1.17557e+1, 1.50103e+1, 1.84654e+1, 2.55565e+1, 3.99407e+1,
      5.43989e+1, 7.61289e+1, 9.78808e+1, 1.26901e+2, 1.48678e+2, 1.85004e+2,
      2.21427e+2, 3.70966e+2, 5.34408e+2, 8.19779e+2, 1.42263e+3, 2.19951e+3,
      4.34427e+3, 7.41870e+3, 1.16070e+4, 1.71056e+4, 2.41444e+4, 3.30098e+4,
      4.40540e+4, 5.76869e+4};
  pf[30] = detail::bc16_interp(qo2, T_K);

  // pf[31]: OH (BC16)
  static constexpr std::array<double, 32> qoh = {
      1.20000e+1, 1.20000e+1, 1.20000e+1, 1.20000e+1, 1.20000e+1, 1.20000e+1,
      1.20000e+1, 1.20005e+1, 1.20163e+1, 1.20966e+1, 1.25753e+1, 1.45073e+1,
      1.70195e+1, 2.11581e+1, 2.54704e+1, 3.13320e+1, 3.57710e+1, 4.32133e+1,
      5.06892e+1, 8.07652e+1, 1.11070e+2, 1.57530e+2, 2.41294e+2, 3.37500e+2,
      5.77700e+2, 8.91168e+2, 1.28807e+3, 1.78078e+3, 2.38461e+3, 3.11647e+3,
      3.99262e+3, 5.02698e+3};
  pf[31] = detail::bc16_interp(qoh, T_K);

  // pf[32]: CO (BC16)
  static constexpr std::array<double, 32> qco = {
      1.01186e+0, 1.04254e+0, 1.11606e+0, 1.18989e+0, 1.49427e+0, 2.18214e+0,
      2.89218e+0, 3.96754e+0, 5.76822e+0, 7.57241e+0, 1.11842e+1, 1.84123e+1,
      2.56425e+1, 3.64899e+1, 4.73391e+1, 6.18073e+1, 7.26603e+1, 9.07524e+1,
      1.08852e+2, 1.81659e+2, 2.56937e+2, 3.80242e+2, 6.25544e+2, 9.28195e+2,
      1.71706e+3, 2.75995e+3, 4.06593e+3, 5.64445e+3, 7.50751e+3, 9.67381e+3,
      1.21758e+4, 1.50689e+4};
  pf[32] = detail::bc16_interp(qco, T_K);

  // pf[33]: H2O (Vidler & Tennyson 2000 polynomial in log10(T))
  {
    double lgT = std::log10(T_K);
    double lgT2 = lgT * lgT, lgT3 = lgT2 * lgT, lgT4 = lgT3 * lgT,
           lgT5 = lgT4 * lgT, lgT6 = lgT5 * lgT;
    double lgq = -14.0874691574179 + 37.9243248539882 * lgT -
                 42.6817978731789 * lgT2 + 25.3302448517916 * lgT3 -
                 8.10851262935532 * lgT4 + 1.33106871720535 * lgT5 -
                 0.0872981051095757 * lgT6;
    pf[33] = std::max(0.25, std::pow(10.0, lgq));
  }

  // pf[34]: HCO (JPL rotational formula for T<=500, I88 for T>500)
  if (T_K <= 500.0) {
    pf[34] = asym_top(670485.83, 44788.0, 41930.4, 2.0);
  } else {
    pf[34] =
        std::pow(10.0, 6.298781639 - 3.85672804 * lg_th + 0.8551678 * lg_th2 +
                           0.321901 * lg_th3 + 0.020274 * lg_th4 +
                           0.15254 * lg_th5 - 0.25298 * lg_th6);
  }

  // pf[35]: HO2 (HITRAN file / 2)
  pf[35] = detail::pf_interp(tbl.pf_table[35], T_K) / 2.0;

  // pf[36]: CO2 (HITRAN file for T<=5000, I88 for T>5000)
  if (T_K <= 5000.0) {
    pf[36] = detail::pf_interp(tbl.pf_table[36], T_K);
  } else {
    pf[36] = std::pow(10.0, 6.01081285 - 4.438833 * lg_th +
                                0.840462 * lg_th2 + 0.2945 * lg_th3 +
                                0.3694 * lg_th4 - 0.273 * lg_th5);
  }

  // pf[37]: H2CO (HITRAN file / 4)
  pf[37] = detail::pf_interp(tbl.pf_table[37], T_K) / 4.0;

  // pf[38]: H2O2 (HITRAN file / 4)
  pf[38] = detail::pf_interp(tbl.pf_table[38], T_K) / 4.0;

  // pf[39]: O+ (BC16)
  static constexpr std::array<double, 32> qop = {
      4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0,
      4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0,
      4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0,
      4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0, 4.00000e+0,
      4.00003e+0, 4.00065e+0, 4.00451e+0, 4.01649e+0, 4.04186e+0, 4.08460e+0,
      4.14679e+0, 4.22885e+0};
  pf[39] = detail::bc16_interp(qop, T_K);

  // pf[40]: O2+ (BC16)
  static constexpr std::array<double, 32> qo2p = {
      3.00032e+0, 3.00294e+0, 3.01694e+0, 3.03979e+0, 3.20139e+0, 3.78014e+0,
      4.49647e+0, 5.65503e+0, 7.66052e+0, 9.69811e+0, 1.38059e+1, 2.21326e+1,
      3.08526e+1, 4.52466e+1, 6.13856e+1, 8.52491e+1, 1.04528e+2, 1.38583e+2,
      1.74381e+2, 3.28202e+2, 4.96506e+2, 7.83125e+2, 1.37289e+3, 2.11772e+3,
      4.10188e+3, 6.77817e+3, 1.01888e+4, 1.44008e+4, 1.95417e+4, 2.58386e+4,
      3.36476e+4, 4.34606e+4};
  pf[40] = detail::bc16_interp(qo2p, T_K);

  // pf[41]: OH+ (BC16)
  static constexpr std::array<double, 32> qohp = {
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00047e+0,
      2.00703e+0, 2.05326e+0, 2.25799e+0, 2.57365e+0, 3.33233e+0, 4.97055e+0,
      6.64415e+0, 9.17264e+0, 1.17099e+1, 1.50998e+1, 1.76456e+1, 2.18936e+1,
      2.61470e+1, 4.32173e+1, 6.04874e+1, 8.74853e+1, 1.38201e+2, 1.98838e+2,
      3.56403e+2, 5.71245e+2, 8.57404e+2, 1.23389e+3, 1.72209e+3, 2.34214e+3,
      3.11051e+3, 4.03909e+3};
  pf[41] = detail::bc16_interp(qohp, T_K);

  // pf[42]: CO+ (BC16)
  static constexpr std::array<double, 32> qcop = {
      2.02087e+0, 2.07708e+0, 2.21516e+0, 2.35588e+0, 2.94397e+0, 4.28449e+0,
      5.67104e+0, 7.77209e+0, 1.12911e+1, 1.48173e+1, 2.18769e+1, 3.60049e+1,
      5.01372e+1, 7.13398e+1, 9.25460e+1, 1.20826e+2, 1.42039e+2, 1.77403e+2,
      2.12779e+2, 3.55006e+2, 5.01729e+2, 7.41213e+2, 1.21592e+3, 1.80037e+3,
      3.32256e+3, 5.34428e+3, 7.92672e+3, 1.11813e+4, 1.52681e+4, 2.03854e+4,
      2.67626e+4, 3.46562e+4};
  pf[42] = detail::bc16_interp(qcop, T_K);

  // pf[43]: H2O+ (asymmetric top)
  pf[43] = asym_top(870581.0, 372365.0, 253880.0);

  // pf[44]: HCO+ (linear rotor for T<=350, I88 for T>350)
  if (T_K <= 350.0) {
    double B = 44594.43 * MHz;
    pf[44] = (k_B * T_K / h_P) / B;
  } else {
    pf[44] = std::pow(10.0, 5.453394934 - 4.14568927 * lg_th +
                                0.8632023 * lg_th2 + 0.482875 * lg_th3 +
                                0.14863 * lg_th4 - 0.281734 * lg_th5);
  }

  // pf[45]: O2H+ (asymmetric top × 3)
  pf[45] = asym_top(663778.0, 38723.0, 36588.0, 3.0);

  // pf[46]: H3O+ (symmetric top / 3 for T<=500, I88 for T>500)
  if (T_K <= 500.0) {
    pf[46] = asym_top(334405.0, 334405.0, 184725.0, 1.0 / 3.0);
  } else {
    pf[46] = std::pow(10.0, 5.6188225644 - 5.599410492 * lg_th +
                                1.8194835 * lg_th2 + 0.774831 * lg_th3 +
                                0.404714 * lg_th4 - 1.53688 * lg_th5 +
                                0.655507 * lg_th6);
  }

  // pf[47]: H2CO+ (asymmetric top)
  pf[47] = asym_top(266035.83, 40232.15, 34416.17);

  // pf[48]: HOCO+ (asymmetric top)
  pf[48] = asym_top(789947.8, 10773.733, 10609.431);

  // pf[49]: H2COH+ (asymmetric top)
  pf[49] = asym_top(197581.56, 34350.55, 29172.65);

  // pf[50]: Li (BC16)
  static constexpr std::array<double, 32> qLi_mg = {
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00013e+0,
      2.00473e+0, 2.02869e+0, 2.08800e+0, 2.19942e+0, 2.39054e+0, 2.70188e+0,
      3.17983e+0, 3.86752e+0};
  pf[50] = detail::bc16_interp(qLi_mg, T_K);

  // pf[51]: LiH (BC16)
  static constexpr std::array<double, 32> qLiH_mg = {
      1.00000e+0, 1.00000e+0, 1.00001e+0, 1.00007e+0, 1.00247e+0, 1.04234e+0,
      1.14353e+0, 1.36475e+0, 1.79698e+0, 2.25093e+0, 3.17620e+0, 5.04593e+0,
      6.92359e+0, 9.74658e+0, 1.25749e+1, 1.63530e+1, 1.91923e+1, 2.39405e+1,
      2.87248e+1, 4.88043e+1, 7.16412e+1, 1.13107e+2, 2.04843e+2, 3.28787e+2,
      6.95286e+2, 1.26125e+3, 2.06412e+3, 3.13695e+3, 4.53795e+3, 6.36625e+3,
      8.75776e+3, 1.18704e+4};
  pf[51] = detail::bc16_interp(qLiH_mg, T_K);

  // pf[52]: Li+  = 1.0 (fill default)
  // pf[53]: Li-  = 1.0 (fill default)

  // pf[54]: LiH+ (ExoMol file, nuclear-spin deg removed: /2)
  pf[54] = detail::pf_interp(tbl.pf_table[54], T_K) / 2.0;

  // pf[55]: Li++
  pf[55] = 2.0;

  // pf[56]: Li+++ = 1.0 (fill default)

  // pf[57]: K (BC16)
  static constexpr std::array<double, 32> qK = {
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00002e+0, 2.00052e+0,
      2.01222e+0, 2.06702e+0, 2.22592e+0, 2.61039e+0, 3.39777e+0, 4.77353e+0,
      6.88524e+0, 9.82105e+0};
  pf[57] = detail::bc16_interp(qK, T_K);

  // pf[58]: K+  = 1.0 (fill default)

  // pf[59]: Na (BC16)
  static constexpr std::array<double, 32> qNa = {
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00003e+0,
      2.00178e+0, 2.01488e+0, 2.06447e+0, 2.21609e+0, 2.60177e+0, 3.40984e+0,
      4.84529e+0, 7.08960e+0};
  pf[59] = detail::bc16_interp(qNa, T_K);

  // pf[60]: Na+ = 1.0 (fill default)

  // pf[61]: Mg (BC16)
  static constexpr std::array<double, 32> qMg = {
      1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0,
      1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0,
      1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0,
      1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0, 1.00000e+0,
      1.00025e+0, 1.00344e+0, 1.01678e+0, 1.04919e+0, 1.11004e+0, 1.21285e+0,
      1.37994e+0, 1.64434e+0};
  pf[61] = detail::bc16_interp(qMg, T_K);

  // pf[62]: Mg+ (BC16)
  static constexpr std::array<double, 32> qMgp = {
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0, 2.00000e+0,
      2.00000e+0, 2.00002e+0, 2.00021e+0, 2.00114e+0, 2.00389e+0, 2.00976e+0,
      2.02002e+0, 2.03571e+0};
  pf[62] = detail::bc16_interp(qMgp, T_K);

  // 63..88: grain/ice pseudo-species = 1.0 (fill default)
}

}  // namespace arche
