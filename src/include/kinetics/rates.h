// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
#include <algorithm>
#include <array>
#include <cmath>

#include "core/state.h"
#include "kinetics/partition_function.h"
#include "kinetics/topology.h"
#include "kinetics/reaction_index.h"

// ---------------------------------------------------------------------------
// Shared reaction-rate primitives (model-independent).
//
// Provides:
//   pow_int(x, n)              — integer-exponent power for small exponents.
//   detail::eval_opacity(...)  — continuum opacity k_prm [cm²/g].
//   PrimRateContext            — species-independent intermediate quantities
//                                (temperature powers, LTE critical densities,
//                                LTE/low-density blend sub-rates).
//   build_prim_context(...)    — fills a PrimRateContext from the gas state.
//   lte_blend / lte_blend_noguard — LTE/low-density blend of two rates.
//   compute_reverse_loop(...)  — detailed-balance reverse rates (shared).
//   make_y_ext(...)            — extended species vector with vacant/photon/CR
//                                sentinels consumed by the reaction loops.
//
// All rate coefficients in CGS (cm³/s or cm⁶/s for 3-body). k_rxn indexing:
// k_rxn[num-1] = forward, k_rxn[num-1+N_react] = reverse (1-based num).
// ---------------------------------------------------------------------------

namespace arche {

// Integer-exponent power for small non-negative exponents.
inline double pow_int(double x, int n) {
  switch (n) {
    case 0:
      return 1.0;
    case 1:
      return x;
    case 2:
      return x * x;
    case 3:
      return x * x * x;
    default:
      return std::pow(x, static_cast<double>(n));
  }
}

namespace detail {

// ---------------------------------------------------------------------------
// eval_opacity — bilinear interpolation of continuum opacity table k_prm
// ---------------------------------------------------------------------------
inline double eval_opacity(double T_K, double rho) {
  static constexpr std::array<double, 57> T6 = {
      0.75e-3, 1.00e-3, 1.25e-3, 1.50e-3, 1.75e-3, 2.00e-3, 2.25e-3, 2.50e-3,
      2.75e-3, 3.00e-3, 3.25e-3, 3.50e-3, 3.75e-3, 4.00e-3, 4.25e-3, 4.50e-3,
      4.75e-3, 5.00e-3, 5.25e-3, 5.50e-3, 5.75e-3, 6.00e-3, 6.25e-3, 6.50e-3,
      6.75e-3, 7.00e-3, 8.00e-3, 9.00e-3, 1.00e-2, 1.10e-2, 1.20e-2, 1.40e-2,
      1.60e-2, 1.80e-2, 2.00e-2, 2.50e-2, 3.00e-2, 3.50e-2, 4.00e-2, 4.50e-2,
      5.00e-2, 5.50e-2, 6.00e-2, 7.00e-2, 8.00e-2, 9.00e-2, 1.00e-1, 1.20e-1,
      1.50e-1, 2.00e-1, 2.50e-1, 3.00e-1, 4.00e-1, 5.00e-1, 6.00e-1, 8.00e-1,
      1.00e+0};

  static constexpr std::array<double, 19> rlgR = {
      -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5,
      0.0,  0.5,  1.0,  1.5,  2.0,  2.5,  3.0,  3.5,  4.0};

  static constexpr std::array<std::array<double, 19>, 57> rlgk = {
      {{-12.43, -11.91, -11.39, -10.86, -10.36, -9.804, -9.350, -8.847, -8.346,
        -7.845, -7.345, -6.845, -6.344, -5.843, -5.343, -4.843, -4.343, -3.843,
        -3.343},
       {-11.85, -11.35, -10.85, -10.35, -9.854, -9.354, -8.854, -8.354, -7.854,
        -7.354, -6.854, -6.354, -5.854, -5.354, -4.854, -4.354, -3.854, -3.354,
        -2.854},
       {-11.46, -10.96, -10.46, -9.964, -9.464, -8.964, -8.464, -7.964, -7.464,
        -6.964, -6.464, -5.964, -5.464, -4.964, -4.464, -3.964, -3.464, -2.964,
        -2.464},
       {-11.26, -10.74, -10.22, -9.704, -9.194, -8.634, -8.180, -7.677, -7.176,
        -6.675, -6.175, -5.675, -5.174, -4.673, -4.173, -3.673, -3.173, -2.673,
        -2.173},
       {-11.86, -11.12, -10.42, -9.757, -9.139, -8.557, -8.018, -7.491, -6.979,
        -6.471, -5.967, -5.465, -4.962, -4.461, -3.961, -3.461, -2.961, -2.461,
        -1.961},
       {-13.67, -12.51, -11.39, -10.43, -9.514, -8.757, -8.039, -7.433, -6.854,
        -6.318, -5.792, -5.281, -4.772, -4.269, -3.767, -3.266, -2.765, -2.264,
        -1.764},
       {-12.58, -12.14, -11.67, -10.99, -10.27, -9.348, -8.440, -7.672, -6.923,
        -6.322, -5.727, -5.193, -4.661, -4.150, -3.640, -3.137, -2.634, -2.133,
        -1.632},
       {-10.62, -10.50, -10.32, -10.07, -9.761, -9.375, -8.827, -8.025, -7.238,
        -6.478, -5.780, -5.178, -4.602, -4.067, -3.542, -3.032, -2.524, -2.020,
        -1.518},
       {-9.597, -9.352, -9.106, -8.857, -8.606, -8.354, -8.075, -7.781, -7.270,
        -6.635, -5.954, -5.245, -4.617, -4.036, -3.489, -2.961, -2.446, -1.937,
        -1.432},
       {-8.520, -8.331, -8.108, -7.879, -7.615, -7.345, -7.093, -6.844, -6.582,
        -6.318, -5.827, -5.300, -4.677, -4.038, -3.472, -2.918, -2.394, -1.875,
        -1.367},
       {-7.730, -7.481, -7.232, -6.982, -6.732, -6.483, -6.233, -5.984, -5.741,
        -5.498, -5.264, -5.016, -4.573, -4.121, -3.559, -3.000, -2.473, -1.947,
        -1.436},
       {-7.011, -6.747, -6.498, -6.255, -6.027, -5.788, -5.516, -5.249, -5.001,
        -4.758, -4.531, -4.309, -4.101, -3.832, -3.395, -2.937, -2.414, -1.894,
        -1.383},
       {-6.359, -6.111, -5.862, -5.612, -5.362, -5.112, -4.862, -4.613, -4.364,
        -4.119, -3.879, -3.661, -3.474, -3.272, -3.051, -2.732, -2.291, -1.825,
        -1.328},
       {-5.798, -5.552, -5.305, -5.057, -4.809, -4.560, -4.312, -4.061, -3.808,
        -3.560, -3.314, -3.089, -2.876, -2.704, -2.559, -2.329, -2.046, -1.667,
        -1.228},
       {-5.288, -5.048, -4.807, -4.563, -4.317, -4.068, -3.819, -3.569, -3.319,
        -3.070, -2.822, -2.587, -2.357, -2.183, -2.027, -1.858, -1.684, -1.385,
        -1.048},
       {-4.810, -4.585, -4.359, -4.119, -3.877, -3.630, -3.382, -3.133, -2.883,
        -2.634, -2.385, -2.145, -1.905, -1.718, -1.536, -1.400, -1.267, -1.034,
        -0.792},
       {-4.344, -4.143, -3.941, -3.710, -3.479, -3.235, -2.991, -2.742, -2.494,
        -2.244, -1.995, -1.751, -1.509, -1.304, -1.104, -0.975, -0.843, -0.661,
        -0.472},
       {-3.875, -3.713, -3.541, -3.328, -3.111, -2.875, -2.636, -2.390, -2.143,
        -1.893, -1.645, -1.399, -1.157, -0.937, -0.732, -0.590, -0.448, -0.300,
        -0.133},
       {-3.410, -3.288, -3.147, -2.961, -2.764, -2.539, -2.309, -2.067, -1.824,
        -1.576, -1.329, -1.081, -0.839, -0.608, -0.400, -0.240, -0.091, 0.037,
        0.191},
       {-2.958, -2.872, -2.757, -2.606, -2.430, -2.224, -2.005, -1.771, -1.532,
        -1.286, -1.039, -0.792, -0.549, -0.312, -0.101, 0.077, 0.229, 0.347,
        0.496},
       {-2.529, -2.467, -2.374, -2.257, -2.103, -1.922, -1.717, -1.493, -1.259,
        -1.017, -0.773, -0.527, -0.283, -0.042, 0.173, 0.366, 0.516, 0.632,
        0.780},
       {-2.130, -2.080, -2.004, -1.916, -1.784, -1.631, -1.441, -1.234, -1.007,
        -0.771, -0.529, -0.284, -0.041, 0.201, 0.420, 0.628, 0.776, 0.896,
        1.042},
       {-1.766, -1.717, -1.653, -1.585, -1.473, -1.349, -1.174, -0.986, -0.768,
        -0.542, -0.303, -0.060, 0.183, 0.425, 0.649, 0.866, 1.014, 1.142,
        1.285},
       {-1.437, -1.380, -1.324, -1.269, -1.174, -1.074, -0.916, -0.751, -0.542,
        -0.327, -0.091, 0.147, 0.390, 0.633, 0.859, 1.083, 1.230, 1.367, 1.507},
       {-1.149, -1.070, -1.019, -0.968, -0.887, -0.806, -0.665, -0.523, -0.324,
        -0.123, 0.108, 0.341, 0.583, 0.825, 1.055, 1.284, 1.433, 1.580, 1.715},
       {-0.900, -0.792, -0.737, -0.681, -0.615, -0.544, -0.424, -0.299, -0.115,
        0.073, 0.298, 0.523, 0.763, 1.002, 1.235, 1.463, 1.620, 1.775, 1.907},
       {-0.224, -0.120, -0.033, 0.037, 0.098, 0.189, 0.286, 0.402, 0.538, 0.684,
        0.842, 1.019, 1.213, 1.405, 1.626, 1.850, 2.077, 2.307, 2.540},
       {-0.096, 0.191, 0.434, 0.632, 0.785, 0.893, 0.987, 1.084, 1.181, 1.290,
        1.410, 1.556, 1.721, 1.875, 2.080, 2.301, 2.510, 2.679, 2.780},
       {-0.107, 0.234, 0.579, 0.929, 1.196, 1.392, 1.527, 1.624, 1.708, 1.799,
        1.901, 2.020, 2.156, 2.270, 2.458, 2.664, 2.861, 3.022, 3.120},
       {-0.145, 0.193, 0.572, 0.991, 1.388, 1.681, 1.885, 2.026, 2.118, 2.202,
        2.294, 2.391, 2.521, 2.606, 2.778, 2.973, 3.160, 3.308, 3.386},
       {-0.180, 0.139, 0.519, 0.960, 1.415, 1.810, 2.106, 2.297, 2.437, 2.537,
        2.632, 2.733, 2.848, 2.905, 3.061, 3.237, 3.422, 3.605, 3.775},
       {-0.205, 0.048, 0.391, 0.825, 1.314, 1.806, 2.241, 2.585, 2.837, 3.012,
        3.139, 3.259, 3.378, 3.386, 3.527, 3.684, 3.853, 4.030, 4.211},
       {-0.229, 0.018, 0.341, 0.741, 1.205, 1.726, 2.218, 2.671, 3.041, 3.307,
        3.499, 3.651, 3.786, 3.762, 3.896, 4.045, 4.197, 4.340, 4.462},
       {-0.251, -0.024, 0.296, 0.709, 1.186, 1.679, 2.175, 2.681, 3.135, 3.490,
        3.755, 3.949, 4.110, 4.076, 4.204, 4.354, 4.495, 4.596, 4.626},
       {-0.258, -0.054, 0.260, 0.684, 1.160, 1.657, 2.170, 2.685, 3.182, 3.612,
        3.939, 4.180, 4.375, 4.351, 4.465, 4.590, 4.727, 4.877, 5.041},
       {-0.229, -0.046, 0.255, 0.676, 1.166, 1.695, 2.231, 2.765, 3.300, 3.808,
        4.238, 4.582, 4.836, 4.831, 4.958, 5.046, 5.097, 5.113, 5.096},
       {-0.168, 0.014, 0.318, 0.743, 1.245, 1.772, 2.317, 2.869, 3.440, 3.971,
        4.442, 4.832, 5.138, 5.101, 5.206, 5.206, 5.206, 5.206, 5.206},
       {-0.078, 0.124, 0.432, 0.845, 1.341, 1.874, 2.420, 2.977, 3.552, 4.100,
        4.589, 4.995, 5.312, 5.250, 5.426, 5.426, 5.426, 5.426, 5.426},
       {-0.064, 0.221, 0.571, 0.985, 1.461, 1.980, 2.516, 3.068, 3.631, 4.175,
        4.666, 5.073, 5.384, 5.392, 5.555, 5.555, 5.555, 5.555, 5.555},
       {-0.101, 0.211, 0.606, 1.084, 1.565, 2.065, 2.583, 3.120, 3.664, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.138, 0.149, 0.541, 1.038, 1.568, 2.091, 2.604, 3.125, 3.649, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.170, 0.099, 0.469, 0.941, 1.476, 2.030, 2.570, 3.091, 3.603, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.198, 0.059, 0.409, 0.854, 1.376, 1.925, 2.489, 3.032, 3.540, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.230, -0.011, 0.311, 0.736, 1.242, 1.761, 2.329, 2.895, 3.431, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.242, -0.044, 0.260, 0.670, 1.169, 1.685, 2.239, 2.799, 3.354, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.240, -0.055, 0.241, 0.649, 1.140, 1.655, 2.202, 2.758, 3.316, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.240, -0.050, 0.243, 0.639, 1.128, 1.652, 2.196, 2.749, 3.305, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.246, -0.045, 0.258, 0.661, 1.145, 1.672, 2.210, 2.754, 3.292, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.286, -0.078, 0.229, 0.634, 1.112, 1.621, 2.135, 2.646, 3.141, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.350, -0.196, 0.054, 0.400, 0.818, 1.280, 1.758, 2.236, 2.700, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.387, -0.281, -0.090, 0.184, 0.539, 0.957, 1.413, 1.881, 2.340, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.409, -0.332, -0.182, 0.041, 0.347, 0.729, 1.163, 1.620, 2.077, 0.000,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.431, -0.384, -0.280, -0.118, 0.123, 0.446, 0.839, 1.275, 1.724,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.442, -0.409, -0.329, -0.200, 0.000, 0.282, 0.641, 1.055, 1.493,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.449, -0.424, -0.357, -0.249, -0.074, 0.179, 0.511, 0.905, 1.333,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.456, -0.439, -0.389, -0.305, -0.164, 0.049, 0.341, 0.704, 1.115,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
       {-0.460, -0.447, -0.405, -0.335, -0.213, -0.026, 0.238, 0.577, 0.972,
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000}}};

  double temp6 = T_K / 1.0e6;
  double rlogR = std::log10(rho / (temp6 * temp6 * temp6));

  bool in_T_lo = (temp6 >= 0.75e-3) && (temp6 <= 2.5e-2);
  bool in_T_mid = (temp6 >= 2.5e-2) && (temp6 <= 4.0e-2);
  bool in_T_hi = (temp6 > 4.0e-2);
  bool in_R_full = (rlogR >= -5.0) && (rlogR <= 4.0);
  bool in_R_mid = (rlogR >= -5.0) && (rlogR <= 3.0);
  bool in_R_hi = (rlogR >= -5.0) && (rlogR <= -1.0);

  bool do_interp =
      (in_T_lo && in_R_full) || (in_T_mid && in_R_mid) || (in_T_hi && in_R_hi);

  if (temp6 < 0.75e-3 || rlogR < -5.0) return 0.0;
  if (!do_interp) return 1.0e5;

  int ms = 56;
  for (int i = 0; i < 57; ++i)
    if (temp6 <= T6[i]) {
      ms = i;
      break;
    }
  if (ms == 0) ms = 1;

  int ns = 18;
  for (int i = 0; i < 19; ++i)
    if (rlogR <= rlgR[i]) {
      ns = i;
      break;
    }
  if (ns == 0) ns = 1;

  double t = (temp6 - T6[ms - 1]) / (T6[ms] - T6[ms - 1]);
  double u = (rlogR - rlgR[ns - 1]) / (rlgR[ns] - rlgR[ns - 1]);
  double rlogk = (1.0 - t) * (1.0 - u) * rlgk[ms - 1][ns - 1] +
                 t * (1.0 - u) * rlgk[ms][ns - 1] + t * u * rlgk[ms][ns] +
                 (1.0 - t) * u * rlgk[ms - 1][ns];
  return std::pow(10.0, rlogk);
}

}  // namespace detail

// ---------------------------------------------------------------------------
// PrimRateContext — intermediate quantities that depend only on temperature,
// density, and the H/H2/He background abundances (not on the per-species
// abundance vector).  Built once per rate evaluation and read by every forward
// rate law and by the reverse / LTE-correction step, so the heavy shared
// pieces (LTE critical densities, the LTE/low-density blend sub-rates) are
// computed a single time.  Sharing these through one context keeps the rate
// laws free of global reaction-id coupling while preserving the exact original
// operation order (bit-for-bit).
// ---------------------------------------------------------------------------
struct PrimRateContext {
  double T_K, nH;
  // Temperature powers used across the network.
  double T_eV, lnT_eV, T300, lnT;
  double lgT, lgT2, lgT3, lgT4, lgT5, lgT6, lT4;
  // LTE critical-density ratios for the collisional H2 (de)excitation blends.
  double ncr, crit_ratio, crit_ratio_HD;
  // LTE/low-density blend sub-rates: each pairs a low-density (v0) and an LTE
  // coefficient consumed by one forward law and its detailed-balance reverse.
  double k12v0, k12LTE;  // H2 + e  -> 2H + e
  double k13v0, k13LTE;  // H2 + H  -> 3H
  double k21v0, k21LTE;  // 2H2     -> 2H + H2
  double k37v0, k37LTE;  // H2 + He -> 2H + He
  double k40v0;          // H2 + e  -> H- + H
  double k47v0, k47LTE;  // HD + e  -> H + D + e
  double k75v0;          // HD + e  -> H + D-
  double delE_22;        // |ΔE| source for k_rxn[21]
};

// ---------------------------------------------------------------------------
// build_prim_context — compute every PrimRateContext field from the gas state.
//
//   T_K     — gas temperature [K].
//   nH      — hydrogen number density [cm^-3].
//   H,H2,He — background abundances relative to nH (drive the ncr mix).
//   delE_22 — |ΔE| for k_rxn[21] (tbl.reactions[21].delE) [erg].
// ---------------------------------------------------------------------------
inline PrimRateContext build_prim_context(double T_K, double nH, double H,
                                          double H2, double He,
                                          double delE_22) {
  PrimRateContext c;
  c.T_K = T_K;
  c.nH = nH;

  c.T_eV = phys::k_B_eV * T_K;
  c.lnT_eV = std::log(c.T_eV);
  c.T300 = T_K / 300.0;
  c.lnT = std::log(T_K);
  c.lgT = std::log10(T_K);
  c.lgT2 = c.lgT * c.lgT;
  c.lgT3 = c.lgT2 * c.lgT;
  c.lgT4 = c.lgT3 * c.lgT;
  c.lgT5 = c.lgT4 * c.lgT;
  c.lgT6 = c.lgT5 * c.lgT;
  c.lT4 = std::log10(T_K / 1.0e4);

  // Critical densities for LTE interpolation (GA08)
  double ncr_H = std::pow(10.0, 3.0 - 0.416 * c.lT4 - 0.327 * c.lT4 * c.lT4);
  double ncr_H2 = std::pow(10.0, 4.845 - 1.3 * c.lT4 + 1.62 * c.lT4 * c.lT4);
  double ncr_He = std::pow(10.0, 5.0792 * (1.0 - 1.23e-5 * (T_K - 2000.0)));
  c.ncr = 1.0 / (H / ncr_H + H2 / ncr_H2 + He / ncr_He);
  c.crit_ratio = nH / c.ncr;
  double ncr_HD = 1.0e2 * c.ncr;
  c.crit_ratio_HD = nH / ncr_HD;

  // LTE/low-density blend sub-rates
  c.k12v0 = 4.49e-9 * std::pow(T_K, 0.11) * std::exp(-101858.0 / T_K);
  c.k12LTE = 1.91e-9 * std::pow(T_K, 0.136) * std::exp(-53407.1 / T_K);
  c.k13v0 = 6.67e-12 * std::sqrt(T_K) * std::exp(-(1.0 + 63593.0 / T_K));
  c.k13LTE = 3.52e-9 * std::exp(-43900.0 / T_K);
  c.k21v0 = (5.996e-30 * std::pow(T_K, 4.1881) /
             std::pow(1.0 + 6.761e-6 * T_K, 5.6881)) *
            std::exp(-54657.4 / T_K);
  c.k21LTE = 1.3e-9 * std::exp(-53300.0 / T_K);
  c.k37v0 = std::pow(10.0, -27.029 + 3.801 * c.lgT - 29487.0 / T_K);
  c.k37LTE = std::pow(10.0, -2.729 - 1.75 * c.lgT - 23474.0 / T_K);
  c.k40v0 = 2.7e-8 * std::pow(T_K, -1.27) * std::exp(-43000.0 / T_K);
  c.k47v0 = 5.09e-9 * std::pow(T_K, 0.128) * std::exp(-103258.0 / T_K);
  c.k47LTE = 1.04e-9 * std::pow(T_K, 0.218) * std::exp(-53070.7 / T_K);
  c.k75v0 = 1.35e-9 * std::pow(T_K, -1.27) * std::exp(-43000.0 / T_K);

  c.delE_22 = delE_22;
  return c;
}

// ---------------------------------------------------------------------------
// LTE / low-density blend of two rate coefficients.
//
// Blends a low-density (v0) and an LTE rate in log10 space, weighted by the
// critical-density ratio crit_ratio:
//   log10(k) = crit_ratio/(1+crit_ratio) * log10(k_LTE) + 1/(1+crit_ratio) *
//   log10(k_v0)
//
// lte_blend         — guarded form: returns 0 if either input is exactly 0
//                     (the log10 would be undefined). Used where v0 or LTE can
//                     legitimately evaluate to 0.
// lte_blend_noguard — unguarded form: the caller guarantees both inputs are
//                     strictly positive (e.g. products of exp()), so no zero
//                     check is performed. Kept separate to preserve the exact
//                     original control flow (bit-for-bit).
// ---------------------------------------------------------------------------
inline double lte_blend(double k_v0, double k_LTE, double crit_ratio) {
  if (k_v0 == 0.0 || k_LTE == 0.0) return 0.0;
  double lgk = (crit_ratio / (1.0 + crit_ratio)) * std::log10(k_LTE) +
               (1.0 / (1.0 + crit_ratio)) * std::log10(k_v0);
  return std::pow(10.0, lgk);
}

inline double lte_blend_noguard(double k_v0, double k_LTE, double crit_ratio) {
  double lgk = (crit_ratio / (1.0 + crit_ratio)) * std::log10(k_LTE) +
               (1.0 / (1.0 + crit_ratio)) * std::log10(k_v0);
  return std::pow(10.0, lgk);
}

// ---------------------------------------------------------------------------
// compute_reverse_loop — main reverse-rate loop via detailed balance.
//
// For each reaction ire in [0, n_loop), computes the equilibrium constant
// lnKeqb[num-1] and the corresponding reverse rate k_rxn[num-1 + N_react].
// Photodissociation reactions (p3 == IDX_PHOTON) include an escape-factor
// correction based on the continuum optical depth tau_cnt.
//
// The Form template parameter preserves expression-form differences:
//   Primordial:  dE/kB/T (two divisions), ternary zero-guard on kfwd
//   MetalGrain:  dE/(kB*T), updates lnKeqb for photo escape factor
//
// Plain `inline` (not force-inlined): benchmarking showed no wall-clock
// benefit from __attribute__((always_inline)) on the rate/cooling kernels
// (forced inlining was marginally slower, within noise), so inlining is left
// to the compiler.
// ---------------------------------------------------------------------------
template <RateForm Form, int N_sp, int N_react>
inline void compute_reverse_loop(std::array<double, 2 * N_react>& k_rxn,
                                 double T_K, double tau_cnt, int n_loop,
                                 const ReactionTable<N_sp, N_react>& tbl,
                                 const std::array<double, N_sp + 3>& pf,
                                 std::array<double, N_react>& lnKeqb) {
  constexpr double k_B = phys::k_B;
  constexpr double h_P = phys::h_P;
  constexpr double pi = phys::pi;

  const double lnC1_base = std::log(2.0 * pi * k_B * T_K / (h_P * h_P));

  for (int ire = 0; ire < n_loop; ++ire) {
    int num = tbl.reactions[ire].num;
    if (num < 1 || num > N_react) continue;
    int r1 = tbl.reactions[ire].reactants[0],
        r2 = tbl.reactions[ire].reactants[1],
        r3 = tbl.reactions[ire].reactants[2];
    int p1 = tbl.reactions[ire].products[0],
        p2 = tbl.reactions[ire].products[1],
        p3 = tbl.reactions[ire].products[2];
    int nr = tbl.reactions[ire].n_reactants, np = tbl.reactions[ire].n_products;
    double Cm = tbl.reactions[ire].Cmass;
    double dE = tbl.reactions[ire].delE;

    double lnC1 = 1.5 * (nr - np) * lnC1_base;
    double lnCm = std::log(Cm);
    double lnCpf = std::log(pf[r1]) + std::log(pf[r2]) + std::log(pf[r3]) -
                   std::log(pf[p1]) - std::log(pf[p2]) - std::log(pf[p3]);

    double lnKeqb_val;
    if constexpr (Form == RateForm::Primordial)
      lnKeqb_val = lnC1 + lnCm + lnCpf - dE / k_B / T_K;
    else
      lnKeqb_val = lnC1 + lnCm + lnCpf - dE / (k_B * T_K);

    lnKeqb[num - 1] = lnKeqb_val;

    int n = num - 1;
    if (p3 == ReactionTable<N_sp, N_react>::IDX_PHOTON) {
      if (tau_cnt <= 1.0e-10) {
        k_rxn[n + N_react] = 0.0;
      } else {
        double esc_fact = (tau_cnt <= 1.0e-3)
                              ? tau_cnt - 0.5 * tau_cnt * tau_cnt
                              : 1.0 - std::exp(-tau_cnt);
        if constexpr (Form == RateForm::MetalGrain) {
          lnKeqb_val += std::log(esc_fact);
          lnKeqb[n] = lnKeqb_val;
          k_rxn[n + N_react] = std::exp(std::log(k_rxn[n]) + lnKeqb_val);
        } else {
          double lnKeqb_esc = lnKeqb_val + std::log(esc_fact);
          k_rxn[n + N_react] = (k_rxn[n] > 0.0)
                                   ? std::exp(std::log(k_rxn[n]) + lnKeqb_esc)
                                   : 0.0;
        }
      }
    } else {
      if constexpr (Form == RateForm::Primordial)
        k_rxn[n + N_react] =
            (k_rxn[n] > 0.0) ? std::exp(std::log(k_rxn[n]) + lnKeqb_val) : 0.0;
      else
        k_rxn[n + N_react] = std::exp(std::log(k_rxn[n]) + lnKeqb_val);
    }

    if (k_rxn[n] == 0.0) k_rxn[n + N_react] = 0.0;
  }
}

// ---------------------------------------------------------------------------
// make_y_ext<N_sp>
//
// Builds the extended species vector consumed by the reaction loops:
//   indices [0, N_sp)       — real species abundances copied from y
//   index   N_sp            — IDX_VACANT  (empty reactant/product slot)
//   index   N_sp + 1        — IDX_PHOTON
//   index   N_sp + 2        — IDX_CR      (cosmic-ray particle)
//
// The three sentinel slots are set to 1.0 so that an "absent" reactant
// (slot index == N_sp/N_sp+1/N_sp+2) contributes a neutral factor under the
// product y_ext[r1] * y_ext[r2] * y_ext[r3] in compute_rates(): a reaction
// with fewer than three reactants points its unused slots at a sentinel and
// the corresponding factor is 1.0 (no contribution).
// ---------------------------------------------------------------------------
template <int N_sp>
inline std::array<double, N_sp + 3> make_y_ext(
    const std::array<double, N_sp>& y) {
  std::array<double, N_sp + 3> y_ext;
  y_ext.fill(0.0);
  for (int i = 0; i < N_sp; ++i) y_ext[i] = y[i];
  y_ext[N_sp] = 1.0;      // IDX_VACANT: neutral element for multiplication
  y_ext[N_sp + 1] = 1.0;  // IDX_PHOTON
  y_ext[N_sp + 2] = 1.0;  // IDX_CR
  return y_ext;
}

}  // namespace arche
