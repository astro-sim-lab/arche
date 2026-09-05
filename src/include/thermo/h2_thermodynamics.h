// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// h2_thermodynamics.h — H2 internal degrees of freedom for the EOS.

#include <algorithm>
#include <cmath>

namespace arche::thermo {

// c_H2 — internal degrees-of-freedom factor for H2: returns 1/(gamma-1).
//
// The rotational term is evaluated from the para/ortho (1:3) partition
// functions through 1e6 K.  Above 1000 K, include enough rotational states
// that the final state's Boltzmann exponent is at least 60.  At still higher
// temperatures use the rigid-rotor expansion c_rot = 1 + (theta_rot/T)^2/45;
// it avoids an impractical O(sqrt(T)) state sum in the public EOS inversion.
// Retaining K_max = 20 through 1000 K preserves the historical
// low-temperature calculation bit-for-bit while removing its discontinuous
// classical switch.
inline double c_H2(double T_K) {
  constexpr double theta_rot_K = 85.4;
  constexpr int k_max_legacy = 20;
  constexpr double asymptotic_start_K = 1.0e6;
  double c_rot;
  if (T_K > asymptotic_start_K) {
    const double theta_over_T = theta_rot_K / T_K;
    c_rot = 1.0 + theta_over_T * theta_over_T / 45.0;
  } else {
    int K_max = k_max_legacy;
    if (T_K > 1.0e3) {
      const double root =
          0.5 * (std::sqrt(1.0 + 240.0 * T_K / theta_rot_K) - 1.0);
      K_max = std::max(k_max_legacy, static_cast<int>(std::ceil(root)));
    }

    constexpr double eps = 1.0e-3;
    const double T_b = (1.0 - eps) * T_K;
    const double T_f = (1.0 + eps) * T_K;
    double Zp_b = 0, Zp_f = 0, Xp_b = 0, Xp_f = 0;
    double Zo_b = 0, Zo_f = 0, Xo_b = 0, Xo_f = 0;
    for (int K = 0; K <= K_max; ++K) {
      double E = theta_rot_K * static_cast<double>(K * (K + 1));
      double wb = static_cast<double>(2 * K + 1) * std::exp(-E / T_b);
      double wf = static_cast<double>(2 * K + 1) * std::exp(-E / T_f);
      if (K % 2 == 0) {
        Zp_b += wb;
        Zp_f += wf;
        Xp_b += E * wb;
        Xp_f += E * wf;
      } else {
        Zo_b += wb;
        Zo_f += wf;
        Xo_b += E * wb;
        Xo_f += E * wf;
      }
    }
    double Erot_f = 0.25 * (Xp_f / Zp_f) + 0.75 * (Xo_f / Zo_f);
    double Erot_b = 0.25 * (Xp_b / Zp_b) + 0.75 * (Xo_b / Zo_b);
    c_rot = (Erot_f - Erot_b) / (2.0 * eps * T_K);
  }
  double x = 6.1e3 / T_K;
  double ex = std::exp(x);
  double c_vib = (x > 1.0e2) ? 0.0 : (x * x * ex) / ((ex - 1.0) * (ex - 1.0));
  return 1.5 + c_rot + c_vib;
}

}  // namespace arche::thermo
