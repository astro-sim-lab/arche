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
// Through 1000 K the sum keeps K_max = 20; the state count is continuous at
// 1000 K, so c_H2 has no classical switch there.
//
// The Boltzmann factors exp(-theta K(K+1)/T) are the powers r^{K(K+1)} of
// r = exp(-theta/T), built by g_K = g_{K-1} s_K with s_K = r^{2K} = s_{K-1} r^2;
// every 16th level takes g_K and s_K from exp directly, so the rounding of the
// products does not accumulate over the long sums above 1000 K.  The vibrational
// term uses e^x - 1 = expm1(x), which keeps its precision as x = 6100 K / T -> 0.
namespace detail {

inline double c_H2_eval(double T_K) {
  constexpr double theta_rot_K = 85.4;
  constexpr int k_max_legacy = 20;
  constexpr int k_anchor = 16;
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
    const double rb = std::exp(-theta_rot_K / T_b);
    const double rf = std::exp(-theta_rot_K / T_f);
    const double rb2 = rb * rb, rf2 = rf * rf;
    double gb = 1.0, gf = 1.0;  // r^{K(K+1)}
    double sb = 1.0, sf = 1.0;  // r^{2K}
    double Zp_b = 0, Zp_f = 0, Xp_b = 0, Xp_f = 0;
    double Zo_b = 0, Zo_f = 0, Xo_b = 0, Xo_f = 0;
    for (int K = 0; K <= K_max; ++K) {
      const double E = theta_rot_K * static_cast<double>(K * (K + 1));
      if (K > 0 && K % k_anchor == 0) {
        const double E2K = theta_rot_K * static_cast<double>(2 * K);
        gb = std::exp(-E / T_b);
        gf = std::exp(-E / T_f);
        sb = std::exp(-E2K / T_b);
        sf = std::exp(-E2K / T_f);
      } else if (K > 0) {
        sb *= rb2;
        sf *= rf2;
        gb *= sb;
        gf *= sf;
      }
      const double wb = static_cast<double>(2 * K + 1) * gb;
      const double wf = static_cast<double>(2 * K + 1) * gf;
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
  const double x = 6.1e3 / T_K;
  double c_vib = 0.0;
  if (x <= 1.0e2) {
    const double em1 = std::expm1(x);
    c_vib = (x * x * (em1 + 1.0)) / (em1 * em1);
  }
  return 1.5 + c_rot + c_vib;
}

}  // namespace detail

// c_H2(T) = detail::c_H2_eval(T).  Evaluated at one temperature many times in a row
// (every abundance column of an energy-row Jacobian shares the temperature), it keeps
// its last value per thread.  Defined in thermo/eos.cc, so every caller shares one
// compiled evaluation and the kept value is exactly what that evaluation returns.
double c_H2(double T_K);

}  // namespace arche::thermo
