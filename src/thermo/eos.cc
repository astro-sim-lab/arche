// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE

#include <cmath>

#include "model_traits.h"
#include "thermo/eos.h"

// Keep this implementation in a translation unit separate from arche_api.cc
// for numerical reproducibility. Adding it to that already-large unit changes
// the compiler's optimization budget and can alter ULP-level code generation
// for inline-template COMDAT definitions emitted there (including metal line
// cooling), which then changes existing applications linked to libarche.
// LTO/IPO removes this isolation and can reintroduce the same coupling.

namespace arche::thermo {

double temperature_from_specific_internal_energy(const EosComposition& c,
                                                 double e_cgs) noexcept {
  if (std::isnan(e_cgs)) return e_cgs;

  double log_lo = std::log(1.0);
  double log_hi = std::log(1.0e12);
  const double e_lo = specific_internal_energy(c, std::exp(log_lo));
  const double e_hi = specific_internal_energy(c, std::exp(log_hi));
  if (e_cgs <= e_lo) return std::exp(log_lo);
  if (e_cgs >= e_hi) return std::exp(log_hi);

  for (int it = 0; it < 80; ++it) {
    const double log_mid = 0.5 * (log_lo + log_hi);
    // The bracket is at full double precision once its midpoint rounds to an
    // endpoint. Further c_H2() evaluations cannot refine either bound.
    if (log_mid == log_lo || log_mid == log_hi) break;
    const double e_mid = specific_internal_energy(c, std::exp(log_mid));
    if (e_mid < e_cgs) {
      log_lo = log_mid;
    } else {
      log_hi = log_mid;
    }
  }
  return std::exp(0.5 * (log_lo + log_hi));
}

template <class Model>
double temperature_from_specific_internal_energy(
    const std::array<double, Model::N_sp>& y, double e_cgs) noexcept {
  return temperature_from_specific_internal_energy(
      make_eos_composition<Model>(y), e_cgs);
}

template double temperature_from_specific_internal_energy<Nakauchi2019>(
    const std::array<double, Nakauchi2019::N_sp>&, double) noexcept;
template double temperature_from_specific_internal_energy<Nakauchi2019_Minimal>(
    const std::array<double, Nakauchi2019_Minimal::N_sp>&, double) noexcept;
template double temperature_from_specific_internal_energy<Nakauchi2021>(
    const std::array<double, Nakauchi2021::N_sp>&, double) noexcept;
template double temperature_from_specific_internal_energy<Nakauchi2021_Minimal>(
    const std::array<double, Nakauchi2021_Minimal::N_sp>&, double) noexcept;

}  // namespace arche::thermo
