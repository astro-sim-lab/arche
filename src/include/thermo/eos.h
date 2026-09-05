// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// eos.h — project-internal, shared forward equation-of-state helpers.

#include <array>
#include <cmath>
#include <limits>

#include "core/state.h"
#include "thermo/h2_thermodynamics.h"

namespace arche::thermo {

// EOS inputs aggregated from a model-specific species abundance vector.
// mass_factor_per_h_nucleus is the mass per H nucleus in proton-mass units.
struct EosComposition {
  double mass_factor_per_h_nucleus;
  double particle_sum;
  double non_h2_particle_sum;
  double y_H2;
};

// Preserve the historical species-addition order used by eos_particle_sum().
// abundance_ref::yHe is the reference total He/H abundance, not y[Sp::He].
template <class Model>
inline EosComposition make_eos_composition(
    const std::array<double, Model::N_sp>& y) noexcept {
  using Sp = typename Model::Sp;
  double non_h2_particle_sum = y[Sp::H];
  non_h2_particle_sum += y[Sp::e];
  non_h2_particle_sum += y[Sp::Hp];
  non_h2_particle_sum += y[Sp::He];
  if constexpr (Model::cooling::has_he_ion) {
    non_h2_particle_sum += y[Sp::Hep];
  }
  if constexpr (Model::cooling::has_he_pp) {
    non_h2_particle_sum += y[Sp::Hepp];
  }
  const double y_H2 = y[Sp::H2];
  double particle_sum = y[Sp::H];
  particle_sum += y_H2;
  particle_sum += y[Sp::e];
  particle_sum += y[Sp::Hp];
  particle_sum += y[Sp::He];
  if constexpr (Model::cooling::has_he_ion) {
    particle_sum += y[Sp::Hep];
  }
  if constexpr (Model::cooling::has_he_pp) {
    particle_sum += y[Sp::Hepp];
  }
  return {1.0 + 4.0 * abundance_ref::yHe, particle_sum, non_h2_particle_sum,
          y_H2};
}

inline double mean_molecular_weight(const EosComposition& c) noexcept {
  return c.mass_factor_per_h_nucleus / c.particle_sum;
}

inline double adiabatic_index(const EosComposition& c, double T_K) noexcept {
  if (!std::isfinite(T_K) || T_K <= 0.0)
    return std::numeric_limits<double>::quiet_NaN();
  const double mu = mean_molecular_weight(c);
  return 1.0 + c.mass_factor_per_h_nucleus /
                   (mu * (1.5 * c.non_h2_particle_sum + c_H2(T_K) * c.y_H2));
}

// Specific internal energy [erg g^-1] for the supplied composition and T_K [K].
inline double specific_internal_energy(const EosComposition& c,
                                       double T_K) noexcept {
  if (!std::isfinite(T_K) || T_K <= 0.0)
    return std::numeric_limits<double>::quiet_NaN();
  const double denom = c.mass_factor_per_h_nucleus * phys::m_p;
  return phys::k_B * T_K * (1.5 * c.non_h2_particle_sum + c_H2(T_K) * c.y_H2) /
         denom;
}

// Defined in thermo/eos.cc; it remains out-of-line to preserve the
// translation-unit isolation required for reproducible builds.
double temperature_from_specific_internal_energy(const EosComposition& c,
                                                 double e_cgs) noexcept;

// Defined and explicitly instantiated in thermo/eos.cc.  Keeping the
// composition extraction out of the facade translation unit preserves its
// code-generation isolation from the large chemistry kernel.
template <class Model>
double temperature_from_specific_internal_energy(
    const std::array<double, Model::N_sp>& y, double e_cgs) noexcept;

}  // namespace arche::thermo
