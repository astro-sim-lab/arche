// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// collapse_driver.h — shared utilities for one-zone collapse apps
//
// Extracts duplicated logic from collapse_primordial and collapse_metal_grain
// into reusable inline functions.  Both apps include this header and call
// these helpers within their model-specific RunCollapse routines.
// ---------------------------------------------------------------------------
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "collapse_defaults.h"

namespace collapse_driver {

using collapse_defaults::ExitReason;
using collapse_defaults::kGGrav;
using collapse_defaults::kKB;
using collapse_defaults::kMp;
using collapse_defaults::kPi;
using collapse_defaults::kTHighStop;

// ---------------------------------------------------------------------------
// LoadFretTable — read a 2-column ASCII table: nH[cm^-3]  f_ret
//
// Format: rows sorted by ascending nH; lines starting with '#' are comments.
// The returned vectors have equal length >= 1.
// The step-function convention: f_ret[i] applies while nH < nH[i+1].
// f_ret[last] applies for all nH >= nH[last].
// ---------------------------------------------------------------------------
inline std::pair<std::vector<double>, std::vector<double>> LoadFretTable(
    const std::string& path, const char* env_label) {
  std::ifstream fin(path);
  if (!fin.is_open()) {
    std::fprintf(stderr, "ERROR: cannot open %s '%s'\n", env_label,
                 path.c_str());
    std::exit(1);
  }
  std::vector<double> nH_tab, fret_tab;
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    double nH, fr;
    if (ss >> nH >> fr) {
      nH_tab.push_back(nH);
      fret_tab.push_back(fr);
    }
  }
  if (nH_tab.empty()) {
    std::fprintf(stderr, "ERROR: %s '%s' contains no valid rows\n", env_label,
                 path.c_str());
    std::exit(1);
  }
  return {nH_tab, fret_tab};
}

// ---------------------------------------------------------------------------
// update_fret — advance f_ret step-function when nH crosses next threshold
// ---------------------------------------------------------------------------
inline void update_fret(int& fret_idx, double& f_ret, double nH,
                        const std::vector<double>& fret_nH,
                        const std::vector<double>& fret_val) {
  while (fret_idx + 1 < static_cast<int>(fret_nH.size()) &&
         nH >= fret_nH[fret_idx + 1]) {
    ++fret_idx;
    f_ret = fret_val[fret_idx];
  }
}

// ---------------------------------------------------------------------------
// compute_tchem — chemistry timescale: min_i( y_i / |Δy_i / Δt| )
// ---------------------------------------------------------------------------
template <int N_sp>
inline double compute_tchem(double dt, const std::array<double, N_sp>& y,
                            const std::array<double, N_sp>& y_prev) {
  double t_chem = 1.0e50;
  for (int i = 0; i < N_sp; ++i) {
    const double dy = std::abs(y[i] - y_prev[i]);
    if (y[i] > 1.0e-30 && dy > 1.0e-40)
      t_chem = std::min(t_chem, dt * y[i] / dy);
  }
  return std::max(t_chem, dt);
}

// ---------------------------------------------------------------------------
// compute_diagnostics — Jeans mass and critical magnetic field
// ---------------------------------------------------------------------------
inline void compute_diagnostics(double rho, double lmbd_J, double& MJ,
                                double& B_cr) {
  MJ = (4.0 * kPi / 3.0) * rho * lmbd_J * lmbd_J * lmbd_J;
  B_cr = std::sqrt(4.0 * kPi * kGGrav * MJ * rho / lmbd_J);
}

// ---------------------------------------------------------------------------
// compute_timestep — dt update for next iteration
// ---------------------------------------------------------------------------
inline double compute_timestep(int it, int n_init_steps, double t_cool,
                               double t_eff, double t_chem, double dt_factor,
                               double dt_factor_chem, double dt_factor_init) {
  if (it <= n_init_steps) return dt_factor_init * t_eff;
  return std::min(dt_factor * std::min(t_cool, t_eff), dt_factor_chem * t_chem);
}

// ---------------------------------------------------------------------------
// update_thermodynamics — advance density, energy, temperature, pressure
//
// Computes drho and de from the current state, then updates rho, e, T_K, p,
// nH in place.  Returns false if e goes non-positive (caller should break).
// When clamp_T is true (metal model), T_K is clamped to >= 1 K.
// ---------------------------------------------------------------------------
inline bool update_thermodynamics(double dt, double t_eff, double Lambda_net,
                                  double mu, double gamma, double yHe_factor,
                                  double& rho, double& e, double& T_K,
                                  double& p, double& nH, double& t,
                                  bool clamp_T = false) {
  double drho = dt * rho / t_eff;
  double de = -Lambda_net * dt + drho * p / (rho * rho);
  rho += drho;
  e += de;
  if (e <= 0.0) return false;
  T_K = e * ((gamma - 1.0) * mu * kMp) / kKB;
  if (clamp_T) T_K = std::max(T_K, 1.0);
  p = rho * kKB * T_K / (mu * kMp);
  nH = rho / (yHe_factor * kMp);
  t += dt;
  return true;
}

// ---------------------------------------------------------------------------
// exit_code / exit_message — map ExitReason to numeric code and string
// ---------------------------------------------------------------------------
inline int to_exit_code(ExitReason reason) {
  switch (reason) {
    case ExitReason::Normal:
      return 0;
    case ExitReason::MaxIter:
      return 1;
    case ExitReason::HighTemp:
      return 2;
    case ExitReason::NegEnergy:
      return 3;
    case ExitReason::NonFinite:
      return 4;
    case ExitReason::SolverFailed:
      return 5;
    default:
      return -1;
  }
}

inline const char* exit_message(ExitReason reason) {
  switch (reason) {
    case ExitReason::Normal:
      return "Normal: nH reached nH_stop";
    case ExitReason::MaxIter:
      return "Abnormal: max_iter reached";
    case ExitReason::HighTemp:
      return "Normal: T_K reached 1e5 K ceiling";
    case ExitReason::NegEnergy:
      return "Abnormal: internal energy e <= 0";
    case ExitReason::NonFinite:
      return "Abnormal: NaN/Inf in nH or T_K";
    case ExitReason::SolverFailed:
      return "Abnormal: NR subcycle max depth exceeded";
    default:
      return "Unknown";
  }
}

// ---------------------------------------------------------------------------
// report_exit — print exit status and return exit code
// ---------------------------------------------------------------------------
inline int report_exit(ExitReason reason, const char* label) {
  const int code = to_exit_code(reason);
  const char* msg = exit_message(reason);
  std::printf("=== EXIT code=%d: %s ===\n", code, msg);
  // Normal (nH ceiling) and HighTemp (T ceiling) are both expected ends of a
  // collapse — reaching the validated T_max = 1e5 K boundary is not an error.
  if (reason != ExitReason::Normal && reason != ExitReason::HighTemp)
    std::fprintf(stderr, "WARNING [%s]: %s\n", label, msg);
  return code;
}

// ---------------------------------------------------------------------------
// check_exit — check exit conditions; returns true + sets reason if should exit
// ---------------------------------------------------------------------------
inline bool check_exit(double nH, double T_K, double e, double nH_stop,
                       ExitReason& reason) {
  if (nH > nH_stop) {
    reason = ExitReason::Normal;
    return true;
  }
  if (T_K > kTHighStop) {
    reason = ExitReason::HighTemp;
    return true;
  }
  if (!std::isfinite(nH) || !std::isfinite(T_K)) {
    reason = ExitReason::NonFinite;
    return true;
  }
  if (e <= 0.0) {
    reason = ExitReason::NegEnergy;
    return true;
  }
  return false;
}

}  // namespace collapse_driver
