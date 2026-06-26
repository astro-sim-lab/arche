// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
//
// main.cc — standalone primordial chemistry demo using chem_full_step()
//
// Demonstrates the chem_full_step() API for the zero-metal (primordial)
// network in isolation, without gravitational-collapse hydrodynamics.
//
// The gas cell is held at fixed density and temperature (isochoric,
// isothermal) so the focus is on the chemistry / cooling API.
//
// In a 3-D fluid code the caller would:
//   - compute zeta from an attenuation model or ray-tracing
//   - supply Nc_H/H2/HD from column integration (TreeCol, Sobolev, etc.)
//   - supply tau_cnt/esc_cnt from the opacity field
//
// Here we approximate column densities as n × L_Jeans (same as the
// 1-D collapse app) purely for illustration.
//
// Environment variables:
//   PRIM_ZETA0       — CR ionization rate [s^-1]                    (required)
//   PRIM_CHEM_TABLE  — path to primordial.h5 chemistry table (optional,
//                       default: compile-time PRIM_CHEM_TABLE)
//   PRIM_XNH      — H number density [cm^-3]         (optional, default: 1e4)
//   PRIM_T_K      — gas temperature [K]               (optional, default:
//   100.0) PRIM_YE0      — initial electron / H+ fraction    (optional,
//   default: 1e-4) PRIM_YH2      — initial H2 fraction               (optional,
//   default: 6e-7) PRIM_YHD      — initial HD fraction               (optional,
//   default: 4e-10) PRIM_NSTEPS   — number of integration steps (optional,
//   default: 200) PRIM_DT       — time step [s]                     (optional,
//   default: 1e10) PRIM_CR_ATTEN_COL_DENS — CR attenuation column density scale
//   [g cm^-2]
//                              (optional, default: 96.0, must be > 0)
//   PRIM_REDSHIFT — cosmological redshift z            (optional, default: 0.0)
//   PRIM_JLW21    — LW intensity J_21                  (optional, default: 0.0)
//   PRIM_ABUNDANCE_SET — abundance preset name         (optional, default:
//   solar)
//                        currently supported: solar, default

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include "api/arche_api.h"  // non-template facade: PrimCell, chem_full_step_prim
#include "solve/chemistry.h"  // kernel helper still used directly: c_H2

int main() {
  constexpr double kCrAttenuColDens = 96.0;

  // ── Chemistry table (HDF5) ────────────────────────────────────────────────
  const char* env_table = std::getenv("PRIM_CHEM_TABLE");
  const std::string chem_table =
      (env_table && env_table[0] != '\0') ? env_table : PRIM_CHEM_TABLE;

  // ── Required: PRIM_ZETA0 ─────────────────────────────────────────────────
  const char* env_zeta0 = std::getenv("PRIM_ZETA0");
  if (!env_zeta0 || env_zeta0[0] == '\0') {
    std::fprintf(stderr,
                 "ERROR: PRIM_ZETA0 is required. Example: PRIM_ZETA0=1e-17\n");
    return 1;
  }
  double zeta0 = 0.0;
  try {
    zeta0 = std::stod(env_zeta0);
  } catch (...) {
    std::fprintf(stderr, "ERROR: PRIM_ZETA0='%s' is not a valid number\n",
                 env_zeta0);
    return 1;
  }

  // ── Optional parameters ───────────────────────────────────────────────────
  double nH = 1.0e4;   // H number density [cm^-3]
  double T_K = 100.0;  // gas temperature [K]
  int nstep = 200;
  double dt = 1.0e10;  // time step [s]

  auto getenv_dbl = [](const char* name, double def) -> double {
    const char* s = std::getenv(name);
    if (!s || s[0] == '\0') return def;
    try {
      return std::stod(s);
    } catch (...) {
      return def;
    }
  };
  auto getenv_int = [](const char* name, int def) -> int {
    const char* s = std::getenv(name);
    if (!s || s[0] == '\0') return def;
    try {
      return std::stoi(s);
    } catch (...) {
      return def;
    }
  };

  nH = getenv_dbl("PRIM_XNH", nH);
  T_K = getenv_dbl("PRIM_T_K", T_K);
  nstep = getenv_int("PRIM_NSTEPS", nstep);
  dt = getenv_dbl("PRIM_DT", dt);
  double zred = getenv_dbl("PRIM_REDSHIFT", 0.0);
  double jlw21 = getenv_dbl("PRIM_JLW21", 0.0);

  double cr_atten_col_dens = kCrAttenuColDens;
  const char* env_cr_col = std::getenv("PRIM_CR_ATTEN_COL_DENS");
  if (env_cr_col && env_cr_col[0] != '\0') {
    try {
      cr_atten_col_dens = std::stod(env_cr_col);
    } catch (...) {
      std::fprintf(stderr,
                   "ERROR: PRIM_CR_ATTEN_COL_DENS='%s' is not a valid number\n",
                   env_cr_col);
      return 1;
    }
    if (!(cr_atten_col_dens > 0.0)) {
      std::fprintf(stderr,
                   "ERROR: PRIM_CR_ATTEN_COL_DENS must be > 0, got %s\n",
                   env_cr_col);
      return 1;
    }
  }
  if (zred < 0.0) {
    std::fprintf(stderr, "ERROR: PRIM_REDSHIFT must be >= 0, got %g\n", zred);
    return 1;
  }
  if (jlw21 < 0.0) {
    std::fprintf(stderr, "ERROR: PRIM_JLW21 must be >= 0, got %g\n", jlw21);
    return 1;
  }

  const char* env_abund = std::getenv("PRIM_ABUNDANCE_SET");
  const std::string abundance_set =
      (env_abund && env_abund[0] != '\0') ? std::string(env_abund) : "solar";
  arche::abundance::PrimordialSet abund{};
  try {
    abund = arche::abundance::get_primordial_set(abundance_set);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }

  // ── Load reaction tables ──────────────────────────────────────────────────
  arche::PrimTablePtr tbl = arche::load_prim_table_owned(chem_table);

  // ── Physical constants ────────────────────────────────────────────────────
  constexpr double kKB = arche::phys::k_B;
  constexpr double kMp = arche::phys::m_p;
  constexpr double kPi = arche::phys::pi;
  constexpr double kGGrav = arche::phys::G;

  // ── Primordial initial conditions ─────────────────────────────────────────
  arche::PrimCellPtr cell = arche::prim_cell_create_owned();
  arche::ChemStateZM& st = arche::prim_cell_state(*cell);
  auto& y = st.y;

  const double y_H2 = getenv_dbl("PRIM_YH2", 6.0e-7);
  const double y_Hp = getenv_dbl("PRIM_YE0", 1.0e-4);
  const double y_Dp = 0.0;
  const double y_HD = getenv_dbl("PRIM_YHD", 4.0e-10);
  const double y_H = 1.0 - y_Hp - 2.0 * y_H2 - y_HD;
  if (y_H <= 0.0) {
    std::fprintf(stderr,
                 "ERROR: Invalid IC: y_H = 1 - y_e0(%.3g) - 2*y_H2(%.3g) - "
                 "y_HD(%.3g) = %.3g <= 0\n",
                 y_Hp, y_H2, y_HD, y_H);
    return 1;
  }
  y[0] = y_H;                      // H
  y[1] = y_H2;                     // H2
  y[3] = y_Hp;                     // H+
  y[7] = abund.yHe;                // He
  y[11] = abund.yD - y_Dp - y_HD;  // D
  y[12] = y_HD;                    // HD
  y[13] = y_Dp;                    // D+
  y[18] = abund.yLi;               // Li+
  y[2] = y_Hp + y_Dp + abund.yLi;  // e-

  const double rho = (1.0 + 4.0 * abund.yHe) * kMp * nH;
  double mu = (1.0 + 4.0 * abund.yHe) /
              (y[0] + y[1] + y[2] + y[3] + y[7] + y[8] + y[9]);
  double gamma =
      1.0 + (1.0 + 4.0 * abund.yHe) /
                (mu * (1.5 * (y[0] + y[2] + y[3] + y[7] + y[8] + y[9]) +
                       arche::c_H2(T_K) * y[1]));

  st.nH = nH;
  st.T_K = T_K;
  st.mu = mu;
  st.gamma = gamma;

  // ── Chemistry parameters ──────────────────────────────────────────────────
  arche::ChemParams params{};
  params.T_rad = 2.725 * (1.0 + zred);  // CMB temperature [K]

  // ── Persistent opacity across steps ──────────────────────────────────────
  double k_gas = 0.0;
  double tau_cnt = 0.0;
  double esc_cnt = 1.0;

  // ── Print header ──────────────────────────────────────────────────────────
  std::printf(
      "# primordial chemistry demo: nH=%g cm^-3  T_K=%g K  zeta0=%g s^-1"
      "  y_e0=%g  y_H2=%g  y_HD=%g  cr_col=%g  zred=%g  T_rad=%g  J_LW21=%g\n",
      nH, T_K, zeta0, y_Hp, y_H2, y_HD, cr_atten_col_dens, zred, params.T_rad,
      jlw21);
  std::printf("# abundance_set=%s\n", abund.name);
  std::printf("# %5s  %10s  %11s  %11s  %11s  %11s  %11s  %11s  %11s\n", "step",
              "t [s]", "y[H]", "y[H2]", "y[e-]", "L_line", "L_cnt", "L_ch",
              "Gam_CR");
  std::printf("# %5s  %10s  %11s  %11s  %11s  %11s  %11s  %11s  %11s\n",
              "-----", "----------", "-----------", "-----------",
              "-----------", "-----------", "-----------", "-----------",
              "-----------");

  double t = 0.0;

  // ── Integration loop ──────────────────────────────────────────────────────
  for (int it = 1; it <= nstep; ++it) {
    // Jeans length (isothermal estimate)
    const double lmbd_J =
        std::sqrt(kPi * kKB * T_K / (kGGrav * mu * kMp * rho));

    // Column densities (Jeans-length approximation)
    // NOTE: in a 3-D fluid code these come from column integration.
    const double Nc_H = y[0] * nH * lmbd_J;
    const double Nc_H2 = y[1] * nH * lmbd_J;
    const double Nc_HD = y[12] * nH * lmbd_J;

    tau_cnt = k_gas * rho * lmbd_J;
    esc_cnt = (tau_cnt > 1.0) ? 1.0 / (tau_cnt * tau_cnt) : 1.0;

    // Build shielding environment (caller-supplied in a real fluid code)
    arche::ChemShielding shield;
    shield.zeta = zeta0 * std::exp(-rho * lmbd_J / cr_atten_col_dens);
    shield.Nc_H = Nc_H;
    shield.Nc_H2 = Nc_H2;
    shield.Nc_HD = Nc_HD;
    shield.J_LW21 = jlw21;
    shield.tau_cnt = tau_cnt;
    shield.esc_cnt = esc_cnt;

    st.nH = nH;  // fixed (isothermal/isochoric demo)
    st.T_K = T_K;
    st.mu = mu;
    st.gamma = gamma;

    const auto rates =
        arche::chem_full_step_prim(*cell, dt, params, shield, *tbl);

    mu = st.mu;
    gamma = st.gamma;
    k_gas = rates.k_gas;

    // Print
    std::printf(
        "  %5d  %10.3E  %11.4E  %11.4E  %11.4E"
        "  %11.4E  %11.4E  %11.4E  %11.4E\n",
        it, t, y[0], y[1], y[2], rates.Lambda_line, rates.Lambda_cnt,
        rates.Lambda_chem, rates.Gamma_CR);

    t += dt;
  }

  return 0;
}
