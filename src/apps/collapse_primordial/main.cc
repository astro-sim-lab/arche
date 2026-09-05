// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
//
// main.cc — one-zone gravitational collapse (zero-metal network)
//
// Single case per run; parameters are read from environment variables:
//   PRIM_ZETA0          — CR ionization rate [s^-1]  (required)
//   PRIM_OUTDIR         — output directory            (optional, default:
//   results/prim/h5) PRIM_FF_RET         — free-fall retardation factor
//   (optional, default: 1.0) PRIM_FRET_TABLE     — 2-column ASCII file:
//   nH[cm^-3]  f_ret  step-function table
//                         Rows sorted by ascending nH; comments with '#'
//                         allowed. If set, PRIM_FF_RET is ignored.  fret_tag =
//                         "-step" → output: collapse_CR<tag>_fret-step.h5
//   PRIM_FF_GAMMA       — gamma-dependent collapse factor (flag; set to 1 to
//   enable)
//                         Uses t_eff = t_ff / sqrt(1−f(γ))  (Higuchi+2018
//                         Eq.5-7) Overrides PRIM_FF_RET and PRIM_FRET_TABLE.
//                         fret_tag = "-gamma"  → output:
//                         collapse_CR<tag>_fret-gamma.h5
//   PRIM_XNH0           — initial H number density [cm^-3] (optional, default:
//   0.1) PRIM_OUTPUT_STRIDE  — write every N-th step to HDF5  (optional,
//   default: 100) PRIM_MAX_ITER       — maximum integration steps (optional,
//   default: 10000000) PRIM_CHEM_TABLE     — path to primordial.h5 chemistry
//   table
//                         (optional, default: compile-time PRIM_CHEM_TABLE)
//   PRIM_JLW21          — Lyman-Werner radiation intensity J_21 [10^-21
//   erg/s/cm^2/Hz/sr]
//                         (optional, default: 0.0 = no LW field)
//                         Activates H2/HD photodissociation and H-
//                         photodetachment.
//   PRIM_ABUNDANCE_SET  — abundance preset name (optional, default: solar)
//                         currently supported: solar, default
//
// HDF5 layout for each collapse_CR<tag>[_fret<fr>].h5
// ─────────────────────────────────────────────────────────────────────────────
//   Attributes (root):
//     description   — human-readable label
//     cr_tag        — CR tag string derived from PRIM_ZETA0
//     zeta0_cgs     — CR ionization rate [s^-1]
//     f_ret         — free-fall retardation factor (initial value; 1.0 =
//     standard free-fall) f_ret_table   — path to f_ret step-function table
//     file (absent if not used) ff_collapse_mode — "gamma" when PRIM_FF_GAMMA
//     mode is active (absent otherwise) network       — "zero_metal N_sp=23
//     N_react=140" units_density — "cm^-3 (nH), g/cm^3 (rho)" units_cooling —
//     "erg g^-1 s^-1" units_time    — "s" units_length  — "cm" units_B       —
//     "G"
//
//   Datasets (all 1-D length N_rows, except y):
//     step          (N_rows,)       int32   — physical step number (multiple of
//     100) y             (N_rows, N_sp)  float64 — species abundances (number
//     fraction / nH)
//                     attr "species" = "H,H2,e-,H+,H2+,H3+,H-,He,He+,He++,HeH+,
//                                       D,HD,D+,HD+,D-,Li,LiH,Li+,Li-,LiH+,Li++,Li+++"
//     nH           (N_rows,)       — H number density [cm^-3]
//     T_K           (N_rows,)       — gas temperature [K]
//     rho           (N_rows,)       — mass density [g cm^-3]
//     Lambda_net     (N_rows,)       — net cooling Λ_line+Λ_cnt+Λ_chem − Γ_CR
//     [erg g^-1 s^-1] Lambda_line    (N_rows,)       — total line cooling
//     (H2+HD+Lya) Lambda_cnt     (N_rows,)       — continuum (dust+H ff+H2 CIA)
//     cooling Lambda_chem      (N_rows,)       — chemical (endothermic
//     reaction) cooling Gamma_cmp      (N_rows,)       — compressional heating
//     p/ρ/t_eff Lambda_gas     (N_rows,)       — gas (H ff + H2 CIA) cooling
//     subset of cnt Lambda_Lya     (N_rows,)       — Lyman-alpha cooling
//     Lambda_H2      (N_rows,)       — H2 line cooling
//     Lambda_HD      (N_rows,)       — HD line cooling
//     Gamma_CR       (N_rows,)       — CR ionization heating
//     t_ff          (N_rows,)       — true free-fall time [s]  (= t_eff /
//     f_ret) t_cool        (N_rows,)       — cooling time e/|Λ_net| [s] t_chem
//     (N_rows,)       — chemistry time scale [s]  min_i(y_i/|Δy_i/Δt|) tau_cnt
//     (N_rows,)       — continuum optical depth lambda_J       (N_rows,) —
//     Jeans length [cm] M_J           (N_rows,)       — Jeans mass [g] B_cr
//     (N_rows,)       — critical (ambipolar) magnetic field [G] y_plus
//     (N_rows,)       — total positive charge fraction y_minus       (N_rows,)
//     — total negative charge fraction charge_imbal  (N_rows,)       — |y+ −
//     y−| / (y+ + y−)

#include <hdf5.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "api/arche_api.h"  // non-template facade: PrimCell / PrimMinimalCell
#include "collapse_defaults.h"
#include "collapse_driver.h"
#include "collapse_dynamics.h"
#include "cooling/cooling.h"         // c_H2 (initial mu/gamma estimate)
#include "cooling/xray_secondary.h"  // xray::sigma_HI_Verner96 (ARCHE_XRAY)
#include "core/hdf5_utils.h"
#include "models/primordial/minimal.h"  // zero_metal_minimal::Sp / N_sp

namespace {

// ─── compile-time network variant (full vs. compact "minimal") ───────────────
// Both builds drive their model through the libarche facade (api/arche_api.h):
// the full build uses Nakauchi2019 (PrimCell), the reduced build
// (-DARCHE_PRIM_MINIMAL) uses the compact 15-species / 33-reaction
// Nakauchi2019_Minimal (PrimMinimalCell).  Each variant defines its species
// count (kNSp), cell/state/step entry points, reaction-table type (AppTable)
// and creator, plus the HDF5 `network` label and an output-filename tag.
#ifdef ARCHE_PRIM_MINIMAL
using AppCell = arche::PrimMinimalCell;
using AppTable = arche::PrimMinimalTable;
using AppSp = arche::zero_metal_minimal::Sp;  // compact species indices
constexpr int kNSp = arche::zero_metal_minimal::N_sp;  // 15
inline double AppChargePlus(const std::array<double, kNSp>& y) {
  return y[AppSp::Hp] + y[AppSp::H2p] + y[AppSp::H3p] + y[AppSp::Hep] +
         y[AppSp::Dp] + y[AppSp::Lip];
}
inline double AppChargeMinus(const std::array<double, kNSp>& y) {
  return y[AppSp::e] + y[AppSp::Hm];
}
inline arche::PrimMinimalCellPtr AppCellCreateOwned() {
  return arche::prim_minimal_cell_create_owned();
}
inline arche::ChemStatePrimMinimal& AppCellState(AppCell& cell) {
  return arche::prim_minimal_cell_state(cell);
}
inline arche::PrimMinimalTablePtr AppTableCreate(const std::string& path) {
  return arche::load_minimal_table_owned(path);
}
inline arche::ChemFullRates AppChemFullStep(AppCell& cell, double dt,
                                            const arche::ChemParams& params,
                                            const arche::ChemShielding& shield,
                                            const AppTable& tbl) {
  return arche::chem_full_step_prim_minimal(cell, dt, params, shield, tbl);
}
inline int AppTableNInvariants(const AppTable& tbl) {
  return arche::prim_minimal_table_n_invariants(tbl);
}
constexpr const char* kNetworkLabel =
    "zero_metal_minimal N_sp=15 N_react=33 "
    "(Nakauchi2019 minimal: compact 15-species / 33-reaction network)";
constexpr const char* kOutputTag = "_min";
constexpr const char* kSpeciesCsv =
    "H,H2,e-,H+,H2+,H3+,H-,He,He+,D,HD,D+,Li,LiH,Li+";
#else
using AppCell = arche::PrimCell;
using AppTable = arche::PrimTable;
using AppSp = arche::zero_metal::Sp;           // full species indices
constexpr int kNSp = arche::zero_metal::N_sp;  // 23
inline double AppChargePlus(const std::array<double, kNSp>& y) {
  return y[AppSp::Hp] + y[AppSp::H2p] + y[AppSp::H3p] + y[AppSp::Hep] +
         y[AppSp::HeHp] + y[AppSp::Dp] + y[AppSp::HDp] + y[AppSp::Lip] +
         y[AppSp::LiHp] + 2.0 * (y[AppSp::Hepp] + y[AppSp::Lipp]) +
         3.0 * y[AppSp::Lippp];
}
inline double AppChargeMinus(const std::array<double, kNSp>& y) {
  return y[AppSp::e] + y[AppSp::Hm] + y[AppSp::Dm] + y[AppSp::Lim];
}
inline arche::PrimCellPtr AppCellCreateOwned() {
  return arche::prim_cell_create_owned();
}
inline arche::ChemStateZM& AppCellState(AppCell& cell) {
  return arche::prim_cell_state(cell);
}
inline arche::PrimTablePtr AppTableCreate(const std::string& path) {
  return arche::load_prim_table_owned(path);
}
inline arche::ChemFullRates AppChemFullStep(AppCell& cell, double dt,
                                            const arche::ChemParams& params,
                                            const arche::ChemShielding& shield,
                                            const AppTable& tbl) {
  return arche::chem_full_step_prim(cell, dt, params, shield, tbl);
}
inline int AppTableNInvariants(const AppTable& tbl) {
  return arche::prim_table_n_invariants(tbl);
}
constexpr const char* kNetworkLabel = "zero_metal N_sp=23 N_react=140";
constexpr const char* kOutputTag = "";
constexpr const char* kSpeciesCsv =
    "H,H2,e-,H+,H2+,H3+,H-,He,He+,He++,HeH+,"
    "D,HD,D+,HD+,D-,Li,LiH,Li+,Li-,LiH+,Li++,Li+++";
#endif

// ─── physical constants (CGS) — from collapse_defaults.h ─────────────────────
using collapse_defaults::kGGrav;
using collapse_defaults::kKB;
using collapse_defaults::kMp;
using collapse_defaults::kPi;

// ─── CMB / integration control / exit — from collapse_defaults.h ─────────────
using collapse_defaults::ExitReason;
using collapse_defaults::kCrAttenuColDens;
using collapse_defaults::kDtFactor;
using collapse_defaults::kDtFactorChem;
using collapse_defaults::kDtFactorInit;
using collapse_defaults::kNInitSteps;
using collapse_defaults::kTCMB0;
using collapse_defaults::kXnHStop;

// ─── Integration control — from collapse_defaults.h ──────────────────────────
using collapse_defaults::kItMax;
using collapse_defaults::kOutputStride;

// ─── Initial conditions / scenario defaults — from collapse_defaults.h
// ────────
using collapse_defaults::kFFRet;
using collapse_defaults::kJLW21;
using collapse_defaults::kTK0;
using collapse_defaults::kXnH0;
using collapse_defaults::kYe0;
using collapse_defaults::kYH2;
using collapse_defaults::kYHD;
using collapse_defaults::kZeta0;
using collapse_defaults::kZred;

// ─────────────────────────────────────────────────────────────────────────────
// OutputRow — one record (every 100 steps) buffered before writing to HDF5
// ─────────────────────────────────────────────────────────────────────────────
struct OutputRow {
  int step;
  std::array<double, kNSp> y;
  double nH, T_K, rho;
  double Lambda_net, Lambda_line, Lambda_cnt, Lambda_chem, Gamma_cmp;
  double Lambda_gas, Lambda_Lya, Lambda_H2, Lambda_HD, Gamma_CR;
  double t_ff, t_cool, t_chem, tau_cnt, lmbd_J, MJ, B_cr;
  double y_plus, y_minus, charge_imbal;
};

// ─── HDF5 helpers (shared implementation in public/include/hdf5_utils.h) ─────
using h5utils::H5Create;
using h5utils::H5Write1d;
using h5utils::H5Write1dInt;
using h5utils::H5Write2d;
using h5utils::H5WriteDblAttr;
using h5utils::H5WriteIntAttr;
using h5utils::H5WriteStrAttr;

// ─────────────────────────────────────────────────────────────────────────────
// WriteHdf5File — write all datasets flat at the root of an open HDF5 file
// ─────────────────────────────────────────────────────────────────────────────
void WriteHdf5File(hid_t fid, const std::string& cr_tag,
                   const std::vector<OutputRow>& rows, double zeta0,
                   double f_ret, double T_rad, double zred, double T_K0_ic,
                   double y_e0_ic, double y_H2_ic, double y_HD_ic,
                   double jlw21 = 0.0, const std::string& fret_table_path = "",
                   bool ff_gamma = false) {
  const hsize_t N = static_cast<hsize_t>(rows.size());
  const hsize_t Nsp = static_cast<hsize_t>(kNSp);

  // ── Root attributes ───────────────────────────────────────────────────────
  {
    char desc[128];
    std::snprintf(desc, sizeof(desc),
                  "zero-metal collapse CR%s (C++ port of Nakauchi+2021)",
                  cr_tag.c_str());
    H5WriteStrAttr(fid, "description", desc);
  }
  // Output schema version (N-B). v1 = legacy x-prefixed dataset names
  // (nH/Lambda_*/Gamma_*/MJ/lmbd_J); v2 = renamed (nH/Lambda_*/Gamma_*/M_J/
  // lambda_J). Post-processing tools branch on this attribute.
  H5WriteIntAttr(fid, "schema_version", 2);
  H5WriteStrAttr(fid, "cr_tag", cr_tag);
  H5WriteDblAttr(fid, "zeta0_cgs", zeta0);
  H5WriteDblAttr(fid, "f_ret", f_ret);
  if (!fret_table_path.empty())
    H5WriteStrAttr(fid, "f_ret_table", fret_table_path);
  if (ff_gamma) H5WriteStrAttr(fid, "ff_collapse_mode", "gamma");
  H5WriteDblAttr(fid, "zred", zred);
  H5WriteDblAttr(fid, "T_rad", T_rad);
  H5WriteDblAttr(fid, "J_LW21", jlw21);
  H5WriteDblAttr(fid, "ic_T_K0", T_K0_ic);
  H5WriteDblAttr(fid, "ic_y_e0", y_e0_ic);
  H5WriteDblAttr(fid, "ic_y_H2", y_H2_ic);
  H5WriteDblAttr(fid, "ic_y_HD", y_HD_ic);
  H5WriteStrAttr(fid, "network", kNetworkLabel);
  H5WriteStrAttr(fid, "units_density", "cm^-3 (nH)  g/cm^3 (rho)");
  H5WriteStrAttr(fid, "units_cooling", "erg g^-1 s^-1");
  H5WriteStrAttr(fid, "units_time", "s");
  H5WriteStrAttr(fid, "units_length", "cm");
  H5WriteStrAttr(fid, "units_B", "G");

  // ── step ─────────────────────────────────────────────────────────────────
  {
    std::vector<int> v(N);
    for (hsize_t i = 0; i < N; ++i) v[i] = rows[i].step;
    H5Write1dInt(fid, "step", v);
  }

  // ── y (N × kNSp) — species abundances ────────────────────────────────────
  {
    std::vector<double> mat(N * Nsp);
    for (hsize_t i = 0; i < N; ++i)
      for (hsize_t j = 0; j < Nsp; ++j)
        mat[i * Nsp + j] = rows[i].y[static_cast<int>(j)];
    H5Write2d(fid, "y", mat, N, Nsp);

    hid_t ds = H5Dopen2(fid, "y", H5P_DEFAULT);
    H5WriteStrAttr(ds, "species", kSpeciesCsv);
    H5WriteStrAttr(ds, "units", "number fraction relative to nH");
    H5Dclose(ds);
  }

  // ── scalar 1-D datasets ───────────────────────────────────────────────────
  auto dump1d = [&](const std::string& name, auto fn) {
    std::vector<double> v(N);
    for (hsize_t i = 0; i < N; ++i) v[i] = fn(rows[i]);
    H5Write1d(fid, name, v);
  };

  dump1d("nH", [](const OutputRow& r) { return r.nH; });
  dump1d("T_K", [](const OutputRow& r) { return r.T_K; });
  dump1d("rho", [](const OutputRow& r) { return r.rho; });
  dump1d("Lambda_net", [](const OutputRow& r) { return r.Lambda_net; });
  dump1d("Lambda_line", [](const OutputRow& r) { return r.Lambda_line; });
  dump1d("Lambda_cnt", [](const OutputRow& r) { return r.Lambda_cnt; });
  dump1d("Lambda_chem", [](const OutputRow& r) { return r.Lambda_chem; });
  dump1d("Gamma_cmp", [](const OutputRow& r) { return r.Gamma_cmp; });
  dump1d("Lambda_gas", [](const OutputRow& r) { return r.Lambda_gas; });
  dump1d("Lambda_Lya", [](const OutputRow& r) { return r.Lambda_Lya; });
  dump1d("Lambda_H2", [](const OutputRow& r) { return r.Lambda_H2; });
  dump1d("Lambda_HD", [](const OutputRow& r) { return r.Lambda_HD; });
  dump1d("Gamma_CR", [](const OutputRow& r) { return r.Gamma_CR; });
  dump1d("t_ff", [](const OutputRow& r) { return r.t_ff; });
  dump1d("t_cool", [](const OutputRow& r) { return r.t_cool; });
  dump1d("t_chem", [](const OutputRow& r) { return r.t_chem; });
  dump1d("tau_cnt", [](const OutputRow& r) { return r.tau_cnt; });
  dump1d("lambda_J", [](const OutputRow& r) { return r.lmbd_J; });
  dump1d("M_J", [](const OutputRow& r) { return r.MJ; });
  dump1d("B_cr", [](const OutputRow& r) { return r.B_cr; });
  dump1d("y_plus", [](const OutputRow& r) { return r.y_plus; });
  dump1d("y_minus", [](const OutputRow& r) { return r.y_minus; });
  dump1d("charge_imbal", [](const OutputRow& r) { return r.charge_imbal; });
}

// ─────────────────────────────────────────────────────────────────────────────
// RunCollapse — drives one one-zone collapse to completion
// Writes:  <out_dir>/collapse_CR<tag>.h5
//          <out_dir>/collapse_CR<tag>_fret<fr>.h5         (f_ret != 1.0)
//          <out_dir>/collapse_CR<tag>_fret<fr>_JLW<jlw>.h5  (+ J_LW21 != 0)
//          <out_dir>/collapse_CR<tag>_fret-gamma.h5       (ff_gamma mode)
//          suffix order: _fret → _JLW → _z
//
// Collapse timescale modes (selected by ff_gamma and fret_tag):
//   (A/B) ff_gamma=false:  t_eff = f_ret * t_ff  (fixed or step-function f_ret)
//   (C)   ff_gamma=true:   t_eff = t_ff / sqrt(1-f(γ))  (Higuchi+2018 Eq.5-7)
// Affects: density update (drho), compressional heating (Γ_cmp),
// shielding lengths (lsh), and the integration timestep (dt).
// ─────────────────────────────────────────────────────────────────────────────
void RunCollapse(
    double y_e0, double T_K0, double y_H2_init, double y_HD_init, double zeta0,
    double nH1, const std::string& cr_tag, double f_ret_init,
    const std::string& fret_tag, double T_rad, double zred,
    const std::string& zred_tag, const arche::abundance::PrimordialSet& abund,
    const AppTable& tbl, const std::string& out_dir, double jlw21 = 0.0,
    const std::string& jlw_tag = "",
    double cr_atten_col_dens = kCrAttenuColDens,
    const std::vector<double>& fret_nH = {},
    const std::vector<double>& fret_val = {},
    const std::string& fret_table_path = "", double dt_factor = kDtFactor,
    double dt_factor_chem = kDtFactorChem,
    double dt_factor_init = kDtFactorInit, int n_init_steps = kNInitSteps,
    double nH_stop = kXnHStop, int output_stride = kOutputStride,
    int max_iter = kItMax, bool ff_gamma = false, bool bench_mode = false
#ifdef ARCHE_XRAY
    ,
    double zeta_X = 0.0, double E_X_eV = 300.0, const std::string& zx_tag = ""
#endif
) {
  // ── Initial conditions ────────────────────────────────────────────────────
  const double kT1 = T_K0;  // initial temperature [K]

  auto cell = AppCellCreateOwned();
  auto& st = AppCellState(*cell);
  auto& y = st.y;  // alias: y IS the cell's species vector

  double y_H2 = y_H2_init;
  double y_Hp = y_e0;
  double y_Dp = 0.0;
  double y_HD = y_HD_init;
  double y_H = 1.0 - y_Hp - 2.0 * y_H2 - y_HD;
  if (y_H <= 0.0) {
    std::fprintf(stderr,
                 "ERROR: Invalid IC: y_H = 1 - y_e0(%.3g) - 2*y_H2(%.3g) - "
                 "y_HD(%.3g) = %.3g <= 0\n",
                 y_Hp, y_H2, y_HD, y_H);
    std::exit(1);
  }
  double y_D = abund.yD - y_Dp - y_HD;
  double y_Lip = abund.yLi;

  y[AppSp::H] = y_H;                  // H
  y[AppSp::H2] = y_H2;                // H2
  y[AppSp::Hp] = y_Hp;                // H+
  y[AppSp::He] = abund.yHe;           // He
  y[AppSp::D] = y_D;                  // D
  y[AppSp::HD] = y_HD;                // HD
  y[AppSp::Dp] = y_Dp;                // D+
  y[AppSp::Li] = 0.0;                 // Li
  y[AppSp::Lip] = abund.yLi;          // Li+
  y[AppSp::e] = y_Hp + y_Dp + y_Lip;  // e-
  y_e0 = y[AppSp::e];

  double rho = ((1.0 + 4.0 * abund.yHe) * kMp) * nH1;
  double nH = nH1;
  double T_K = kT1;
  double p = (1.0 + abund.yHe) * nH1 * kKB * kT1;

  // Free-particle / thermal-particle sums for the initial mu, gamma estimate
  // (the kernel recomputes these each step).  He ions only for variants that
  // carry them; the term order matches the historical full-model expression.
  double mu_den =
      y[AppSp::H] + y[AppSp::H2] + y[AppSp::e] + y[AppSp::Hp] + y[AppSp::He];
  double g_den = y[AppSp::H] + y[AppSp::e] + y[AppSp::Hp] + y[AppSp::He];
#ifdef ARCHE_PRIM_MINIMAL
  mu_den += y[AppSp::Hep];  // compact minimal carries He+ (but not He++)
  g_den += y[AppSp::Hep];
#else
  mu_den += y[AppSp::Hep];
  mu_den += y[AppSp::Hepp];
  g_den += y[AppSp::Hep];
  g_den += y[AppSp::Hepp];
#endif
  double mu = (1.0 + 4.0 * abund.yHe) / mu_den;
  double gamma =
      1.0 + (1.0 + 4.0 * abund.yHe) /
                (mu * (1.5 * g_den + arche::c_H2(T_K) * y[AppSp::H2]));
  double e = kKB * T_K / ((gamma - 1.0) * mu * kMp);

  arche::ChemParams params{};
  params.T_rad = T_rad;
  // params.zeta is set each step via ChemShielding

  double t = 0.0;
  double dt = 1.0e-1;
  double t_chem = 1.0e-1;
  double t_cool = 0.0;
  double esc_cnt = 1.0;
  double tau_cnt = 0.0;
  double k_gas = 0.0;  // gas opacity [cm^2/g]; persists across steps
  double Nc_H = 0.0;
  double Nc_H2 = 0.0;
  double Nc_HD = 0.0;

  // ── f_ret step-function table state ──────────────────────────────────────
  double f_ret = f_ret_init;
  int fret_idx = 0;
  const bool has_fret_tab = !fret_nH.empty();

  // ── HDF5 row buffer ───────────────────────────────────────────────────────
  std::vector<OutputRow> h5rows;
  h5rows.reserve(600);

  // ── Bench timing CSV ─────────────────────────────────────────────────────
  using Clock = std::chrono::high_resolution_clock;
  std::FILE* bench_fp = nullptr;
  if (bench_mode) {
    std::string bench_path =
        out_dir + "/bench_CR" + cr_tag + kOutputTag + ".csv";
    bench_fp = std::fopen(bench_path.c_str(), "w");
    if (!bench_fp) {
      std::fprintf(stderr, "WARNING: cannot open bench file %s\n",
                   bench_path.c_str());
    } else {
      std::setvbuf(bench_fp, nullptr, _IOFBF, 1 << 20);
      std::fprintf(bench_fp, "step,nH_cm3,wall_us\n");
      std::printf("  bench file: %s\n", bench_path.c_str());
    }
  }

  // ── Time integration ──────────────────────────────────────────────────────
  ExitReason exit_reason = ExitReason::MaxIter;
  collapse_driver::ConservationTally cons_tally;
  // Read once, before the loop: this is a property of the loaded table,
  // not of the run, and reading it through the facade the app actually
  // uses is the point (see arche_api.h).
  const int cons_rows = AppTableNInvariants(tbl);
  for (int it = 1; it <= max_iter; ++it) {
    if (has_fret_tab)
      collapse_driver::update_fret(fret_idx, f_ret, nH, fret_nH, fret_val);

    double t_ff = std::sqrt(3.0 * kPi / (32.0 * kGGrav * rho));
    double t_eff = t_eff_collapse(t_ff, f_ret, gamma, ff_gamma);

    double vD_H = std::sqrt(2.0 * kKB * T_K / kMp);
    double vD_H2 = std::sqrt(2.0 * kKB * T_K / (kMp * 2.0));
    double vD_HD = std::sqrt(2.0 * kKB * T_K / (kMp * 3.0));

    double lmbd_J = std::sqrt(kPi * kKB * T_K / (kGGrav * mu * kMp * rho));
    double lsh_H = std::min(lmbd_J, 6.0 * vD_H * t_eff);
    double lsh_H2 = std::min(lmbd_J, 6.0 * vD_H2 * t_eff);
    double lsh_HD = std::min(lmbd_J, 6.0 * vD_HD * t_eff);

    Nc_H = y[AppSp::H] * nH * lsh_H;
    Nc_H2 = y[AppSp::H2] * nH * lsh_H2;
    Nc_HD = y[AppSp::HD] * nH * lsh_HD;

    tau_cnt = k_gas * rho * lmbd_J;
    esc_cnt = (tau_cnt > 1.0) ? 1.0 / (tau_cnt * tau_cnt) : 1.0;

    // ── Build shielding environment ───────────────────────────────────────
    arche::ChemShielding shield;
    shield.zeta = zeta0 * std::exp(-rho * lmbd_J / cr_atten_col_dens);
    shield.Nc_H = Nc_H;
    shield.Nc_H2 = Nc_H2;
    shield.Nc_HD = Nc_HD;
    shield.tau_cnt = tau_cnt;
    shield.esc_cnt = esc_cnt;
    shield.J_LW21 = jlw21;
#ifdef ARCHE_XRAY
    {
      // Beer-Lambert X-ray shielding:
      // tau_X = N_H^eff * sigma_HI(E_X)   [Inayoshi & Omukai 2011, MNRAS 416,
      // 2748] N_H^eff = (y_HI + 2*y_H2) * nH * lambda_J
      //   H2 contribution included: sigma(H2) ~ 2*sigma(HI) at 300 eV
      const double Nc_H_eff = (y[AppSp::H] + 2.0 * y[AppSp::H2]) * nH * lsh_H;
      const double tau_X = Nc_H_eff * arche::xray::sigma_HI_Verner96(E_X_eV);
      shield.zeta_X = zeta_X * std::exp(-tau_X);
    }
    shield.E_X_eV = E_X_eV;
#endif

    st.nH = nH;
    st.T_K = T_K;
    st.mu = mu;
    st.gamma = gamma;

    const auto y_prev = y;
    Clock::time_point t_bench;
    if (bench_fp) t_bench = Clock::now();
    const auto rates = AppChemFullStep(*cell, dt, params, shield, tbl);
    if (rates.solver_failed) {
      exit_reason = ExitReason::SolverFailed;
      break;
    }
    cons_tally.record(rates.conservation_projected, it, nH);
    if (bench_fp) {
      double wall_us =
          std::chrono::duration<double, std::micro>(Clock::now() - t_bench)
              .count();
      std::fprintf(bench_fp, "%d,%.6e,%.4f\n", it, nH, wall_us);
    }

    mu = st.mu;
    gamma = st.gamma;
    k_gas = rates.k_gas;

    const double Lambda_line = rates.Lambda_line;
    const double Lambda_H2 = rates.Lambda_H2;
    const double Lambda_HD = rates.Lambda_HD;
    const double Lambda_Lya = rates.Lambda_Lya;
    const double Lambda_gas = rates.Lambda_gas;
    const double Lambda_cnt = rates.Lambda_cnt;
    const double Lambda_ch = rates.Lambda_chem;
    const double Gamma_CR = rates.Gamma_CR;
    const double Lambda_net = rates.Lambda_net;

    t_chem = collapse_driver::compute_tchem<kNSp>(dt, y, y_prev);
    double Gamma_cmp = p / rho / t_eff;

    double MJ, B_cr;
    collapse_driver::compute_diagnostics(rho, lmbd_J, MJ, B_cr);

    double y_plus = AppChargePlus(y);
    double y_minus = AppChargeMinus(y);
    double y_sum = y_plus + y_minus;
    double charge_imbal =
        (y_sum > 0.0) ? std::abs(y_plus - y_minus) / y_sum : 0.0;

    // ── Output (every output_stride steps) ───────────────────────────────
    if (it % output_stride == 0) {
      OutputRow row{};
      row.step = it;
      row.y = y;
      row.nH = nH;
      row.T_K = T_K;
      row.rho = rho;
      row.Lambda_net = Lambda_net;
      row.Lambda_line = Lambda_line;
      row.Lambda_cnt = Lambda_cnt;
      row.Lambda_chem = Lambda_ch;
      row.Gamma_cmp = Gamma_cmp;
      row.Lambda_gas = Lambda_gas;
      row.Lambda_Lya = Lambda_Lya;
      row.Lambda_H2 = Lambda_H2;
      row.Lambda_HD = Lambda_HD;
      row.Gamma_CR = Gamma_CR;
      row.t_ff = t_ff;
      row.t_cool = t_cool;
      row.t_chem = t_chem;
      row.tau_cnt = tau_cnt;
      row.lmbd_J = lmbd_J;
      row.MJ = MJ;
      row.B_cr = B_cr;
      row.y_plus = y_plus;
      row.y_minus = y_minus;
      row.charge_imbal = charge_imbal;
      h5rows.push_back(row);
    }

    // ── Update thermodynamic state ────────────────────────────────────────
    if (!collapse_driver::update_thermodynamics(dt, t_eff, Lambda_net, mu,
                                                gamma, 1.0 + 4.0 * abund.yHe,
                                                rho, e, T_K, p, nH, t)) {
      exit_reason = ExitReason::NegEnergy;
      break;
    }

    t_cool = (Lambda_net != 0.0) ? e / std::abs(Lambda_net) : 1.0e50;
    dt = collapse_driver::compute_timestep(it, n_init_steps, t_cool, t_eff,
                                           t_chem, dt_factor, dt_factor_chem,
                                           dt_factor_init);

    // Stdout progress
    std::printf(
        "%7d %11.3E %11.3E %11.3E %11.3E %11.3E %11.3E"
        " %11.3E %11.3E %11.3E\n",
        it, y_e0, zeta0 / 1.0e-17, nH, T_K, y[AppSp::e], y[AppSp::H2], y_plus,
        y_minus, charge_imbal);

    if (collapse_driver::check_exit(nH, T_K, e, nH_stop, exit_reason)) break;
  }

  if (bench_fp) {
    std::fclose(bench_fp);
    bench_fp = nullptr;
  }

  // ── Exit status report ────────────────────────────────────────────────────
  std::string exit_label = "prim collapse CR" + cr_tag;
  int exit_code = collapse_driver::report_exit(exit_reason, exit_label.c_str());
  const char* exit_msg = collapse_driver::exit_message(exit_reason);
  collapse_driver::report_conservation(cons_tally, cons_rows);

  // ── Write HDF5 file ───────────────────────────────────────────────────────
  std::string h5_path = out_dir + "/collapse_CR" + cr_tag + kOutputTag;
  if (!fret_tag.empty()) h5_path += "_fret" + fret_tag;
  if (!jlw_tag.empty()) h5_path += "_JLW" + jlw_tag;
  if (!zred_tag.empty()) h5_path += "_z" + zred_tag;
#ifdef ARCHE_XRAY
  if (!zx_tag.empty()) h5_path += "_ZX" + zx_tag;
#endif
  h5_path += ".h5";
  hid_t fid = H5Create(h5_path);
  if (fid < 0) {
    std::fprintf(stderr, "ERROR: cannot create %s\n", h5_path.c_str());
    return;
  }
  WriteHdf5File(fid, cr_tag, h5rows, zeta0, f_ret_init, T_rad, zred, T_K0, y_e0,
                y_H2_init, y_HD_init, jlw21, fret_table_path, ff_gamma);
  H5WriteIntAttr(fid, "exit_code", exit_code);
  H5WriteStrAttr(fid, "exit_message", std::string(exit_msg));
  if (!collapse_driver::write_conservation_attrs(fid, cons_tally, cons_rows)) {
    std::fprintf(stderr,
                 "ERROR: could not write the conservation attributes to %s\n",
                 h5_path.c_str());
  }
#ifdef ARCHE_XRAY
  H5WriteDblAttr(fid, "zeta_X", zeta_X);
  H5WriteDblAttr(fid, "E_X_eV", E_X_eV);
#endif
  H5Fclose(fid);

  std::printf("  -> %s  (%zu rows)\n", h5_path.c_str(), h5rows.size());
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
int main() {
  // ── Chemistry table (compile-time default; overridable at runtime) ────────
  const char* env_chem_table = std::getenv("PRIM_CHEM_TABLE");
  const std::string chem_table = (env_chem_table && env_chem_table[0] != '\0')
                                     ? env_chem_table
                                     : PRIM_CHEM_TABLE;

  // ── Output directory ──────────────────────────────────────────────────────
  const char* env_out = std::getenv("PRIM_OUTDIR");
  std::string out_dir = env_out ? env_out : "results/prim/h5";
  std::filesystem::create_directories(out_dir);
  std::printf("Output directory: %s/\n", out_dir.c_str());

  // ── Build the reaction table (topology + partition functions) ─────────────
  auto tbl_owned =
      AppTableCreate(chem_table);  // PrimTablePtr/PrimMinimalTablePtr
  const AppTable& tbl = *tbl_owned;

  // ── Required: PRIM_ZETA0 [s^-1] ──────────────────────────────────────────
  const char* env_zeta0 = std::getenv("PRIM_ZETA0");
  if (!env_zeta0 || env_zeta0[0] == '\0') {
    std::fprintf(
        stderr,
        "ERROR: environment variable PRIM_ZETA0 is required\n"
        "       Set the CR ionization rate [s^-1], e.g.: PRIM_ZETA0=1e-17\n");
    return 1;
  }
  double zeta0 = kZeta0;
  try {
    zeta0 = std::stod(env_zeta0);
  } catch (const std::exception&) {
    std::fprintf(stderr, "ERROR: PRIM_ZETA0='%s' is not a valid number\n",
                 env_zeta0);
    return 1;
  }
  if (zeta0 < 0.0) {
    std::fprintf(stderr, "ERROR: PRIM_ZETA0 must be >= 0, got %s\n", env_zeta0);
    return 1;
  }

  std::string cr_tag = collapse_defaults::make_tag(env_zeta0, zeta0);

  // ── Optional: PRIM_REDSHIFT (cosmological redshift, default 0.0) ──────────
  // T_rad = kTCMB0 * (1 + z).  Affects H2/HD line cooling, continuum
  // cooling (T^4 - T_rad^4), molecular cooling CMB correction, and grain
  // temperature equilibrium.
  double zred = kZred;
  std::string zred_tag;  // empty → no suffix in filename
  const char* env_zred = std::getenv("PRIM_REDSHIFT");
  if (env_zred && env_zred[0] != '\0') {
    try {
      zred = std::stod(env_zred);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: PRIM_REDSHIFT='%s' is not a valid number\n",
                   env_zred);
      return 1;
    }
    if (zred < 0.0) {
      std::fprintf(stderr, "ERROR: PRIM_REDSHIFT must be >= 0, got %s\n",
                   env_zred);
      return 1;
    }
    if (zred != 0.0) {
      zred_tag = collapse_defaults::make_tag(env_zred, zred);
    }
  }
  const double T_rad = kTCMB0 * (1.0 + zred);

  // ── Optional: PRIM_FRET_TABLE (step-function f_ret table; overrides
  // PRIM_FF_RET) ── ASCII 2-col file: nH[cm^-3]  f_ret  (sorted by ascending
  // nH, '#' comments allowed)
  const char* env_fret_table = std::getenv("PRIM_FRET_TABLE");
  std::string fret_table_path;
  std::vector<double> fret_nH, fret_val;
  if (env_fret_table && env_fret_table[0] != '\0') {
    fret_table_path = env_fret_table;
    auto [nh, fv] =
        collapse_driver::LoadFretTable(fret_table_path, "PRIM_FRET_TABLE");
    fret_nH = std::move(nh);
    fret_val = std::move(fv);
  }

  // ── Optional: PRIM_FF_GAMMA (gamma-dependent collapse factor; highest
  // priority) ── When set (non-empty, non-"0"), uses t_eff = t_ff /
  // sqrt(1-f(γ))  [Higuchi+2018 Eq.5-7]. Overrides PRIM_FF_RET and
  // PRIM_FRET_TABLE.
  const char* env_ff_gamma = std::getenv("PRIM_FF_GAMMA");
  const bool ff_gamma =
      (env_ff_gamma && env_ff_gamma[0] != '\0' && env_ff_gamma[0] != '0');

  // ── Optional: PRIM_FF_RET (free-fall retardation factor, default 1.0) ────
  // f_ret > 1 slows the collapse; f_ret = 1 gives standard free-fall.
  // Ignored when PRIM_FRET_TABLE or PRIM_FF_GAMMA is set.
  const char* env_fret = std::getenv("PRIM_FF_RET");
  double f_ret = kFFRet;
  std::string fret_tag;  // empty → no suffix in filename
  if (ff_gamma) {
    // Gamma mode: tag "-gamma"; warn if other fret options were also supplied
    fret_tag = "-gamma";
    if (!fret_table_path.empty())
      std::fprintf(
          stderr,
          "WARNING: PRIM_FF_GAMMA is set; PRIM_FRET_TABLE='%s' is ignored\n",
          env_fret_table);
    if (env_fret && env_fret[0] != '\0')
      std::fprintf(
          stderr,
          "WARNING: PRIM_FF_GAMMA is set; PRIM_FF_RET='%s' is ignored\n",
          env_fret);
  } else if (!fret_table_path.empty()) {
    // Table mode: initial f_ret from first table entry; fixed tag "-step"
    // → output filename: collapse_CR<tag>_fret-step.h5
    f_ret = fret_val[0];
    fret_tag = "-step";
    if (env_fret && env_fret[0] != '\0')
      std::fprintf(
          stderr,
          "WARNING: PRIM_FRET_TABLE is set; PRIM_FF_RET='%s' is ignored\n",
          env_fret);
  } else if (env_fret && env_fret[0] != '\0') {
    try {
      f_ret = std::stod(env_fret);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: PRIM_FF_RET='%s' is not a valid number\n",
                   env_fret);
      return 1;
    }
    if (f_ret <= 0.0) {
      std::fprintf(stderr, "ERROR: PRIM_FF_RET must be > 0, got %s\n",
                   env_fret);
      return 1;
    }
    if (f_ret != 1.0) {
      fret_tag = collapse_defaults::make_tag(env_fret, f_ret);
    }
  }

  // ── Optional: PRIM_JLW21 (Lyman-Werner intensity J_21, default 0.0) ──────
  // J_21 units: 10^-21 erg/s/cm^2/Hz/sr.  0.0 = no LW field (default).
  double jlw21 = kJLW21;
  const char* env_jlw21 = std::getenv("PRIM_JLW21");
  if (env_jlw21 && env_jlw21[0] != '\0') {
    try {
      jlw21 = std::stod(env_jlw21);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: PRIM_JLW21='%s' is not a valid number\n",
                   env_jlw21);
      return 1;
    }
    if (jlw21 < 0.0) {
      std::fprintf(stderr, "ERROR: PRIM_JLW21 must be >= 0, got %s\n",
                   env_jlw21);
      return 1;
    }
  }
  // Tag: "." → "p"; zero → empty (no suffix)
  std::string jlw_tag;
  if (jlw21 != 0.0) {
    jlw_tag = collapse_defaults::make_tag(env_jlw21, jlw21);
  }

  // ── Optional: PRIM_ABUNDANCE_SET (default: solar) ───────────────────────
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

  // ── Optional: PRIM_XNH0 [cm^-3] (initial H number density) ──────────────
  double nH0 = kXnH0;
  const char* env_xnH0 = std::getenv("PRIM_XNH0");
  if (env_xnH0 && env_xnH0[0] != '\0') {
    try {
      nH0 = std::stod(env_xnH0);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: PRIM_XNH0='%s' is not a valid number\n",
                   env_xnH0);
      return 1;
    }
    if (nH0 <= 0.0) {
      std::fprintf(stderr, "ERROR: PRIM_XNH0 must be > 0, got %s\n", env_xnH0);
      return 1;
    }
  }

  // ── Optional: PRIM_TK0 [K] (initial gas temperature) ─────────────────────
  double T_K0 = kTK0;
  const char* env_T_K0 = std::getenv("PRIM_TK0");
  if (env_T_K0 && env_T_K0[0] != '\0') {
    try {
      T_K0 = std::stod(env_T_K0);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: PRIM_TK0='%s' is not a valid number\n",
                   env_T_K0);
      return 1;
    }
    if (T_K0 <= 0.0) {
      std::fprintf(stderr, "ERROR: PRIM_TK0 must be > 0, got %s\n", env_T_K0);
      return 1;
    }
  }

  // ── Optional: PRIM_YE0 (initial electron / H+ fraction) ──────────────────
  double ye0 = kYe0;
  const char* env_ye0 = std::getenv("PRIM_YE0");
  if (env_ye0 && env_ye0[0] != '\0') {
    try {
      ye0 = std::stod(env_ye0);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: PRIM_YE0='%s' is not a valid number\n",
                   env_ye0);
      return 1;
    }
    if (ye0 < 0.0 || ye0 >= 1.0) {
      std::fprintf(stderr, "ERROR: PRIM_YE0 must be in [0, 1), got %s\n",
                   env_ye0);
      return 1;
    }
  }

  // ── Optional: PRIM_YH2 (initial H2 fraction) ─────────────────────────────
  double y_H2_init = kYH2;
  const char* env_yH2 = std::getenv("PRIM_YH2");
  if (env_yH2 && env_yH2[0] != '\0') {
    try {
      y_H2_init = std::stod(env_yH2);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: PRIM_YH2='%s' is not a valid number\n",
                   env_yH2);
      return 1;
    }
    if (y_H2_init < 0.0 || y_H2_init >= 0.5) {
      std::fprintf(stderr, "ERROR: PRIM_YH2 must be in [0, 0.5), got %s\n",
                   env_yH2);
      return 1;
    }
  }

  // ── Optional: PRIM_YHD (initial HD fraction) ─────────────────────────────
  double y_HD_init = kYHD;
  const char* env_yHD = std::getenv("PRIM_YHD");
  if (env_yHD && env_yHD[0] != '\0') {
    try {
      y_HD_init = std::stod(env_yHD);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: PRIM_YHD='%s' is not a valid number\n",
                   env_yHD);
      return 1;
    }
    if (y_HD_init < 0.0) {
      std::fprintf(stderr, "ERROR: PRIM_YHD must be >= 0, got %s\n", env_yHD);
      return 1;
    }
  }

  // ── Optional: PRIM_OUTPUT_STRIDE (write every N-th step to HDF5) ─────────
  int output_stride = kOutputStride;  // default: 100
  const char* env_stride = std::getenv("PRIM_OUTPUT_STRIDE");
  if (env_stride && env_stride[0] != '\0') {
    try {
      output_stride = std::stoi(env_stride);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: PRIM_OUTPUT_STRIDE='%s' is not a valid integer\n",
                   env_stride);
      return 1;
    }
    if (output_stride <= 0) {
      std::fprintf(stderr, "ERROR: PRIM_OUTPUT_STRIDE must be > 0, got %s\n",
                   env_stride);
      return 1;
    }
  }

  // ── Optional: PRIM_MAX_ITER (maximum integration steps) ──────────────────
  int max_iter = kItMax;  // default: 10000000
  const char* env_max_iter = std::getenv("PRIM_MAX_ITER");
  if (env_max_iter && env_max_iter[0] != '\0') {
    try {
      max_iter = std::stoi(env_max_iter);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: PRIM_MAX_ITER='%s' is not a valid integer\n",
                   env_max_iter);
      return 1;
    }
    if (max_iter <= 0) {
      std::fprintf(stderr, "ERROR: PRIM_MAX_ITER must be > 0, got %s\n",
                   env_max_iter);
      return 1;
    }
  }

  // ── Optional: PRIM_CR_ATTEN_COL_DENS [g cm^-2] ──────────────────────────
  double cr_atten_col_dens = kCrAttenuColDens;
  const char* env_cr_col = std::getenv("PRIM_CR_ATTEN_COL_DENS");
  if (env_cr_col && env_cr_col[0] != '\0') {
    try {
      cr_atten_col_dens = std::stod(env_cr_col);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: PRIM_CR_ATTEN_COL_DENS='%s' is not a valid number\n",
                   env_cr_col);
      return 1;
    }
    if (cr_atten_col_dens <= 0.0) {
      std::fprintf(stderr,
                   "ERROR: PRIM_CR_ATTEN_COL_DENS must be > 0, got %s\n",
                   env_cr_col);
      return 1;
    }
  }

  // ── Optional: PRIM_DT_FACTOR ─────────────────────────────────────────────
  double dt_factor = kDtFactor;
  const char* env_dt_factor = std::getenv("PRIM_DT_FACTOR");
  if (env_dt_factor && env_dt_factor[0] != '\0') {
    try {
      dt_factor = std::stod(env_dt_factor);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: PRIM_DT_FACTOR='%s' is not a valid number\n",
                   env_dt_factor);
      return 1;
    }
    if (dt_factor <= 0.0) {
      std::fprintf(stderr, "ERROR: PRIM_DT_FACTOR must be > 0, got %s\n",
                   env_dt_factor);
      return 1;
    }
  }

  // ── Optional: PRIM_DT_FACTOR_CHEM ────────────────────────────────────────
  double dt_factor_chem = kDtFactorChem;
  const char* env_dt_factor_chem = std::getenv("PRIM_DT_FACTOR_CHEM");
  if (env_dt_factor_chem && env_dt_factor_chem[0] != '\0') {
    try {
      dt_factor_chem = std::stod(env_dt_factor_chem);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: PRIM_DT_FACTOR_CHEM='%s' is not a valid number\n",
                   env_dt_factor_chem);
      return 1;
    }
    if (dt_factor_chem <= 0.0) {
      std::fprintf(stderr, "ERROR: PRIM_DT_FACTOR_CHEM must be > 0, got %s\n",
                   env_dt_factor_chem);
      return 1;
    }
  }

  // ── Optional: PRIM_DT_FACTOR_INIT ────────────────────────────────────────
  double dt_factor_init = kDtFactorInit;
  const char* env_dt_factor_init = std::getenv("PRIM_DT_FACTOR_INIT");
  if (env_dt_factor_init && env_dt_factor_init[0] != '\0') {
    try {
      dt_factor_init = std::stod(env_dt_factor_init);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: PRIM_DT_FACTOR_INIT='%s' is not a valid number\n",
                   env_dt_factor_init);
      return 1;
    }
    if (dt_factor_init <= 0.0) {
      std::fprintf(stderr, "ERROR: PRIM_DT_FACTOR_INIT must be > 0, got %s\n",
                   env_dt_factor_init);
      return 1;
    }
  }

  // ── Optional: PRIM_N_INIT_STEPS ──────────────────────────────────────────
  int n_init_steps = kNInitSteps;
  const char* env_n_init = std::getenv("PRIM_N_INIT_STEPS");
  if (env_n_init && env_n_init[0] != '\0') {
    try {
      n_init_steps = std::stoi(env_n_init);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: PRIM_N_INIT_STEPS='%s' is not a valid integer\n",
                   env_n_init);
      return 1;
    }
    if (n_init_steps < 0) {
      std::fprintf(stderr, "ERROR: PRIM_N_INIT_STEPS must be >= 0, got %s\n",
                   env_n_init);
      return 1;
    }
  }

  // ── Optional: PRIM_XNH_STOP [cm^-3] ──────────────────────────────────────
  double nH_stop = kXnHStop;
  const char* env_xnh_stop = std::getenv("PRIM_XNH_STOP");
  if (env_xnh_stop && env_xnh_stop[0] != '\0') {
    try {
      nH_stop = std::stod(env_xnh_stop);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: PRIM_XNH_STOP='%s' is not a valid number\n",
                   env_xnh_stop);
      return 1;
    }
    if (nH_stop <= 0.0) {
      std::fprintf(stderr, "ERROR: PRIM_XNH_STOP must be > 0, got %s\n",
                   env_xnh_stop);
      return 1;
    }
  }

#ifdef ARCHE_XRAY
  // ── Optional: PRIM_ZETA_X [s^-1] (X-ray primary photoionization rate) ─────
  // 0.0 = no X-ray field (default).  Optically thin; no self-shielding.
  double zeta_X = 0.0;
  const char* env_zeta_X = std::getenv("PRIM_ZETA_X");
  if (env_zeta_X && env_zeta_X[0] != '\0') {
    try {
      zeta_X = std::stod(env_zeta_X);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: PRIM_ZETA_X='%s' is not a valid number\n",
                   env_zeta_X);
      return 1;
    }
    if (zeta_X < 0.0) {
      std::fprintf(stderr, "ERROR: PRIM_ZETA_X must be >= 0, got %s\n",
                   env_zeta_X);
      return 1;
    }
  }

  // ── Optional: PRIM_E_X_EV [eV] (representative X-ray photon energy) ──────
  double E_X_eV = 300.0;
  const char* env_E_X_eV = std::getenv("PRIM_E_X_EV");
  if (env_E_X_eV && env_E_X_eV[0] != '\0') {
    try {
      E_X_eV = std::stod(env_E_X_eV);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: PRIM_E_X_EV='%s' is not a valid number\n",
                   env_E_X_eV);
      return 1;
    }
    if (E_X_eV <= 13.6) {
      std::fprintf(
          stderr,
          "ERROR: PRIM_E_X_EV must be > 13.6 eV (HI ionization threshold),"
          " got %s\n",
          env_E_X_eV);
      return 1;
    }
  }

  // Tag: "." → "p"; zero → empty (no suffix)
  std::string zx_tag;
  if (zeta_X != 0.0) {
    zx_tag = collapse_defaults::make_tag(env_zeta_X, zeta_X);
  }
#endif  // ARCHE_XRAY

  std::printf(
      "=== PRIM_ZETA0=%s  tag=%s  f_ret=%g  fret_tag=%s  ff_gamma=%d"
      "  PRIM_REDSHIFT=%g  T_rad=%g K"
      "  abund=%s"
      "  J_LW21=%g  nH0=%g"
      "  T_K0=%g  y_e0=%g  y_H2=%g  y_HD=%g"
      "  stride=%d  max_iter=%d"
      "  cr_col=%g  dt_factor=%g  dt_factor_chem=%g  dt_factor_init=%g"
      "  n_init=%d  nH_stop=%g"
#ifdef ARCHE_XRAY
      "  ZETA_X=%g  E_X_EV=%g"
#endif
      " ===\n",
      env_zeta0, cr_tag.c_str(), f_ret, fret_tag.c_str(),
      static_cast<int>(ff_gamma), zred, T_rad, abund.name, jlw21, nH0, T_K0,
      ye0, y_H2_init, y_HD_init, output_stride, max_iter, cr_atten_col_dens,
      dt_factor, dt_factor_chem, dt_factor_init, n_init_steps, nH_stop
#ifdef ARCHE_XRAY
      ,
      zeta_X, E_X_eV
#endif
  );
  if (!fret_table_path.empty() && !ff_gamma)
    std::printf("  fret_table: %s  (%zu rows)\n", fret_table_path.c_str(),
                fret_nH.size());
  const char* env_bench = std::getenv("PRIM_BENCH");
  const bool bench_mode =
      (env_bench && env_bench[0] != '\0' && env_bench[0] != '0');

  RunCollapse(ye0, T_K0, y_H2_init, y_HD_init, zeta0, nH0, cr_tag, f_ret,
              fret_tag, T_rad, zred, zred_tag, abund, tbl, out_dir, jlw21,
              jlw_tag, cr_atten_col_dens, fret_nH, fret_val, fret_table_path,
              dt_factor, dt_factor_chem, dt_factor_init, n_init_steps, nH_stop,
              output_stride, max_iter, ff_gamma, bench_mode
#ifdef ARCHE_XRAY
              ,
              zeta_X, E_X_eV, zx_tag
#endif
  );

  return 0;
}
