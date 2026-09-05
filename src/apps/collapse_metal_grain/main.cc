// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
//
// main.cc — one-zone gravitational collapse (metal_grain network)
//
// Single case per run; parameters are read from environment variables:
//   METAL_ZETA0          — CR ionization rate [s^-1]       (required)
//   METAL_Z_METAL        — metallicity [Z_sun]              (required)
//   METAL_OUTDIR         — output directory                 (optional, default:
//   results/metal/h5) METAL_FF_RET         — free-fall retardation factor
//   (optional, default: 1.0) METAL_FRET_TABLE     — 2-column ASCII file:
//   nH[cm^-3]  f_ret  step-function table
//                          Rows sorted by ascending nH; comments with '#'
//                          allowed. If set, METAL_FF_RET is ignored.  fret_tag
//                          = "-step" → output:
//                          collapse_CR<cr>_Z<z>_fret-step.h5
//   METAL_FF_GAMMA       — gamma-dependent collapse factor (flag; set to 1 to
//   enable)
//                          Uses t_eff = t_ff / sqrt(1−f(γ))  (Higuchi+2018
//                          Eq.5-7) Overrides METAL_FF_RET and METAL_FRET_TABLE.
//                          fret_tag = "-gamma"  → output:
//                          collapse_CR<cr>_Z<z>_fret-gamma.h5
//   METAL_XNH0           — initial H number density [cm^-3] (optional,
//   default: 1.0) METAL_OUTPUT_STRIDE  — write every N-th step to HDF5
//   (optional, default: 10) METAL_MAX_ITER       — maximum integration steps
//   (optional, default: 1000000) METAL_CHEM_TABLE     — path to metal_grain.h5
//   chemistry table
//                          (optional, default: compile-time METAL_CHEM_TABLE)
//   METAL_JLW21          — Lyman-Werner radiation intensity J_21 [10^-21
//   erg/s/cm^2/Hz/sr]
//                          (optional, default: 0.0 = no LW field)
//                          Activates H2/HD photodissociation and H-
//                          photodetachment.
//   METAL_SRA_RATE       — short-lived radionuclide ionization scaling
//                          (optional, default: 0.0; zeta_short = rate * 7.6e-19
//                          [s^-1])
//   METAL_LRA_RATE       — long-lived radionuclide ionization scaling
//                          (optional, default: 0.0; zeta_long  = rate * 1.4e-22
//                          * Z_metal [s^-1])
//   METAL_ABUNDANCE_SET  — abundance preset name (optional, default: solar)
//                          currently supported: solar, default
//
// HDF5 layout (root-level datasets):
//   y          (N_rows, N_sp=89) float64 — species abundances / nH
//   nH        (N_rows,)   — H number density [cm^-3]
//   T_K        (N_rows,)   — gas temperature [K]
//   T_gr_K     (N_rows,)   — grain temperature [K]
//   rho        (N_rows,)   — mass density [g/cm^3]
//   Lambda_net  (N_rows,)   — net cooling [erg g^-1 s^-1]
//   Lambda_line (N_rows,)   — total line cooling
//   Lambda_cnt  (N_rows,)   — continuum cooling
//   Lambda_gr   (N_rows,)   — grain continuum cooling
//   Lambda_gas  (N_rows,)   — gas continuum cooling (ff + CIA)
//   Lambda_chem   (N_rows,)   — chemical cooling
//   Gamma_cmp   (N_rows,)   — compressional heating
//   Lambda_H2   (N_rows,)   — H2 line cooling
//   Lambda_HD   (N_rows,)   — HD line cooling
//   Lambda_Lya  (N_rows,)   — Lyman-alpha cooling
//   Lambda_CO   (N_rows,)   — CO line cooling
//   Lambda_OH   (N_rows,)   — OH line cooling
//   Lambda_H2O  (N_rows,)   — H2O line cooling
//   Lambda_CII  (N_rows,)   — CII line cooling
//   Lambda_CI   (N_rows,)   — CI line cooling
//   Lambda_OI   (N_rows,)   — OI line cooling
//   Gamma_CR    (N_rows,)   — CR ionization heating
//   t_ff       (N_rows,)   — true free-fall time [s]  (= t_eff / f_ret)
//   t_cool     (N_rows,)   — cooling time [s]
//   t_chem     (N_rows,)   — chemistry timescale [s]
//   tau_cnt    (N_rows,)   — continuum optical depth
//   lambda_J    (N_rows,)   — Jeans length [cm]
//   M_J        (N_rows,)   — Jeans mass [g]
//   B_cr       (N_rows,)   — critical magnetic field [G]
//   y_plus     (N_rows,)   — total positive charge fraction
//   y_minus    (N_rows,)   — total negative charge fraction
//   charge_imbal (N_rows,) — |y+−y−|/(y++y−)
//   step       (N_rows,)   int32 — step number (every 10 steps)
//
//   Attributes (root):
//     f_ret        — free-fall retardation factor (initial value; 1.0 =
//     standard free-fall) f_ret_table  — path to f_ret step-function table file
//     (absent if not used) ff_collapse_mode — "gamma" when METAL_FF_GAMMA mode
//     is active (absent otherwise)

#include <hdf5.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "api/arche_api.h"  // non-template facade: MetalCell, chem_full_step_metal
#include "collapse_defaults.h"
#include "collapse_driver.h"
#include "collapse_dynamics.h"
#include "core/hdf5_utils.h"
#include "solve/chemistry.h"  // kernel helpers still used: c_H2, detail::, xray::

namespace {

// ─── compile-time network variant (full vs. compact "minimal") ───────────────
// Both builds drive their model through the libarche facade (api/arche_api.h):
// the full build uses Nakauchi2021 (MetalCell), the reduced build
// (-DARCHE_METAL_MINIMAL) uses the compact 40-species / 103+10-reaction
// Nakauchi2021_Minimal (MetalMinimalCell).  Each variant defines its species
// count (kNSp), species-index enum (AppSp), cell/state/step entry points,
// reaction-table type (AppTable) and creator, the charge/EOS reductions that
// differ by species set, plus the HDF5 `network` label, species list and an
// output-filename tag.
#ifdef ARCHE_METAL_MINIMAL
using AppCell = arche::MetalMinimalCell;
using AppTable = arche::MetalMinimalTable;
using AppCellPtr = arche::MetalMinimalCellPtr;
using AppTablePtr = arche::MetalMinimalTablePtr;
using AppSp = arche::metal_grain_minimal::Sp;
constexpr int kNSp = arche::metal_grain_minimal::N_sp;  // 40
inline AppCellPtr AppCellCreateOwned() {
  return arche::metal_minimal_cell_create_owned();
}
inline arche::ChemStateMetalMinimal& AppCellState(AppCell& cell) {
  return arche::metal_minimal_cell_state(cell);
}
inline AppTablePtr AppTableCreate(const std::string& path) {
  return arche::load_metal_minimal_table_owned(path);
}
inline arche::ChemFullRates AppChemFullStep(AppCell& cell, double dt,
                                            const arche::ChemParams& params,
                                            const arche::ChemShielding& shield,
                                            const AppTable& tbl) {
  return arche::chem_full_step_metal_minimal(cell, dt, params, shield, tbl);
}
inline int AppTableNInvariants(const AppTable& tbl) {
  return arche::metal_minimal_table_n_invariants(tbl);
}
// Net positive / negative charge carried by the compact ion set (the same ions
// the compact Saha balances; no He++, higher C/O ions, K, Na or Gr2+).
inline double AppChargePlus(const std::array<double, kNSp>& y) {
  return y[AppSp::Hp] + y[AppSp::H2p] + y[AppSp::H3p] + y[AppSp::Hep] +
         y[AppSp::Dp] + y[AppSp::Cp] + y[AppSp::Op] + y[AppSp::OHp] +
         y[AppSp::H2Op] + y[AppSp::HCOp] + y[AppSp::H3Op] + y[AppSp::Lip] +
         y[AppSp::Mgp] + y[AppSp::Grp];
}
inline double AppChargeMinus(const std::array<double, kNSp>& y) {
  return y[AppSp::e] + y[AppSp::Hm] + y[AppSp::Grm] + 2.0 * y[AppSp::Gr2m];
}
constexpr const char* kNetworkLabel =
    "metal_grain_minimal N_sp=40 N_react=113 "
    "(Nakauchi2021 minimal: compact 40-species / 103+10-reaction network)";
constexpr const char* kOutputTag = "_min";
constexpr const char* kSpeciesCsv =
    "H,H2,e-,H+,H2+,H3+,H-,He,He+,D,HD,D+,"
    "C,CH,CH2,C+,O,O2,OH,CO,H2O,HCO,O+,OH+,H2O+,HCO+,H3O+,"
    "Li,Li+,Mg,Mg+,Gr,Gr+,Gr-,Gr2-,"
    "O(gr),OH(gr),CO(gr),H2O(gr),C(gr)";
#else
using AppCell = arche::MetalCell;
using AppTable = arche::MetalTable;
using AppCellPtr = arche::MetalCellPtr;
using AppTablePtr = arche::MetalTablePtr;
using AppSp = arche::metal_grain::Sp;
constexpr int kNSp = arche::metal_grain::N_sp;  // 89
inline AppCellPtr AppCellCreateOwned() {
  return arche::metal_cell_create_owned();
}
inline arche::ChemStateMG& AppCellState(AppCell& cell) {
  return arche::metal_cell_state(cell);
}
inline AppTablePtr AppTableCreate(const std::string& path) {
  return arche::load_metal_table_owned(path);
}
inline arche::ChemFullRates AppChemFullStep(AppCell& cell, double dt,
                                            const arche::ChemParams& params,
                                            const arche::ChemShielding& shield,
                                            const AppTable& tbl) {
  return arche::chem_full_step_metal(cell, dt, params, shield, tbl);
}
inline int AppTableNInvariants(const AppTable& tbl) {
  return arche::metal_table_n_invariants(tbl);
}
inline double AppChargePlus(const std::array<double, kNSp>& y) {
  return y[3] + y[4] + y[5] + y[8] + y[10] + y[13] + y[14] + y[22] + y[23] +
         y[24] + y[25] + y[26] + y[27] + y[28] + y[39] + y[40] + y[41] + y[42] +
         y[43] + y[44] + y[45] + y[46] + y[47] + y[48] + y[49] + y[52] + y[54] +
         y[58] + y[60] + y[62] + y[64] + 2.0 * (y[9] + y[55] + y[65]) +
         3.0 * y[56];
}
inline double AppChargeMinus(const std::array<double, kNSp>& y) {
  return y[2] + y[6] + y[15] + y[53] + y[66] + 2.0 * y[67];
}
constexpr const char* kNetworkLabel = "metal_grain N_sp=89 N_react=1200";
constexpr const char* kOutputTag = "";
constexpr const char* kSpeciesCsv =
    "H,H2,e-,H+,H2+,H3+,H-,He,He+,He++,HeH+,"
    "D,HD,D+,HD+,D-,"
    "C,C2,CH,CH2,CH3,CH4,C+,C2+,CH+,CH2+,CH3+,CH4+,CH5+,"
    "O,O2,OH,CO,H2O,HCO,HO2,CO2,H2CO,H2O2,"
    "O+,O2+,OH+,CO+,H2O+,HCO+,O2H+,H3O+,H2CO+,HOCO+,H2COH+,"
    "Li,LiH,Li+,Li-,LiH+,Li++,Li+++,"
    "K,K+,Na,Na+,Mg,Mg+,"
    "Gr,Gr+,Gr2+,Gr-,Gr2-,"
    "H(gr),H(ch),H2(gr),D(gr),D(ch),HD(gr),"
    "O(gr),O2(gr),OH(gr),CO(gr),CO2(gr),H2O(gr),"
    "HO2(gr),H2O2(gr),HCO(gr),H2CO(gr),"
    "C(gr),CH(gr),CH2(gr),CH3(gr),CH4(gr)";
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
using collapse_defaults::kTHighStop;
using collapse_defaults::kXnHStop;

// ─── CR attenuation constants (metal-specific)
// ──────────────────────────────── Secondary CR fraction (Padovani et al.)
constexpr double kCrAttenuSecondFrac = 7.6e-2;
// Metal background CR floor [erg s^-1 g^-1 / Z_sun]
constexpr double kCrMetalBkgnd = 1.4e-22;
// Short-lived radionuclide ionization scaling (default off)
constexpr double kSraRateDefault = 0.0;
// Long-lived radionuclide ionization scaling (default off)
constexpr double kLraRateDefault = 0.0;
// Effective CR desorption spike temperature [K]
constexpr double kTcrDesorp = 70.0;

// ─── MRN grain size distribution parameters ──────────────────────────────────
// Bipartite power-law: n(a) ∝ a^{-3.5} (a_min..a_mid), a^{-5.5} (a_mid..a_max)
constexpr double kGrAMin = arche::metal_grain::mrn_a_min;
constexpr double kGrAMid = arche::metal_grain::mrn_a_mid;
constexpr double kGrAMax = arche::metal_grain::mrn_a_max;
constexpr double kGrInd1 = arche::metal_grain::mrn_ind1;
constexpr double kGrInd2 = arche::metal_grain::mrn_ind2;
constexpr double kGrNorm = arche::metal_grain::mrn_norm_zsun;

// ─── Initial gas-phase element fractions ─────────────────────────────────────
// (1 - fraction) is initially locked in grain mantles
constexpr double kCgasFrac = arche::metal_grain::c_gas_frac_default;
constexpr double kOgasFrac = arche::metal_grain::o_gas_frac_default;
constexpr double kMggasFrac = arche::metal_grain::mg_gas_frac_default;

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
using collapse_defaults::kZmetal;
using collapse_defaults::kZred;

// ─────────────────────────────────────────────────────────────────────────────
// OutputRow — one record (every 10 steps) buffered before HDF5 write
// ─────────────────────────────────────────────────────────────────────────────
struct OutputRow {
  int step;
  std::array<double, kNSp> y;
  double nH, T_K, T_gr_K, rho;
  double Lambda_net, Lambda_line, Lambda_cnt, Lambda_gr, Lambda_gas;
  double Lambda_chem, Gamma_cmp;
  double Lambda_H2, Lambda_HD, Lambda_Lya;
  double Lambda_CO, Lambda_OH, Lambda_H2O;
  double Lambda_CII, Lambda_CI, Lambda_OI;
  double Gamma_CR;
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
                   const std::string& z_tag, const std::vector<OutputRow>& rows,
                   double zeta0, double Z_metal, double f_ret, double sra_rate,
                   double lra_rate, double T_rad, double zred, double T_K0_ic,
                   double y_e0_ic, double y_H2_ic, double y_HD_ic,
                   double jlw21 = 0.0, const std::string& fret_table_path = "",
                   bool ff_gamma = false) {
  const hsize_t N = static_cast<hsize_t>(rows.size());
  const hsize_t Nsp = static_cast<hsize_t>(kNSp);

  // ── Root attributes ───────────────────────────────────────────────────────
  {
    char desc[160];
    std::snprintf(desc, sizeof(desc),
                  "metal_grain collapse CR%s Z%s (C++ port)", cr_tag.c_str(),
                  z_tag.c_str());
    H5WriteStrAttr(fid, "description", std::string(desc));
  }
  // Output schema version (N-B). v1 = legacy x-prefixed dataset names
  // (nH/Lambda_*/Gamma_*/MJ/lmbd_J); v2 = renamed (nH/Lambda_*/Gamma_*/M_J/
  // lambda_J). Post-processing tools branch on this attribute.
  H5WriteIntAttr(fid, "schema_version", 2);
  H5WriteStrAttr(fid, "cr_tag", cr_tag);
  H5WriteStrAttr(fid, "z_tag", z_tag);
  H5WriteDblAttr(fid, "zeta0_cgs", zeta0);
  H5WriteDblAttr(fid, "Z_metal", Z_metal);
  H5WriteDblAttr(fid, "f_ret", f_ret);
  H5WriteDblAttr(fid, "sra_rate", sra_rate);
  H5WriteDblAttr(fid, "lra_rate", lra_rate);
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
  dump1d("T_gr_K", [](const OutputRow& r) { return r.T_gr_K; });
  dump1d("rho", [](const OutputRow& r) { return r.rho; });
  dump1d("Lambda_net", [](const OutputRow& r) { return r.Lambda_net; });
  dump1d("Lambda_line", [](const OutputRow& r) { return r.Lambda_line; });
  dump1d("Lambda_cnt", [](const OutputRow& r) { return r.Lambda_cnt; });
  dump1d("Lambda_gr", [](const OutputRow& r) { return r.Lambda_gr; });
  dump1d("Lambda_gas", [](const OutputRow& r) { return r.Lambda_gas; });
  dump1d("Lambda_chem", [](const OutputRow& r) { return r.Lambda_chem; });
  dump1d("Gamma_cmp", [](const OutputRow& r) { return r.Gamma_cmp; });
  dump1d("Lambda_H2", [](const OutputRow& r) { return r.Lambda_H2; });
  dump1d("Lambda_HD", [](const OutputRow& r) { return r.Lambda_HD; });
  dump1d("Lambda_Lya", [](const OutputRow& r) { return r.Lambda_Lya; });
  dump1d("Lambda_CO", [](const OutputRow& r) { return r.Lambda_CO; });
  dump1d("Lambda_OH", [](const OutputRow& r) { return r.Lambda_OH; });
  dump1d("Lambda_H2O", [](const OutputRow& r) { return r.Lambda_H2O; });
  dump1d("Lambda_CII", [](const OutputRow& r) { return r.Lambda_CII; });
  dump1d("Lambda_CI", [](const OutputRow& r) { return r.Lambda_CI; });
  dump1d("Lambda_OI", [](const OutputRow& r) { return r.Lambda_OI; });
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
// Writes:  <out_dir>/collapse_CR<cr>_Z<z>.h5
//          <out_dir>/collapse_CR<cr>_Z<z>_fret<fr>.h5          (f_ret != 1.0)
//          <out_dir>/collapse_CR<cr>_Z<z>_fret<fr>_JLW<jlw>.h5  (+ J_LW21 != 0)
//          <out_dir>/collapse_CR<cr>_Z<z>_fret-gamma.h5        (ff_gamma mode)
//          suffix order: _fret → _JLW → _z
//
// Collapse timescale modes (selected by ff_gamma and fret_tag):
//   (A/B) ff_gamma=false:  t_eff = f_ret * t_ff  (fixed or step-function f_ret)
//   (C)   ff_gamma=true:   t_eff = t_ff / sqrt(1-f(γ))  (Higuchi+2018 Eq.5-7)
// Affects: density update (drho), compressional heating (Γ_cmp),
// shielding lengths (lsh), and the integration timestep (dt).
// ─────────────────────────────────────────────────────────────────────────────
void RunCollapse(
    double zeta0, const std::string& cr_tag, double Z_metal,
    const std::string& z_tag, double f_ret_init, const std::string& fret_tag,
    double T_rad, double zred, const std::string& zred_tag,
    const arche::abundance::MetalSet& abund, const AppTable& tbl,
    const std::string& out_dir, double jlw21 = 0.0,
    const std::string& jlw_tag = "",
    double cr_atten_col_dens = kCrAttenuColDens,
    double cr_atten_second_frac = kCrAttenuSecondFrac,
    double cr_metal_bkgnd = kCrMetalBkgnd, double sra_rate = kSraRateDefault,
    double lra_rate = kLraRateDefault, double t_cr_desorp = kTcrDesorp,
    double nH0 = 1.0, double T_K0 = 1.0e2, double y_e0 = 1.0e-4,
    double y_H2_init = 6.0e-7, double y_HD_init = 4.0e-10,
    double c_gas_frac = kCgasFrac, double o_gas_frac = kOgasFrac,
    double mg_gas_frac = kMggasFrac, const std::vector<double>& fret_nH = {},
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
  // ── Metallicity-scaled abundances ─────────────────────────────────────────
  const double XC = abund.yC * Z_metal;
  const double XO = abund.yO * Z_metal;
#ifndef ARCHE_METAL_MINIMAL
  const double XKa = abund.yK * Z_metal;   // K / Na are absent from the compact
  const double XNa = abund.yNa * Z_metal;  // network
#endif
  const double XMg = abund.yMg * Z_metal;

  // ── Initial conditions ────────────────────────────────────────────────────
  const double nH1 = nH0;   // initial H number density [cm^-3]
  const double kT1 = T_K0;  // initial temperature [K]

  AppCellPtr cell = AppCellCreateOwned();
  auto& st = AppCellState(*cell);
  auto& y = st.y;

  const double y_H2 = y_H2_init;
  const double y_Dp = 0.0;
  const double y_HD = y_HD_init;
  const double y_Hp = y_e0;
  const double y_H = 1.0 - y_Hp - 2.0 * y_H2 - y_HD;
  if (y_H <= 0.0) {
    std::fprintf(stderr,
                 "ERROR: Invalid IC: y_H = 1 - y_e0(%.3g) - 2*y_H2(%.3g) - "
                 "y_HD(%.3g) = %.3g <= 0\n",
                 y_Hp, y_H2, y_HD, y_H);
    std::exit(1);
  }
  const double y_D = abund.yD - y_Dp - y_HD;

  y[AppSp::H] = y_H;              // H
  y[AppSp::H2] = y_H2;            // H2
  y[AppSp::Hp] = y_Hp;            // H+
  y[AppSp::He] = abund.yHe;       // He
  y[AppSp::D] = y_D;              // D
  y[AppSp::HD] = y_HD;            // HD
  y[AppSp::Dp] = y_Dp;            // D+
  y[AppSp::Lip] = abund.yLi;      // Li+
  y[AppSp::C] = XC * c_gas_frac;  // CI
  y[AppSp::O] = XO * o_gas_frac;  // OI
#ifndef ARCHE_METAL_MINIMAL
  y[AppSp::Kp] = 0.0;   // K+  (initially zero; absent from the compact network)
  y[AppSp::Nap] = 0.0;  // Na+ (initially zero; absent from the compact network)
#endif
  y[AppSp::Mgp] = XMg * mg_gas_frac;  // Mg+
  // e- = sum of positive ions
#ifdef ARCHE_METAL_MINIMAL
  y[AppSp::e] = y_Hp + y_Dp + abund.yLi + y[AppSp::Mgp];
#else
  y[AppSp::e] =
      y_Hp + y_Dp + abund.yLi + y[AppSp::Kp] + y[AppSp::Nap] + y[AppSp::Mgp];
#endif

  // MRN grain abundance: yGr = Z_metal * kGrNorm * (1+4yHe)*mp / <vol>
  // Bipartite power-law: a^{-kGrInd1} (kGrAMin..kGrAMid), a^{-kGrInd2}
  // (kGrAMid..kGrAMax)
  double yGr = 0.0;
  if (Z_metal > 0.0) {
    double sum_ana = std::pow(kGrAMid, kGrInd1) *
                         (std::pow(kGrAMid, 1.0 - kGrInd1) -
                          std::pow(kGrAMin, 1.0 - kGrInd1)) /
                         (1.0 - kGrInd1) +
                     std::pow(kGrAMid, kGrInd2) *
                         (std::pow(kGrAMax, 1.0 - kGrInd2) -
                          std::pow(kGrAMid, 1.0 - kGrInd2)) /
                         (1.0 - kGrInd2);
    double ave_gr_vol = std::pow(kGrAMid, kGrInd1) *
                            (std::pow(kGrAMid, 4.0 - kGrInd1) -
                             std::pow(kGrAMin, 4.0 - kGrInd1)) /
                            (4.0 - kGrInd1) +
                        std::pow(kGrAMid, kGrInd2) *
                            (std::pow(kGrAMax, 4.0 - kGrInd2) -
                             std::pow(kGrAMid, 4.0 - kGrInd2)) /
                            (4.0 - kGrInd2);
    ave_gr_vol = (ave_gr_vol / sum_ana) * (4.0 * kPi / 3.0);
    yGr = Z_metal * kGrNorm * (1.0 + 4.0 * abund.yHe) * kMp / ave_gr_vol;
  }
  y[AppSp::Gr] = yGr;  // Gr (neutral grain)

  // Grain surface depletion reservoirs (elements locked in grains initially)
  double yCs = XC * (1.0 - c_gas_frac);
  double yOs = XO * (1.0 - o_gas_frac);
  double yMgs = XMg * (1.0 - mg_gas_frac);
#ifndef ARCHE_METAL_MINIMAL
  double yKs = XKa;  // K / Na reservoirs are absent from the compact network
  double yNas = XNa;
#endif
  int switch_gr = 1;  // 1 = grains present, 2 = fully evaporated

  // ── Thermodynamic state ───────────────────────────────────────────────────
  double rho = ((1.0 + 4.0 * abund.yHe) * kMp) * nH1;
  double nH = nH1;
  double T_K = kT1;
  double p = (1.0 + abund.yHe) * nH1 * kKB * kT1;

  // Free-particle sums for the EOS.  The full network adds He++ (compact has
  // no He++ species); the += keeps the full sum byte-identical to the original
  // left-to-right fold.
  double mu_denom = y[AppSp::H] + y[AppSp::H2] + y[AppSp::e] + y[AppSp::Hp] +
                    y[AppSp::He] + y[AppSp::Hep];
  double mono =
      y[AppSp::H] + y[AppSp::e] + y[AppSp::Hp] + y[AppSp::He] + y[AppSp::Hep];
#ifndef ARCHE_METAL_MINIMAL
  mu_denom += y[AppSp::Hepp];
  mono += y[AppSp::Hepp];
#endif
  double mu = (1.0 + 4.0 * abund.yHe) / mu_denom;
  double gamma =
      1.0 + (1.0 + 4.0 * abund.yHe) /
                (mu * (1.5 * mono + arche::c_H2(T_K) * y[AppSp::H2]));
  double e = kKB * T_K / ((gamma - 1.0) * mu * kMp);

  double T_gr_K = T_rad;
  double rho_old = rho;
  double T_gr_old = T_gr_K;

  // ── Chemistry parameters ──────────────────────────────────────────────────
  arche::ChemParams params{};
  params.T_rad = T_rad;
  params.zeta = zeta0;
  params.Z_metal = Z_metal;
  params.T_gr_K = T_gr_K;
  params.T_cr_desorp = t_cr_desorp;

  // ── Scalar loop state ─────────────────────────────────────────────────────
  double t = 0.0;
  double dt = 1.0e-1;
  double t_chem = 1.0e-1;
  double t_cool = 0.0;
  double k_gr = 0.0;
  double k_gas = 0.0;
  double tau_cnt = 0.0;
  double esc_cnt = 1.0;

  // Column densities
  double Nc_H = 0.0;
  double Nc_H2 = 0.0;
  double Nc_HD = 0.0;
  double Nc_CII = 0.0;
  double Nc_CI = 0.0;
  double Nc_OI = 0.0;
  double Nc_CO = 0.0;
  double Nc_H2O = 0.0;
  double Nc_OH = 0.0;

  // ── f_ret step-function table state ──────────────────────────────────────
  double f_ret = f_ret_init;
  int fret_idx = 0;
  const bool has_fret_tab = !fret_nH.empty();

  // ── HDF5 row buffer ───────────────────────────────────────────────────────
  std::vector<OutputRow> h5rows;
  h5rows.reserve(6000);

  std::printf("CR%s Z%s: zeta0=%g z=%g yGr=%g\n", cr_tag.c_str(), z_tag.c_str(),
              zeta0 / 1.0e-17, Z_metal, yGr);

  // ── Bench timing CSV ─────────────────────────────────────────────────────
  using Clock = std::chrono::high_resolution_clock;
  std::FILE* bench_fp = nullptr;
  if (bench_mode) {
    std::string bench_path =
        out_dir + "/bench_CR" + cr_tag + "_Z" + z_tag + ".csv";
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
    // Save state from previous step for grain evaporation ratio
    rho_old = rho;
    T_gr_old = T_gr_K;

    if (has_fret_tab)
      collapse_driver::update_fret(fret_idx, f_ret, nH, fret_nH, fret_val);

    // ── Free-fall time & Jeans length ─────────────────────────────────────
    double t_ff = std::sqrt(3.0 * kPi / (32.0 * kGGrav * rho));
    double t_eff = t_eff_collapse(t_ff, f_ret, gamma, ff_gamma);
    double lmbd_J = std::sqrt(kPi * kKB * T_K / (kGGrav * mu * kMp * rho));

    // ── Thermal velocities for column density ─────────────────────────────
    auto vD = [&](double mass_amu) {
      return std::sqrt(2.0 * kKB * T_K / (kMp * mass_amu));
    };
    double lsh_H = std::min(lmbd_J, 6.0 * vD(1.0) * t_eff);
    double lsh_H2 = std::min(lmbd_J, 6.0 * vD(2.0) * t_eff);
    double lsh_HD = std::min(lmbd_J, 6.0 * vD(3.0) * t_eff);
    double lsh_C = std::min(lmbd_J, 6.0 * vD(12.0) * t_eff);
    double lsh_O = std::min(lmbd_J, 6.0 * vD(16.0) * t_eff);
    double lsh_CO = std::min(lmbd_J, 6.0 * vD(28.0) * t_eff);
    double lsh_H2O = std::min(lmbd_J, 6.0 * vD(18.0) * t_eff);
    double lsh_OH = std::min(lmbd_J, 6.0 * vD(17.0) * t_eff);

    Nc_H = y[AppSp::H] * nH * lsh_H;
    Nc_H2 = y[AppSp::H2] * nH * lsh_H2;
    Nc_HD = y[AppSp::HD] * nH * lsh_HD;
    Nc_CII = y[AppSp::Cp] * nH * lsh_C;  // C+
    Nc_CI = y[AppSp::C] * nH * lsh_C;    // C
    Nc_OI = y[AppSp::O] * nH * lsh_O;    // O
    Nc_CO = y[AppSp::CO] * nH * lsh_CO;
    Nc_H2O = y[AppSp::H2O] * nH * lsh_H2O;
    Nc_OH = y[AppSp::OH] * nH * lsh_OH;

    // ── Continuum optical depth ───────────────────────────────────────────
    tau_cnt = (k_gr + k_gas) * rho * lmbd_J;
    if (tau_cnt > 1.0)
      esc_cnt = 1.0 / (tau_cnt * tau_cnt);
    else
      esc_cnt = 1.0;

    // ── Build shielding environment ───────────────────────────────────────
    arche::ChemShielding shield;
    const double zeta_r_short = sra_rate * 7.6e-19;
    const double zeta_r_long = lra_rate * 1.4e-22 * Z_metal;
    if (zeta0 > 0.0)
      shield.zeta = zeta0 * (std::exp(-rho * lmbd_J / cr_atten_col_dens) +
                             cr_atten_second_frac) +
                    cr_metal_bkgnd * Z_metal + zeta_r_short + zeta_r_long;
    else
      shield.zeta = zeta_r_short + zeta_r_long;
    shield.Nc_H = Nc_H;
    shield.Nc_H2 = Nc_H2;
    shield.Nc_HD = Nc_HD;
    shield.Nc_CO = Nc_CO;
    shield.Nc_OH = Nc_OH;
    shield.Nc_H2O = Nc_H2O;
    shield.Nc_CII = Nc_CII;
    shield.Nc_CI = Nc_CI;
    shield.Nc_OI = Nc_OI;
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
    // Grain temperature is now returned; feed it back as the next step's
    // warm-start seed and keep the local copy for diagnostics/output.
    params.T_gr_K = rates.T_gr_K;
    T_gr_K = params.T_gr_K;
    k_gr = rates.k_gr;
    k_gas = rates.k_gas;

    const double Lambda_line = rates.Lambda_line;
    const double Lambda_gr = rates.Lambda_gr;
    const double Lambda_gas_out = rates.Lambda_gas;
    const double Lambda_cnt = rates.Lambda_cnt;
    const double Lambda_ch = rates.Lambda_chem;
    const double Gamma_CR = rates.Gamma_CR;
    const double Lambda_net = rates.Lambda_net;

    t_chem = collapse_driver::compute_tchem<kNSp>(dt, y, y_prev);
    double Gamma_cmp = p / rho / t_eff;

    double MJ, B_cr;
    collapse_driver::compute_diagnostics(rho, lmbd_J, MJ, B_cr);

    // ── Grain evaporation ─────────────────────────────────────────────────
    if (switch_gr == 1) {
      double T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol;
      arche::detail::vaptemp(rho, T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol);

      double vol_cur = arche::detail::vol_gr(rho, T_gr_K);
      double vol_old = arche::detail::vol_gr(rho_old, T_gr_old);

      if (vol_cur < 1.0e-60) {
        // Complete evaporation: release all grain charge and the locked-in
        // metals.  The compact network has no K / Na (and no Gr2+), so its net
        // grain charge goes to the electrons directly and only Mg is released.
#ifdef ARCHE_METAL_MINIMAL
        double del_yp = y[AppSp::Grp];                         // Gr+
        double del_ym = y[AppSp::Grm] + 2.0 * y[AppSp::Gr2m];  // Gr- + 2*Gr2-
        y[AppSp::e] += (del_ym - del_yp);  // net grain charge -> electrons
        y[AppSp::Mg] += yMgs;              // Mg neutral released
        yMgs = 0.0;
        yCs = 0.0;
        yOs = 0.0;
        y[AppSp::Gr] = 0.0;
        y[AppSp::Grp] = 0.0;
        y[AppSp::Grm] = 0.0;
        y[AppSp::Gr2m] = 0.0;
#else
        double del_yp = y[64] + 2.0 * y[65];  // Gr+ + 2*Gr2+
        double del_ym = y[66] + 2.0 * y[67];  // Gr- + 2*Gr2-
        if (del_ym >= del_yp) {
          y[2] += (del_ym - del_yp);  // e-
        } else {
          y[58] += (del_yp - del_ym);  // K+
          y[57] -= (del_yp - del_ym);  // K (neutral)
        }
        y[57] += yKs;   // K neutral released
        y[59] += yNas;  // Na neutral released
        y[61] += yMgs;  // Mg neutral released
        yKs = 0.0;
        yNas = 0.0;
        yMgs = 0.0;
        yCs = 0.0;
        yOs = 0.0;
        y[63] = 0.0;
        y[64] = 0.0;
        y[65] = 0.0;
        y[66] = 0.0;
        y[67] = 0.0;
#endif
        switch_gr = 2;

      } else {
        double gr_frac = (vol_old > 0.0) ? vol_cur / vol_old : 1.0;

        if (T_K >= 0.975 * T_pyr) {
          // Partial evaporation of grain charges and locked-in metals.
#ifdef ARCHE_METAL_MINIMAL
          double del_yp = y[AppSp::Grp] * (1.0 - gr_frac);
          double del_ym =
              (y[AppSp::Grm] + 2.0 * y[AppSp::Gr2m]) * (1.0 - gr_frac);
          y[AppSp::e] += (del_ym - del_yp);
          y[AppSp::Mg] += yMgs * (1.0 - gr_frac);
          yCs *= gr_frac;
          yOs *= gr_frac;
          yMgs *= gr_frac;
          y[AppSp::Grp] *= gr_frac;
          y[AppSp::Grm] *= gr_frac;
          y[AppSp::Gr2m] *= gr_frac;
#else
          double del_yp = (y[64] + 2.0 * y[65]) * (1.0 - gr_frac);
          double del_ym = (y[66] + 2.0 * y[67]) * (1.0 - gr_frac);
          if (del_ym >= del_yp) {
            y[2] += (del_ym - del_yp);
          } else {
            y[58] += (del_yp - del_ym);
            y[57] -= (del_yp - del_ym);
          }
          y[57] += yKs * (1.0 - gr_frac);
          y[59] += yNas * (1.0 - gr_frac);
          y[61] += yMgs * (1.0 - gr_frac);
          yCs *= gr_frac;
          yOs *= gr_frac;
          yKs *= gr_frac;
          yNas *= gr_frac;
          yMgs *= gr_frac;
          y[64] *= gr_frac;
          y[65] *= gr_frac;
          y[66] *= gr_frac;
          y[67] *= gr_frac;
#endif
        }
        y[AppSp::Gr] *= gr_frac;  // Gr neutral scales with grain volume
      }
    }

    // ── Charge accounting ─────────────────────────────────────────────────
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
      row.T_gr_K = T_gr_K;
      row.rho = rho;
      row.Lambda_net = Lambda_net;
      row.Lambda_line = Lambda_line;
      row.Lambda_cnt = Lambda_cnt;
      row.Lambda_gr = Lambda_gr;
      row.Lambda_gas = Lambda_gas_out;
      row.Lambda_chem = Lambda_ch;
      row.Gamma_cmp = Gamma_cmp;
      row.Lambda_H2 = rates.Lambda_H2;
      row.Lambda_HD = rates.Lambda_HD;
      row.Lambda_Lya = rates.Lambda_Lya;
      row.Lambda_CO = rates.Lambda_CO;
      row.Lambda_OH = rates.Lambda_OH;
      row.Lambda_H2O = rates.Lambda_H2O;
      row.Lambda_CII = rates.Lambda_CII;
      row.Lambda_CI = rates.Lambda_CI;
      row.Lambda_OI = rates.Lambda_OI;
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
                                                rho, e, T_K, p, nH, t, true)) {
      exit_reason = ExitReason::NegEnergy;
      break;
    }

    t_cool = (Lambda_net != 0.0) ? e / std::abs(Lambda_net) : 1.0e50;
    dt = collapse_driver::compute_timestep(it, n_init_steps, t_cool, t_eff,
                                           t_chem, dt_factor, dt_factor_chem,
                                           dt_factor_init);

    // ── Progress ──────────────────────────────────────────────────────────
    std::printf(
        "%7d %11.3E %11.3E %11.3E %11.3E %11.3E %11.3E"
        " %11.3E %11.3E %11.3E\n",
        it, zeta0 / 1.0e-17, Z_metal, nH, T_K, y[AppSp::e], y[AppSp::H2],
        y_plus, y_minus, charge_imbal);

    if (collapse_driver::check_exit(nH, T_K, e, nH_stop, exit_reason)) break;
  }

  if (bench_fp) {
    std::fclose(bench_fp);
    bench_fp = nullptr;
  }

  // ── Exit status report ────────────────────────────────────────────────────
  std::string exit_label = "metal collapse CR" + cr_tag + " Z" + z_tag;
  int exit_code = collapse_driver::report_exit(exit_reason, exit_label.c_str());
  const char* exit_msg = collapse_driver::exit_message(exit_reason);
  collapse_driver::report_conservation(cons_tally, cons_rows);

  // ── Write HDF5 file ───────────────────────────────────────────────────────
  std::string h5_path = out_dir + "/collapse_CR" + cr_tag + "_Z" + z_tag;
  if (!fret_tag.empty()) h5_path += "_fret" + fret_tag;
  if (!jlw_tag.empty()) h5_path += "_JLW" + jlw_tag;
#ifdef ARCHE_XRAY
  if (!zx_tag.empty()) h5_path += "_ZX" + zx_tag;
#endif
  if (!zred_tag.empty()) h5_path += "_z" + zred_tag;
  h5_path += kOutputTag;  // "_min" for the compact metal variant
  h5_path += ".h5";
  hid_t fid = H5Create(h5_path);
  if (fid < 0) {
    std::fprintf(stderr, "ERROR: cannot create %s\n", h5_path.c_str());
    return;
  }
  WriteHdf5File(fid, cr_tag, z_tag, h5rows, zeta0, Z_metal, f_ret_init,
                sra_rate, lra_rate, T_rad, zred, T_K0, y_e0, y_H2_init,
                y_HD_init, jlw21, fret_table_path, ff_gamma);
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
  const char* env_chem_table = std::getenv("METAL_CHEM_TABLE");
  std::string chem_table = (env_chem_table && env_chem_table[0] != '\0')
                               ? env_chem_table
                               : METAL_CHEM_TABLE;

  // ── Output directory ──────────────────────────────────────────────────────
  const char* env_out = std::getenv("METAL_OUTDIR");
  std::string out_dir = env_out ? env_out : "results/metal/h5";
  std::filesystem::create_directories(out_dir);
  std::printf("Output directory: %s/\n", out_dir.c_str());

  // ── Load reaction table from HDF5 (via facade; opaque handle) ─────────────
  AppTablePtr tbl = AppTableCreate(chem_table);

  // ── Required: METAL_ZETA0 [s^-1] ─────────────────────────────────────────
  const char* env_zeta0 = std::getenv("METAL_ZETA0");
  if (!env_zeta0 || env_zeta0[0] == '\0') {
    std::fprintf(
        stderr,
        "ERROR: environment variable METAL_ZETA0 is required\n"
        "       Set the CR ionization rate [s^-1], e.g.: METAL_ZETA0=1e-17\n");
    return 1;
  }
  double zeta0 = kZeta0;
  try {
    zeta0 = std::stod(env_zeta0);
  } catch (const std::exception&) {
    std::fprintf(stderr, "ERROR: METAL_ZETA0='%s' is not a valid number\n",
                 env_zeta0);
    return 1;
  }
  if (zeta0 < 0.0) {
    std::fprintf(stderr, "ERROR: METAL_ZETA0 must be >= 0, got %s\n",
                 env_zeta0);
    return 1;
  }

  // ── Required: METAL_Z_METAL [Z_sun] ──────────────────────────────────────
  const char* env_z = std::getenv("METAL_Z_METAL");
  if (!env_z || env_z[0] == '\0') {
    std::fprintf(
        stderr,
        "ERROR: environment variable METAL_Z_METAL is required\n"
        "       Set the metallicity [Z_sun], e.g.: METAL_Z_METAL=1e-3\n");
    return 1;
  }
  double Z_metal = kZmetal;
  try {
    Z_metal = std::stod(env_z);
  } catch (const std::exception&) {
    std::fprintf(stderr, "ERROR: METAL_Z_METAL='%s' is not a valid number\n",
                 env_z);
    return 1;
  }
  if (Z_metal < 0.0) {
    std::fprintf(stderr, "ERROR: METAL_Z_METAL must be >= 0, got %s\n", env_z);
    return 1;
  }

  std::string cr_tag = collapse_defaults::make_tag(env_zeta0, zeta0);
  std::string z_tag = collapse_defaults::make_tag(env_z, Z_metal);

  // ── Optional: METAL_REDSHIFT (cosmological redshift, default 0.0) ─────────
  // T_rad = kTCMB0 * (1 + z).  Affects line/continuum/molecular cooling,
  // grain temperature equilibrium, and initial grain temperature.
  double zred = kZred;
  std::string zred_tag;  // empty → no suffix in filename
  const char* env_zred = std::getenv("METAL_REDSHIFT");
  if (env_zred && env_zred[0] != '\0') {
    try {
      zred = std::stod(env_zred);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_REDSHIFT='%s' is not a valid number\n",
                   env_zred);
      return 1;
    }
    if (zred < 0.0) {
      std::fprintf(stderr, "ERROR: METAL_REDSHIFT must be >= 0, got %s\n",
                   env_zred);
      return 1;
    }
    if (zred != 0.0) {
      zred_tag = collapse_defaults::make_tag(env_zred, zred);
    }
  }
  const double T_rad = kTCMB0 * (1.0 + zred);

  // ── Optional: METAL_FRET_TABLE (step-function f_ret table; overrides
  // METAL_FF_RET) ─ ASCII 2-col file: nH[cm^-3]  f_ret  (sorted by ascending
  // nH, '#' comments allowed)
  const char* env_fret_table = std::getenv("METAL_FRET_TABLE");
  std::string fret_table_path;
  std::vector<double> fret_nH, fret_val;
  if (env_fret_table && env_fret_table[0] != '\0') {
    fret_table_path = env_fret_table;
    auto [nh, fv] =
        collapse_driver::LoadFretTable(fret_table_path, "METAL_FRET_TABLE");
    fret_nH = std::move(nh);
    fret_val = std::move(fv);
  }

  // ── Optional: METAL_FF_GAMMA (gamma-dependent collapse factor; highest
  // priority) ── When set (non-empty, non-"0"), uses t_eff = t_ff /
  // sqrt(1-f(γ))  [Higuchi+2018 Eq.5-7]. Overrides METAL_FF_RET and
  // METAL_FRET_TABLE.
  const char* env_ff_gamma = std::getenv("METAL_FF_GAMMA");
  const bool ff_gamma =
      (env_ff_gamma && env_ff_gamma[0] != '\0' && env_ff_gamma[0] != '0');

  // ── Optional: METAL_FF_RET (free-fall retardation factor, default 1.0) ───
  // f_ret > 1 slows the collapse; f_ret = 1 gives standard free-fall.
  // Ignored when METAL_FRET_TABLE or METAL_FF_GAMMA is set.
  const char* env_fret = std::getenv("METAL_FF_RET");
  double f_ret = kFFRet;
  std::string fret_tag;  // empty → no suffix in filename
  if (ff_gamma) {
    // Gamma mode: tag "-gamma"; warn if other fret options were also supplied
    fret_tag = "-gamma";
    if (!fret_table_path.empty())
      std::fprintf(
          stderr,
          "WARNING: METAL_FF_GAMMA is set; METAL_FRET_TABLE='%s' is ignored\n",
          env_fret_table);
    if (env_fret && env_fret[0] != '\0')
      std::fprintf(
          stderr,
          "WARNING: METAL_FF_GAMMA is set; METAL_FF_RET='%s' is ignored\n",
          env_fret);
  } else if (!fret_table_path.empty()) {
    // Table mode: initial f_ret from first table entry; fixed tag "-step"
    // → output filename: collapse_CR<cr>_Z<z>_fret-step.h5
    f_ret = fret_val[0];
    fret_tag = "-step";
    if (env_fret && env_fret[0] != '\0')
      std::fprintf(
          stderr,
          "WARNING: METAL_FRET_TABLE is set; METAL_FF_RET='%s' is ignored\n",
          env_fret);
  } else if (env_fret && env_fret[0] != '\0') {
    try {
      f_ret = std::stod(env_fret);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_FF_RET='%s' is not a valid number\n",
                   env_fret);
      return 1;
    }
    if (f_ret <= 0.0) {
      std::fprintf(stderr, "ERROR: METAL_FF_RET must be > 0, got %s\n",
                   env_fret);
      return 1;
    }
    if (f_ret != 1.0) {
      fret_tag = collapse_defaults::make_tag(env_fret, f_ret);
    }
  }

  // ── Optional: METAL_JLW21 (Lyman-Werner intensity J_21, default 0.0) ─────
  // J_21 units: 10^-21 erg/s/cm^2/Hz/sr.  0.0 = no LW field (default).
  double jlw21 = kJLW21;
  const char* env_jlw21 = std::getenv("METAL_JLW21");
  if (env_jlw21 && env_jlw21[0] != '\0') {
    try {
      jlw21 = std::stod(env_jlw21);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_JLW21='%s' is not a valid number\n",
                   env_jlw21);
      return 1;
    }
    if (jlw21 < 0.0) {
      std::fprintf(stderr, "ERROR: METAL_JLW21 must be >= 0, got %s\n",
                   env_jlw21);
      return 1;
    }
  }
  // Tag: "." → "p"; zero → empty (no suffix)
  std::string jlw_tag;
  if (jlw21 != 0.0) {
    jlw_tag = collapse_defaults::make_tag(env_jlw21, jlw21);
  }

  // ── Optional: METAL_ABUNDANCE_SET (default: solar) ──────────────────────
  const char* env_abund = std::getenv("METAL_ABUNDANCE_SET");
  const std::string abundance_set =
      (env_abund && env_abund[0] != '\0') ? std::string(env_abund) : "solar";
  arche::abundance::MetalSet abund{};
  try {
    abund = arche::abundance::get_metal_set(abundance_set);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }

  // ── Optional: METAL_XNH0 [cm^-3] (initial H number density) ─────────────
  double nH0 = kXnH0;
  const char* env_xnH0 = std::getenv("METAL_XNH0");
  if (env_xnH0 && env_xnH0[0] != '\0') {
    try {
      nH0 = std::stod(env_xnH0);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_XNH0='%s' is not a valid number\n",
                   env_xnH0);
      return 1;
    }
    if (nH0 <= 0.0) {
      std::fprintf(stderr, "ERROR: METAL_XNH0 must be > 0, got %s\n", env_xnH0);
      return 1;
    }
  }

  // ── Optional: METAL_TK0 [K] (initial gas temperature) ────────────────────
  double T_K0 = kTK0;
  const char* env_T_K0 = std::getenv("METAL_TK0");
  if (env_T_K0 && env_T_K0[0] != '\0') {
    try {
      T_K0 = std::stod(env_T_K0);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_TK0='%s' is not a valid number\n",
                   env_T_K0);
      return 1;
    }
    if (T_K0 <= 0.0) {
      std::fprintf(stderr, "ERROR: METAL_TK0 must be > 0, got %s\n", env_T_K0);
      return 1;
    }
  }

  // ── Optional: METAL_YE0 (initial electron / H+ fraction) ─────────────────
  double ye0 = kYe0;
  const char* env_ye0 = std::getenv("METAL_YE0");
  if (env_ye0 && env_ye0[0] != '\0') {
    try {
      ye0 = std::stod(env_ye0);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_YE0='%s' is not a valid number\n",
                   env_ye0);
      return 1;
    }
    if (ye0 < 0.0 || ye0 >= 1.0) {
      std::fprintf(stderr, "ERROR: METAL_YE0 must be in [0, 1), got %s\n",
                   env_ye0);
      return 1;
    }
  }

  // ── Optional: METAL_YH2 (initial H2 fraction) ────────────────────────────
  double y_H2_init = kYH2;
  const char* env_yH2 = std::getenv("METAL_YH2");
  if (env_yH2 && env_yH2[0] != '\0') {
    try {
      y_H2_init = std::stod(env_yH2);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_YH2='%s' is not a valid number\n",
                   env_yH2);
      return 1;
    }
    if (y_H2_init < 0.0 || y_H2_init >= 0.5) {
      std::fprintf(stderr, "ERROR: METAL_YH2 must be in [0, 0.5), got %s\n",
                   env_yH2);
      return 1;
    }
  }

  // ── Optional: METAL_YHD (initial HD fraction) ────────────────────────────
  double y_HD_init = kYHD;
  const char* env_yHD = std::getenv("METAL_YHD");
  if (env_yHD && env_yHD[0] != '\0') {
    try {
      y_HD_init = std::stod(env_yHD);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_YHD='%s' is not a valid number\n",
                   env_yHD);
      return 1;
    }
    if (y_HD_init < 0.0) {
      std::fprintf(stderr, "ERROR: METAL_YHD must be >= 0, got %s\n", env_yHD);
      return 1;
    }
  }

  // ── Optional: METAL_OUTPUT_STRIDE (write every N-th step to HDF5) ────────
  int output_stride = kOutputStride;  // default: 10
  const char* env_stride = std::getenv("METAL_OUTPUT_STRIDE");
  if (env_stride && env_stride[0] != '\0') {
    try {
      output_stride = std::stoi(env_stride);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: METAL_OUTPUT_STRIDE='%s' is not a valid integer\n",
                   env_stride);
      return 1;
    }
    if (output_stride <= 0) {
      std::fprintf(stderr, "ERROR: METAL_OUTPUT_STRIDE must be > 0, got %s\n",
                   env_stride);
      return 1;
    }
  }

  // ── Optional: METAL_MAX_ITER (maximum integration steps) ─────────────────
  int max_iter = kItMax;  // default: 1000000
  const char* env_max_iter = std::getenv("METAL_MAX_ITER");
  if (env_max_iter && env_max_iter[0] != '\0') {
    try {
      max_iter = std::stoi(env_max_iter);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: METAL_MAX_ITER='%s' is not a valid integer\n",
                   env_max_iter);
      return 1;
    }
    if (max_iter <= 0) {
      std::fprintf(stderr, "ERROR: METAL_MAX_ITER must be > 0, got %s\n",
                   env_max_iter);
      return 1;
    }
  }

  // ── Optional: METAL_CR_ATTEN_COL_DENS [g cm^-2] ─────────────────────────
  double cr_atten_col_dens = kCrAttenuColDens;
  const char* env_cr_col = std::getenv("METAL_CR_ATTEN_COL_DENS");
  if (env_cr_col && env_cr_col[0] != '\0') {
    try {
      cr_atten_col_dens = std::stod(env_cr_col);
    } catch (const std::exception&) {
      std::fprintf(
          stderr, "ERROR: METAL_CR_ATTEN_COL_DENS='%s' is not a valid number\n",
          env_cr_col);
      return 1;
    }
    if (cr_atten_col_dens <= 0.0) {
      std::fprintf(stderr,
                   "ERROR: METAL_CR_ATTEN_COL_DENS must be > 0, got %s\n",
                   env_cr_col);
      return 1;
    }
  }

  // ── Optional: METAL_CR_ATTEN_SECOND_FRAC ────────────────────────────────
  double cr_atten_second_frac = kCrAttenuSecondFrac;
  const char* env_cr_second = std::getenv("METAL_CR_ATTEN_SECOND_FRAC");
  if (env_cr_second && env_cr_second[0] != '\0') {
    try {
      cr_atten_second_frac = std::stod(env_cr_second);
    } catch (const std::exception&) {
      std::fprintf(
          stderr,
          "ERROR: METAL_CR_ATTEN_SECOND_FRAC='%s' is not a valid number\n",
          env_cr_second);
      return 1;
    }
    if (cr_atten_second_frac < 0.0) {
      std::fprintf(stderr,
                   "ERROR: METAL_CR_ATTEN_SECOND_FRAC must be >= 0, got %s\n",
                   env_cr_second);
      return 1;
    }
  }

  // ── Optional: METAL_CR_METAL_BKGND [erg s^-1 g^-1 / Z_sun] ──────────────
  double cr_metal_bkgnd = kCrMetalBkgnd;
  const char* env_cr_bkg = std::getenv("METAL_CR_METAL_BKGND");
  if (env_cr_bkg && env_cr_bkg[0] != '\0') {
    try {
      cr_metal_bkgnd = std::stod(env_cr_bkg);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: METAL_CR_METAL_BKGND='%s' is not a valid number\n",
                   env_cr_bkg);
      return 1;
    }
    if (cr_metal_bkgnd < 0.0) {
      std::fprintf(stderr, "ERROR: METAL_CR_METAL_BKGND must be >= 0, got %s\n",
                   env_cr_bkg);
      return 1;
    }
  }

  // ── Optional: METAL_SRA_RATE / METAL_LRA_RATE (radionuclides) ───────────
  double sra_rate = kSraRateDefault;
  double lra_rate = kLraRateDefault;
  const char* env_sra_rate = std::getenv("METAL_SRA_RATE");
  if (env_sra_rate && env_sra_rate[0] != '\0') {
    try {
      sra_rate = std::stod(env_sra_rate);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_SRA_RATE='%s' is not a valid number\n",
                   env_sra_rate);
      return 1;
    }
    if (sra_rate < 0.0) {
      std::fprintf(stderr, "ERROR: METAL_SRA_RATE must be >= 0, got %s\n",
                   env_sra_rate);
      return 1;
    }
  }
  const char* env_lra_rate = std::getenv("METAL_LRA_RATE");
  if (env_lra_rate && env_lra_rate[0] != '\0') {
    try {
      lra_rate = std::stod(env_lra_rate);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_LRA_RATE='%s' is not a valid number\n",
                   env_lra_rate);
      return 1;
    }
    if (lra_rate < 0.0) {
      std::fprintf(stderr, "ERROR: METAL_LRA_RATE must be >= 0, got %s\n",
                   env_lra_rate);
      return 1;
    }
  }

  // ── Optional: METAL_T_CR_DES [K] (advanced) ─────────────────────────────
  double t_cr_desorp = kTcrDesorp;
  const char* env_tcr = std::getenv("METAL_T_CR_DES");
  if (env_tcr && env_tcr[0] != '\0') {
    try {
      t_cr_desorp = std::stod(env_tcr);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_T_CR_DES='%s' is not a valid number\n",
                   env_tcr);
      return 1;
    }
    if (t_cr_desorp <= 0.0) {
      std::fprintf(stderr, "ERROR: METAL_T_CR_DES must be > 0, got %s\n",
                   env_tcr);
      return 1;
    }
  }

  // ── Optional: METAL_C_GAS_FRAC / METAL_O_GAS_FRAC / METAL_MG_GAS_FRAC ───
  double c_gas_frac = kCgasFrac;
  double o_gas_frac = kOgasFrac;
  double mg_gas_frac = kMggasFrac;

  const char* env_cfrac = std::getenv("METAL_C_GAS_FRAC");
  if (env_cfrac && env_cfrac[0] != '\0') {
    try {
      c_gas_frac = std::stod(env_cfrac);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: METAL_C_GAS_FRAC='%s' is not a valid number\n",
                   env_cfrac);
      return 1;
    }
  }
  const char* env_ofrac = std::getenv("METAL_O_GAS_FRAC");
  if (env_ofrac && env_ofrac[0] != '\0') {
    try {
      o_gas_frac = std::stod(env_ofrac);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: METAL_O_GAS_FRAC='%s' is not a valid number\n",
                   env_ofrac);
      return 1;
    }
  }
  const char* env_mgfrac = std::getenv("METAL_MG_GAS_FRAC");
  if (env_mgfrac && env_mgfrac[0] != '\0') {
    try {
      mg_gas_frac = std::stod(env_mgfrac);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: METAL_MG_GAS_FRAC='%s' is not a valid number\n",
                   env_mgfrac);
      return 1;
    }
  }
  if (c_gas_frac < 0.0 || c_gas_frac > 1.0 || o_gas_frac < 0.0 ||
      o_gas_frac > 1.0 || mg_gas_frac < 0.0 || mg_gas_frac > 1.0) {
    std::fprintf(stderr,
                 "ERROR: gas fractions must be in [0,1]: C=%g O=%g Mg=%g\n",
                 c_gas_frac, o_gas_frac, mg_gas_frac);
    return 1;
  }

  // ── Optional: METAL_DT_FACTOR / METAL_DT_FACTOR_CHEM / METAL_DT_FACTOR_INIT
  // / METAL_N_INIT_STEPS ─
  double dt_factor = kDtFactor;
  double dt_factor_chem = kDtFactorChem;
  double dt_factor_init = kDtFactorInit;
  int n_init_steps = kNInitSteps;

  const char* env_dt_factor = std::getenv("METAL_DT_FACTOR");
  if (env_dt_factor && env_dt_factor[0] != '\0') {
    try {
      dt_factor = std::stod(env_dt_factor);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: METAL_DT_FACTOR='%s' is not a valid number\n",
                   env_dt_factor);
      return 1;
    }
    if (dt_factor <= 0.0) {
      std::fprintf(stderr, "ERROR: METAL_DT_FACTOR must be > 0, got %s\n",
                   env_dt_factor);
      return 1;
    }
  }
  const char* env_dt_factor_chem = std::getenv("METAL_DT_FACTOR_CHEM");
  if (env_dt_factor_chem && env_dt_factor_chem[0] != '\0') {
    try {
      dt_factor_chem = std::stod(env_dt_factor_chem);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: METAL_DT_FACTOR_CHEM='%s' is not a valid number\n",
                   env_dt_factor_chem);
      return 1;
    }
    if (dt_factor_chem <= 0.0) {
      std::fprintf(stderr, "ERROR: METAL_DT_FACTOR_CHEM must be > 0, got %s\n",
                   env_dt_factor_chem);
      return 1;
    }
  }
  const char* env_dt_factor_init = std::getenv("METAL_DT_FACTOR_INIT");
  if (env_dt_factor_init && env_dt_factor_init[0] != '\0') {
    try {
      dt_factor_init = std::stod(env_dt_factor_init);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: METAL_DT_FACTOR_INIT='%s' is not a valid number\n",
                   env_dt_factor_init);
      return 1;
    }
    if (dt_factor_init <= 0.0) {
      std::fprintf(stderr, "ERROR: METAL_DT_FACTOR_INIT must be > 0, got %s\n",
                   env_dt_factor_init);
      return 1;
    }
  }
  const char* env_n_init = std::getenv("METAL_N_INIT_STEPS");
  if (env_n_init && env_n_init[0] != '\0') {
    try {
      n_init_steps = std::stoi(env_n_init);
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "ERROR: METAL_N_INIT_STEPS='%s' is not a valid integer\n",
                   env_n_init);
      return 1;
    }
    if (n_init_steps < 0) {
      std::fprintf(stderr, "ERROR: METAL_N_INIT_STEPS must be >= 0, got %s\n",
                   env_n_init);
      return 1;
    }
  }

  // ── Optional: METAL_XNH_STOP [cm^-3] ─────────────────────────────────────
  double nH_stop = kXnHStop;
  const char* env_xnh_stop = std::getenv("METAL_XNH_STOP");
  if (env_xnh_stop && env_xnh_stop[0] != '\0') {
    try {
      nH_stop = std::stod(env_xnh_stop);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_XNH_STOP='%s' is not a valid number\n",
                   env_xnh_stop);
      return 1;
    }
    if (nH_stop <= 0.0) {
      std::fprintf(stderr, "ERROR: METAL_XNH_STOP must be > 0, got %s\n",
                   env_xnh_stop);
      return 1;
    }
  }

#ifdef ARCHE_XRAY
  // ── Optional: METAL_ZETA_X [s^-1] (X-ray primary photoionization rate) ────
  double zeta_X = 0.0;
  const char* env_zeta_X = std::getenv("METAL_ZETA_X");
  if (env_zeta_X && env_zeta_X[0] != '\0') {
    try {
      zeta_X = std::stod(env_zeta_X);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_ZETA_X='%s' is not a valid number\n",
                   env_zeta_X);
      return 1;
    }
    if (zeta_X < 0.0) {
      std::fprintf(stderr, "ERROR: METAL_ZETA_X must be >= 0, got %s\n",
                   env_zeta_X);
      return 1;
    }
  }

  // ── Optional: METAL_E_X_EV [eV] (representative X-ray photon energy) ──────
  double E_X_eV = 300.0;
  const char* env_E_X_eV = std::getenv("METAL_E_X_EV");
  if (env_E_X_eV && env_E_X_eV[0] != '\0') {
    try {
      E_X_eV = std::stod(env_E_X_eV);
    } catch (const std::exception&) {
      std::fprintf(stderr, "ERROR: METAL_E_X_EV='%s' is not a valid number\n",
                   env_E_X_eV);
      return 1;
    }
    if (E_X_eV <= 13.6) {
      std::fprintf(stderr, "ERROR: METAL_E_X_EV must be > 13.6 eV, got %s\n",
                   env_E_X_eV);
      return 1;
    }
  }

  std::string zx_tag;
  if (zeta_X != 0.0) {
    zx_tag = collapse_defaults::make_tag(env_zeta_X, zeta_X);
  }
#endif  // ARCHE_XRAY

  std::printf(
      "=== METAL_ZETA0=%s  METAL_Z_METAL=%s"
      "  cr_tag=%s  z_tag=%s  f_ret=%g  fret_tag=%s  ff_gamma=%d"
      "  METAL_REDSHIFT=%g  T_rad=%g K"
      "  abund=%s"
      "  J_LW21=%g  nH0=%g"
      "  T_K0=%g  y_e0=%g  y_H2=%g  y_HD=%g"
      "  stride=%d  max_iter=%d"
      "  cr_col=%g  cr_sec=%g  cr_bkg=%g"
      "  sra_rate=%g lra_rate=%g"
      "  Cgas=%g Ogas=%g Mggas=%g  T_cr_des=%g"
      "  dt_factor=%g dt_factor_chem=%g dt_factor_init=%g n_init=%d nH_stop=%g"
#ifdef ARCHE_XRAY
      "  ZETA_X=%g  E_X_EV=%g"
#endif
      " ===\n",
      env_zeta0, env_z, cr_tag.c_str(), z_tag.c_str(), f_ret, fret_tag.c_str(),
      static_cast<int>(ff_gamma), zred, T_rad, abund.name, jlw21, nH0, T_K0,
      ye0, y_H2_init, y_HD_init, output_stride, max_iter, cr_atten_col_dens,
      cr_atten_second_frac, cr_metal_bkgnd, sra_rate, lra_rate, c_gas_frac,
      o_gas_frac, mg_gas_frac, t_cr_desorp, dt_factor, dt_factor_chem,
      dt_factor_init, n_init_steps, nH_stop
#ifdef ARCHE_XRAY
      ,
      zeta_X, E_X_eV
#endif
  );
  if (!fret_table_path.empty() && !ff_gamma)
    std::printf("  fret_table: %s  (%zu rows)\n", fret_table_path.c_str(),
                fret_nH.size());
  const char* env_bench = std::getenv("METAL_BENCH");
  const bool bench_mode =
      (env_bench && env_bench[0] != '\0' && env_bench[0] != '0');

  RunCollapse(zeta0, cr_tag, Z_metal, z_tag, f_ret, fret_tag, T_rad, zred,
              zred_tag, abund, *tbl, out_dir, jlw21, jlw_tag, cr_atten_col_dens,
              cr_atten_second_frac, cr_metal_bkgnd, sra_rate, lra_rate,
              t_cr_desorp, nH0, T_K0, ye0, y_H2_init, y_HD_init, c_gas_frac,
              o_gas_frac, mg_gas_frac, fret_nH, fret_val, fret_table_path,
              dt_factor, dt_factor_chem, dt_factor_init, n_init_steps, nH_stop,
              output_stride, max_iter, ff_gamma, bench_mode
#ifdef ARCHE_XRAY
              ,
              zeta_X, E_X_eV, zx_tag
#endif
  );

  return 0;
}
