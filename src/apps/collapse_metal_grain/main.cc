// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
//
// main.cc — one-zone gravitational collapse (metal_grain network)
//
// Port of:
//   fortran/metal_grain/collapse_metal_grain.f
//
// Single case per run; parameters are read from environment variables:
//   METAL_ZETA0          — CR ionization rate [s^-1]       (required)
//   METAL_Z_METAL        — metallicity [Z_sun]              (required)
//   METAL_OUTDIR         — output directory                 (optional, default: results/metal/h5)
//   METAL_FF_RET         — free-fall retardation factor     (optional, default: 1.0)
//   METAL_FRET_TABLE     — 2-column ASCII file: nH[cm^-3]  f_ret  step-function table
//                          Rows sorted by ascending nH; comments with '#' allowed.
//                          If set, METAL_FF_RET is ignored.  fret_tag = "-step"
//                          → output: collapse_CR<cr>_Z<z>_fret-step.h5
//   METAL_FF_GAMMA       — gamma-dependent collapse factor (flag; set to 1 to enable)
//                          Uses t_eff = t_ff / sqrt(1−f(γ))  (Higuchi+2018 Eq.5-7)
//                          Overrides METAL_FF_RET and METAL_FRET_TABLE.
//                          fret_tag = "-gamma"  → output: collapse_CR<cr>_Z<z>_fret-gamma.h5
//   METAL_XNH0           — initial H number density [cm^-3] (optional, default: 1.0)
//   METAL_OUTPUT_STRIDE  — write every N-th step to HDF5    (optional, default: 10)
//   METAL_MAX_ITER       — maximum integration steps        (optional, default: 1000000)
//   METAL_DATA_DIR       — path to metal reaction data      (optional, default: compile-time METAL_DATA_DIR)
//   PRIM_DATA_DIR        — path to shared primordial data   (optional, default: compile-time DATA_DIR)
//   METAL_JLW21          — Lyman-Werner radiation intensity J_21 [10^-21 erg/s/cm^2/Hz/sr]
//                          (optional, default: 0.0 = no LW field)
//                          Activates H2/HD photodissociation and H- photodetachment.
//   METAL_SRA_RATE       — short-lived radionuclide ionization scaling
//                          (optional, default: 0.0; zeta_short = rate * 7.6e-19 [s^-1])
//   METAL_LRA_RATE       — long-lived radionuclide ionization scaling
//                          (optional, default: 0.0; zeta_long  = rate * 1.4e-22 * Z_metal [s^-1])
//   METAL_ABUNDANCE_SET  — abundance preset name (optional, default: solar)
//                          currently supported: solar, default
//
// HDF5 layout (root-level datasets):
//   y          (N_rows, N_sp=89) float64 — species abundances / xnH
//   xnH        (N_rows,)   — H number density [cm^-3]
//   T_K        (N_rows,)   — gas temperature [K]
//   T_gr_K     (N_rows,)   — grain temperature [K]
//   rho        (N_rows,)   — mass density [g/cm^3]
//   xLmbd_net  (N_rows,)   — net cooling [erg g^-1 s^-1]
//   xLmbd_line (N_rows,)   — total line cooling
//   xLmbd_cnt  (N_rows,)   — continuum cooling
//   xLmbd_gr   (N_rows,)   — grain continuum cooling
//   xLmbd_gas  (N_rows,)   — gas continuum cooling (ff + CIA)
//   xLmbd_ch   (N_rows,)   — chemical cooling
//   xGam_cmp   (N_rows,)   — compressional heating
//   xLmbd_H2   (N_rows,)   — H2 line cooling
//   xLmbd_HD   (N_rows,)   — HD line cooling
//   xLmbd_Lya  (N_rows,)   — Lyman-alpha cooling
//   xLmbd_CO   (N_rows,)   — CO line cooling
//   xLmbd_OH   (N_rows,)   — OH line cooling
//   xLmbd_H2O  (N_rows,)   — H2O line cooling
//   xLmbd_CII  (N_rows,)   — CII line cooling
//   xLmbd_CI   (N_rows,)   — CI line cooling
//   xLmbd_OI   (N_rows,)   — OI line cooling
//   xGam_CR    (N_rows,)   — CR ionization heating
//   t_ff       (N_rows,)   — true free-fall time [s]  (= t_eff / f_ret)
//   t_cool     (N_rows,)   — cooling time [s]
//   t_chem     (N_rows,)   — chemistry timescale [s]
//   tau_cnt    (N_rows,)   — continuum optical depth
//   xlmbd_J    (N_rows,)   — Jeans length [cm]
//   xMJ        (N_rows,)   — Jeans mass [g]
//   B_cr       (N_rows,)   — critical magnetic field [G]
//   y_plus     (N_rows,)   — total positive charge fraction
//   y_minus    (N_rows,)   — total negative charge fraction
//   charge_imbal (N_rows,) — |y+−y−|/(y++y−)
//   step       (N_rows,)   int32 — step number (every 10 steps)
//
//   Attributes (root):
//     f_ret        — free-fall retardation factor (initial value; 1.0 = standard free-fall)
//     f_ret_table  — path to f_ret step-function table file (absent if not used)
//     ff_collapse_mode — "gamma" when METAL_FF_GAMMA mode is active (absent otherwise)

#include <cmath>
#include <cstdio>
#include <cstring>
#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <hdf5.h>

#include "hdf5_utils.h"

#include "chemistry.h"
#include "collapse_dynamics.h"

namespace {

// ─── compile-time network sizes ──────────────────────────────────────────────
constexpr int kNSp    = chemistry::metal_grain::N_sp;    // 89

// ─── physical constants (CGS) ─────────────────────────────────────────────────
constexpr double kKB    = chemistry::phys::xk_B;    // 1.380662e-16  [erg/K]
constexpr double kMp    = chemistry::phys::xm_p;    // 1.67262e-24   [g]
constexpr double kPi    = chemistry::phys::pi;
constexpr double kGGrav = chemistry::phys::G;        // 6.6720e-8     [cm^3 g^-1 s^-2]

// ─── abundance baseline / CMB ────────────────────────────────────────────────
constexpr double kTCMB0 = 2.725;                         // CMB temperature at z=0 [K]

// ─── CR attenuation constants ─────────────────────────────────────────────────
// CR attenuation column density [g cm^-2]
constexpr double kCrAttenuColDens   = 96.0;
// Secondary CR fraction (Padovani et al.)
constexpr double kCrAttenuSecondFrac = 7.6e-2;
// Metal background CR floor [erg s^-1 g^-1 / Z_sun]
constexpr double kCrMetalBkgnd      = 1.4e-22;
// Short-lived radionuclide ionization scaling (1_make-compatible, default off)
constexpr double kSraRateDefault    = 0.0;
// Long-lived radionuclide ionization scaling (1_make-compatible, default off)
constexpr double kLraRateDefault    = 0.0;
// Effective CR desorption spike temperature [K]
constexpr double kTcrDesorp         = 70.0;

// ─── MRN grain size distribution parameters ──────────────────────────────────
// Bipartite power-law: n(a) ∝ a^{-3.5} (a_min..a_mid), a^{-5.5} (a_mid..a_max)
constexpr double kGrAMin    = chemistry::metal_grain::mrn_a_min;
constexpr double kGrAMid    = chemistry::metal_grain::mrn_a_mid;
constexpr double kGrAMax    = chemistry::metal_grain::mrn_a_max;
constexpr double kGrInd1    = chemistry::metal_grain::mrn_ind1;
constexpr double kGrInd2    = chemistry::metal_grain::mrn_ind2;
constexpr double kGrNorm    = chemistry::metal_grain::mrn_norm_zsun;

// ─── Initial gas-phase element fractions ─────────────────────────────────────
// (1 - fraction) is initially locked in grain mantles
constexpr double kCgasFrac  = chemistry::metal_grain::c_gas_frac_default;
constexpr double kOgasFrac  = chemistry::metal_grain::o_gas_frac_default;
constexpr double kMggasFrac = chemistry::metal_grain::mg_gas_frac_default;

// ─── Integration control ─────────────────────────────────────────────────────
// Timestep as fraction of min(t_cool, t_eff) after initial phase
constexpr double kDtFactor     = 1.0e-3;
// Timestep factor during the first kNInitSteps steps (short initial kick)
constexpr double kDtFactorInit = 1.0e-8;
constexpr int    kNInitSteps   = 10;
// Default HDF5 output stride (write every N-th step); overridable via METAL_OUTPUT_STRIDE
constexpr int    kOutputStride = 10;
// Maximum number of integration steps; overridable via METAL_MAX_ITER
constexpr int    kItMax        = 1000000;
// Stop when xnH exceeds this threshold [cm^-3]
constexpr double kXnHStop      = 1.0e23;

// ─────────────────────────────────────────────────────────────────────────────
// OutputRow — one record (every 10 steps) buffered before HDF5 write
// ─────────────────────────────────────────────────────────────────────────────
struct OutputRow {
    int    step;
    std::array<double, kNSp> y;
    double xnH, T_K, T_gr_K, rho;
    double xLmbd_net, xLmbd_line, xLmbd_cnt, xLmbd_gr, xLmbd_gas;
    double xLmbd_ch, xGam_cmp;
    double xLmbd_H2, xLmbd_HD, xLmbd_Lya;
    double xLmbd_CO, xLmbd_OH, xLmbd_H2O;
    double xLmbd_CII, xLmbd_CI, xLmbd_OI;
    double xGam_CR;
    double t_ff, t_cool, t_chem, tau_cnt, xlmbd_J, xMJ, B_cr;
    double y_plus, y_minus, charge_imbal;
};

// ─── HDF5 helpers (shared implementation in public/include/hdf5_utils.h) ─────
using h5utils::H5Write1dInt;
using h5utils::H5Write1d;
using h5utils::H5Write2d;
using h5utils::H5WriteStrAttr;
using h5utils::H5WriteDblAttr;
using h5utils::H5Create;

// ─────────────────────────────────────────────────────────────────────────────
// LoadFretTable — read a 2-column ASCII table: nH[cm^-3]  f_ret
//
// Format: rows sorted by ascending nH; lines starting with '#' are comments.
// The returned vectors have equal length >= 1.
// The step-function convention: f_ret[i] applies while nH < nH[i+1].
// f_ret[last] applies for all nH >= nH[last].
// ─────────────────────────────────────────────────────────────────────────────
std::pair<std::vector<double>, std::vector<double>>
LoadFretTable(const std::string& path)
{
    std::ifstream fin(path);
    if (!fin.is_open()) {
        std::fprintf(stderr, "ERROR: cannot open METAL_FRET_TABLE '%s'\n",
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
        std::fprintf(stderr,
            "ERROR: METAL_FRET_TABLE '%s' contains no valid rows\n",
            path.c_str());
        std::exit(1);
    }
    return {nH_tab, fret_tab};
}

// ─────────────────────────────────────────────────────────────────────────────
// WriteHdf5File — write all datasets flat at the root of an open HDF5 file
// ─────────────────────────────────────────────────────────────────────────────
void WriteHdf5File(hid_t fid,
                   const std::string& cr_tag, const std::string& z_tag,
                   const std::vector<OutputRow>& rows,
                   double zeta0, double Z_metal, double f_ret,
                   double sra_rate, double lra_rate,
                   double T_rad, double zred,
                   double T_K0_ic, double y_e0_ic,
                   double y_H2_ic, double y_HD_ic,
                   double jlw21 = 0.0,
                   const std::string& fret_table_path = "",
                   bool ff_gamma = false)
{
    const hsize_t N   = static_cast<hsize_t>(rows.size());
    const hsize_t Nsp = static_cast<hsize_t>(kNSp);

    // ── Root attributes ───────────────────────────────────────────────────────
    {
        char desc[160];
        std::snprintf(desc, sizeof(desc),
            "metal_grain collapse CR%s Z%s (C++ port)",
            cr_tag.c_str(), z_tag.c_str());
        H5WriteStrAttr(fid, "description", std::string(desc));
    }
    H5WriteStrAttr(fid, "cr_tag",        cr_tag);
    H5WriteStrAttr(fid, "z_tag",         z_tag);
    H5WriteDblAttr(fid, "zeta0_cgs",     zeta0);
    H5WriteDblAttr(fid, "Z_metal",       Z_metal);
    H5WriteDblAttr(fid, "f_ret",         f_ret);
    H5WriteDblAttr(fid, "sra_rate",      sra_rate);
    H5WriteDblAttr(fid, "lra_rate",      lra_rate);
    if (!fret_table_path.empty())
        H5WriteStrAttr(fid, "f_ret_table", fret_table_path);
    if (ff_gamma)
        H5WriteStrAttr(fid, "ff_collapse_mode", "gamma");
    H5WriteDblAttr(fid, "zred",          zred);
    H5WriteDblAttr(fid, "T_rad",         T_rad);
    H5WriteDblAttr(fid, "J_LW21",        jlw21);
    H5WriteDblAttr(fid, "ic_T_K0",       T_K0_ic);
    H5WriteDblAttr(fid, "ic_y_e0",       y_e0_ic);
    H5WriteDblAttr(fid, "ic_y_H2",       y_H2_ic);
    H5WriteDblAttr(fid, "ic_y_HD",       y_HD_ic);
    H5WriteStrAttr(fid, "network",       "metal_grain N_sp=89 N_react=1200");
    H5WriteStrAttr(fid, "units_density", "cm^-3 (xnH)  g/cm^3 (rho)");
    H5WriteStrAttr(fid, "units_cooling", "erg g^-1 s^-1");
    H5WriteStrAttr(fid, "units_time",    "s");
    H5WriteStrAttr(fid, "units_length",  "cm");
    H5WriteStrAttr(fid, "units_B",       "G");

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
                mat[i*Nsp + j] = rows[i].y[static_cast<int>(j)];
        H5Write2d(fid, "y", mat, N, Nsp);

        hid_t ds = H5Dopen2(fid, "y", H5P_DEFAULT);
        H5WriteStrAttr(ds, "species",
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
            "C(gr),CH(gr),CH2(gr),CH3(gr),CH4(gr)");
        H5WriteStrAttr(ds, "units", "number fraction relative to xnH");
        H5Dclose(ds);
    }

    // ── scalar 1-D datasets ───────────────────────────────────────────────────
    auto dump1d = [&](std::string_view name, auto fn) {
        std::vector<double> v(N);
        for (hsize_t i = 0; i < N; ++i) v[i] = fn(rows[i]);
        H5Write1d(fid, name, v);
    };

    dump1d("xnH",          [](const OutputRow& r){ return r.xnH;         });
    dump1d("T_K",          [](const OutputRow& r){ return r.T_K;         });
    dump1d("T_gr_K",       [](const OutputRow& r){ return r.T_gr_K;      });
    dump1d("rho",          [](const OutputRow& r){ return r.rho;         });
    dump1d("xLmbd_net",    [](const OutputRow& r){ return r.xLmbd_net;   });
    dump1d("xLmbd_line",   [](const OutputRow& r){ return r.xLmbd_line;  });
    dump1d("xLmbd_cnt",    [](const OutputRow& r){ return r.xLmbd_cnt;   });
    dump1d("xLmbd_gr",     [](const OutputRow& r){ return r.xLmbd_gr;    });
    dump1d("xLmbd_gas",    [](const OutputRow& r){ return r.xLmbd_gas;   });
    dump1d("xLmbd_ch",     [](const OutputRow& r){ return r.xLmbd_ch;    });
    dump1d("xGam_cmp",     [](const OutputRow& r){ return r.xGam_cmp;    });
    dump1d("xLmbd_H2",     [](const OutputRow& r){ return r.xLmbd_H2;    });
    dump1d("xLmbd_HD",     [](const OutputRow& r){ return r.xLmbd_HD;    });
    dump1d("xLmbd_Lya",    [](const OutputRow& r){ return r.xLmbd_Lya;   });
    dump1d("xLmbd_CO",     [](const OutputRow& r){ return r.xLmbd_CO;    });
    dump1d("xLmbd_OH",     [](const OutputRow& r){ return r.xLmbd_OH;    });
    dump1d("xLmbd_H2O",    [](const OutputRow& r){ return r.xLmbd_H2O;   });
    dump1d("xLmbd_CII",    [](const OutputRow& r){ return r.xLmbd_CII;   });
    dump1d("xLmbd_CI",     [](const OutputRow& r){ return r.xLmbd_CI;    });
    dump1d("xLmbd_OI",     [](const OutputRow& r){ return r.xLmbd_OI;    });
    dump1d("xGam_CR",      [](const OutputRow& r){ return r.xGam_CR;     });
    dump1d("t_ff",         [](const OutputRow& r){ return r.t_ff;        });
    dump1d("t_cool",       [](const OutputRow& r){ return r.t_cool;      });
    dump1d("t_chem",       [](const OutputRow& r){ return r.t_chem;      });
    dump1d("tau_cnt",      [](const OutputRow& r){ return r.tau_cnt;     });
    dump1d("xlmbd_J",      [](const OutputRow& r){ return r.xlmbd_J;     });
    dump1d("xMJ",          [](const OutputRow& r){ return r.xMJ;         });
    dump1d("B_cr",         [](const OutputRow& r){ return r.B_cr;        });
    dump1d("y_plus",       [](const OutputRow& r){ return r.y_plus;      });
    dump1d("y_minus",      [](const OutputRow& r){ return r.y_minus;     });
    dump1d("charge_imbal", [](const OutputRow& r){ return r.charge_imbal;});
}

// ─────────────────────────────────────────────────────────────────────────────
// RunCollapse — port of Fortran subroutine collapse()
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
// shielding lengths (xlsh), and the integration timestep (dt).
// ─────────────────────────────────────────────────────────────────────────────
void RunCollapse(double zeta0, const std::string& cr_tag,
                 double Z_metal, const std::string& z_tag,
                 double f_ret_init, const std::string& fret_tag,
                 double T_rad, double zred, const std::string& zred_tag,
                 const chemistry::abundance::MetalSet& abund,
                 const chemistry::MetalGrainTable& tbl,
                 const std::string& out_dir,
                 double jlw21 = 0.0,
                 const std::string& jlw_tag = "",
                 double cr_atten_col_dens = kCrAttenuColDens,
                 double cr_atten_second_frac = kCrAttenuSecondFrac,
                 double cr_metal_bkgnd = kCrMetalBkgnd,
                 double sra_rate = kSraRateDefault,
                 double lra_rate = kLraRateDefault,
                 double t_cr_desorp = kTcrDesorp,
                 double xnH0       = 1.0,
                 double T_K0       = 1.0e2,
                 double y_e0       = 1.0e-4,
                 double y_H2_init  = 6.0e-7,
                 double y_HD_init  = 4.0e-10,
                 double c_gas_frac = kCgasFrac,
                 double o_gas_frac = kOgasFrac,
                 double mg_gas_frac = kMggasFrac,
                 const std::vector<double>& fret_nH  = {},
                 const std::vector<double>& fret_val = {},
                 const std::string& fret_table_path  = "",
                 double dt_factor = kDtFactor,
                 double dt_factor_init = kDtFactorInit,
                 int n_init_steps = kNInitSteps,
                 double xnH_stop = kXnHStop,
                 int output_stride = kOutputStride,
                 int max_iter      = kItMax,
                 bool ff_gamma     = false)
{
    // ── Metallicity-scaled abundances ─────────────────────────────────────────
    const double XC  = abund.yC  * Z_metal;
    const double XO  = abund.yO  * Z_metal;
    const double XKa = abund.yK  * Z_metal;
    const double XNa = abund.yNa * Z_metal;
    const double XMg = abund.yMg * Z_metal;

    // ── Initial conditions ────────────────────────────────────────────────────
    const double     xnH1      = xnH0;    // initial H number density [cm^-3]
    const double     kT1       = T_K0;    // initial temperature [K]

    chemistry::MetalGrainCell cell;
    auto& y = cell.state.y;

    const double y_H2 = y_H2_init;
    const double y_Dp = 0.0;
    const double y_HD = y_HD_init;
    const double y_Hp = y_e0;
    const double y_H  = 1.0 - y_Hp - 2.0*y_H2 - y_HD;
    if (y_H <= 0.0) {
        std::fprintf(stderr,
            "ERROR: Invalid IC: y_H = 1 - y_e0(%.3g) - 2*y_H2(%.3g) - y_HD(%.3g) = %.3g <= 0\n",
            y_Hp, y_H2, y_HD, y_H);
        std::exit(1);
    }
    const double y_D  = abund.yD - y_Dp - y_HD;

    y[0]  = y_H;          // H
    y[1]  = y_H2;         // H2
    y[3]  = y_Hp;         // H+
    y[7]  = abund.yHe;    // He
    y[11] = y_D;          // D
    y[12] = y_HD;         // HD
    y[13] = y_Dp;         // D+
    y[52] = abund.yLi;    // Li+ (Fortran y(53) = Li+)
    y[16] = XC  * c_gas_frac;    // CI
    y[29] = XO  * o_gas_frac;    // OI
    y[58] = 0.0;                // K+  (initially zero)
    y[60] = 0.0;                // Na+ (initially zero)
    y[62] = XMg * mg_gas_frac;   // Mg+
    // e- = sum of positive ions
    y[2]  = y_Hp + y_Dp + abund.yLi + y[58] + y[60] + y[62];

    // MRN grain abundance: yGr = Z_metal * kGrNorm * (1+4yHe)*mp / <vol>
    // Bipartite power-law: a^{-kGrInd1} (kGrAMin..kGrAMid), a^{-kGrInd2} (kGrAMid..kGrAMax)
    double yGr = 0.0;
    if (Z_metal > 0.0) {
        double sum_ana
            = std::pow(kGrAMid, kGrInd1)
              * (std::pow(kGrAMid, 1.0-kGrInd1) - std::pow(kGrAMin, 1.0-kGrInd1)) / (1.0-kGrInd1)
            + std::pow(kGrAMid, kGrInd2)
              * (std::pow(kGrAMax, 1.0-kGrInd2) - std::pow(kGrAMid, 1.0-kGrInd2)) / (1.0-kGrInd2);
        double ave_gr_vol
            = std::pow(kGrAMid, kGrInd1)
              * (std::pow(kGrAMid, 4.0-kGrInd1) - std::pow(kGrAMin, 4.0-kGrInd1)) / (4.0-kGrInd1)
            + std::pow(kGrAMid, kGrInd2)
              * (std::pow(kGrAMax, 4.0-kGrInd2) - std::pow(kGrAMid, 4.0-kGrInd2)) / (4.0-kGrInd2);
        ave_gr_vol = (ave_gr_vol / sum_ana) * (4.0 * kPi / 3.0);
        yGr = Z_metal * kGrNorm * (1.0 + 4.0*abund.yHe) * kMp / ave_gr_vol;
    }
    y[63] = yGr;   // Gr (neutral grain)

    // Grain surface depletion reservoirs (elements locked in grains initially)
    double yCs   = XC  * (1.0 - c_gas_frac);
    double yOs   = XO  * (1.0 - o_gas_frac);
    double yKs   = XKa;
    double yNas  = XNa;
    double yMgs  = XMg * (1.0 - mg_gas_frac);
    int switch_gr = 1;   // 1 = grains present, 2 = fully evaporated

    // ── Thermodynamic state ───────────────────────────────────────────────────
    double rho = ((1.0 + 4.0*abund.yHe) * kMp) * xnH1;
    double xnH = xnH1;
    double T_K = kT1;
    double p   = (1.0 + abund.yHe) * xnH1 * kKB * kT1;

    double xmu   = (1.0 + 4.0*abund.yHe)
                 / (y[0] + y[1] + y[2] + y[3] + y[7] + y[8] + y[9]);
    double gamma = 1.0 + (1.0 + 4.0*abund.yHe)
                 / (xmu * (1.5*(y[0] + y[2] + y[3] + y[7] + y[8] + y[9])
                          + chemistry::c_H2(T_K) * y[1]));
    double e     = kKB * T_K / ((gamma - 1.0) * xmu * kMp);

    double T_gr_K   = T_rad;
    double rho_old  = rho;
    double T_gr_old = T_gr_K;

    // ── Chemistry parameters ──────────────────────────────────────────────────
    chemistry::ChemParams params{};
    params.T_rad   = T_rad;
    params.zeta    = zeta0;
    params.Z_metal = Z_metal;
    params.T_gr_K  = T_gr_K;
    params.T_cr_desorp = t_cr_desorp;

    // ── Scalar loop state ─────────────────────────────────────────────────────
    double t        = 0.0;
    double dt       = 1.0e-1;
    double t_chem   = 1.0e-1;
    double t_cool   = 0.0;
    double xk_gr    = 0.0;
    double xk_gas   = 0.0;
    double tau_cnt  = 0.0;
    double esc_cnt  = 1.0;

    // Column densities
    double xNc_H    = 0.0;
    double xNc_H2   = 0.0;
    double xNc_HD   = 0.0;
    double xNc_CII  = 0.0;
    double xNc_CI   = 0.0;
    double xNc_OI   = 0.0;
    double xNc_CO   = 0.0;
    double xNc_H2O  = 0.0;
    double xNc_OH   = 0.0;

    // ── f_ret step-function table state ──────────────────────────────────────
    double f_ret = f_ret_init;
    int    fret_idx         = 0;
    const bool has_fret_tab = !fret_nH.empty();

    // ── HDF5 row buffer ───────────────────────────────────────────────────────
    std::vector<OutputRow> h5rows;
    h5rows.reserve(6000);

    std::printf("CR%s Z%s: zeta0=%g z=%g yGr=%g\n",
                cr_tag.c_str(), z_tag.c_str(), zeta0/1.0e-17, Z_metal, yGr);

    // ── Time integration ──────────────────────────────────────────────────────
    for (int it = 1; it <= max_iter; ++it) {

        // Save state from previous step for grain evaporation ratio
        rho_old  = rho;
        T_gr_old = T_gr_K;

        // ── f_ret ratchet: advance step-function when xnH crosses next threshold ─
        if (has_fret_tab) {
            while (fret_idx + 1 < static_cast<int>(fret_nH.size()) &&
                   xnH >= fret_nH[fret_idx + 1]) {
                ++fret_idx;
                f_ret = fret_val[fret_idx];
            }
        }

        // ── Free-fall time & Jeans length ─────────────────────────────────────
        double t_ff    = std::sqrt(3.0*kPi / (32.0*kGGrav*rho));
        double t_eff   = t_eff_collapse(t_ff, f_ret, gamma, ff_gamma);
        double xlmbd_J = std::sqrt(kPi * kKB * T_K / (kGGrav * xmu * kMp * rho));

        // ── Thermal velocities for column density ─────────────────────────────
        auto vD = [&](double mass_amu) {
            return std::sqrt(2.0 * kKB * T_K / (kMp * mass_amu));
        };
        double xlsh_H   = std::min(xlmbd_J, 6.0 * vD(1.0)  * t_eff);
        double xlsh_H2  = std::min(xlmbd_J, 6.0 * vD(2.0)  * t_eff);
        double xlsh_HD  = std::min(xlmbd_J, 6.0 * vD(3.0)  * t_eff);
        double xlsh_C   = std::min(xlmbd_J, 6.0 * vD(12.0) * t_eff);
        double xlsh_O   = std::min(xlmbd_J, 6.0 * vD(16.0) * t_eff);
        double xlsh_CO  = std::min(xlmbd_J, 6.0 * vD(28.0) * t_eff);
        double xlsh_H2O = std::min(xlmbd_J, 6.0 * vD(18.0) * t_eff);
        double xlsh_OH  = std::min(xlmbd_J, 6.0 * vD(17.0) * t_eff);

        xNc_H   = y[0]  * xnH * xlsh_H;
        xNc_H2  = y[1]  * xnH * xlsh_H2;
        xNc_HD  = y[12] * xnH * xlsh_HD;
        xNc_CII = y[22] * xnH * xlsh_C;   // C+
        xNc_CI  = y[16] * xnH * xlsh_C;   // C
        xNc_OI  = y[29] * xnH * xlsh_O;   // O
        xNc_CO  = y[32] * xnH * xlsh_CO;
        xNc_H2O = y[33] * xnH * xlsh_H2O;
        xNc_OH  = y[31] * xnH * xlsh_OH;

        // ── Continuum optical depth ───────────────────────────────────────────
        tau_cnt = (xk_gr + xk_gas) * rho * xlmbd_J;
        if (tau_cnt > 1.0)
            esc_cnt = 1.0 / (tau_cnt * tau_cnt);
        else
            esc_cnt = 1.0;

        // ── Build shielding environment ───────────────────────────────────────
        chemistry::ChemShielding shield;
        const double zeta_r_short = sra_rate * 7.6e-19;
        const double zeta_r_long  = lra_rate * 1.4e-22 * Z_metal;
        if (zeta0 > 0.0)
            shield.zeta = zeta0 * (std::exp(-rho * xlmbd_J / cr_atten_col_dens)
                                   + cr_atten_second_frac)
                        + cr_metal_bkgnd * Z_metal
                        + zeta_r_short + zeta_r_long;
        else
            shield.zeta = zeta_r_short + zeta_r_long;
        shield.xNc_H   = xNc_H;
        shield.xNc_H2  = xNc_H2;
        shield.xNc_HD  = xNc_HD;
        shield.xNc_CO  = xNc_CO;
        shield.xNc_OH  = xNc_OH;
        shield.xNc_H2O = xNc_H2O;
        shield.xNc_CII = xNc_CII;
        shield.xNc_CI  = xNc_CI;
        shield.xNc_OI  = xNc_OI;
        shield.tau_cnt = tau_cnt;
        shield.esc_cnt = esc_cnt;
        shield.J_LW21  = jlw21;

        cell.state.xnH   = xnH;
        cell.state.T_K   = T_K;
        cell.state.xmu   = xmu;
        cell.state.gamma = gamma;

        const auto y_prev = y;
        const auto rates  = chemistry::chem_full_step(cell, dt, params, shield, tbl);

        xmu    = cell.state.xmu;
        gamma  = cell.state.gamma;
        T_gr_K = params.T_gr_K;
        xk_gr  = rates.xk_gr;
        xk_gas = rates.xk_gas;

        const double xLmbd_line    = rates.xLmbd_line;
        const double xLmbd_gr      = rates.xLmbd_gr;
        const double xLmbd_gas_out = rates.xLmbd_gas;
        const double xLmbd_cnt     = rates.xLmbd_cnt;
        const double xLmbd_ch      = rates.xLmbd_ch;
        const double xGam_CR       = rates.xGam_CR;
        const double xLmbd_net     = rates.xLmbd_net;

        // ── Chemistry timescale: min_i( y_i / |Δy_i / Δt| ) ─────────────────
        t_chem = 1.0e50;
        for (int i = 0; i < kNSp; ++i) {
            const double dy = std::abs(y[i] - y_prev[i]);
            if (y[i] > 1.0e-30 && dy > 1.0e-40)
                t_chem = std::min(t_chem, dt * y[i] / dy);
        }
        t_chem = std::max(t_chem, dt);
        double xGam_cmp = p / rho / t_eff;

        double xMJ = (4.0*kPi/3.0) * rho * xlmbd_J*xlmbd_J*xlmbd_J;
        double B_cr = std::sqrt(4.0*kPi * kGGrav * xMJ * rho / xlmbd_J);

        // ── Grain evaporation ─────────────────────────────────────────────────
        if (switch_gr == 1) {
            double T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol;
            chemistry::detail::vaptemp(rho, T_ice, T_vo, T_ro, T_tr, T_ir, T_pyr, T_ol);

            double vol_cur = chemistry::detail::vol_gr(rho,     T_gr_K);
            double vol_old = chemistry::detail::vol_gr(rho_old, T_gr_old);

            if (vol_cur < 1.0e-60) {
                // Complete evaporation: release all grain charge and alkali elements
                double del_yp = y[64] + 2.0*y[65];   // Gr+ + 2*Gr2+
                double del_ym = y[66] + 2.0*y[67];   // Gr- + 2*Gr2-
                if (del_ym >= del_yp) {
                    y[2]  += (del_ym - del_yp);       // e-
                } else {
                    y[58] += (del_yp - del_ym);       // K+
                    y[57] -= (del_yp - del_ym);       // K (neutral)
                }
                y[57] += yKs;    // K neutral released
                y[59] += yNas;   // Na neutral released
                y[61] += yMgs;   // Mg neutral released
                yKs  = 0.0; yNas = 0.0; yMgs = 0.0;
                yCs  = 0.0; yOs  = 0.0;
                y[63] = 0.0; y[64] = 0.0; y[65] = 0.0;
                y[66] = 0.0; y[67] = 0.0;
                switch_gr = 2;

            } else {
                double gr_frac = (vol_old > 1.0e-60) ? std::min(vol_cur / vol_old, 1.0) : 1.0;

                if (T_K >= 0.975 * T_pyr) {
                    // Partial evaporation of grain charges and alkali elements
                    double del_yp = (y[64] + 2.0*y[65]) * (1.0 - gr_frac);
                    double del_ym = (y[66] + 2.0*y[67]) * (1.0 - gr_frac);
                    if (del_ym >= del_yp) {
                        y[2]  += (del_ym - del_yp);
                    } else {
                        y[58] += (del_yp - del_ym);
                        y[57] -= (del_yp - del_ym);
                    }
                    y[57] += yKs  * (1.0 - gr_frac);
                    y[59] += yNas * (1.0 - gr_frac);
                    y[61] += yMgs * (1.0 - gr_frac);
                    yCs  *= gr_frac;  yOs  *= gr_frac;
                    yKs  *= gr_frac;  yNas *= gr_frac;  yMgs *= gr_frac;
                    y[64] *= gr_frac; y[65] *= gr_frac;
                    y[66] *= gr_frac; y[67] *= gr_frac;
                }
                y[63] *= gr_frac;   // Gr neutral scales with grain volume
            }
        }

        // ── Charge accounting ─────────────────────────────────────────────────
        double y_plus
            = y[3] + y[4] + y[5] + y[8] + y[10] + y[13] + y[14]
            + y[22] + y[23] + y[24] + y[25] + y[26] + y[27] + y[28]
            + y[39] + y[40] + y[41] + y[42] + y[43] + y[44] + y[45]
            + y[46] + y[47] + y[48] + y[49]
            + y[52] + y[54] + y[58] + y[60] + y[62] + y[64]
            + 2.0*(y[9] + y[55] + y[65]) + 3.0*y[56];
        double y_minus
            = y[2] + y[6] + y[15] + y[53] + y[66] + 2.0*y[67];
        double y_sum = y_plus + y_minus;
        double charge_imbal = (y_sum > 0.0)
            ? std::abs(y_plus - y_minus) / y_sum : 0.0;

        // ── Output (every output_stride steps) ───────────────────────────────
        if (it % output_stride == 0) {
            OutputRow row{};
            row.step        = it;
            row.y           = y;
            row.xnH         = xnH;
            row.T_K         = T_K;
            row.T_gr_K      = T_gr_K;
            row.rho         = rho;
            row.xLmbd_net   = xLmbd_net;
            row.xLmbd_line  = xLmbd_line;
            row.xLmbd_cnt   = xLmbd_cnt;
            row.xLmbd_gr    = xLmbd_gr;
            row.xLmbd_gas   = xLmbd_gas_out;
            row.xLmbd_ch    = xLmbd_ch;
            row.xGam_cmp    = xGam_cmp;
            row.xLmbd_H2    = rates.xLmbd_H2;
            row.xLmbd_HD    = rates.xLmbd_HD;
            row.xLmbd_Lya   = rates.xLmbd_Lya;
            row.xLmbd_CO    = rates.xLmbd_CO;
            row.xLmbd_OH    = rates.xLmbd_OH;
            row.xLmbd_H2O   = rates.xLmbd_H2O;
            row.xLmbd_CII   = rates.xLmbd_CII;
            row.xLmbd_CI    = rates.xLmbd_CI;
            row.xLmbd_OI    = rates.xLmbd_OI;
            row.xGam_CR     = xGam_CR;
            row.t_ff        = t_ff;
            row.t_cool      = t_cool;
            row.t_chem      = t_chem;
            row.tau_cnt     = tau_cnt;
            row.xlmbd_J     = xlmbd_J;
            row.xMJ         = xMJ;
            row.B_cr        = B_cr;
            row.y_plus      = y_plus;
            row.y_minus     = y_minus;
            row.charge_imbal = charge_imbal;
            h5rows.push_back(row);
        }

        // ── Update thermodynamic state ────────────────────────────────────────
        double drho = dt * rho / t_eff;
        double de   = -xLmbd_net * dt + drho * p / (rho * rho);
        rho  += drho;
        e    += de;
        if (e <= 0.0) break;  // catch before T_K goes negative
        T_K   = e * ((gamma - 1.0) * xmu * kMp) / kKB;
        T_K   = std::max(T_K, 1.0);
        p     = rho * kKB * T_K / (xmu * kMp);
        xnH   = rho / ((1.0 + 4.0*abund.yHe) * kMp);
        t    += dt;

        // ── Timestep ──────────────────────────────────────────────────────────
        t_cool = (xLmbd_net != 0.0) ? e / std::abs(xLmbd_net) : 1.0e50;
        if (it <= n_init_steps)
            dt = dt_factor_init * t_eff;
        else
            //dt = dt_factor * std::min({t_cool, t_eff, t_chem});
            dt = dt_factor * std::min({t_cool, t_eff});

        // ── Progress ──────────────────────────────────────────────────────────
        std::printf("%7d %11.3E %11.3E %11.3E %11.3E %11.3E %11.3E"
                    " %11.3E %11.3E %11.3E\n",
                    it, zeta0/1.0e-17, Z_metal, xnH, T_K, y[2], y[1],
                    y_plus, y_minus, charge_imbal);

        if (xnH > xnH_stop || !std::isfinite(xnH) || !std::isfinite(T_K) || e <= 0.0) break;
    }

    // ── Write HDF5 file ───────────────────────────────────────────────────────
    std::string h5_path = out_dir + "/collapse_CR" + cr_tag + "_Z" + z_tag;
    if (!fret_tag.empty()) h5_path += "_fret" + fret_tag;
    if (!jlw_tag.empty())  h5_path += "_JLW"  + jlw_tag;
    if (!zred_tag.empty()) h5_path += "_z"    + zred_tag;
    h5_path += ".h5";
    hid_t fid = H5Create(h5_path);
    if (fid < 0) {
        std::fprintf(stderr, "ERROR: cannot create %s\n", h5_path.c_str());
        return;
    }
    WriteHdf5File(fid, cr_tag, z_tag, h5rows, zeta0, Z_metal, f_ret_init,
                  sra_rate, lra_rate,
                  T_rad, zred, T_K0, y_e0, y_H2_init, y_HD_init,
                  jlw21, fret_table_path, ff_gamma);
    H5Fclose(fid);

    std::printf("  -> %s  (%zu rows)\n", h5_path.c_str(), h5rows.size());
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
int main()
{
    // ── Data directories (compile-time defaults; overridable at runtime) ──────
    const char* env_prim_data  = std::getenv("PRIM_DATA_DIR");
    const char* env_metal_data = std::getenv("METAL_DATA_DIR");
    std::string data_dir       = (env_prim_data  && env_prim_data[0]  != '\0')
                                 ? env_prim_data  : DATA_DIR;
    std::string metal_data_dir = (env_metal_data && env_metal_data[0] != '\0')
                                 ? env_metal_data : METAL_DATA_DIR;

    // ── Output directory ──────────────────────────────────────────────────────
    const char* env_out = std::getenv("METAL_OUTDIR");
    std::string out_dir = env_out ? env_out : "results/metal/h5";
    std::filesystem::create_directories(out_dir);
    std::printf("Output directory: %s/\n", out_dir.c_str());

    // ── Load reaction tables ──────────────────────────────────────────────────
    chemistry::MetalGrainTable tbl;
    chemistry::load_reaction_table     (tbl, metal_data_dir + "/react_metal_grain.dat");
    chemistry::load_grain_surface_table(tbl, metal_data_dir + "/react_grain_surface.dat");
    chemistry::load_mass_table         (tbl, data_dir       + "/mass_metal.dat");
    chemistry::load_saha_table         (tbl, metal_data_dir + "/react_metal_saha.dat");

    // Indices are 1-based Fortran species numbers (C++ index + 1).
    // Shared files from DATA_DIR; metal-specific from METAL_DATA_DIR.
    chemistry::load_partition_functions(tbl, {
        {  6,  data_dir       + "/pf_H3p.dat"       },  // H3+   C++ y[5]
        { 13,  data_dir       + "/pf_HD.dat"         },  // HD    C++ y[12]
        { 15,  data_dir       + "/pf_HDp.dat"        },  // HD+   C++ y[14]
        { 21,  metal_data_dir + "/pf_CH3.dat"        },  // CH3   C++ y[20]
        { 22,  metal_data_dir + "/pf_CH4.dat"        },  // CH4   C++ y[21]
        { 36,  metal_data_dir + "/pf_HO2.dat"        },  // HO2   C++ y[35]
        { 37,  metal_data_dir + "/pf_CO2.dat"        },  // CO2   C++ y[36]
        { 38,  metal_data_dir + "/pf_H2CO.dat"       },  // H2CO  C++ y[37]
        { 39,  metal_data_dir + "/pf_H2O2.dat"       },  // H2O2  C++ y[38]
        { 55,  metal_data_dir + "/pf_LiHp.dat"       },  // LiH+  C++ y[54]
    });

    // ── Required: METAL_ZETA0 [s^-1] ─────────────────────────────────────────
    const char* env_zeta0 = std::getenv("METAL_ZETA0");
    if (!env_zeta0 || env_zeta0[0] == '\0') {
        std::fprintf(stderr,
            "ERROR: environment variable METAL_ZETA0 is required\n"
            "       Set the CR ionization rate [s^-1], e.g.: METAL_ZETA0=1e-17\n");
        return 1;
    }
    double zeta0 = 0.0;
    try {
        zeta0 = std::stod(env_zeta0);
    } catch (const std::exception&) {
        std::fprintf(stderr,
            "ERROR: METAL_ZETA0='%s' is not a valid number\n", env_zeta0);
        return 1;
    }
    if (zeta0 < 0.0) {
        std::fprintf(stderr,
            "ERROR: METAL_ZETA0 must be >= 0, got %s\n", env_zeta0);
        return 1;
    }

    // ── Required: METAL_Z_METAL [Z_sun] ──────────────────────────────────────
    const char* env_z = std::getenv("METAL_Z_METAL");
    if (!env_z || env_z[0] == '\0') {
        std::fprintf(stderr,
            "ERROR: environment variable METAL_Z_METAL is required\n"
            "       Set the metallicity [Z_sun], e.g.: METAL_Z_METAL=1e-3\n");
        return 1;
    }
    double Z_metal = 0.0;
    try {
        Z_metal = std::stod(env_z);
    } catch (const std::exception&) {
        std::fprintf(stderr,
            "ERROR: METAL_Z_METAL='%s' is not a valid number\n", env_z);
        return 1;
    }
    if (Z_metal < 0.0) {
        std::fprintf(stderr,
            "ERROR: METAL_Z_METAL must be >= 0, got %s\n", env_z);
        return 1;
    }

    // ── Tags: input string with '.' -> 'p'; any zero value -> "0" ────────────
    auto make_tag = [](const char* s, double val) -> std::string {
        if (val == 0.0) return "0";
        std::string t(s);
        for (char& c : t) if (c == '.') c = 'p';
        return t;
    };
    std::string cr_tag = make_tag(env_zeta0, zeta0);
    std::string z_tag  = make_tag(env_z,     Z_metal);

    // ── Optional: METAL_REDSHIFT (cosmological redshift, default 0.0) ─────────
    // T_rad = kTCMB0 * (1 + z).  Affects line/continuum/molecular cooling,
    // grain temperature equilibrium, and initial grain temperature.
    double zred = 0.0;
    std::string zred_tag;   // empty → no suffix in filename
    const char* env_zred = std::getenv("METAL_REDSHIFT");
    if (env_zred && env_zred[0] != '\0') {
        try {
            zred = std::stod(env_zred);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: METAL_REDSHIFT='%s' is not a valid number\n", env_zred);
            return 1;
        }
        if (zred < 0.0) {
            std::fprintf(stderr,
                "ERROR: METAL_REDSHIFT must be >= 0, got %s\n", env_zred);
            return 1;
        }
        if (zred != 0.0) {
            zred_tag = env_zred;
            for (char& c : zred_tag) if (c == '.') c = 'p';
        }
    }
    const double T_rad = kTCMB0 * (1.0 + zred);

    // ── Optional: METAL_FRET_TABLE (step-function f_ret table; overrides METAL_FF_RET) ─
    // ASCII 2-col file: nH[cm^-3]  f_ret  (sorted by ascending nH, '#' comments allowed)
    const char* env_fret_table = std::getenv("METAL_FRET_TABLE");
    std::string fret_table_path;
    std::vector<double> fret_nH, fret_val;
    if (env_fret_table && env_fret_table[0] != '\0') {
        fret_table_path = env_fret_table;
        auto [nh, fv] = LoadFretTable(fret_table_path);
        fret_nH  = std::move(nh);
        fret_val = std::move(fv);
    }

    // ── Optional: METAL_FF_GAMMA (gamma-dependent collapse factor; highest priority) ──
    // When set (non-empty, non-"0"), uses t_eff = t_ff / sqrt(1-f(γ))  [Higuchi+2018 Eq.5-7].
    // Overrides METAL_FF_RET and METAL_FRET_TABLE.
    const char* env_ff_gamma = std::getenv("METAL_FF_GAMMA");
    const bool  ff_gamma = (env_ff_gamma && env_ff_gamma[0] != '\0'
                             && env_ff_gamma[0] != '0');

    // ── Optional: METAL_FF_RET (free-fall retardation factor, default 1.0) ───
    // f_ret > 1 slows the collapse; f_ret = 1 gives standard free-fall.
    // Ignored when METAL_FRET_TABLE or METAL_FF_GAMMA is set.
    const char* env_fret = std::getenv("METAL_FF_RET");
    double f_ret = 1.0;
    std::string fret_tag;   // empty → no suffix in filename
    if (ff_gamma) {
        // Gamma mode: tag "-gamma"; warn if other fret options were also supplied
        fret_tag = "-gamma";
        if (!fret_table_path.empty())
            std::fprintf(stderr,
                "WARNING: METAL_FF_GAMMA is set; METAL_FRET_TABLE='%s' is ignored\n",
                env_fret_table);
        if (env_fret && env_fret[0] != '\0')
            std::fprintf(stderr,
                "WARNING: METAL_FF_GAMMA is set; METAL_FF_RET='%s' is ignored\n",
                env_fret);
    } else if (!fret_table_path.empty()) {
        // Table mode: initial f_ret from first table entry; fixed tag "-step"
        // → output filename: collapse_CR<cr>_Z<z>_fret-step.h5
        f_ret    = fret_val[0];
        fret_tag = "-step";
        if (env_fret && env_fret[0] != '\0')
            std::fprintf(stderr,
                "WARNING: METAL_FRET_TABLE is set; METAL_FF_RET='%s' is ignored\n",
                env_fret);
    } else if (env_fret && env_fret[0] != '\0') {
        try {
            f_ret = std::stod(env_fret);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: METAL_FF_RET='%s' is not a valid number\n", env_fret);
            return 1;
        }
        if (f_ret <= 0.0) {
            std::fprintf(stderr,
                "ERROR: METAL_FF_RET must be > 0, got %s\n", env_fret);
            return 1;
        }
        if (f_ret != 1.0) {
            fret_tag = env_fret;
            for (char& c : fret_tag) if (c == '.') c = 'p';
        }
    }

    // ── Optional: METAL_JLW21 (Lyman-Werner intensity J_21, default 0.0) ─────
    // J_21 units: 10^-21 erg/s/cm^2/Hz/sr.  0.0 = no LW field (default).
    double jlw21 = 0.0;
    const char* env_jlw21 = std::getenv("METAL_JLW21");
    if (env_jlw21 && env_jlw21[0] != '\0') {
        try {
            jlw21 = std::stod(env_jlw21);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: METAL_JLW21='%s' is not a valid number\n", env_jlw21);
            return 1;
        }
        if (jlw21 < 0.0) {
            std::fprintf(stderr,
                "ERROR: METAL_JLW21 must be >= 0, got %s\n", env_jlw21);
            return 1;
        }
    }
    // Tag: "." → "p"; zero → empty (no suffix)
    std::string jlw_tag;
    if (jlw21 != 0.0) {
        jlw_tag = env_jlw21;
        for (char& c : jlw_tag) if (c == '.') c = 'p';
    }

    // ── Optional: METAL_ABUNDANCE_SET (default: solar) ──────────────────────
    const char* env_abund = std::getenv("METAL_ABUNDANCE_SET");
    const std::string abundance_set = (env_abund && env_abund[0] != '\0')
                                      ? std::string(env_abund) : "solar";
    chemistry::abundance::MetalSet abund{};
    try {
        abund = chemistry::abundance::get_metal_set(abundance_set);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "ERROR: %s\n", e.what());
        return 1;
    }

    // ── Optional: METAL_XNH0 [cm^-3] (initial H number density) ─────────────
    double xnH0 = 1.0;  // default: 1.0 cm^-3
    const char* env_xnH0 = std::getenv("METAL_XNH0");
    if (env_xnH0 && env_xnH0[0] != '\0') {
        try {
            xnH0 = std::stod(env_xnH0);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: METAL_XNH0='%s' is not a valid number\n", env_xnH0);
            return 1;
        }
        if (xnH0 <= 0.0) {
            std::fprintf(stderr,
                "ERROR: METAL_XNH0 must be > 0, got %s\n", env_xnH0);
            return 1;
        }
    }

    // ── Optional: METAL_TK0 [K] (initial gas temperature) ────────────────────
    double T_K0 = 1.0e2;  // default: 100 K
    const char* env_T_K0 = std::getenv("METAL_TK0");
    if (env_T_K0 && env_T_K0[0] != '\0') {
        try {
            T_K0 = std::stod(env_T_K0);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: METAL_TK0='%s' is not a valid number\n", env_T_K0);
            return 1;
        }
        if (T_K0 <= 0.0) {
            std::fprintf(stderr,
                "ERROR: METAL_TK0 must be > 0, got %s\n", env_T_K0);
            return 1;
        }
    }

    // ── Optional: METAL_YE0 (initial electron / H+ fraction) ─────────────────
    double ye0 = 1.0e-4;  // default: 1e-4
    const char* env_ye0 = std::getenv("METAL_YE0");
    if (env_ye0 && env_ye0[0] != '\0') {
        try {
            ye0 = std::stod(env_ye0);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: METAL_YE0='%s' is not a valid number\n", env_ye0);
            return 1;
        }
        if (ye0 < 0.0 || ye0 >= 1.0) {
            std::fprintf(stderr,
                "ERROR: METAL_YE0 must be in [0, 1), got %s\n", env_ye0);
            return 1;
        }
    }

    // ── Optional: METAL_YH2 (initial H2 fraction) ────────────────────────────
    double y_H2_init = 6.0e-7;  // default: 6e-7
    const char* env_yH2 = std::getenv("METAL_YH2");
    if (env_yH2 && env_yH2[0] != '\0') {
        try {
            y_H2_init = std::stod(env_yH2);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: METAL_YH2='%s' is not a valid number\n", env_yH2);
            return 1;
        }
        if (y_H2_init < 0.0 || y_H2_init >= 0.5) {
            std::fprintf(stderr,
                "ERROR: METAL_YH2 must be in [0, 0.5), got %s\n", env_yH2);
            return 1;
        }
    }

    // ── Optional: METAL_YHD (initial HD fraction) ────────────────────────────
    double y_HD_init = 4.0e-10;  // default: 4e-10
    const char* env_yHD = std::getenv("METAL_YHD");
    if (env_yHD && env_yHD[0] != '\0') {
        try {
            y_HD_init = std::stod(env_yHD);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: METAL_YHD='%s' is not a valid number\n", env_yHD);
            return 1;
        }
        if (y_HD_init < 0.0) {
            std::fprintf(stderr,
                "ERROR: METAL_YHD must be >= 0, got %s\n", env_yHD);
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
                "ERROR: METAL_OUTPUT_STRIDE='%s' is not a valid integer\n", env_stride);
            return 1;
        }
        if (output_stride <= 0) {
            std::fprintf(stderr,
                "ERROR: METAL_OUTPUT_STRIDE must be > 0, got %s\n", env_stride);
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
                "ERROR: METAL_MAX_ITER='%s' is not a valid integer\n", env_max_iter);
            return 1;
        }
        if (max_iter <= 0) {
            std::fprintf(stderr,
                "ERROR: METAL_MAX_ITER must be > 0, got %s\n", env_max_iter);
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
            std::fprintf(stderr,
                "ERROR: METAL_CR_ATTEN_COL_DENS='%s' is not a valid number\n", env_cr_col);
            return 1;
        }
        if (cr_atten_col_dens <= 0.0) {
            std::fprintf(stderr,
                "ERROR: METAL_CR_ATTEN_COL_DENS must be > 0, got %s\n", env_cr_col);
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
            std::fprintf(stderr,
                "ERROR: METAL_CR_ATTEN_SECOND_FRAC='%s' is not a valid number\n", env_cr_second);
            return 1;
        }
        if (cr_atten_second_frac < 0.0) {
            std::fprintf(stderr,
                "ERROR: METAL_CR_ATTEN_SECOND_FRAC must be >= 0, got %s\n", env_cr_second);
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
                "ERROR: METAL_CR_METAL_BKGND='%s' is not a valid number\n", env_cr_bkg);
            return 1;
        }
        if (cr_metal_bkgnd < 0.0) {
            std::fprintf(stderr,
                "ERROR: METAL_CR_METAL_BKGND must be >= 0, got %s\n", env_cr_bkg);
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
            std::fprintf(stderr,
                "ERROR: METAL_SRA_RATE='%s' is not a valid number\n", env_sra_rate);
            return 1;
        }
        if (sra_rate < 0.0) {
            std::fprintf(stderr,
                "ERROR: METAL_SRA_RATE must be >= 0, got %s\n", env_sra_rate);
            return 1;
        }
    }
    const char* env_lra_rate = std::getenv("METAL_LRA_RATE");
    if (env_lra_rate && env_lra_rate[0] != '\0') {
        try {
            lra_rate = std::stod(env_lra_rate);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: METAL_LRA_RATE='%s' is not a valid number\n", env_lra_rate);
            return 1;
        }
        if (lra_rate < 0.0) {
            std::fprintf(stderr,
                "ERROR: METAL_LRA_RATE must be >= 0, got %s\n", env_lra_rate);
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
            std::fprintf(stderr,
                "ERROR: METAL_T_CR_DES='%s' is not a valid number\n", env_tcr);
            return 1;
        }
        if (t_cr_desorp <= 0.0) {
            std::fprintf(stderr,
                "ERROR: METAL_T_CR_DES must be > 0, got %s\n", env_tcr);
            return 1;
        }
    }

    // ── Optional: METAL_C_GAS_FRAC / METAL_O_GAS_FRAC / METAL_MG_GAS_FRAC ───
    double c_gas_frac = kCgasFrac;
    double o_gas_frac = kOgasFrac;
    double mg_gas_frac = kMggasFrac;

    const char* env_cfrac = std::getenv("METAL_C_GAS_FRAC");
    if (env_cfrac && env_cfrac[0] != '\0') {
        try { c_gas_frac = std::stod(env_cfrac); }
        catch (const std::exception&) {
            std::fprintf(stderr, "ERROR: METAL_C_GAS_FRAC='%s' is not a valid number\n", env_cfrac);
            return 1;
        }
    }
    const char* env_ofrac = std::getenv("METAL_O_GAS_FRAC");
    if (env_ofrac && env_ofrac[0] != '\0') {
        try { o_gas_frac = std::stod(env_ofrac); }
        catch (const std::exception&) {
            std::fprintf(stderr, "ERROR: METAL_O_GAS_FRAC='%s' is not a valid number\n", env_ofrac);
            return 1;
        }
    }
    const char* env_mgfrac = std::getenv("METAL_MG_GAS_FRAC");
    if (env_mgfrac && env_mgfrac[0] != '\0') {
        try { mg_gas_frac = std::stod(env_mgfrac); }
        catch (const std::exception&) {
            std::fprintf(stderr, "ERROR: METAL_MG_GAS_FRAC='%s' is not a valid number\n", env_mgfrac);
            return 1;
        }
    }
    if (c_gas_frac < 0.0 || c_gas_frac > 1.0 ||
        o_gas_frac < 0.0 || o_gas_frac > 1.0 ||
        mg_gas_frac < 0.0 || mg_gas_frac > 1.0) {
        std::fprintf(stderr,
            "ERROR: gas fractions must be in [0,1]: C=%g O=%g Mg=%g\n",
            c_gas_frac, o_gas_frac, mg_gas_frac);
        return 1;
    }

    // ── Optional: METAL_DT_FACTOR / METAL_DT_FACTOR_INIT / METAL_N_INIT_STEPS ─
    double dt_factor = kDtFactor;
    double dt_factor_init = kDtFactorInit;
    int n_init_steps = kNInitSteps;

    const char* env_dt_factor = std::getenv("METAL_DT_FACTOR");
    if (env_dt_factor && env_dt_factor[0] != '\0') {
        try { dt_factor = std::stod(env_dt_factor); }
        catch (const std::exception&) {
            std::fprintf(stderr, "ERROR: METAL_DT_FACTOR='%s' is not a valid number\n", env_dt_factor);
            return 1;
        }
        if (dt_factor <= 0.0) {
            std::fprintf(stderr, "ERROR: METAL_DT_FACTOR must be > 0, got %s\n", env_dt_factor);
            return 1;
        }
    }
    const char* env_dt_factor_init = std::getenv("METAL_DT_FACTOR_INIT");
    if (env_dt_factor_init && env_dt_factor_init[0] != '\0') {
        try { dt_factor_init = std::stod(env_dt_factor_init); }
        catch (const std::exception&) {
            std::fprintf(stderr, "ERROR: METAL_DT_FACTOR_INIT='%s' is not a valid number\n", env_dt_factor_init);
            return 1;
        }
        if (dt_factor_init <= 0.0) {
            std::fprintf(stderr, "ERROR: METAL_DT_FACTOR_INIT must be > 0, got %s\n", env_dt_factor_init);
            return 1;
        }
    }
    const char* env_n_init = std::getenv("METAL_N_INIT_STEPS");
    if (env_n_init && env_n_init[0] != '\0') {
        try { n_init_steps = std::stoi(env_n_init); }
        catch (const std::exception&) {
            std::fprintf(stderr, "ERROR: METAL_N_INIT_STEPS='%s' is not a valid integer\n", env_n_init);
            return 1;
        }
        if (n_init_steps < 0) {
            std::fprintf(stderr, "ERROR: METAL_N_INIT_STEPS must be >= 0, got %s\n", env_n_init);
            return 1;
        }
    }

    // ── Optional: METAL_XNH_STOP [cm^-3] ─────────────────────────────────────
    double xnH_stop = kXnHStop;
    const char* env_xnh_stop = std::getenv("METAL_XNH_STOP");
    if (env_xnh_stop && env_xnh_stop[0] != '\0') {
        try { xnH_stop = std::stod(env_xnh_stop); }
        catch (const std::exception&) {
            std::fprintf(stderr, "ERROR: METAL_XNH_STOP='%s' is not a valid number\n", env_xnh_stop);
            return 1;
        }
        if (xnH_stop <= 0.0) {
            std::fprintf(stderr, "ERROR: METAL_XNH_STOP must be > 0, got %s\n", env_xnh_stop);
            return 1;
        }
    }

    std::printf("=== METAL_ZETA0=%s  METAL_Z_METAL=%s"
                "  cr_tag=%s  z_tag=%s  f_ret=%g  fret_tag=%s  ff_gamma=%d"
                "  METAL_REDSHIFT=%g  T_rad=%g K"
                "  abund=%s"
                "  J_LW21=%g  xnH0=%g"
                "  T_K0=%g  y_e0=%g  y_H2=%g  y_HD=%g"
                "  stride=%d  max_iter=%d"
                "  cr_col=%g  cr_sec=%g  cr_bkg=%g"
                "  sra_rate=%g lra_rate=%g"
                "  Cgas=%g Ogas=%g Mggas=%g  T_cr_des=%g"
                "  dt_factor=%g dt_factor_init=%g n_init=%d xnH_stop=%g ===\n",
                env_zeta0, env_z, cr_tag.c_str(), z_tag.c_str(), f_ret,
                fret_tag.c_str(), (int)ff_gamma, zred, T_rad, abund.name, jlw21, xnH0,
                T_K0, ye0, y_H2_init, y_HD_init,
                output_stride, max_iter,
                cr_atten_col_dens, cr_atten_second_frac, cr_metal_bkgnd,
                sra_rate, lra_rate,
                c_gas_frac, o_gas_frac, mg_gas_frac, t_cr_desorp,
                dt_factor, dt_factor_init, n_init_steps, xnH_stop);
    if (!fret_table_path.empty() && !ff_gamma)
        std::printf("  fret_table: %s  (%zu rows)\n",
                    fret_table_path.c_str(), fret_nH.size());
    RunCollapse(zeta0, cr_tag, Z_metal, z_tag, f_ret, fret_tag,
                T_rad, zred, zred_tag, abund, tbl, out_dir,
                jlw21, jlw_tag,
                cr_atten_col_dens, cr_atten_second_frac, cr_metal_bkgnd,
                sra_rate, lra_rate, t_cr_desorp,
                xnH0, T_K0, ye0, y_H2_init, y_HD_init,
                c_gas_frac, o_gas_frac, mg_gas_frac,
                fret_nH, fret_val, fret_table_path,
                dt_factor, dt_factor_init, n_init_steps, xnH_stop,
                output_stride, max_iter, ff_gamma);

    return 0;
}
