// main.cc — one-zone gravitational collapse (zero-metal network)
//
// Port of:
//   fortran/zero_metal/collapse_prm.f  — PROGRAM collapse_prm + SUBROUTINE collapse
//
// Single case per run; parameters are read from environment variables:
//   PRIM_ZETA0          — CR ionization rate [s^-1]  (required)
//   PRIM_OUTDIR         — output directory            (optional, default: results/prim/h5)
//   PRIM_FF_RET         — free-fall retardation factor (optional, default: 1.0)
//   PRIM_FRET_TABLE     — 2-column ASCII file: nH[cm^-3]  f_ret  step-function table
//                         Rows sorted by ascending nH; comments with '#' allowed.
//                         If set, PRIM_FF_RET is ignored.  fret_tag = "-step"
//                         → output: collapse_CR<tag>_fret-step.h5
//   PRIM_XNH0           — initial H number density [cm^-3] (optional, default: 0.1)
//   PRIM_OUTPUT_STRIDE  — write every N-th step to HDF5  (optional, default: 100)
//   PRIM_MAX_ITER       — maximum integration steps       (optional, default: 10000000)
//   PRIM_DATA_DIR       — path to reaction data directory (optional, default: compile-time DATA_DIR)
//   PRIM_JLW21          — Lyman-Werner radiation intensity J_21 [10^-21 erg/s/cm^2/Hz/sr]
//                         (optional, default: 0.0 = no LW field)
//                         Activates H2/HD photodissociation and H- photodetachment.
//   PRIM_ABUNDANCE_SET  — abundance preset name (optional, default: solar)
//                         currently supported: solar, default
//
// HDF5 layout for each collapse_CR<tag>[_fret<fr>].h5
// ─────────────────────────────────────────────────────────────────────────────
//   Attributes (root):
//     description   — human-readable label
//     cr_tag        — CR tag string derived from PRIM_ZETA0
//     zeta0_cgs     — CR ionization rate [s^-1]
//     f_ret         — free-fall retardation factor (initial value; 1.0 = standard free-fall)
//     f_ret_table   — path to f_ret step-function table file (absent if not used)
//     network       — "zero_metal N_sp=23 N_react=140"
//     units_density — "cm^-3 (xnH), g/cm^3 (rho)"
//     units_cooling — "erg g^-1 s^-1"
//     units_time    — "s"
//     units_length  — "cm"
//     units_B       — "G"
//
//   Datasets (all 1-D length N_rows, except y):
//     step          (N_rows,)       int32   — physical step number (multiple of 100)
//     y             (N_rows, N_sp)  float64 — species abundances (number fraction / xnH)
//                     attr "species" = "H,H2,e-,H+,H2+,H3+,H-,He,He+,He++,HeH+,
//                                       D,HD,D+,HD+,D-,Li,LiH,Li+,Li-,LiH+,Li++,Li+++"
//     xnH           (N_rows,)       — H number density [cm^-3]
//     T_K           (N_rows,)       — gas temperature [K]
//     rho           (N_rows,)       — mass density [g cm^-3]
//     xLmbd_net     (N_rows,)       — net cooling Λ_line+Λ_cnt+Λ_chem − Γ_CR [erg g^-1 s^-1]
//     xLmbd_line    (N_rows,)       — total line cooling (H2+HD+Lya)
//     xLmbd_cnt     (N_rows,)       — continuum (dust+H ff+H2 CIA) cooling
//     xLmbd_ch      (N_rows,)       — chemical (endothermic reaction) cooling
//     xGam_cmp      (N_rows,)       — compressional heating p/ρ/t_eff
//     xLmbd_gas     (N_rows,)       — gas (H ff + H2 CIA) cooling subset of cnt
//     xLmbd_Lya     (N_rows,)       — Lyman-alpha cooling
//     xLmbd_H2      (N_rows,)       — H2 line cooling
//     xLmbd_HD      (N_rows,)       — HD line cooling
//     xGam_CR       (N_rows,)       — CR ionization heating
//     t_ff          (N_rows,)       — true free-fall time [s]  (= t_eff / f_ret)
//     t_cool        (N_rows,)       — cooling time e/|Λ_net| [s]
//     t_chem        (N_rows,)       — chemistry time scale [s]  min_i(y_i/|Δy_i/Δt|)
//     tau_cnt       (N_rows,)       — continuum optical depth
//     xlmbd_J       (N_rows,)       — Jeans length [cm]
//     xMJ           (N_rows,)       — Jeans mass [g]
//     B_cr          (N_rows,)       — critical (ambipolar) magnetic field [G]
//     y_plus        (N_rows,)       — total positive charge fraction
//     y_minus       (N_rows,)       — total negative charge fraction
//     charge_imbal  (N_rows,)       — |y+ − y−| / (y+ + y−)

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

namespace {

// ─── compile-time network sizes ──────────────────────────────────────────────
constexpr int kNSp    = chemistry::zero_metal::N_sp;    // 23
constexpr int kNReact = chemistry::zero_metal::N_react; // 140

// ─── physical constants (CGS) ─────────────────────────────────────────────────
constexpr double kKB     = chemistry::phys::xk_B;    // 1.380662e-16  [erg/K]
constexpr double kMp     = chemistry::phys::xm_p;    // 1.67262e-24   [g]
constexpr double kPi     = chemistry::phys::pi;
constexpr double kGGrav  = chemistry::phys::G;        // 6.6720e-8     [cm^3 g^-1 s^-2]
constexpr double kSigmaB = chemistry::phys::sigma_B;  // 5.67e-5       [erg cm^-2 s^-1 K^-4]

// ─── primordial composition ───────────────────────────────────────────────────
constexpr double kTCMB0 = 2.725;  // CMB temperature at z=0 [K]

// ─── CR reaction constants ────────────────────────────────────────────────────
// Index range in var[] for CR ionization reactions (Fortran reactions 131–136)
constexpr int    kCrReactFirst     = 130;
constexpr int    kCrReactLast      = 135;
// CR ionization energy deposited as heat: 3.4 eV × elementary charge [erg]
constexpr double kCrIonizEnergyErg = 3.4 * 1.6022e-12;
// CR attenuation column density [g cm^-2]
constexpr double kCrAttenuColDens  = 96.0;

// ─── Integration control ─────────────────────────────────────────────────────
// Timestep as fraction of min(t_cool, t_eff) after initial phase
constexpr double kDtFactor     = 1.0e-3;
// Timestep factor during the first kNInitSteps steps (short initial kick)
constexpr double kDtFactorInit = 1.0e-8;
constexpr int    kNInitSteps   = 10;
// Default HDF5 output stride (write every N-th step); overridable via PRIM_OUTPUT_STRIDE
constexpr int    kOutputStride = 100;
// Maximum number of integration steps; overridable via PRIM_MAX_ITER
constexpr int    kItMax        = 10000000;
// Stop when xnH exceeds this threshold [cm^-3]
constexpr double kXnHStop      = 1.0e23;

// ─────────────────────────────────────────────────────────────────────────────
// OutputRow — one record (every 100 steps) buffered before writing to HDF5
// ─────────────────────────────────────────────────────────────────────────────
struct OutputRow {
    int    step;
    std::array<double, kNSp> y;
    double xnH, T_K, rho;
    double xLmbd_net, xLmbd_line, xLmbd_cnt, xLmbd_ch, xGam_cmp;
    double xLmbd_gas, xLmbd_Lya, xLmbd_H2, xLmbd_HD, xGam_CR;
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
        std::fprintf(stderr, "ERROR: cannot open PRIM_FRET_TABLE '%s'\n",
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
            "ERROR: PRIM_FRET_TABLE '%s' contains no valid rows\n",
            path.c_str());
        std::exit(1);
    }
    return {nH_tab, fret_tab};
}

// ─────────────────────────────────────────────────────────────────────────────
// WriteHdf5File — write all datasets flat at the root of an open HDF5 file
// ─────────────────────────────────────────────────────────────────────────────
void WriteHdf5File(hid_t fid, const std::string& cr_tag,
                   const std::vector<OutputRow>& rows,
                   double zeta0, double f_ret,
                   double T_rad, double zred,
                   double T_K0_ic, double y_e0_ic,
                   double y_H2_ic, double y_HD_ic,
                   double jlw21 = 0.0,
                   const std::string& fret_table_path = "")
{
    const hsize_t N   = static_cast<hsize_t>(rows.size());
    const hsize_t Nsp = static_cast<hsize_t>(kNSp);

    // ── Root attributes ───────────────────────────────────────────────────────
    {
        char desc[128];
        std::snprintf(desc, sizeof(desc),
                      "zero-metal collapse CR%s (C++ port of Nakauchi+2021)",
                      cr_tag.c_str());
        H5WriteStrAttr(fid, "description",  desc);
    }
    H5WriteStrAttr(fid, "cr_tag",        cr_tag);
    H5WriteDblAttr(fid, "zeta0_cgs",     zeta0);
    H5WriteDblAttr(fid, "f_ret",         f_ret);
    if (!fret_table_path.empty())
        H5WriteStrAttr(fid, "f_ret_table", fret_table_path);
    H5WriteDblAttr(fid, "zred",          zred);
    H5WriteDblAttr(fid, "T_rad",         T_rad);
    H5WriteDblAttr(fid, "J_LW21",        jlw21);
    H5WriteDblAttr(fid, "ic_T_K0",       T_K0_ic);
    H5WriteDblAttr(fid, "ic_y_e0",       y_e0_ic);
    H5WriteDblAttr(fid, "ic_y_H2",       y_H2_ic);
    H5WriteDblAttr(fid, "ic_y_HD",       y_HD_ic);
    H5WriteStrAttr(fid, "network",       "zero_metal N_sp=23 N_react=140");
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
            "D,HD,D+,HD+,D-,Li,LiH,Li+,Li-,LiH+,Li++,Li+++");
        H5WriteStrAttr(ds, "units", "number fraction relative to xnH");
        H5Dclose(ds);
    }

    // ── scalar 1-D datasets ───────────────────────────────────────────────────
    auto dump1d = [&](std::string_view name, auto fn) {
        std::vector<double> v(N);
        for (hsize_t i = 0; i < N; ++i) v[i] = fn(rows[i]);
        H5Write1d(fid, name, v);
    };

    dump1d("xnH",          [](const OutputRow& r){ return r.xnH;          });
    dump1d("T_K",          [](const OutputRow& r){ return r.T_K;          });
    dump1d("rho",          [](const OutputRow& r){ return r.rho;          });
    dump1d("xLmbd_net",    [](const OutputRow& r){ return r.xLmbd_net;    });
    dump1d("xLmbd_line",   [](const OutputRow& r){ return r.xLmbd_line;   });
    dump1d("xLmbd_cnt",    [](const OutputRow& r){ return r.xLmbd_cnt;    });
    dump1d("xLmbd_ch",     [](const OutputRow& r){ return r.xLmbd_ch;     });
    dump1d("xGam_cmp",     [](const OutputRow& r){ return r.xGam_cmp;     });
    dump1d("xLmbd_gas",    [](const OutputRow& r){ return r.xLmbd_gas;    });
    dump1d("xLmbd_Lya",    [](const OutputRow& r){ return r.xLmbd_Lya;    });
    dump1d("xLmbd_H2",     [](const OutputRow& r){ return r.xLmbd_H2;     });
    dump1d("xLmbd_HD",     [](const OutputRow& r){ return r.xLmbd_HD;     });
    dump1d("xGam_CR",      [](const OutputRow& r){ return r.xGam_CR;      });
    dump1d("t_ff",         [](const OutputRow& r){ return r.t_ff;         });
    dump1d("t_cool",       [](const OutputRow& r){ return r.t_cool;       });
    dump1d("t_chem",       [](const OutputRow& r){ return r.t_chem;       });
    dump1d("tau_cnt",      [](const OutputRow& r){ return r.tau_cnt;      });
    dump1d("xlmbd_J",      [](const OutputRow& r){ return r.xlmbd_J;      });
    dump1d("xMJ",          [](const OutputRow& r){ return r.xMJ;          });
    dump1d("B_cr",         [](const OutputRow& r){ return r.B_cr;         });
    dump1d("y_plus",       [](const OutputRow& r){ return r.y_plus;       });
    dump1d("y_minus",      [](const OutputRow& r){ return r.y_minus;      });
    dump1d("charge_imbal", [](const OutputRow& r){ return r.charge_imbal; });
}

// ─────────────────────────────────────────────────────────────────────────────
// RunCollapse — port of Fortran subroutine collapse()
// Writes:  <out_dir>/collapse_CR<tag>.h5
//          <out_dir>/collapse_CR<tag>_fret<fr>.h5         (f_ret != 1.0)
//          <out_dir>/collapse_CR<tag>_fret<fr>_JLW<jlw>.h5  (+ J_LW21 != 0)
//          suffix order: _fret → _JLW → _z
//
// f_ret (free-fall retardation factor): scales the effective collapse
// timescale as t_eff = f_ret * t_ff.  f_ret > 1 slows the collapse;
// f_ret = 1 (default) gives standard free-fall.  Affects: density
// update (drho), compressional heating (Γ_cmp), shielding lengths
// (xlsh) and the integration timestep (dt).
// ─────────────────────────────────────────────────────────────────────────────
void RunCollapse(double y_e0, double T_K0, double y_H2_init, double y_HD_init,
                 double zeta0, double xnH1,
                 const std::string& cr_tag,
                 double f_ret_init, const std::string& fret_tag,
                 double T_rad, double zred, const std::string& zred_tag,
                 const chemistry::abundance::PrimordialSet& abund,
                 const chemistry::ZeroMetalTable& tbl,
                 const std::string& out_dir,
                 double jlw21 = 0.0,
                 const std::string& jlw_tag = "",
                 double cr_atten_col_dens = kCrAttenuColDens,
                 const std::vector<double>& fret_nH  = {},
                 const std::vector<double>& fret_val = {},
                 const std::string& fret_table_path  = "",
                 double dt_factor   = kDtFactor,
                 double dt_factor_init = kDtFactorInit,
                 int n_init_steps   = kNInitSteps,
                 double xnH_stop    = kXnHStop,
                 int output_stride = kOutputStride,
                 int max_iter      = kItMax)
{
    // ── Initial conditions ────────────────────────────────────────────────────
    const double kT1 = T_K0;        // initial temperature [K]

    chemistry::ZeroMetalCell cell;
    auto& y = cell.state.y;   // alias: y IS cell.state.y

    double y_H2  = y_H2_init;
    double y_Hp  = y_e0;
    double y_Dp  = 0.0;
    double y_HD  = y_HD_init;
    double y_H   = 1.0 - y_Hp - 2.0*y_H2 - y_HD;
    if (y_H <= 0.0) {
        std::fprintf(stderr,
            "ERROR: Invalid IC: y_H = 1 - y_e0(%.3g) - 2*y_H2(%.3g) - y_HD(%.3g) = %.3g <= 0\n",
            y_Hp, y_H2, y_HD, y_H);
        std::exit(1);
    }
    double y_D   = abund.yD - y_Dp - y_HD;
    double y_Lip = abund.yLi;

    y[0]  = y_H;    // H
    y[1]  = y_H2;   // H2
    y[3]  = y_Hp;   // H+
    y[7]  = abund.yHe;   // He
    y[11] = y_D;    // D
    y[12] = y_HD;   // HD
    y[13] = y_Dp;   // D+
    y[16] = 0.0;    // Li
    y[18] = abund.yLi;   // Li+
    y[2]  = y_Hp + y_Dp + y_Lip;   // e-
    y_e0  = y[2];

    double rho = ((1.0 + 4.0*abund.yHe) * kMp) * xnH1;
    double xnH = xnH1;
    double T_K = kT1;
    double p   = (1.0 + abund.yHe) * xnH1 * kKB * kT1;

    double xmu = (1.0 + 4.0*abund.yHe)
               / (y[0] + y[1] + y[2] + y[3] + y[7] + y[8] + y[9]);
    double gamma = 1.0 + (1.0 + 4.0*abund.yHe)
                 / (xmu * (1.5*(y[0] + y[2] + y[3] + y[7] + y[8] + y[9])
                          + chemistry::c_H2(T_K) * y[1]));
    double e = kKB * T_K / ((gamma - 1.0) * xmu * kMp);

    chemistry::ChemParams params{};
    params.T_rad = T_rad;
    // params.zeta is set each step via ChemShielding

    double t        = 0.0;
    double dt       = 1.0e-1;
    double t_chem   = 1.0e-1;
    double t_cool   = 0.0;
    double esc_cnt  = 1.0;
    double tau_cnt  = 0.0;
    double xk_gas   = 0.0;  // gas opacity [cm^2/g]; persists across steps
    double xNc_H    = 0.0;
    double xNc_H2   = 0.0;
    double xNc_HD   = 0.0;

    // ── f_ret step-function table state ──────────────────────────────────────
    double f_ret = f_ret_init;
    int    fret_idx          = 0;
    const bool has_fret_tab  = !fret_nH.empty();

    // ── HDF5 row buffer ───────────────────────────────────────────────────────
    std::vector<OutputRow> h5rows;
    h5rows.reserve(600);

    // ── Time integration ──────────────────────────────────────────────────────
    for (int it = 1; it <= max_iter; ++it) {

        // ── f_ret ratchet: advance step-function when xnH crosses next threshold ─
        if (has_fret_tab) {
            while (fret_idx + 1 < static_cast<int>(fret_nH.size()) &&
                   xnH >= fret_nH[fret_idx + 1]) {
                ++fret_idx;
                f_ret = fret_val[fret_idx];
            }
        }

        double t_ff  = std::sqrt(3.0*kPi / (32.0*kGGrav*rho));
        double t_eff = f_ret * t_ff;   // effective collapse timescale [s]

        double vD_H  = std::sqrt(2.0*kKB*T_K / kMp);
        double vD_H2 = std::sqrt(2.0*kKB*T_K / (kMp * 2.0));
        double vD_HD = std::sqrt(2.0*kKB*T_K / (kMp * 3.0));

        double xlmbd_J = std::sqrt(kPi * kKB * T_K / (kGGrav * xmu * kMp * rho));
        double xlsh_H  = std::min(xlmbd_J, 6.0*vD_H  * t_eff);
        double xlsh_H2 = std::min(xlmbd_J, 6.0*vD_H2 * t_eff);
        double xlsh_HD = std::min(xlmbd_J, 6.0*vD_HD * t_eff);

        xNc_H  = y[0]  * xnH * xlsh_H;
        xNc_H2 = y[1]  * xnH * xlsh_H2;
        xNc_HD = y[12] * xnH * xlsh_HD;

        tau_cnt = xk_gas * rho * xlmbd_J;
        esc_cnt = (tau_cnt > 1.0) ? 1.0 / (tau_cnt * tau_cnt) : 1.0;

        // ── Build shielding environment ───────────────────────────────────────
        chemistry::ChemShielding shield;
        shield.zeta    = zeta0 * std::exp(-rho * xlmbd_J / cr_atten_col_dens);
        shield.xNc_H   = xNc_H;
        shield.xNc_H2  = xNc_H2;
        shield.xNc_HD  = xNc_HD;
        shield.tau_cnt = tau_cnt;
        shield.esc_cnt = esc_cnt;
        shield.J_LW21  = jlw21;

        cell.state.xnH   = xnH;
        cell.state.T_K   = T_K;
        cell.state.xmu   = xmu;
        cell.state.gamma = gamma;

        const auto y_prev = y;
        const auto rates  = chemistry::chem_full_step(cell, dt, params, shield, tbl);

        xmu   = cell.state.xmu;
        gamma = cell.state.gamma;
        xk_gas = rates.xk_gas;

        const double xLmbd_line = rates.xLmbd_line;
        const double xLmbd_H2   = rates.xLmbd_H2;
        const double xLmbd_HD   = rates.xLmbd_HD;
        const double xLmbd_Lya  = rates.xLmbd_Lya;
        const double xLmbd_gas  = rates.xLmbd_gas;
        const double xLmbd_cnt  = rates.xLmbd_cnt;
        const double xLmbd_ch   = rates.xLmbd_ch;
        const double xGam_CR    = rates.xGam_CR;
        const double xLmbd_net  = rates.xLmbd_net;

        // ── Chemistry timescale: min_i( y_i / |Δy_i / Δt| ) ─────────────────
        t_chem = 1.0e50;
        for (int i = 0; i < kNSp; ++i) {
            const double dy = std::abs(y[i] - y_prev[i]);
            if (y[i] > 1.0e-30 && dy > 1.0e-40)
                t_chem = std::min(t_chem, dt * y[i] / dy);
        }
        t_chem = std::max(t_chem, dt);
        double xGam_cmp  = p / rho / t_eff;

        double xMJ  = (4.0*kPi/3.0) * rho * xlmbd_J*xlmbd_J*xlmbd_J;
        double B_cr = std::sqrt(4.0*kPi * kGGrav * xMJ * rho / xlmbd_J);

        double y_plus  = y[3] + y[4] + y[5] + y[8] + y[10] + y[13] + y[14]
                       + y[18] + y[20]
                       + 2.0*(y[9] + y[21]) + 3.0*y[22];
        double y_minus = y[2] + y[6] + y[15] + y[19];
        double y_sum   = y_plus + y_minus;
        double charge_imbal = (y_sum > 0.0)
            ? std::abs(y_plus - y_minus) / y_sum : 0.0;

        // ── Output (every output_stride steps) ───────────────────────────────
        if (it % output_stride == 0) {
            OutputRow row{};
            row.step         = it;
            row.y            = y;
            row.xnH          = xnH;
            row.T_K          = T_K;
            row.rho          = rho;
            row.xLmbd_net    = xLmbd_net;
            row.xLmbd_line   = xLmbd_line;
            row.xLmbd_cnt    = xLmbd_cnt;
            row.xLmbd_ch     = xLmbd_ch;
            row.xGam_cmp     = xGam_cmp;
            row.xLmbd_gas    = xLmbd_gas;
            row.xLmbd_Lya    = xLmbd_Lya;
            row.xLmbd_H2     = xLmbd_H2;
            row.xLmbd_HD     = xLmbd_HD;
            row.xGam_CR      = xGam_CR;
            row.t_ff         = t_ff;
            row.t_cool       = t_cool;
            row.t_chem       = t_chem;
            row.tau_cnt      = tau_cnt;
            row.xlmbd_J      = xlmbd_J;
            row.xMJ          = xMJ;
            row.B_cr         = B_cr;
            row.y_plus       = y_plus;
            row.y_minus      = y_minus;
            row.charge_imbal = charge_imbal;
            h5rows.push_back(row);
        }

        // ── Update thermodynamic state ────────────────────────────────────────
        double drho = dt * rho / t_eff;
        double de   = -xLmbd_net * dt + drho * p / (rho * rho);
        rho  += drho;
        e    += de;
        T_K   = e * ((gamma - 1.0) * xmu * kMp) / kKB;
        p     = rho * kKB * T_K / (xmu * kMp);
        xnH   = rho / ((1.0 + 4.0*abund.yHe) * kMp);
        t    += dt;

        t_cool = (xLmbd_net != 0.0) ? e / std::abs(xLmbd_net) : 1.0e50;
        if (it <= n_init_steps)
            dt = dt_factor_init * t_eff;
        else
            dt = dt_factor * std::min(t_cool, t_eff);

        // Stdout progress
        std::printf("%7d %11.3E %11.3E %11.3E %11.3E %11.3E %11.3E"
                    " %11.3E %11.3E %11.3E\n",
                    it, y_e0, zeta0 / 1.0e-17, xnH, T_K, y[2], y[1],
                    y_plus, y_minus, charge_imbal);

        if (xnH > xnH_stop || T_K > 1.0e5 ||
            !std::isfinite(xnH) || !std::isfinite(T_K) || e <= 0.0) break;
    }

    // ── Write HDF5 file ───────────────────────────────────────────────────────
    std::string h5_path = out_dir + "/collapse_CR" + cr_tag;
    if (!fret_tag.empty()) h5_path += "_fret" + fret_tag;
    if (!jlw_tag.empty())  h5_path += "_JLW"  + jlw_tag;
    if (!zred_tag.empty()) h5_path += "_z"    + zred_tag;
    h5_path += ".h5";
    hid_t fid = H5Create(h5_path);
    if (fid < 0) {
        std::fprintf(stderr, "ERROR: cannot create %s\n", h5_path.c_str());
        return;
    }
    WriteHdf5File(fid, cr_tag, h5rows, zeta0, f_ret_init, T_rad, zred,
                  T_K0, y_e0, y_H2_init, y_HD_init,
                  jlw21, fret_table_path);
    H5Fclose(fid);

    std::printf("  -> %s  (%zu rows)\n", h5_path.c_str(), h5rows.size());
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
int main()
{
    // ── Data directory (compile-time default; overridable at runtime) ─────────
    const char* env_data_dir = std::getenv("PRIM_DATA_DIR");
    const std::string data_dir = (env_data_dir && env_data_dir[0] != '\0')
                                 ? env_data_dir : DATA_DIR;

    // ── Output directory ──────────────────────────────────────────────────────
    const char* env_out = std::getenv("PRIM_OUTDIR");
    std::string out_dir = env_out ? env_out : "results/prim/h5";
    std::filesystem::create_directories(out_dir);
    std::printf("Output directory: %s/\n", out_dir.c_str());

    // ── Load reaction tables ──────────────────────────────────────────────────
    chemistry::ZeroMetalTable tbl;
    chemistry::load_reaction_table(tbl, data_dir + "/react_prm.dat");
    chemistry::load_mass_table    (tbl, data_dir + "/mass.dat");
    chemistry::load_saha_table    (tbl, data_dir + "/react_prm_saha.dat");

    chemistry::load_partition_functions(tbl, {
        { 6,  data_dir + "/pf_H3p.dat"  },
        {13,  data_dir + "/pf_HD.dat"   },
        {15,  data_dir + "/pf_HDp.dat"  },
        {18,  data_dir + "/pf_LiH.dat"  },
        {21,  data_dir + "/pf_LiHp.dat" }
    });

    // ── Required: PRIM_ZETA0 [s^-1] ──────────────────────────────────────────
    const char* env_zeta0 = std::getenv("PRIM_ZETA0");
    if (!env_zeta0 || env_zeta0[0] == '\0') {
        std::fprintf(stderr,
            "ERROR: environment variable PRIM_ZETA0 is required\n"
            "       Set the CR ionization rate [s^-1], e.g.: PRIM_ZETA0=1e-17\n");
        return 1;
    }
    double zeta0 = 0.0;
    try {
        zeta0 = std::stod(env_zeta0);
    } catch (const std::exception&) {
        std::fprintf(stderr,
            "ERROR: PRIM_ZETA0='%s' is not a valid number\n", env_zeta0);
        return 1;
    }
    if (zeta0 < 0.0) {
        std::fprintf(stderr,
            "ERROR: PRIM_ZETA0 must be >= 0, got %s\n", env_zeta0);
        return 1;
    }

    // ── Tag: input string with '.' -> 'p'; any zero value -> "0" ─────────────
    std::string cr_tag = (zeta0 == 0.0) ? "0" : env_zeta0;
    for (char& c : cr_tag) if (c == '.') c = 'p';

    // ── Optional: PRIM_REDSHIFT (cosmological redshift, default 0.0) ──────────
    // T_rad = kTCMB0 * (1 + z).  Affects H2/HD line cooling, continuum
    // cooling (T^4 - T_rad^4), molecular cooling CMB correction, and grain
    // temperature equilibrium.
    double zred = 0.0;
    std::string zred_tag;   // empty → no suffix in filename
    const char* env_zred = std::getenv("PRIM_REDSHIFT");
    if (env_zred && env_zred[0] != '\0') {
        try {
            zred = std::stod(env_zred);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: PRIM_REDSHIFT='%s' is not a valid number\n", env_zred);
            return 1;
        }
        if (zred < 0.0) {
            std::fprintf(stderr,
                "ERROR: PRIM_REDSHIFT must be >= 0, got %s\n", env_zred);
            return 1;
        }
        if (zred != 0.0) {
            zred_tag = env_zred;
            for (char& c : zred_tag) if (c == '.') c = 'p';
        }
    }
    const double T_rad = kTCMB0 * (1.0 + zred);

    // ── Optional: PRIM_FRET_TABLE (step-function f_ret table; overrides PRIM_FF_RET) ──
    // ASCII 2-col file: nH[cm^-3]  f_ret  (sorted by ascending nH, '#' comments allowed)
    const char* env_fret_table = std::getenv("PRIM_FRET_TABLE");
    std::string fret_table_path;
    std::vector<double> fret_nH, fret_val;
    if (env_fret_table && env_fret_table[0] != '\0') {
        fret_table_path = env_fret_table;
        auto [nh, fv] = LoadFretTable(fret_table_path);
        fret_nH  = std::move(nh);
        fret_val = std::move(fv);
    }

    // ── Optional: PRIM_FF_RET (free-fall retardation factor, default 1.0) ────
    // f_ret > 1 slows the collapse; f_ret = 1 gives standard free-fall.
    // Ignored when PRIM_FRET_TABLE is set.
    const char* env_fret = std::getenv("PRIM_FF_RET");
    double f_ret = 1.0;
    std::string fret_tag;   // empty → no suffix in filename
    if (!fret_table_path.empty()) {
        // Table mode: initial f_ret from first table entry; fixed tag "-step"
        // → output filename: collapse_CR<tag>_fret-step.h5
        f_ret    = fret_val[0];
        fret_tag = "-step";
        if (env_fret && env_fret[0] != '\0')
            std::fprintf(stderr,
                "WARNING: PRIM_FRET_TABLE is set; PRIM_FF_RET='%s' is ignored\n",
                env_fret);
    } else if (env_fret && env_fret[0] != '\0') {
        try {
            f_ret = std::stod(env_fret);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: PRIM_FF_RET='%s' is not a valid number\n", env_fret);
            return 1;
        }
        if (f_ret <= 0.0) {
            std::fprintf(stderr,
                "ERROR: PRIM_FF_RET must be > 0, got %s\n", env_fret);
            return 1;
        }
        if (f_ret != 1.0) {
            fret_tag = env_fret;
            for (char& c : fret_tag) if (c == '.') c = 'p';
        }
    }

    // ── Optional: PRIM_JLW21 (Lyman-Werner intensity J_21, default 0.0) ──────
    // J_21 units: 10^-21 erg/s/cm^2/Hz/sr.  0.0 = no LW field (default).
    double jlw21 = 0.0;
    const char* env_jlw21 = std::getenv("PRIM_JLW21");
    if (env_jlw21 && env_jlw21[0] != '\0') {
        try {
            jlw21 = std::stod(env_jlw21);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: PRIM_JLW21='%s' is not a valid number\n", env_jlw21);
            return 1;
        }
        if (jlw21 < 0.0) {
            std::fprintf(stderr,
                "ERROR: PRIM_JLW21 must be >= 0, got %s\n", env_jlw21);
            return 1;
        }
    }
    // Tag: "." → "p"; zero → empty (no suffix)
    std::string jlw_tag;
    if (jlw21 != 0.0) {
        jlw_tag = env_jlw21;
        for (char& c : jlw_tag) if (c == '.') c = 'p';
    }

    // ── Optional: PRIM_ABUNDANCE_SET (default: solar) ───────────────────────
    const char* env_abund = std::getenv("PRIM_ABUNDANCE_SET");
    const std::string abundance_set = (env_abund && env_abund[0] != '\0')
                                      ? std::string(env_abund) : "solar";
    chemistry::abundance::PrimordialSet abund{};
    try {
        abund = chemistry::abundance::get_primordial_set(abundance_set);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "ERROR: %s\n", e.what());
        return 1;
    }

    // ── Optional: PRIM_XNH0 [cm^-3] (initial H number density) ──────────────
    double xnH0 = 0.1;  // default: 0.1 cm^-3
    const char* env_xnH0 = std::getenv("PRIM_XNH0");
    if (env_xnH0 && env_xnH0[0] != '\0') {
        try {
            xnH0 = std::stod(env_xnH0);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: PRIM_XNH0='%s' is not a valid number\n", env_xnH0);
            return 1;
        }
        if (xnH0 <= 0.0) {
            std::fprintf(stderr,
                "ERROR: PRIM_XNH0 must be > 0, got %s\n", env_xnH0);
            return 1;
        }
    }

    // ── Optional: PRIM_TK0 [K] (initial gas temperature) ─────────────────────
    double T_K0 = 1.0e2;  // default: 100 K
    const char* env_T_K0 = std::getenv("PRIM_TK0");
    if (env_T_K0 && env_T_K0[0] != '\0') {
        try {
            T_K0 = std::stod(env_T_K0);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: PRIM_TK0='%s' is not a valid number\n", env_T_K0);
            return 1;
        }
        if (T_K0 <= 0.0) {
            std::fprintf(stderr,
                "ERROR: PRIM_TK0 must be > 0, got %s\n", env_T_K0);
            return 1;
        }
    }

    // ── Optional: PRIM_YE0 (initial electron / H+ fraction) ──────────────────
    double ye0 = 1.0e-4;  // default: 1e-4
    const char* env_ye0 = std::getenv("PRIM_YE0");
    if (env_ye0 && env_ye0[0] != '\0') {
        try {
            ye0 = std::stod(env_ye0);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: PRIM_YE0='%s' is not a valid number\n", env_ye0);
            return 1;
        }
        if (ye0 < 0.0 || ye0 >= 1.0) {
            std::fprintf(stderr,
                "ERROR: PRIM_YE0 must be in [0, 1), got %s\n", env_ye0);
            return 1;
        }
    }

    // ── Optional: PRIM_YH2 (initial H2 fraction) ─────────────────────────────
    double y_H2_init = 6.0e-7;  // default: 6e-7
    const char* env_yH2 = std::getenv("PRIM_YH2");
    if (env_yH2 && env_yH2[0] != '\0') {
        try {
            y_H2_init = std::stod(env_yH2);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: PRIM_YH2='%s' is not a valid number\n", env_yH2);
            return 1;
        }
        if (y_H2_init < 0.0 || y_H2_init >= 0.5) {
            std::fprintf(stderr,
                "ERROR: PRIM_YH2 must be in [0, 0.5), got %s\n", env_yH2);
            return 1;
        }
    }

    // ── Optional: PRIM_YHD (initial HD fraction) ─────────────────────────────
    double y_HD_init = 4.0e-10;  // default: 4e-10
    const char* env_yHD = std::getenv("PRIM_YHD");
    if (env_yHD && env_yHD[0] != '\0') {
        try {
            y_HD_init = std::stod(env_yHD);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: PRIM_YHD='%s' is not a valid number\n", env_yHD);
            return 1;
        }
        if (y_HD_init < 0.0) {
            std::fprintf(stderr,
                "ERROR: PRIM_YHD must be >= 0, got %s\n", env_yHD);
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
                "ERROR: PRIM_OUTPUT_STRIDE='%s' is not a valid integer\n", env_stride);
            return 1;
        }
        if (output_stride <= 0) {
            std::fprintf(stderr,
                "ERROR: PRIM_OUTPUT_STRIDE must be > 0, got %s\n", env_stride);
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
            std::fprintf(stderr,
                "ERROR: PRIM_MAX_ITER='%s' is not a valid integer\n", env_max_iter);
            return 1;
        }
        if (max_iter <= 0) {
            std::fprintf(stderr,
                "ERROR: PRIM_MAX_ITER must be > 0, got %s\n", env_max_iter);
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
                "ERROR: PRIM_CR_ATTEN_COL_DENS='%s' is not a valid number\n", env_cr_col);
            return 1;
        }
        if (cr_atten_col_dens <= 0.0) {
            std::fprintf(stderr,
                "ERROR: PRIM_CR_ATTEN_COL_DENS must be > 0, got %s\n", env_cr_col);
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
            std::fprintf(stderr,
                "ERROR: PRIM_DT_FACTOR='%s' is not a valid number\n", env_dt_factor);
            return 1;
        }
        if (dt_factor <= 0.0) {
            std::fprintf(stderr,
                "ERROR: PRIM_DT_FACTOR must be > 0, got %s\n", env_dt_factor);
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
                "ERROR: PRIM_DT_FACTOR_INIT='%s' is not a valid number\n", env_dt_factor_init);
            return 1;
        }
        if (dt_factor_init <= 0.0) {
            std::fprintf(stderr,
                "ERROR: PRIM_DT_FACTOR_INIT must be > 0, got %s\n", env_dt_factor_init);
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
                "ERROR: PRIM_N_INIT_STEPS='%s' is not a valid integer\n", env_n_init);
            return 1;
        }
        if (n_init_steps < 0) {
            std::fprintf(stderr,
                "ERROR: PRIM_N_INIT_STEPS must be >= 0, got %s\n", env_n_init);
            return 1;
        }
    }

    // ── Optional: PRIM_XNH_STOP [cm^-3] ──────────────────────────────────────
    double xnH_stop = kXnHStop;
    const char* env_xnh_stop = std::getenv("PRIM_XNH_STOP");
    if (env_xnh_stop && env_xnh_stop[0] != '\0') {
        try {
            xnH_stop = std::stod(env_xnh_stop);
        } catch (const std::exception&) {
            std::fprintf(stderr,
                "ERROR: PRIM_XNH_STOP='%s' is not a valid number\n", env_xnh_stop);
            return 1;
        }
        if (xnH_stop <= 0.0) {
            std::fprintf(stderr,
                "ERROR: PRIM_XNH_STOP must be > 0, got %s\n", env_xnh_stop);
            return 1;
        }
    }

    std::printf("=== PRIM_ZETA0=%s  tag=%s  f_ret=%g  fret_tag=%s"
                "  PRIM_REDSHIFT=%g  T_rad=%g K"
                "  abund=%s"
                "  J_LW21=%g  xnH0=%g"
                "  T_K0=%g  y_e0=%g  y_H2=%g  y_HD=%g"
                "  stride=%d  max_iter=%d"
                "  cr_col=%g  dt_factor=%g  dt_factor_init=%g"
                "  n_init=%d  xnH_stop=%g ===\n",
                env_zeta0, cr_tag.c_str(), f_ret, fret_tag.c_str(),
                zred, T_rad, abund.name, jlw21, xnH0,
                T_K0, ye0, y_H2_init, y_HD_init,
                output_stride, max_iter,
                cr_atten_col_dens, dt_factor, dt_factor_init,
                n_init_steps, xnH_stop);
    if (!fret_table_path.empty())
        std::printf("  fret_table: %s  (%zu rows)\n",
                    fret_table_path.c_str(), fret_nH.size());
    RunCollapse(ye0, T_K0, y_H2_init, y_HD_init,
                zeta0, xnH0, cr_tag, f_ret, fret_tag,
                T_rad, zred, zred_tag, abund, tbl, out_dir,
                jlw21, jlw_tag, cr_atten_col_dens,
                fret_nH, fret_val, fret_table_path,
                dt_factor, dt_factor_init, n_init_steps, xnH_stop,
                output_stride, max_iter);

    return 0;
}
