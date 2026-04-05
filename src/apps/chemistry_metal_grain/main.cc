// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
//
// main.cc — standalone metal-grain chemistry demo using chem_full_step()
//
// Demonstrates the chem_full_step() API for the metal-grain (89-species)
// network in isolation, without gravitational-collapse hydrodynamics.
//
// The gas cell is held at fixed density and temperature (isochoric,
// isothermal) so the focus is on the chemistry / cooling API.
//
// In a 3-D fluid code the caller would:
//   - compute zeta from an attenuation model or ray-tracing
//   - supply all column densities from column integration (TreeCol, Sobolev, etc.)
//   - supply tau_cnt/esc_cnt from the opacity field
//   - update T_gr_K from the previous step's params.T_gr_K
//
// Here we approximate column densities as n × L_Jeans (same as the
// 1-D collapse app) purely for illustration.
//
// Environment variables:
//   METAL_ZETA0     — CR ionization rate [s^-1]      (required)
//   METAL_Z_METAL   — metallicity [Z_sun]             (required)
//   PRIM_DATA_DIR   — path to shared primordial data  (optional, default: DATA_DIR)
//   METAL_DATA_DIR  — path to metal-grain data        (optional, default: METAL_DATA_DIR_DEF)
//   METAL_XNH       — H number density [cm^-3]        (optional, default: 1e4)
//   METAL_T_K       — gas temperature [K]              (optional, default: 100.0)
//   METAL_YE0       — initial electron / H+ fraction  (optional, default: 1e-4)
//   METAL_YH2       — initial H2 fraction              (optional, default: 6e-7)
//   METAL_YHD       — initial HD fraction              (optional, default: 4e-10)
//   METAL_NSTEPS    — number of integration steps      (optional, default: 200)
//   METAL_DT        — time step [s]                    (optional, default: 1e10)
//   METAL_CR_ATTEN_COL_DENS   — CR attenuation column density [g cm^-2]
//                                (optional, default: 96.0, must be > 0)
//   METAL_CR_ATTEN_SECOND_FRAC — secondary CR attenuation fraction
//                                (optional, default: 7.6e-2, must be >= 0)
//   METAL_CR_METAL_BKGND      — CR metal floor coefficient
//                                (optional, default: 1.4e-22, must be >= 0)
//   METAL_SRA_RATE            — short-lived radionuclide ionization scale
//                                (optional, default: 0.0, must be >= 0)
//   METAL_LRA_RATE            — long-lived radionuclide ionization scale
//                                (optional, default: 0.0, must be >= 0)
//   METAL_C_GAS_FRAC          — initial C gas fraction (optional, default: 0.28, [0,1])
//   METAL_O_GAS_FRAC          — initial O gas fraction (optional, default: 0.54, [0,1])
//   METAL_MG_GAS_FRAC         — initial Mg gas fraction(optional, default: 0.02, [0,1])
//   METAL_T_CR_DES            — effective CR desorption spike temperature [K]
//                                (optional, default: 70.0, must be > 0)
//   METAL_REDSHIFT            — cosmological redshift z (optional, default: 0.0)
//   METAL_JLW21               — LW intensity J_21       (optional, default: 0.0)
//   METAL_ABUNDANCE_SET       — abundance preset name   (optional, default: solar)
//                                currently supported: solar, default

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include "chemistry.h"

int main()
{
    constexpr double kCrAttenuColDens   = 96.0;
    constexpr double kCrAttenuSecondFrac = 7.6e-2;
    constexpr double kCrMetalBkgnd      = 1.4e-22;
    constexpr double kSraRateDefault    = 0.0;
    constexpr double kLraRateDefault    = 0.0;
    constexpr double kTcrDesorp         = 70.0;
    constexpr double kCgasFrac          = chemistry::metal_grain::c_gas_frac_default;
    constexpr double kOgasFrac          = chemistry::metal_grain::o_gas_frac_default;
    constexpr double kMggasFrac         = chemistry::metal_grain::mg_gas_frac_default;

    // ── Data directories ──────────────────────────────────────────────────────
    const char* env_prim  = std::getenv("PRIM_DATA_DIR");
    const char* env_metal = std::getenv("METAL_DATA_DIR");
    const std::string data_dir       = (env_prim  && env_prim[0]  != '\0')
                                       ? env_prim  : DATA_DIR;
    const std::string metal_data_dir = (env_metal && env_metal[0] != '\0')
                                       ? env_metal : METAL_DATA_DIR_DEF;

    // ── Required: METAL_ZETA0 ─────────────────────────────────────────────────
    const char* env_zeta0 = std::getenv("METAL_ZETA0");
    if (!env_zeta0 || env_zeta0[0] == '\0') {
        std::fprintf(stderr,
            "ERROR: METAL_ZETA0 is required. Example: METAL_ZETA0=1e-17\n");
        return 1;
    }
    double zeta0 = 0.0;
    try { zeta0 = std::stod(env_zeta0); }
    catch (...) {
        std::fprintf(stderr, "ERROR: METAL_ZETA0='%s' is not a valid number\n", env_zeta0);
        return 1;
    }

    // ── Required: METAL_Z_METAL ───────────────────────────────────────────────
    const char* env_z = std::getenv("METAL_Z_METAL");
    if (!env_z || env_z[0] == '\0') {
        std::fprintf(stderr,
            "ERROR: METAL_Z_METAL is required. Example: METAL_Z_METAL=1e-3\n");
        return 1;
    }
    double Z_metal = 0.0;
    try { Z_metal = std::stod(env_z); }
    catch (...) {
        std::fprintf(stderr, "ERROR: METAL_Z_METAL='%s' is not a valid number\n", env_z);
        return 1;
    }

    // ── Optional parameters ───────────────────────────────────────────────────
    double xnH   = 1.0e4;
    double T_K   = 100.0;
    int    nstep = 200;
    double dt    = 1.0e10;

    auto getenv_dbl = [](const char* name, double def) -> double {
        const char* s = std::getenv(name);
        if (!s || s[0] == '\0') return def;
        try { return std::stod(s); } catch (...) { return def; }
    };
    auto getenv_int = [](const char* name, int def) -> int {
        const char* s = std::getenv(name);
        if (!s || s[0] == '\0') return def;
        try { return std::stoi(s); } catch (...) { return def; }
    };

    xnH   = getenv_dbl("METAL_XNH",    xnH);
    T_K   = getenv_dbl("METAL_T_K",    T_K);
    nstep = getenv_int("METAL_NSTEPS", nstep);
    dt    = getenv_dbl("METAL_DT",     dt);
    double zred = getenv_dbl("METAL_REDSHIFT", 0.0);
    double jlw21 = getenv_dbl("METAL_JLW21", 0.0);

    double cr_atten_col_dens = kCrAttenuColDens;
    double cr_atten_second_frac = kCrAttenuSecondFrac;
    double cr_metal_bkgnd = kCrMetalBkgnd;
    double sra_rate = kSraRateDefault;
    double lra_rate = kLraRateDefault;
    double t_cr_desorp = kTcrDesorp;
    double c_gas_frac = kCgasFrac;
    double o_gas_frac = kOgasFrac;
    double mg_gas_frac = kMggasFrac;

    const char* env_cr_col = std::getenv("METAL_CR_ATTEN_COL_DENS");
    if (env_cr_col && env_cr_col[0] != '\0') {
        try { cr_atten_col_dens = std::stod(env_cr_col); }
        catch (...) {
            std::fprintf(stderr,
                "ERROR: METAL_CR_ATTEN_COL_DENS='%s' is not a valid number\n", env_cr_col);
            return 1;
        }
        if (!(cr_atten_col_dens > 0.0)) {
            std::fprintf(stderr,
                "ERROR: METAL_CR_ATTEN_COL_DENS must be > 0, got %s\n", env_cr_col);
            return 1;
        }
    }

    const char* env_cr_second = std::getenv("METAL_CR_ATTEN_SECOND_FRAC");
    if (env_cr_second && env_cr_second[0] != '\0') {
        try { cr_atten_second_frac = std::stod(env_cr_second); }
        catch (...) {
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

    const char* env_cr_bkg = std::getenv("METAL_CR_METAL_BKGND");
    if (env_cr_bkg && env_cr_bkg[0] != '\0') {
        try { cr_metal_bkgnd = std::stod(env_cr_bkg); }
        catch (...) {
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

    const char* env_sra_rate = std::getenv("METAL_SRA_RATE");
    if (env_sra_rate && env_sra_rate[0] != '\0') {
        try { sra_rate = std::stod(env_sra_rate); }
        catch (...) {
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
        try { lra_rate = std::stod(env_lra_rate); }
        catch (...) {
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

    const char* env_tcr = std::getenv("METAL_T_CR_DES");
    if (env_tcr && env_tcr[0] != '\0') {
        try { t_cr_desorp = std::stod(env_tcr); }
        catch (...) {
            std::fprintf(stderr,
                "ERROR: METAL_T_CR_DES='%s' is not a valid number\n", env_tcr);
            return 1;
        }
        if (!(t_cr_desorp > 0.0)) {
            std::fprintf(stderr,
                "ERROR: METAL_T_CR_DES must be > 0, got %s\n", env_tcr);
            return 1;
        }
    }

    auto parse_fraction = [](const char* key, double& v) -> bool {
        const char* s = std::getenv(key);
        if (!s || s[0] == '\0') return true;
        try { v = std::stod(s); }
        catch (...) {
            std::fprintf(stderr, "ERROR: %s='%s' is not a valid number\n", key, s);
            return false;
        }
        if (v < 0.0 || v > 1.0) {
            std::fprintf(stderr, "ERROR: %s must be in [0,1], got %s\n", key, s);
            return false;
        }
        return true;
    };
    if (!parse_fraction("METAL_C_GAS_FRAC", c_gas_frac)) return 1;
    if (!parse_fraction("METAL_O_GAS_FRAC", o_gas_frac)) return 1;
    if (!parse_fraction("METAL_MG_GAS_FRAC", mg_gas_frac)) return 1;
    if (zred < 0.0) {
        std::fprintf(stderr, "ERROR: METAL_REDSHIFT must be >= 0, got %g\n", zred);
        return 1;
    }
    if (jlw21 < 0.0) {
        std::fprintf(stderr, "ERROR: METAL_JLW21 must be >= 0, got %g\n", jlw21);
        return 1;
    }

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

    // ── Load reaction tables ──────────────────────────────────────────────────
    chemistry::MetalGrainTable tbl =
        chemistry::make_metal_grain_table(data_dir, metal_data_dir);

    // ── Physical constants ────────────────────────────────────────────────────
    constexpr double kKB    = chemistry::phys::xk_B;
    constexpr double kMp    = chemistry::phys::xm_p;
    constexpr double kPi    = chemistry::phys::pi;
    constexpr double kGGrav = chemistry::phys::G;

    // ── Initial conditions (metal-grain network, 89 species) ─────────────────
    chemistry::MetalGrainCell cell;
    auto& y = cell.state.y;

    const double XC  = abund.yC  * Z_metal;
    const double XO  = abund.yO  * Z_metal;
    const double XMg = abund.yMg * Z_metal;

    const double y_H2 = getenv_dbl("METAL_YH2", 6.0e-7);
    const double y_Hp = getenv_dbl("METAL_YE0", 1.0e-4);
    const double y_Dp = 0.0;
    const double y_HD = getenv_dbl("METAL_YHD", 4.0e-10);
    const double y_H  = 1.0 - y_Hp - 2.0*y_H2 - y_HD;
    if (y_H <= 0.0) {
        std::fprintf(stderr,
            "ERROR: Invalid IC: y_H = 1 - y_e0(%.3g) - 2*y_H2(%.3g) - y_HD(%.3g) = %.3g <= 0\n",
            y_Hp, y_H2, y_HD, y_H);
        return 1;
    }
    y[0]  = y_H;                             // H
    y[1]  = y_H2;                            // H2
    y[3]  = y_Hp;                            // H+
    y[7]  = abund.yHe;                       // He
    y[11] = abund.yD - y_Dp - y_HD;         // D
    y[12] = y_HD;                            // HD
    y[13] = y_Dp;                            // D+
    y[52] = abund.yLi;                       // Li+
    y[16] = XC  * c_gas_frac;               // CI
    y[29] = XO  * o_gas_frac;               // OI
    y[62] = XMg * mg_gas_frac;              // Mg+
    y[2]  = y_Hp + y_Dp + abund.yLi + y[62]; // e-

    const double rho = (1.0 + 4.0*abund.yHe) * kMp * xnH;
    double xmu  = (1.0 + 4.0*abund.yHe) / (y[0]+y[1]+y[2]+y[3]+y[7]+y[8]+y[9]);
    double gamma = 1.0 + (1.0 + 4.0*abund.yHe)
                 / (xmu * (1.5*(y[0]+y[2]+y[3]+y[7]+y[8]+y[9])
                          + chemistry::c_H2(T_K)*y[1]));

    cell.state.xnH   = xnH;
    cell.state.T_K   = T_K;
    cell.state.xmu   = xmu;
    cell.state.gamma = gamma;

    // ── Chemistry parameters ──────────────────────────────────────────────────
    chemistry::ChemParams params{};
    params.T_rad   = 2.725 * (1.0 + zred);
    params.Z_metal = Z_metal;
    params.T_gr_K  = params.T_rad;  // initial grain temperature = CMB temperature
    params.T_cr_desorp = t_cr_desorp;

    // ── Persistent opacities across steps ────────────────────────────────────
    double xk_gr   = 0.0;
    double xk_gas  = 0.0;
    double tau_cnt = 0.0;
    double esc_cnt = 1.0;

    // ── Print header ──────────────────────────────────────────────────────────
    std::printf("# metal-grain chemistry demo:"
                " xnH=%g cm^-3  T_K=%g K  Z=%g Zsun  zeta0=%g s^-1"
                "  y_e0=%g  y_H2=%g  y_HD=%g  cr_col=%g  cr_sec=%g  cr_bkg=%g"
                "  sra_rate=%g lra_rate=%g"
                "  Cgas=%g Ogas=%g Mggas=%g  T_cr_des=%g  abund=%s"
                "  zred=%g  T_rad=%g  J_LW21=%g\n",
                xnH, T_K, Z_metal, zeta0, y_Hp, y_H2, y_HD,
                cr_atten_col_dens, cr_atten_second_frac, cr_metal_bkgnd,
                sra_rate, lra_rate,
                c_gas_frac, o_gas_frac, mg_gas_frac, t_cr_desorp, abund.name,
                zred, params.T_rad, jlw21);
    std::printf("# %5s  %10s  %11s  %11s  %11s  %10s"
                "  %11s  %11s  %11s  %11s\n",
                "step", "t [s]",
                "y[H]", "y[H2]", "y[CO]", "T_gr[K]",
                "L_line", "L_cnt", "L_gr", "Gam_CR");
    std::printf("# %5s  %10s  %11s  %11s  %11s  %10s"
                "  %11s  %11s  %11s  %11s\n",
                "-----", "----------",
                "-----------", "-----------", "-----------", "----------",
                "-----------", "-----------", "-----------", "-----------");

    double t = 0.0;

    // ── Integration loop ──────────────────────────────────────────────────────
    for (int it = 1; it <= nstep; ++it) {

        // Jeans length (isothermal estimate)
        const double xlmbd_J = std::sqrt(kPi * kKB * T_K
                                         / (kGGrav * xmu * kMp * rho));

        // Column densities (Jeans-length approximation).
        // NOTE: in a 3-D fluid code these come from column integration.
        auto Nc = [&](int idx, double mass_amu) -> double {
            const double vD  = std::sqrt(2.0 * kKB * T_K / (kMp * mass_amu));
            const double lsh = std::min(xlmbd_J, 6.0 * vD * 1.0e13); // 1e13 s ~ t_Hubble
            return y[idx] * xnH * lsh;
        };

        tau_cnt = (xk_gr + xk_gas) * rho * xlmbd_J;
        esc_cnt = (tau_cnt > 1.0) ? 1.0 / (tau_cnt * tau_cnt) : 1.0;

        chemistry::ChemShielding shield;
        const double zeta_r_short = sra_rate * 7.6e-19;
        const double zeta_r_long  = lra_rate * 1.4e-22 * Z_metal;
        shield.zeta    = zeta0 * (std::exp(-rho * xlmbd_J / cr_atten_col_dens)
                       + cr_atten_second_frac)
                       + cr_metal_bkgnd * Z_metal
                       + zeta_r_short + zeta_r_long;
        shield.xNc_H   = Nc(0,  1.0);
        shield.xNc_H2  = Nc(1,  2.0);
        shield.xNc_HD  = Nc(12, 3.0);
        shield.xNc_CO  = Nc(32, 28.0);
        shield.xNc_OH  = Nc(31, 17.0);
        shield.xNc_H2O = Nc(33, 18.0);
        shield.xNc_CII = Nc(22, 12.0);
        shield.xNc_CI  = Nc(16, 12.0);
        shield.xNc_OI  = Nc(29, 16.0);
        shield.J_LW21  = jlw21;
        shield.tau_cnt = tau_cnt;
        shield.esc_cnt = esc_cnt;

        cell.state.xnH   = xnH;   // fixed (isothermal/isochoric demo)
        cell.state.T_K   = T_K;
        cell.state.xmu   = xmu;
        cell.state.gamma = gamma;

        const auto rates = chemistry::chem_full_step(cell, dt, params, shield, tbl);

        xmu    = cell.state.xmu;
        gamma  = cell.state.gamma;
        xk_gr  = rates.xk_gr;
        xk_gas = rates.xk_gas;
        // params.T_gr_K is updated in-place by chem_full_step

        // Print
        std::printf("  %5d  %10.3E  %11.4E  %11.4E  %11.4E  %10.3E"
                    "  %11.4E  %11.4E  %11.4E  %11.4E\n",
                    it, t,
                    y[0], y[1], y[32], params.T_gr_K,
                    rates.xLmbd_line, rates.xLmbd_cnt,
                    rates.xLmbd_gr,   rates.xGam_CR);

        t += dt;
    }

    return 0;
}
