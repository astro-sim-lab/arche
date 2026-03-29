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
//   - supply xNc_H/H2/HD from column integration (TreeCol, Sobolev, etc.)
//   - supply tau_cnt/esc_cnt from the opacity field
//
// Here we approximate column densities as n × L_Jeans (same as the
// 1-D collapse app) purely for illustration.
//
// Environment variables:
//   PRIM_ZETA0    — CR ionization rate [s^-1]       (required)
//   PRIM_DATA_DIR — path to reaction data directory  (optional, default: DATA_DIR)
//   PRIM_XNH      — H number density [cm^-3]         (optional, default: 1e4)
//   PRIM_T_K      — gas temperature [K]               (optional, default: 100.0)
//   PRIM_YE0      — initial electron / H+ fraction    (optional, default: 1e-4)
//   PRIM_YH2      — initial H2 fraction               (optional, default: 6e-7)
//   PRIM_YHD      — initial HD fraction               (optional, default: 4e-10)
//   PRIM_NSTEPS   — number of integration steps       (optional, default: 200)
//   PRIM_DT       — time step [s]                     (optional, default: 1e10)
//   PRIM_CR_ATTEN_COL_DENS — CR attenuation column density scale [g cm^-2]
//                              (optional, default: 96.0, must be > 0)
//   PRIM_REDSHIFT — cosmological redshift z            (optional, default: 0.0)
//   PRIM_JLW21    — LW intensity J_21                  (optional, default: 0.0)
//   PRIM_ABUNDANCE_SET — abundance preset name         (optional, default: solar)
//                        currently supported: solar, default

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include "chemistry.h"

int main()
{
    constexpr double kCrAttenuColDens = 96.0;

    // ── Data directory ────────────────────────────────────────────────────────
    const char* env_data = std::getenv("PRIM_DATA_DIR");
    const std::string data_dir = (env_data && env_data[0] != '\0')
                                 ? env_data : DATA_DIR;

    // ── Required: PRIM_ZETA0 ─────────────────────────────────────────────────
    const char* env_zeta0 = std::getenv("PRIM_ZETA0");
    if (!env_zeta0 || env_zeta0[0] == '\0') {
        std::fprintf(stderr,
            "ERROR: PRIM_ZETA0 is required. Example: PRIM_ZETA0=1e-17\n");
        return 1;
    }
    double zeta0 = 0.0;
    try { zeta0 = std::stod(env_zeta0); }
    catch (...) {
        std::fprintf(stderr, "ERROR: PRIM_ZETA0='%s' is not a valid number\n", env_zeta0);
        return 1;
    }

    // ── Optional parameters ───────────────────────────────────────────────────
    double xnH   = 1.0e4;   // H number density [cm^-3]
    double T_K   = 100.0;   // gas temperature [K]
    int    nstep = 200;
    double dt    = 1.0e10;  // time step [s]

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

    xnH   = getenv_dbl("PRIM_XNH",    xnH);
    T_K   = getenv_dbl("PRIM_T_K",    T_K);
    nstep = getenv_int("PRIM_NSTEPS", nstep);
    dt    = getenv_dbl("PRIM_DT",     dt);
    double zred = getenv_dbl("PRIM_REDSHIFT", 0.0);
    double jlw21 = getenv_dbl("PRIM_JLW21", 0.0);

    double cr_atten_col_dens = kCrAttenuColDens;
    const char* env_cr_col = std::getenv("PRIM_CR_ATTEN_COL_DENS");
    if (env_cr_col && env_cr_col[0] != '\0') {
        try { cr_atten_col_dens = std::stod(env_cr_col); }
        catch (...) {
            std::fprintf(stderr,
                "ERROR: PRIM_CR_ATTEN_COL_DENS='%s' is not a valid number\n", env_cr_col);
            return 1;
        }
        if (!(cr_atten_col_dens > 0.0)) {
            std::fprintf(stderr,
                "ERROR: PRIM_CR_ATTEN_COL_DENS must be > 0, got %s\n", env_cr_col);
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
    const std::string abundance_set = (env_abund && env_abund[0] != '\0')
                                      ? std::string(env_abund) : "solar";
    chemistry::abundance::PrimordialSet abund{};
    try {
        abund = chemistry::abundance::get_primordial_set(abundance_set);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "ERROR: %s\n", e.what());
        return 1;
    }

    // ── Load reaction tables ──────────────────────────────────────────────────
    chemistry::ZeroMetalTable tbl = chemistry::make_zero_metal_table(data_dir);

    // ── Physical constants ────────────────────────────────────────────────────
    constexpr double kKB    = chemistry::phys::xk_B;
    constexpr double kMp    = chemistry::phys::xm_p;
    constexpr double kPi    = chemistry::phys::pi;
    constexpr double kGGrav = chemistry::phys::G;

    // ── Primordial initial conditions ─────────────────────────────────────────
    chemistry::ZeroMetalCell cell;
    auto& y = cell.state.y;

    const double y_H2 = getenv_dbl("PRIM_YH2", 6.0e-7);
    const double y_Hp = getenv_dbl("PRIM_YE0", 1.0e-4);
    const double y_Dp = 0.0;
    const double y_HD = getenv_dbl("PRIM_YHD", 4.0e-10);
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
    y[18] = abund.yLi;                       // Li+
    y[2]  = y_Hp + y_Dp + abund.yLi;        // e-

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
    params.T_rad = 2.725 * (1.0 + zred);  // CMB temperature [K]

    // ── Persistent opacity across steps ──────────────────────────────────────
    double xk_gas  = 0.0;
    double tau_cnt = 0.0;
    double esc_cnt = 1.0;

    // ── Print header ──────────────────────────────────────────────────────────
    std::printf("# primordial chemistry demo: xnH=%g cm^-3  T_K=%g K  zeta0=%g s^-1"
                "  y_e0=%g  y_H2=%g  y_HD=%g  cr_col=%g  zred=%g  T_rad=%g  J_LW21=%g\n",
                xnH, T_K, zeta0, y_Hp, y_H2, y_HD, cr_atten_col_dens, zred, params.T_rad, jlw21);
    std::printf("# abundance_set=%s\n", abund.name);
    std::printf("# %5s  %10s  %11s  %11s  %11s  %11s  %11s  %11s  %11s\n",
                "step", "t [s]",
                "y[H]", "y[H2]", "y[e-]",
                "L_line", "L_cnt", "L_ch", "Gam_CR");
    std::printf("# %5s  %10s  %11s  %11s  %11s  %11s  %11s  %11s  %11s\n",
                "-----", "----------",
                "-----------", "-----------", "-----------",
                "-----------", "-----------", "-----------", "-----------");

    double t = 0.0;

    // ── Integration loop ──────────────────────────────────────────────────────
    for (int it = 1; it <= nstep; ++it) {

        // Jeans length (isothermal estimate)
        const double xlmbd_J = std::sqrt(kPi * kKB * T_K
                                         / (kGGrav * xmu * kMp * rho));

        // Column densities (Jeans-length approximation)
        // NOTE: in a 3-D fluid code these come from column integration.
        const double xNc_H  = y[0]  * xnH * xlmbd_J;
        const double xNc_H2 = y[1]  * xnH * xlmbd_J;
        const double xNc_HD = y[12] * xnH * xlmbd_J;

        tau_cnt = xk_gas * rho * xlmbd_J;
        esc_cnt = (tau_cnt > 1.0) ? 1.0 / (tau_cnt * tau_cnt) : 1.0;

        // Build shielding environment (caller-supplied in a real fluid code)
        chemistry::ChemShielding shield;
        shield.zeta    = zeta0 * std::exp(-rho * xlmbd_J / cr_atten_col_dens);
        shield.xNc_H   = xNc_H;
        shield.xNc_H2  = xNc_H2;
        shield.xNc_HD  = xNc_HD;
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
        xk_gas = rates.xk_gas;

        // Print
        std::printf("  %5d  %10.3E  %11.4E  %11.4E  %11.4E"
                    "  %11.4E  %11.4E  %11.4E  %11.4E\n",
                    it, t,
                    y[0], y[1], y[2],
                    rates.xLmbd_line, rates.xLmbd_cnt,
                    rates.xLmbd_ch,   rates.xGam_CR);

        t += dt;
    }

    return 0;
}
