// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// chemistry.h — portable chemistry module umbrella header
//
// Provides a simplified interface for embedding the chemistry network
// into a fluid/hydrodynamics code. Key components:
//
//   ChemCell<N_sp, N_react>   — per-cell state (ChemState + var cache)
//   ZeroMetalCell             — alias for ChemCell<23, 140>
//   MetalGrainCell            — alias for ChemCell<89, 1200>
//
//   chem_step()               — single time-step chemistry wrapper
//
//   make_zero_metal_table()   — load ZeroMetalTable from data files
//   make_metal_grain_table()  — load MetalGrainTable from data files
//
// Minimal usage (zero-metal network):
// ─────────────────────────────────────────────────────────────────────────
//   #include "chemistry.h"
//   using namespace chemistry;
//
//   // Load tables once at startup (shared across cells/threads)
//   ZeroMetalTable tbl = make_zero_metal_table("/path/to/data");
//
//   // Allocate per-cell state (heap for OpenMP/GPU parallel loops)
//   auto cell = std::make_unique<ZeroMetalCell>();
//   cell->state.xnH = 1.0e2;  // H number density [cm^-3]
//   cell->state.T_K = 100.0;  // temperature [K]
//   // ... initialise cell->state.y[0..22] ...
//   cell->state.xmu   = 1.22;   // mean molecular weight [m_p]
//   cell->state.gamma = 5.0/3.0;
//
//   ChemParams params;
//   params.zeta  = 1.0e-17;   // CR ionization rate [s^-1]
//   params.T_rad = 2.725;     // CMB temperature [K]
//
//   // Per time step:
//   ChemRates rates = chem_step(*cell, dt, params, tbl);
//   // rates.xLmbdch — chemistry cooling [erg g^-1 s^-1]
//   // rates.xGam_CR — CR heating        [erg g^-1 s^-1]
// ─────────────────────────────────────────────────────────────────────────
//
// Note on stack vs heap:
//   MetalGrainCell contains var[2400] (~19 KB) + state (~0.8 KB).
//   Inside chem_step, chemreact<89,1200> also allocates A_mat[89×89]
//   and dr_fdy[89×89] (~63 KB each) on the stack.
//   For OpenMP: declare MetalGrainCell with thread-private storage or heap.
//   For GPU:    allocate cells in device memory; see Phase 6 notes.
// ---------------------------------------------------------------------------

#include "solver_metal.h"      // metal_grain specializations; pulls in solver.h,
                                 // reaction.h, reaction_metal.h, partition_function.h
#include "cooling.h"           // line_cool, cnt_cool
#include "cooling_molecular.h" // line_cool_metal, LineCoolRates, EscapeState
#include "cooling_grain.h"     // cnt_cool_metal
#include "reaction_table.h"    // ReactionTable + load_* helpers
#include "state.h"             // ChemState, ChemParams, ChemRates, ChemShielding,
                                 // ChemFullRates, phys::
#include "species.h"           // zero_metal::, metal_grain:: constants
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace chemistry {

// ---------------------------------------------------------------------------
// ChemCell<N_sp, N_react>
//
// Bundles a per-cell chemical state with the inter-step reaction rate cache
// (var[2*N_react]).  The solver uses var[] to carry forward rate estimates
// between consecutive time steps (predictor-corrector warm start).
//
// Memory layout:
//   state  — ChemState<N_sp>  (y[N_sp], xnH, T_K, xmu, gamma)
//   var    — double[2*N_react] (inter-step rate cache; zero-init on creation)
//
// Thread safety: each thread/cell must own its own ChemCell.
// ReactionTable and ChemParams are read-only and can be shared.
// ---------------------------------------------------------------------------
template<int N_sp, int N_react>
struct ChemCell {
    ChemState<N_sp>               state{};
    std::array<double, 2*N_react> var{};

    // Zero the inter-step rate cache.
    // Call once after initialising state (or when copying a cell to a new
    // location where the previous rate history is not applicable).
    void reset_var() noexcept { var.fill(0.0); }
};

// Convenience aliases matching the two supported networks
using ZeroMetalCell  = ChemCell<zero_metal::N_sp,  zero_metal::N_react>;

// ---------------------------------------------------------------------------
// ChemCell<89, 1200>  full specialization — metal_grain network
//
// Extends the primary template with:
//   es      — per-cell escape-probability arrays (warm-start NR solver).
//             Must be initialised to 1.0 before first use; default-init is
//             correct (EscapeState members are initialised to {1.0, …}).
//
// T_gr_K, xk_gr, xk_gas are stored in ChemParams (passed to chem_full_step)
// rather than here, so they remain visible to the calling code for
// tau_cnt / esc_cnt computation between steps.
// ---------------------------------------------------------------------------
template<>
struct ChemCell<metal_grain::N_sp, metal_grain::N_react> {
    ChemState<metal_grain::N_sp>               state{};
    std::array<double, 2*metal_grain::N_react> var{};
    EscapeState                                es{};   // escape fractions (warm-start NR)

    void reset_var() noexcept { var.fill(0.0); }
};

using MetalGrainCell = ChemCell<metal_grain::N_sp, metal_grain::N_react>;

// ---------------------------------------------------------------------------
// chem_step()
//
// Advances the chemistry network by one time step dt for a single cell.
//
// On entry:
//   cell.state.{xnH, T_K, y, xmu, gamma}  — current cell state
//   cell.var                               — rate cache from previous step
//                                            (zero-filled on first call)
//   dt      — time step [s]
//   params  — external parameters (zeta, T_rad, Z_metal, ...)
//   tbl     — reaction table (read-only, shared)
//
// On exit (updated in-place):
//   cell.state.y      — species abundances [dimensionless, / xnH]
//   cell.state.xmu    — updated mean molecular weight [m_p]
//   cell.state.gamma  — updated adiabatic index
//   cell.var          — rate cache for the next time step
//
// Returns:
//   ChemRates.xLmbdch — net chemistry cooling rate [erg g^-1 s^-1]
//   ChemRates.xGam_CR — CR ionization heating rate  [erg g^-1 s^-1]
//
// The returned rates should be added to the fluid energy equation:
//   de/dt = - xLmbdch + xGam_CR + (line cooling) + (continuum cooling) + ...
//
// Thread safety: cell is per-thread; tbl and params are read-only.
// ---------------------------------------------------------------------------
template<int N_sp, int N_react>
ChemRates chem_step(ChemCell<N_sp, N_react>& cell,
                    double dt,
                    const ChemParams& params,
                    const ReactionTable<N_sp, N_react>& tbl)
{
    if (!std::isfinite(params.T_rad) || params.T_rad <= 0.0) {
        throw std::invalid_argument(
            "chem_step: params.T_rad must be set to a finite positive value [K]");
    }
    if constexpr (N_react == metal_grain::N_react) {
        if (!std::isfinite(params.T_gr_K) || params.T_gr_K <= 0.0) {
            throw std::invalid_argument(
                "chem_step: params.T_gr_K must be set to a finite positive value [K] for metal_grain");
        }
    }

    double xLmbdch = 0.0;

    chemcool<N_sp, N_react>(cell.state.xnH, cell.state.T_K,
                             cell.state.y, dt,
                             cell.state.xmu, cell.state.gamma,
                             xLmbdch, cell.var, tbl, params);

    // CR ionization heating: sum per-reaction CR rates stored in var[],
    // then convert to specific energy rate [erg g^-1 s^-1].
    //
    // Index ranges differ by network (Fortran CR reaction block positions):
    //   zero_metal  (N_react=140):  var[130..135]  (6 terms: CR[r1..r6])
    //   metal_grain (N_react=1200): var[543..551]  (9 terms: CR[r1..r9])
    double rate_CR = 0.0;
    if constexpr (N_react == zero_metal::N_react) {
        for (int i = 130; i <= 135; ++i) rate_CR += cell.var[i];
    } else {
        for (int i = 543; i <= 551; ++i) rate_CR += cell.var[i];
    }

    // 3.4 eV per CR ionization event; yHe = 8.33e-2 in both networks.
    // phys::xm_p is used here (app-layer precision, 1.67262e-24 g).
    // Note: inside chemcool, xm_p = 1.67e-24 (Fortran update mimicry).
    constexpr double yHe = abundance_ref::yHe;
    double xGam_CR = 3.4 * phys::eV_to_erg * rate_CR
                   / ((1.0 + 4.0 * yHe) * phys::xm_p);

    return ChemRates{ xLmbdch, xGam_CR };
}

// ---------------------------------------------------------------------------
// chem_full_step()
//
// Advances the chemistry network by one time step AND computes all
// cooling/heating rate components for a single cell.  This is the
// recommended entry point for fluid/hydro code integration.
//
// On entry:
//   cell.state.{xnH, T_K, y, xmu, gamma}  — current cell state
//   cell.var                               — rate cache from previous step
//   cell.es                                — escape fractions (metal only)
//   dt      — time step [s]
//   params  — external params (T_rad, Z_metal, T_gr_K, ...)
//             params.zeta is IGNORED; shield.zeta is used instead.
//             params.T_gr_K is UPDATED in-place (metal only).
//   shield  — shielding environment pre-computed by the caller:
//             zeta (pre-attenuated), column densities, tau_cnt, esc_cnt
//   tbl     — reaction table (read-only, shared)
//
// On exit (updated in-place):
//   cell.state.y, xmu, gamma  — updated chemical state
//   cell.var                  — rate cache for the next step
//   params.T_gr_K             — updated grain temperature (metal only)
//
// Returns:
//   ChemFullRates with all cooling/heating rate components.
//   rates.xk_gas (and rates.xk_gr for metal) are the continuum opacities
//   needed to compute tau_cnt in the next step:
//     tau_cnt = (xk_gr + xk_gas) * rho * L_shield
//
// Thread safety: cell and params are per-thread; tbl is read-only.
// ---------------------------------------------------------------------------
template<int N_sp, int N_react>
ChemFullRates chem_full_step(ChemCell<N_sp, N_react>& cell,
                             double dt,
                             ChemParams& params,
                             const ChemShielding& shield,
                             const ReactionTable<N_sp, N_react>& tbl)
{
    if (!std::isfinite(params.T_rad) || params.T_rad <= 0.0) {
        throw std::invalid_argument(
            "chem_full_step: params.T_rad must be set to a finite positive value [K]");
    }
    if constexpr (N_react == metal_grain::N_react) {
        if (!std::isfinite(params.T_gr_K) || params.T_gr_K <= 0.0) {
            throw std::invalid_argument(
                "chem_full_step: params.T_gr_K must be set to a finite positive value [K] for metal_grain");
        }
    }

    params.zeta = shield.zeta;

    ChemFullRates rates;

    // ── Line cooling ──────────────────────────────────────────────────────────
    if constexpr (N_react == metal_grain::N_react) {
        // Metal: H2/HD/Lya primitive part from zero-metal line_cool
        double xLmbd_line_prim;
        LineCoolRates mcool{};
        line_cool<N_sp>(cell.state.y, shield.xNc_H2, shield.xNc_HD,
                        cell.state.xnH, cell.state.T_K, params.T_rad, shield.tau_cnt,
                        xLmbd_line_prim, mcool.H2, mcool.HD, mcool.Lya);
        // Molecular + atomic cooling (CO, OH, H2O, CII, CI, OI); updates mcool
        const auto& y = cell.state.y;
        line_cool_metal(cell.state.xnH, cell.state.T_K, params.T_rad,
                        y[0], y[1], y[2], y[3], y[7],
                        y[16], y[22], y[29], y[31], y[32], y[33],
                        shield.tau_cnt,
                        shield.xNc_CO, shield.xNc_OH, shield.xNc_H2O,
                        shield.xNc_CII, shield.xNc_CI, shield.xNc_OI,
                        cell.es, mcool);
        rates.xLmbd_H2  = mcool.H2;
        rates.xLmbd_HD  = mcool.HD;
        rates.xLmbd_Lya = mcool.Lya;
        rates.xLmbd_CO  = mcool.CO;
        rates.xLmbd_OH  = mcool.OH;
        rates.xLmbd_H2O = mcool.H2O;
        rates.xLmbd_CII = mcool.CII;
        rates.xLmbd_CI  = mcool.CI;
        rates.xLmbd_OI  = mcool.OI;
        rates.xLmbd_line = mcool.total();
    } else {
        // Zero-metal: H2 / HD / Ly-alpha only
        double xLmbd_line;
        line_cool<N_sp>(cell.state.y, shield.xNc_H2, shield.xNc_HD,
                        cell.state.xnH, cell.state.T_K, params.T_rad, shield.tau_cnt,
                        xLmbd_line, rates.xLmbd_H2, rates.xLmbd_HD, rates.xLmbd_Lya);
        rates.xLmbd_line = xLmbd_line;
    }

    // ── Continuum cooling ─────────────────────────────────────────────────────
    if constexpr (N_react == metal_grain::N_react) {
        cnt_cool_metal(cell.state.xnH, cell.state.T_K, params.T_rad,
                       shield.tau_cnt, shield.esc_cnt, params.Z_metal,
                       params.T_gr_K,
                       rates.xk_gr, rates.xk_gas,
                       rates.xLmbd_gr, rates.xLmbd_gas, rates.xLmbd_cnt);
    } else {
        cnt_cool(cell.state.xnH, cell.state.T_K, params.T_rad, shield.esc_cnt,
                 rates.xk_gas, rates.xLmbd_gas, rates.xLmbd_cnt);
    }

    // ── Chemistry ─────────────────────────────────────────────────────────────
    double xLmbd_ch = 0.0;
    chemcool<N_sp, N_react>(cell.state.xnH, cell.state.T_K,
                             cell.state.y, dt,
                             cell.state.xmu, cell.state.gamma,
                             xLmbd_ch, cell.var, tbl, params);
    rates.xLmbd_ch = xLmbd_ch;

    // ── Lyman-Werner photodissociation (operator-split) ───────────────────────
    // Applied as first-order exponential decay after the chemistry solver,
    // using the LW intensity J_LW21 stored in shield (units of J_21 =
    // 10^-21 erg/s/cm^2/Hz/sr).  Reactions:
    //   H2  + hν(LW) → H  + H   (Abel et al. 1997; self-shielding WG2019)
    //   HD  + hν(LW) → H  + D   (Wolcott-Green & Haiman 2011)
    //   H-  + hν     → H  + e-  (Tegmark et al. 1997; no self-shielding)
    // Species indices (same for both zero_metal and metal_grain networks):
    //   y[0]=H  y[1]=H2  y[2]=e-  y[6]=H-  y[11]=D  y[12]=HD
    if (shield.J_LW21 > 0.0) {
        auto& y          = cell.state.y;
        const double T_K = cell.state.T_K;
        const double nH  = cell.state.xnH;

        // Thermal Doppler b-parameters [units of 10^5 cm/s]
        const double b5_H2 = std::sqrt(phys::xk_B * T_K / phys::xm_p) / 1.0e5; // sqrt(2kT/2mp) = sqrt(kT/mp)
        const double b5_HD = std::sqrt(phys::xk_B * T_K / (1.5 * phys::xm_p)) / 1.0e5; // sqrt(2kT/3mp)

        // H2 self-shielding: Wolcott-Green & Haiman (2019) MNRAS 484, 2467, eq. 7–8
        // The density/temperature-dependent exponent alpha(n,T) captures the
        // weakening of self-shielding at high n and T, where excited H2
        // rovibrational states become significantly populated (LTE approach).
        //   fsh = 0.965/(1+x/b5)^alpha + 0.035/sqrt(1+x)*exp(-8.5e-4*sqrt(1+x))
        //   x = N_H2 / 5e14 cm^-2
        //   alpha(n,T) = A1(T)*10^(-c1*log10(n)) + A2(T)
        //   A1(T) = c2*log10(T/K) - c3,  A2(T) = -c4*log10(T/K) + c5
        //   c1=0.2856, c2=0.8711, c3=1.928, c4=0.9639, c5=3.892 (best-fit, eq. 8)
        // Valid range: n <= 1e7 cm^-3, T <= 8000 K, N_H2 <= 1e17 cm^-2.
        // Reduces to WG2011 (alpha≈2) in the cold-gas limit (T~100 K, high n).
        const double log10_nH  = std::log10(nH  > 0.0 ? nH  : 1.0);
        const double log10_TK  = std::log10(T_K > 0.0 ? T_K : 1.0);
        const double A1_H2     = 0.8711 * log10_TK - 1.928;
        const double A2_H2     = -0.9639 * log10_TK + 3.892;
        const double alpha_H2  = A1_H2 * std::pow(10.0, -0.2856 * log10_nH) + A2_H2;
        const double x_H2      = shield.xNc_H2 / 5.0e14;
        const double fsh_H2    = 0.965 / std::pow(1.0 + x_H2 / b5_H2, alpha_H2)
                                + 0.035 / std::sqrt(1.0 + x_H2)
                                * std::exp(-8.5e-4 * std::sqrt(1.0 + x_H2));
        const double k_H2_LW   = 1.38e-12 * shield.J_LW21 * fsh_H2;  // [s^-1]

        // HD self-shielding: WG2011 functional form (alpha=2 fixed).
        // WG2019 calibrated for H2 only; HD uses the same fitting function
        // as an approximation.
        const double x_HD  = shield.xNc_HD / 5.0e14;
        const double fsh_HD = 0.965 / std::pow(1.0 + x_HD / b5_HD, 2.0)
                            + 0.035 / std::sqrt(1.0 + x_HD)
                            * std::exp(-8.5e-4 * std::sqrt(1.0 + x_HD));
        const double k_HD_LW = 1.38e-12 * shield.J_LW21 * fsh_HD;   // [s^-1]

        // H- photodetachment (Tegmark et al. 1997).
        // H- lacks a discrete band structure, so self-shielding does not apply.
        const double k_Hm_LW = 1.1e-10 * shield.J_LW21;              // [s^-1]

        // H2 + hν(LW) → H + H
        const double dy_H2 = y[1] * (1.0 - std::exp(-k_H2_LW * dt));
        y[1]  -= dy_H2;
        y[0]  += 2.0 * dy_H2;

        // HD + hν(LW) → H + D
        const double dy_HD = y[12] * (1.0 - std::exp(-k_HD_LW * dt));
        y[12] -= dy_HD;
        y[0]  += dy_HD;
        y[11] += dy_HD;

        // H- + hν → H + e-
        const double dy_Hm = y[6] * (1.0 - std::exp(-k_Hm_LW * dt));
        y[6]  -= dy_Hm;
        y[0]  += dy_Hm;
        y[2]  += dy_Hm;
    }

    // ── CR heating ────────────────────────────────────────────────────────────
    double rate_CR = 0.0;
    if constexpr (N_react == zero_metal::N_react) {
        for (int i = 130; i <= 135; ++i) rate_CR += cell.var[i];
    } else {
        for (int i = 543; i <= 551; ++i) rate_CR += cell.var[i];
    }
    constexpr double yHe = abundance_ref::yHe;
    rates.xGam_CR = 3.4 * phys::eV_to_erg * rate_CR
                  / ((1.0 + 4.0 * yHe) * phys::xm_p);

    // ── Net cooling ───────────────────────────────────────────────────────────
    rates.xLmbd_net = rates.xLmbd_line + rates.xLmbd_cnt
                    + rates.xLmbd_ch   - rates.xGam_CR;

    return rates;
}

// ---------------------------------------------------------------------------
// make_zero_metal_table()
//
// Constructs and returns a ZeroMetalTable loaded from data files.
// Call once at startup; the returned table is read-only and can be
// shared across threads.
//
// data_dir must contain:
//   react_prm.dat         — reaction network (140 entries)
//   mass.dat              — species masses [g]
//   react_prm_saha.dat    — Saha equilibrium entries
//   pf_H3p.dat            — partition function for H3+  (species 6, 1-based)
//   pf_HD.dat             — partition function for HD   (species 13)
//   pf_HDp.dat            — partition function for HD+  (species 15)
//   pf_LiH.dat            — partition function for LiH  (species 18)
//   pf_LiHp.dat           — partition function for LiH+ (species 21)
// ---------------------------------------------------------------------------
inline ZeroMetalTable make_zero_metal_table(const std::string& data_dir)
{
    ZeroMetalTable tbl;
    load_reaction_table(tbl, data_dir + "/react_prm.dat");
    load_mass_table    (tbl, data_dir + "/mass.dat");
    load_saha_table    (tbl, data_dir + "/react_prm_saha.dat");
    load_partition_functions(tbl, {
        {  6, data_dir + "/pf_H3p.dat"  },   // H3+  (C++ y[5])
        { 13, data_dir + "/pf_HD.dat"   },   // HD   (C++ y[12])
        { 15, data_dir + "/pf_HDp.dat"  },   // HD+  (C++ y[14])
        { 18, data_dir + "/pf_LiH.dat"  },   // LiH  (C++ y[17])
        { 21, data_dir + "/pf_LiHp.dat" }    // LiH+ (C++ y[20])
    });
    return tbl;
}

// ---------------------------------------------------------------------------
// make_metal_grain_table()
//
// Constructs and returns a MetalGrainTable loaded from two data directories.
// Call once at startup; the returned table is read-only and can be
// shared across threads.
//
// data_dir must contain:
//   mass_metal.dat        — species masses [g]  (89 species)
//   pf_H3p.dat            — H3+  (species 6,  1-based)
//   pf_HD.dat             — HD   (species 13)
//   pf_HDp.dat            — HD+  (species 15)
//
// metal_data_dir must contain:
//   react_metal_grain.dat — reaction network (1200 entries)
//   react_grain_surface.dat — grain surface reactions (up to 200 entries)
//   react_metal_saha.dat  — Saha equilibrium entries
//   pf_CH3.dat            — CH3   (species 21)
//   pf_CH4.dat            — CH4   (species 22)
//   pf_HO2.dat            — HO2   (species 36)
//   pf_CO2.dat            — CO2   (species 37)
//   pf_H2CO.dat           — H2CO  (species 38)
//   pf_H2O2.dat           — H2O2  (species 39)
//   pf_LiHp.dat           — LiH+  (species 55)
// ---------------------------------------------------------------------------
inline MetalGrainTable make_metal_grain_table(const std::string& data_dir,
                                               const std::string& metal_data_dir)
{
    MetalGrainTable tbl;
    load_reaction_table     (tbl, metal_data_dir + "/react_metal_grain.dat");
    load_grain_surface_table(tbl, metal_data_dir + "/react_grain_surface.dat");
    load_mass_table         (tbl, data_dir       + "/mass_metal.dat");
    load_saha_table         (tbl, metal_data_dir + "/react_metal_saha.dat");
    load_partition_functions(tbl, {
        {  6, data_dir       + "/pf_H3p.dat"    },   // H3+  (C++ y[5])
        { 13, data_dir       + "/pf_HD.dat"     },   // HD   (C++ y[12])
        { 15, data_dir       + "/pf_HDp.dat"    },   // HD+  (C++ y[14])
        { 21, metal_data_dir + "/pf_CH3.dat"    },   // CH3  (C++ y[20])
        { 22, metal_data_dir + "/pf_CH4.dat"    },   // CH4  (C++ y[21])
        { 36, metal_data_dir + "/pf_HO2.dat"    },   // HO2  (C++ y[35])
        { 37, metal_data_dir + "/pf_CO2.dat"    },   // CO2  (C++ y[36])
        { 38, metal_data_dir + "/pf_H2CO.dat"   },   // H2CO (C++ y[37])
        { 39, metal_data_dir + "/pf_H2O2.dat"   },   // H2O2 (C++ y[38])
        { 55, metal_data_dir + "/pf_LiHp.dat"   }    // LiH+ (C++ y[54])
    });
    return tbl;
}

} // namespace chemistry
