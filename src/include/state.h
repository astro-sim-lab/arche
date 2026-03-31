#pragma once
#include <array>
#include <cmath>
#include <limits>
#include "species.h"

namespace chemistry {

// ---------------------------------------------------------------------------
// ChemState  — single-cell state updated by the chemistry kernel.
//
// Holds only what the kernel needs to read/write per time step.
// Thermodynamic variables that can be derived (rho, p, e) are kept
// in CollapseState (app layer) to avoid redundancy.
// ---------------------------------------------------------------------------
template<int N_sp>
struct ChemState {
    std::array<double, N_sp> y{};  // species abundances (dimensionless, relative to xnH)
    double xnH   = 0.0;            // H number density [cm^-3]
    double T_K   = 0.0;            // gas temperature [K]
    double xmu   = 1.0;            // mean molecular weight [m_p]
    double gamma = 5.0/3.0;        // adiabatic index
};

// Convenience aliases
using ChemStateZM = ChemState<zero_metal::N_sp>;
using ChemStateMG = ChemState<metal_grain::N_sp>;

// ---------------------------------------------------------------------------
// ChemRates  — scalar rates returned by the kernel to the app layer.
//
// The full var[2*N_react] reaction rate array is NOT exposed.
// CR heating (var[131:136] in Fortran) is summed inside the kernel.
// ---------------------------------------------------------------------------
struct ChemRates {
    double xLmbdch = 0.0;  // net chemistry cooling rate [erg g^-1 s^-1]
    double xGam_CR = 0.0;  // CR heating rate            [erg g^-1 s^-1]
};

// ---------------------------------------------------------------------------
// ChemShielding  — shielding environment pre-computed by the calling code.
//
// For a 1-D collapse app: derived from local Jeans / shielding lengths.
// For a 3-D fluid code:   from RT, TreeCol, or Sobolev approximation.
//
// zeta already carries the full attenuation factor; the chemistry kernel
// uses this value directly without further reduction.
//
// Metal-only column densities (xNc_CO, xNc_OH, ...) are zero-initialised
// and ignored by the zero_metal network.
// ---------------------------------------------------------------------------
struct ChemShielding {
    double zeta;           // effective CR ionization rate [s^-1]
    double xNc_H   = 0.0; // H  column density [cm^-2]
    double xNc_H2  = 0.0; // H2 column density [cm^-2]
    double xNc_HD  = 0.0; // HD column density [cm^-2]
    double tau_cnt = 0.0; // continuum optical depth
    double esc_cnt = 1.0; // continuum escape fraction  (1 = optically thin)
    double J_LW21  = 0.0; // Lyman-Werner intensity [10^-21 erg/s/cm^2/Hz/sr]
                           // Photodissociates H2/HD and photodetaches H-.
                           // 0.0 = no LW field (default).
    // metal_grain only (leave at 0 for zero_metal network):
    double xNc_CO  = 0.0;
    double xNc_OH  = 0.0;
    double xNc_H2O = 0.0;
    double xNc_CII = 0.0;
    double xNc_CI  = 0.0;
    double xNc_OI  = 0.0;
};

// ---------------------------------------------------------------------------
// ChemFullRates  — complete cooling/heating breakdown returned by
//                  chem_full_step().
//
// All rates in [erg g^-1 s^-1].
// Metal-only fields are 0 for the zero_metal network.
// ---------------------------------------------------------------------------
struct ChemFullRates {
    // Aggregate rates
    double xLmbd_net  = 0.0;  // net cooling  = line + cnt + ch − CR
    double xLmbd_line = 0.0;  // total line cooling
    double xLmbd_cnt  = 0.0;  // total continuum cooling  (grain + gas)
    double xLmbd_ch   = 0.0;  // chemistry cooling
    double xGam_CR    = 0.0;  // CR ionization heating
    // Per-line rates
    double xLmbd_H2   = 0.0;
    double xLmbd_HD   = 0.0;
    double xLmbd_Lya  = 0.0;
    // metal_grain only (0 for zero_metal):
    double xLmbd_CO   = 0.0;
    double xLmbd_OH   = 0.0;
    double xLmbd_H2O  = 0.0;
    double xLmbd_CII  = 0.0;
    double xLmbd_CI   = 0.0;
    double xLmbd_OI   = 0.0;
    double xLmbd_gr   = 0.0;  // grain continuum cooling
    double xLmbd_gas  = 0.0;  // gas (ff+CIA) continuum cooling
    // Opacities — needed to compute tau_cnt in the next step
    double xk_gas     = 0.0;  // gas opacity [cm^2/g]
    double xk_gr      = 0.0;  // grain opacity × Z_metal [cm^2/g]  (metal only)
};

// ---------------------------------------------------------------------------
// ChemParams  — external parameters that are read-only for the kernel.
//
// Contains both zero_metal and metal_grain parameters.
// For zero_metal runs: T_gr_K = 0, Z_metal = 0 (unused).
// ---------------------------------------------------------------------------
struct ChemParams {
    double zeta    = 0.0;   // cosmic-ray ionization rate [s^-1]
    // NaN by default so missing caller assignment is caught explicitly.
    double T_rad   = std::numeric_limits<double>::quiet_NaN(); // CMB radiation temperature [K]
    double T_gr_K  = 0.0;   // grain temperature [K]  (metal_grain only)
    double Z_metal = 0.0;   // metallicity [Z_sun]    (metal_grain only)
    double T_cr_desorp = 70.0; // effective CR desorption spike temperature [K]

    // critical densities for LTE/low-density interpolation (filled by kernel)
    double xH  = 0.0;  // H abundance (copy of y[H],  used in xcrit COMMON)
    double xH2 = 0.0;  // H2 abundance (2*y[H2])
    double xHe = 0.0;  // He abundance (copy of y[He])

    // metal_grain extra: UV radiation field strengths (xJH2, xJH2O, xJtot)
    double xJH2  = 0.0;
    double xJH2O = 0.0;
    double xJtot = 0.0;
};

// ---------------------------------------------------------------------------
// Physical constants  (CGS)
// ---------------------------------------------------------------------------
namespace phys {
    constexpr double xk_B = 1.380662e-16;   // Boltzmann [erg/K]
    constexpr double h_P  = 6.626176e-27;   // Planck [erg·s]
    constexpr double hbar = 1.0545919e-27;  // reduced Planck [erg·s]
    constexpr double xm_p = 1.67262e-24;    // proton mass [g]
    constexpr double xm_e = 9.19941e-28;    // electron mass [g]
    constexpr double qe   = 4.80653e-10;    // elementary charge [statC]
    constexpr double pi   = 3.14159265358979;
    constexpr double G    = 6.6720e-8;      // gravitational [cgs]
    constexpr double sigma_B = 5.67e-5;     // Stefan-Boltzmann [erg/cm^2/s/K^4]
    constexpr double eV_to_erg = 1.6022e-12;
    constexpr double k_B_eV = 8.61735e-5;   // Boltzmann [eV/K]

    // Legacy precision constants used in Fortran-port fitting formulas.
    // Keep this set to preserve historical numerics where those reduced
    // significant digits were used intentionally.
    namespace legacy {
        constexpr double xk_B = 1.38e-16;
        constexpr double h_P  = 6.63e-27;
        constexpr double xm_p = 1.67e-24;
    } // namespace legacy
} // namespace phys

namespace numerics {
    constexpr double eps_gaussj = 1.0e-13;
    constexpr double eps_it_prim = 1.0e-13;
    constexpr double eps_it_metal = 1.0e-5;
    constexpr double eps_y = 1.0e-10;
    constexpr double xnH_eq = 1.0e18;

} // namespace numerics

namespace model {
    constexpr double cr_photo_albedo = 0.6;
} // namespace model

} // namespace chemistry
