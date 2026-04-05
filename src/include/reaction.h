// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
#include <array>
#include <cmath>
#include <algorithm>
#include "reaction_table.h"
#include "state.h"
#include "partition_function.h"

// ---------------------------------------------------------------------------
// Chemistry reaction-rate kernel.
//
// Provides:
//   eval_opacity(T_K, rho)        — continuum opacity xk_prm [cm²/g]
//   compute_base_rates(...)       — forward + reverse rate coefficients xk[2*N_react]
//   compute_rates(...)            — dy/dt and Jacobian from xk and y
//
// Mirrors Fortran react_coef() and react_rat() in reaction_prm.f.
// All rate coefficients in CGS (cm³/s or cm⁶/s for 3-body).
// xk indexing: xk[num-1] = forward, xk[num-1+N_react] = reverse (1-based num).
// ---------------------------------------------------------------------------

namespace chemistry {
namespace detail {

// ---------------------------------------------------------------------------
// eval_opacity — bilinear interpolation of continuum opacity table xk_prm
// Ported from xk_prm.f (Fortran function xk_prm(T, rho)).
// Returns mean opacity [cm²/g].  Returns 0 outside low-T or low-density range,
// 1e5 outside high-density range (optically thick limit).
// ---------------------------------------------------------------------------
inline double eval_opacity(double T_K, double rho)
{
    // T6 grid (57 points) and R grid (19 points)
    static constexpr std::array<double,57> T6 = {
        0.75e-3,1.00e-3,1.25e-3,1.50e-3,1.75e-3,2.00e-3,
        2.25e-3,2.50e-3,2.75e-3,3.00e-3,3.25e-3,
        3.50e-3,3.75e-3,4.00e-3,4.25e-3,4.50e-3,
        4.75e-3,5.00e-3,5.25e-3,5.50e-3,5.75e-3,
        6.00e-3,6.25e-3,6.50e-3,6.75e-3,7.00e-3,
        8.00e-3,9.00e-3,1.00e-2,1.10e-2,1.20e-2,
        1.40e-2,1.60e-2,1.80e-2,2.00e-2,2.50e-2,
        3.00e-2,3.50e-2,4.00e-2,4.50e-2,5.00e-2,
        5.50e-2,6.00e-2,7.00e-2,8.00e-2,9.00e-2,
        1.00e-1,1.20e-1,1.50e-1,2.00e-1,2.50e-1,
        3.00e-1,4.00e-1,5.00e-1,6.00e-1,8.00e-1,1.00e+0};

    static constexpr std::array<double,19> rlgR = {
        -5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,
         0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};

    // rlgk[i][j] = log10 opacity at T6[i], rlgR[j]  (57 x 19)
    static constexpr std::array<std::array<double,19>,57> rlgk = {{
        {-12.43,-11.91,-11.39,-10.86,-10.36,-9.804,-9.350,-8.847,-8.346,-7.845,-7.345,-6.845,-6.344,-5.843,-5.343,-4.843,-4.343,-3.843,-3.343},
        {-11.85,-11.35,-10.85,-10.35,-9.854,-9.354,-8.854,-8.354,-7.854,-7.354,-6.854,-6.354,-5.854,-5.354,-4.854,-4.354,-3.854,-3.354,-2.854},
        {-11.46,-10.96,-10.46,-9.964,-9.464,-8.964,-8.464,-7.964,-7.464,-6.964,-6.464,-5.964,-5.464,-4.964,-4.464,-3.964,-3.464,-2.964,-2.464},
        {-11.26,-10.74,-10.22,-9.704,-9.194,-8.634,-8.180,-7.677,-7.176,-6.675,-6.175,-5.675,-5.174,-4.673,-4.173,-3.673,-3.173,-2.673,-2.173},
        {-11.86,-11.12,-10.42,-9.757,-9.139,-8.557,-8.018,-7.491,-6.979,-6.471,-5.967,-5.465,-4.962,-4.461,-3.961,-3.461,-2.961,-2.461,-1.961},
        {-13.67,-12.51,-11.39,-10.43,-9.514,-8.757,-8.039,-7.433,-6.854,-6.318,-5.792,-5.281,-4.772,-4.269,-3.767,-3.266,-2.765,-2.264,-1.764},
        {-12.58,-12.14,-11.67,-10.99,-10.27,-9.348,-8.440,-7.672,-6.923,-6.322,-5.727,-5.193,-4.661,-4.150,-3.640,-3.137,-2.634,-2.133,-1.632},
        {-10.62,-10.50,-10.32,-10.07,-9.761,-9.375,-8.827,-8.025,-7.238,-6.478,-5.780,-5.178,-4.602,-4.067,-3.542,-3.032,-2.524,-2.020,-1.518},
        {-9.597,-9.352,-9.106,-8.857,-8.606,-8.354,-8.075,-7.781,-7.270,-6.635,-5.954,-5.245,-4.617,-4.036,-3.489,-2.961,-2.446,-1.937,-1.432},
        {-8.520,-8.331,-8.108,-7.879,-7.615,-7.345,-7.093,-6.844,-6.582,-6.318,-5.827,-5.300,-4.677,-4.038,-3.472,-2.918,-2.394,-1.875,-1.367},
        {-7.730,-7.481,-7.232,-6.982,-6.732,-6.483,-6.233,-5.984,-5.741,-5.498,-5.264,-5.016,-4.573,-4.121,-3.559,-3.000,-2.473,-1.947,-1.436},
        {-7.011,-6.747,-6.498,-6.255,-6.027,-5.788,-5.516,-5.249,-5.001,-4.758,-4.531,-4.309,-4.101,-3.832,-3.395,-2.937,-2.414,-1.894,-1.383},
        {-6.359,-6.111,-5.862,-5.612,-5.362,-5.112,-4.862,-4.613,-4.364,-4.119,-3.879,-3.661,-3.474,-3.272,-3.051,-2.732,-2.291,-1.825,-1.328},
        {-5.798,-5.552,-5.305,-5.057,-4.809,-4.560,-4.312,-4.061,-3.808,-3.560,-3.314,-3.089,-2.876,-2.704,-2.559,-2.329,-2.046,-1.667,-1.228},
        {-5.288,-5.048,-4.807,-4.563,-4.317,-4.068,-3.819,-3.569,-3.319,-3.070,-2.822,-2.587,-2.357,-2.183,-2.027,-1.858,-1.684,-1.385,-1.048},
        {-4.810,-4.585,-4.359,-4.119,-3.877,-3.630,-3.382,-3.133,-2.883,-2.634,-2.385,-2.145,-1.905,-1.718,-1.536,-1.400,-1.267,-1.034,-0.792},
        {-4.344,-4.143,-3.941,-3.710,-3.479,-3.235,-2.991,-2.742,-2.494,-2.244,-1.995,-1.751,-1.509,-1.304,-1.104,-0.975,-0.843,-0.661,-0.472},
        {-3.875,-3.713,-3.541,-3.328,-3.111,-2.875,-2.636,-2.390,-2.143,-1.893,-1.645,-1.399,-1.157,-0.937,-0.732,-0.590,-0.448,-0.300,-0.133},
        {-3.410,-3.288,-3.147,-2.961,-2.764,-2.539,-2.309,-2.067,-1.824,-1.576,-1.329,-1.081,-0.839,-0.608,-0.400,-0.240,-0.091, 0.037, 0.191},
        {-2.958,-2.872,-2.757,-2.606,-2.430,-2.224,-2.005,-1.771,-1.532,-1.286,-1.039,-0.792,-0.549,-0.312,-0.101, 0.077, 0.229, 0.347, 0.496},
        {-2.529,-2.467,-2.374,-2.257,-2.103,-1.922,-1.717,-1.493,-1.259,-1.017,-0.773,-0.527,-0.283,-0.042, 0.173, 0.366, 0.516, 0.632, 0.780},
        {-2.130,-2.080,-2.004,-1.916,-1.784,-1.631,-1.441,-1.234,-1.007,-0.771,-0.529,-0.284,-0.041, 0.201, 0.420, 0.628, 0.776, 0.896, 1.042},
        {-1.766,-1.717,-1.653,-1.585,-1.473,-1.349,-1.174,-0.986,-0.768,-0.542,-0.303,-0.060, 0.183, 0.425, 0.649, 0.866, 1.014, 1.142, 1.285},
        {-1.437,-1.380,-1.324,-1.269,-1.174,-1.074,-0.916,-0.751,-0.542,-0.327,-0.091, 0.147, 0.390, 0.633, 0.859, 1.083, 1.230, 1.367, 1.507},
        {-1.149,-1.070,-1.019,-0.968,-0.887,-0.806,-0.665,-0.523,-0.324,-0.123, 0.108, 0.341, 0.583, 0.825, 1.055, 1.284, 1.433, 1.580, 1.715},
        {-0.900,-0.792,-0.737,-0.681,-0.615,-0.544,-0.424,-0.299,-0.115, 0.073, 0.298, 0.523, 0.763, 1.002, 1.235, 1.463, 1.620, 1.775, 1.907},
        {-0.224,-0.120,-0.033, 0.037, 0.098, 0.189, 0.286, 0.402, 0.538, 0.684, 0.842, 1.019, 1.213, 1.405, 1.626, 1.850, 2.077, 2.307, 2.540},
        {-0.096, 0.191, 0.434, 0.632, 0.785, 0.893, 0.987, 1.084, 1.181, 1.290, 1.410, 1.556, 1.721, 1.875, 2.080, 2.301, 2.510, 2.679, 2.780},
        {-0.107, 0.234, 0.579, 0.929, 1.196, 1.392, 1.527, 1.624, 1.708, 1.799, 1.901, 2.020, 2.156, 2.270, 2.458, 2.664, 2.861, 3.022, 3.120},
        {-0.145, 0.193, 0.572, 0.991, 1.388, 1.681, 1.885, 2.026, 2.118, 2.202, 2.294, 2.391, 2.521, 2.606, 2.778, 2.973, 3.160, 3.308, 3.386},
        {-0.180, 0.139, 0.519, 0.960, 1.415, 1.810, 2.106, 2.297, 2.437, 2.537, 2.632, 2.733, 2.848, 2.905, 3.061, 3.237, 3.422, 3.605, 3.775},
        {-0.205, 0.048, 0.391, 0.825, 1.314, 1.806, 2.241, 2.585, 2.837, 3.012, 3.139, 3.259, 3.378, 3.386, 3.527, 3.684, 3.853, 4.030, 4.211},
        {-0.229, 0.018, 0.341, 0.741, 1.205, 1.726, 2.218, 2.671, 3.041, 3.307, 3.499, 3.651, 3.786, 3.762, 3.896, 4.045, 4.197, 4.340, 4.462},
        {-0.251,-0.024, 0.296, 0.709, 1.186, 1.679, 2.175, 2.681, 3.135, 3.490, 3.755, 3.949, 4.110, 4.076, 4.204, 4.354, 4.495, 4.596, 4.626},
        {-0.258,-0.054, 0.260, 0.684, 1.160, 1.657, 2.170, 2.685, 3.182, 3.612, 3.939, 4.180, 4.375, 4.351, 4.465, 4.590, 4.727, 4.877, 5.041},
        {-0.229,-0.046, 0.255, 0.676, 1.166, 1.695, 2.231, 2.765, 3.300, 3.808, 4.238, 4.582, 4.836, 4.831, 4.958, 5.046, 5.097, 5.113, 5.096},
        {-0.168, 0.014, 0.318, 0.743, 1.245, 1.772, 2.317, 2.869, 3.440, 3.971, 4.442, 4.832, 5.138, 5.101, 5.206, 5.206, 5.206, 5.206, 5.206},
        {-0.078, 0.124, 0.432, 0.845, 1.341, 1.874, 2.420, 2.977, 3.552, 4.100, 4.589, 4.995, 5.312, 5.250, 5.426, 5.426, 5.426, 5.426, 5.426},
        {-0.064, 0.221, 0.571, 0.985, 1.461, 1.980, 2.516, 3.068, 3.631, 4.175, 4.666, 5.073, 5.384, 5.392, 5.555, 5.555, 5.555, 5.555, 5.555},
        {-0.101, 0.211, 0.606, 1.084, 1.565, 2.065, 2.583, 3.120, 3.664, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.138, 0.149, 0.541, 1.038, 1.568, 2.091, 2.604, 3.125, 3.649, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.170, 0.099, 0.469, 0.941, 1.476, 2.030, 2.570, 3.091, 3.603, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.198, 0.059, 0.409, 0.854, 1.376, 1.925, 2.489, 3.032, 3.540, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.230,-0.011, 0.311, 0.736, 1.242, 1.761, 2.329, 2.895, 3.431, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.242,-0.044, 0.260, 0.670, 1.169, 1.685, 2.239, 2.799, 3.354, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.240,-0.055, 0.241, 0.649, 1.140, 1.655, 2.202, 2.758, 3.316, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.240,-0.050, 0.243, 0.639, 1.128, 1.652, 2.196, 2.749, 3.305, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.246,-0.045, 0.258, 0.661, 1.145, 1.672, 2.210, 2.754, 3.292, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.286,-0.078, 0.229, 0.634, 1.112, 1.621, 2.135, 2.646, 3.141, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.350,-0.196, 0.054, 0.400, 0.818, 1.280, 1.758, 2.236, 2.700, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.387,-0.281,-0.090, 0.184, 0.539, 0.957, 1.413, 1.881, 2.340, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.409,-0.332,-0.182, 0.041, 0.347, 0.729, 1.163, 1.620, 2.077, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.431,-0.384,-0.280,-0.118, 0.123, 0.446, 0.839, 1.275, 1.724, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.442,-0.409,-0.329,-0.200, 0.000, 0.282, 0.641, 1.055, 1.493, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.449,-0.424,-0.357,-0.249,-0.074, 0.179, 0.511, 0.905, 1.333, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.456,-0.439,-0.389,-0.305,-0.164, 0.049, 0.341, 0.704, 1.115, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
        {-0.460,-0.447,-0.405,-0.335,-0.213,-0.026, 0.238, 0.577, 0.972, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000}
    }};

    double temp6 = T_K / 1.0e6;
    double rlogR  = std::log10(rho / (temp6*temp6*temp6));

    // Range checks (mirrors Fortran if-else chain)
    bool in_T_lo  = (temp6 >= 0.75e-3) && (temp6 <= 2.5e-2);
    bool in_T_mid = (temp6 >= 2.5e-2)  && (temp6 <= 4.0e-2);
    bool in_T_hi  = (temp6 >  4.0e-2);
    bool in_R_full = (rlogR >= -5.0) && (rlogR <= 4.0);
    bool in_R_mid  = (rlogR >= -5.0) && (rlogR <= 3.0);
    bool in_R_hi   = (rlogR >= -5.0) && (rlogR <= -1.0);

    bool do_interp = (in_T_lo && in_R_full)
                  || (in_T_mid && in_R_mid)
                  || (in_T_hi  && in_R_hi);

    if (temp6 < 0.75e-3 || rlogR < -5.0) return 0.0;
    if (!do_interp) return 1.0e5;

    // Bilinear interpolation on (temp6, rlogR)
    int ms = 56;
    for (int i = 0; i < 57; ++i)
        if (temp6 <= T6[i]) { ms = i; break; }
    if (ms == 0) ms = 1;

    int ns = 18;
    for (int i = 0; i < 19; ++i)
        if (rlogR <= rlgR[i]) { ns = i; break; }
    if (ns == 0) ns = 1;

    double t = (temp6 - T6[ms-1]) / (T6[ms] - T6[ms-1]);
    double u = (rlogR - rlgR[ns-1]) / (rlgR[ns] - rlgR[ns-1]);
    double rlogk = (1.0-t)*(1.0-u)*rlgk[ms-1][ns-1]
                 +      t *(1.0-u)*rlgk[ms  ][ns-1]
                 +      t *     u *rlgk[ms  ][ns  ]
                 + (1.0-t)*     u *rlgk[ms-1][ns  ];
    return std::pow(10.0, rlogk);
}

} // namespace detail


// ---------------------------------------------------------------------------
// compute_base_rates<N_sp, N_react>
//
// Computes all rate coefficients xk[0..2*N_react-1].
//   xk[num-1]          = forward rate for reaction num (1-based)
//   xk[num-1 + N_react] = reverse rate (computed from detailed balance)
//
// Mirrors react_coef() in reaction_prm.f.
// Reactions 1-130: H/He/D/Li forward rates (hardcoded analytic fits).
// Reactions 131-139: cosmic-ray ionization rates.
// Reverse rates 1-130 computed via detailed balance (partition functions).
//
// Inputs:
//   xnH    — H number density [cm⁻³]
//   T_K    — temperature [K]
//   xmu    — mean molecular weight
//   params — ChemParams (zeta, xH, xH2, xHe)
//   tbl    — reaction table (for Cmass, delE, and pf_table)
// Output:
//   xk     — rate coefficient array [2*N_react]
// ---------------------------------------------------------------------------
template<int N_sp, int N_react>
void compute_base_rates(double xnH, double T_K, double xmu,
                        const ChemParams& params,
                        const ReactionTable<N_sp, N_react>& tbl,
                        std::array<double, 2*N_react>& xk)
{
    xk.fill(0.0);

    // Physical constants (CGS)
    constexpr double xk_B = phys::xk_B;
    constexpr double h_P  = phys::h_P;
    constexpr double pi   = phys::pi;
    constexpr double G    = phys::G;
    constexpr double xm_p = phys::xm_p;
    constexpr double yHe  = abundance_ref::yHe;

    // Temperature variables
    double T_eV    = phys::k_B_eV * T_K;
    double xlnT_eV = std::log(T_eV);
    double T300    = T_K / 300.0;
    double xlnT    = std::log(T_K);
    double xlgT    = std::log10(T_K);
    double xlgT2   = xlgT*xlgT;
    double xlgT3   = xlgT2*xlgT;
    double xlgT4   = xlgT3*xlgT;
    double xlgT5   = xlgT4*xlgT;
    double xlgT6   = xlgT5*xlgT;
    double xlT4    = std::log10(T_K / 1.0e4);

    // Critical densities for LTE interpolation (GA08 eqs.)
    double xncr_H  = std::pow(10.0,  3.0 - 0.416*xlT4 - 0.327*xlT4*xlT4);
    double xncr_H2 = std::pow(10.0,  4.845 - 1.3*xlT4 + 1.62*xlT4*xlT4);
    double xncr_He = std::pow(10.0,  5.0792*(1.0 - 1.23e-5*(T_K - 2000.0)));

    double xncr = 1.0 / (params.xH/xncr_H + params.xH2/xncr_H2 + params.xHe/xncr_He);
    double xcr  = xnH / xncr;
    double xncr_HD = 1.0e2 * xncr;
    double xcr_HD  = xnH / xncr_HD;

    // Jeans length and continuum optical depth
    double rho     = xnH * ((1.0 + 4.0*yHe) * xm_p);
    double xlmbd_J = std::sqrt(pi * xk_B * T_K / (G * xmu * xm_p * rho));
    double kap     = detail::eval_opacity(T_K, rho);
    double tau_cnt = kap * rho * xlmbd_J;

    // Partition functions
    std::array<double, 102> pf;
    eval_partition_functions(T_K, tbl, pf);

    // -----------------------------------------------------------------------
    // H, He related reactions (1-46)  [GA08 = Glover & Abel 2008]
    // -----------------------------------------------------------------------
    // 1) H + e -> H+ + 2e
    if (T_K < 1000.0)
        xk[0] = 0.0;
    else
        xk[0] = std::exp(-32.71396786
            + xlnT_eV*(13.536556
            + xlnT_eV*(-5.73932875
            + xlnT_eV*(1.56315498
            + xlnT_eV*(-0.2877056
            + xlnT_eV*(3.48255977e-2
            + xlnT_eV*(-2.63197617e-3
            + xlnT_eV*(1.11954395e-4
            + xlnT_eV*(-2.03914985e-6)))))))));

    // 2) H+ + e -> H + ph.  (case B)
    xk[1] = 2.753e-14 * std::pow(315614.0/T_K, 1.500)
            * std::pow(1.0 + std::pow(115188.0/T_K, 0.407), -2.242);

    // 3) He + e -> He+ + 2e
    if (T_K < 1000.0)
        xk[2] = 0.0;
    else
        xk[2] = std::exp(-44.09864886
            + xlnT_eV*(23.91596563
            + xlnT_eV*(-10.7532302
            + xlnT_eV*(3.05803875
            + xlnT_eV*(-0.56851189
            + xlnT_eV*(6.79539123e-2
            + xlnT_eV*(-5.00905610e-3
            + xlnT_eV*(2.06723616e-4
            + xlnT_eV*(-3.64916141e-6)))))))));

    // 4) He+ + e -> He + ph.  (optically thick: 0.68*caseA + 0.32*caseB + dielectric)
    {
        double xkrrA = 1.0e-11/std::sqrt(T_K)*(12.72 - 1.615*xlgT - 0.3162*xlgT2 + 0.0493*xlgT3);
        double xkrrB = 1.0e-11/std::sqrt(T_K)*(11.19 - 1.676*xlgT - 0.2852*xlgT2 + 0.04433*xlgT3);
        double xkdi  = 1.9e-3*std::pow(T_K,-1.5)*std::exp(-473421.0/T_K)
                      *(1.0 + 0.3*std::exp(-94684.0/T_K));
        xk[3] = 0.68*xkrrA + 0.32*xkrrB + xkdi;
    }

    // 5) He+ + e -> He++ + 2e
    if (T_K < 1000.0)
        xk[4] = 0.0;
    else
        xk[4] = std::exp(-68.71040990
            + xlnT_eV*(43.93347633
            + xlnT_eV*(-18.4806699
            + xlnT_eV*(4.70162649
            + xlnT_eV*(-0.76924663
            + xlnT_eV*(8.113042e-2
            + xlnT_eV*(-5.32402063e-3
            + xlnT_eV*(1.97570531e-4
            + xlnT_eV*(-3.16558106e-6)))))))));

    // 6) He++ + e -> He+ + ph.  (case B)
    xk[5] = 5.506e-14 * std::pow(1262456.0/T_K, 1.500)
            * std::pow(1.0 + std::pow(460752.0/T_K, 0.407), -2.242);

    // 7) H + e -> H- + ph.
    if (T_K < 6000.0)
        xk[6] = std::pow(10.0, -17.845 + 0.762*xlgT + 0.1523*xlgT2 - 0.03274*xlgT3);
    else
        xk[6] = std::pow(10.0, -16.4199 + 0.1998*xlgT2 - 5.447e-3*xlgT4 + 4.0415e-5*xlgT6);

    // 8) H- + H -> H2 + e  (Kreckel+2010)
    xk[7] = 1.35e-9 * (std::pow(T_K, 9.8493e-2) + 3.2852e-1*std::pow(T_K, 5.561e-1)
            + 2.771e-7*std::pow(T_K, 2.1826))
            / (1.0 + 6.191e-3*std::pow(T_K, 1.0461) + 8.9712e-11*std::pow(T_K, 3.0424)
            + 3.2576e-14*std::pow(T_K, 3.7741));

    // 9) H + H+ -> H2+ + ph.  (Coppola2011, Glover2015)
    xk[8] = std::pow(10.0, -18.2 - 3.194*xlgT + 1.786*xlgT2 - 0.2072*xlgT3);

    // 10) H2+ + H -> H2 + H+
    xk[9] = 6.4e-10;

    // 11) H2 + H+ -> H2+ + H  (!! zeroed below !!)
    if (T_K < 100.0)
        xk[10] = 0.0;
    else
        xk[10] = (-3.3232183e-7 + 3.3735382e-7*xlnT - 1.4491368e-7*xlnT*xlnT
                  + 3.4172805e-8*std::pow(xlnT,3) - 4.7813720e-9*std::pow(xlnT,4)
                  + 3.9731542e-10*std::pow(xlnT,5) - 1.8171411e-11*std::pow(xlnT,6)
                  + 3.5311932e-13*std::pow(xlnT,7))
                 * std::exp(-21237.15/T_K);

    // 12) H2 + e -> 2H + e  (LTE/v0 interpolation)
    double xk12v0  = 4.49e-9*std::pow(T_K,0.11)*std::exp(-101858.0/T_K);
    double xk12LTE = 1.91e-9*std::pow(T_K,0.136)*std::exp(-53407.1/T_K);
    double xlgk12  = (xcr/(1.0+xcr))*std::log10(xk12LTE)
                   + (1.0/(1.0+xcr))*std::log10(xk12v0);
    xk[11] = std::pow(10.0, xlgk12);

    // 13) H2 + H -> 3H  (!! zeroed below; reverse computed via 19 !!)
    double xk13v0  = 6.67e-12*std::sqrt(T_K)*std::exp(-(1.0 + 63593.0/T_K));
    double xk13LTE_init = 3.52e-9*std::exp(-43900.0/T_K);
    double xlgk13  = (xcr/(1.0+xcr))*std::log10(xk13LTE_init)
                   + (1.0/(1.0+xcr))*std::log10(xk13v0);
    xk[12] = std::pow(10.0, xlgk13);

    // 14) H- + e -> H + 2e
    if (T_eV <= 0.04)  // T_K <= ~464 K
        xk[13] = 0.0;
    else
        xk[13] = std::exp(-18.01849334
            + xlnT_eV*(2.3608522
            + xlnT_eV*(-0.28274430
            + xlnT_eV*(1.62331664e-2
            + xlnT_eV*(-3.36501203e-2
            + xlnT_eV*(1.17832978e-2
            + xlnT_eV*(-1.65619470e-3
            + xlnT_eV*(1.06827520e-4
            + xlnT_eV*(-2.63128581e-6)))))))));

    // 15) H- + H+ -> 2H
    xk[14] = 2.4e-6 / std::sqrt(T_K) * (1.0 + T_K/20000.0);

    // 16) H- + H+ -> H2+ + e
    if (T_K <= 8.0e3)
        xk[15] = 6.9e-9 * std::pow(T_K, -0.35);
    else
        xk[15] = 9.6e-7 * std::pow(T_K, -0.90);

    // 17) H2+ + e -> 2H
    if (T_K < 617.0)
        xk[16] = 1.0e-8;
    else
        xk[16] = 1.32e-6 * std::pow(T_K, -0.76);

    // 18) H2+ + H- -> H2 + H
    xk[17] = 1.4e-7 * std::pow(T300, -0.5);

    // 19) 3H -> H2 + H  (Forrey2013, Glover2015)
    double xk19v0 = 6.0e-32*std::pow(T_K,-0.25) + 2.0e-31*std::pow(T_K,-0.5);
    xk[18] = xk19v0;

    // 20) 2H + H2 -> 2H2  (!! zeroed below !!)
    xk[19] = xk[18] / 8.0;

    // 21) 2H2 -> 2H + H2  (LTE/v0 interpolation)
    double xk21v0  = (5.996e-30*std::pow(T_K,4.1881)/std::pow(1.0+6.761e-6*T_K,5.6881))
                   * std::exp(-54657.4/T_K);
    double xk21LTE = 1.3e-9 * std::exp(-53300.0/T_K);
    double xlgk21  = (xcr/(1.0+xcr))*std::log10(xk21LTE)
                   + (1.0/(1.0+xcr))*std::log10(xk21v0);
    xk[20] = std::pow(10.0, xlgk21);

    // 22) H + H -> H + e + H+  (Lenzuni+1991)
    xk[21] = 1.2e-17 * std::pow(T_K, 1.2)
             * std::exp(-std::abs(tbl.delE[21]) / xk_B / T_K);

    // 23) 2H + grain -> H2  (disabled)
    xk[22] = 0.0;

    // 24) He+ + H2 -> H+ + H + He
    xk[23] = 3.70e-14 * std::exp(-35.0/T_K);

    // 25) H2+ + He -> HeH+ + H  (!! zeroed below !!)
    xk[24] = 3.0e-10 * std::exp(-6717.0/T_K);

    // 26) H2+ + H2 -> H3+ + H
    xk[25] = 2.24e-9 * std::pow(T_K/300.0, 0.042) * std::exp(-T_K/4.66e4);

    // 27) H3+ + H- -> 2H2
    xk[26] = 2.30e-7 * std::pow(T300, -0.50);

    // 28) He+ + H -> H+ + He
    xk[27] = 1.2e-15 * std::pow(T300, 0.25);

    // 29) He+ + H- -> H + He
    xk[28] = 2.32e-7 * std::pow(T300, -0.52) * std::exp(T_K/22400.0);

    // 30) He+ + H2 -> H2+ + He
    xk[29] = 7.20e-15;

    // 31) HeH+ + H -> H2+ + He  (!! zeroed below !!)
    xk[30] = 1.04e-9 * std::pow(T300, 0.13) * std::exp(-T_K/3.31e4);

    // 32) HeH+ + H2 -> H3+ + He
    xk[31] = 1.53e-9 * std::pow(T300, 0.24) * std::exp(-T_K/1.48e4);

    // 33) H2+ + H- -> H3+ + e
    xk[32] = 2.7e-10 * std::pow(T300, -0.485) * std::exp(T_K/3.12e4);

    // 34) H3+ + e -> H2 + H
    xk[33] = 2.34e-8 * std::pow(T300, -0.52);

    // 35) H3+ + e -> 3H
    xk[34] = 4.36e-8 * std::pow(T300, -0.52);

    // 36) HeH+ + e -> H + He
    xk[35] = 3.0e-8 * std::pow(T300, -0.47);

    // 37) H2 + He -> 2H + He  (LTE/v0 interpolation)
    double xk37v0  = std::pow(10.0, -27.029 + 3.801*xlgT - 29487.0/T_K);
    double xk37LTE = std::pow(10.0,  -2.729 - 1.75*xlgT  - 23474.0/T_K);
    double xlgk37  = (xcr/(1.0+xcr))*std::log10(xk37LTE)
                   + (1.0/(1.0+xcr))*std::log10(xk37v0);
    xk[36] = std::pow(10.0, xlgk37);

    // 38) H- + H -> 2H + e
    if (T_eV < 0.1)
        xk[37] = 0.0;
    else
        xk[37] = std::exp(-20.372609
            + xlnT_eV*(1.13944933
            + xlnT_eV*(-1.4210135e-1
            + xlnT_eV*(8.4644554e-3
            + xlnT_eV*(-1.4327641e-3
            + xlnT_eV*(2.0122503e-4
            + xlnT_eV*(8.6639632e-5
            + xlnT_eV*(-2.5850097e-5
            + xlnT_eV*(2.4555012e-6
            + xlnT_eV*(-8.0683825e-8))))))))));

    // 39) H- + H2+ -> 3H
    xk[38] = 1.4e-7 * std::pow(T300, -0.5);

    // 40) H2 + e -> H- + H  (!! zeroed below; reverse computed via 8 !!)
    double xk40v0 = 2.7e-8 * std::pow(T_K, -1.27) * std::exp(-43000.0/T_K);
    xk[39] = xk40v0;

    // 41) He + H+ -> He+ + H  (!! zeroed below !!)
    if (T_K < 1.0e4)
        xk[40] = 1.26e-9 * std::pow(T_K, -0.75) * std::exp(-127500.0/T_K);
    else
        xk[40] = 4.0e-37 * std::pow(T_K, 4.74);

    // 42) He + H- -> He + H + e
    xk[41] = 4.1e-17 * T_K*T_K * std::exp(-19870.0/T_K);

    // 43) 2H + He -> H2 + He  (!! zeroed below !!)
    xk[42] = 6.9e-32 * std::pow(T_K, -0.4);

    // 44) H + H2+ -> H3+ + ph.
    xk[43] = 1.5e-17 * std::pow(T300, 1.8) * std::exp(20.0/T_K);

    // 45) H2 + H+ -> H3+ + ph.
    xk[44] = 1.0e-20;

    // 46) He + H+ -> HeH+ + ph.
    xk[45] = 8.0e-20 * std::pow(T300, -0.24) * std::exp(-T_K/4.0e3);

    // Zero out disabled forward reactions (Fortran lines 466-475)
    xk[10] = 0.0;   // 11
    xk[12] = 0.0;   // 13
    xk[19] = 0.0;   // 20
    xk[30] = 0.0;   // 31
    xk[39] = 0.0;   // 40
    xk[40] = 0.0;   // 41
    xk[42] = 0.0;   // 43

    // -----------------------------------------------------------------------
    // D related reactions (47-100)
    // -----------------------------------------------------------------------
    // 47) HD + e -> H + D + e  (LTE/v0, HD critical density)
    double xk47v0  = 5.09e-9*std::pow(T_K,0.128)*std::exp(-103258.0/T_K);
    double xk47LTE = 1.04e-9*std::pow(T_K,0.218)*std::exp(-53070.7/T_K);
    double xlgk47  = (xcr_HD/(1.0+xcr_HD))*std::log10(xk47LTE)
                   + (1.0/(1.0+xcr_HD))*std::log10(xk47v0);
    xk[46] = std::pow(10.0, xlgk47);

    // 48) HD + He -> H + D + He  (same as 37 but with HD critical density)
    double xlgk48 = (xcr_HD/(1.0+xcr_HD))*std::log10(xk37LTE)
                  + (1.0/(1.0+xcr_HD))*std::log10(xk37v0);
    xk[47] = std::pow(10.0, xlgk48);

    // 49) HD + H2 -> H + D + H2  (same as 21 but HD critical density)
    double xlgk49 = (xcr_HD/(1.0+xcr_HD))*std::log10(xk21LTE)
                  + (1.0/(1.0+xcr_HD))*std::log10(xk21v0);
    xk[48] = std::pow(10.0, xlgk49);

    // 50) HD + H -> 2H + D  (same as 13 but HD critical density)
    double xlgk50 = (xcr_HD/(1.0+xcr_HD))*std::log10(xk13LTE_init)
                  + (1.0/(1.0+xcr_HD))*std::log10(xk13v0);
    xk[49] = std::pow(10.0, xlgk50);

    // 51) D+ + e -> D + ph.
    xk[50] = xk[1];

    // 52) D + H+ -> D+ + H
    if (T_K < 2.0e5)
        xk[51] = 2.0e-10*std::pow(T_K,0.402)*std::exp(-37.1/T_K)
                - 3.31e-17*std::pow(T_K,1.48);
    else
        xk[51] = 3.44e-10 * std::pow(T_K, 0.35);

    // 53) D+ + H -> D + H+  (!! zeroed below !!)
    xk[52] = 2.06e-10*std::pow(T_K,0.396)*std::exp(-33.0/T_K)
            + 2.03e-9*std::pow(T_K,-0.332);

    // 54) D + H -> HD + ph.
    if (T_K < 200.0)
        xk[53] = 1.0e-25*(2.80202 - 6.63697*xlnT + 4.75619*xlnT*xlnT
                - 1.39325*std::pow(xlnT,3) + 0.178259*std::pow(xlnT,4)
                - 0.00817097*std::pow(xlnT,5));
    else
        xk[53] = 1.0e-25*std::exp(507.207 - 370.889*xlnT + 104.854*xlnT*xlnT
                - 14.4192*std::pow(xlnT,3) + 0.971469*std::pow(xlnT,4)
                - 0.0258076*std::pow(xlnT,5));

    // 55) D + H2 -> H + HD
    if (T_K < 2000.0)
        xk[54] = std::pow(10.0, -56.4737 + 5.88886*xlgT + 7.19692*xlgT2
                + 2.25069*xlgT3 - 2.16903*xlgT4 + 0.317887*xlgT5);
    else
        xk[54] = 3.17e-10 * std::exp(-5207.0/T_K);

    // 56) HD+ + H -> H+ + HD
    xk[55] = xk[9];  // = xk(10) in Fortran = 6.4e-10

    // 57) D+ + H2 -> H+ + HD
    xk[56] = (0.417 + 0.846*xlgT - 0.137*xlgT2) * 1.0e-9;

    // 58) HD + H -> H2 + D  (!! zeroed below !!)
    if (T_K < 200.0)
        xk[57] = 5.25e-11 * std::exp(-4430.0/T_K);
    else
        xk[57] = 5.25e-11 * std::exp(-4430.0/T_K + 173900.0/(T_K*T_K));

    // 59) HD + H+ -> H2 + D+  (!! zeroed below !!)
    xk[58] = 1.1e-9 * std::exp(-488.0/T_K);

    // 60) D + H+ -> HD+ + ph.
    xk[59] = 3.9e-19 * std::pow(T300, 1.8) * std::exp(20.0/T_K);

    // 61) D+ + H -> HD+ + ph.
    xk[60] = 3.9e-19 * std::pow(T300, 1.8) * std::exp(20.0/T_K);

    // 62) HD+ + e -> H + D
    xk[61] = 7.2e-8 * std::pow(T_K, -0.5);

    // 63) D + e -> D- + ph.
    xk[62] = xk[6];

    // 64) D+ + D- -> 2D
    xk[63] = xk[14];

    // 65) H+ + D- -> D + H
    xk[64] = xk[14];

    // 66) H- + D -> H + D-
    xk[65] = 6.4e-9 * std::pow(T300, 0.41);

    // 67) D- + H -> D + H-  (!! zeroed below !!)
    xk[66] = 6.4e-9 * std::pow(T300, 0.41);

    // 68) D- + H -> HD + e
    xk[67] = 0.5 * xk[7];

    // 69) D + e -> D+ + 2e
    xk[68] = xk[0];

    // 70) He+ + D -> D+ + He
    xk[69] = 1.1e-15 * std::pow(T300, 0.25);

    // 71) He + D+ -> D + He+  (!! zeroed below !!)
    if (T_K < 10000.0)
        xk[70] = 1.85e-9 * std::pow(T_K, -0.75) * std::exp(-127500.0/T_K);
    else
        xk[70] = 5.9e-37 * std::pow(T_K, 4.74);

    // 72) H2+ + D -> HD+ + H
    xk[71] = 1.07e-9 * std::pow(T300, 0.062) * std::exp(-T_K/41400.0);

    // 73) HD+ + D -> HD + D+
    xk[72] = xk[9];  // 6.4e-10

    // 74) HD+ + H -> H2+ + D  (!! zeroed below !!)
    xk[73] = 1.0e-9 * std::exp(-154.0/T_K);

    // 75) HD + e -> H + D-  (LTE/v0, HD critical density)
    double xk75v0  = 1.35e-9 * std::pow(T_K,-1.27) * std::exp(-43000.0/T_K);
    xk[74] = xk75v0;  // initial value; special override computed in reverse loop

    // 76) HD + e -> D + H-
    xk[75] = 1.35e-9 * std::pow(T_K,-1.27) * std::exp(-43000.0/T_K);

    // 77) H+ + D- -> HD+ + e
    xk[76] = 1.1e-9 * std::pow(T300, -0.4);

    // 78) D+ + H- -> HD+ + e
    xk[77] = 1.1e-9 * std::pow(T300, -0.4);

    // 79) D- + e -> D + 2e
    xk[78] = xk[13];

    // 80) D- + H -> D + H + e
    xk[79] = xk[37];

    // 81) D- + He -> D + He + e
    xk[80] = 1.5e-17 * T_K*T_K * std::exp(-19870.0/T_K);

    // 82) D+ + H- -> D + H
    xk[81] = xk[14];

    // 83) H2+ + D- -> H2 + D
    xk[82] = 1.7e-7 * std::pow(T300, -0.5);

    // 84) H2+ + D- -> 2H + D
    xk[83] = 1.7e-7 * std::pow(T300, -0.5);

    // 85) HD+ + H- -> HD + H
    xk[84] = 1.5e-7 * std::pow(T300, -0.5);

    // 86) HD+ + H- -> D + 2H
    xk[85] = 1.5e-7 * std::pow(T300, -0.5);

    // 87) HD+ + D- -> HD + D
    xk[86] = 1.9e-7 * std::pow(T300, -0.5);

    // 88) HD+ + D- -> 2D + H
    xk[87] = 1.9e-7 * std::pow(T300, -0.5);

    // 89) He+ + D- -> He + D
    xk[88] = 3.03e-7 * std::pow(T300, -0.52) * std::exp(T_K/22400.0);

    // 90) D + H2+ -> H2 + D+
    xk[89] = xk[9];

    // 91) H2+ + D -> HD + H+
    xk[90] = 1.0e-9;

    // 92) HD+ + H -> H2 + D+
    xk[91] = 1.0e-9;

    // 93) H2 + D+ -> H2+ + D  (!! zeroed below !!)
    xk[92] = xk[10];  // currently 0 since xk[10] was zeroed above

    // 94) H2 + D+ -> HD+ + H  (!! zeroed below !!)
    xk[93] = (1.04e-9 + 9.52e-9*(T_K/10000.0)
              - 1.81e-9*(T_K/10000.0)*(T_K/10000.0))
             * std::exp(-21000.0/T_K);

    // 95) HD + H+ -> HD+ + H  (!! zeroed below !!)
    xk[94] = xk[10];  // 0

    // 96) HD + H+ -> H2+ + D  (!! zeroed below !!)
    xk[95] = 1.0e-9 * std::exp(-21600.0/T_K);

    // 97) HD + D+ -> HD+ + D  (!! zeroed below !!)
    xk[96] = xk[10];  // 0

    // 98) HD + He+ -> HD+ + He
    xk[97] = xk[29];  // xk(30) = 7.20e-15

    // 99) HD + He+ -> He + H+ + D
    xk[98] = 1.85e-14 * std::exp(-35.0/T_K);

    // 100) HD + He+ -> He + H + D+
    xk[99] = 1.85e-14 * std::exp(-35.0/T_K);

    // Zero out disabled D reactions (Fortran lines 748-761)
    xk[52] = 0.0;   // 53
    xk[57] = 0.0;   // 58
    xk[58] = 0.0;   // 59
    xk[66] = 0.0;   // 67
    xk[70] = 0.0;   // 71
    xk[73] = 0.0;   // 74
    xk[74] = 0.0;   // 75 (forward; reverse set in detail balance loop)
    xk[92] = 0.0;   // 93
    xk[93] = 0.0;   // 94
    xk[94] = 0.0;   // 95
    xk[95] = 0.0;   // 96
    xk[96] = 0.0;   // 97

    // -----------------------------------------------------------------------
    // Li related reactions (101-130)  [Bovino+2011: 101-120, Lepp+2002: 121-127]
    // -----------------------------------------------------------------------
    // 101) Li+ + e -> Li + ph.
    xk[100] = 1.036e-11 / (std::sqrt(T_K/107.7)
                          * std::pow(1.0 + std::sqrt(T_K/107.7), 0.6612)
                          * std::pow(1.0 + std::sqrt(T_K/1.177e7), 1.3388));

    // 102) Li+ + H- -> Li + H
    xk[101] = 6.3e-9 / std::sqrt(T_K) * (1.0 + T_K/14000.0);

    // 103) Li- + H+ -> Li + H
    xk[102] = 2.3e-6 / std::sqrt(T_K);

    // 104) Li + e -> Li- + ph.
    xk[103] = 6.1e-17 * std::pow(T_K, 0.58) * std::exp(-T_K/1.72e4);

    // 105) Li + H+ -> Li+ + H
    xk[104] = 2.5e-40 * std::pow(T_K, 7.9) * std::exp(-T_K/1210.0);

    // 106) Li + H+ -> Li+ + H + ph.
    xk[105] = 1.7e-13 * std::pow(T_K, -0.051) * std::exp(-T_K/282000.0);

    // 107) Li + H- -> LiH + e
    xk[106] = 4.0e-10;

    // 108) Li- + H -> LiH + e
    xk[107] = 4.0e-10;

    // 109) LiH+ + H -> LiH + H+  (!! zeroed below !!)
    xk[108] = 1.0e-11 * std::exp(-67900.0/T_K);

    // 110) LiH + H+ -> LiH+ + H
    xk[109] = 1.0e-9;

    // 111) LiH + H -> Li + H2
    xk[110] = 2.0e-12 * T_K * std::exp(-T_K/1200.0);

    // 112) Li+ + H -> LiH+ + ph.
    xk[111] = 1.4e-20 * std::pow(T_K, -0.9) * std::exp(-T_K/7.0e3);

    // 113) Li + H+ -> LiH+ + ph.
    xk[112] = 5.3e-14 * std::pow(T_K, -0.49);

    // 114) LiH + H+ -> Li+ + H2
    xk[113] = 1.0e-9;

    // 115) LiH+ + e -> Li + H
    xk[114] = 3.9e-6 * std::pow(T_K, -0.70) * std::exp(-T_K/1200.0);

    // 116) LiH+ + H -> Li + H2+
    xk[115] = 9.0e-10 * std::exp(-66400.0/T_K);

    // 117) LiH+ + H -> Li+ + H2
    xk[116] = 8.7e-10 * std::pow(T_K, 0.040) * std::exp(T_K/5.92e8);

    // 118) Li + H -> LiH + ph.
    xk[117] = 4.0e-20 * std::exp(-T_K/4065.0 + std::pow(T_K/13193.0, 3.0));

    // 119) Li + H2+ -> LiH + H+
    if (T_K < 500.0)
        xk[118] = 6.3e-10 * std::exp(-2553.0/T_K);
    else
        xk[118] = 7.2e-14 * std::pow(T_K, 1.18) * std::exp(-1470.0/T_K);

    // 120) LiH + H+ -> Li + H2+  (!! zeroed below !!)
    xk[119] = 2.9e-10*std::pow(T_K,0.59) - 2.6e-10*std::pow(T_K,0.6)*std::exp(-400.0/T_K);

    // 121) Li++ + e -> Li+ + ph.
    xk[120] = 5.34e-8 * std::pow(T300, -1.23) * std::exp(-T_K/9.23e5);

    // 122) Li+++ + e -> Li++ + ph.
    xk[121] = 4.83e-11 * std::pow(T300, -0.621) * std::exp(-T_K/1.67e6);

    // 123) Li+ + D- -> Li + D
    xk[122] = 3.71e-7 * std::pow(T300, -0.51) * std::exp(T_K/4.41e4);

    // 124) Li- + D+ -> Li + D
    xk[123] = 2.28e-7 * std::pow(T300, -0.51) * std::exp(T_K/4.41e4);

    // 125) Li + e -> Li+ + 2e
    xk[124] = 3.11e-8 * std::pow(T300, 0.163) * std::exp(-6.27e4/T_K);

    // 126) Li+ + e -> Li++ + 2e
    xk[125] = 5.67e-12 * std::pow(T300, 0.715) * std::exp(-8.77e5/T_K);

    // 127) Li++ + e -> Li+++ + 2e
    xk[126] = 1.70e-12 * std::pow(T300, 0.709) * std::exp(-1.42e6/T_K);

    // 128) Li + 2H -> LiH + H
    xk[127] = 2.5e-29 / T_K;

    // 129) Li + H + H2 -> LiH + H2
    xk[128] = 4.1e-30 / T_K;

    // 130) H2 + Li -> H2 + Li+ + e  (positron-like ionisation)
    xk[129] = 9.9e-9 * std::sqrt(T_K) * std::exp(tbl.delE[129] / xk_B / T_K);

    // Zero out disabled Li reactions (Fortran lines 868-870)
    xk[108] = 0.0;   // 109
    xk[119] = 0.0;   // 120

    // -----------------------------------------------------------------------
    // Reverse rates via detailed balance (reactions 1-130)
    // -----------------------------------------------------------------------
    std::array<double, N_react> xlnKeqb{};
    xlnKeqb.fill(0.0);

    for (int ire = 0; ire < 130; ++ire) {
        int num = tbl.rean[ire];  // 1-based reaction number
        if (num < 1 || num > N_react) continue;
        int r1 = tbl.rea1[ire], r2 = tbl.rea2[ire], r3 = tbl.rea3[ire];
        int p1 = tbl.pro1[ire], p2 = tbl.pro2[ire], p3 = tbl.pro3[ire];
        int nr = tbl.nrea[ire], np = tbl.npro[ire];
        double Cm = tbl.Cmass[ire], dE = tbl.delE[ire];

        // ln(K_eq) = 1.5*(nr-np)*ln(2*pi*kB*T/h^2) + ln(Cm) + ln(pf_ratio) - dE/kB/T
        double xlnC1 = 1.5 * double(nr - np)
                      * std::log(2.0*pi*xk_B*T_K / (h_P*h_P));
        double xlnCm = std::log(Cm);

        // pf: use sentinel 1.0 for index 0 (vacant) and 101 (photon)
        auto get_pf = [&](int idx) -> double {
            if (idx == 0 || idx >= 101) return 1.0;
            return pf[idx];
        };
        double xlnCpf = std::log(get_pf(r1)) + std::log(get_pf(r2)) + std::log(get_pf(r3))
                      - std::log(get_pf(p1)) - std::log(get_pf(p2)) - std::log(get_pf(p3));

        xlnKeqb[num-1] = xlnC1 + xlnCm + xlnCpf - dE / xk_B / T_K;

        int rev_idx = num - 1 + N_react;  // 0-based index for reverse rate
        double kfwd = xk[num-1];

        if (p3 == 101) {
            // Photodissociation: reverse only if tau_cnt > 0
            if (tau_cnt <= 1.0e-10) {
                xk[rev_idx] = 0.0;
            } else {
                double esc_fact;
                if (tau_cnt <= 1.0e-3)
                    esc_fact = tau_cnt - tau_cnt*tau_cnt/2.0;
                else
                    esc_fact = 1.0 - std::exp(-tau_cnt);
                double lnKeqb_esc = xlnKeqb[num-1] + std::log(esc_fact);
                xk[rev_idx] = (kfwd > 0.0) ? std::exp(std::log(kfwd) + lnKeqb_esc) : 0.0;
            }
        } else {
            xk[rev_idx] = (kfwd > 0.0) ? std::exp(std::log(kfwd) + xlnKeqb[num-1]) : 0.0;
        }

        if (kfwd == 0.0) xk[rev_idx] = 0.0;
    }

    // -----------------------------------------------------------------------
    // Special LTE/v0 overrides for reactions with both v=0 and LTE limits
    // -----------------------------------------------------------------------
    // 8+N_react: H2 + e -> H- + H  (reverse of 8)
    {
        double xlnk40LTE = std::log(xk[7]) + xlnKeqb[7];  // ln(k8) + lnKeqb(8)
        double xk40LTE   = std::exp(xlnk40LTE);
        double xlgk40 = (xcr/(1.0+xcr))*std::log10(xk40LTE)
                      + (1.0/(1.0+xcr))*std::log10(xk40v0);
        xk[7 + N_react] = std::pow(10.0, xlgk40);
    }

    // 12+N_react: 2H + e -> H2 + e  (reverse of 12)
    {
        double xlnk12rev = std::log(xk12LTE) + xlnKeqb[11];
        xk[11 + N_react] = std::exp(xlnk12rev);
    }

    // 19+N_react: H2 + H -> 3H  (reverse of 19, uses LTE k13 derived from k19 Keqb)
    // xk13LTE declared at outer scope so reaction 50 reverse can reuse it (mirrors
    // Fortran where xk13LTE is overwritten at line 969 and reused at line 1006).
    double xk13LTE;
    {
        double xlnk13LTE = std::log(xk[18]) + xlnKeqb[18];  // ln(k19) + lnKeqb(19)
        xk13LTE          = std::exp(xlnk13LTE);
        double xlgk13 = (xcr/(1.0+xcr))*std::log10(xk13LTE)
                      + (1.0/(1.0+xcr))*std::log10(xk13v0);
        xk[18 + N_react] = std::pow(10.0, xlgk13);
    }

    // 21+N_react: 2H + H2 -> 2H2  (reverse of 21)
    {
        double xlnk21rev = std::log(xk21LTE) + xlnKeqb[20];
        xk[20 + N_react] = std::exp(xlnk21rev);
    }

    // 37+N_react: 2H + He -> H2 + He  (reverse of 37)
    {
        double xlnk37rev = std::log(xk37LTE) + xlnKeqb[36];
        xk[36 + N_react] = std::exp(xlnk37rev);
    }

    // 47+N_react: H + D + e -> HD + e  (reverse of 47)
    {
        double xlnk47rev = std::log(xk47LTE) + xlnKeqb[46];
        xk[46 + N_react] = std::exp(xlnk47rev);
    }

    // 48+N_react: H + D + He -> HD + He  (reverse of 48, uses k37LTE)
    {
        double xlnk48rev = std::log(xk37LTE) + xlnKeqb[47];
        xk[47 + N_react] = std::exp(xlnk48rev);
    }

    // 49+N_react: H + D + H2 -> HD + H2  (reverse of 49, uses k21LTE)
    {
        double xlnk49rev = std::log(xk21LTE) + xlnKeqb[48];
        xk[48 + N_react] = std::exp(xlnk49rev);
    }

    // 50+N_react: 2H + D -> HD + H  (reverse of 50, uses corrected xk13LTE)
    {
        double xlnk50rev = std::log(xk13LTE) + xlnKeqb[49];
        xk[49 + N_react] = std::exp(xlnk50rev);
    }

    // 68+N_react: HD + e -> D- + H  (reverse of 68, LTE/v0 interp, HD xcr)
    {
        double xlnk75LTE = std::log(xk[67]) + xlnKeqb[67];  // ln(k68) + lnKeqb(68)
        double xk75LTE   = std::exp(xlnk75LTE);
        double xlgk75 = (xcr_HD/(1.0+xcr_HD))*std::log10(xk75LTE)
                      + (1.0/(1.0+xcr_HD))*std::log10(xk75v0);
        xk[67 + N_react] = std::pow(10.0, xlgk75);
    }

    // Suppress He ionization reversal at low xcr (xcr < 1)
    if (xcr < 1.0) {
        xk[0 + N_react] = 0.0;   // reverse of 1
        xk[2 + N_react] = 0.0;   // reverse of 3
    }

    // -----------------------------------------------------------------------
    // CR (cosmic-ray) reactions 131-139
    // -----------------------------------------------------------------------
    constexpr double zeta_ref = 1.36e-17;
    double zeta_ratio = params.zeta / zeta_ref;
    constexpr double omega = model::cr_photo_albedo;   // albedo for CR photo-reactions

    xk[130] = 5.98e-18 * zeta_ratio;        // 131: H + CR -> H+ + e
    xk[131] = 2.86e-19 * zeta_ratio;        // 132: H2 + CR -> H+ + H + e
    xk[132] = 1.20e-17 * zeta_ratio;        // 133: H2 + CR -> H2+ + e
    xk[133] = 1.30e-18 * zeta_ratio;        // 134: H2 + CR -> 2H
    xk[134] = 3.90e-21 * zeta_ratio;        // 135: H2 + CR -> H+ + H-
    xk[135] = 6.50e-18 * zeta_ratio;        // 136: He + CR -> He+ + e
    xk[136] = 1.3e-17 * zeta_ratio * 0.2  / (1.0 - omega);  // 137: H + CRph -> H+ + e
    xk[137] = 1.3e-17 * zeta_ratio * 250.0 / (1.0 - omega); // 138: H- + CRph -> H + e
    xk[138] = 1.3e-17 * zeta_ratio * 0.2  / (1.0 - omega);  // 139: He + CRph -> He+ + e
}


// ---------------------------------------------------------------------------
// compute_rates<N_sp, N_react>
//
// Computes dy/dt = r_f[i] and Jacobian dr_fdy[i][j] for the NR solver.
// Also records per-reaction forward/reverse rates in var[0..2*N_react-1].
//
// Mirrors react_rat() in reaction_prm.f.
// Indexing: y[0..N_sp-1] is C++ 0-based; species map to y_ext[1..N_sp]
// to match 1-based Fortran species indices stored in the reaction table.
// Sentinel: y_ext[0] = y_ext[101] = 1.0 (vacant / photon slots).
//
// Inputs:
//   xk     — rate coefficients [2*N_react]
//   xnH    — H number density [cm⁻³]
//   y      — species abundances [N_sp]  (0-based)
//   tbl    — reaction table
// Outputs:
//   r_f    — net rate of change for each species [N_sp]  (0-based)
//   dr_fdy — Jacobian dr_f[i]/dy[j], stored as [N_sp][N_sp] row-major
//   var    — per-reaction rates [2*N_react] (var[num-1]=fwd, var[num-1+N_react]=rev)
// ---------------------------------------------------------------------------
template<int N_sp, int N_react>
void compute_rates(const std::array<double, 2*N_react>& xk,
                   double xnH,
                   const std::array<double, N_sp>& y,
                   const ReactionTable<N_sp, N_react>& tbl,
                   std::array<double, N_sp>& r_f,
                   std::array<double, N_sp*N_sp>& dr_fdy,
                   std::array<double, 2*N_react>& var)
{
    r_f.fill(0.0);
    dr_fdy.fill(0.0);
    var.fill(0.0);

    // Extended species array with sentinels at index 0 and 101
    std::array<double, 102> y_ext;
    y_ext.fill(0.0);
    y_ext[0]   = 1.0;   // vacant slot sentinel
    y_ext[101] = 1.0;   // photon sentinel
    for (int i = 0; i < N_sp; ++i) y_ext[i+1] = y[i];  // 1-based mapping

    // Accumulator arrays with extended sentinel indices
    std::array<double, 102> r_f_dum{};
    r_f_dum.fill(0.0);
    std::array<double, 102*102> dr_fdy_dum{};
    dr_fdy_dum.fill(0.0);

    auto dJ = [&](int i, int j) -> double& {
        return dr_fdy_dum[i * 102 + j];
    };

    // -----------------------------------------------------------------------
    // Reactions 1-130: forward + reverse, all species
    // -----------------------------------------------------------------------
    for (int ire = 0; ire < 130; ++ire) {
        int num = tbl.rean[ire];
        if (num < 1 || num > N_react) continue;
        int r1 = tbl.rea1[ire], r2 = tbl.rea2[ire], r3 = tbl.rea3[ire];
        int p1 = tbl.pro1[ire], p2 = tbl.pro2[ire], p3 = tbl.pro3[ire];
        int nr = tbl.nrea[ire], np = tbl.npro[ire];

        double kf = xk[num-1];
        double kr = xk[num-1 + N_react];

        double rate_fwd = kf * y_ext[r1] * y_ext[r2] * y_ext[r3]
                          * std::pow(xnH, double(nr-1));
        double rate_rev = kr * y_ext[p1] * y_ext[p2] * y_ext[p3]
                          * std::pow(xnH, double(np-1));
        double rate = rate_fwd - rate_rev;

        // Accumulate r_f (net rate contribution for each participating index)
        r_f_dum[r1] -= rate;
        r_f_dum[r2] -= rate;
        r_f_dum[r3] -= rate;
        r_f_dum[p1] += rate;
        r_f_dum[p2] += rate;
        r_f_dum[p3] += rate;

        var[num-1]          = rate_fwd;
        var[num-1 + N_react] = rate_rev;

        // Jacobian contributions: d(rate_fwd)/dy[r1], etc.
        double nH_fwd = std::pow(xnH, double(nr-1));
        double nH_rev = std::pow(xnH, double(np-1));

        // d/dr1
        double dr1 = kf * y_ext[r2] * y_ext[r3] * nH_fwd;
        dJ(r1,r1) -= dr1; dJ(r2,r1) -= dr1; dJ(r3,r1) -= dr1;
        dJ(p1,r1) += dr1; dJ(p2,r1) += dr1; dJ(p3,r1) += dr1;

        // d/dr2
        double dr2 = kf * y_ext[r1] * y_ext[r3] * nH_fwd;
        dJ(r1,r2) -= dr2; dJ(r2,r2) -= dr2; dJ(r3,r2) -= dr2;
        dJ(p1,r2) += dr2; dJ(p2,r2) += dr2; dJ(p3,r2) += dr2;

        // d/dr3
        double dr3 = kf * y_ext[r1] * y_ext[r2] * nH_fwd;
        dJ(r1,r3) -= dr3; dJ(r2,r3) -= dr3; dJ(r3,r3) -= dr3;
        dJ(p1,r3) += dr3; dJ(p2,r3) += dr3; dJ(p3,r3) += dr3;

        // d/dp1
        double dp1 = kr * y_ext[p2] * y_ext[p3] * nH_rev;
        dJ(r1,p1) += dp1; dJ(r2,p1) += dp1; dJ(r3,p1) += dp1;
        dJ(p1,p1) -= dp1; dJ(p2,p1) -= dp1; dJ(p3,p1) -= dp1;

        // d/dp2
        double dp2 = kr * y_ext[p1] * y_ext[p3] * nH_rev;
        dJ(r1,p2) += dp2; dJ(r2,p2) += dp2; dJ(r3,p2) += dp2;
        dJ(p1,p2) -= dp2; dJ(p2,p2) -= dp2; dJ(p3,p2) -= dp2;

        // d/dp3
        double dp3 = kr * y_ext[p1] * y_ext[p2] * nH_rev;
        dJ(r1,p3) += dp3; dJ(r2,p3) += dp3; dJ(r3,p3) += dp3;
        dJ(p1,p3) -= dp3; dJ(p2,p3) -= dp3; dJ(p3,p3) -= dp3;
    }

    // -----------------------------------------------------------------------
    // CR reactions 131-139: first-order (rate = xk * y[r2], no reverse)
    // -----------------------------------------------------------------------
    for (int ire = 130; ire < 139; ++ire) {
        int num = tbl.rean[ire];
        if (num < 1 || num > N_react) continue;
        int r2 = tbl.rea2[ire];  // r1 = 0 (CR particle, sentinel)
        int p1 = tbl.pro1[ire], p2 = tbl.pro2[ire], p3 = tbl.pro3[ire];

        double rate = xk[num-1] * y_ext[r2];

        r_f_dum[r2] -= rate;
        r_f_dum[p1] += rate;
        r_f_dum[p2] += rate;
        r_f_dum[p3] += rate;

        var[num-1] = rate;

        // Jacobian d/dr2
        double d2 = xk[num-1];
        dJ(r2,r2) -= d2;
        dJ(p1,r2) += d2;
        dJ(p2,r2) += d2;
        dJ(p3,r2) += d2;
    }

    // -----------------------------------------------------------------------
    // Copy extended arrays [1..N_sp] back to compact 0-based output arrays
    // -----------------------------------------------------------------------
    for (int i = 0; i < N_sp; ++i)
        r_f[i] = r_f_dum[i+1];

    for (int i = 0; i < N_sp; ++i)
        for (int j = 0; j < N_sp; ++j)
            dr_fdy[i*N_sp + j] = dr_fdy_dum[(i+1)*102 + (j+1)];
}

} // namespace chemistry
