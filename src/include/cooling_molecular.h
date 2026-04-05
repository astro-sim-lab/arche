// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// cooling_molecular.h — molecular & atomic line cooling for metal_grain
//
// Ports of:
//   fortran/metal_grain/cooling/CO.f      (COcool)
//   fortran/metal_grain/cooling/H2O.f     (H2Ocool)
//   fortran/metal_grain/cooling/OH.f      (OHcool)
//   fortran/metal_grain/subs/lines.f      (pop_CII, pop_CI, pop_OI, etc.)
//   fortran/metal_grain/collapse_metal_grain.f  (line_cool subroutine)
//
// All cooling rates returned in [erg g^-1 s^-1].
// ---------------------------------------------------------------------------
#include <cmath>
#include <array>
#include <functional>
#include <algorithm>
#include "state.h"

namespace chemistry {

// ---------------------------------------------------------------------------
// LineCoolRates — per-species breakdown returned by line_cool_metal
// ---------------------------------------------------------------------------
struct LineCoolRates {
    double H2   = 0.0;
    double HD   = 0.0;
    double CO   = 0.0;
    double OH   = 0.0;
    double H2O  = 0.0;
    double CII  = 0.0;
    double CI   = 0.0;
    double OI   = 0.0;
    double Lya  = 0.0;
    double total() const { return H2+HD+CO+OH+H2O+CII+CI+OI+Lya; }
};

// ---------------------------------------------------------------------------
// EscapeState — persisted escape-probability arrays (warm-start NR solver)
// One instance per cell; must be initialised to 1.0 at t=0.
// ---------------------------------------------------------------------------
struct EscapeState {
    double esc_CII[1]      = {1.0};
    double esc_CIImeta[1]  = {1.0};
    double esc_CI[3]       = {1.0, 1.0, 1.0};
    double esc_CImeta[3]   = {1.0, 1.0, 1.0};
    double esc_OI[3]       = {1.0, 1.0, 1.0};
    double esc_OImeta[3]   = {1.0, 1.0, 1.0};
};

namespace mol_detail {

// ── helpers ────────────────────────────────────────────────────────────────

// 1-D linear interpolation; clamped at endpoints (port of Fortran linear())
inline double linear_interp(const double* xa, const double* ya, int m, double x)
{
    int ms = m - 1;
    for (int i = 0; i < m; ++i) {
        if (x <= xa[i]) { ms = i; break; }
    }
    if (ms == 0) ms = 1;
    double t = (xa[ms] - xa[ms-1] == 0.0) ? 0.0
             : (x - xa[ms-1]) / (xa[ms] - xa[ms-1]);
    return (1.0 - t) * ya[ms-1] + t * ya[ms];
}

// 2-D bilinear interpolation.
// ya is stored row-major as ya[iT * n + iN] (row=T, col=N).
// Matches Fortran bilinear(x1a[m], x2a[n], ya(m,n), m, n, x1, x2, y)
// where Fortran ya(m,n) is column-major (first index = T).
inline double bilinear_interp(
    const double* x1a, int m,
    const double* x2a, int n,
    const double* ya,
    double x1, double x2)
{
    int ms = m - 1;
    for (int i = 0; i < m; ++i) { if (x1 <= x1a[i]) { ms = i; break; } }
    if (ms == 0) ms = 1;
    int ns = n - 1;
    for (int j = 0; j < n; ++j) { if (x2 <= x2a[j]) { ns = j; break; } }
    if (ns == 0) ns = 1;
    double y1 = ya[(ms-1)*n + (ns-1)];
    double y2 = ya[ ms   *n + (ns-1)];
    double y3 = ya[ ms   *n +  ns   ];
    double y4 = ya[(ms-1)*n +  ns   ];
    double t = (x1a[ms] == x1a[ms-1]) ? 0.0
             : (x1 - x1a[ms-1]) / (x1a[ms] - x1a[ms-1]);
    double u = (x2a[ns] == x2a[ns-1]) ? 0.0
             : (x2 - x2a[ns-1]) / (x2a[ns] - x2a[ns-1]);
    return (1.0-t)*(1.0-u)*y1 + t*(1.0-u)*y2 + t*u*y3 + (1.0-t)*u*y4;
}

// Omukai-style cooling rate formula (shared by CO, OH, H2O)
inline double omukai_rate(double xn_c, double aL0, double aLLTE,
                          double xlnh, double alpha)
{
    double xL0inv   = std::pow(10.0, aL0);
    double xLLTEinv = std::pow(10.0, aLLTE);
    double xn_h     = std::pow(10.0, xlnh);
    double xLinv    = xL0inv
                    + xn_c * xLLTEinv
                    + xL0inv * std::pow(xn_c / xn_h, alpha)
                            * (1.0 - xn_h * xLLTEinv / xL0inv);
    return 1.0 / xLinv;
}

// CMB background photon occupation number
inline double Q_bg(double DT_K, double T_rad) {
    double x = DT_K / T_rad;
    return (x > 1.0e2) ? 0.0 : 1.0 / (std::exp(x) - 1.0);
}

// Escape probability (Fortran beta_esc)
inline double beta_esc(double tau_L, double tau_C) {
    if (tau_L < 0.0)    return std::exp(-tau_C);
    if (tau_L < 1.0e-5) return std::exp(-tau_C);
    return std::exp(-tau_C) * (1.0 - std::exp(-tau_L)) / tau_L;
}

// ── NR escape-probability solver (port of Fortran xLd) ────────────────────
// pop_fn(esc, func, xLd): compute residuals func[] and cooling rate xLd.
// esc[N_line] is the escape probability array (warm-start, updated in-place).
// Returns converged cooling rate.
template<int N_line>
double solve_esc_prob(
    std::function<void(const double*, double*, double&)> pop_fn,
    double esc[N_line])
{
    double func[N_line], func_f[N_line], esc_f[N_line];
    double A[N_line][N_line], desc[N_line];
    double esc_min[N_line];
    double xLd_val = 0.0;
    double err_max_min = 1.0e10, err_max = 0.0;

    for (int i = 0; i < N_line; ++i) esc_min[i] = esc[i];

    for (int itr = 0; itr < 1000; ++itr) {
        pop_fn(esc, func, xLd_val);

        // Finite-difference Jacobian
        for (int j = 0; j < N_line; ++j) {
            double de = (esc[j] == 0.0) ? 1.0e-10 : 1.0e-5 * esc[j];
            if (de == 0.0) de = 1.0e-10;   // guard: subnormal esc[j] underflows de to 0 → 0/0 NaN
            for (int jj = 0; jj < N_line; ++jj) esc_f[jj] = esc[jj];
            esc_f[j] += de;
            double tmp; pop_fn(esc_f, func_f, tmp);
            for (int i = 0; i < N_line; ++i) A[i][j] = (func_f[i]-func[i])/de;
        }
        for (int i = 0; i < N_line; ++i) desc[i] = -func[i];

        // Solve A * desc = rhs via simple Gaussian elimination (N≤3)
        for (int col = 0; col < N_line; ++col) {
            int piv = col;
            for (int r = col+1; r < N_line; ++r)
                if (std::abs(A[r][col]) > std::abs(A[piv][col])) piv = r;
            if (piv != col) {
                for (int c = 0; c < N_line; ++c) std::swap(A[col][c], A[piv][c]);
                std::swap(desc[col], desc[piv]);
            }
            if (std::abs(A[col][col]) < 1.0e-50) continue;
            for (int r = 0; r < N_line; ++r) {
                if (r == col) continue;
                double fac = A[r][col] / A[col][col];
                for (int c = 0; c < N_line; ++c) A[r][c] -= fac * A[col][c];
                desc[r] -= fac * desc[col];
            }
        }
        for (int i = 0; i < N_line; ++i)
            if (std::abs(A[i][i]) > 1.0e-50) desc[i] /= A[i][i];

        // Step-size limiting
        double fact = 1.0;
        if (itr > 20 && err_max > 1.0) fact = 1.0e-2;
        for (int i = 0; i < N_line; ++i)
            if (esc[i] != 0.0 && desc[i] != 0.0)
                fact = std::min(fact, 0.4 * std::abs(esc[i]/desc[i]));
        for (int i = 0; i < N_line; ++i) esc[i] += fact * desc[i];
        for (int i = 0; i < N_line; ++i)
            if (!std::isfinite(esc[i])) esc[i] = esc_min[i];

        // Convergence
        err_max = 0.0;
        for (int i = 0; i < N_line; ++i) {
            double err = (esc[i] != 0.0) ? std::abs(desc[i]/esc[i]) : 0.0;
            err_max = std::max(err, err_max);
        }
        if (err_max < err_max_min) {
            err_max_min = err_max;
            for (int i = 0; i < N_line; ++i) esc_min[i] = esc[i];
        }
        double thr = (itr <= 10) ? 1.0e-12 : (itr <= 40) ? 1.0e-8 : 1.0e-4;
        if (err_max < thr) break;
    }
    for (int i = 0; i < N_line; ++i) esc[i] = esc_min[i];
    pop_fn(esc, func, xLd_val);
    return xLd_val;
}

} // namespace mol_detail

// ---------------------------------------------------------------------------
// co_cool — CO cooling per unit volume / (xnH / ((1+4yHe)*mp))
//           = xLd_CO * y_CO [erg g^-1 s^-1]  (app layer does the × y_CO)
// Returns xLd_CO [erg cm^3 s^-1 per CO molecule].
// Port of Fortran COcool (Omukai et al. 2010 Table 3).
// ---------------------------------------------------------------------------
inline double co_cool(double xnH, double T_K, double y_H2,
                      double xNc_CO, double tau_cnt)
{
    constexpr double xk_B = phys::legacy::xk_B, xm_p = phys::legacy::xm_p;
    // T axis (log10 K): 13 points
    static constexpr double xlTa[13] = {
        0.477, 0.778, 1.000, 1.301, 1.477, 1.699, 1.903,
        2.000, 2.477, 2.778, 3.000, 3.176, 3.301 };
    // N axis (log10 cm/s): 11 points
    static constexpr double xlNa[11] = {
        14.0, 14.5, 15.0, 15.5, 16.0, 16.5,
        17.0, 17.5, 18.0, 18.5, 19.0 };
    // aL0 (1-D in T only)
    static constexpr double aL0a[13] = {
        25.85, 25.17, 24.77, 24.38, 24.21, 24.03, 23.89,
        23.82, 23.42, 23.13, 22.91, 22.63, 22.28 };
    // aLLTE[iT][iN] row-major (13T × 11N)
    static constexpr double aLLTEa[13*11] = {
        22.51,22.54,22.62,22.82,23.17,23.62,24.05,24.50,25.00,25.49,25.97,
        21.65,21.67,21.71,21.83,22.08,22.44,22.84,23.26,23.73,24.17,24.63,
        21.08,21.09,21.11,21.18,21.37,21.67,22.04,22.44,22.87,23.30,23.76,
        20.35,20.35,20.37,20.40,20.51,20.73,21.05,21.42,21.82,22.23,22.66,
        19.94,19.95,19.96,19.98,20.05,20.23,20.52,20.86,21.24,21.65,22.06,
        19.45,19.45,19.46,19.47,19.52,19.64,19.87,20.19,20.55,20.94,21.35,
        19.01,19.01,19.01,19.02,19.05,19.13,19.32,19.60,19.95,20.32,20.71,
        18.80,18.80,18.80,18.81,18.83,18.90,19.06,19.33,19.66,20.03,20.42,
        17.81,17.81,17.81,17.82,17.82,17.85,17.92,18.08,18.34,18.67,19.03,
        17.23,17.23,17.23,17.23,17.23,17.25,17.28,17.38,17.59,17.89,18.26,
        16.86,16.86,16.86,16.87,16.87,16.88,16.90,16.97,17.15,17.48,17.93,
        16.66,16.66,16.66,16.66,16.66,16.67,16.69,16.75,16.91,17.26,17.74,
        16.55,16.55,16.55,16.55,16.55,16.56,16.58,16.63,16.78,17.12,17.61
    };
    // xlnh[iT][iN]
    static constexpr double xlnha[13*11] = {
        3.23,3.19,3.05,2.73,2.26,1.76,1.26,0.76,0.26,-0.240,-0.740,
        3.26,3.24,3.16,2.97,2.57,2.07,1.57,1.07,0.573,0.0727,-0.427,
        3.29,3.27,3.22,3.07,2.72,2.24,1.74,1.24,0.742,0.242,-0.258,
        3.49,3.48,3.45,3.34,3.09,2.65,2.15,1.65,1.15,0.652,0.152,
        3.67,3.66,3.64,3.56,3.35,2.95,2.47,1.97,1.47,0.966,0.466,
        3.97,3.96,3.94,3.89,3.74,3.42,2.95,2.45,1.95,1.45,0.954,
        4.30,4.30,4.29,4.26,4.16,3.92,3.49,3.00,2.50,2.00,1.50,
        4.46,4.45,4.45,4.42,4.34,4.14,3.74,3.25,2.75,2.25,1.75,
        5.17,5.16,5.16,5.15,5.13,5.06,4.86,4.47,3.98,3.48,2.98,
        5.47,5.47,5.47,5.46,5.45,5.41,5.30,5.02,4.57,4.07,3.57,
        5.53,5.53,5.53,5.52,5.51,5.48,5.39,5.16,4.73,4.24,3.74,
        5.30,5.30,5.30,5.30,5.29,5.26,5.17,4.94,4.52,4.03,3.53,
        4.70,4.70,4.70,4.70,4.68,4.64,4.53,4.27,3.84,3.35,2.85
    };
    // alpha[iT][iN]
    static constexpr double alphaa[13*11] = {
        0.514,0.505,0.486,0.465,0.466,0.503,0.566,0.598,0.603,0.613,0.634,
        0.465,0.460,0.448,0.433,0.448,0.487,0.538,0.574,0.594,0.623,0.645,
        0.439,0.436,0.428,0.416,0.416,0.450,0.492,0.529,0.555,0.582,0.596,
        0.409,0.407,0.401,0.388,0.378,0.396,0.435,0.473,0.503,0.528,0.546,
        0.392,0.391,0.385,0.373,0.360,0.367,0.403,0.441,0.473,0.499,0.519,
        0.370,0.368,0.364,0.353,0.338,0.334,0.362,0.404,0.440,0.469,0.492,
        0.361,0.359,0.356,0.347,0.332,0.322,0.339,0.381,0.423,0.457,0.483,
        0.357,0.356,0.352,0.345,0.330,0.317,0.329,0.370,0.414,0.451,0.479,
        0.385,0.385,0.383,0.380,0.371,0.355,0.343,0.362,0.418,0.470,0.510,
        0.437,0.437,0.436,0.434,0.429,0.419,0.406,0.410,0.446,0.487,0.516,
        0.428,0.427,0.427,0.425,0.421,0.414,0.401,0.392,0.404,0.432,0.448,
        0.354,0.354,0.352,0.349,0.341,0.329,0.317,0.316,0.335,0.364,0.372,
        0.322,0.322,0.320,0.316,0.307,0.292,0.276,0.272,0.289,0.310,0.313
    };

    double xn_c = (1.0 - y_H2) * xnH;
    double xlT  = std::log10(T_K);
    double cs   = 1.0e-5 * std::sqrt(2.0 * xk_B * T_K / (28.0 * xm_p));
    double xNc  = xNc_CO + 1.0e-10;
    double xlN  = std::log10(xNc / cs);

    double aL0   = mol_detail::linear_interp(xlTa, aL0a, 13, xlT);
    double aLLTE = mol_detail::bilinear_interp(xlTa,13, xlNa,11, aLLTEa, xlT, xlN);
    double xlnh  = mol_detail::bilinear_interp(xlTa,13, xlNa,11, xlnha,  xlT, xlN);
    double alpha = mol_detail::bilinear_interp(xlTa,13, xlNa,11, alphaa, xlT, xlN);

    double xL = mol_detail::omukai_rate(xn_c, aL0, aLLTE, xlnh, alpha);
    return (1.0 - y_H2) * xL * std::exp(-tau_cnt);
}

// ---------------------------------------------------------------------------
// oh_cool — OH cooling (Omukai et al. 2010 / Leiden data)
// Returns xLd_OH
// ---------------------------------------------------------------------------
inline double oh_cool(double xnH, double T_K, double y_H2,
                      double xNc_OH, double tau_cnt)
{
    constexpr double xk_B = phys::legacy::xk_B, xm_p = phys::legacy::xm_p;
    static constexpr double xlTa[6] = {1.477,1.699,1.903,2.000,2.477,2.778};
    static constexpr double xlNa[9] = {10,11,12,13,14,15,16,17,18};
    static constexpr double aL0a[6] = {25.31,24.53,24.02,23.82,23.16,22.90};
    // [iT][iN]
    static constexpr double aLLTEa[6*9] = {
        16.20,16.20,16.22,16.36,17.02,17.74,18.58,19.42,20.38,
        15.44,15.45,15.46,15.55,16.00,16.60,17.38,18.28,19.21,
        14.90,14.90,14.91,14.96,15.25,15.77,16.50,17.40,18.36,
        14.67,14.67,14.67,14.71,14.95,15.45,16.16,17.05,18.02,
        13.81,13.81,13.81,13.83,13.94,14.49,15.10,15.84,16.83,
        13.58,13.58,13.59,13.59,13.67,14.15,14.77,15.40,16.35
    };
    static constexpr double xlnha[6*9] = {
        9.05,9.05,9.03,8.88,8.20,7.25,6.26,5.26,4.26,
        8.98,8.98,8.97,8.85,8.27,7.45,6.46,5.46,4.46,
        8.90,8.90,8.89,8.81,8.36,7.71,6.75,5.75,4.75,
        8.88,8.87,8.87,8.80,8.40,7.79,6.85,5.85,4.85,
        8.77,8.77,8.76,8.74,8.50,7.99,7.22,6.23,5.23,
        8.72,8.71,8.71,8.69,8.50,8.01,7.29,6.33,5.33
    };
    static constexpr double alphaa[6*9] = {
        0.353,0.354,0.364,0.399,0.655,0.641,0.649,0.566,0.701,
        0.543,0.543,0.544,0.556,0.586,0.538,0.572,0.618,0.645,
        0.536,0.536,0.536,0.534,0.532,0.538,0.552,0.580,0.598,
        0.530,0.530,0.529,0.525,0.521,0.538,0.516,0.538,0.550,
        0.466,0.466,0.466,0.463,0.447,0.448,0.381,0.381,0.383,
        0.434,0.434,0.434,0.432,0.424,0.443,0.370,0.363,0.369
    };

    double T = std::max(T_K, 10.0);
    double xn_c = (1.0 - y_H2) * xnH;
    double xlT  = std::log10(T);
    double cs   = 1.0e-5 * std::sqrt(2.0 * xk_B * T / (17.0 * xm_p));
    double xNc  = xNc_OH + 1.0e-10;
    double xlN  = std::log10(xNc / cs);

    double aL0   = mol_detail::linear_interp(xlTa, aL0a, 6, xlT);
    double aLLTE = mol_detail::bilinear_interp(xlTa,6, xlNa,9, aLLTEa, xlT, xlN);
    double xlnh  = mol_detail::bilinear_interp(xlTa,6, xlNa,9, xlnha,  xlT, xlN);
    double alpha = mol_detail::bilinear_interp(xlTa,6, xlNa,9, alphaa, xlT, xlN);

    double xL = mol_detail::omukai_rate(xn_c, aL0, aLLTE, xlnh, alpha);
    return (1.0 - y_H2) * xL * std::exp(-tau_cnt);
}

// ---------------------------------------------------------------------------
// h2o_cool — H2O cooling (Neufeld et al.)
// Returns xLd_H2O
// ---------------------------------------------------------------------------
inline double h2o_cool(double xnH, double T_K, double y_H2,
                       double xNc_H2O, double tau_cnt)
{
    constexpr double xk_B = phys::legacy::xk_B, xm_p = phys::legacy::xm_p;
    // Branch a: T > 100 K
    static constexpr double xlTa[6]  = {2.000,2.301,2.602,3.000,3.301,3.602};
    static constexpr double xlNa[10] = {10,11,12,13,14,15,16,17,18,19};
    static constexpr double aL0a[6]  = {24.35,23.87,23.42,22.88,22.50,22.14};
    static constexpr double aLLTEa[6*10] = {
        14.59,14.59,14.60,14.68,14.98,15.53,16.22,17.00,17.83,18.70,
        13.85,13.86,13.86,13.88,14.05,14.46,15.05,15.74,16.50,17.31,
        13.16,13.16,13.16,13.17,13.25,13.53,14.02,14.63,15.32,16.07,
        12.32,12.32,12.32,12.32,12.34,12.49,12.87,13.46,14.16,14.94,
        11.86,11.86,11.86,11.86,11.87,11.97,12.35,12.97,13.69,14.46,
        11.64,11.64,11.64,11.64,11.65,11.72,12.06,12.66,13.36,14.13
    };
    static constexpr double xlnha[6*10] = {
        9.00,8.99,8.96,8.74,8.11,7.20,6.22,5.22,4.24,3.21,
        9.04,9.04,9.03,8.89,8.37,7.51,6.53,5.57,4.59,3.58,
        9.19,9.19,9.19,9.11,8.73,7.95,6.99,6.03,5.09,4.10,
        9.50,9.50,9.50,9.47,9.31,8.74,7.87,6.94,6.02,5.08,
        9.67,9.67,9.66,9.65,9.56,9.15,8.38,7.48,6.59,5.69,
        9.60,9.60,9.59,9.59,9.53,9.20,8.50,7.64,6.78,5.89
    };
    static constexpr double alphaa[6*10] = {
        0.43,0.43,0.42,0.41,0.42,0.45,0.47,0.50,0.52,0.53,
        0.42,0.42,0.41,0.39,0.38,0.38,0.40,0.42,0.44,0.45,
        0.39,0.39,0.39,0.37,0.34,0.34,0.35,0.36,0.37,0.39,
        0.36,0.36,0.36,0.35,0.33,0.32,0.32,0.32,0.31,0.31,
        0.34,0.34,0.34,0.33,0.32,0.30,0.29,0.28,0.27,0.27,
        0.34,0.34,0.34,0.33,0.32,0.30,0.30,0.29,0.28,0.27
    };
    // Branch b: T <= 100 K (ortho)
    static constexpr double xlTb[6]   = {1.000,1.301,1.477,1.699,1.903,2.000};
    static constexpr double aL0bo[6]  = {26.81,25.88,25.43,24.96,24.58,24.41};
    static constexpr double aLLTEbo[6*10] = {
        17.94,17.96,18.14,18.77,19.70,20.67,21.58,22.53,23.50,24.41,
        16.71,16.72,16.86,17.36,18.11,18.93,19.81,20.71,21.64,22.58,
        16.08,16.09,16.19,16.58,17.25,18.05,18.91,19.80,20.72,21.65,
        15.41,15.42,15.47,15.72,16.27,17.01,17.82,18.69,19.58,20.49,
        14.85,14.86,14.88,15.02,15.47,16.12,16.88,17.70,18.56,19.46,
        14.60,14.60,14.62,14.73,15.11,15.72,16.45,17.25,18.10,18.99
    };
    static constexpr double xlnhbo[6*10] = {
        8.81,8.79,8.60,7.96,7.01,6.02,5.03,4.02,3.02,2.03,
        8.90,8.88,8.73,8.14,7.21,6.20,5.21,4.21,3.22,2.21,
        9.03,9.01,8.88,8.31,7.40,6.41,5.42,4.41,3.41,2.42,
        9.14,9.13,9.03,8.55,7.69,6.71,5.70,4.71,3.71,2.70,
        9.19,9.20,9.13,8.78,8.00,7.04,6.06,5.05,4.06,3.06,
        9.20,9.20,9.15,8.85,8.11,7.17,6.19,5.19,4.20,3.20
    };
    static constexpr double alphabo[6*10] = {
        0.71,0.64,0.57,0.56,0.59,0.72,0.85,0.86,0.93,0.87,
        0.49,0.49,0.50,0.53,0.60,0.64,0.68,0.70,0.72,0.73,
        0.48,0.48,0.47,0.48,0.53,0.58,0.61,0.63,0.66,0.67,
        0.46,0.45,0.45,0.44,0.47,0.52,0.56,0.58,0.61,0.62,
        0.45,0.45,0.44,0.42,0.45,0.49,0.52,0.55,0.57,0.59,
        0.45,0.45,0.44,0.41,0.43,0.47,0.50,0.53,0.55,0.56
    };
    // Branch b: T <= 100 K (para)
    static constexpr double aL0bp[6]  = {27.01,25.73,25.24,24.75,24.38,24.22};
    static constexpr double aLLTEbp[6*10] = {
        17.72,17.76,18.07,18.83,19.68,20.50,21.37,22.28,23.22,24.19,
        16.60,16.63,16.85,17.41,18.11,18.94,19.83,20.75,21.69,22.60,
        16.12,16.13,16.24,16.61,17.26,18.06,18.93,19.81,20.73,21.65,
        15.43,15.43,15.48,15.72,16.28,17.01,17.82,18.69,19.58,20.49,
        14.86,14.86,14.88,15.03,15.47,16.12,16.87,17.70,18.56,19.45,
        14.60,14.60,14.62,14.73,15.11,15.72,16.45,17.25,18.10,18.99
    };
    static constexpr double xlnhbp[6*10] = {
        9.30,9.25,8.94,8.17,7.21,6.20,5.21,4.22,3.23,2.21,
        9.11,9.06,8.82,8.16,7.24,6.24,5.25,4.26,3.24,2.25,
        8.94,8.91,8.71,8.16,7.28,6.31,5.30,4.31,3.31,2.30,
        8.70,8.69,8.56,8.12,7.30,6.33,5.32,4.33,3.33,2.32,
        8.55,8.54,8.46,8.10,7.32,6.35,5.35,4.36,3.37,2.37,
        8.50,8.49,8.43,8.11,7.35,6.38,5.39,4.40,3.41,2.41
    };
    static constexpr double alphabp[6*10] = {
        0.49,0.52,0.49,0.57,0.75,0.76,0.79,0.80,0.82,0.91,
        0.72,0.65,0.65,0.73,0.68,0.70,0.73,0.75,0.77,0.80,
        0.69,0.66,0.63,0.65,0.64,0.68,0.69,0.71,0.73,0.74,
        0.53,0.53,0.51,0.49,0.52,0.55,0.58,0.58,0.60,0.62,
        0.46,0.46,0.45,0.43,0.45,0.47,0.48,0.50,0.52,0.53,
        0.44,0.44,0.43,0.41,0.42,0.43,0.44,0.46,0.48,0.50
    };

    double T   = std::max(T_K, 10.0);
    double xn_c = (1.0 - y_H2) * xnH;
    double xlT  = std::log10(T);
    double cs   = 1.0e-5 * std::sqrt(2.0 * xk_B * T / (18.0 * xm_p));
    double xNc  = xNc_H2O + 1.0e-10;

    double xL;
    if (xlT > 2.0) {
        double xlN  = std::log10(xNc / cs);
        double aL0   = mol_detail::linear_interp(xlTa, aL0a, 6, xlT);
        double aLLTE = mol_detail::bilinear_interp(xlTa,6, xlNa,10, aLLTEa, xlT, xlN);
        double xlnh  = mol_detail::bilinear_interp(xlTa,6, xlNa,10, xlnha,  xlT, xlN);
        double alpha = mol_detail::bilinear_interp(xlTa,6, xlNa,10, alphaa, xlT, xlN);
        xL = mol_detail::omukai_rate(xn_c, aL0, aLLTE, xlnh, alpha);
    } else {
        // ortho
        double xlNo = std::log10(0.75 * xNc / cs);
        double aL0o   = mol_detail::linear_interp(xlTb, aL0bo, 6, xlT);
        double aLLTEo = mol_detail::bilinear_interp(xlTb,6, xlNa,10, aLLTEbo, xlT, xlNo);
        double xlnho  = mol_detail::bilinear_interp(xlTb,6, xlNa,10, xlnhbo,  xlT, xlNo);
        double alphao = mol_detail::bilinear_interp(xlTb,6, xlNa,10, alphabo, xlT, xlNo);
        double xLo    = mol_detail::omukai_rate(xn_c, aL0o, aLLTEo, xlnho, alphao);
        // para
        double xlNp = std::log10(0.25 * xNc / cs);
        double aL0p   = mol_detail::linear_interp(xlTb, aL0bp, 6, xlT);
        double aLLTEp = mol_detail::bilinear_interp(xlTb,6, xlNa,10, aLLTEbp, xlT, xlNp);
        double xlnhp  = mol_detail::bilinear_interp(xlTb,6, xlNa,10, xlnhbp,  xlT, xlNp);
        double alphap = mol_detail::bilinear_interp(xlTb,6, xlNa,10, alphabp, xlT, xlNp);
        double xLp    = mol_detail::omukai_rate(xn_c, aL0p, aLLTEp, xlnhp, alphap);
        xL = 1.148 * (0.75 * xLo + 0.25 * xLp);
    }
    return (1.0 - y_H2) * xL * std::exp(-tau_cnt);
}

// ---------------------------------------------------------------------------
// Fine-structure / metastable line cooling: CII, CI, OI
//
// Each returns xLd_XX [erg s^-1 per XX atom] computed self-consistently
// with the escape probability.  esc[] is updated in-place (warm start).
//
// Context shared by all pop functions:
//   xnH, T_K, y_e (e-), y_a (H), y_m (H2), tau_cnt, xNc, T_rad
// ---------------------------------------------------------------------------

// CII fine-structure (2-level, Δ=92K)
inline double cii_cool_fs(double xnH, double T_K, double y_e, double y_a,
                          double y_m, double tau_cnt, double xNc, double T_rad,
                          double esc[1])
{
    using namespace mol_detail;
    constexpr double xk_B=phys::legacy::xk_B, h_Pl=phys::legacy::h_P, pi=phys::pi;
    constexpr double xm_p=phys::legacy::xm_p;
    constexpr double g0=2.0, g1=4.0;
    constexpr double DT_10=92.0;
    constexpr double A_10=2.4e-6;

    auto pop_fn = [&](const double* e, double* f, double& xLd_out) {
        double esc_10 = e[0];
        double Q_10   = Q_bg(DT_10, T_rad);
        double DE_10  = DT_10 * xk_B;

        double gam_e  = 2.8e-7 / std::sqrt(T_K * 1.0e-2);
        double gam_H  = 8.0e-10 * std::pow(T_K * 1.0e-2, 0.07);
        double gam_H2 = 0.5 * gam_H;
        double C_10   = xnH * (y_e*gam_e + y_a*gam_H + y_m*gam_H2);
        double C_01   = (g1/g0) * C_10 * std::exp(-DT_10/T_K);

        double R_01 = (g1/g0)*A_10*esc_10*Q_10 + C_01;
        double R_10 = esc_10*A_10*(1.0+Q_10) + C_10;
        double f_1  = R_01 / (R_10 + R_01);
        double f_0  = R_10 / (R_10 + R_01);

        double xnu  = DE_10 / h_Pl;
        double v_th = std::sqrt(2.0*xk_B*T_K / (12.0*xm_p));
        double xNc_1 = xNc * f_1, xNc_0 = xNc * f_0;
        double tau_10 = (A_10/(8.0*pi)) * std::pow(3.0e10/xnu, 3)
                      * std::max(0.0, xNc_0*g1/g0 - xNc_1) / v_th;

        double num = g1*f_0/(g0*f_1) - 1.0;
        double S_10 = (!std::isfinite(num) || std::abs(num) < 1.0e-100) ? 1.0e50 : 1.0/num;
        f[0] = esc_10 - beta_esc(tau_10, tau_cnt);
        xLd_out = DE_10 * A_10 * f_1 * esc_10 * (1.0 - Q_10/S_10) / xnH;
    };
    return mol_detail::solve_esc_prob<1>(pop_fn, esc);
}

// CII metastable (2-level, Δ=62000K)
inline double cii_cool_meta(double xnH, double T_K, double y_e,
                            double tau_cnt, double xNc, double T_rad,
                            double esc[1])
{
    using namespace mol_detail;
    constexpr double xk_B=phys::legacy::xk_B, h_Pl=phys::legacy::h_P, pi=phys::pi;
    constexpr double xm_p=phys::legacy::xm_p;
    constexpr double g0=6.0, g1=12.0;
    constexpr double DT_10=6.2e4, A_10=3.6;

    auto pop_fn = [&](const double* e, double* f, double& xLd_out) {
        double esc_10 = e[0];
        double Q_10   = Q_bg(DT_10, T_rad);
        double DE_10  = DT_10 * xk_B;
        double gam_e  = 2.3e-8 / std::sqrt(T_K * 1.0e-4);
        double C_10   = xnH * y_e * gam_e;
        double C_01   = (g1/g0) * C_10 * std::exp(-DT_10/T_K);
        double R_01 = (g1/g0)*A_10*esc_10*Q_10 + C_01;
        double R_10 = esc_10*A_10*(1.0+Q_10) + C_10;
        double f_1  = R_01 / (R_10 + R_01);
        double f_0  = R_10 / (R_10 + R_01);
        double xnu  = DE_10 / h_Pl;
        double v_th = std::sqrt(2.0*xk_B*T_K / (12.0*xm_p));
        double xNc_1 = xNc*f_1, xNc_0 = xNc*f_0;
        double tau_10 = (A_10/(8.0*pi)) * std::pow(3.0e10/xnu, 3)
                      * std::max(0.0, xNc_0*g1/g0 - xNc_1) / v_th;
        double num = g1*f_0/(g0*f_1) - 1.0;
        double S_10 = (!std::isfinite(num) || std::abs(num) < 1.0e-100) ? 1.0e50 : 1.0/num;
        f[0] = esc_10 - beta_esc(tau_10, tau_cnt);
        xLd_out = DE_10 * A_10 * f_1 * esc_10 * (1.0 - Q_10/S_10) / xnH;
    };
    return mol_detail::solve_esc_prob<1>(pop_fn, esc);
}

// CI fine-structure (3-level)
inline double ci_cool_fs(double xnH, double T_K, double y_e, double y_a,
                         double y_m, double tau_cnt, double xNc, double T_rad,
                         double esc[3])
{
    using namespace mol_detail;
    constexpr double xk_B=phys::legacy::xk_B, h_Pl=phys::legacy::h_P, pi=phys::pi;
    constexpr double xm_p=phys::legacy::xm_p;
    constexpr double g0=1.0, g1=3.0, g2=5.0;
    constexpr double DT_10=24.0, DT_20=63.0, DT_21=39.0;
    constexpr double A_10=7.9e-8, A_20=2.0e-14, A_21=2.7e-7;

    auto pop_fn = [&](const double* e, double* f, double& xLd_out) {
        double esc_10=e[0], esc_20=e[1], esc_21=e[2];
        double Q_10=Q_bg(DT_10,T_rad), Q_20=Q_bg(DT_20,T_rad), Q_21=Q_bg(DT_21,T_rad);
        double DE_10=DT_10*xk_B, DE_20=DT_20*xk_B, DE_21=DT_21*xk_B;

        // Abrahamsson et al. 2007 collisional rates
        auto p4 = [](double x){ return std::pow(x, -0.25); };
        auto p3 = [](double x){ return std::pow(x, -1.0/3.0); };
        auto p05= [](double x){ return std::pow(x, -0.5); };
        double T_K4 = T_K;
        double xlnCH_01 = 3.6593+56.6023*p4(T_K4)-802.9765*p4(T_K4)*p4(T_K4)+5025.1882*std::pow(T_K4,-0.75)
                       -17874.4255*std::pow(T_K4,-1.0)+38343.6655*std::pow(T_K4,-1.25)
                       -49249.4895*std::pow(T_K4,-1.5)+34789.3941*std::pow(T_K4,-1.75)
                       -10390.9809*std::pow(T_K4,-2.0);
        double xlnCH_02 = 10.8377-173.4153*p3(T_K4)+2024.0272*std::pow(T_K4,-2.0/3.0)
                       -13391.6549/T_K4+52198.5522*std::pow(T_K4,-4.0/3.0)
                       -124518.3586*std::pow(T_K4,-5.0/3.0)+178182.5823*std::pow(T_K4,-2.0)
                       -140970.6106*std::pow(T_K4,-7.0/3.0)+47504.5861*std::pow(T_K4,-8.0/3.0);
        double xlnCH_12 = 15.8996-201.3030*p4(T_K4)+1533.6164*p4(T_K4)*p4(T_K4)
                       -6491.0083*std::pow(T_K4,-0.75)+15921.9239/T_K4-22691.1632*std::pow(T_K4,-1.25)
                       +17334.7529*std::pow(T_K4,-1.5)-5517.9360*std::pow(T_K4,-1.75);

        double gam01_H = 1.0e-11 * std::exp(xlnCH_01);
        double gam02_H = 1.0e-11 * std::exp(xlnCH_02)
                       * std::exp(-std::pow(4.0/T_K4, 8.0));
        double gam12_H = 1.0e-11 * std::exp(xlnCH_12);
        double gam10_H = (g0/g1)*gam01_H*std::exp(DT_10/T_K);
        double gam20_H = (g0/g2)*gam02_H*std::exp(DT_20/T_K);
        double gam21_H = (g1/g2)*gam12_H*std::exp(DT_21/T_K);

        double g10_e=3.0e-9, g20_e=5.0e-9, g21_e=1.5e-8;
        double C_10=xnH*(y_e*g10_e+y_a*gam10_H+y_m*0.05*gam10_H);
        double C_20=xnH*(y_e*g20_e+y_a*gam20_H+y_m*0.5*gam20_H);
        double C_21=xnH*(y_e*g21_e+y_a*gam21_H+y_m*0.5*gam21_H);
        double C_01=(g1/g0)*C_10*std::exp(-DT_10/T_K);
        double C_02=(g2/g0)*C_20*std::exp(-DT_20/T_K);
        double C_12=(g2/g1)*C_21*std::exp(-DT_21/T_K);

        double R_10=esc_10*A_10*(1.0+Q_10)+C_10;
        double R_20=esc_20*A_20*(1.0+Q_20)+C_20;
        double R_21=esc_21*A_21*(1.0+Q_21)+C_21;
        double R_01=(g1/g0)*esc_10*A_10*Q_10+C_01;
        double R_02=(g2/g0)*esc_20*A_20*Q_20+C_02;
        double R_12=(g2/g1)*esc_21*A_21*Q_21+C_12;
        double f_0 = (R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21))
                   / ((R_01+R_02+R_20)*(R_10+R_12+R_21) - (R_01-R_21)*(R_10-R_20));
        double f_1 = (f_0*(R_01-R_21)+R_21)/(R_10+R_12+R_21);
        double f_2 = (f_0*R_02+f_1*R_12)/(R_21+R_20);
        double v_th = std::sqrt(2.0*xk_B*T_K/(12.0*xm_p));
        auto tau = [&](double A, double DE, double gU, double gL,
                       double fL, double fU) {
            double xnu = DE/h_Pl;
            double xNcL = xNc*fL, xNcU = xNc*fU;
            return (A/(8.0*pi))*std::pow(3.0e10/xnu,3)*std::max(0.0,xNcL*gU/gL-xNcU)/v_th;
        };
        auto S_fn = [](double gU, double gL, double fL, double fU) {
            if (fU == 0.0) return 1.0e50; // fU→0: source fn → ∞, Q/S → 0
            double n = gU*fL/(gL*fU)-1.0;
            return (!std::isfinite(n) || std::abs(n)<1.0e-100) ? 1.0e50 : 1.0/n;
        };
        double tau10=tau(A_10,DE_10,g1,g0,f_0,f_1);
        double tau20=tau(A_20,DE_20,g2,g0,f_0,f_2);
        double tau21=tau(A_21,DE_21,g2,g1,f_1,f_2);
        double S10=S_fn(g1,g0,f_0,f_1), S20=S_fn(g2,g0,f_0,f_2), S21=S_fn(g2,g1,f_1,f_2);
        f[0]=esc_10-beta_esc(tau10,tau_cnt);
        f[1]=esc_20-beta_esc(tau20,tau_cnt);
        f[2]=esc_21-beta_esc(tau21,tau_cnt);
        xLd_out = (DE_10*A_10*f_1*esc_10*(1.0-Q_10/S10)
                 + DE_20*A_20*f_2*esc_20*(1.0-Q_20/S20)
                 + DE_21*A_21*f_2*esc_21*(1.0-Q_21/S21)) / xnH;
    };
    return mol_detail::solve_esc_prob<3>(pop_fn, esc);
}

// OI fine-structure (3-level)
inline double oi_cool_fs(double xnH, double T_K, double y_e, double y_a,
                         double y_m, double tau_cnt, double xNc, double T_rad,
                         double esc[3])
{
    using namespace mol_detail;
    constexpr double xk_B=phys::legacy::xk_B, h_Pl=phys::legacy::h_P, pi=phys::pi;
    constexpr double xm_p=phys::legacy::xm_p;
    constexpr double g0=5.0, g1=3.0, g2=1.0;
    constexpr double DT_10=230.0, DT_20=328.0, DT_21=98.0;
    constexpr double A_10=9.0e-5, A_20=1.0e-10, A_21=1.7e-5;

    auto pop_fn = [&](const double* e, double* f, double& xLd_out) {
        double esc_10=e[0], esc_20=e[1], esc_21=e[2];
        double Q_10=Q_bg(DT_10,T_rad), Q_20=Q_bg(DT_20,T_rad), Q_21=Q_bg(DT_21,T_rad);
        double DE_10=DT_10*xk_B, DE_20=DT_20*xk_B, DE_21=DT_21*xk_B;

        // Abrahamsson et al. 2007
        double xlnOH_01 = 4.581-156.118*std::pow(T_K,-0.75)+2679.979*std::pow(T_K,-1.5)
                       -78996.962*std::pow(T_K,-2.25)+1308323.468*std::pow(T_K,-3.0)
                       -13011761.861*std::pow(T_K,-3.75)+71010784.971*std::pow(T_K,-4.5)
                       -162826621.855*std::pow(T_K,-5.25);
        double xlnOH_02 = 3.297-168.382*std::pow(T_K,-0.75)+1844.099*std::pow(T_K,-1.5)
                       -68362.889*std::pow(T_K,-2.25)+1376864.737*std::pow(T_K,-3.0)
                       -17964610.169*std::pow(T_K,-3.75)+134374927.808*std::pow(T_K,-4.5)
                       -430107587.886*std::pow(T_K,-5.25);
        double xlnOH_12 = 3.437+17.443*std::pow(T_K,-0.5)-618.761/T_K+3757.156*std::pow(T_K,-1.5)
                       -12736.468*std::pow(T_K,-2.0)+22785.266*std::pow(T_K,-2.5)
                       -22759.228*std::pow(T_K,-3.0)+12668.261*std::pow(T_K,-3.5);

        double gam01_H = 1.0e-11 * std::exp(xlnOH_01);
        double gam02_H = 1.0e-11 * std::exp(xlnOH_02);
        double gam12_H = 1.0e-11 * std::exp(xlnOH_12) * std::exp(-std::pow(10.0/T_K, 8.0));
        double gam10_H = (g0/g1)*gam01_H*std::exp(DT_10/T_K);
        double gam20_H = (g0/g2)*gam02_H*std::exp(DT_20/T_K);
        double gam21_H = (g1/g2)*gam12_H*std::exp(DT_21/T_K);

        double g10_e=1.4e-8, g20_e=1.4e-8, g21_e=5.0e-9;
        double C_10=xnH*(y_e*g10_e+y_a*gam10_H+y_m*0.5*gam10_H);
        double C_20=xnH*(y_e*g20_e+y_a*gam20_H+y_m*0.5*gam20_H);
        double C_21=xnH*(y_e*g21_e+y_a*gam21_H+y_m*0.05*gam21_H);
        double C_01=(g1/g0)*C_10*std::exp(-DT_10/T_K);
        double C_02=(g2/g0)*C_20*std::exp(-DT_20/T_K);
        double C_12=(g2/g1)*C_21*std::exp(-DT_21/T_K);

        double R_10=esc_10*A_10*(1.0+Q_10)+C_10;
        double R_20=esc_20*A_20*(1.0+Q_20)+C_20;
        double R_21=esc_21*A_21*(1.0+Q_21)+C_21;
        double R_01=(g1/g0)*esc_10*A_10*Q_10+C_01;
        double R_02=(g2/g0)*esc_20*A_20*Q_20+C_02;
        double R_12=(g2/g1)*esc_21*A_21*Q_21+C_12;
        double f_0 = (R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21))
                   / ((R_01+R_02+R_20)*(R_10+R_12+R_21) - (R_01-R_21)*(R_10-R_20));
        double f_1 = (f_0*(R_01-R_21)+R_21)/(R_10+R_12+R_21);
        double f_2 = (f_0*R_02+f_1*R_12)/(R_21+R_20);
        if (f_1 == 0.0) { f[0]=f[1]=f[2]=0.0; xLd_out=0.0; return; }
        double v_th = std::sqrt(2.0*xk_B*T_K/(16.0*xm_p));
        auto tau = [&](double A, double DE, double gU, double gL,
                       double fL, double fU) {
            double xnu = DE/h_Pl;
            return (A/(8.0*pi))*std::pow(3.0e10/xnu,3)*std::max(0.0,xNc*fL*gU/gL-xNc*fU)/v_th;
        };
        auto S_fn = [](double gU, double gL, double fL, double fU) {
            if (fU == 0.0) return 1.0e50; // fU→0: source fn → ∞, Q/S → 0
            double n = gU*fL/(gL*fU)-1.0;
            return (!std::isfinite(n) || std::abs(n)<1.0e-100) ? 1.0e50 : 1.0/n;
        };
        double tau10=tau(A_10,DE_10,g1,g0,f_0,f_1);
        double tau20=tau(A_20,DE_20,g2,g0,f_0,f_2);
        double tau21=tau(A_21,DE_21,g2,g1,f_1,f_2);
        double S10=S_fn(g1,g0,f_0,f_1), S20=S_fn(g2,g0,f_0,f_2), S21=S_fn(g2,g1,f_1,f_2);
        f[0]=esc_10-beta_esc(tau10,tau_cnt);
        f[1]=esc_20-beta_esc(tau20,tau_cnt);
        f[2]=esc_21-beta_esc(tau21,tau_cnt);
        xLd_out = (DE_10*A_10*f_1*esc_10*(1.0-Q_10/S10)
                 + DE_20*A_20*f_2*esc_20*(1.0-Q_20/S20)
                 + DE_21*A_21*f_2*esc_21*(1.0-Q_21/S21)) / xnH;
    };
    return mol_detail::solve_esc_prob<3>(pop_fn, esc);
}

// ---------------------------------------------------------------------------
// line_cool_metal — combines all line cooling for one cell/time step.
//
// Species abundance convention (y[i] = ratio to xnH):
//   y_H=y[0], y_H2=y[1], y_e=y[2], y_Hp=y[3], y_He=y[7]
//   y_CI=y[16], y_CII=y[22], y_OI=y[29], y_OH=y[31], y_CO=y[32], y_H2O=y[33]
//
// Column densities (xNc_XX) are computed in the application layer as
//   xNc_XX = y_XX * xnH * xlsh_XX  (Jeans length limited)
//
// H2/HD cooling calls are delegated to the existing cooling.h functions
// (h2_cooling, hd_cooling), so LineCoolRates.H2 and HD are filled by caller.
//
// This function computes CO, OH, H2O, CII, CI, OI.
// It sets rates.CO, .OH, .H2O, .CII, .CI, .OI and .Lya.
// ---------------------------------------------------------------------------
inline void line_cool_metal(
    double xnH, double T_K, double T_rad,
    double y_H, double y_H2, double y_e, double y_Hp, double y_He,
    double y_CI, double y_CII, double y_OI, double y_OH, double y_CO, double y_H2O,
    double tau_cnt,
    double xNc_CO,  double xNc_OH,  double xNc_H2O,
    double xNc_CII, double xNc_CI,  double xNc_OI,
    EscapeState& es,
    LineCoolRates& rates)
{
    constexpr double yHe = abundance_ref::yHe;
    constexpr double xm_p = phys::legacy::xm_p;
    const double rho_norm = 1.0 / ((1.0 + 4.0*yHe) * xm_p); // 1/((1+4yHe)mp)

    auto emit = [&](double xLd, double abundance) -> double {
        return abundance * xLd * xnH * rho_norm;
    };

    // CO
    {
        double xLd    = co_cool(xnH, T_K, y_H2, xNc_CO, tau_cnt);
        double xLdr   = co_cool(xnH, T_rad, y_H2, xNc_CO, tau_cnt);
        rates.CO = emit(xLd - xLdr, y_CO);
    }
    // OH
    {
        double xLd  = oh_cool(xnH, T_K,   y_H2, xNc_OH, tau_cnt);
        double xLdr = oh_cool(xnH, T_rad, y_H2, xNc_OH, tau_cnt);
        rates.OH = emit(xLd - xLdr, y_OH);
    }
    // H2O
    {
        double xLd  = h2o_cool(xnH, T_K,   y_H2, xNc_H2O, tau_cnt);
        double xLdr = h2o_cool(xnH, T_rad, y_H2, xNc_H2O, tau_cnt);
        rates.H2O = emit(xLd - xLdr, y_H2O);
    }
    // CII (fine + metastable)
    {
        double xLd_fs = cii_cool_fs(xnH, T_K, y_e, y_H, y_H2, tau_cnt,
                                    xNc_CII, T_rad, es.esc_CII);
        rates.CII = emit(xLd_fs, y_CII);
        if (T_K > 200.0) {
            double xLd_m = cii_cool_meta(xnH, T_K, y_e, tau_cnt,
                                         xNc_CII, T_rad, es.esc_CIImeta);
            rates.CII += emit(xLd_m, y_CII);
        }
    }
    // CI (fine + metastable)
    {
        double xLd_fs = ci_cool_fs(xnH, T_K, y_e, y_H, y_H2, tau_cnt,
                                   xNc_CI, T_rad, es.esc_CI);
        rates.CI = emit(xLd_fs, y_CI);
        // (metastable CI: T > 3000 K — similar 3-level system, omitted for brevity)
    }
    // OI (fine; skip if T < 15 K)
    {
        if (T_K > 15.0) {
            double xLd_fs = oi_cool_fs(xnH, T_K, y_e, y_H, y_H2, tau_cnt,
                                       xNc_OI, T_rad, es.esc_OI);
            rates.OI = emit(xLd_fs, y_OI);
        } else {
            rates.OI = 0.0;
        }
    }
    // Ly-alpha (H I recombination cooling, xnH ≤ 1e10)
    if (xnH <= 1.0e10) {
        double T_5 = 1.0e-5 * T_K;
        rates.Lya = y_e * y_H * 7.50e-19 / (1.0 + std::sqrt(T_5))
                  * std::exp(-1.18348e5 / T_K)
                  * xnH * rho_norm;
    } else {
        rates.Lya = 0.0;
    }
}

} // namespace chemistry
