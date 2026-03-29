#pragma once
// solver.h — implicit chemistry ODE integrator (zero-metal network)
//
// Ported from:
//   fortran/zero_metal/subs/chemistry_prm.f  — chemcool (secant), chemreact (NR)
//   fortran/zero_metal/subs/saha_prm.f       — equichem (Saha equilibrium)
//   fortran/zero_metal/subs/math.f           — gaussj
//
// Optional Eigen support: define CHEMISTRY_USE_EIGEN before including, or
// add -DCHEMISTRY_USE_EIGEN to the compiler flags.

#include "cooling.h"             // c_H2
#include "partition_function.h"  // eval_partition_functions
#include "reaction.h"            // compute_base_rates, compute_rates
#include "reaction_table.h"
#include <array>
#include <cmath>
#include <algorithm>
#include <cstring>

#ifdef CHEMISTRY_USE_EIGEN
#  include <Eigen/Dense>
#endif

namespace chemistry {

// ─────────────────────────────────────────────────────────────────────────────
// gaussj_solve<N> — solve A·x = b in-place; b receives the solution.
//   a[N×N] row-major (destroyed), b[N].
//   Port of Numerical Recipes gaussj (math.f) adapted to C++ row-major.
// ─────────────────────────────────────────────────────────────────────────────
template<int N>
void gaussj_solve(std::array<double, N*N>& a, std::array<double, N>& b)
{
#ifdef CHEMISTRY_USE_EIGEN
    Eigen::Map<Eigen::Matrix<double, N, N, Eigen::RowMajor>> A(a.data());
    Eigen::Map<Eigen::Matrix<double, N, 1>> x(b.data());
    x = A.partialPivLu().solve(x);
#else
    constexpr double eps = numerics::eps_gaussj;
    std::array<int, N> indxc{}, indxr{}, ipiv{};
    ipiv.fill(0);

    // Track running max of |b[icol]*dum| for post-solve zero-flush
    std::array<double, N> b_max{};
    b_max.fill(0.0);

    int irow = 0, icol = 0;
    bool singular = false;

    for (int i = 0; i < N && !singular; ++i) {
        // Find pivot
        double big = 0.0;
        for (int j = 0; j < N; ++j) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < N; ++k) {
                    if (ipiv[k] == 0) {
                        double aabs = std::abs(a[j*N + k]);
                        if (aabs >= big) { big = aabs; irow = j; icol = k; }
                    } else if (ipiv[k] > 1) {
                        singular = true; break;
                    }
                }
            }
            if (singular) break;
        }
        if (singular) break;

        ++ipiv[icol];
        if (irow != icol) {
            for (int l = 0; l < N; ++l) std::swap(a[irow*N+l], a[icol*N+l]);
            std::swap(b[irow], b[icol]);
        }
        indxr[i] = irow;
        indxc[i] = icol;

        if (a[icol*N + icol] == 0.0) { singular = true; break; }

        double pivinv = 1.0 / a[icol*N + icol];
        a[icol*N + icol] = 1.0;
        for (int l = 0; l < N; ++l) a[icol*N + l] *= pivinv;
        b[icol] *= pivinv;

        for (int ll = 0; ll < N; ++ll) {
            if (ll != icol) {
                double dum = a[ll*N + icol];
                a[ll*N + icol] = 0.0;
                for (int l = 0; l < N; ++l) {
                    a[ll*N + l] -= a[icol*N + l] * dum;
                    // Flush near-zero fill-in
                    double ref = std::abs(a[icol*N+l] * dum);
                    if (ref != 0.0 && std::abs(a[ll*N+l]) / ref < eps)
                        a[ll*N+l] = 0.0;
                }
                double ref_b = std::abs(b[icol] * dum);
                b_max[ll] = std::max(b_max[ll], ref_b);
                b[ll] -= b[icol] * dum;
                if (ref_b != 0.0 && std::abs(b[ll]) / ref_b < eps)
                    b[ll] = 0.0;
            }
        }
    }

    // Flush residual noise in solution vector
    for (int ll = 0; ll < N; ++ll)
        if (b_max[ll] != 0.0 && std::abs(b[ll] / b_max[ll]) <= eps)
            b[ll] = 0.0;

    // Undo column permutation
    for (int l = N-1; l >= 0; --l)
        if (indxr[l] != indxc[l])
            for (int k = 0; k < N; ++k)
                std::swap(a[k*N + indxr[l]], a[k*N + indxc[l]]);
#endif
}

// ─────────────────────────────────────────────────────────────────────────────
// equichem — Saha equilibrium abundances (high-density branch xnH > 1e18)
//
//   Computes equilibrium constants K_eq for each Saha reaction using
//   eval_partition_functions, then iterates 2D Newton-Raphson on
//   charge neutrality (F_cha) and hydrogen conservation (F_hyd).
// ─────────────────────────────────────────────────────────────────────────────
template<int N_sp, int N_react>
void equichem(double xnH, double T_K,
              std::array<double, N_sp>& y,
              const ReactionTable<N_sp, N_react>& tbl)
{
    constexpr double xk_B = phys::xk_B;
    constexpr double h_P  = phys::h_P;
    constexpr double pi   = phys::pi;
    constexpr double yHe  = abundance_ref::yHe;
    constexpr double yD   = abundance_ref::yD;
    constexpr double yLi  = abundance_ref::yLi;
    constexpr double eps_it = numerics::eps_it_prim;

    // Compute partition functions
    std::array<double, 102> pf{};
    eval_partition_functions<N_sp, N_react>(T_K, tbl, pf);

    // Build equilibrium constants for each Saha entry
    // Keqb indexed by reaction number (1-based → 0-based stored at saha_num-1)
    const int N_keqb = N_react;
    std::array<double, N_react> Keqb{};
    Keqb.fill(0.0);

    double lnT32 = 1.5 * std::log(2.0*pi*xk_B*T_K / (h_P*h_P));

    for (int i = 0; i < tbl.n_saha; ++i) {
        int  num = tbl.saha_num[i];
        int  r1  = tbl.saha_rs1[i];   // 1-based species index
        int  r2  = tbl.saha_rs2[i];
        int  p1  = tbl.saha_ps1[i];
        int  p2  = tbl.saha_ps2[i];   // 101 = photon
        int  nr  = tbl.saha_nsr[i];
        int  np  = tbl.saha_nsp[i];
        double Cm = tbl.saha_Cm[i];
        double dE = tbl.saha_dE[i];

        double xlnC1  = double(nr - np) * lnT32;
        double xlnCpf = std::log(pf[r1]) + std::log(pf[r2]) - std::log(pf[p1]);
        if (np == 2 && p2 != 101)
            xlnCpf -= std::log(pf[p2]);

        double xlnKeqb = xlnC1 + std::log(Cm) + xlnCpf - dE / (xk_B * T_K);
        if (num >= 1 && num <= N_keqb)
            Keqb[num - 1] = std::exp(xlnKeqb);   // 0-based storage
    }

    // Convenience: Fortran Keqb(k) → Keqb[k-1]
    auto K = [&](int k) -> double { return Keqb[k-1]; };

    // Build abundance ratios (independent of y_e, y_H except where noted)
    // H sequence
    double K_Hp  = K(2)  / xnH;
    double K_Hm  = xnH   / K(7);
    double K_H2p = K(2)  / K(9);
    double K_H2  = xnH   / K(7) / K(8);
    double K_H3p = xnH * K(2) / K(7) / K(8) / K(9) / K(26);
    // He sequence
    double K_Hep  = K(4)  / xnH;
    double K_He2p = K(4) * K(6) / xnH / xnH;
    double K_HeHp = xnH   / K(44);
    // D sequence
    double K_Dp  = K(51) / xnH;
    double K_HD  = xnH   / K(54);
    double K_HDp = xnH   / K(60);
    double K_Dm  = xnH   / K(63);
    // Li sequence
    double K_Lip  = K(101) / xnH;
    double K_Li2p = K_Lip * K(121) / xnH;
    double K_Li3p = K_Li2p * K(122) / xnH;
    double K_LiH  = xnH / K(118);
    double K_LiHp = xnH / K(113);
    double K_Lim  = xnH / K(104);

    // Initial guess from current y
    double y_e = y[2];   // e-   (0-based: index 2)
    double y_H = y[0];   // H    (0-based: index 0)

    // Helper: fill all species from y_e, y_H
    auto fill_y = [&](double ye, double yh,
                       std::array<double, N_sp>& yy) {
        double y_Hp  = K_Hp  * yh / ye;
        double y_Hm  = K_Hm  * yh * ye;
        double y_H2p = K_H2p * yh*yh / ye;
        double y_H2  = K_H2  * yh*yh;
        double y_H3p = K_H3p * yh*yh*yh / ye;

        yy[0] = yh;            // H
        yy[1] = y_H2;          // H2
        yy[2] = ye;            // e-
        yy[3] = y_Hp;          // H+
        yy[4] = y_H2p;         // H2+
        yy[5] = y_H3p;         // H3+
        yy[6] = y_Hm;          // H-

        // He sequence
        {
            double den = 1.0 + K_Hep/ye + K_He2p/ye/ye + K_HeHp*y_Hp;
            yy[7]  = yHe / den;            // He
            yy[8]  = K_Hep  * yy[7] / ye; // He+
            yy[9]  = K_He2p * yy[7] / ye / ye; // He++
            yy[10] = K_HeHp * y_Hp * yy[7];    // HeH+
        }
        // D sequence
        {
            double den = 1.0 + K_HD*yh + K_Dp/ye + K_HDp*y_Hp + K_Dm*ye;
            yy[11] = yD / den;             // D
            yy[12] = K_HD  * yy[11] * yh; // HD
            yy[13] = K_Dp  * yy[11] / ye; // D+
            yy[14] = K_HDp * yy[11] * y_Hp; // HD+
            yy[15] = K_Dm  * yy[11] * ye; // D-
        }
        // Li sequence (only if N_sp >= 23)
        if constexpr (N_sp >= 23) {
            double den = 1.0 + K_LiH*yh + K_Lim*ye + K_LiHp*y_Hp
                        + K_Lip/ye + K_Li2p/ye/ye + K_Li3p/(ye*ye*ye);
            yy[16] = yLi / den;              // Li
            yy[17] = K_LiH  * yy[16] * yh;  // LiH
            yy[18] = K_Lip  * yy[16] / ye;  // Li+
            yy[19] = K_Lim  * yy[16] * ye;  // Li-
            yy[20] = K_LiHp * yy[16] * y_Hp;// LiH+
            yy[21] = K_Li2p * yy[16] / ye / ye; // Li++
            yy[22] = K_Li3p * yy[16] / (ye*ye*ye); // Li3+
        }
    };

    // Charge conservation:  Σ(z_i * y_i) = 0
    auto F_cha = [&](const std::array<double, N_sp>& yy) -> double {
        double val = yy[3] + yy[4] + yy[5] + yy[8] + yy[10] + yy[13] + yy[14]
                   + 2.0*(yy[9] + yy[21]) + 3.0*yy[22]
                   - (yy[2] + yy[6] + yy[15] + yy[19]);
        if constexpr (N_sp >= 23)
            val += yy[18] + yy[20]; // Li+, LiH+
        return val;
    };

    // Hydrogen conservation:  Σ(H-content_i * y_i) = 1
    auto F_hyd = [&](const std::array<double, N_sp>& yy) -> double {
        double val = yy[0] + yy[3] + yy[6] + yy[10] + yy[12] + yy[14]
                   + 2.0*(yy[1] + yy[4]) + 3.0*yy[5] - 1.0;
        if constexpr (N_sp >= 23)
            val += yy[17] + yy[20]; // LiH, LiH+
        return val;
    };

    std::array<double, N_sp> ytmp{};

    // 2D Newton-Raphson on (y_e, y_H)
    for (int itr = 0; itr < 100; ++itr) {
        fill_y(y_e, y_H, ytmp);
        double fc0 = F_cha(ytmp);
        double fh0 = F_hyd(ytmp);

        double conv = std::abs(fc0/y_e) + std::abs(fh0/y_H);
        if (conv <= 1.0e-10) break;

        // Numerical Jacobian (forward differences; base values are fc0, fh0)
        // Perturb y_e
        fill_y(y_e*(1.0+eps_it), y_H, ytmp);
        double dFc_dye = (F_cha(ytmp) - fc0) / (eps_it * y_e);
        double dFh_dye = (F_hyd(ytmp) - fh0) / (eps_it * y_e);

        // Perturb y_H
        fill_y(y_e, y_H*(1.0+eps_it), ytmp);
        double dFc_dyH = (F_cha(ytmp) - fc0) / (eps_it * y_H);
        double dFh_dyH = (F_hyd(ytmp) - fh0) / (eps_it * y_H);

        // Cramer's rule for 2×2 system
        double detA = dFc_dye * dFh_dyH - dFc_dyH * dFh_dye;
        double dy_e = -(dFh_dyH * fc0 - dFc_dyH * fh0) / detA;
        double dy_H =  (dFh_dye * fc0 - dFc_dye * fh0) / detA;

        y_e += dy_e;
        y_H += dy_H;
    }

    // Final fill
    fill_y(y_e, y_H, y);
}

// ─────────────────────────────────────────────────────────────────────────────
// chemreact — one timestep NR chemistry solver
//   Fortran chemreact in chemistry_prm.f.
//   Solves dy/dt = r_f(y) implicitly, optionally via Saha for xnH > 1e18.
//   Returns updated y, xmu, gamma, xLmbdch (chemistry cooling), var.
// ─────────────────────────────────────────────────────────────────────────────
template<int N_sp, int N_react>
void chemreact(double xnH, double T_K,
               std::array<double, N_sp>& y,
               double dt,
               double& xmu, double& gamma, double& xLmbdch,
               std::array<double, 2*N_react>& var,
               const ReactionTable<N_sp, N_react>& tbl,
               const ChemParams& params)
{
    constexpr double yHe  = abundance_ref::yHe;
    constexpr double xm_p = phys::xm_p;  // 1.67262e-24 g (high-precision; unified in Phase 5)
    constexpr double eps_y   = numerics::eps_y;
    constexpr double xnH_eq  = numerics::xnH_eq;

    std::array<double, N_sp>     y_init = y;
    std::array<double, N_sp>     dy{};    dy.fill(0.0);
    std::array<double, N_sp>     ddy{};
    std::array<double, N_sp>     r_f{};
    std::array<double, N_sp*N_sp> dr_fdy{};
    std::array<double, N_sp*N_sp> A_mat{};
    std::array<double, 2*N_react> xk{};

    if (xnH <= xnH_eq) {
        // ─── Non-Equilibrium (low-density) NR solver ───────────────────────
        // Set xcrit from current y, matching Fortran COMMON /xcrit/:
        //   xH = y(1), xH2 = 2*y(2), xHe = y(8)  (Fortran 1-based)
        ChemParams p_loc = params;
        p_loc.xH  = y[0];
        p_loc.xH2 = 2.0 * y[1];
        p_loc.xHe = y[7];

        // Compute reaction rate coefficients
        xmu = (1.0 + 4.0*yHe)
            / (y[0] + y[1] + y[2] + y[3] + y[7] + y[8] + y[9]);
        compute_base_rates<N_sp, N_react>(xnH, T_K, xmu, p_loc, tbl, xk);

        for (int itr = 0; itr < 30; ++itr) {
            r_f.fill(0.0);
            dr_fdy.fill(0.0);

            compute_rates<N_sp, N_react>(xk, xnH, y, tbl, r_f, dr_fdy, var);

            // Build A = I - dt * J  (J = dr_fdy, row-major)
            for (int isp = 0; isp < N_sp; ++isp)
                for (int jsp = 0; jsp < N_sp; ++jsp)
                    A_mat[isp*N_sp + jsp] = (isp == jsp ? 1.0 : 0.0)
                                          - dt * dr_fdy[isp*N_sp + jsp];

            // rhs: ddy = r_f*dt - dy
            double err_fnc = 0.0;
            for (int isp = 0; isp < N_sp; ++isp) {
                ddy[isp]  = r_f[isp]*dt - dy[isp];
                err_fnc  += std::abs(ddy[isp]);
            }

            // Solve A * ddy = ddy  (in-place)
            gaussj_solve<N_sp>(A_mat, ddy);

            // Update abundances
            double err_y = 0.0;
            for (int isp = 0; isp < N_sp; ++isp) {
                if (y[isp] + ddy[isp] < 0.0) {
                    ddy[isp] = -0.1 * y[isp];
                    err_y += 1.0;
                }
                dy[isp] += ddy[isp];
                y[isp]  += ddy[isp];
                err_y   += std::abs(ddy[isp]);
            }

            if (err_y <= eps_y && err_fnc <= 1.0e-4) break;
        }
    } else {
        // ─── Saha equilibrium (high-density) ───────────────────────────────
        equichem<N_sp, N_react>(xnH, T_K, y, tbl);
        for (int isp = 0; isp < N_sp; ++isp)
            dy[isp] = y[isp] - y_init[isp];
    }

    // Update thermodynamic variables
    xmu = (1.0 + 4.0*yHe)
        / (y[0] + y[1] + y[2] + y[3] + y[7] + y[8] + y[9]);

    gamma = 1.0 + (1.0 + 4.0*yHe)
          / (xmu * (1.5*(y[0] + y[2] + y[3] + y[7] + y[8] + y[9])
                   + c_H2(T_K)*y[1]));

    // ── Chemistry cooling xLmbdch ──────────────────────────────────────────
    double xL_H2   = 0.0;
    double dyHp    = dy[3];   // d(H+)
    double dyHep   = dy[8];   // d(He+)
    double dyHepp  = dy[9];   // d(He++)

    if (xnH < 1.0e13) {
        // H2 formation/dissociation cooling
        double xn_cr = 1.0e6 / std::sqrt(T_K)
                     / (1.6*y[0]*std::exp(-std::pow(400.0/T_K, 2.0))
                       + 1.4*y[1]*std::exp(-12000.0/(T_K + 1200.0)));
        double crit  = 1.0 / (1.0 + xn_cr/xnH);

        // H- + H → H2 + e  (reaction 8, 0-based: xk[7])
        double rtHm  = xk[7]  * y[6] * y[0];
        // H2+ + H → H2 + H+ (reaction 10, 0-based: xk[9])
        double rtH2p = xk[9]  * y[4] * y[0];
        // 3-body formation: 3H→H2+H, 2H+H2→2H2, 2H+He→H2+He
        //   uses: xk[18] (forward r19), xk[20+N_react] (reverse r21), xk[36+N_react] (reverse r37)
        double rt3b  = (xk[18]          * y[0]*y[0]*y[0]
                       + xk[20+N_react] * y[0]*y[0]*y[1]
                       + xk[36+N_react] * y[0]*y[0]*y[7]) * xnH;
        // Dissociation: H2+H→3H, H2+H2→2H+H2, H2+He→2H+He
        //   uses: xk[18+N_react] (reverse r19), xk[20] (forward r21), xk[36] (forward r37)
        double rtdis = xk[18+N_react] * y[1]*y[0]
                     + xk[20]          * y[1]*y[1]
                     + xk[36]          * y[1]*y[7];

        xL_H2 = -(rtHm  * 3.53 * crit
               +  rtH2p * 1.83 * crit
               +  (rt3b*crit - rtdis) * 4.48)
               * xnH * 1.60219e-12;

        // Chemistry ionization contributions
        dyHp  = dy[3]
              + (xk[1]*y[3]*y[2]*xnH
                - (xk[130]+xk[136])*y[0]
                - (xk[131]+xk[134])*y[1]) * dt;
        dyHep = dy[8]
              + (xk[3]*y[8]*y[2]*xnH - (xk[135]+xk[138])*y[7]) * dt;
    } else {
        double dyH2 = dy[1];
        xL_H2  = -7.18e-12 * dyH2 / dt;
    }

    // Specific cooling rate [erg g^-1 s^-1]
    xLmbdch = ((2.18e-11*dyHp + 3.94e-11*dyHep + 12.66e-11*dyHepp) / dt
               + xL_H2) / ((1.0 + 4.0*yHe) * xm_p);
}

// ─────────────────────────────────────────────────────────────────────────────
// chemcool — secant method wrapper for temperature self-consistent integration
//   Fortran chemcool in chemistry_prm.f.
//   Finds T* such that T_K - T* - (γ-1)*μ*mp*ΛΔt/kB = 0,
//   then updates y, xmu, gamma, xLmbdch at converged T*.
//   Skip secant loop when xnH > 1e16 and T < 1650 K (Fortran condition).
// ─────────────────────────────────────────────────────────────────────────────
template<int N_sp, int N_react>
void chemcool(double xnH, double Tp,
              std::array<double, N_sp>& y,
              double dt,
              double& xmu, double& gamma, double& xLmbdch,
              std::array<double, 2*N_react>& var,
              const ReactionTable<N_sp, N_react>& tbl,
              const ChemParams& params)
{
    constexpr double xm_p  = phys::xm_p;  // 1.67262e-24 g (high-precision; unified in Phase 5)
    constexpr double xk_B  = phys::xk_B;  // 1.380662e-16 erg/K (high-precision; unified in Phase 5)
    constexpr int    maxit  = 100;

    // Residual: f(T) = Tp - T - (γ-1)*μ*mp*Λ*dt/kB
    auto func = [&](double t, double xLc) -> double {
        return Tp - t - (gamma - 1.0) * xmu * xm_p * xLc * dt / xk_B;
    };

    // First evaluation at Tp1 = Tp
    double Tp1 = Tp;
    std::array<double, N_sp> ytmp = y;
    double xmu1 = xmu, gamma1 = gamma, xLmbdch1 = 0.0;
    chemreact<N_sp, N_react>(xnH, Tp1, ytmp, dt, xmu1, gamma1, xLmbdch1, var, tbl, params);
    double fL = func(Tp1, xLmbdch1);

    // Skip secant if low-temperature dense condition
    if (xnH > 1.0e16 && Tp <= 1.65e3) {
        y       = ytmp;
        xmu     = xmu1;
        gamma   = gamma1;
        xLmbdch = xLmbdch1;
        return;
    }

    // Second evaluation at Tp2 = 0.999*Tp1
    double Tp2 = 0.999 * Tp1;
    ytmp = y;
    double xmu2 = xmu, gamma2 = gamma, xLmbdch2 = 0.0;
    chemreact<N_sp, N_react>(xnH, Tp2, ytmp, dt, xmu2, gamma2, xLmbdch2, var, tbl, params);
    double f = func(Tp2, xLmbdch2);

    // Ensure |fL| >= |f|  (larger residual at TpL)
    double TpL, Tpsec;
    if (std::abs(fL) < std::abs(f)) {
        TpL = Tp2; Tpsec = Tp1;
        std::swap(fL, f);
    } else {
        TpL = Tp1; Tpsec = Tp2;
    }

    // Secant iterations
    for (int i = 0; i < maxit; ++i) {
        double dTp = (TpL - Tpsec) * f / (f - fL);
        TpL  = Tpsec;
        fL   = f;
        Tpsec += dTp;

        ytmp = y;
        chemreact<N_sp, N_react>(xnH, Tpsec, ytmp, dt,
                                 xmu1, gamma1, xLmbdch1, var, tbl, params);
        f = func(Tpsec, xLmbdch1);

        if (std::abs(f / Tpsec) <= 1.0e-7) break;
    }

    y       = ytmp;
    xmu     = xmu1;
    gamma   = gamma1;
    xLmbdch = xLmbdch1;
}

}  // namespace chemistry
