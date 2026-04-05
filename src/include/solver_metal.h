// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// solver_metal.h — metal_grain solver specializations (89 species, 1200 reactions)
//
// Template specializations of chemreact and chemcool for <89, 1200>, plus
// a detail::equichem_metal helper (receives Z_metal that the generic
// equichem signature lacks).
//
// Port of:
//   fortran/metal_grain/subs/saha_metal.f      — equichem
//   fortran/metal_grain/subs/chemistry_metal.f — chemreact, chemcool
//
// Must be included AFTER solver.h and reaction_metal.h.
// Species index: Fortran y(N) → C++ y[N-1]  (0-based).
// ---------------------------------------------------------------------------
#include "solver.h"
#include "reaction_metal.h"

namespace chemistry {

// ─────────────────────────────────────────────────────────────────────────────
// detail::equichem_metal
//
// Saha equilibrium for the metal_grain network.
// 4D Newton–Raphson on (y_e, y_H, y_C, y_O).
// Fills y[0..62]; leaves y[63..88] unchanged.
//
// Port of equichem() in saha_metal.f.
// ─────────────────────────────────────────────────────────────────────────────
namespace detail {

inline void equichem_metal(double xnH, double T_K, double Z_metal,
                            std::array<double, 89>& y,
                            const ReactionTable<89, 1200>& tbl)
{
    constexpr double xk_B  = phys::xk_B;
    constexpr double h_P   = phys::h_P;
    constexpr double pi    = phys::pi;
    // Solar reference abundances
    constexpr double yHe_s  = abundance_ref::yHe;
    constexpr double yD_s   = abundance_ref::yD;
    constexpr double yLi_s  = abundance_ref::yLi;
    constexpr double yC_s   = abundance_ref::yC;
    constexpr double yO_s   = abundance_ref::yO;
    constexpr double yNa_s  = abundance_ref::yNa;
    constexpr double yMg_s  = abundance_ref::yMg;
    constexpr double yK_s   = abundance_ref::yK;
    constexpr double eps_it = numerics::eps_it_metal;

    const double yHe = yHe_s;
    const double yD  = yD_s;
    const double yLi = yLi_s;
    const double XC  = yC_s  * Z_metal;
    const double XO  = yO_s  * Z_metal;
    const double XKa = yK_s  * Z_metal;
    const double XNa = yNa_s * Z_metal;
    const double XMg = yMg_s * Z_metal;

    // ── Partition functions ──────────────────────────────────────────────────
    std::array<double, 102> pf{};
    eval_partition_functions<89, 1200>(T_K, tbl, pf);

    // ── Build Keqb[0..1199] indexed by (reaction_number - 1) ─────────────────
    // Keqb[num-1] = exp( 1.5*(nr-np)*ln(2πkT/h²) + ln(Cm)
    //                   + Σln(pf_reactants) - Σln(pf_products) - dE/(kT) )
    std::array<double, 1200> Keqb{};
    Keqb.fill(0.0);
    {
        const double lnT32 = 1.5 * std::log(2.0*pi*xk_B*T_K / (h_P*h_P));
        for (int i = 0; i < tbl.n_saha; ++i) {
            int    num = tbl.saha_num[i];
            int    r1  = tbl.saha_rs1[i];
            int    r2  = tbl.saha_rs2[i];
            int    p1  = tbl.saha_ps1[i];
            int    p2  = tbl.saha_ps2[i];
            int    nr  = tbl.saha_nsr[i];
            int    np  = tbl.saha_nsp[i];
            double Cm  = tbl.saha_Cm[i];
            double dE  = tbl.saha_dE[i];

            double xlnCpf = std::log(pf[r1]) + std::log(pf[r2])
                          - std::log(pf[p1]);
            if (np == 2 && p2 != 101)
                xlnCpf -= std::log(pf[p2]);

            double xlnKeqb = double(nr - np)*lnT32 + std::log(Cm)
                           + xlnCpf - dE / (xk_B * T_K);

            if (num >= 1 && num <= 1200)
                Keqb[num - 1] = std::exp(xlnKeqb);
        }
    }

    // ── Equilibrium ratio constants (Fortran Keqb(N) → C++ Keqb[N-1]) ────────
    // H sequence
    double K_Hp  = Keqb[1]  / xnH;
    double K_Hm  = xnH      / Keqb[6];
    double K_H2p = Keqb[1]  / Keqb[8];
    double K_H2  = xnH      / Keqb[6] / Keqb[7];
    double K_H3p = xnH * Keqb[1] / (Keqb[6]*Keqb[7]*Keqb[8]*Keqb[25]);
    // He sequence
    double K_Hep  = Keqb[3]  / xnH;
    double K_He2p = Keqb[3]  * Keqb[5] / (xnH*xnH);
    double K_HeHp = xnH      / Keqb[640];
    // D sequence
    double K_Dp   = Keqb[50] / xnH;
    double K_HD   = xnH      / Keqb[53];
    double K_HDp  = xnH      / Keqb[59];
    double K_Dm   = xnH      / Keqb[62];
    // Li sequence
    double K_Lip  = Keqb[800] / xnH;
    double K_Li2p = K_Lip     * Keqb[820] / xnH;
    double K_Li3p = K_Li2p    * Keqb[821] / xnH;
    double K_LiH  = xnH       / Keqb[817];
    double K_LiHp = xnH       / Keqb[812];
    double K_Lim  = xnH       / Keqb[803];

    // ── fill_y: write y[0..62] from (y_e, y_H, y_C, y_O) ───────────────────
    // Mirrors the "0 th" block repeated in each Fortran perturbation step.
    auto fill_y = [&](double ye, double yh, double yc, double yo) {
        double yHp  = K_Hp  * yh / ye;
        double yHm  = K_Hm  * yh * ye;
        double yH2p = K_H2p * yh * yh / ye;
        double yH2  = K_H2  * yh * yh;
        double yH3p = K_H3p * yh * yh * yh / ye;

        y[0] = yh;    // H
        y[1] = yH2;   // H2
        y[2] = ye;    // e-
        y[3] = yHp;   // H+
        y[4] = yH2p;  // H2+
        y[5] = yH3p;  // H3+
        y[6] = yHm;   // H-

        // He sequence
        y[7]  = yHe / (1.0 + K_Hep/ye + K_He2p/(ye*ye) + K_HeHp*yHp);
        y[8]  = K_Hep  * y[7] / ye;
        y[9]  = K_He2p * y[7] / (ye*ye);
        y[10] = K_HeHp * yHp  * y[7];

        // D sequence
        y[11] = yD / (1.0 + K_HD*yh + K_Dp/ye + K_HDp*yHp + K_Dm*ye);
        y[12] = K_HD  * y[11] * yh;
        y[13] = K_Dp  * y[11] / ye;
        y[14] = K_HDp * y[11] * yHp;
        y[15] = K_Dm  * y[11] * ye;

        // Li sequence  (Fortran y(51..57) → C++ y[50..56])
        y[50] = yLi / (1.0 + K_LiH*yh + K_Lim*ye + K_LiHp*yHp
                      + K_Lip/ye + K_Li2p/(ye*ye) + K_Li3p/(ye*ye*ye));
        y[51] = K_LiH  * y[50] * yh;
        y[52] = K_Lip  * y[50] / ye;
        y[53] = K_Lim  * y[50] * ye;
        y[54] = K_LiHp * y[50] * yHp;
        y[55] = K_Li2p * y[50] / (ye*ye);
        y[56] = K_Li3p * y[50] / (ye*ye*ye);

        // Alkali ion sequences — y(59),y(61),y(63) are K+, Na+, Mg+ (appear in F_cha)
        // Fortran y(59)→C++ y[58]=K+, y(61)→y[60]=Na+, y(63)→y[62]=Mg+
        y[58] = XKa / (1.0 + xnH*ye / Keqb[700]);  // K+  ← Keqb(701)
        y[60] = XNa / (1.0 + xnH*ye / Keqb[702]);  // Na+ ← Keqb(703)
        y[62] = XMg / (1.0 + xnH*ye / Keqb[717]);  // Mg+ ← Keqb(718)

        // Carbon sequence  (Fortran y(17..29) → C++ y[16..28])
        y[16] = yc;
        y[17] = yc*yc*xnH / Keqb[185];                // C2 ← Keqb(186)
        y[18] = yh*yc*xnH  / Keqb[184];               // CH ← Keqb(185)
        y[19] = y[1]*yc*xnH / Keqb[537];              // CH2 ← Keqb(538)
        y[20] = y[1]*y[18]*xnH / Keqb[642];           // CH3 ← Keqb(643)
        y[21] = y[1]*y[20] / (yh * Keqb[135]);        // CH4 ← Keqb(136)
        y[22] = (yc/ye)*Keqb[506] / xnH;              // C+  ← Keqb(507)
        y[23] = yc*y[22]*xnH / Keqb[643];             // C2+ ← Keqb(644)
        y[24] = yh*y[22]*xnH / Keqb[294];             // CH+ ← Keqb(295)
        y[25] = y[1]*y[22]*xnH / Keqb[297];           // CH2+ ← Keqb(298)
        y[26] = (y[20]/ye)*Keqb[513] / xnH;           // CH3+ ← Keqb(514)
        y[27] = y[3]*y[21] / (yh * Keqb[212]);        // CH4+ ← Keqb(213)
        y[28] = y[1]*y[26]*xnH / Keqb[339];           // CH5+ ← Keqb(340)

        // Oxygen sequence  (Fortran y(30..50) → C++ y[29..49])
        y[29] = yo;
        y[30] = yo*yo*xnH / Keqb[193];                // O2 ← Keqb(194)
        y[31] = yh*yo*xnH / Keqb[192];                // OH ← Keqb(193)
        y[32] = yc*yo*xnH / Keqb[186];                // CO ← Keqb(187)
        y[33] = yh*y[31]*xnH / Keqb[641];             // H2O ← Keqb(642)
        y[34] = Keqb[183]*y[1]*y[32] / yh;            // HCO ← Keqb(184)
        y[35] = Keqb[114]*yo*y[33] / yh;              // H2CO ← Keqb(115)
        y[36] = yo*y[34] / (Keqb[200]*yh);            // CO2 ← Keqb(201)
        y[37] = Keqb[111]*y[1]*y[34] / yh;            // ← Keqb(112)
        y[38] = Keqb[601]*y[31]*y[33] / yh;           // ← Keqb(602)
        y[39] = (yo/ye)*Keqb[514] / xnH;              // O+  ← Keqb(515)
        y[40] = Keqb[533]*yo*yo / ye;                 // O2+ ← Keqb(534)
        y[41] = Keqb[517]*yh*yo / ye;                 // OH+ ← Keqb(518)
        y[42] = yc*y[39]*xnH / Keqb[644];             // CO+ ← Keqb(645)
        y[43] = Keqb[520]*yh*y[31] / ye;              // H2O+ ← Keqb(521)
        y[44] = Keqb[526]*yh*y[32] / ye;              // HCO+ ← Keqb(527)
        y[45] = Keqb[534]*yh*y[30] / ye;              // O2H+ ← Keqb(535)
        y[46] = Keqb[522]*yh*y[33] / ye;              // H3O+ ← Keqb(523)
        y[47] = (y[37]/ye)*Keqb[529] / xnH;           // ← Keqb(530)
        y[48] = Keqb[535]*yh*y[36] / ye;              // ← Keqb(536)
        y[49] = Keqb[532]*yh*y[37] / ye;              // ← Keqb(533)
    };

    // ── Conservation residuals (lambdas read from y[]) ────────────────────────
    // F_cha: charge conservation Σ(z_i * y_i) = 0
    auto F_cha_fn = [&]() -> double {
        return y[3]+y[4]+y[5]+y[8]+y[10]+y[13]+y[14]
              +y[22]+y[23]+y[24]+y[25]+y[26]+y[27]+y[28]
              +y[39]+y[40]+y[41]+y[42]+y[43]+y[44]+y[45]
              +y[46]+y[47]+y[49]
              +y[52]+y[54]+y[58]+y[60]+y[62]
              +2.0*(y[9]+y[55])+3.0*y[56]
              -(y[2]+y[6]+y[15]+y[53]);
    };
    // F_hyd: hydrogen conservation Σ(n_H(i) * y_i) = 1
    auto F_hyd_fn = [&]() -> double {
        return y[0]+y[3]+y[6]+y[10]+y[12]+y[14]
              +y[18]+y[24]+y[31]+y[34]+y[35]+y[41]+y[44]+y[45]+y[48]
              +y[51]+y[54]
              +2.0*(y[1]+y[4]+y[19]+y[25]+y[33]+y[37]+y[38]+y[43]+y[47])
              +3.0*(y[5]+y[20]+y[26]+y[46]+y[49])
              +4.0*(y[21]+y[27])+5.0*y[28]-1.0;
    };
    // F_car: carbon conservation Σ(n_C(i) * y_i) = XC
    auto F_car_fn = [&]() -> double {
        return y[16]+y[18]+y[19]+y[20]+y[21]+y[22]
              +y[24]+y[25]+y[26]+y[27]+y[28]+y[32]+y[34]+y[36]+y[37]
              +y[42]+y[44]+y[47]+y[48]+y[49]
              +2.0*(y[17]+y[23])-XC;
    };
    // F_oxy: oxygen conservation Σ(n_O(i) * y_i) = XO
    auto F_oxy_fn = [&]() -> double {
        return y[29]+y[31]+y[32]+y[33]+y[34]+y[37]
              +y[39]+y[41]+y[42]+y[43]+y[44]+y[46]+y[47]+y[49]
              +2.0*(y[30]+y[35]+y[36]+y[38]+y[40]+y[45]+y[48])-XO;
    };

    // ── Initial guess from current y ──────────────────────────────────────────
    double y_e = y[2];
    double y_H = y[0];
    double y_C = y[16];
    double y_O = y[29];
    if (y_e <= 0.0) y_e = 1.0e-20;
    if (y_H <= 0.0) y_H = 1.0e-20;
    if (y_C <= 0.0) y_C = (Z_metal > 0.0) ? 1.0e-20 : 0.0;
    if (y_O <= 0.0) y_O = (Z_metal > 0.0) ? 1.0e-20 : 0.0;

    // ── 4D Newton–Raphson loop ─────────────────────────────────────────────────
    for (int itr = 0; itr < 100; ++itr) {

        // ── 0. Base evaluation ─────────────────────────────────────────────────
        fill_y(y_e, y_H, y_C, y_O);
        double Fc0   = F_cha_fn();
        double Fh0   = F_hyd_fn();
        double Fcar0 = F_car_fn();
        double Fo0   = F_oxy_fn();

        // ── 1. Convergence check ───────────────────────────────────────────────
        double conv = std::abs(Fc0/y_e) + std::abs(Fh0/y_H);
        if (y_C > 0.0 && XC > 0.0) conv += std::abs(Fcar0/y_C);
        if (y_O > 0.0 && XO > 0.0) conv += std::abs(Fo0  /y_O);
        if (conv <= 1.0e-10) break;

        // ── 2. Numerical Jacobian (forward differences) ────────────────────────
        // Column 0: perturb y_e
        fill_y(y_e*(1.0+eps_it), y_H, y_C, y_O);
        double dFc_de   = (F_cha_fn() - Fc0)   / (eps_it * y_e);
        double dFh_de   = (F_hyd_fn() - Fh0)   / (eps_it * y_e);
        double dFcar_de = (F_car_fn() - Fcar0) / (eps_it * y_e);
        double dFo_de   = (F_oxy_fn() - Fo0)   / (eps_it * y_e);

        // Column 1: perturb y_H
        fill_y(y_e, y_H*(1.0+eps_it), y_C, y_O);
        double dFc_dH   = (F_cha_fn() - Fc0)   / (eps_it * y_H);
        double dFh_dH   = (F_hyd_fn() - Fh0)   / (eps_it * y_H);
        double dFcar_dH = (F_car_fn() - Fcar0) / (eps_it * y_H);
        double dFo_dH   = (F_oxy_fn() - Fo0)   / (eps_it * y_H);

        // Column 2: perturb y_C
        double dFc_dC = 0.0, dFh_dC = 0.0, dFcar_dC = 0.0, dFo_dC = 0.0;
        if (y_C > 0.0 && XC > 0.0) {
            fill_y(y_e, y_H, y_C*(1.0+eps_it), y_O);
            dFc_dC   = (F_cha_fn() - Fc0)   / (eps_it * y_C);
            dFh_dC   = (F_hyd_fn() - Fh0)   / (eps_it * y_C);
            dFcar_dC = (F_car_fn() - Fcar0) / (eps_it * y_C);
            dFo_dC   = (F_oxy_fn() - Fo0)   / (eps_it * y_C);
        }

        // Column 3: perturb y_O
        double dFc_dO = 0.0, dFh_dO = 0.0, dFcar_dO = 0.0, dFo_dO = 0.0;
        if (y_O > 0.0 && XO > 0.0) {
            fill_y(y_e, y_H, y_C, y_O*(1.0+eps_it));
            dFc_dO   = (F_cha_fn() - Fc0)   / (eps_it * y_O);
            dFh_dO   = (F_hyd_fn() - Fh0)   / (eps_it * y_O);
            dFcar_dO = (F_car_fn() - Fcar0) / (eps_it * y_O);
            dFo_dO   = (F_oxy_fn() - Fo0)   / (eps_it * y_O);
        }

        // ── 3. Solve  J · δ = −F  (row-major: [F_cha, F_hyd, F_car, F_oxy]) ──
        // J[i][j] = dF_i/dx_j,  x = (y_e, y_H, y_C, y_O)
        std::array<double, 16> J_mat = {
            dFc_de,   dFc_dH,   dFc_dC,   dFc_dO,
            dFh_de,   dFh_dH,   dFh_dC,   dFh_dO,
            dFcar_de, dFcar_dH, dFcar_dC, dFcar_dO,
            dFo_de,   dFo_dH,   dFo_dC,   dFo_dO
        };
        std::array<double, 4> rhs = { -Fc0, -Fh0, -Fcar0, -Fo0 };
        gaussj_solve<4>(J_mat, rhs);

        // ── 4. Update and clamp positive ──────────────────────────────────────
        y_e = std::max(y_e + rhs[0], 1.0e-30);
        y_H = std::max(y_H + rhs[1], 1.0e-30);
        if (XC > 0.0) y_C = std::max(y_C + rhs[2], 1.0e-30);
        if (XO > 0.0) y_O = std::max(y_O + rhs[3], 1.0e-30);
    }

    // ── Final fill at converged (y_e, y_H, y_C, y_O) ─────────────────────────
    fill_y(y_e, y_H, y_C, y_O);

    // ── Neutral alkali abundances (Fortran: filled after NR loop) ────────────
    // Fortran y(58)=K_neutral, y(60)=Na_neutral, y(62)=Mg_neutral
    // C++ y[57]=K, y[59]=Na, y[61]=Mg  (species.h: K=57, Na=59, Mg=61)
    // Formula: K_neutral = K+ * y_e * xnH / Keqb(K_ion)
    // Guard: Keqb can underflow to 0 at low T → division by zero
    y[57] = (Keqb[700] > 1.0e-300) ? xnH * y_e * y[58] / Keqb[700] : 0.0;
    y[59] = (Keqb[702] > 1.0e-300) ? xnH * y_e * y[60] / Keqb[702] : 0.0;
    y[61] = (Keqb[717] > 1.0e-300) ? xnH * y_e * y[62] / Keqb[717] : 0.0;
}

} // namespace detail


// ─────────────────────────────────────────────────────────────────────────────
// chemreact<89, 1200>
//
// One-timestep implicit NR chemistry solver for the metal_grain network.
// Key differences from the generic template:
//   - xJH2, xJH2O, xJtot are read from species array (y[70], y[79], y[69]+y[72])
//   - Equilibrium branch calls equichem_metal (passes Z_metal)
//   - rtgrain includes grain-catalysed H2 formation channels
//   - CR ionisation uses metal_grain indices (xk[543], xk[676..677] etc.)
//
// Port of chemreact() in chemistry_metal.f.
// ─────────────────────────────────────────────────────────────────────────────

// TODO: define number of species and reactions
template<>
inline void chemreact<89, 1200>(
    double xnH, double T_K,
    std::array<double, 89>& y,
    double dt,
    double& xmu, double& gamma, double& xLmbdch,
    std::array<double, 2400>& var,
    const ReactionTable<89, 1200>& tbl,
    const ChemParams& params)
{
    constexpr double yHe   = abundance_ref::yHe;
    constexpr double xm_p  = phys::xm_p;  // 1.67262e-24 g (high-precision; unified in Phase 5)
    constexpr double eps_y = numerics::eps_y;
    constexpr double xnH_eq = numerics::xnH_eq;

    std::array<double, 89>   y_init = y;
    std::array<double, 89>   dy{};   dy.fill(0.0);
    std::array<double, 89>   ddy{};
    std::array<double, 89>   r_f{};
    std::array<double, 89*89> dr_fdy{};
    std::array<double, 89*89> A_mat{};
    std::array<double, 2400>  xk{};

    if (xnH <= xnH_eq) {
        // ── Non-equilibrium (low-density) NR solver ──────────────────────────
        // xcrit abundances; xJH2/xJH2O/xJtot come from UV-shield species.
        ChemParams p_loc = params;
        p_loc.xH    = y[0];
        p_loc.xH2   = 2.0 * y[1];
        p_loc.xHe   = y[7];
        p_loc.xJH2  = y[70];          // Fortran y(71)
        p_loc.xJH2O = y[79];          // Fortran y(80)
        p_loc.xJtot = y[69] + y[72];  // Fortran y(70)+y(73)

        xmu = (1.0 + 4.0*yHe)
            / (y[0]+y[1]+y[2]+y[3]+y[7]+y[8]+y[9]);

        compute_base_rates<89, 1200>(xnH, T_K, xmu, p_loc, tbl, xk);

        // TODO: define max iterations
        for (int itr = 0; itr < 60; ++itr) {
            r_f.fill(0.0);
            dr_fdy.fill(0.0);

            compute_rates<89, 1200>(xk, xnH, y, tbl, r_f, dr_fdy, var);

            // Build A = I − dt·J
            for (int isp = 0; isp < 89; ++isp)
                for (int jsp = 0; jsp < 89; ++jsp)
                    A_mat[isp*89 + jsp] = (isp == jsp ? 1.0 : 0.0)
                                        - dt * dr_fdy[isp*89 + jsp];

            // rhs: ddy = r_f*dt − dy
            double err_fnc = 0.0;
            for (int isp = 0; isp < 89; ++isp) {
                ddy[isp]  = r_f[isp]*dt - dy[isp];
                err_fnc  += std::abs(ddy[isp]);
            }

            gaussj_solve<89>(A_mat, ddy);

            double err_y = 0.0;
            for (int isp = 0; isp < 89; ++isp) {
                if (y[isp] + ddy[isp] < 0.0) {
                    ddy[isp] = -0.1 * y[isp];
                    err_y += 1.0;
                }
                dy[isp] += ddy[isp];
                y[isp]  += ddy[isp];
                err_y   += std::abs(ddy[isp]);
            }
            // Floor clamp: prevent cumulative negative drift
            for (int isp = 0; isp < 89; ++isp)
                y[isp] = std::max(y[isp], 0.0);

            // Guard: NR divergence from near-singular Jacobian (specific Z/density
            // combinations can make det(A)≈0, producing huge or non-finite ddy).
            // Reset to initial state to prevent NaN propagation into xmu/gamma.
            bool nr_ok = true;
            for (int isp = 0; isp < 89 && nr_ok; ++isp)
                if (!std::isfinite(y[isp]) || y[isp] > 1.0e50)
                    nr_ok = false;
            if (!nr_ok) { y = y_init; dy.fill(0.0); break; }

            if (err_y <= eps_y && err_fnc <= 1.0e-4) break;
        }

    } else {
        // ── Saha equilibrium (high-density) — only first 63 species ──────────
        detail::equichem_metal(xnH, T_K, params.Z_metal, y, tbl);
        for (int i = 0; i < 63; ++i)
            dy[i] = y[i] - y_init[i];
        // dy[63..88] = 0 (grain/UV-shield species unchanged)
    }

    // ── Thermodynamic update ──────────────────────────────────────────────────
    xmu = (1.0 + 4.0*yHe)
        / (y[0]+y[1]+y[2]+y[3]+y[7]+y[8]+y[9]);

    gamma = 1.0 + (1.0 + 4.0*yHe)
          / (xmu * (1.5*(y[0]+y[2]+y[3]+y[7]+y[8]+y[9])
                   + c_H2(T_K)*y[1]));

    // ── Chemistry cooling xLmbdch ─────────────────────────────────────────────
    double xL_H2   = 0.0;
    double dyHp    = dy[3];
    double dyHep   = dy[8];
    double dyHepp  = dy[9];

    if (xnH < 1.0e13) {
        double denom_cr = 1.6*y[0]*std::exp(-std::pow(400.0/T_K, 2.0))
                        + 1.4*y[1]*std::exp(-12000.0/(T_K + 1200.0));
        double xn_cr = (denom_cr > 1.0e-30)
                     ? 1.0e6 / std::sqrt(T_K) / denom_cr : 1.0e50;
        double crit  = 1.0 / (1.0 + xn_cr/xnH);

        // Grain-catalysed H2 formation (5 channels, Fortran L206-209)
        // xk[22]   = xk(23)   3-body on grain
        // xk[1048] = xk(1049) surface:  gr + gr
        // xk[1108] = xk(1109) surface:  gr + dust
        // xk[1107] = xk(1108) surface:  H + dust + H2
        // xk[1114] = xk(1115) surface:  H + gr + H2
        double rtgrain = xk[22]*y[0]
                       + xk[1048]*y[68]*y[68]
                       + xk[1108]*y[68]*y[69]
                       + xk[1107]*y[0]*y[69]*xnH
                       + xk[1114]*y[0]*y[68]*xnH;

        // H- + H → H2 + e-  (reaction 8, 0-based xk[7])
        double rtHm  = xk[7]  * y[6] * y[0];
        // H2+ + H → H2 + H+ (reaction 10, 0-based xk[9])
        double rtH2p = xk[9]  * y[4] * y[0];
        // 3H, 2H+H2, 2H+He three-body formation
        double rt3b  = (xk[18]          * y[0]*y[0]*y[0]
                       + xk[20+1200]    * y[0]*y[0]*y[1]
                       + xk[36+1200]    * y[0]*y[0]*y[7]) * xnH;
        // H2 collisional dissociation
        double rtdis = xk[18+1200]*y[1]*y[0]
                     + xk[20]         *y[1]*y[1]
                     + xk[36]         *y[1]*y[7];

        xL_H2 = -(rtgrain*(0.2 + 4.2*crit)
               +  rtHm   * 3.53 * crit
               +  rtH2p  * 1.83 * crit
               +  (rt3b*crit - rtdis) * 4.48)
               * xnH * 1.60219e-12;

        // CR/photo ionisation contributions to cooling
        // Fortran: xk(2)→xk[1], xk(544)→xk[543], xk(677)→xk[676]
        //          xk(548)→xk[547], xk(551)→xk[550]
        dyHp = dy[3]
             + (xk[1]*y[3]*y[2]*xnH
               - (xk[543]+xk[676])*y[0]
               - (xk[547]+xk[550])*y[1]) * dt;
        // Fortran: xk(4)→xk[3], xk(545)→xk[544], xk(678)→xk[677]
        dyHep = dy[8]
              + (xk[3]*y[8]*y[2]*xnH
                - (xk[544]+xk[677])*y[7]) * dt;

    } else {
        double dyH2 = dy[1];
        xL_H2  = -7.18e-12 * dyH2 / dt;
        // Clamp to prevent extreme cooling from NR noise
        constexpr double xL_H2_max = 1.0e20;
        xL_H2 = std::max(-xL_H2_max, std::min(xL_H2, xL_H2_max));
        // dyHp, dyHep, dyHepp stay as dy[3], dy[8], dy[9]
    }

    xLmbdch = ((2.18e-11*dyHp + 3.94e-11*dyHep + 12.66e-11*dyHepp) / dt
               + xL_H2) / ((1.0 + 4.0*yHe) * xm_p);
}


// ─────────────────────────────────────────────────────────────────────────────
// chemcool<89, 1200>
//
// Secant method for temperature self-consistent integration.
// Differs from generic template only in the early-exit condition:
//   generic:      skip secant when  xnH  > 1e16  AND  T <= 1650 K
//   metal_grain:  skip secant when  xnH <= 1e16  AND  T <= 1650 K
//
// Port of chemcool() in chemistry_metal.f (line 26).
// ─────────────────────────────────────────────────────────────────────────────
template<>
inline void chemcool<89, 1200>(
    double xnH, double Tp,
    std::array<double, 89>& y,
    double dt,
    double& xmu, double& gamma, double& xLmbdch,
    std::array<double, 2400>& var,
    const ReactionTable<89, 1200>& tbl,
    const ChemParams& params)
{
    constexpr double xm_p  = phys::xm_p;  // 1.67262e-24 g (high-precision; unified in Phase 5)
    constexpr double xk_B  = phys::xk_B;  // 1.380662e-16 erg/K (high-precision; unified in Phase 5)
    constexpr int    maxit  = 100;

    // Residual using updated gamma/xmu from chemreact output
    auto func = [&](double t, double g, double mu, double xLc) -> double {
        return Tp - t - (g - 1.0) * mu * xm_p * xLc * dt / xk_B;
    };

    // First evaluation at Tp1 = Tp
    double Tp1 = Tp;
    std::array<double, 89> ytmp = y;
    double xmu1 = xmu, gamma1 = gamma, xLmbdch1 = 0.0;
    chemreact<89, 1200>(xnH, Tp1, ytmp, dt, xmu1, gamma1, xLmbdch1, var, tbl, params);
    double fL = func(Tp1, gamma1, xmu1, xLmbdch1);

    // Skip secant when xnH <= 1e16 AND Tp <= 1650 K
    // (metal_grain condition; inverse of zero_metal which skips at xnH > 1e16)
    if (xnH <= 1.0e16 && Tp <= 1.65e3) {
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
    chemreact<89, 1200>(xnH, Tp2, ytmp, dt, xmu2, gamma2, xLmbdch2, var, tbl, params);
    double f = func(Tp2, gamma2, xmu2, xLmbdch2);

    // Ensure |fL| >= |f|
    double TpL, Tpsec;
    if (std::abs(fL) < std::abs(f)) {
        TpL = Tp2; Tpsec = Tp1;
        std::swap(fL, f);
    } else {
        TpL = Tp1; Tpsec = Tp2;
    }

    // Secant iterations
    for (int i = 0; i < maxit; ++i) {
        // Guard against f ≈ fL which causes division by zero → NaN
        double df = f - fL;
        if (std::abs(df) < 1.0e-30 * (std::abs(f) + std::abs(fL) + 1.0e-30))
            break;
        double dTp = (TpL - Tpsec) * f / df;
        TpL   = Tpsec;
        fL    = f;
        Tpsec += dTp;

        ytmp = y;
        chemreact<89, 1200>(xnH, Tpsec, ytmp, dt,
                            xmu1, gamma1, xLmbdch1, var, tbl, params);
        f = func(Tpsec, gamma1, xmu1, xLmbdch1);

        if (std::abs(f / Tpsec) <= 1.0e-7) break;
    }

    y       = ytmp;
    xmu     = xmu1;
    gamma   = gamma1;
    xLmbdch = xLmbdch1;
}

} // namespace chemistry
