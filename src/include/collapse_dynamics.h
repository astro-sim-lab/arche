// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// collapse_dynamics.h  Eshared helpers for one-zone collapse dynamics
//
// Provides the gamma-dependent free-fall correction factor from:
//   Higuchi, Machida & Susa (2018) MNRAS 475, 3331  Eequations (5)-(7)
//   (original formulation: Omukai et al. 2005)

#include <cmath>
#include <algorithm>


// ─────────────────────────────────────────────────────────────────────────────
// fgamma_collapse — pressure-gradient correction factor f(γ)
//
// γ (≡ dlnP/dlnρ) is the effective ratio of specific heats.
// f represents the ratio of pressure gradient force to gravity at the cloud
// centre (Omukai et al. 2005).  Used in Higuchi+2018 Eq.(7).
//
// Returns f ∈ [0, ~1):  f = 0 → pure free-fall;  f → 1 → collapse halts.
// ─────────────────────────────────────────────────────────────────────────────
inline double fgamma_collapse(double gamma) {
    if (gamma < 0.83)
        return 0.0;
    else if (gamma < 1.0)
        return 0.6 + 2.5*(gamma - 1.0) - 6.0*(gamma - 1.0)*(gamma - 1.0);
    else
        return 1.0 + 0.2*(gamma - 4.0/3.0) - 2.9*(gamma - 4.0/3.0)*(gamma - 4.0/3.0);
}

// ─────────────────────────────────────────────────────────────────────────────
// t_eff_collapse — effective collapse timescale [s]
//
// Two modes selected by ff_gamma:
//
//   ff_gamma = false  (modes A and B — fret):
//     t_eff = f_ret * t_ff
//
//   ff_gamma = true   (mode C — gamma):
//     t_eff = t_ff / sqrt(1 - f(γ))     (Higuchi+2018 Eq.5-7)
//     where f = fgamma_collapse(gamma).
//     f_ret is unused in this mode.
//
// In both cases: drho/dt = rho / t_eff
// ─────────────────────────────────────────────────────────────────────────────
inline double t_eff_collapse(double t_ff, double f_ret,
                              double gamma, bool ff_gamma) {
    if (ff_gamma) {
        const double f     = fgamma_collapse(gamma);
        const double denom = std::sqrt(std::max(1.0 - f, 1.0e-10));
        return t_ff / denom;
    }
    return f_ret * t_ff;
}
