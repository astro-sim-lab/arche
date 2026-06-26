// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// newton.h — model-independent numerical solvers.
//
// A dense Gauss-Jordan linear solve, the linear-solve policies, and a
// dimension-generic Newton-Raphson driver.  These carry no chemistry or model
// knowledge, so they live in the foundation layer: both the implicit chemistry
// integrator (solve/chemreact.h) and the per-model Saha-equilibrium solvers
// (models/*/equilibrium.h) depend on them without depending on each other.
//
// Optional Eigen support: define ARCHE_USE_EIGEN before including, or add
// -DARCHE_USE_EIGEN to the compiler flags.
// ---------------------------------------------------------------------------
#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

#include "core/state.h"  // numerics::eps_gaussj

#ifdef ARCHE_USE_EIGEN
#include <Eigen/Dense>
#endif

namespace arche {

// ─────────────────────────────────────────────────────────────────────────────
// gaussj_solve<N> — solve A·x = b in-place; b receives the solution.
//   a[N×N] row-major (destroyed), b[N].
//   Gauss-Jordan elimination (Numerical Recipes gaussj), C++ row-major.
//   Returns true on success; false if A is singular or the Eigen solve is
//   non-finite (b is zeroed in that case to prevent NaN propagation).
// ─────────────────────────────────────────────────────────────────────────────
template <int N>
bool gaussj_solve(std::array<double, N * N>& a, std::array<double, N>& b) {
#ifdef ARCHE_USE_EIGEN
  Eigen::Map<Eigen::Matrix<double, N, N, Eigen::RowMajor>> A(a.data());
  Eigen::Map<Eigen::Matrix<double, N, 1>> x(b.data());
  x = A.partialPivLu().solve(x);
  // Guard: near-singular A produces Inf/NaN in x; zero out to prevent
  // NaN propagation through the chemistry network.
  if (!x.allFinite()) {
    x.setZero();
    return false;
  }
  return true;
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
            double aabs = std::abs(a[j * N + k]);
            if (aabs >= big) {
              big = aabs;
              irow = j;
              icol = k;
            }
          } else if (ipiv[k] > 1) {
            singular = true;
            break;
          }
        }
      }
      if (singular) break;
    }
    if (singular) break;

    ++ipiv[icol];
    if (irow != icol) {
      for (int l = 0; l < N; ++l) std::swap(a[irow * N + l], a[icol * N + l]);
      std::swap(b[irow], b[icol]);
    }
    indxr[i] = irow;
    indxc[i] = icol;

    if (a[icol * N + icol] == 0.0) {
      // Pivot is zero (matrix is singular): return zero update to prevent
      // NaN/Inf propagation from a partial-elimination state in b[].
      b.fill(0.0);
      return false;
    }

    double pivinv = 1.0 / a[icol * N + icol];
    a[icol * N + icol] = 1.0;
    for (int l = 0; l < N; ++l) a[icol * N + l] *= pivinv;
    b[icol] *= pivinv;

    for (int ll = 0; ll < N; ++ll) {
      if (ll != icol) {
        double dum = a[ll * N + icol];
        a[ll * N + icol] = 0.0;
        for (int l = 0; l < N; ++l) a[ll * N + l] -= a[icol * N + l] * dum;
        double ref_b = std::abs(b[icol] * dum);
        b_max[ll] = std::max(b_max[ll], ref_b);
        b[ll] -= b[icol] * dum;
        if (ref_b != 0.0 && std::abs(b[ll]) / ref_b < eps) b[ll] = 0.0;
      }
    }
  }

  if (singular) {
    // ipiv[k] > 1 path: pivot-tracking logic error.  Return zero update.
    b.fill(0.0);
    return false;
  }

  // Flush residual noise in solution vector
  for (int ll = 0; ll < N; ++ll)
    if (b_max[ll] != 0.0 && std::abs(b[ll] / b_max[ll]) <= eps) b[ll] = 0.0;

  // Undo column permutation
  for (int l = N - 1; l >= 0; --l)
    if (indxr[l] != indxc[l]) {
      for (int k = 0; k < N; ++k)
        std::swap(a[k * N + indxr[l]], a[k * N + indxc[l]]);
    }
  return true;
#endif
}

// ─────────────────────────────────────────────────────────────────────────────
// newton_solve<N, LinSolve> — dimension-generic Newton–Raphson driver.
//
//   Solves F(x) = 0 for x ∈ R^N with a forward-difference Jacobian.  Each
//   iteration evaluates the residual, tests the weighted convergence measure
//   Σ|f_j / x_j| over the active unknowns, assembles the Jacobian column by
//   column by perturbing x_j → x_j(1+eps), solves J·δ = −F through the LinSolve
//   policy, and applies the caller's (optionally clamped) update.  The
//   Saha-equilibrium solvers equichem / equichem_minimal (2D) and
//   equichem_metal (4D) share this driver: the only per-call variation is the
//   residual, the active mask, the update rule, and the linear-solve policy, so
//   the arithmetic of every common step is identical across callers.
//
//   active[j]  selects whether unknown j is varied — its column is perturbed,
//   residual j enters the convergence sum, and the update may touch it; an
//   inactive unknown contributes a zero Jacobian column and is left to the
//   update rule.
//
//   LinSolve policy:  static void solve(std::array<double,N*N>& J,
//                                       std::array<double,N>& b)
//   with b entering as the residual F and leaving as δ solving J·δ = −F.
//   CramerSolve2 preserves the 2×2 Cramer evaluation order; GaussjLinSolve<N>
//   defers to gaussj_solve<N> (Eigen partialPivLu when ARCHE_USE_EIGEN).
// ─────────────────────────────────────────────────────────────────────────────
struct NewtonOpts {
  double eps;    // forward-difference step (relative)
  int max_iter;  // iteration cap
  double tol;    // convergence threshold on Σ_active |f_j / x_j|
};

// 2×2 Cramer's rule.  J row-major [a00 a01 a10 a11]; b = (f0, f1) → δ = −J⁻¹ f.
struct CramerSolve2 {
  static void solve(std::array<double, 4>& J, std::array<double, 2>& b) {
    double det = J[0] * J[3] - J[1] * J[2];
    double f0 = b[0], f1 = b[1];
    b[0] = -(J[3] * f0 - J[1] * f1) / det;
    b[1] = (J[2] * f0 - J[0] * f1) / det;
  }
};

// N×N dense solve via gaussj_solve.
template <int N>
struct GaussjLinSolve {
  static void solve(std::array<double, N * N>& J, std::array<double, N>& b) {
    for (int i = 0; i < N; ++i) b[i] = -b[i];  // rhs = −F
    gaussj_solve<N>(J, b);                     // b ← δ
  }
};

template <int N, class LinSolve, class Residual, class Update>
void newton_solve(std::array<double, N>& x, const std::array<bool, N>& active,
                  Residual&& residual, Update&& update, const NewtonOpts& opt) {
  std::array<double, N> f{};      // residual F(x)
  std::array<double, N> fp{};     // residual at the perturbed point
  std::array<double, N * N> J{};  // Jacobian, row-major (residual-major)

  for (int itr = 0; itr < opt.max_iter; ++itr) {
    residual(x, f);

    double conv = 0.0;
    for (int j = 0; j < N; ++j)
      if (active[j]) conv += std::abs(f[j] / x[j]);
    if (conv <= opt.tol) break;

    // Forward-difference Jacobian, one unknown (column) at a time.
    for (int j = 0; j < N; ++j) {
      if (active[j]) {
        std::array<double, N> xp = x;
        xp[j] = x[j] * (1.0 + opt.eps);
        residual(xp, fp);
        double h = opt.eps * x[j];
        for (int i = 0; i < N; ++i) J[i * N + j] = (fp[i] - f[i]) / h;
      } else {
        for (int i = 0; i < N; ++i) J[i * N + j] = 0.0;
      }
    }

    std::array<double, N> delta = f;  // residual; LinSolve overwrites with δ
    LinSolve::solve(J, delta);
    update(x, delta);
  }
}

}  // namespace arche
