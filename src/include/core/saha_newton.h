// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// saha_newton.h — the Newton iteration shared by every Saha-equilibrium solver
//
// Each Saha solver reduces the high-density equilibrium to a few unknown free
// abundances -- (y_e, y_H) for the primordial models, (y_e, y_H, y_C, y_O) for
// metal_grain, where C and O also bind to each other and to H -- and to as many
// balance residuals (charge, H, C, O). The unknowns differ with the chemistry;
// the iteration must not.
//
// A plain step `x += d` from a state far from equilibrium can land on the root with
// y_e < 0 that the charge balance also has, and rounding decides whether it does. A
// floored step `max(x + d, 1e-30)` can pin y_C or y_O at the floor; the iteration then
// never leaves it and the state violates the C or O budget (by up to 20 orders for the
// metal_grain full network below ~2800 K). Neither detects the failure unless the
// convergence is checked.
//
// saha_newton_log iterates on u = ln y, which keeps every unknown positive with no
// floor, limits a step to |du| <= 2, halves a step whose residual is not finite,
// and reports whether it converged so the caller can try another start or fail.
// ---------------------------------------------------------------------------
#include <algorithm>
#include <array>
#include <cmath>

namespace arche {

// Solve A d = b in place for an m x m system (m <= N) held in the leading block
// of a row-major N x N array, with partial pivoting. False if singular.
template <int N>
bool saha_lu_solve(std::array<double, N * N>& A, std::array<double, N>& b,
                   int m) {
  for (int k = 0; k < m; ++k) {
    int piv = k;
    double big = std::abs(A[k * N + k]);
    for (int i = k + 1; i < m; ++i) {
      const double v = std::abs(A[i * N + k]);
      if (v > big) {
        big = v;
        piv = i;
      }
    }
    if (!(big > 0.0) || !std::isfinite(big)) return false;
    if (piv != k) {
      for (int j = 0; j < m; ++j) std::swap(A[k * N + j], A[piv * N + j]);
      std::swap(b[k], b[piv]);
    }
    for (int i = k + 1; i < m; ++i) {
      const double l = A[i * N + k] / A[k * N + k];
      for (int j = k; j < m; ++j) A[i * N + j] -= l * A[k * N + j];
      b[i] -= l * b[k];
    }
  }
  for (int k = m - 1; k >= 0; --k) {
    double s = b[k];
    for (int j = k + 1; j < m; ++j) s -= A[k * N + j] * b[j];
    b[k] = s / A[k * N + k];
    if (!std::isfinite(b[k])) return false;
  }
  return true;
}

// x: the start on entry (active entries must be positive and finite), the root on
// a true return (unspecified otherwise). Inactive entries are held fixed.
// residual(y, f): the balances at the linear abundances y (N of each; the rows of
// inactive unknowns are ignored). Convergence: sum over active j of |f_j / y_j| <=
// tol, or a step that moves no active
// unknown by more than 1e-12 relative, which is where an iteration sits once a
// balance's rounding floor exceeds tol * y_j (an H balance with y_H << 1).
template <int N, class Residual>
bool saha_newton_log(std::array<double, N>& x,
                     const std::array<bool, N>& active, Residual&& residual,
                     int max_iter = 100, double tol = 1.0e-10) {
  constexpr double kMaxStep = 2.0;    // max |d ln y| per iteration
  constexpr double kFdStep = 1.0e-7;  // forward-difference step in ln y
  constexpr double kSettled = 1.0e-12;
  std::array<int, N> idx{};  // active unknowns, in order
  int m = 0;
  for (int j = 0; j < N; ++j) {
    if (!active[j]) continue;
    if (!(x[j] > 0.0) || !std::isfinite(x[j])) return false;
    idx[m++] = j;
  }
  if (m == 0) return true;

  std::array<double, N> u{};
  for (int j = 0; j < N; ++j) u[j] = active[j] ? std::log(x[j]) : x[j];
  auto to_y = [&](const std::array<double, N>& uu) {
    std::array<double, N> y{};
    for (int j = 0; j < N; ++j) y[j] = active[j] ? std::exp(uu[j]) : uu[j];
    return y;
  };
  auto eval = [&](const std::array<double, N>& uu, std::array<double, N>& f) {
    residual(to_y(uu), f);
    for (int k = 0; k < m; ++k)
      if (!std::isfinite(f[idx[k]])) return false;
    return true;
  };

  std::array<double, N> f{};
  if (!eval(u, f)) return false;
  for (int itr = 0; itr < max_iter; ++itr) {
    const std::array<double, N> y = to_y(u);
    double conv = 0.0;
    for (int k = 0; k < m; ++k) conv += std::abs(f[idx[k]] / y[idx[k]]);
    if (conv <= tol) {
      x = y;
      return true;
    }
    // Jacobian of the active rows with respect to the active ln y, column by column.
    std::array<double, N * N> J{};
    for (int c = 0; c < m; ++c) {
      std::array<double, N> up = u, fp{};
      up[idx[c]] += kFdStep;
      if (!eval(up, fp)) return false;
      for (int r = 0; r < m; ++r)
        J[r * N + c] = (fp[idx[r]] - f[idx[r]]) / kFdStep;
    }
    std::array<double, N> d{};
    for (int r = 0; r < m; ++r) d[r] = -f[idx[r]];
    if (!saha_lu_solve<N>(J, d, m)) return false;
    double big = 0.0;
    for (int k = 0; k < m; ++k) big = std::max(big, std::abs(d[k]));
    if (big > kMaxStep)
      for (int k = 0; k < m; ++k) d[k] *= kMaxStep / big;

    std::array<double, N> un = u, fn{};
    bool ok = false;
    for (int h = 0; h < 40 && !ok; ++h) {
      for (int k = 0; k < m; ++k) un[idx[k]] = u[idx[k]] + d[k];
      ok = eval(un, fn);
      if (!ok)
        for (int k = 0; k < m; ++k) d[k] *= 0.5;
    }
    if (!ok) return false;
    u = un;
    f = fn;
    double moved = 0.0;
    for (int k = 0; k < m; ++k) moved = std::max(moved, std::abs(d[k]));
    if (moved <= kSettled) {
      x = to_y(u);
      return true;
    }
  }
  return false;
}

}  // namespace arche
