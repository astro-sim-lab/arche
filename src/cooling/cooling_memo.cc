// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE

// The line-cooling and continuum entry points that keep their last results per thread
// (core/exact_memo.h).  They are compiled here once: every caller in a program reaches
// the same compiled evaluation, so a kept result is exactly what that evaluation
// returns for the same input bits, whatever the callers' own compiler flags.

#include <array>

#include "cooling/cooling_grain.h"
#include "cooling/cooling_molecular.h"
#include "core/exact_memo.h"

namespace arche {

namespace {

// A molecular line's cooling (co_cool_eval, ...) is a pure function of its five
// arguments.
template <class Eval>
double memo_line(ExactMemo<5, 1, 8>& memo, double nH, double T_K, double y_H2,
                 double Nc, double tau_cnt, Eval eval) {
  const std::array<double, 5> in{nH, T_K, y_H2, Nc, tau_cnt};
  std::array<double, 1> out{};
  if (memo.find(in, out)) return out[0];
  out[0] = eval(nH, T_K, y_H2, Nc, tau_cnt);
  memo.insert(in, out);
  return out[0];
}

}  // namespace

double co_cool(double nH, double T_K, double y_H2, double Nc_CO,
               double tau_cnt) {
  thread_local ExactMemo<5, 1, 8> memo;
  return memo_line(memo, nH, T_K, y_H2, Nc_CO, tau_cnt, co_cool_eval);
}

double oh_cool(double nH, double T_K, double y_H2, double Nc_OH,
               double tau_cnt) {
  thread_local ExactMemo<5, 1, 8> memo;
  return memo_line(memo, nH, T_K, y_H2, Nc_OH, tau_cnt, oh_cool_eval);
}

double h2o_cool(double nH, double T_K, double y_H2, double Nc_H2O,
                double tau_cnt) {
  thread_local ExactMemo<5, 1, 8> memo;
  return memo_line(memo, nH, T_K, y_H2, Nc_H2O, tau_cnt, h2o_cool_eval);
}

double cii_cool_fs(double nH, double T_K, double y_e, double y_a, double y_m,
                   double tau_cnt, double Nc, double T_rad, double esc[1]) {
  return mol_detail::memo_esc_solve<1>(
      std::array<double, 8>{nH, T_K, y_e, y_a, y_m, tau_cnt, Nc, T_rad}, esc,
      [&] {
        return cii_cool_fs_solve(nH, T_K, y_e, y_a, y_m, tau_cnt, Nc, T_rad,
                                 esc);
      });
}

double cii_cool_meta(double nH, double T_K, double y_e, double tau_cnt,
                     double Nc, double T_rad, double esc[1]) {
  return mol_detail::memo_esc_solve<1>(
      std::array<double, 6>{nH, T_K, y_e, tau_cnt, Nc, T_rad}, esc, [&] {
        return cii_cool_meta_solve(nH, T_K, y_e, tau_cnt, Nc, T_rad, esc);
      });
}

double ci_cool_fs(double nH, double T_K, double y_e, double y_a, double y_m,
                  double tau_cnt, double Nc, double T_rad, double esc[3]) {
  return mol_detail::memo_esc_solve<3>(
      std::array<double, 8>{nH, T_K, y_e, y_a, y_m, tau_cnt, Nc, T_rad}, esc,
      [&] {
        return ci_cool_fs_solve(nH, T_K, y_e, y_a, y_m, tau_cnt, Nc, T_rad,
                                esc);
      });
}

double oi_cool_fs(double nH, double T_K, double y_e, double y_a, double y_m,
                  double tau_cnt, double Nc, double T_rad, double esc[3]) {
  return mol_detail::memo_esc_solve<3>(
      std::array<double, 8>{nH, T_K, y_e, y_a, y_m, tau_cnt, Nc, T_rad}, esc,
      [&] {
        return oi_cool_fs_solve(nH, T_K, y_e, y_a, y_m, tau_cnt, Nc, T_rad,
                                esc);
      });
}

void cnt_cool_metal_consistent(double nH, double T_K, double T_rad,
                               double column, double Z_metal, double& T_gr_K,
                               double& k_gr, double& k_gas, double& Lambda_gr,
                               double& Lambda_gas, double& Lambda_cnt,
                               double& tau_cnt, double& esc_cnt) {
  thread_local ExactMemo<6, 8, 8> memo;
  const std::array<double, 6> in{nH, T_K, T_rad, column, Z_metal, T_gr_K};
  std::array<double, 8> out{};
  if (!memo.find(in, out)) {
    cnt_cool_metal_consistent_solve(nH, T_K, T_rad, column, Z_metal, T_gr_K,
                                    k_gr, k_gas, Lambda_gr, Lambda_gas,
                                    Lambda_cnt, tau_cnt, esc_cnt);
    memo.insert(in, {T_gr_K, k_gr, k_gas, Lambda_gr, Lambda_gas, Lambda_cnt,
                     tau_cnt, esc_cnt});
    return;
  }
  T_gr_K = out[0];
  k_gr = out[1];
  k_gas = out[2];
  Lambda_gr = out[3];
  Lambda_gas = out[4];
  Lambda_cnt = out[5];
  tau_cnt = out[6];
  esc_cnt = out[7];
}

}  // namespace arche
