// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
//
// main.cc — magnetic Reynolds number over a one-zone collapse track
//
// Reads the output of a collapse application (metal_collapse,
// arche_collapse_metal_minimal, prim_collapse, arche_collapse_prim_minimal)
// and, for every output row and every field strength B of a logarithmic grid,
// evaluates the Ohmic and ambipolar diffusivities of the row's composition
// (arche::model_resistivity) and the magnetic Reynolds number at the Jeans
// scale (Nakauchi, Omukai & Susa 2021, eq. 41 with v = v_ff, L_B = lambda_J):
//
//   Rm(lambda_J) = v_ff lambda_J / (eta_O + eta_A),  v_ff = lambda_J / (3 t_ff).
//
// The collapse is computed without magnetic fields, so one track serves every
// B.  The y columns of the input are in the model's species order, as the
// collapse applications write them; the column count must equal the model's
// species count.
//
// Usage:
//   magnetic_reynolds <collapse.h5> <output.h5>
//
// Environment variables:
//   MAGRE_MODEL       — registered model name (default: Nakauchi2021)
//   MAGRE_CHEM_TABLE  — reaction table of that model
//                       (default: compile-time MAGRE_CHEM_TABLE)
//   MAGRE_LOG10B_MIN  — lowest log10(B/G) of the grid (default: -10)
//   MAGRE_LOG10B_MAX  — highest log10(B/G) of the grid (default: 10)
//   MAGRE_NB          — number of grid points, >= 2 (default: 401)
//
// Output datasets (N rows of the input, N_B field strengths):
//   nH, T_K, rho, lambda_J, t_ff, B_cr  (N)       copied from the input
//   log10_B                             (N_B)
//   eta_O                               (N)       [cm^2 s^-1], B-independent
//   eta_A, Rm                           (N, N_B)  [cm^2 s^-1], dimensionless

#include <hdf5.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

#include "api/arche_api.h"
#include "core/hdf5_utils.h"

#ifndef MAGRE_CHEM_TABLE
#define MAGRE_CHEM_TABLE "data/metal_grain.h5"
#endif

namespace {

using h5utils::H5Create;
using h5utils::H5DatasetSize;
using h5utils::H5FileGuard;
using h5utils::H5OpenRead;
using h5utils::H5Read1d;
using h5utils::H5Write1d;
using h5utils::H5Write2d;
using h5utils::H5WriteIntAttr;
using h5utils::H5WriteStrAttr;

std::string env_str(const char* name, const char* def) {
  const char* v = std::getenv(name);
  return (v && *v) ? std::string(v) : std::string(def);
}

double env_double(const char* name, double def) {
  const char* v = std::getenv(name);
  if (!v || !*v) return def;
  char* end = nullptr;
  const double x = std::strtod(v, &end);
  if (end == v || *end != '\0' || !std::isfinite(x))
    throw std::invalid_argument(std::string(name) + ": not a finite number");
  return x;
}

// Row-major (rows x cols) double dataset.
std::vector<double> read_2d(hid_t loc, const std::string& name, hsize_t& rows,
                            hsize_t& cols) {
  hid_t ds = H5Dopen2(loc, name.c_str(), H5P_DEFAULT);
  if (ds < 0)
    throw std::runtime_error("HDF5: cannot open dataset '" + name + "'");
  hid_t sp = H5Dget_space(ds);
  hsize_t dims[2] = {0, 0};
  const int rank = H5Sget_simple_extent_ndims(sp);
  if (rank != 2) {
    H5Sclose(sp);
    H5Dclose(ds);
    throw std::runtime_error("HDF5: dataset '" + name + "' is not 2-D");
  }
  H5Sget_simple_extent_dims(sp, dims, nullptr);
  rows = dims[0];
  cols = dims[1];
  std::vector<double> v(static_cast<size_t>(rows * cols));
  const herr_t st =
      H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
  H5Sclose(sp);
  H5Dclose(ds);
  if (st < 0) throw std::runtime_error("HDF5: cannot read '" + name + "'");
  return v;
}

int run(const std::string& in_path, const std::string& out_path) {
  const std::string model_name = env_str("MAGRE_MODEL", "Nakauchi2021");
  const std::string table = env_str("MAGRE_CHEM_TABLE", MAGRE_CHEM_TABLE);
  const double lb_min = env_double("MAGRE_LOG10B_MIN", -10.0);
  const double lb_max = env_double("MAGRE_LOG10B_MAX", 10.0);
  const double nb_d = env_double("MAGRE_NB", 401.0);
  if (!(nb_d >= 2.0) || nb_d != std::floor(nb_d) || nb_d > 1.0e6)
    throw std::invalid_argument("MAGRE_NB: an integer in [2, 1e6]");
  if (!(lb_max > lb_min))
    throw std::invalid_argument(
        "MAGRE_LOG10B_MAX must exceed MAGRE_LOG10B_MIN");
  const int nb = static_cast<int>(nb_d);

  // ── input track ──────────────────────────────────────────────────────────
  hsize_t n = 0, nsp = 0;
  std::vector<double> y, nH, T_K, rho, lambda_J, t_ff, B_cr;
  {
    H5FileGuard fin(H5OpenRead(in_path));
    if (fin.id < 0) throw std::runtime_error("cannot open " + in_path);
    y = read_2d(fin, "y", n, nsp);
    nH = H5Read1d(fin, "nH");
    T_K = H5Read1d(fin, "T_K");
    rho = H5Read1d(fin, "rho");
    lambda_J = H5Read1d(fin, "lambda_J");
    t_ff = H5Read1d(fin, "t_ff");
    B_cr = H5Read1d(fin, "B_cr");
    for (const char* name : {"nH", "T_K", "rho", "lambda_J", "t_ff", "B_cr"})
      if (H5DatasetSize(fin, name) != n)
        throw std::runtime_error(std::string("dataset '") + name +
                                 "' length differs from the rows of y");
  }

  arche::ModelRuntimePtr model = arche::model_create_owned(model_name, table);
  const int n_model = arche::model_n_species(*model);
  if (static_cast<hsize_t>(n_model) != nsp)
    throw std::runtime_error(in_path + ": y has " + std::to_string(nsp) +
                             " columns, model " + model_name + " has " +
                             std::to_string(n_model) + " species");
  arche::ModelCellPtr cell = arche::model_cell_create_owned(*model);

  // ── grid ─────────────────────────────────────────────────────────────────
  std::vector<double> log10_B(nb);
  for (int k = 0; k < nb; ++k)
    log10_B[k] = lb_min + (lb_max - lb_min) * k / (nb - 1);

  std::vector<double> eta_O(n), eta_A(n * nb), Rm(n * nb);
  for (hsize_t i = 0; i < n; ++i) {
    double* yc = arche::model_cell_y(*cell);
    for (hsize_t j = 0; j < nsp; ++j) yc[j] = y[i * nsp + j];
    arche::model_cell_nH(*cell) = nH[i];
    arche::model_cell_T_K(*cell) = T_K[i];
    // v_ff lambda_J = lambda_J^2 / (3 t_ff)
    const double vl = lambda_J[i] * lambda_J[i] / (3.0 * t_ff[i]);
    for (int k = 0; k < nb; ++k) {
      const arche::MagneticDiffusivities d =
          arche::model_resistivity(*cell, T_K[i], std::pow(10.0, log10_B[k]));
      if (k == 0) eta_O[i] = d.eta_O;
      eta_A[i * nb + k] = d.eta_A;
      Rm[i * nb + k] = vl / (d.eta_O + d.eta_A);
    }
  }

  // ── output ───────────────────────────────────────────────────────────────
  H5FileGuard fout(H5Create(out_path));
  if (fout.id < 0) throw std::runtime_error("cannot create " + out_path);
  H5WriteStrAttr(fout, "description",
                 "Magnetic Reynolds number Rm(lambda_J) = lambda_J^2 / "
                 "(3 t_ff (eta_O + eta_A)) over a collapse track");
  H5WriteStrAttr(fout, "source", in_path);
  H5WriteStrAttr(fout, "model", model_name);
  H5WriteIntAttr(fout, "n_B", nb);
  H5WriteStrAttr(fout, "units_eta", "cm^2 s^-1");
  H5WriteStrAttr(fout, "units_B", "G");
  H5Write1d(fout, "nH", nH);
  H5Write1d(fout, "T_K", T_K);
  H5Write1d(fout, "rho", rho);
  H5Write1d(fout, "lambda_J", lambda_J);
  H5Write1d(fout, "t_ff", t_ff);
  H5Write1d(fout, "B_cr", B_cr);
  H5Write1d(fout, "log10_B", log10_B);
  H5Write1d(fout, "eta_O", eta_O);
  H5Write2d(fout, "eta_A", eta_A, n, static_cast<hsize_t>(nb));
  H5Write2d(fout, "Rm", Rm, n, static_cast<hsize_t>(nb));
  std::printf("%s -> %s  (%llu rows x %d field strengths, model %s)\n",
              in_path.c_str(), out_path.c_str(),
              static_cast<unsigned long long>(n), nb, model_name.c_str());
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::fprintf(stderr, "usage: %s <collapse.h5> <output.h5>\n", argv[0]);
    return 2;
  }
  try {
    return run(argv[1], argv[2]);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "magnetic_reynolds: %s\n", e.what());
    return 1;
  }
}
