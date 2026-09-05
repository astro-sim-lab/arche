// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
// hdf5_utils.h — thin inline wrappers around the HDF5 C API
//
// Shared by primordial_collapse and metal_grain_collapse.
// Include this header and bring names into scope with:
//
//   using h5utils::H5Write1dInt;
//   using h5utils::H5Write1d;
//   using h5utils::H5Write2d;
//   using h5utils::H5WriteStrAttr;
//   using h5utils::H5WriteDblAttr;
//   using h5utils::H5Create;
//
#pragma once

#include <hdf5.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace h5utils {

// ─── RAII scope guards ───────────────────────────────────────────────────────
// Close an HDF5 handle automatically at scope exit, so an exception thrown
// mid-load (e.g. a dimension-mismatch error) cannot leak the file/group.
// Non-copyable; implicitly convertible to hid_t for use with the HDF5 C API.

struct H5FileGuard {
  hid_t id;
  explicit H5FileGuard(hid_t h) : id(h) {}
  ~H5FileGuard() {
    if (id >= 0) H5Fclose(id);
  }
  H5FileGuard(const H5FileGuard&) = delete;
  H5FileGuard& operator=(const H5FileGuard&) = delete;
  operator hid_t() const { return id; }
};

struct H5GroupGuard {
  hid_t id;
  explicit H5GroupGuard(hid_t h) : id(h) {}
  ~H5GroupGuard() {
    if (id >= 0) H5Gclose(id);
  }
  H5GroupGuard(const H5GroupGuard&) = delete;
  H5GroupGuard& operator=(const H5GroupGuard&) = delete;
  operator hid_t() const { return id; }
};

inline void H5Write1dInt(hid_t loc, const std::string& name,
                         const std::vector<int>& v) {
  hsize_t n = static_cast<hsize_t>(v.size());
  hid_t sp = H5Screate_simple(1, &n, nullptr);
  hid_t ds = H5Dcreate2(loc, name.c_str(), H5T_NATIVE_INT, sp, H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
  H5Dclose(ds);
  H5Sclose(sp);
}

inline void H5Write1d(hid_t loc, const std::string& name,
                      const std::vector<double>& v) {
  hsize_t n = static_cast<hsize_t>(v.size());
  hid_t sp = H5Screate_simple(1, &n, nullptr);
  hid_t ds = H5Dcreate2(loc, name.c_str(), H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
  H5Dclose(ds);
  H5Sclose(sp);
}

// 2-D dataset: row-major layout [nrows × ncols]
inline void H5Write2d(hid_t loc, const std::string& name,
                      const std::vector<double>& data, hsize_t nrows,
                      hsize_t ncols) {
  hsize_t dims[2] = {nrows, ncols};
  hid_t sp = H5Screate_simple(2, dims, nullptr);
  hid_t ds = H5Dcreate2(loc, name.c_str(), H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
  H5Dclose(ds);
  H5Sclose(sp);
}

// Scalar fixed-length string attribute on any HDF5 object
// Returns false if the attribute could not be created or written (a full disk,
// a broken file handle, a duplicate name).  HDF5's default handler prints a
// diagnostic to stderr, but that is not something a caller can act on, so the
// status is returned as well.  Existing callers may ignore it.
inline bool H5WriteStrAttr(hid_t obj, const std::string& name,
                           const std::string& val) {
  hid_t sp = H5Screate(H5S_SCALAR);
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, val.size() + 1);
  hid_t at = H5Acreate2(obj, name.c_str(), type, sp, H5P_DEFAULT, H5P_DEFAULT);
  bool ok = (at >= 0);
  if (ok) {
    ok = (H5Awrite(at, type, val.c_str()) >= 0);
    H5Aclose(at);
  }
  H5Tclose(type);
  H5Sclose(sp);
  return ok;
}

// Scalar double attribute
// Returns false on failure; see H5WriteStrAttr.
inline bool H5WriteDblAttr(hid_t obj, const std::string& name, double val) {
  hid_t sp = H5Screate(H5S_SCALAR);
  hid_t at = H5Acreate2(obj, name.c_str(), H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT,
                        H5P_DEFAULT);
  bool ok = (at >= 0);
  if (ok) {
    ok = (H5Awrite(at, H5T_NATIVE_DOUBLE, &val) >= 0);
    H5Aclose(at);
  }
  H5Sclose(sp);
  return ok;
}

// Scalar int attribute
// Returns false on failure; see H5WriteStrAttr.
inline bool H5WriteIntAttr(hid_t obj, const std::string& name, int val) {
  hid_t sp = H5Screate(H5S_SCALAR);
  hid_t at = H5Acreate2(obj, name.c_str(), H5T_NATIVE_INT, sp, H5P_DEFAULT,
                        H5P_DEFAULT);
  bool ok = (at >= 0);
  if (ok) {
    ok = (H5Awrite(at, H5T_NATIVE_INT, &val) >= 0);
    H5Aclose(at);
  }
  H5Sclose(sp);
  return ok;
}

// Open (create/truncate) HDF5 file with POSIX locking disabled
// (required on WSL / some NFS mounts)
inline hid_t H5Create(const std::string& path) {
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_file_locking(fapl, /*use_file_locking=*/0, /*ignore_when_disabled=*/1);
  hid_t fid = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  H5Pclose(fapl);
  return fid;
}

// ─── Read helpers ──────────────────────────────────────────────────────────

// Open an existing HDF5 file read-only with POSIX locking disabled.
inline hid_t H5OpenRead(const std::string& path) {
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_file_locking(fapl, /*use_file_locking=*/0, /*ignore_when_disabled=*/1);
  hid_t fid = H5Fopen(path.c_str(), H5F_ACC_RDONLY, fapl);
  H5Pclose(fapl);
  return fid;
}

inline hsize_t H5DatasetSize(hid_t loc, const std::string& name) {
  hid_t ds = H5Dopen2(loc, name.c_str(), H5P_DEFAULT);
  if (ds < 0)
    throw std::runtime_error("HDF5: cannot open dataset '" + name + "'");
  hid_t sp = H5Dget_space(ds);
  hsize_t n = 0;
  H5Sget_simple_extent_dims(sp, &n, nullptr);
  H5Sclose(sp);
  H5Dclose(ds);
  return n;
}

inline std::vector<int> H5Read1dInt(hid_t loc, const std::string& name) {
  hid_t ds = H5Dopen2(loc, name.c_str(), H5P_DEFAULT);
  if (ds < 0)
    throw std::runtime_error("HDF5: cannot open dataset '" + name + "'");
  hid_t sp = H5Dget_space(ds);
  hsize_t n = 0;
  H5Sget_simple_extent_dims(sp, &n, nullptr);
  std::vector<int> v(static_cast<size_t>(n));
  H5Dread(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
  H5Sclose(sp);
  H5Dclose(ds);
  return v;
}

inline std::vector<double> H5Read1d(hid_t loc, const std::string& name) {
  hid_t ds = H5Dopen2(loc, name.c_str(), H5P_DEFAULT);
  if (ds < 0)
    throw std::runtime_error("HDF5: cannot open dataset '" + name + "'");
  hid_t sp = H5Dget_space(ds);
  hsize_t n = 0;
  H5Sget_simple_extent_dims(sp, &n, nullptr);
  std::vector<double> v(static_cast<size_t>(n));
  H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
  H5Sclose(sp);
  H5Dclose(ds);
  return v;
}

inline int H5ReadIntAttr(hid_t obj, const std::string& name) {
  hid_t at = H5Aopen(obj, name.c_str(), H5P_DEFAULT);
  if (at < 0)
    throw std::runtime_error("HDF5: cannot open attribute '" + name + "'");
  int val = 0;
  H5Aread(at, H5T_NATIVE_INT, &val);
  H5Aclose(at);
  return val;
}

inline bool H5GroupExists(hid_t loc, const std::string& name) {
  return H5Lexists(loc, name.c_str(), H5P_DEFAULT) > 0;
}

}  // namespace h5utils
