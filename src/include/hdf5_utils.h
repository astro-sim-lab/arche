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

#include <string>
#include <string_view>
#include <vector>

#include <hdf5.h>

namespace h5utils {

inline void H5Write1dInt(hid_t loc, std::string_view name,
                         const std::vector<int>& v)
{
    hsize_t n = static_cast<hsize_t>(v.size());
    hid_t sp  = H5Screate_simple(1, &n, nullptr);
    hid_t ds  = H5Dcreate2(loc, name.data(), H5T_NATIVE_INT, sp,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
    H5Dclose(ds);
    H5Sclose(sp);
}

inline void H5Write1d(hid_t loc, std::string_view name,
                      const std::vector<double>& v)
{
    hsize_t n = static_cast<hsize_t>(v.size());
    hid_t sp  = H5Screate_simple(1, &n, nullptr);
    hid_t ds  = H5Dcreate2(loc, name.data(), H5T_NATIVE_DOUBLE, sp,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
    H5Dclose(ds);
    H5Sclose(sp);
}

// 2-D dataset: row-major layout [nrows × ncols]
inline void H5Write2d(hid_t loc, std::string_view name,
                      const std::vector<double>& data,
                      hsize_t nrows, hsize_t ncols)
{
    hsize_t dims[2] = {nrows, ncols};
    hid_t sp = H5Screate_simple(2, dims, nullptr);
    hid_t ds = H5Dcreate2(loc, name.data(), H5T_NATIVE_DOUBLE, sp,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    H5Dclose(ds);
    H5Sclose(sp);
}

// Scalar fixed-length string attribute on any HDF5 object
inline void H5WriteStrAttr(hid_t obj, std::string_view name, const std::string& val)
{
    hid_t sp   = H5Screate(H5S_SCALAR);
    hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, val.size() + 1);
    hid_t at   = H5Acreate2(obj, name.data(), type, sp, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(at, type, val.c_str());
    H5Aclose(at);
    H5Tclose(type);
    H5Sclose(sp);
}

// Scalar double attribute
inline void H5WriteDblAttr(hid_t obj, std::string_view name, double val)
{
    hid_t sp = H5Screate(H5S_SCALAR);
    hid_t at = H5Acreate2(obj, name.data(), H5T_NATIVE_DOUBLE, sp,
                           H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(at, H5T_NATIVE_DOUBLE, &val);
    H5Aclose(at);
    H5Sclose(sp);
}

// Open (create/truncate) HDF5 file with POSIX locking disabled
// (required on WSL / some NFS mounts)
inline hid_t H5Create(const std::string& path)
{
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_file_locking(fapl, /*use_file_locking=*/0, /*ignore_when_disabled=*/1);
    hid_t fid  = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    H5Pclose(fapl);
    return fid;
}

}  // namespace h5utils
