// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
#include <array>
#include <vector>
#include <utility>
#include <string>
#include <fstream>
#include <stdexcept>
#include "species.h"

namespace chemistry {

// ---------------------------------------------------------------------------
// ReactionTable<N_sp, N_react>
//
// Read-only after initialization. Shared across threads/cells.
// Loaded once at startup from data files.
//
// Index convention:
//   Fortran uses 1-based species indices; 0 and 101 are sentinel values.
//   C++ stores them as-is (int). Access helpers map to 0-based arrays
//   via y_ext[] (see solver.h).
// ---------------------------------------------------------------------------
template<int N_sp_, int N_react_>
struct ReactionTable {
    static constexpr int N_sp    = N_sp_;
    static constexpr int N_react = N_react_;

    // Reaction indices (Fortran 1-based; 0 = vacant slot, 101 = photon)
    std::array<int, N_react> rean{};  // reaction number (= 1..N_react)
    std::array<int, N_react> rea1{};  // reactant 1 species index
    std::array<int, N_react> rea2{};  // reactant 2 species index
    std::array<int, N_react> rea3{};  // reactant 3 species index
    std::array<int, N_react> pro1{};  // product 1 species index
    std::array<int, N_react> pro2{};  // product 2 species index
    std::array<int, N_react> pro3{};  // product 3 species index
    std::array<int, N_react> nrea{};  // number of reactants
    std::array<int, N_react> npro{};  // number of products

    std::array<double, N_react> Cmass{};  // mass factor for equilibrium constant
    std::array<double, N_react> delE{};   // reaction enthalpy [erg]

    // Partition function tables: pf_table[isp] = vector of (T_K, pf_value) pairs
    // Loaded from pf_*.dat files (two-column: T[K]  pf).
    // isp is 1-based (0 unused) to match Fortran indexing.
    // Only species loaded from files have non-empty vectors; analytic species
    // are computed directly in eval_partition_functions().
    std::array<std::vector<std::pair<double,double>>, N_sp + 2> pf_table{};

    // Species masses [g], indexed 1-based (element 0 unused)
    std::array<double, N_sp + 1> mass{};

    // Number of reactions actually loaded (may be < N_react for CR reactions
    // which are hardcoded in react_coef rather than from the data file)
    int n_loaded = 0;

    // ── Grain surface reaction table (metal_grain only) ───────────────────────
    // Loaded from react_grain_surface.dat via load_grain_surface_table().
    // Indexed 0-based sequentially (entry 0 = first file line).
    // All species indices are 1-based (Fortran convention); 0 = vacant sentinel.
    static constexpr int N_grain = 200;
    std::array<int, N_grain> grean{};  // sequential index (= ire + 1)
    std::array<int, N_grain> grea1{};  // reactant 1 species (0=sentinel)
    std::array<int, N_grain> grea2{};  // reactant 2 species (1-based)
    std::array<int, N_grain> gpro1{};  // product  1 species (0=sentinel)
    std::array<int, N_grain> gpro2{};  // product  2 species (1-based)
    std::array<int, N_grain> gnrea{};  // exponent: rate ∝ xnH^gnrea
    int n_grain = 0;                   // number of grain surface reactions loaded

    // ── Saha equilibrium data (high-density branch, xnH > 1e18) ─────────────
    // Read from react_prm_saha.dat. Up to 100 entries.
    // nsp == 1 → r1+r2 → p1+photon  (type 1, entries 0-15)
    // nsp == 2 → r1+r2 → p1+p2      (type 2, entries 16-17)
    static constexpr int N_saha = 100;
    std::array<int, N_saha>    saha_num{};   // reaction number (index into Keqb)
    std::array<int, N_saha>    saha_rs1{};   // reactant 1 species (1-based)
    std::array<int, N_saha>    saha_rs2{};   // reactant 2 species (1-based)
    std::array<int, N_saha>    saha_ps1{};   // product 1 species  (1-based)
    std::array<int, N_saha>    saha_ps2{};   // product 2 species  (1-based; 101=photon)
    std::array<int, N_saha>    saha_nsr{};   // number of reactants
    std::array<int, N_saha>    saha_nsp{};   // number of products
    std::array<double, N_saha> saha_Cm{};    // mass factor
    std::array<double, N_saha> saha_dE{};    // enthalpy [erg]
    int n_saha = 0;
};

// Convenience aliases
using ZeroMetalTable  = ReactionTable<zero_metal::N_sp,  zero_metal::N_react>;
using MetalGrainTable = ReactionTable<metal_grain::N_sp, metal_grain::N_react>;

// ---------------------------------------------------------------------------
// load_reaction_table()
//
// Reads react_prm.dat (or react_metal_grain.dat) into a ReactionTable.
// File format (space-delimited):
//   rean  rea1  rea2  rea3  pro1  pro2  pro3  nrea  npro  Cmass  delE
// Fortran equivalent: read_reaction subroutine.
// ---------------------------------------------------------------------------
template<int N_sp, int N_react>
void load_reaction_table(ReactionTable<N_sp, N_react>& tbl,
                         const std::string& react_dat_path)
{
    std::ifstream f(react_dat_path);
    if (!f) throw std::runtime_error("Cannot open: " + react_dat_path);

    int i = 0;
    while (i < N_react) {
        int rn, r1, r2, r3, p1, p2, p3, nr, np;
        double cm, de;
        if (!(f >> rn >> r1 >> r2 >> r3 >> p1 >> p2 >> p3 >> nr >> np >> cm >> de))
            break;
        // Store at sequential file-position index i (0-based).
        // tbl.rean[i] = rn  (the 1-based reaction number of the i-th file line).
        // This matches the Fortran sequential array access pattern used in
        // react_rat / compute_rates: rean(ire) → tbl.rean[ire-1].
        int idx = i;
        tbl.rean[idx] = rn;
        tbl.rea1[idx] = r1;
        tbl.rea2[idx] = r2;
        tbl.rea3[idx] = r3;
        tbl.pro1[idx] = p1;
        tbl.pro2[idx] = p2;
        tbl.pro3[idx] = p3;
        tbl.nrea[idx] = nr;
        tbl.npro[idx] = np;
        tbl.Cmass[idx] = cm;
        tbl.delE[idx]  = de;
        ++i;
    }
    tbl.n_loaded = i;
}

// ---------------------------------------------------------------------------
// load_partition_functions()
//
// Reads pf_*.dat files for species that have tabulated partition functions.
// Each file has 102 lines: line 0 = temperature [K], lines 1-101 = log10(pf).
// Fortran equivalent: read_pf (part_fnc_prm.f / part_fnc_metal.f).
//
// pf_files: list of (species_index_1based, filename) pairs.
// For species without a pf file, pf_table[isp] remains empty
// (the partition function evaluator will return 1.0 as default).
// ---------------------------------------------------------------------------
template<int N_sp, int N_react>
void load_partition_functions(ReactionTable<N_sp, N_react>& tbl,
    const std::vector<std::pair<int, std::string>>& pf_files)
{
    for (auto& [isp, path] : pf_files) {
        if (isp < 1 || isp > N_sp)
            throw std::runtime_error("PF species index out of range: " + std::to_string(isp));
        std::ifstream f(path);
        if (!f) throw std::runtime_error("Cannot open PF file: " + path);
        tbl.pf_table[isp].clear();
        double T_val, pf_val;
        while (f >> T_val >> pf_val)
            tbl.pf_table[isp].emplace_back(T_val, pf_val);
        if (tbl.pf_table[isp].empty())
            throw std::runtime_error("Empty PF file: " + path);
    }
}

// ---------------------------------------------------------------------------
// load_mass_table()
//
// Reads mass.dat into tbl.mass[1..N_sp] (1-based).
// File format: one mass [g] per line, ordered by species index.
// ---------------------------------------------------------------------------
template<int N_sp, int N_react>
void load_mass_table(ReactionTable<N_sp, N_react>& tbl,
                     const std::string& mass_dat_path)
{
    std::ifstream f(mass_dat_path);
    if (!f) throw std::runtime_error("Cannot open: " + mass_dat_path);
    for (int isp = 1; isp <= N_sp; ++isp) {
        if (!(f >> tbl.mass[isp]))
            throw std::runtime_error("Short read in mass.dat at species " + std::to_string(isp));
    }
}

// ---------------------------------------------------------------------------
// load_saha_table()
//
// Reads react_prm_saha.dat into tbl.saha_* arrays.
// File format (space-delimited, list-directed Fortran):
//   num rs1 rs2 ps1 ps2 nsr nsp Cm dE
// Entries with nsp==1 are type-1 (r1+r2→p1+photon).
// Entries with nsp==2 are type-2 (r1+r2→p1+p2).
// ---------------------------------------------------------------------------
template<int N_sp, int N_react>
void load_saha_table(ReactionTable<N_sp, N_react>& tbl,
                     const std::string& saha_dat_path)
{
    std::ifstream f(saha_dat_path);
    if (!f) throw std::runtime_error("Cannot open saha file: " + saha_dat_path);

    int i = 0;
    while (i < ReactionTable<N_sp, N_react>::N_saha) {
        int num, r1, r2, p1, p2, nr, np;
        double cm, de;
        if (!(f >> num >> r1 >> r2 >> p1 >> p2 >> nr >> np >> cm >> de))
            break;
        tbl.saha_num[i] = num;
        tbl.saha_rs1[i] = r1;
        tbl.saha_rs2[i] = r2;
        tbl.saha_ps1[i] = p1;
        tbl.saha_ps2[i] = p2;
        tbl.saha_nsr[i] = nr;
        tbl.saha_nsp[i] = np;
        tbl.saha_Cm[i]  = cm;
        tbl.saha_dE[i]  = de;
        ++i;
    }
    tbl.n_saha = i;
}

// ---------------------------------------------------------------------------
// load_grain_surface_table()
//
// Reads react_grain_surface.dat into tbl.grean/grea1/grea2/gpro1/gpro2/gnrea.
// File format (space-delimited, up to 200 entries):
//   grean  grea1  grea2  gpro1  gpro2  gnrea
// Fortran equivalent: read_grain_reaction subroutine.
// ---------------------------------------------------------------------------
template<int N_sp, int N_react>
void load_grain_surface_table(ReactionTable<N_sp, N_react>& tbl,
                               const std::string& grain_dat_path)
{
    std::ifstream f(grain_dat_path);
    if (!f) throw std::runtime_error("Cannot open: " + grain_dat_path);

    int i = 0;
    static constexpr int N_grain = ReactionTable<N_sp, N_react>::N_grain;
    while (i < N_grain) {
        int ge, g1, g2, gp1, gp2, gnr;
        if (!(f >> ge >> g1 >> g2 >> gp1 >> gp2 >> gnr))
            break;
        tbl.grean[i] = ge;
        tbl.grea1[i] = g1;
        tbl.grea2[i] = g2;
        tbl.gpro1[i] = gp1;
        tbl.gpro2[i] = gp2;
        tbl.gnrea[i] = gnr;
        ++i;
    }
    tbl.n_grain = i;
}

} // namespace chemistry
