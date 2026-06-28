add_executable(prim_collapse main.cc)
target_link_libraries(prim_collapse PRIVATE arche::arche collapse_common)

target_compile_definitions(prim_collapse PRIVATE
    PRIM_CHEM_TABLE="data/primordial.h5"
)

# ── Compact "minimal" variant ──────────────────────────────────────────────────
# Second target built from the same main.cc with -DARCHE_PRIM_MINIMAL, selecting
# the Nakauchi2019_Minimal network (compact 15-species / 33-reaction primordial
# model — a first-class composed network, not a mask over the full one).  It
# reuses the same primordial.h5 table for the partition functions.  Same
# pattern as the ARCHE_XRAY / ARCHE_SUBCYCLE compile-time variants.
add_executable(arche_collapse_prim_minimal main.cc)
target_link_libraries(arche_collapse_prim_minimal PRIVATE arche::arche collapse_common)

target_compile_definitions(arche_collapse_prim_minimal PRIVATE
    ARCHE_PRIM_MINIMAL
    PRIM_CHEM_TABLE="data/primordial.h5"
)

# ── HDF5 output ────────────────────────────────────────────────────────────────
set(HDF5_PREFER_PARALLEL FALSE)
find_package(HDF5 REQUIRED COMPONENTS C)

target_include_directories(prim_collapse PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(prim_collapse PRIVATE ${HDF5_C_LIBRARIES})

target_include_directories(arche_collapse_prim_minimal PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(arche_collapse_prim_minimal PRIVATE ${HDF5_C_LIBRARIES})
