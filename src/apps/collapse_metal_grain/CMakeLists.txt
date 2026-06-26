add_executable(metal_collapse main.cc)
target_link_libraries(metal_collapse PRIVATE arche::arche collapse_common)

target_compile_definitions(metal_collapse PRIVATE
    METAL_CHEM_TABLE="data/metal_grain.h5"
)

# ── Compact "minimal" variant ──────────────────────────────────────────────────
# Second target built from the same main.cc with -DARCHE_METAL_MINIMAL, selecting
# the Nakauchi2021_Minimal network (compact 40-species / 103+10-reaction metal
# model — a first-class composed network, not a mask over the full one).  Same
# pattern as collapse_primordial's ARCHE_PRIM_MINIMAL; it reuses metal_grain.h5
# for the partition functions (the compact table retains the full PF-loaded table
# for the strategy-1 rate path and the compact Saha).  main.cc selects the
# variant through the AppCell/AppSp abstraction at its head; the compact build
# omits the K/Na/Gr2+/He++ species the network drops.
add_executable(arche_collapse_metal_minimal main.cc)
target_link_libraries(arche_collapse_metal_minimal PRIVATE arche::arche collapse_common)

target_compile_definitions(arche_collapse_metal_minimal PRIVATE
    ARCHE_METAL_MINIMAL
    METAL_CHEM_TABLE="data/metal_grain.h5"
)

# ── HDF5 output ────────────────────────────────────────────────────────────────
set(HDF5_PREFER_PARALLEL FALSE)
find_package(HDF5 REQUIRED COMPONENTS C)

target_include_directories(metal_collapse PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(metal_collapse PRIVATE ${HDF5_C_LIBRARIES})

target_include_directories(arche_collapse_metal_minimal PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(arche_collapse_metal_minimal PRIVATE ${HDF5_C_LIBRARIES})
