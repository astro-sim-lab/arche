add_executable(prim_collapse main.cc)
target_link_libraries(prim_collapse PRIVATE chemistry)

target_compile_definitions(prim_collapse PRIVATE
    DATA_DIR="data/primordial"
)

# ── HDF5 output ────────────────────────────────────────────────────────────────
set(HDF5_PREFER_PARALLEL FALSE)
find_package(HDF5 REQUIRED COMPONENTS C)

target_include_directories(prim_collapse PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(prim_collapse PRIVATE ${HDF5_C_LIBRARIES})
