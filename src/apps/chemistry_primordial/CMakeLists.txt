add_executable(prim_chem main.cc)
target_link_libraries(prim_chem PRIVATE arche::arche)

target_compile_definitions(prim_chem PRIVATE
    PRIM_CHEM_TABLE="data/primordial.h5"
)

# ── HDF5 input ─────────────────────────────────────────────────────────────────
set(HDF5_PREFER_PARALLEL FALSE)
find_package(HDF5 REQUIRED COMPONENTS C)
target_include_directories(prim_chem PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(prim_chem PRIVATE ${HDF5_C_LIBRARIES})
