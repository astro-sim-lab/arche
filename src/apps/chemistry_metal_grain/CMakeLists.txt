add_executable(metal_chem main.cc)
target_link_libraries(metal_chem PRIVATE arche::arche)

target_compile_definitions(metal_chem PRIVATE
    METAL_CHEM_TABLE="data/metal_grain.h5"
)

# ── HDF5 input ─────────────────────────────────────────────────────────────────
set(HDF5_PREFER_PARALLEL FALSE)
find_package(HDF5 REQUIRED COMPONENTS C)
target_include_directories(metal_chem PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(metal_chem PRIVATE ${HDF5_C_LIBRARIES})
