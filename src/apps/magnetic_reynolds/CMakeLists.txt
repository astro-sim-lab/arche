add_executable(magnetic_reynolds main.cc)
target_link_libraries(magnetic_reynolds PRIVATE arche::arche)

target_compile_definitions(magnetic_reynolds PRIVATE
    MAGRE_CHEM_TABLE="data/metal_grain.h5"
)

# ── HDF5 input/output ──────────────────────────────────────────────────────────
set(HDF5_PREFER_PARALLEL FALSE)
find_package(HDF5 REQUIRED COMPONENTS C)
target_include_directories(magnetic_reynolds PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(magnetic_reynolds PRIVATE ${HDF5_C_LIBRARIES})
