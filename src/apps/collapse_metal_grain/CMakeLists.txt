add_executable(metal_collapse main.cc)
target_link_libraries(metal_collapse PRIVATE chemistry)

target_compile_definitions(metal_collapse PRIVATE
    DATA_DIR="data/metal_grain"
    METAL_DATA_DIR="data/metal_grain"
)

# ── HDF5 output ────────────────────────────────────────────────────────────────
set(HDF5_PREFER_PARALLEL FALSE)
find_package(HDF5 REQUIRED COMPONENTS C)

target_include_directories(metal_collapse PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(metal_collapse PRIVATE ${HDF5_C_LIBRARIES})
