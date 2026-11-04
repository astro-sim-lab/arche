# The bundled applications instantiate the kernel templates in their own
# translation units, so they are compiled with the same CPU and floating-point
# flags as libarche (ARCHE_CXX_CPU_FLAGS, set in the top-level CMakeLists.txt).
add_compile_options(${ARCHE_CXX_CPU_FLAGS})

add_subdirectory(apps/common)
add_subdirectory(apps/collapse_primordial)
add_subdirectory(apps/collapse_metal_grain)
add_subdirectory(apps/chemistry_primordial)
add_subdirectory(apps/chemistry_metal_grain)
add_subdirectory(apps/magnetic_reynolds)
