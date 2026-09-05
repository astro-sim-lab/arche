cmake_minimum_required(VERSION 3.18)
project(arche VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

include(GNUInstallDirs)

# --- Eigen ---
find_package(Eigen3 REQUIRED NO_MODULE)

# --- HDF5 (used by chemistry table loader; must be visible to all apps) ---
set(HDF5_PREFER_PARALLEL FALSE)
find_package(HDF5 REQUIRED COMPONENTS C)

# --- arche chemistry kernel (STATIC library: libarche) ---
# The kernel is templated and header-only, but the non-template public facade
# (src/include/api/arche_api.h) is compiled and explicitly instantiated for both
# models in src/lib/arche_api.cc.  Building it into a STATIC library means the
# kernel object code lives in libarche once: consumer translation units link it
# instead of re-instantiating the templates, and an external app can drive the
# network via api/arche_api.h alone (no Eigen/HDF5 in its own sources).
#
# Thermo inversion implementations intentionally remain separate translation
# units for numerical reproducibility. Moving them into arche_api.cc changes
# the large unit's optimization budget and can alter ULP-level code generation
# for its inline-template COMDAT definitions (including metal line cooling),
# affecting existing linked applications.  LTO/IPO removes this isolation and
# restores that code-generation coupling.
#
# The target is named `arche` (artefact libarche.a) and is also exposed as the
# namespaced alias `arche::arche`, which is the name external consumers should
# link (via add_subdirectory / FetchContent). The `::` makes CMake treat it as a
# target and error early on a typo. All usage requirements are PUBLIC so they
# apply both to libarche's own compilation and to every consumer.
add_library(arche STATIC
    src/lib/arche_api.cc
    src/lib/arche_capi.cc
    src/thermo/eos.cc
)
add_library(arche::arche ALIAS arche)
# Compile the static archive as position-independent code.  Modern toolchains
# link executables as PIE by default, and consumers may embed libarche into a
# shared library; without -fPIC the archive's absolute (R_X86_64_32)
# relocations cannot go into a PIE/shared consumer ("recompile with -fPIE").
set_target_properties(arche PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(arche PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/arche>
    ${HDF5_INCLUDE_DIRS}
)
target_link_libraries(arche PUBLIC Eigen3::Eigen ${HDF5_C_LIBRARIES})
option(ARCHE_XRAY     "Enable X-ray background radiation (primary + secondary ionization + heating)" OFF)  # experimental; disabled by default until validated
option(ARCHE_SUBCYCLE "Enable NR subcycling: halve dt and retry when chemcool NR does not converge" OFF)

target_compile_definitions(arche PUBLIC
    ARCHE_USE_EIGEN
    $<$<BOOL:${ARCHE_XRAY}>:ARCHE_XRAY>
    $<$<BOOL:${ARCHE_SUBCYCLE}>:ARCHE_SUBCYCLE>
)
target_compile_options(arche PUBLIC
    $<$<CONFIG:Release>:
        $<$<CXX_COMPILER_ID:GNU>: -march=native>
        $<$<CXX_COMPILER_ID:IntelLLVM>: -xHost -fp-model=precise>
    >
)

# --- install / export (find_package(arche)) ---
# Only the lightweight public-facade headers are installed (the closure reachable
# from api/arche_api.h / api/arche_capi.h): api/arche_api.h -> core/state.h ->
# core/species_index.h -> core/species_catalog.h, plus the standalone C ABI
# api/arche_capi.h. The Eigen/HDF5-
# dependent kernel headers (solve/chemistry.h, solve/chemreact.h,
# kinetics/topology.h, ...) and the internal rate-kernel bookkeeping
# (kinetics/reaction_index.h) are deliberately NOT shipped, so an installed
# consumer keeps the same "no Eigen/HDF5 in your own sources" contract.
option(ARCHE_INSTALL "Generate install/export rules for libarche (find_package(arche))" ON)
if(ARCHE_INSTALL)
    include(CMakePackageConfigHelpers)

    install(TARGETS arche
        EXPORT archeTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
    # The subdirectory layout under src/include is preserved on install so that
    # the facade headers' root-relative includes (api/arche_api.h -> core/state.h
    # -> core/species_index.h) resolve against the installed arche/ include root.
    install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/include/api/arche_api.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/include/api/arche_capi.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arche/api
    )
    install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/include/core/state.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/include/core/species_index.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/include/core/species_catalog.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arche/core
    )
    install(EXPORT archeTargets
        FILE archeTargets.cmake
        NAMESPACE arche::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/arche
    )

    configure_package_config_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/archeConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/archeConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/arche
    )
    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/archeConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/archeConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/archeConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/arche
    )
endif()

# --- applications ---
add_subdirectory(src)
