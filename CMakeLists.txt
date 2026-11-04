cmake_minimum_required(VERSION 3.18)
project(arche VERSION 2.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

# Default to an optimized build when ARCHE is the top-level project and no
# build type is given (single-config generators only).  An embedding project
# keeps its own CMAKE_BUILD_TYPE.
get_property(ARCHE_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR AND NOT ARCHE_MULTI_CONFIG AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Build type (Debug, Release, RelWithDebInfo, MinSizeRel)" FORCE)
    message(STATUS "arche: CMAKE_BUILD_TYPE not set; using Release")
endif()

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
    src/cooling/cooling_memo.cc
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

target_compile_definitions(arche PUBLIC
    ARCHE_USE_EIGEN
    $<$<BOOL:${ARCHE_XRAY}>:ARCHE_XRAY>
)

# --- target CPU and floating-point model ---
# ARCHE_MARCH is passed as -march=<value> to every C++ source of libarche and
# of the bundled applications, in every build type.  "native" (default) targets
# the CPU of the build machine.  When the machines that run ARCHE differ from
# the build machine (e.g. a login node and compute nodes of different CPUs),
# name the run machines' instruction set instead, e.g. x86-64-v3 (AVX2) or
# x86-64-v4 (AVX-512).  An empty value adds no -march flag (compiler default).
#
# Intel oneAPI icpx additionally gets -fp-model=precise, because its default
# floating-point model allows value-changing optimizations.  -march is used for
# icpx as well, not -x<code>/-xHost: the -x options embed an Intel-only
# processor check, and the resulting binaries abort on AMD processors.
#
# These flags are not exported through arche::arche: a project that links
# libarche keeps its own CPU and floating-point settings.  The bundled
# applications receive them as directory options in src/CMakeLists.txt, so
# that the kernel templates they instantiate are compiled exactly as libarche.
set(ARCHE_MARCH "native" CACHE STRING
    "Target CPU passed as -march=<value> (empty: no -march flag)")
set(ARCHE_CPU_FLAGS "")
if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|IntelLLVM)$")
    if(NOT ARCHE_MARCH STREQUAL "")
        list(APPEND ARCHE_CPU_FLAGS "-march=${ARCHE_MARCH}")
    endif()
    if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        list(APPEND ARCHE_CPU_FLAGS "-fp-model=precise")
    endif()
else()
    message(WARNING "arche: no CPU or floating-point flags are set for the C++ "
                    "compiler '${CMAKE_CXX_COMPILER_ID}' (supported: GNU, Clang, "
                    "IntelLLVM). The classic Intel compiler (icpc) is not supported; use icpx.")
endif()
# C++ sources only: C sources are compiled by the C compiler, which need not
# accept these flags.
set(ARCHE_CXX_CPU_FLAGS "")
foreach(flag IN LISTS ARCHE_CPU_FLAGS)
    list(APPEND ARCHE_CXX_CPU_FLAGS "$<$<COMPILE_LANGUAGE:CXX>:${flag}>")
endforeach()
message(STATUS "arche: C++ compiler ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}, "
               "CPU/floating-point flags: ${ARCHE_CPU_FLAGS}")
target_compile_options(arche PRIVATE ${ARCHE_CXX_CPU_FLAGS})
# Warnings are PRIVATE: they apply to libarche's own translation units and are
# not forced on projects that link arche::arche.
target_compile_options(arche PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang,IntelLLVM>: -Wall -Wextra>
)

# --- install / export (find_package(arche)) ---
# Only the lightweight public-facade headers are installed (the closure
# reachable from api/arche_api.h / api/arche_capi.h): api/arche_api.h ->
# core/state.h -> core/species_index.h -> core/species_catalog.h, plus the
# standalone C ABI api/arche_capi.h.  The Eigen/HDF5-dependent kernel headers
# (solve/chemistry.h, solve/chemreact.h, kinetics/topology.h, ...), internal
# thermo headers (thermo/eos.h, thermo/h2_thermodynamics.h), and internal
# rate-kernel bookkeeping (kinetics/reaction_index.h) are deliberately NOT
# shipped.  Installed consumers use only the facade thermo queries and retain
# the "no Eigen/HDF5 in your own sources" contract.
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
