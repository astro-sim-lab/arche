# Standalone external C consumer of libarche.
#
# Build it AGAINST AN INSTALLED arche, not from inside the arche tree:
#   cmake -B build -S . -DCMAKE_PREFIX_PATH=<arche-install-prefix>
#   cmake --build build
#   ./build/c_minimal <path-to>/primordial.h5
#
# The link language is forced to CXX so the C++ runtime that libarche needs is
# pulled in, even though the only source here is C.
cmake_minimum_required(VERSION 3.18)
project(arche_c_minimal LANGUAGES C CXX)

find_package(arche REQUIRED)

add_executable(c_minimal main.c)
target_link_libraries(c_minimal PRIVATE arche::arche)
set_target_properties(c_minimal PROPERTIES LINKER_LANGUAGE CXX)
