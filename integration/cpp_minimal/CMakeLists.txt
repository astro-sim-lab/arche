# Standalone external C++ consumer of libarche.
#
# Build it AGAINST AN INSTALLED arche, not from inside the arche tree:
#   cmake -B build -S . -DCMAKE_PREFIX_PATH=<arche-install-prefix>
#   cmake --build build
#   ./build/cpp_minimal <path-to>/primordial.h5
cmake_minimum_required(VERSION 3.18)
project(arche_cpp_minimal LANGUAGES CXX)

find_package(arche REQUIRED)

add_executable(cpp_minimal main.cc)
target_link_libraries(cpp_minimal PRIVATE arche::arche)
