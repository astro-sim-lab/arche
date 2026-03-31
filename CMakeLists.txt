cmake_minimum_required(VERSION 3.18)
project(chemistry_cpp VERSION 0.1.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

# --- Eigen ---
find_package(Eigen3 REQUIRED NO_MODULE)

# --- chemistry kernel (header-only interface library) ---
add_library(chemistry INTERFACE)
target_include_directories(chemistry INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/include>
)
target_link_libraries(chemistry INTERFACE Eigen3::Eigen)
target_compile_definitions(chemistry INTERFACE CHEMISTRY_USE_EIGEN)
target_compile_options(chemistry INTERFACE
    $<$<CONFIG:Release>:
        $<$<CXX_COMPILER_ID:GNU>: -march=native>
        $<$<CXX_COMPILER_ID:IntelLLVM>: -xHost>
    >
)

# --- applications ---
add_subdirectory(src)
