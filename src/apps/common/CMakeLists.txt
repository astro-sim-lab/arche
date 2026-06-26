# collapse_common — shared application-layer headers for the one-zone collapse apps
#
# Header-only interface library. Carries the include path for the collapse_*
# helpers (collapse_defaults.h / collapse_dynamics.h / collapse_driver.h) that
# both collapse_primordial and collapse_metal_grain depend on. These are
# application-layer code; the chemistry kernel (src/include) does not reference
# them (dependency direction is app -> lib only). collapse_defaults.h pulls in
# solve/chemistry.h, so the arche library interface is propagated transitively.
add_library(collapse_common INTERFACE)
target_include_directories(collapse_common INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(collapse_common INTERFACE arche::arche)
