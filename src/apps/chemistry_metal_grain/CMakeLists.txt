add_executable(metal_chem main.cc)
target_link_libraries(metal_chem PRIVATE chemistry)

target_compile_definitions(metal_chem PRIVATE
    DATA_DIR="data/metal_grain"
    METAL_DATA_DIR_DEF="data/metal_grain"
)
