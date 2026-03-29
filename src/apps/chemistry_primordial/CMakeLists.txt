add_executable(prim_chem main.cc)
target_link_libraries(prim_chem PRIVATE chemistry)

target_compile_definitions(prim_chem PRIVATE
    DATA_DIR="data/primordial"
)
