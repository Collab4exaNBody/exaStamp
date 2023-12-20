
add_executable(snap-compute-cg-cpp snap-compute-cg.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyCG.cpp)
target_include_directories(snap-compute-cg-cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../include)
target_compile_options(snap-compute-cg-cpp PRIVATE -DDEBUG_SNAP_CG)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/snap-compute-cg.sh" "${CMAKE_CURRENT_BINARY_DIR}/snap-compute-cg.sh")
xstamp_add_test(unit-snap-compute-cg ${CMAKE_CURRENT_BINARY_DIR}/snap-compute-cg.sh)

