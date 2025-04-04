
add_executable(snap-compute-gsh-cpp snap-compute-gsh.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyGSH.cpp)
target_include_directories(snap-compute-gsh-cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../include)
target_compile_options(snap-compute-gsh-cpp PRIVATE -DDEBUG_SNAP_GSH)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/snap-compute-gsh.sh" "${CMAKE_CURRENT_BINARY_DIR}/snap-compute-gsh.sh")
onika_add_test(unit-snap-compute-gsh ${CMAKE_CURRENT_BINARY_DIR}/snap-compute-gsh.sh)

