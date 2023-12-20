
add_executable(snap-compute-bs-cpp snap-compute-bs.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyCG.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyGSH.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyBS.cpp)
target_include_directories(snap-compute-bs-cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../include)
target_compile_options(snap-compute-bs-cpp PRIVATE -DDEBUG_SNAP_BS)

add_executable(snap-compute-bs-cpp-lammps snap-compute-bs-lammps.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyCG.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyGSH.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyBS.cpp)
target_include_directories(snap-compute-bs-cpp-lammps PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../include)
target_compile_options(snap-compute-bs-cpp-lammps PRIVATE -DDEBUG_SNAP_BS -DLAMMPS=1)

file(COPY atom_positions.ref dbs.ref bs.ref2 Mo_conf.xyz bs_lammps.ref DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/snap-compute-bs.sh" "${CMAKE_CURRENT_BINARY_DIR}/snap-compute-bs.sh")
xstamp_add_test(utests-snap-compute-bs ${CMAKE_CURRENT_BINARY_DIR}/snap-compute-bs.sh)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/snap-compute-bs-lammps.sh" "${CMAKE_CURRENT_BINARY_DIR}/snap-compute-bs-lammps.sh")
xstamp_add_test(utests-snap-compute-bs-lammps ${CMAKE_CURRENT_BINARY_DIR}/snap-compute-bs-lammps.sh)

