# ExaStamp specific packages

add_subdirectory(debug)
add_subdirectory(analytics)
add_subdirectory(particle_species)
add_subdirectory(io)
add_subdirectory(compute)
add_subdirectory(mpi)
add_subdirectory(parrinellorahman)
add_subdirectory(ttm)
add_subdirectory(potential_factory)
add_subdirectory(potential)
add_subdirectory(particle_charge)
add_subdirectory(molecule)
add_subdirectory(mechanical)
add_subdirectory(rigidmol)
add_subdirectory(tutorial)
add_subdirectory(flat_arrays)

if(XSTAMP_ENABLE_SLIPLINK)
  add_subdirectory(sliplink)
endif()

