# ==========================================
# === Rigid molecule droplet generator   ===
# ==========================================

add_executable(genDropletMolRigid genDropletMolRigid.cpp StampV3LegacyIOStructures.cpp)
target_include_directories(genDropletMolRigid PUBLIC ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(genDropletMolRigid ${MPI_CXX_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
target_compile_options(genDropletMolRigid PUBLIC ${OpenMP_CXX_FLAGS})

install(TARGETS genDropletMolRigid DESTINATION bin)



