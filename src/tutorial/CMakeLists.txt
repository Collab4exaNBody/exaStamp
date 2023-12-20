
set(exaStampTutorial_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors exanbMPI)
xstamp_add_plugin(exaStampTutorial ${CMAKE_CURRENT_SOURCE_DIR})

# copy tutorial exemples right into build dir to ease exercices
#configure_file(${PROJECT_SOURCE_DIR}/data/tutorial/sample1.msp ${PROJECT_BINARY_DIR}/tutorial/sample1.msp COPYONLY)
#configure_file(${PROJECT_SOURCE_DIR}/data/tutorial/sample2.msp ${PROJECT_BINARY_DIR}/tutorial/sample2.msp COPYONLY)
#configure_file(${PROJECT_SOURCE_DIR}/data/tutorial/sample3.msp ${PROJECT_BINARY_DIR}/tutorial/sample3.msp COPYONLY)

