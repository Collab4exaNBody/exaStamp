# pyexastamp: Python package for exaStamp, on top of pyexanbody and pyonika.
# _exastamp_data registers extractors for exaStamp types (ThermodynamicState,
# ParticleSpecies, MoleculeSpeciesVector) into pyonika's slot_as_array registry.
# They only use inline code from the headers, so no plugin library is linked,
# but the compile definitions that size these structures must match the plugins.
xnb_add_python_package(
  NAME pyexastamp
  PACKAGE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pyexastamp
  EXTENSION _exastamp_data
  SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/bind_exastamp_data.cpp
  INCLUDE_DIRECTORIES
    ${XNB_ROOT_DIR}/include
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src/thermo_state/include
    ${PROJECT_SOURCE_DIR}/src/particle_species/include
    ${PROJECT_SOURCE_DIR}/src/molecule/include
  COMPILE_DEFINITIONS
    XSTAMP_MAX_RIGID_MOLECULE_ATOMS=${XSTAMP_MAX_RIGID_MOLECULE_ATOMS}
    XSTAMP_MAX_MOLECULE_ATOMS=${XSTAMP_MAX_MOLECULE_ATOMS}
    XSTAMP_MAX_MOLECULE_BONDS=${XSTAMP_MAX_MOLECULE_BONDS}
    XSTAMP_MAX_MOLECULE_BENDS=${XSTAMP_MAX_MOLECULE_BENDS}
    XSTAMP_MAX_MOLECULE_TORSIONS=${XSTAMP_MAX_MOLECULE_TORSIONS}
    XSTAMP_MAX_MOLECULE_IMPROPERS=${XSTAMP_MAX_MOLECULE_IMPROPERS}
    XSTAMP_MAX_MOLECULE_PAIRS=${XSTAMP_MAX_MOLECULE_PAIRS}
  LINK_LIBRARIES onika
  VERSION ${EXASTAMP_VERSION}
  DEPENDENCIES pyonika pyexanbody
  DESCRIPTION "Python interface to the exaStamp molecular dynamics code")

file(GLOB PYEXASTAMP_EXAMPLES ${CMAKE_CURRENT_SOURCE_DIR}/examples/*.py)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/README.md DESTINATION python)
install(FILES ${PYEXASTAMP_EXAMPLES} DESTINATION python/examples)

if(ONIKA_REGRESSION_TEST_ENABLE)
  # multi-process checks, run under mpirun
  if(NOT MPIEXEC_NUMPROC_FLAG)
    set(MPIEXEC_NUMPROC_FLAG -n)
  endif()
  set(PYEXASTAMP_MPIRUN ${ONIKA_MPIRUN_CMD} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${Python3_EXECUTABLE})
  add_test(NAME pyexastamp_mpi COMMAND ${PYEXASTAMP_MPIRUN} ${CMAKE_CURRENT_SOURCE_DIR}/tests/mpi/check_mpi.py)
  add_test(NAME pyexastamp_mpi_abort COMMAND ${PYEXASTAMP_MPIRUN} ${CMAKE_CURRENT_SOURCE_DIR}/tests/mpi/check_abort.py)
  set_tests_properties(pyexastamp_mpi_abort PROPERTIES WILL_FAIL TRUE)
  set_tests_properties(pyexastamp_mpi pyexastamp_mpi_abort PROPERTIES TIMEOUT 300 ENVIRONMENT
    "PYTHONPATH=${PROJECT_BINARY_DIR}/python:${XNB_ROOT_DIR}/python:${onika_DIR}/python:$ENV{PYTHONPATH};ONIKA_PLUGIN_PATH=${onika_DIR}/plugins:${XNB_ROOT_DIR}/plugins:${PROJECT_BINARY_DIR}/plugins;ONIKA_CONFIG_PATH=${PROJECT_SOURCE_DIR}/data/config;OMP_NUM_THREADS=1")

  add_test(NAME pyexastamp_tests
           COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  set_tests_properties(pyexastamp_tests PROPERTIES ENVIRONMENT
    "PYTHONPATH=${PROJECT_BINARY_DIR}/python:${XNB_ROOT_DIR}/python:${onika_DIR}/python:$ENV{PYTHONPATH};ONIKA_PLUGIN_PATH=${onika_DIR}/plugins:${XNB_ROOT_DIR}/plugins:${PROJECT_BINARY_DIR}/plugins;ONIKA_CONFIG_PATH=${PROJECT_SOURCE_DIR}/data/config")
endif()
