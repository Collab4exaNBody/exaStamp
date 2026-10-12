# ---------------------------------------------------------------------------
# AddSingleTest(TEST_NAME CMDVAR NUMPROCS NUMCORES [INPUT_DIR])
#
# Registers a single test with CTest and creates a matching custom target
# for running it through the build system (e.g. make <TEST_NAME>).
#
# Arguments:
#   TEST_NAME  - Unique CTest name for this test
#   CMDVAR     - Name of the CMake variable holding the command list
#   NUMPROCS   - Number of MPI processes
#   NUMCORES   - Number of OpenMP threads / CPU cores
#   INPUT_DIR  - Optional folder holding the test input files. The test then
#                runs in <build>/regression_runs/<INPUT_DIR relative to data/regression_new>,
#                shared by the tests of that folder and kept between runs, holding a copy
#                of INPUT_DIR and a 'data' link to data/potentials (cmake/run_in_test_copy.sh)
# ---------------------------------------------------------------------------
function(AddSingleTest TEST_NAME CMDVAR NUMPROCS NUMCORES)
    MakeRunCommand(${CMDVAR} ${NUMPROCS} ${NUMCORES} FULL_COMMAND)
    MakeDebugRunCommand(${CMDVAR} ${NUMPROCS} ${NUMCORES} FULL_COMMAND_DBG)

    set(_wrapper "")
    if(${ARGC} GREATER 4)
        file(RELATIVE_PATH _run_subdir "${CMAKE_SOURCE_DIR}/data/regression_new" "${ARGV4}")
        if(_run_subdir STREQUAL "" OR _run_subdir MATCHES "^\\.\\.")
            set(_run_subdir ${TEST_NAME})
        endif()
        set(_wrapper bash "${CMAKE_SOURCE_DIR}/cmake/run_in_test_copy.sh"
            "${CMAKE_BINARY_DIR}/regression_runs/${_run_subdir}" "${ARGV4}" "${CMAKE_SOURCE_DIR}/data/potentials")
    endif()

    add_test(NAME ${TEST_NAME} COMMAND ${_wrapper} ${FULL_COMMAND})
    add_custom_target(${TEST_NAME} COMMAND ${_wrapper} ${FULL_COMMAND_DBG})
    set_tests_properties(${TEST_NAME} PROPERTIES
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    )
endfunction()

# function(AddSingleTest TEST_NAME CMDVAR NUMPROCS NUMCORES)
# #  message(STATUS "AddTestWithDebugTarget ${TEST_NAME} | ${CMDVAR} (${${CMDVAR}}) | ${NUMPROCS} | ${NUMCORES}")
#   MakeRunCommand(${CMDVAR} ${NUMPROCS} ${NUMCORES} FULL_COMMAND)
#   MakeDebugRunCommand(${CMDVAR} ${NUMPROCS} ${NUMCORES} FULL_COMMAND_DBG)
#   # message(STATUS "${TEST_NAME} | ${FULL_COMMAND} | ${FULL_COMMAND_DBG}")
#   add_test(NAME ${TEST_NAME} COMMAND ${FULL_COMMAND})
#   add_custom_target(${TEST_NAME} COMMAND ${FULL_COMMAND_DBG})
#   set_tests_properties(${TEST_NAME} PROPERTIES WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
# endfunction()
