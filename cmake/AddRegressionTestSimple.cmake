# ---------------------------------------------------------------------------
# _ParseMspHeader(MSP_FILE)
#
# Internal helper — reads the first line of an .msp file and populates the
# following variables in the PARENT scope:
#
#   ENABLE_TEST_SEQ  / ENABLE_TEST_MT  / ENABLE_TEST_MPI  (ON/OFF)
#   TEST_ENABLED     — set to OFF if any enable-if-<VAR>[=<VAL>] condition fails
#   MT_NTHREADS      — thread count override (default: 4)
# ---------------------------------------------------------------------------
function(_ParseMspHeader MSP_FILE)
    # ---- Read first line ------------------------------------------------
    file(STRINGS "${MSP_FILE}" _all_lines)
    list(GET _all_lines 0 _header)

    # ---- no-{mpi,seq,mt} flags ------------------------------------------
    foreach(_mode MPI SEQ MT)
        string(TOLOWER "${_mode}" _mode_lc)
        string(FIND "${_header}" "no-${_mode_lc}" _pos)
        if(_pos EQUAL -1)
            set(ENABLE_TEST_${_mode} ON  PARENT_SCOPE)
        else()
            set(ENABLE_TEST_${_mode} OFF PARENT_SCOPE)
        endif()
    endforeach()

    # ---- enable-if-<VAR>[=<VALUE>] conditions ---------------------------
    set(_test_enabled ON)
    string(REGEX MATCHALL "enable-if-[^ ]+" _conditions "${_header}")

    foreach(_cond ${_conditions})
        string(REPLACE "enable-if-" "" _cond "${_cond}")

        # Check for optional =<value> suffix
        string(REGEX MATCH "=.*" _expected_value "${_cond}")
        if(_expected_value)
            string(REGEX REPLACE "=.*" "" _var_name   "${_cond}")
            string(REPLACE "="    ""   _expected_value "${_expected_value}")
            if(NOT "${${_var_name}}" STREQUAL "${_expected_value}")
                set(_test_enabled OFF)
            endif()
        else()
            # No value specified — variable must be truthy
            if(NOT ${_cond})
                set(_test_enabled OFF)
            endif()
        endif()
    endforeach()

    if(NOT _test_enabled)
        set(ENABLE_TEST_MPI OFF PARENT_SCOPE)
        set(ENABLE_TEST_SEQ OFF PARENT_SCOPE)
        set(ENABLE_TEST_MT  OFF PARENT_SCOPE)
    endif()
    set(TEST_ENABLED ${_test_enabled} PARENT_SCOPE)

    # ---- nthreads override ----------------------------------------------
    set(_nthreads 4)
    string(REGEX MATCH "nthreads=[^ ]+" _nthreads_kv "${_header}")
    if(_nthreads_kv)
        string(REPLACE "nthreads=" "" _nthreads_kv "${_nthreads_kv}")
        if("${_nthreads_kv}" STREQUAL "max")
            set(_nthreads ${ONIKA_HOST_HW_CORES})
        elseif(${_nthreads_kv} GREATER 1)
            set(_nthreads ${_nthreads_kv})
        endif()
    endif()
    set(MT_NTHREADS ${_nthreads} PARENT_SCOPE)
endfunction()


# ---------------------------------------------------------------------------
# _RegisterMspTest(PREFIX TESTCASE TESTVARIANT MSP_FILE)
#
# Internal helper — registers the seq / mt / par variants of a single .msp
# test, using the flags populated by _ParseMspHeader().
#
# The resulting CTest names follow the pattern:
#   <APP>_<PREFIX>_<TESTCASE>_<TESTVARIANT>_{seq|mt|par}
#
# where PREFIX is typically the top-level folder name (e.g. "flexible_molecule")
# so that `ctest -R flexible_molecule` selects all tests from that folder.
# ---------------------------------------------------------------------------
function(_RegisterMspTest PREFIX TESTCASE TESTVARIANT MSP_FILE)
    _ParseMspHeader("${MSP_FILE}")

    if(NOT ONIKA_RUN_WRAPPER)
        set(ONIKA_RUN_WRAPPER ${ONIKA_RUN})
    endif()

    set(_registered OFF)
    set(_base_name ${XNB_APP_NAME}_${PREFIX}_${TESTCASE}_${TESTVARIANT})

    if(ONIKA_REGRESSION_TEST_ENABLE_SEQ AND ENABLE_TEST_SEQ)
        set(_registered ON)
        set(_test_name ${_base_name}_seq)
        set(${_test_name}Command ${ONIKA_RUN_WRAPPER} "${MSP_FILE}" ${EXASTAMP_TEST_ADDITIONAL_ARGS})
        AddSingleTest(${_test_name} ${_test_name}Command 1 1)
    endif()

    if(ONIKA_REGRESSION_TEST_ENABLE_MT AND ENABLE_TEST_MT)
        set(_registered ON)
        set(_test_name ${_base_name}_mt)
        set(${_test_name}Command ${ONIKA_RUN_WRAPPER} "${MSP_FILE}" ${EXASTAMP_TEST_ADDITIONAL_ARGS})
        AddSingleTest(${_test_name} ${_test_name}Command 1 ${MT_NTHREADS})
    endif()

    if(ONIKA_REGRESSION_TEST_ENABLE_MPI AND ENABLE_TEST_MPI)
        set(_registered ON)
        set(_test_name ${_base_name}_par)
        set(${_test_name}Command ${ONIKA_RUN_WRAPPER} "${MSP_FILE}" ${EXASTAMP_TEST_ADDITIONAL_ARGS})
        AddSingleTest(${_test_name} ${_test_name}Command 4 ${MT_NTHREADS})
    endif()

    if(TEST_ENABLED AND NOT _registered)
        message(STATUS "Warning: ${_base_name} is enabled but not registered \
(check ONIKA_REGRESSION_TEST_ENABLE_* flags)")
    endif()
endfunction()

# ---------------------------------------------------------------------------
# _RegisterMspFlatTest(TESTCASE TESTVARIANT MSP_FILE)
#
# Internal helper — registers the seq / mt / par variants of a single .msp
# test, using the flags populated by _ParseMspHeader().
#
# The resulting CTest names follow the pattern:
#   <APP>_<TESTCASE>_<TESTVARIANT>_{seq|mt|par}
#
# ---------------------------------------------------------------------------
function(_RegisterMspFlatTest TESTCASE TESTVARIANT MSP_FILE)
    _ParseMspHeader("${MSP_FILE}")

    if(NOT ONIKA_RUN_WRAPPER)
        set(ONIKA_RUN_WRAPPER ${ONIKA_RUN})
    endif()

    set(_registered OFF)
    set(_base_name ${XNB_APP_NAME}_${TESTCASE}_${TESTVARIANT})

    if(ONIKA_REGRESSION_TEST_ENABLE_SEQ AND ENABLE_TEST_SEQ)
        set(_registered ON)
        set(_test_name ${_base_name}_seq)
        set(${_test_name}Command ${ONIKA_RUN_WRAPPER} "${MSP_FILE}" ${EXASTAMP_TEST_ADDITIONAL_ARGS})
        AddSingleTest(${_test_name} ${_test_name}Command 1 1)
    endif()

    if(ONIKA_REGRESSION_TEST_ENABLE_MT AND ENABLE_TEST_MT)
        set(_registered ON)
        set(_test_name ${_base_name}_mt)
        set(${_test_name}Command ${ONIKA_RUN_WRAPPER} "${MSP_FILE}" ${EXASTAMP_TEST_ADDITIONAL_ARGS})
        AddSingleTest(${_test_name} ${_test_name}Command 1 ${MT_NTHREADS})
    endif()

    if(ONIKA_REGRESSION_TEST_ENABLE_MPI AND ENABLE_TEST_MPI)
        set(_registered ON)
        set(_test_name ${_base_name}_par)
        set(${_test_name}Command ${ONIKA_RUN_WRAPPER} "${MSP_FILE}" ${EXASTAMP_TEST_ADDITIONAL_ARGS})
        AddSingleTest(${_test_name} ${_test_name}Command 4 ${MT_NTHREADS})
    endif()

    if(TEST_ENABLED AND NOT _registered)
        message(STATUS "Warning: ${_base_name} is enabled but not registered \
(check ONIKA_REGRESSION_TEST_ENABLE_* flags)")
    endif()
endfunction()


# ---------------------------------------------------------------------------
# AddRegressionTestSimple(REGRESSION_DIR)
#
# Scans REGRESSION_DIR for sub-directories (test cases). Inside each
# sub-directory it looks for *.msp files (test variants) and registers them.
#
# Expected layout:
#   REGRESSION_DIR/
#     <testcase>/
#       <variant>.msp
#       ...
#
# CTest name pattern:  <APP>_<folder>_<testcase>_<variant>_{seq|mt|par}
# Filter example:      ctest -R numerical_schemes
# ---------------------------------------------------------------------------
function(AddRegressionTestSimple REGRESSION_DIR)
    # Derive a short prefix from the leaf folder name for test naming/filtering
    get_filename_component(_prefix "${REGRESSION_DIR}" NAME)

    file(GLOB _testcases RELATIVE "${REGRESSION_DIR}" "${REGRESSION_DIR}/*")
    set(_all_variants)

    foreach(_testcase ${_testcases})
        set(_testcase_dir "${REGRESSION_DIR}/${_testcase}")
        if(NOT IS_DIRECTORY "${_testcase_dir}")
            continue()
        endif()

        install(DIRECTORY "${_testcase_dir}" DESTINATION share/examples)

        file(GLOB _msp_files RELATIVE "${_testcase_dir}" "${_testcase_dir}/*.msp")
        foreach(_msp_file ${_msp_files})
            string(REGEX REPLACE "\\.msp$" "" _testvariant "${_msp_file}")
            _RegisterMspTest(
                "${_prefix}"
                "${_testcase}"
                "${_testvariant}"
                "${_testcase_dir}/${_msp_file}"
            )
            list(APPEND _all_variants "${_testcase}_${_testvariant}")
        endforeach()
    endforeach()

    list(LENGTH _all_variants _count)
    message(STATUS "AddRegressionTestSimple: ${_count} test variant(s) found in ${REGRESSION_DIR} [prefix: ${_prefix}]")
endfunction()


# ---------------------------------------------------------------------------
# AddRegressionTestFlat(REGRESSION_DIR)
#
# Like AddRegressionTestSimple() but expects *.msp files directly inside
# REGRESSION_DIR — no sub-directory structure required.
#
# Expected layout:
#   REGRESSION_DIR/
#     <variant>.msp
#     ...
#
# The test-case name is set to the leaf folder name itself so that filtering
# still works naturally:
#   ctest -R flexible_molecule
#
# CTest name pattern:  <APP>_<folder>_<folder>_<variant>_{seq|mt|par}
# ---------------------------------------------------------------------------
function(AddRegressionTestFlat REGRESSION_DIR)
    get_filename_component(_prefix "${REGRESSION_DIR}" NAME)

    install(DIRECTORY "${REGRESSION_DIR}" DESTINATION share/examples)

    file(GLOB _msp_files RELATIVE "${REGRESSION_DIR}" "${REGRESSION_DIR}/*.msp")
    set(_all_variants)

    foreach(_msp_file ${_msp_files})
        string(REGEX REPLACE "\\.msp$" "" _testvariant "${_msp_file}")
        # Use the folder name as both prefix and testcase so the naming stays
        # consistent and `ctest -R <folder>` still selects all tests here.
        _RegisterMspFlatTest(
            "${_prefix}"
            "${_testvariant}"
            "${REGRESSION_DIR}/${_msp_file}"
        )
        list(APPEND _all_variants "${_testvariant}")
    endforeach()

    list(LENGTH _all_variants _count)
    message(STATUS "AddRegressionTestFlat: ${_count} test variant(s) found in ${REGRESSION_DIR} [prefix: ${_prefix}]")
endfunction()

# # Licensed to the Apache Software Foundation (ASF) under one
# # or more contributor license agreements.  See the NOTICE file
# # distributed with this work for additional information
# # regarding copyright ownership.  The ASF licenses this file
# # to you under the Apache License, Version 2.0 (the
# # "License"); you may not use this file except in compliance
# # with the License.  You may obtain a copy of the License at
# # 
# #   http://www.apache.org/licenses/LICENSE-2.0
# # 
# # Unless required by applicable law or agreed to in writing,
# # software distributed under the License is distributed on an
# # "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# # KIND, either express or implied.  See the License for the
# # specific language governing permissions and limitations
# # under the License.
# #
# function(AddRegressionTestSimple REGRESSION_DIR)
#   file(GLOB REGRESSION_TESTS RELATIVE ${REGRESSION_DIR} ${REGRESSION_DIR}/*)
#   set(USTAMP_REGRESSION_TEST_LIST)
#   foreach(testcase ${REGRESSION_TESTS})
#     install(DIRECTORY ${REGRESSION_DIR}/${testcase} DESTINATION share/examples)

#     file(GLOB testruns RELATIVE ${REGRESSION_DIR}/${testcase} ${REGRESSION_DIR}/${testcase}/*.msp)
#     foreach(testrun ${testruns})
    
#       string(REGEX REPLACE ".msp$" "" testvariant ${testrun})

#       # scan file for special keywords on the first line : no-mpi no-seq no-mt
#       file(STRINGS ${REGRESSION_DIR}/${testcase}/${testrun} test_lines)
#       list(GET test_lines 0 test_header)
#       string(FIND "${test_header}" "no-mpi" NO_MPI)
#       string(FIND "${test_header}" "no-seq" NO_SEQ)
#       string(FIND "${test_header}" "no-mt" NO_MT)    
#       if(NO_MPI EQUAL -1)
#         set(ENABLE_TEST_MPI ON)
#       else()
#         set(ENABLE_TEST_MPI OFF)
#       endif()
#       if(NO_SEQ EQUAL -1)
#         set(ENABLE_TEST_SEQ ON)
#       else()
#         set(ENABLE_TEST_SEQ OFF)
#       endif()
#       if(NO_MT EQUAL -1)
#         set(ENABLE_TEST_MT ON)
#       else()
#         set(ENABLE_TEST_MT OFF)
#       endif()

#       if(NOT ONIKA_RUN_WRAPPER)
#         set(ONIKA_RUN_WRAPPER ${ONIKA_RUN})
#       endif()

#       string(REGEX MATCHALL "enable-if-[^ ]*" ALL_COND_VAR "${test_header}")
#       set(${testcase}_${testvariant}_ENABLED ON)
      
#       foreach(COND_VAR ${ALL_COND_VAR})
#         if(COND_VAR)
#           string(REPLACE "enable-if-" "" COND_VAR "${COND_VAR}")
#           string(REGEX MATCH "=.*" COND_VAR_VALUE "${COND_VAR}")
#           if(COND_VAR_VALUE)
#             string(REGEX REPLACE "=.*" "" COND_VAR "${COND_VAR}")
#             string(REPLACE "=" "" COND_VAR_VALUE "${COND_VAR_VALUE}")
#             if(NOT "${${COND_VAR}}" STREQUAL "${COND_VAR_VALUE}" )
#               # message(STATUS "Variable ${COND_VAR} (${${COND_VAR}}) doesn't match '${COND_VAR_VALUE}'")
#               set(XSTAMP_COND_FALSE OFF)
#               set(COND_VAR XSTAMP_COND_FALSE)
#             endif()
#           endif()
          
#           # message(STATUS "${testcase}/${testrun} condition ${COND_VAR} = ${${COND_VAR}}")
#           if(NOT ${COND_VAR})
#             if(${testcase}_${testvariant}_ENABLED)
#               # message(STATUS "disable test ${testcase}_${testvariant} (${COND_VAR}=${${COND_VAR}})")
#               set(ENABLE_TEST_MPI OFF)
#               set(ENABLE_TEST_SEQ OFF)
#               set(ENABLE_TEST_MT OFF)
#               set(${testcase}_${testvariant}_ENABLED OFF)
#             endif()
#           endif()
#         endif()
#       endforeach()

#       set(MT_NTHREADS 4)
#       string(REGEX MATCH "nthreads=[^ ]*" NTHREADS_VAR "${test_header}")
#       if(NTHREADS_VAR)
#         string(REPLACE "nthreads=" "" NTHREADS_VAR "${NTHREADS_VAR}")
#         if("${NTHREADS_VAR}" STREQUAL "max")
#           set(MT_NTHREADS ${ONIKA_HOST_HW_CORES})
#         elseif(${NTHREADS_VAR} GREATER 1)
#           set(MT_NTHREADS ${NTHREADS_VAR})
#         endif()
#         # message(STATUS "${testcase}_${testvariant} : NTHREADS overriden to ${MT_NTHREADS}")
#       endif()

#       # message(STATUS "${testcase}_${testvariant} : SEQ=${ENABLE_TEST_SEQ} MT=${ENABLE_TEST_MT} MPI=${ENABLE_TEST_MPI} ONIKA_ALWAYS_USE_MPIRUN=${ONIKA_ALWAYS_USE_MPIRUN}")
#       list(APPEND USTAMP_REGRESSION_TEST_LIST "${testcase}_${testvariant}")

#       set(${testcase}_${testvariant}_REGISTERED OFF)

#       if(ONIKA_REGRESSION_TEST_ENABLE_SEQ AND ENABLE_TEST_SEQ)
#         set(${testcase}_${testvariant}_REGISTERED ON)
#         set(TestName ${XNB_APP_NAME}_${testcase}_${testvariant}_seq)
#         set(${TestName}Command ${ONIKA_RUN_WRAPPER} ${REGRESSION_DIR}/${testcase}/${testrun} ${EXASTAMP_TEST_ADDITIONAL_ARGS})
#         AddSingleTest(${TestName} ${TestName}Command 1 1)
#       endif()

#       if(ONIKA_REGRESSION_TEST_ENABLE_MT AND ENABLE_TEST_MT)
#         set(${testcase}_${testvariant}_REGISTERED ON)
#         set(TestName ${XNB_APP_NAME}_${testcase}_${testvariant}_mt)
#         set(${TestName}Command ${ONIKA_RUN_WRAPPER} ${REGRESSION_DIR}/${testcase}/${testrun} ${EXASTAMP_TEST_ADDITIONAL_ARGS})
#         AddSingleTest(${TestName} ${TestName}Command 1 ${MT_NTHREADS})
#       endif()

#       if(ONIKA_REGRESSION_TEST_ENABLE_MPI AND ENABLE_TEST_MPI)
#         set(${testcase}_${testvariant}_REGISTERED ON)
#         set(TestName ${XNB_APP_NAME}_${testcase}_${testvariant}_par)
#         set(${TestName}Command ${ONIKA_RUN_WRAPPER} ${REGRESSION_DIR}/${testcase}/${testrun} ${EXASTAMP_TEST_ADDITIONAL_ARGS})
#         AddSingleTest(${TestName} ${TestName}Command 4 ${MT_NTHREADS})
#       endif()
      
#       if(${testcase}_${testvariant}_ENABLED AND NOT ${testcase}_${testvariant}_REGISTERED)
#         message(STATUS "Warning: ${testcase}_${testvariant} enabled but not registered")
#       endif()

#     endforeach()
#   endforeach()

#   list(LENGTH USTAMP_REGRESSION_TEST_LIST USTAMP_REGRESSION_TEST_COUNT)
#   message(STATUS "found ${USTAMP_REGRESSION_TEST_COUNT} regression tests in ${REGRESSION_DIR}")
# endfunction()

