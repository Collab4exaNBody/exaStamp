# licensed to the apache software foundation (asf) under one
# or more contributor license agreements.  see the notice file
# distributed with this work for additional information
# regarding copyright ownership.  the asf licenses this file
# to you under the apache license, version 2.0 (the
# "license"); you may not use this file except in compliance
# with the license.  you may obtain a copy of the license at
# 
#   http://www.apache.org/licenses/license-2.0
# 
# unless required by applicable law or agreed to in writing,
# software distributed under the license is distributed on an
# "as is" basis, without warranties or conditions of any
# kind, either express or implied.  see the license for the
# specific language governing permissions and limitations
# under the license.

# Pair potential generation template

set(PAIR_POTENTIAL_TEMPLATE ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_template)
file(GLOB PAIR_POTENTIAL_TEMPLATE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_template/*)

set_property(GLOBAL PROPERTY GLOBAL_XSTAMP_PAIR_POTENTIALS "")

macro(AddPairPotential PPDIR PluginName)
#  message(STATUS "AddPairPotential(${PPDIR} ${PluginName})")
  get_property(XSTAMP_PAIR_POTENTIALS GLOBAL PROPERTY GLOBAL_XSTAMP_PAIR_POTENTIALS)
  list(APPEND XSTAMP_PAIR_POTENTIALS "${PPDIR}:${PluginName}")
  set_property(GLOBAL PROPERTY GLOBAL_XSTAMP_PAIR_POTENTIALS ${XSTAMP_PAIR_POTENTIALS})
endmacro()

macro(GeneratePairPotential PPDIR PluginName)
#  message(STATUS "GeneratePairPotential(${PPDIR} ${PluginName})")
  if(IS_DIRECTORY ${PPDIR})
    set(PairPotGenDir ${CMAKE_CURRENT_BINARY_DIR}/gen/pair_potentials/${PluginName})
    file(MAKE_DIRECTORY ${PairPotGenDir})
    file(GLOB PP_FILES ${PPDIR}/*)
    foreach(PP_FILE ${PAIR_POTENTIAL_TEMPLATE_SRCS} ${PP_FILES})
      get_filename_component(bname ${PP_FILE} NAME)
      file(CREATE_LINK ${PP_FILE} ${PairPotGenDir}/${bname} SYMBOLIC)
    endforeach()       
#    string(STRIP "${XSTAMP_PAIR_POTENTIAL_DEFINITIONS}" ${PluginName}_COMPILE_DEFINITIONS)
#    message("${PluginName}_COMPILE_DEFINITIONS = '${${PluginName}_COMPILE_DEFINITIONS}'")
    set(${PluginName}_LINK_LIBRARIES ${${PluginName}_LINK_LIBRARIES} exaStampPotential exanbParticleNeighbors)
    onika_add_plugin(${PluginName} ${PairPotGenDir})
  else()
    message(FATAL_ERROR "${PPDIR} is not a directory")
  endif()
endmacro()

macro(GenerateAllPairPotentials)
#  message(STATUS "GenerateAllPairPotentials")
  get_property(XSTAMP_PAIR_POTENTIALS GLOBAL PROPERTY GLOBAL_XSTAMP_PAIR_POTENTIALS)
  foreach(ppdata ${XSTAMP_PAIR_POTENTIALS})
    string(REPLACE ":" ";" PPLIST ${ppdata})
    list(GET PPLIST 0 PPDIR)
    list(GET PPLIST 1 PluginName)
#    message(STATUS "${PPDIR} ${PluginName}")
    GeneratePairPotential(${PPDIR} ${PluginName})
    set(XSTAMP_PAIR_POTENTIAL_NAMES "${XSTAMP_PAIR_POTENTIAL_NAMES} ${PluginName}")
  endforeach()
  if(XSTAMP_PAIR_POTENTIAL_NAMES)
    message(STATUS "found pair potentials :${XSTAMP_PAIR_POTENTIAL_NAMES}")
  else()
    message(STATUS "no pair potential found")
  endif()
endmacro()

