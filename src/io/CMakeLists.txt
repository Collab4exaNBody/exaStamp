#
# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements. See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership. The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.
#
# =================================
# === micro-stamp core library  ===
# =================================
# Function to find a package and define a compile macro if found. EXASTAMP_<macro> (cache)
# tells the regression decks whether the codec is there (# enable-if-EXASTAMP_<macro>).
function(find_and_enable_compression target package macro enabled)
    set(EXASTAMP_${macro} OFF CACHE INTERNAL "")
    if(${enabled})
        find_package(${package})
        # Dynamically build the found variable name, e.g., ZLIB_FOUND
        string(TOUPPER "${package}" PACKAGE_UPPER)
        set(found_var "${PACKAGE_UPPER}_FOUND")

        if(${${found_var}})
            # Dynamically get _LIBRARIES and _INCLUDE_DIRS
            set(lib_var     "${PACKAGE_UPPER}_LIBRARIES")
            set(include_var "${PACKAGE_UPPER}_INCLUDE_DIRS")

            message(STATUS "${package} found — enabling ${macro}")

            target_link_libraries(${target} INTERFACE "${${lib_var}}")
            target_include_directories(${target} INTERFACE "${${include_var}}")
            target_compile_definitions(${target} INTERFACE "${macro}=1")
            set(EXASTAMP_${macro} ON CACHE INTERNAL "")
        else()
            message(STATUS "${package} not found — ${macro} disabled")
        endif()
    else()
        message(STATUS "${package} support disabled by user (${enabled}=OFF)")
    endif()
endfunction()

# interface for handling compression libraries
add_library(libZIP INTERFACE)

option(WITH_ZLIB   "Enable ZLIB compression support"    ON)
option(WITH_BZIP2  "Enable BZip2 compression support"   ON)
option(WITH_LZMA   "Enable LZMA/XZ compression support" ON)

# aio.hpp reads AIO_WITH_*: a codec not found gives a clear error on such files
find_and_enable_compression(libZIP ZLIB    AIO_WITH_ZLIB  WITH_ZLIB )
find_and_enable_compression(libZIP BZip2   AIO_WITH_BZIP2 WITH_BZIP2)
find_and_enable_compression(libZIP LibLZMA AIO_WITH_LZMA  WITH_LZMA )


# HDF5 input (read_hdf5, aio ADR 0016): aio.hpp reads AIO_WITH_HDF5. A parallel
# HDF5 (H5_HAVE_PARALLEL) enables MPI-IO reading; pick one with -DHDF5_ROOT=<prefix>.
option(WITH_HDF5 "Enable HDF5 input (read_hdf5)" ON)
set(EXASTAMP_AIO_WITH_HDF5 OFF CACHE INTERNAL "")
if(WITH_HDF5)
    find_package(HDF5 COMPONENTS C)
    if(HDF5_FOUND)
        target_link_libraries(libZIP INTERFACE HDF5::HDF5)
        target_compile_definitions(libZIP INTERFACE AIO_WITH_HDF5=1)
        set(EXASTAMP_AIO_WITH_HDF5 ON CACHE INTERNAL "")
        # h5pcc -show (HDF5 1.14 autotools) gives no -I: without this, another hdf5.h
        # on the system path (e.g. a serial /usr/include one) would be used
        if(NOT HDF5_INCLUDE_DIRS AND HDF5_C_LIBRARY_hdf5)
            get_filename_component(xsp_hdf5_libdir "${HDF5_C_LIBRARY_hdf5}" DIRECTORY)
            target_include_directories(libZIP INTERFACE "${xsp_hdf5_libdir}/../include")
        endif()
        message(STATUS "HDF5 ${HDF5_VERSION} found (parallel: ${HDF5_IS_PARALLEL}) — enabling AIO_WITH_HDF5")
    else()
        message(STATUS "HDF5 not found — read_hdf5 disabled")
    endif()
else()
    message(STATUS "HDF5 support disabled by user (WITH_HDF5=OFF)")
endif()

set(exaStampIO_LINK_LIBRARIES exaStampCompute exaStampThermoState exanbMPI exanbIO exanbDefBox exaStampParticleSpecies exanbGridCellParticles exanbParticleNeighbors libZIP stdc++fs)
set(exaStampIO_INCLUDE_DIRS ${BASEN_INCLUDE_DIRS})

onika_add_plugin(exaStampIO ${CMAKE_CURRENT_SOURCE_DIR})
