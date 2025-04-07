# =================================
# === micro-stamp core library  ===
# =================================
# find_package(ZLIB REQUIRED)

# set(exaStampIO_LINK_LIBRARIES exaStampMechanical exaStampCompute exanbMPI exanbIO exanbDefBox exaStampParticleSpecies exanbGridCellParticles exanbParticleNeighbors ${ZLIB_LIBRARIES} stdc++fs)
# set(exaStampIO_INCLUDE_DIRS ${ZLIB_INCLUDE_DIR} ${BASEN_INCLUDE_DIRS}) 
# onika_add_plugin(exaStampIO ${CMAKE_CURRENT_SOURCE_DIR})

# Function to find a package and define a compile macro if found
function(find_and_enable_compression target package macro enabled)
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
            target_compile_definitions(${target} INTERFACE "${macro}")
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
option(WITH_BZIP2  "Enable BZip2 compression support"   OFF)
option(WITH_LZMA   "Enable LZMA/XZ compression support" OFF)

find_and_enable_compression(libZIP ZLIB    USE_ZLIB  WITH_ZLIB )
find_and_enable_compression(libZIP BZip2   USE_BZIP2 WITH_BZIP2)
find_and_enable_compression(libZIP LibLZMA USE_LZMA  WITH_LZMA )

set(exaStampIO_LINK_LIBRARIES exaStampMechanical exaStampCompute exanbMPI exanbIO exanbDefBox exaStampParticleSpecies exanbGridCellParticles exanbParticleNeighbors libZIP stdc++fs)
set(exaStampIO_INCLUDE_DIRS ${BASEN_INCLUDE_DIRS}) 
onika_add_plugin(exaStampIO ${CMAKE_CURRENT_SOURCE_DIR})
