# FindMetis.cmake
#
#  METIS_FOUND - system has METIS
#  METIS_INCLUDE_DIRS - the METIS include directories
#  METIS_LIBRARIES - Link these to use METIS

if(NOT METIS_FOUND)

  if(NOT METIS_DIR)
    set(METIS_DIR "/usr/local/metis" CACHE PATH "METIS install path")
  else()
    if(NOT METIS_INCLUDE_DIRS)
        find_path(METIS_INCLUDE_DIRS parmetis.h PATHS ${METIS_DIR}/include)
    endif()
    if(NOT METIS_LIBRARIES)
      find_library(METIS_LIBRARY NAMES parmetis PATHS ${METIS_DIR}/lib)
      if(METIS_LIBRARY)
        set(METIS_LIBRARIES ${METIS_LIBRARY})
      endif()
    endif()
  endif()

  if(METIS_INCLUDE_DIRS AND METIS_LIBRARIES)
    set(METIS_FOUND TRUE)
  endif()

# message(STATUS "METIS_INCLUDE_DIRS=${METIS_INCLUDE_DIRS}")
# message(STATUS "METIS_LIBRARIES=${METIS_LIBRARIES}")

endif()
