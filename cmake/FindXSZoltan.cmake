# FindXSZoltan.cmake
#
#  ZOLTAN_FOUND - system has ZOLTAN
#  ZOLTAN_INCLUDE_DIRS - the ZOLTAN include directories
#  ZOLTAN_LIBRARIES - Link these to use ZOLTAN

if(NOT XSZOLTAN_FOUND)

  if(NOT ZOLTAN_DIR)
    #message(STATUS "ZOLTAN_ROOT='$ENV{ZOLTAN_ROOT}'")
    if(DEFINED ENV{ZOLTAN_ROOT})
      #message(STATUS "found ZOLTAN_ROOT env var, using as ZOLTAN_DIR default value")
      set(ZOLTAN_DIR_DEFAULT "$ENV{ZOLTAN_ROOT}")
    else()
      set(ZOLTAN_DIR_DEFAULT "/usr/local/zoltan")
    endif()
    set(ZOLTAN_DIR ${ZOLTAN_DIR_DEFAULT} CACHE PATH "ZOLTAN install path")
  else()
    if(NOT ZOLTAN_INCLUDE_DIRS)
        find_path(ZOLTAN_INCLUDE_DIRS zoltan.h PATHS ${ZOLTAN_DIR}/include)
    endif()
    if(NOT ZOLTAN_LIBRARIES)
      find_library(ZOLTAN_LIBRARY NAMES zoltan PATHS ${ZOLTAN_DIR}/lib)
      if(ZOLTAN_LIBRARY)
        set(ZOLTAN_LIBRARIES ${ZOLTAN_LIBRARY})
      endif()
    endif()
  endif()

  if(ZOLTAN_INCLUDE_DIRS AND ZOLTAN_LIBRARIES)
    set(XSZOLTAN_FOUND TRUE)
  endif()

 #message("ZOLTAN_INCLUDE_DIRS=${ZOLTAN_INCLUDE_DIRS}")
 #message("ZOLTAN_LIBRARIES=${ZOLTAN_LIBRARIES}")

endif()

