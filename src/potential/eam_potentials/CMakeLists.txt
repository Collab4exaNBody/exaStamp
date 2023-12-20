# ===========================================
# === micro-stamp EAM potential library  ===
# ===========================================

# Pair potentials generation from template sources
file(GLOB USTAMP_POT_SUBDIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*)
foreach(dir ${USTAMP_POT_SUBDIRS})
  if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${dir})
    AddEAMPotential(${CMAKE_CURRENT_SOURCE_DIR}/${dir} ${dir})
  endif()
endforeach()


