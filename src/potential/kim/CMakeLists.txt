# ===================================
# === micro-stamp meam potential  ===
# ===================================

#find_package(KIM-API 2.2 REQUIRED CONFIG)

set(ONIKA_BUILD_exaStampKIM_DEFAULT OFF)

#set(exaStampKIM_LINK_LIBRARIES exaStampPotential KIM-API::kim-api)
set(exaStampKIM_LINK_LIBRARIES exaStampPotential)
onika_add_plugin(exaStampKIM ${CMAKE_CURRENT_SOURCE_DIR})

