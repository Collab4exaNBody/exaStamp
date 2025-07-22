# ===================================
# === micro-stamp meam potential  ===
# ===================================


set(ONIKA_BUILD_exaStampKIM_DEFAULT OFF)

find_package(KIM-API 2.2 REQUIRED CONFIG)
set(exaStampKIM_LINK_LIBRARIES exaStampPotential KIM-API::kim-api)
#set(exaStampKIM_LINK_LIBRARIES exaStampPotential)
onika_add_plugin(exaStampKIM ${CMAKE_CURRENT_SOURCE_DIR})
