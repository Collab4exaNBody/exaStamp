# exaStamp snap is available only if contribs/md/snap has been built in exaNBody

if(XNB_BUILD_mdSnap)
  message(STATUS "exaNBody has snap")
  set(exaStampSnap_LINK_LIBRARIES mdSnap)
  onika_add_plugin(exaStampSnap ${CMAKE_CURRENT_SOURCE_DIR})
endif()
