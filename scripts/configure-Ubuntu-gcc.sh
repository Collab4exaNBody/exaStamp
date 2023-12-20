#!/bin/bash

source $(dirname $0)/Ubuntu.common

read -r -d '' CMAKE_CONFIG_VARS <<- EOF
  ${CMAKE_CONFIG_VARS}
  exaNBody_DIR=/usr/local/exaNBody-xsp-release-3.3
  CMAKE_CXX_FLAGS = -march=native -mtune=native
  CMAKE_CXX_FLAGS_DEBUG = -g -fsignaling-nans
EOF

EXASTAMP_INSTALL_SUFFIX="-gcc"

source $(dirname $0)/configure.sh

