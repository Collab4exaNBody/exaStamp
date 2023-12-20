#!/bin/bash

source $(dirname $0)/Ubuntu.common

read -r -d '' ENV_SETUP_COMMANDS <<- EOF
${ENV_SETUP_COMMANDS}
export COMPILERVARS_ARCHITECTURE=intel64
export COMPILERVARS_PLATFORM=linux
EOF

read -r -d '' CMAKE_CONFIG_VARS <<- EOF
  ${CMAKE_CONFIG_VARS}
  CMAKE_CXX_COMPILER        = /opt/intel/bin/icpc
  CMAKE_Fortran_COMPILER    = /opt/intel/bin/ifort
  CMAKE_EXE_LINKER_FLAGS    = -Wl,-rpath,/opt/intel/lib/intel64_lin
  CMAKE_MODULE_LINKER_FLAGS = -Wl,-rpath,/opt/intel/lib/intel64_lin
  CMAKE_SHARED_LINKER_FLAGS = -Wl,-rpath,/opt/intel/lib/intel64_lin
  CMAKE_CXX_FLAGS           = -gcc-name=gcc-8 -gxx-name=g++-8 -xHOST
EOF

EXASTAMP_INSTALL_SUFFIX="-intel"

source $(dirname $0)/configure.sh

