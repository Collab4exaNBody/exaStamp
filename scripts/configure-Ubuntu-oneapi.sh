#!/bin/bash

source $(dirname $0)/Ubuntu.common

read -r -d '' ENV_SETUP_COMMANDS <<- EOF
${ENV_SETUP_COMMANDS}
export COMPILERVARS_ARCHITECTURE=intel64
export COMPILERVARS_PLATFORM=linux
source /opt/intel/oneapi/setvars.sh
EOF

read -r -d '' CMAKE_CONFIG_VARS <<- EOF
  ${CMAKE_CONFIG_VARS}
  CMAKE_CXX_COMPILER        = /opt/intel/oneapi/compiler/latest/linux/bin/icpx
  CMAKE_EXE_LINKER_FLAGS    = -Wl,-rpath,/opt/intel/oneapi/compiler/latest/linux/lib -Wl,-rpath,/opt/intel/lib/intel64_lin
  CMAKE_MODULE_LINKER_FLAGS = -Wl,-rpath,/opt/intel/oneapi/compiler/latest/linux/lib -Wl,-rpath,/opt/intel/lib/intel64_lin
  CMAKE_SHARED_LINKER_FLAGS = -Wl,-rpath,/opt/intel/oneapi/compiler/latest/linux/lib -Wl,-rpath,/opt/intel/lib/intel64_lin
EOF

EXASTAMP_INSTALL_SUFFIX="-oneapi"

source $(dirname $0)/configure.sh

