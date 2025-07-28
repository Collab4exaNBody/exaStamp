#
# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements. See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership. The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.
#
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

