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

REQUIRED_PACKAGES="<none>"
source $(dirname $0)/Ubuntu.common

# Global product configuration
EXASTAMP_INSTALL_SUFFIX=""
CMAKE_ROOT="/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/cmake-3.19.6"
GIT="/cea/home/pmcdev/pmcdev/RedHat-6-x86_64/git-2.9.3/bin/git"
DEFAULT_RELEASE_DIR="/cea/home/materio/materio/exaStamp-release"
# FORCE_INTERACTIVE=YES

read -r -d '' ENV_SETUP_COMMANDS <<- EOF
    export PATH=/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/gcc-10.3/bin:/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/gdb-8.3/bin:/usr/cea/ssh/bin:/usr/cea/krb/bin:/usr/lib/jvm/java/bin:/bin:/usr/bin:/sbin:/usr/sbin:/usr/local/sr/bin:/usr/local/libre/bin:/usr/lib64/mpich/bin
    export LD_LIBRARY_PATH=/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/gcc-10.3/lib64
EOF

# pre-run environment setup commands to use defined env vars in cmake vars definition
PRE_RUN_ENV_COMMANDS=`tr "\n" ";"<<<"${ENV_SETUP_COMMANDS}"`
echo "Bootstrap environment ..."
eval ${PRE_RUN_ENV_COMMANDS} 2>&1 >/dev/null

read -r -d '' CMAKE_CONFIG_VARS <<- EOF
  ${CMAKE_CONFIG_VARS}
  CMAKE_C_COMPILER          = /cea/home/pmcdev/pmcdev/CentOS-7-x86_64/gcc-10.3/bin/gcc
  CMAKE_CXX_COMPILER        = /cea/home/pmcdev/pmcdev/CentOS-7-x86_64/gcc-10.3/bin/g++
  CMAKE_EXE_LINKER_FLAGS    = -Wl,-rpath,/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/gcc-10.3/lib64
  GCC_PSTL_TBB_INCDIR       = /cea/home/pmcdev/pmcdev/CentOS-7-x86_64/intel-2019.4.243/tbb/include
  CMAKE_CXX_FLAGS_RELEASE   = -O3 -march=native -mtune=native -DNDEBUG
  EXASTAMP_TEST_DATA_DIR    = /cea/home/materio/materio/exaStamp-release/data
  MPIEXEC_PREFLAGS          =
  MPIEXEC_PREFLAGS_DBG      = xterm;-e
  GIT_EXECUTABLE            = /cea/home/pmcdev/pmcdev/RedHat-6-x86_64/git-2.9.3/bin/git
  yaml-cpp_DIR              = /cea/home/materio/materio/exaStamp-release/tools/yaml/lib/cmake/yaml-cpp
EOF

source $(dirname $0)/configure.sh

