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

# Global product configuration
EXASTAMP_INSTALL_SUFFIX=""
CMAKE_ROOT="/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/cmake-3.19.6"
GIT="/cea/home/pmcdev/pmcdev/RedHat-6-x86_64/git-2.9.3/bin/git"
DEFAULT_RELEASE_DIR="/cea/home/materio/materio/exaStamp-release"
EXASTAMP_INSTALL_SUFFIX="-clang"

# FORCE_INTERACTIVE=YES

read -r -d '' ENV_SETUP_COMMANDS <<- EOF
    export CXX_LLVM_ROOT=/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/llvm-12.0.0
    export CXX_GNU_ROOT=/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/gcc-9.1.0
    export PATH=/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/llvm-12.0.0/bin:/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/gcc-9.1.0/bin:/usr/cea/ssh/bin:/usr/cea/krb/bin:/usr/lib/jvm/java/bin:/bin:/usr/bin:/sbin:/usr/sbin:/usr/local/sr/bin:/usr/local/libre/bin:/usr/lib64/mpich/bin
    export LD_LIBRARY_PATH=/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/llvm-12.0.0/lib:/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/gcc-9.1.0/lib64
EOF

# pre-run environment setup commands to use defined env vars in cmake vars definition
PRE_RUN_ENV_COMMANDS=`tr "\n" ";"<<<"${ENV_SETUP_COMMANDS}"`
echo "Bootstrap environment ..."
eval ${PRE_RUN_ENV_COMMANDS} 2>&1 >/dev/null

echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"

read -r -d '' CMAKE_CONFIG_VARS <<- EOF
  CMAKE_CXX_COMPILER        = ${CXX_LLVM_ROOT}/bin/clang++
  CMAKE_C_COMPILER          = ${CXX_LLVM_ROOT}/bin/clang
  CMAKE_CXX_FLAGS           = --target=x86_64-pc-linux-gnu --gcc-toolchain=/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/gcc-9.1.0
  CMAKE_C_FLAGS             = --target=x86_64-pc-linux-gnu --gcc-toolchain=/cea/home/pmcdev/pmcdev/CentOS-7-x86_64/gcc-9.1.0
  CMAKE_EXE_LINKER_FLAGS    = -Wl,-rpath,${CXX_LLVM_ROOT}/lib -Wl,-rpath,${CXX_GNU_ROOT}/lib64
  CMAKE_SHARED_LINKER_FLAGS = -Wl,-rpath,${CXX_LLVM_ROOT}/lib -Wl,-rpath,${CXX_GNU_ROOT}/lib64
  CMAKE_CXX_FLAGS_RELEASE   = -O3 -march=native -mtune=native -DNDEBUG
  EXASTAMP_TEST_DATA_DIR    = /cea/home/materio/materio/exaStamp-release/data
  MPIEXEC_PREFLAGS          =
  MPIEXEC_PREFLAGS_DBG      = xterm;-e
  GIT_EXECUTABLE            = /cea/home/pmcdev/pmcdev/RedHat-6-x86_64/git-2.9.3/bin/git
  yaml-cpp_DIR              = /cea/home/materio/materio/exaStamp-release/tools/yaml-centos7-clang/lib/cmake/yaml-cpp
  ParaView_DIR              = /cea/home/materio/materio/exaStamp-release/tools/paraview-centos7/lib/cmake/paraview-5.6
  Qt5_DIR                   = /cea/home/pmcdev/pmcdev/RedHat-6-x86_64/qt5/lib/cmake/Qt5
  Qt5Core_DIR               = /cea/home/pmcdev/pmcdev/RedHat-6-x86_64/qt5/lib/cmake/Qt5Core
  Qt5Gui_DIR                = /cea/home/pmcdev/pmcdev/RedHat-6-x86_64/qt5/lib/cmake/Qt5Gui
  Qt5Help_DIR               = /cea/home/pmcdev/pmcdev/RedHat-6-x86_64/qt-5.1.1/lib/cmake/Qt5Help
  Qt5Network_DIR            = /cea/home/pmcdev/pmcdev/RedHat-6-x86_64/qt5/lib/cmake/Qt5Network
  Qt5Widgets_DIR            = /cea/home/pmcdev/pmcdev/RedHat-6-x86_64/qt5/lib/cmake/Qt5Widgets
  Eigen3_DIR                = /cea/home/pmcdev/pmcdev/CentOS-7-x86_64/eigen-3.3.7/share/eigen3/cmake 
EOF

source $(dirname $0)/configure.sh

