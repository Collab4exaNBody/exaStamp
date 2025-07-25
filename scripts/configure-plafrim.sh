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
#!/bin/sh

# Global product configuration
PRODUCT_NAME="ExaStamp"
BIN_NAME="xstamp"
TARGET_OS="x86_64"
TARGET_NETWORK="formation.cluster"
DEFAULT_RELEASE_DIR="$HOME/releases"
TEMP_DIR="$HOME/tmp"
DEFAULT_BUILD_DIR="$HOME/build"
#HOST_ALWAYS_USE_MPIRUN=ON

CMAKE_ROOT="/usr"
CMAKE_EXE="cmake3"
CCMAKE_EXE="ccmake3"

#MPI_LAUNCHER_CMD="ccc_mprun -p${TARGET_PARTITION} -n\${NUMPROCS} -c\${NUMTHREADS}"
#MPI_LAUNCHER_CMD="ccc_mprun -n\${NUMPROCS} -p${TARGET_PARTITION} -c\${NUMTHREADS} /bin/env OMP_NUM_THREADS=\${NUMTHREADS}"

read -r -d '' ENV_SETUP_COMMANDS <<- EOF
module purge
module load compiler/gcc/8.3.0 compiler/intel/2019_update4 mpi/openmpi/3.1.4-all
unset TBB_ROOT MPI_ROOT
EOF


read -r -d '' CMAKE_CONFIG_VARS <<- EOF
    CMAKE_CXX_FLAGS_RELEASE     = -O3 -xHOST -DNDEBUG
    EXASTAMP_PARALLEL_MODEL     = OpenMP
    EXASTAMP_USE_HWLOC          = OFF
    EXASTAMP_TEST_DATA_DIR      = /home/cisd-colombe/data
    SOATL_ENABLE_BENCHMARKS     = OFF
    yaml-cpp_DIR                = /home/cisd-colombe/local/yaml-cpp/share/cmake/yaml-cpp
    EXASTAMP_BUILD_SLIPLINK     = OFF
    EXASTAMP_BUILD_SNAP         = OFF
    EXASTAMP_BUILD_TAZ          = OFF
    SOATL_ENABLE_TESTS          = OFF 
EOF

source $(dirname $0)/configure.sh


#    MPIEXEC_PREFLAGS            = -p${TARGET_PARTITION}
#    MPIEXEC_EXECUTABLE          = /usr/bin/ccc_mprun
#    MPIEXEC_NUMCORE_FLAG        = -c
#    MPIEXEC_NUMPROC_FLAG        = -n
#    MPIEXEC_PREFLAGS_DBG        = -p${TARGET_PARTITION};-Xall;xterm;-e
 
