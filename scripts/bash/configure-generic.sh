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
######################## Ubuntu GCC+Cuda version ##########################
export PROJECT_SETUP_ENV_COMMANDS=""
#eval ${PROJECT_SETUP_ENV_COMMANDS}

export CUDA_SDK=/usr/local/cuda
BUILD_DIR=$HOME/build/exaStamp-cuda
SRC_DIR=$HOME/dev/exaStamp-main
XNB_DIR=$HOME/dev/exaNBody

mkdir -p $BUILD_DIR
rm -rf $BUILD_DIR/*
cd $BUILD_DIR

ccmake -DXNB_PRODUCT_VARIANT=rigidmol \
       -DXNB_BUILD_CUDA=ON \
       -DCMAKE_CUDA_ARCHITECTURES=86 \
       -DONIKA_HAVE_OPENMP_DETACH=OFF \
       -DONIKA_HAVE_OPENMP_TOOLS=OFF \
       -DCUDA_SDK_ROOT_DIR=${CUDA_SDK} \
       -DCMAKE_CUDA_COMPILER=${CUDA_SDK}/bin/nvcc \
       -DCMAKE_BUILD_TYPE=Release \
       -DXSTAMP_BUILD_exaStampLCHBOP=OFF \
       -DEXASTAMP_TEST_DATA_DIR=${HOME}/data \
       -DPROJECT_SETUP_ENV_COMMANDS="${PROJECT_SETUP_ENV_COMMANDS}" \
       -DexaNBody_DIR=${XNB_DIR} \
       -DONIKA_CU_MAX_THREADS_PER_BLOCK=256 \
       -DONIKA_CU_MIN_BLOCKS_PER_SM=4 \
       -DONIKA_DEFAULT_ALIGNMENT=64 \
       -DONIKA_DEFAULT_CHUNK_SIZE=8 \
       -DONIKA_MINIMUM_CUDA_ALIGNMENT=64 \
       -DXNB_CHUNK_NBH_DELAYED_COMPUTE_BUFER_SIZE=4 \
       -DXNB_CHUNK_NBH_DELAYED_COMPUTE_CS1=OFF \
       -DXNB_CHUNK_NBH_DELAYED_COMPUTE_MAX_BLOCK_SIZE=64 \
       -DXNB_CHUNK_NEIGHBORS_CS_LIST="4,1" \
       ${SRC_DIR}

exit 0


######################## Ubuntu GCC+HIP version ##########################
PROJECT_SETUP_ENV_COMMANDS=""
#eval ${PROJECT_SETUP_ENV_COMMANDS}

BUILD_DIR=$HOME/build/exaStamp-hip
SRC_DIR=$HOME/dev/exaStamp-main
XNB_DIR=$HOME/dev/exaNBody

mkdir -p $BUILD_DIR
rm -rf $BUILD_DIR/*
cd $BUILD_DIR

ccmake -DXNB_PRODUCT_VARIANT=rigidmol \
       -DCMAKE_C_COMPILER=hipcc \
       -DCMAKE_CXX_COMPILER=hipcc \
       -DCMAKE_CXX_FLAGS="-Wpass-failed" \
       -DXNB_BUILD_CUDA=ON \
       -DXNB_ENABLE_HIP=ON \
       -DCMAKE_HIP_PLATFORM=amd \
       -DONIKA_HAVE_OPENMP_DETACH=OFF \
       -DONIKA_HAVE_OPENMP_TOOLS=OFF \
       -DCMAKE_BUILD_TYPE=Release \
       -DXSTAMP_BUILD_exaStampLCHBOP=OFF \
       -DEXASTAMP_TEST_DATA_DIR=${HOME}/data \
       -DPROJECT_SETUP_ENV_COMMANDS="${PROJECT_SETUP_ENV_COMMANDS}" \
       -DexaNBody_DIR=${XNB_DIR} \
       ${SRC_DIR}

exit 0



# OpenMPI (bxi only)
#export PROJECT_SETUP_ENV_COMMANDS="module purge ; module load inteloneapi/22.1.2 gnu/8.4.0 nvhpc/23.1 mpi/openmpi/4.1.4 texlive cmake/3.20.3"
export PROJECT_SETUP_ENV_COMMANDS="module purge ; module load inteloneapi/22.1.2 gnu/8.4.0 nvhpc/24.3 mpi/openmpi/4.1.4 texlive cmake/3.26.4"
eval ${PROJECT_SETUP_ENV_COMMANDS}

BUILD_DIR=$CCCSCRATCHDIR/build/exaStamp-3.5-mol
SRC_DIR=$HOME/dev/ExaStamp-main
XNB_DIR=$HOME/dev/exaNBody

mkdir -p $BUILD_DIR
rm -rf $BUILD_DIR/*
cd $BUILD_DIR

# Note:
# remove -Dyaml-cpp_DIR if yamlcpp-dev is installed as a standard package
# remove compiler specification if default compiler used
# optional flags to reproduce the same optimization as LAMMPS : (-xHost) -fp-model fast=2 -no-prec-div -qoverride-limits
cmake -DXNB_PRODUCT_VARIANT=rigidmol \
       -DXNB_BUILD_CUDA=ON \
       -DCMAKE_CUDA_ARCHITECTURES=80 \
       -DONIKA_HAVE_OPENMP_DETACH=OFF \
       -DONIKA_HAVE_OPENMP_TOOLS=OFF \
       -DCMAKE_C_COMPILER=icc \
       -DCMAKE_CXX_COMPILER=icpc \
       -DCMAKE_CXX_FLAGS="-diag-disable=15518 -diag-disable=15552" \
       -DCMAKE_Fortran_COMPILER=ifort \
       -DCMAKE_CUDA_COMPILER=/ccc/products/cuda-12.4/system/default/bin/nvcc \
       -DCMAKE_CUDA_FLAGS="-ccbin /ccc/products/icx-22.1.2/system/default/22.1.2/bin/intel64/icpc" \
       -DMPIEXEC_MAX_NUMPROCS=32 \
       -DMPIEXEC_PREFLAGS="-pa100-bxi" \
       -DMPIEXEC_EXECUTABLE="/usr/bin/ccc_mprun" \
       -DMPIEXEC_NUMCORE_FLAG="-c" \
       -DMPIEXEC_NUMPROC_FLAG="-n" \
       -DMPIEXEC_PREFLAGS_DBG="-pa100-bxi;-Xall;xterm;-e" \
       -DMPIEXEC_PREALLOC_FLAG="-K" \
       -DHOST_ALWAYS_USE_MPIRUN=ON \
       -DHOST_HW_CORES=128 \
       -DHOST_HW_THREADS=256 \
       -Dyaml-cpp_DIR=/ccc/home/cont001/xstampdev/xstampdev/tools/yaml-cpp-gcc8.4/lib/cmake/yaml-cpp \
       -DCMAKE_BUILD_TYPE=Release \
       -DXSTAMP_BUILD_exaStampLCHBOP=OFF \
       -DXSTAMP_BUILD_exaStampSnapLMP=ON \
       -DXSTAMP_THIRD_PARTY_TOOLS_ROOT=/ccc/home/cont001/xstampdev/xstampdev/tools \
       -DEXASTAMP_TEST_DATA_DIR=/ccc/home/cont001/xstampdev/xstampdev/data \
       -DPROJECT_SETUP_ENV_COMMANDS="${PROJECT_SETUP_ENV_COMMANDS}" \
       -DexaNBody_DIR=${XNB_DIR} \
       ${SRC_DIR}

cat README.txt

exit 0


######################## GCC+Cuda version ##########################
export PROJECT_SETUP_ENV_COMMANDS="module purge ; module load gnu/11.2.0 nvhpc/24.3 mpi/openmpi texlive cmake/3.26.4"
eval ${PROJECT_SETUP_ENV_COMMANDS}

BUILD_DIR=$CCCSCRATCHDIR/build/exaStamp-gcc
SRC_DIR=$HOME/dev/ExaStamp-main
XNB_DIR=$HOME/dev/exaNBody

mkdir -p $BUILD_DIR
rm -rf $BUILD_DIR/*
cd $BUILD_DIR

cmake -DXNB_PRODUCT_VARIANT=rigidmol \
       -DXNB_BUILD_CUDA=ON \
       -DCMAKE_CUDA_ARCHITECTURES=80 \
       -DONIKA_HAVE_OPENMP_DETACH=OFF \
       -DONIKA_HAVE_OPENMP_TOOLS=OFF \
       -DCMAKE_CUDA_COMPILER=/ccc/products/cuda-12.4/system/default/bin/nvcc \
       -DMPIEXEC_MAX_NUMPROCS=32 \
       -DMPIEXEC_PREFLAGS="-pa100-bxi" \
       -DMPIEXEC_EXECUTABLE="/usr/bin/ccc_mprun" \
       -DMPIEXEC_NUMCORE_FLAG="-c" \
       -DMPIEXEC_NUMPROC_FLAG="-n" \
       -DMPIEXEC_PREFLAGS_DBG="-pa100-bxi;-Xall;xterm;-e" \
       -DMPIEXEC_PREALLOC_FLAG="-K" \
       -DHOST_ALWAYS_USE_MPIRUN=ON \
       -DHOST_HW_CORES=128 \
       -DHOST_HW_THREADS=256 \
       -Dyaml-cpp_DIR=/ccc/home/cont001/xstampdev/xstampdev/tools/yaml-cpp-gcc8.4/lib/cmake/yaml-cpp \
       -DCMAKE_BUILD_TYPE=Release \
       -DXSTAMP_BUILD_exaStampLCHBOP=OFF \
       -DXSTAMP_BUILD_exaStampSnapLMP=ON \
       -DXSTAMP_THIRD_PARTY_TOOLS_ROOT=/ccc/home/cont001/xstampdev/xstampdev/tools \
       -DEXASTAMP_TEST_DATA_DIR=/ccc/home/cont001/xstampdev/xstampdev/data \
       -DPROJECT_SETUP_ENV_COMMANDS="${PROJECT_SETUP_ENV_COMMANDS}" \
       -DexaNBody_DIR=${XNB_DIR} \
       ${SRC_DIR}

cat README.txt

exit 0

######################## GCC+HIP version ##########################
export PROJECT_SETUP_ENV_COMMANDS="module purge ; module load rocm/5.6.0 mpi/intelmpi/21.4.0 cmake/3.26.4"
eval ${PROJECT_SETUP_ENV_COMMANDS}

unset TBB_ROOT MPI_ROOT
BUILD_DIR=$CCCSCRATCHDIR/build/generic-hip
SRC_DIR=$HOME/dev/ExaStamp-main
XNB_DIR=$HOME/dev/exaNBody

mkdir -p $BUILD_DIR
rm -rf $BUILD_DIR/*
cd $BUILD_DIR

ccmake -DXNB_PRODUCT_VARIANT=rigidmol \
       -DCMAKE_C_COMPILER=hipcc \
       -DCMAKE_CXX_COMPILER=hipcc \
       -DCMAKE_CXX_FLAGS="-Wpass-failed" \
       -DXNB_BUILD_CUDA=ON \
       -DXNB_ENABLE_HIP=ON \
       -DCMAKE_HIP_PLATFORM=amd \
       -DCMAKE_HIP_ARCHITECTURES=gfx90a \
       -DOpenMP_CXX_LIB_NAMES="gomp;pthread" \
       -DOpenMP_gomp_LIBRARY=/opt/rocm-5.6.0/llvm/lib/libgomp.so \
       -DHOST_HW_CORES=128 \
       -DHOST_HW_THREADS=256 \
       -DONIKA_HAVE_OPENMP_DETACH=OFF \
       -DONIKA_HAVE_OPENMP_TOOLS=OFF \
       -Dyaml-cpp_DIR=/ccc/home/cont001/xstampdev/xstampdev/tools/yaml-cpp-gcc8.4/lib/cmake/yaml-cpp \
       -DCMAKE_BUILD_TYPE=Release \
       -DXSTAMP_BUILD_exaStampLCHBOP=OFF \
       -DEXASTAMP_TEST_DATA_DIR=${HOME}/data \
       -DPROJECT_SETUP_ENV_COMMANDS="${PROJECT_SETUP_ENV_COMMANDS}" \
       -DexaNBody_DIR=${XNB_DIR} \
       ${SRC_DIR}

# -DCMAKE_C_FLAGS="-march=znver3" -DCMAKE_CXX_FLAGS="-march=znver3" 

exit 0


