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
export PROJECT_SETUP_ENV_COMMANDS="module purge ; module load inteloneapi/22.1.2 gnu/8.4.0 nvhpc/23.1 mpi/openmpi/4.1.4 texlive cmake/3.20.3"
eval ${PROJECT_SETUP_ENV_COMMANDS}

BUILD_DIR=$CCCSCRATCHDIR/build/exaStamp-cuda
SRC_DIR=$HOME/dev/ExaStamp-main
XNB_DIR=$HOME/dev/exaNBody

mkdir -p $BUILD_DIR
rm -rf $BUILD_DIR/*
cd $BUILD_DIR

# Note:
# remove -Dyaml-cpp_DIR if yamlcpp-dev is installed as a standard package
# # remove compiler specification if default compiler used
cmake -DXNB_PRODUCT_VARIANT=rigidmol \
       -DXNB_BUILD_CUDA=ON \
       -DCMAKE_CUDA_ARCHITECTURES=80 \
       -DONIKA_HAVE_OPENMP_DETACH=OFF \
       -DONIKA_HAVE_OPENMP_TOOLS=OFF \
       -DCMAKE_C_COMPILER=icc \
       -DCMAKE_CXX_COMPILER=icpc \
       -DCMAKE_Fortran_COMPILER=ifort \
       -DCMAKE_CUDA_COMPILER=/ccc/products/cuda-12.0/system/nvhpc-231/bin/nvcc \
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
export PROJECT_SETUP_ENV_COMMANDS="module purge ; module load gnu/11.2.0 nvhpc/23.11 mpi/openmpi texlive cmake"
eval ${PROJECT_SETUP_ENV_COMMANDS}

export CUDA_SDK=/ccc/products/nvhpc-23.11/system/default/Linux_x86_64/23.11/cuda/12.3
unset TBB_ROOT MPI_ROOT
BUILD_DIR=$CCCSCRATCHDIR/build/generic-gcc
SRC_DIR=$HOME/dev/ExaStamp-main
XNB_DIR=$HOME/dev/exaNBody

mkdir -p $BUILD_DIR
rm -rf $BUILD_DIR/*
cd $BUILD_DIR

ccmake -DXNB_PRODUCT_VARIANT=rigidmol \
       -DXNB_BUILD_CUDA=ON \
       -DCMAKE_CUDA_ARCHITECTURES=80 \
       -DONIKA_HAVE_OPENMP_DETACH=OFF \
       -DONIKA_HAVE_OPENMP_TOOLS=OFF \
       -DCUDA_SDK_ROOT_DIR=${CUDA_SDK} \
       -DCMAKE_CUDA_COMPILER=${CUDA_SDK}/bin/nvcc \
       -Dyaml-cpp_DIR=/ccc/home/cont001/xstampdev/xstampdev/tools/yaml-cpp-gcc8.4/lib/cmake/yaml-cpp \
       -DCMAKE_BUILD_TYPE=Release \
       -DXSTAMP_BUILD_exaStampLCHBOP=OFF \
       -DEXASTAMP_TEST_DATA_DIR=${HOME}/data \
       -DPROJECT_SETUP_ENV_COMMANDS="${PROJECT_SETUP_ENV_COMMANDS}" \
       -DexaNBody_DIR=${XNB_DIR} \
       ${SRC_DIR}

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


