######################## Ubuntu GCC version ##########################
export PROJECT_SETUP_ENV_COMMANDS=""
#eval ${PROJECT_SETUP_ENV_COMMANDS}

export XSNVCC=/usr/local/cuda/bin/nvcc
export XSNVARCH=86
BUILD_DIR=$CCCSCRATCHDIR/build/generic-gcc
SRC_DIR=$HOME/dev/exaStamp-main
XNB_DIR=$HOME/dev/exaNBody

mkdir -p $BUILD_DIR
rm -rf $BUILD_DIR/*
cd $BUILD_DIR

ccmake -DXNB_PRODUCT_VARIANT=rigidmol \
       -DXSTAMP_BUILD_CUDA=ON \
       -DCMAKE_CUDA_ARCHITECTURES=86 \
       -DONIKA_HAVE_OPENMP_DETACH=OFF \
       -DONIKA_HAVE_OPENMP_TOOLS=OFF \
       -DCMAKE_CUDA_COMPILER=${XSNVCC} \
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
       -DXNB_ENABLE_HIP=ON \
       -DONIKA_HAVE_OPENMP_DETACH=OFF \
       -DONIKA_HAVE_OPENMP_TOOLS=OFF \
       -DCMAKE_BUILD_TYPE=Release \
       -DXSTAMP_BUILD_exaStampLCHBOP=OFF \
       -DEXASTAMP_TEST_DATA_DIR=${HOME}/data \
       -DPROJECT_SETUP_ENV_COMMANDS="${PROJECT_SETUP_ENV_COMMANDS}" \
       -DexaNBody_DIR=${XNB_DIR} \
       ${SRC_DIR}

#       -DCMAKE_CXX_FLAGS="-architecture=sm_86" \

exit 0



# OpenMPI (bxi only)
export PROJECT_SETUP_ENV_COMMANDS="module purge ; module load inteloneapi/22.1.2 gnu/8.4.0 nvhpc/23.1 mpi/openmpi/4.1.4 texlive cmake/3.20.3"
eval ${PROJECT_SETUP_ENV_COMMANDS}

export COMPILER_LIBDIRS=/ccc/products/icx-22.1.2/system/default/22.1.2/compiler/lib/intel64_lin:/ccc/products/icx-22.1.2/system/default/22.1.2/lib
export COMPILER_FORTRAN_EXE=/ccc/products/ifx-22.1.2/system/default/22.1.2/bin/intel64/ifort
export COMPILER_CC_EXE=/ccc/products/icx-22.1.2/system/default/22.1.2/bin/intel64/icc
export COMPILER_CXX_EXE=/ccc/products/icx-22.1.2/system/default/22.1.2/bin/intel64/icpc
export XSNVCC=/ccc/products/cuda-12.0/system/nvhpc-231/bin/nvcc
export XSNVARCH=80
export YAML_CPP_ROOT=/ccc/home/cont001/xstampdev/xstampdev/tools/yaml-cpp-intl21.4-gcc8.4
unset TBB_ROOT MPI_ROOT
BUILD_DIR=$CCCSCRATCHDIR/build/generic
SRC_DIR=$HOME/dev/ExaStamp-xsp-release-3.3
XNB_DIR=$HOME/dev/exaNBody

mkdir -p $BUILD_DIR
rm -rf $BUILD_DIR/*
cd $BUILD_DIR

# Note:
# remove -Dyaml-cpp_DIR if yamlcpp-dev is installed as a standard package
# # remove compiler specification if default compiler used
ccmake -DXNB_PRODUCT_VARIANT=rigidmol \
       -DXSTAMP_BUILD_CUDA=ON \
       -DXSTAMP_CUDA_ARCH=80 \
       -DONIKA_HAVE_OPENMP_DETACH=OFF \
       -DONIKA_HAVE_OPENMP_TOOLS=OFF \
       -DCMAKE_C_COMPILER=icc \
       -DCMAKE_CXX_COMPILER=icpc \
       -DCMAKE_Fortran_COMPILER=ifort \
       -DCMAKE_CUDA_COMPILER=/ccc/products/cuda-12.0/system/nvhpc-231/bin/nvcc \
       -Dyaml-cpp_DIR=/ccc/home/cont001/xstampdev/xstampdev/tools/yaml-cpp-gcc8.4/lib/cmake/yaml-cpp \
       -DCMAKE_BUILD_TYPE=Release \
       -DXSTAMP_BUILD_exaStampLCHBOP=OFF \
       -DEXASTAMP_TEST_DATA_DIR=/ccc/home/cont001/xstampdev/xstampdev/data \
       -DPROJECT_SETUP_ENV_COMMANDS="${PROJECT_SETUP_ENV_COMMANDS}" \
       -DexaNBody_DIR=${XNB_DIR} \
       $HOME/dev/ExaStamp-xsp-release-3.3
# compile
echo "Compile with : cd $BUILD_DIR && make -j8 && make UpdatePluginDataBase"

# execute
# OMP_NUM_THREADS=64 ccc_mprun -n4 -c32 -pa100 ./exaStamp myInput.msp
echo "In build dir, run with : OMP_NUM_THREADS=XX mpirun <mpi options> ./exaStamp myInput.msp"
echo "Note: use one mpi process per GPU card, i.e. 4 GPU per node => 4 mpi process per node"
exit 0


######################## GCC version ##########################
export PROJECT_SETUP_ENV_COMMANDS="module purge ; module load gnu/11.2.0 nvhpc/23.11 mpi/openmpi texlive cmake"
eval ${PROJECT_SETUP_ENV_COMMANDS}

export XSNVCC=/ccc/products/nvhpc-23.11/system/default/Linux_x86_64/23.11/cuda/12.3/bin/nvcc
export XSNVARCH=80
export YAML_CPP_ROOT=/ccc/home/cont001/xstampdev/xstampdev/tools/yaml-cpp-intl21.4-gcc8.4
unset TBB_ROOT MPI_ROOT
BUILD_DIR=$CCCSCRATCHDIR/build/generic-gcc
SRC_DIR=$HOME/dev/ExaStamp-xsp-release-3.3
XNB_DIR=$HOME/dev/exaNBody

mkdir -p $BUILD_DIR
rm -rf $BUILD_DIR/*
cd $BUILD_DIR

ccmake -DXNB_PRODUCT_VARIANT=rigidmol \
       -DXSTAMP_BUILD_CUDA=ON \
       -DXSTAMP_CUDA_ARCH=80 \
       -DONIKA_HAVE_OPENMP_DETACH=OFF \
       -DONIKA_HAVE_OPENMP_TOOLS=OFF \
       -DCMAKE_CUDA_COMPILER=${XSNVCC} \
       -Dyaml-cpp_DIR=/ccc/home/cont001/xstampdev/xstampdev/tools/yaml-cpp-gcc8.4/lib/cmake/yaml-cpp \
       -DCMAKE_BUILD_TYPE=Release \
       -DXSTAMP_BUILD_exaStampLCHBOP=OFF \
       -DEXASTAMP_TEST_DATA_DIR=/ccc/home/cont001/xstampdev/xstampdev/data \
       -DPROJECT_SETUP_ENV_COMMANDS="${PROJECT_SETUP_ENV_COMMANDS}" \
       -DexaNBody_DIR=${XNB_DIR} \
       ${SRC_DIR}

exit 0


