#!/bin/bash

# Global product configuration
PRODUCT_NAME="ExaStamp"
BIN_NAME="xstamp"

set +x

### Possible compilers on topaze :
# rocm/4.2.0
# llvm/11.0.0
# aocc-compiler/3.0.0
# inteloneapi/21.4.0 (default)
DEFAULT_TARGET_COMPILER="intel"
if [ "$TARGET_COMPILER" == "" ]
then
  echo -n "Compiler (intel,llvm,rocm,aocc) ? [$DEFAULT_TARGET_COMPILER] "
  read TARGET_COMPILER
  [ $TARGET_COMPILER ] || TARGET_COMPILER=${DEFAULT_TARGET_COMPILER}
fi

# Available MPI implementations
# mpi/openmpi/4.0.5.2 , mpi/intelmpi/21.4.0
DEFAULT_MPI_MODULE="intelmpi"
if [ "$MPI_MODULE" == "" ]
then
  echo -n "MPI (intelmpi,openmpi) ? [$DEFAULT_MPI_MODULE] "
  read MPI_MODULE
  [ $MPI_MODULE ] || MPI_MODULE=${DEFAULT_MPI_MODULE}
fi
case $MPI_MODULE in
  intelmpi) MPI_MODULE="mpi/intelmpi/21.4.0" ;;
  openmpi) MPI_MODULE="mpi/openmpi/4.0.5.2" ;;
  *) echo "Unknown '${MPI_MODULE}' MPI module"; exit 1 ;;
esac

#echo "TARGET_COMPILER=$TARGET_COMPILER"
TARGET_OS=`ccc_os`
TARGET_NETWORK="c-topaze.mg1.ccrt.ccc.cea.fr"
TARGET_LOGIN_NODE="topaze.ccc.cea.fr"
XSNVCC=/ccc/products/cuda-11.5/system/nvhpc-221/bin/nvcc
XSNVARCH=80
XSNVMODULE="nvhpc/22.1"
COMPILER_GNU_BACKEND="gnu/11.1.0"
DEFAULT_TARGET_PARTITION="a100"
case $TARGET_COMPILER in
  intel)
    COMPILER_MODULE="inteloneapi/21.4.0" 
    COMPILER_LIBDIRS="/ccc/products/icx-21.4.0/system/default/21.4.0/compiler/lib/intel64_lin:/ccc/products/icx-21.4.0/system/default/21.4.0/lib"
    COMPILER_FORTRAN_EXE="/ccc/products/ifx-21.4.0/system/default/21.4.0/bin/intel64/ifort"
    COMPILER_CC_EXE="/ccc/products/icx-21.4.0/system/default/21.4.0/bin/intel64/icc"
    COMPILER_CXX_EXE="/ccc/products/icx-21.4.0/system/default/21.4.0/bin/intel64/icpc"
    COMPILER_ARCH_FLAGS="-xHOST"
    YAML_CPP_ROOT=/ccc/cont002/home/exadem22/exadem22/tools/yaml-cpp
    ;;
  llvm)
    COMPILER_MODULE="llvm/12.0.0"
    COMPILER_LIBDIRS="/ccc/products/llvm-12.0.0/system/default/lib"
    COMPILER_FORTRAN_EXE="/ccc/products/gcc-11.1.0/system/default/bin/gfortran"
    COMPILER_CC_EXE="/ccc/products/llvm-12.0.0/system/default/bin/clang"
    COMPILER_CXX_EXE="/ccc/products/llvm-12.0.0/system/default/bin/clang++"
    COMPILER_ARCH_FLAGS="-march=native -mtune=native"
    YAML_CPP_ROOT=/ccc/cont002/home/exadem22/exadem22/tools/yaml-cpp-llvm12
    ;;
  *) echo "${TARGET_COMPILER} not supported yet" ; exit 1 ;;
esac

DEFAULT_RELEASE_DIR="/ccc/cont002/home/exadem22/exadem22/releases"
HOST_ALWAYS_USE_MPIRUN=ON

#PRINT_ENV_COMMANDS=1
#CMAKE_VERBOSE_MAKEFILE=ON

if [ "$TARGET_PARTITION" == "" ]
then
    echo -n "Partition ? [$DEFAULT_TARGET_PARTITION] "
    read TARGET_PARTITION
fi
[ $TARGET_PARTITION ] || TARGET_PARTITION=$DEFAULT_TARGET_PARTITION
[ "$TARGET_PARTITION" == "knl" ] && TARGET_OS="Atos_7__mic"

echo "********* Host configuration *********"
echo "TARGET_OS         ${TARGET_OS}"
echo "TARGET_NETWORK    ${TARGET_NETWORK}"
echo "TARGET_LOGIN_NODE ${TARGET_LOGIN_NODE}"
echo "XSNVCC            ${XSNVCC}"
echo "TARGET_PARTITION  ${TARGET_PARTITION}"
echo "TARGET_COMPILER   ${TARGET_COMPILER}"
echo "YAML_CPP_ROOT     ${YAML_CPP_ROOT}"
echo "**************************************"

MPI_LAUNCHER_CMD="ccc_mprun -p${TARGET_PARTITION} -n\${NUMPROCS} -c\${NUMTHREADS}"
MPI_PREALLOCATE_CMD="ccc_mprun -p${TARGET_PARTITION} -n\${NUMPROCS} -c\${NUMTHREADS} -K"
EXASTAMP_INSTALL_SUFFIX="-${TARGET_COMPILER}"

read -r -d '' ENV_SETUP_COMMANDS <<- EOF
module purge
module load ${COMPILER_MODULE} ${COMPILER_GNU_BACKEND} ${MPI_MODULE} ${XSNVMODULE} texlive cmake/3.20.3
export COMPILER_LIBDIRS=${COMPILER_LIBDIRS}
export COMPILER_FORTRAN_EXE=${COMPILER_FORTRAN_EXE}
export COMPILER_CC_EXE=${COMPILER_CC_EXE}
export COMPILER_CXX_EXE=${COMPILER_CXX_EXE}
export XSNVCC=${XSNVCC}
export XSNVARCH=${XSNVARCH}
export YAML_CPP_ROOT=${YAML_CPP_ROOT}
export PRODUCT_NAME=exaDEM
unset TBB_ROOT MPI_ROOT
EOF

# pre-run environment setup commands to use defined env vars in cmake vars definition
PRE_RUN_ENV_COMMANDS=`tr "\n" ";"<<<"${ENV_SETUP_COMMANDS}"`
echo "Bootstrap environment ..."
eval ${PRE_RUN_ENV_COMMANDS} 2>&1 >/dev/null

COMPILER_LIBS_OPT=`echo ":${COMPILER_LIBDIRS}:${CXX_GNU_ROOT}/lib64"|sed "s/:/ -Wl,-rpath,/g"`
echo "COMPILER_LIBS_OPT=${COMPILER_LIBS_OPT}"

read -r -d '' CMAKE_CONFIG_VARS <<- EOF
    CMAKE_CXX_FLAGS_RELEASE     = -O3 ${COMPILER_ARCH_FLAGS} -DNDEBUG
    CMAKE_CXX_COMPILER          = ${COMPILER_CXX_EXE}
    CMAKE_C_COMPILER            = ${COMPILER_CC_EXE}
    CMAKE_Fortran_COMPILER      = ${COMPILER_FORTRAN_EXE}
    CMAKE_INSTALL_RPATH         = ${COMPILER_LIBDIRS}
    CMAKE_EXE_LINKER_FLAGS      = ${COMPILER_LIBS_OPT}
    CMAKE_SHARED_LINKER_FLAGS   = ${COMPILER_LIBS_OPT}
    MPIEXEC_MAX_NUMPROCS        = 256
    MPIEXEC_PREFLAGS            = -p${TARGET_PARTITION}
    MPIEXEC_EXECUTABLE          = /usr/bin/ccc_mprun
    MPIEXEC_NUMCORE_FLAG        = -c
    MPIEXEC_NUMPROC_FLAG        = -n
    MPIEXEC_PREFLAGS_DBG        = -p${TARGET_PARTITION};-Xall;xterm;-e
    yaml-cpp_DIR                = ${YAML_CPP_ROOT}/share/cmake/yaml-cpp
    EXASTAMP_TEST_DATA_DIR      = /ccc/store/cont002/exadem22/exadem22/data
    SOATL_ENABLE_BENCHMARKS     = OFF
    CMAKE_CUDA_COMPILER         = ${XSNVCC}
    XNB_BUILD_CUDA           = ON
    CMAKE_CUDA_ARCHITECTURES    = ${XSNVARCH}
    EXASTAMP_BUILD_TAZ          = OFF
    XSTAMP_ENABLE_SLIPLINK      = OFF
    XSTAMP_BUILD_microStampLCHBOP = OFF
    XSTAMP_ENABLE_SlipLinkFieldSet = OFF
    ONIKA_HAVE_OPENMP_DETACH    = OFF
    ONIKA_HAVE_OPENMP_TOOLS     = OFF
EOF

#    ONIKA_ENABLE_TASK_PROFILING = OFF

source $(dirname $0)/configure.sh

