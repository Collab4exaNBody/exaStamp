#!/bin/sh

### Possible compilers on inti :
# rocm/4.2.0
# llvm/11.0.0
# aocc-compiler/3.0.0
# inteloneapi/21.4.0 (default)
# gnu/11.2.0
DEFAULT_TARGET_COMPILER="intel"
if [ "$TARGET_COMPILER" == "" ]
then
  echo -n "Compiler (intel,gcc) ? [$DEFAULT_TARGET_COMPILER] "
  read TARGET_COMPILER
  [ $TARGET_COMPILER ] || TARGET_COMPILER=${DEFAULT_TARGET_COMPILER}
fi

# Available MPI implementations
# mpi/openmpi/4.0.5.2 , mpi/intelmpi/21.4.0, mpi/openmpi/4.1.4
DEFAULT_MPI_MODULE="openmpi"
if [ "$MPI_MODULE" == "" ]
then
  echo -n "MPI (intelmpi,openmpi,latest) ? [$DEFAULT_MPI_MODULE] "
  read MPI_MODULE
  [ $MPI_MODULE ] || MPI_MODULE=${DEFAULT_MPI_MODULE}
fi
case $MPI_MODULE in
  intelmpi) MPI_MODULE="mpi/intelmpi/21.4.0" ;;
  openmpi) MPI_MODULE="mpi/openmpi/4.1.4" ;;
  latest) MPI_MODULE="mpi/openmpi/4.1.4" ;;
  *) echo "Unknown '${MPI_MODULE}' MPI module"; exit 1 ;;
esac

DEFAULT_TOOLS_DIR="/ccc/home/cont001/xstampdev/xstampdev/tools"
if [ "$TOOLS_DIR" == "" ]
then
  echo -n "Tools directory ? [$DEFAULT_TOOLS_DIR] "
  read TOOLS_DIR
  [ $TOOLS_DIR ] || TOOLS_DIR=${DEFAULT_TOOLS_DIR}
fi

DEFAULT_DATA_DIR="/ccc/home/cont001/xstampdev/xstampdev/data"
if [ "$DATA_DIR" == "" ]
then
  echo -n "Data directory ? [$DEFAULT_DATA_DIR] "
  read DATA_DIR
  [ $DATA_DIR ] || DATA_DIR=${DEFAULT_DATA_DIR}
fi

#echo "TARGET_COMPILER=$TARGET_COMPILER"
TARGET_OS=`ccc_os`
case $TARGET_OS in
Rhel_8__x86_64)
  TARGET_NETWORK="c-inti.mg1.ccc.ocre.cea.fr"
  TARGET_LOGIN_NODE="inti-amd.ocre.cea.fr"

#/ccc/products/cuda-11.8/system/nvhpc-2211/bin/
  XSNVMODULE="nvhpc/22.11"
  XSNVCC=/ccc/products/cuda-11.8/system/nvhpc-2211/bin/nvcc
  XSNVARCH=80
  DEFAULT_TARGET_PARTITION="a100-bxi"
  case $TARGET_COMPILER in
    intel)
      COMPILER_MODULE="inteloneapi/22.1.2"
      COMPILER_LIBDIRS="/ccc/products/icx-22.1.2/system/default/22.1.2/compiler/lib/intel64_lin:/ccc/products/icx-22.1.2/system/default/22.1.2/lib"
      COMPILER_FORTRAN_EXE="/ccc/products/ifx-22.1.2/system/default/22.1.2/bin/intel64/ifort"
      COMPILER_CC_EXE="/ccc/products/icx-22.1.2/system/default/22.1.2/bin/intel64/icc"
      COMPILER_CXX_EXE="/ccc/products/icx-22.1.2/system/default/22.1.2/bin/intel64/icpc"
      COMPILER_ARCH_FLAGS="-xHOST"
      COMPILER_GNU_BACKEND="gnu/8.4.0"
      YAML_CPP_ROOT=${TOOLS_DIR}/yaml-cpp-intl21.4-gcc8.4
      ;;
    gcc)
      COMPILER_MODULE="gnu/11"
      COMPILER_LIBDIRS="/ccc/products/gcc-11.2.0/system/default/lib64:/ccc/products/gcc-11.2.0/system/default/lib/gcc/x86_64-pc-linux-gnu/11.2.0"
      COMPILER_FORTRAN_EXE="/ccc/products/gcc-11.2.0/system/default/bin/gfortran"
      COMPILER_CC_EXE="/ccc/products/gcc-11.2.0/system/default/bin/gcc"
      COMPILER_CXX_EXE="/ccc/products/gcc-11.2.0/system/default/bin/g++"
      COMPILER_ARCH_FLAGS="-march=native"
      COMPILER_GNU_BACKEND="gnu/11.2.0"
      YAML_CPP_ROOT=${TOOLS_DIR}/yaml-cpp-gcc8.4
      ;;
    *) echo "${TARGET_COMPILER} not supported yet" ; exit 1 ;;
  esac
  ;;
*)
  echo "Unsupported OS ${TARGET_OS}"
  exit 1
  ;;
esac

DEFAULT_RELEASE_DIR="/ccc/home/cont001/xstampdev/xstampdev/releases"
HOST_ALWAYS_USE_MPIRUN=ON

#PRINT_ENV_COMMANDS=1
#CMAKE_VERBOSE_MAKEFILE=ON

if [ "$TARGET_PARTITION" == "" ]
then
    echo -n "Partition ? [$DEFAULT_TARGET_PARTITION] "
    read TARGET_PARTITION
    [ $TARGET_PARTITION ] || TARGET_PARTITION=$DEFAULT_TARGET_PARTITION
fi

# replace the -xHost intel compiler flag with the appropriate architecture flag for AMD zen2/3 arch (at least +30% perf...)   
if [ "$TARGET_PARTITION" != "sklc" ] && [ "$TARGET_PARTITION" != "sklb" ] ;
then
    case $TARGET_COMPILER in
	intel)
	    COMPILER_ARCH_FLAGS="-march=core-avx2"
	    ;;
	*)
	    COMPILER_ARCH_FLAGS=${COMPILER_ARCH_FLAGS}
	    ;;
    esac
fi

echo "********* Host configuration *********"
echo "TARGET_OS         ${TARGET_OS}"
echo "TARGET_NETWORK    ${TARGET_NETWORK}"
echo "TARGET_LOGIN_NODE ${TARGET_LOGIN_NODE}"
echo "XSNVCC            ${XSNVCC}"
echo "TARGET_PARTITION  ${TARGET_PARTITION}"
echo "TARGET_COMPILER   ${TARGET_COMPILER}"
echo "**************************************"

MPI_LAUNCHER_CMD="ccc_mprun -p${TARGET_PARTITION} -n\${NUMPROCS} -c\${NUMTHREADS}"
EXASTAMP_INSTALL_SUFFIX="-${TARGET_COMPILER}"

read -r -d '' ENV_SETUP_COMMANDS <<- EOF
module purge
module load ${COMPILER_MODULE} ${COMPILER_GNU_BACKEND} ${XSNVMODULE} ${MPI_MODULE} texlive cmake/3.20.3
export COMPILER_LIBDIRS=${COMPILER_LIBDIRS}
export COMPILER_FORTRAN_EXE=${COMPILER_FORTRAN_EXE}
export COMPILER_CC_EXE=${COMPILER_CC_EXE}
export COMPILER_CXX_EXE=${COMPILER_CXX_EXE}
export XSNVCC=${XSNVCC}
export XSNVARCH=${XSNVARCH}
export YAML_CPP_ROOT=${YAML_CPP_ROOT}
unset TBB_ROOT MPI_ROOT
EOF

# pre-run environment setup commands to use defined env vars in cmake vars definition
PRE_RUN_ENV_COMMANDS=`tr "\n" ";"<<<"${ENV_SETUP_COMMANDS}"`
echo "Bootstrap environment ..."
eval ${PRE_RUN_ENV_COMMANDS} 2>&1 >/dev/null

COMPILER_LIBS_OPT=`echo ":${COMPILER_LIBDIRS}:${CXX_GNU_ROOT}/lib64"|sed "s/:/ -Wl,-rpath,/g"`
#echo "COMPILER_LIBDIRS=${COMPILER_LIBDIRS} , CXX_GNU_ROOT=${CXX_GNU_ROOT} , COMPILER_LIBS_OPT=${COMPILER_LIBS_OPT}"

read -r -d '' CMAKE_CONFIG_VARS <<- EOF
    CMAKE_CXX_FLAGS_RELEASE     = -O3 ${COMPILER_ARCH_FLAGS} -DNDEBUG
    CMAKE_CXX_COMPILER          = ${COMPILER_CXX_EXE}
    CMAKE_C_COMPILER            = ${COMPILER_CC_EXE}
    CMAKE_Fortran_COMPILER      = ${COMPILER_FORTRAN_EXE}
    CMAKE_INSTALL_RPATH         = ${COMPILER_LIBDIRS}
    CMAKE_EXE_LINKER_FLAGS      = ${COMPILER_LIBS_OPT}
    CMAKE_SHARED_LINKER_FLAGS   = ${COMPILER_LIBS_OPT}
    MPIEXEC_MAX_NUMPROCS        = 32
    MPIEXEC_PREFLAGS            = -p${TARGET_PARTITION}
    MPIEXEC_EXECUTABLE          = /usr/bin/ccc_mprun
    MPIEXEC_NUMCORE_FLAG        = -c
    MPIEXEC_NUMPROC_FLAG        = -n
    MPIEXEC_PREFLAGS_DBG        = -p${TARGET_PARTITION};-Xall;xterm;-e
    MPIEXEC_PREALLOC_FLAG       = -K
    yaml-cpp_DIR                = ${YAML_CPP_ROOT}/lib/cmake/yaml-cpp
    EXASTAMP_TEST_DATA_DIR      = ${DATA_DIR}
    SOATL_ENABLE_BENCHMARKS     = OFF
    CMAKE_CUDA_COMPILER         = ${XSNVCC}
    XNB_BUILD_CUDA           = ON
    CMAKE_CUDA_ARCHITECTURES    = ${XSNVARCH}
    ONIKA_HAVE_OPENMP_DETACH    = OFF
    ONIKA_HAVE_OPENMP_TOOLS     = OFF
    exaNBody_DIR                = /ccc/home/cont001/xstampdev/xstampdev/releases/exaNBody-v1.1.0
EOF

source $(dirname $0)/configure.sh

