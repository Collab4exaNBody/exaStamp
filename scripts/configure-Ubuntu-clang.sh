#!/bin/bash

source $(dirname $0)/Ubuntu.common

[ $DEFAULT_LLVM_ROOT ] || DEFAULT_LLVM_ROOT=/usr
if [ "$LLVM_ROOT" == "" ]
then
    echo -n "llvm install root ? [$DEFAULT_LLVM_ROOT] "
    read LLVM_ROOT
fi
[ $LLVM_ROOT ] || LLVM_ROOT=$DEFAULT_LLVM_ROOT

[ $DEFAULT_LLVM_LIBDIR ] || [ -d ${LLVM_ROOT}/lib/llvm-11 ] && DEFAULT_LLVM_LIBDIR=${LLVM_ROOT}/lib/llvm-11 || DEFAULT_LLVM_LIBDIR=${LLVM_ROOT}/lib

if [ "$LLVM_LIBDIR" == "" ]
then
    echo -n "llvm lib dir ? [$DEFAULT_LLVM_LIBDIR] "
    read LLVM_LIBDIR
fi
[ $LLVM_LIBDIR ] || LLVM_LIBDIR=$DEFAULT_LLVM_LIBDIR

# Tested option : LLVM-12 , GCC-9 , CUDA-11
GCC_ROOT=`g++ -v 2>&1 | grep "\-\-prefix=" | sed "s/.*--prefix=//g" | cut -d' ' -f1`
GCC_VERSION=`g++ -v 2>&1 | grep "gcc version" | sed "s/gcc version //g" | cut -d' ' -f1`
GCC_VER_MAJOR=`echo $GCC_VERSION|cut -d'.' -f1`
if [ "$GCC_ROOT" == "/usr" ]
then
  GCC_LIBDIR=/usr/lib/gcc/x86_64-linux-gnu/$GCC_VER_MAJOR
else
  GCC_LIBDIR=$GCC_ROOT/lib64/lib
  LLVM_GCC_TOOLCHAIN_OPT="--gcc-toolchain=$GCC_ROOT"
fi

#echo "Using gcc $GCC_VERSION, root=$GCC_ROOT, libdir=$GCC_LIBDIR"

read -r -d '' ENV_SETUP_COMMANDS <<- EOF
export PATH=${LLVM_ROOT}/bin:${PATH}
export LD_LIBRARY_PATH=${LLVM_LIBDIR}:${GCC_LIBDIR}
export GCC_ROOT=$GCC_ROOT
export GCC_LIBDIR=$GCC_LIBDIR
export LLVM_GCC_TOOLCHAIN_OPT=$LLVM_GCC_TOOLCHAIN_OPT
EOF

# pre-run environment setup commands to use defined env vars in cmake vars definition
PRE_RUN_ENV_COMMANDS=`tr "\n" ";"<<<"${ENV_SETUP_COMMANDS}"`
echo "Bootstrap environment ..."
eval ${PRE_RUN_ENV_COMMANDS} 2>&1 >/dev/null

echo "--- clang configuration ---"
echo "PATH            : ${PATH}"
echo "LD_LIBRARY_PATH : ${LD_LIBRARY_PATH}"
echo "LLVM_ROOT       : ${LLVM_ROOT}"
echo "LLVM_LIBDIR     : ${LLVM_LIBDIR}"
echo "GCC_ROOT        : ${GCC_ROOT}"
echo "GCC_LIBDIR      : ${GCC_LIBDIR}"
echo "---------------------------"

read -r -d '' CMAKE_CONFIG_VARS <<- EOF
  ${CMAKE_CONFIG_VARS}
  CMAKE_C_COMPILER          = ${LLVM_ROOT}/bin/clang
  CMAKE_CXX_COMPILER        = ${LLVM_ROOT}/bin/clang++
  CMAKE_CXX_FLAGS           = -march=native -mtune=native $LLVM_GCC_TOOLCHAIN_OPT
  CMAKE_CUDA_FLAGS          = -allow-unsupported-compiler
  CMAKE_INSTALL_RPATH       = ${LLVM_LIBDIR}:${GCC_LIBDIR}
  CMAKE_EXE_LINKER_FLAGS    = -Wl,-rpath,${LLVM_LIBDIR} -Wl,-rpath,${GCC_LIBDIR}
  CMAKE_SHARED_LINKER_FLAGS = -Wl,-rpath,${LLVM_LIBDIR} -Wl,-rpath,${GCC_LIBDIR}
EOF

EXASTAMP_INSTALL_SUFFIX="-clang"

source $(dirname $0)/configure.sh

