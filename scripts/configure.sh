# get source tree from script location
[ $PROJECT_SRC_DIR ] || PROJECT_SRC_DIR="$(readlink -f $(dirname $0)/..)"
[ $DEFAULT_SOURCE_DIR ] || DEFAULT_SOURCE_DIR="${PROJECT_SRC_DIR}"

# sub scripts files
COMMON_FUNC_SCRIPT=${PROJECT_SRC_DIR}/scripts/functions.sh
[ $HOST_INFO_SCRIPT ] || HOST_INFO_SCRIPT=${PROJECT_SRC_DIR}/scripts/host_info.sh
[ $HOST_INFO_CMD ] || HOST_INFO_CMD="run_mpi_job 1 1 ${HOST_INFO_SCRIPT}"
[ -f ${COMMON_FUNC_SCRIPT} ] || ( echo "Error, ${COMMON_FUNC_SCRIPT} not found" && exit 1 )
[ -f ${HOST_INFO_SCRIPT} ] || ( echo "Error, ${HOST_INFO_SCRIPT} not found" && exit 1 )

# load common functions
source ${COMMON_FUNC_SCRIPT}

# general options
[ $CMAKE_MIN_VERSION ] || CMAKE_MIN_VERSION=3.10.0
PRODUCT_NAME="exaStamp"

if [ -d ${PROJECT_SRC_DIR}/addons ]
then
  VARIANTS=`ls ${PROJECT_SRC_DIR}/addons/*.cmake|sed "s/\/.*\/\(.*\).cmake\$/\1/g"|tr "\n" " "`
  DEFVARIANT=`echo $VARIANTS|cut -d' ' -f1|tr -d "\n"`
  echo -n "Build variant (${VARIANTS}) ? [${DEFVARIANT}] "
  read PRODUCT_VARIANT
  [ "${PRODUCT_VARIANT}" == "" ] && PRODUCT_VARIANT=${DEFVARIANT}
fi

# host command execution properties
[ ${HOST_ALWAYS_USE_MPIRUN} ] || HOST_ALWAYS_USE_MPIRUN=OFF

# check for cmake
[ $CMAKE_EXE ] || CMAKE_EXE=cmake
[ $CCMAKE_EXE ] || CCMAKE_EXE=ccmake

if [ "$CMAKE_ROOT" == "" ]
then
        CMAKE_PATH=`which ${CMAKE_EXE}`
        [ $CMAKE_PATH ] && CMAKE_VERSION=`$CMAKE_PATH --version | head -1 | cut -d' ' -f3`
        [ $CMAKE_PATH ] && [ $(version $CMAKE_VERSION) -ge $(version $CMAKE_MIN_VERSION) ] && CMAKE_ROOT=$(dirname $(dirname $CMAKE_PATH))
fi
if [ "$CMAKE_ROOT" == "" ]
then
        echo "${CMAKE_EXE} (version >= $CMAKE_MIN_VERSION) not found"
	echo "please set CMAKE_ROOT to the install prefix of cmake, or add cmake binary dir to your PATH env"
        exit 1
fi

echo "Using cmake here : ${CMAKE_EXE}"

# find Git command
[ $GIT ] || GIT="git"

# set temp directory
CCC_HOME_CMD=`which ccc_home 2>/dev/null`
HOME_CMD=`which home 2>/dev/null`
[ $CCC_HOME_CMD ] && TEMP_DIR=`$CCC_HOME_CMD -s 2>/dev/null`
[ $CCC_HOME_CMD ] && [ "$TEMP_DIR" == "" ] && TEMP_DIR=`$CCC_HOME_CMD -t 2>/dev/null`
[ $HOME_CMD ] && [ "$TEMP_DIR" == "" ] && TEMP_DIR=`$HOME_CMD -t 2>/dev/null`
[ $TEMP_DIR ] || TEMP_DIR="/tmp/exaStampBuild" # TEMP_DIR=`mktemp -d` 

# host machine informations
# run_mpi_job 1 1 DEBUG ${HOST_INFO_SCRIPT}
echo "Scanning host configuration ..."
HOST_INFO=`${HOST_INFO_CMD}`
read HOST_OS HOST_PARTITION NCORES NTHREADS <<<"${HOST_INFO}"
echo "Host configuration : OS=${HOST_OS}, partition=${HOST_PARTITION}, cores=${NCORES}, threads=${NTHREADS}"

# if no target os has been specified, default to host OS
[ ${TARGET_OS} ] || TARGET_OS=${HOST_OS}

# default path to tools and directories
CMAKE="${CMAKE_ROOT}/bin/${CMAKE_EXE}"
CCMAKE="${CMAKE_ROOT}/bin/${CCMAKE_EXE}"
[ $DEFAULT_SOURCE_DIR ] || DEFAULT_SOURCE_DIR=`readlink -f "$(dirname $(dirname $0))"`
PARTITION_PATH_SUFFIX=""
[ ${TARGET_PARTITION} ] && PARTITION_PATH_SUFFIX="/${TARGET_PARTITION}"

if [ "${TARGET_NETWORK}" != "" ]
then
    HOST_NETWORK=`hostname -f|cut -d'.' -f2-`
    if [ "${TARGET_NETWORK}" != "${HOST_NETWORK}" ]
    then
        echo "Error: Must be logged on ${TARGET_LOGIN_NODE} (found network ${HOST_NETWORK} instead). Use ssh ${TARGET_LOGIN_NODE} -l <login>"
        exit 1
    fi
fi

if [ "${TARGET_PARTITION}" != "" ] && [ "${HOST_PARTITION}" != "${TARGET_PARTITION}" ]
then
  echo "Error: Host partition '${HOST_PARTITION}' doest not match target partition '${TARGET_PARTITION}'"
  exit 1
fi

if [ "${TARGET_OS}" != "" ] && [ "${TARGET_OS}" != "${HOST_OS}" ]
then
    echo "Error: Target OS '${HOST_OS}' noes not match '${TARGET_OS}'"
    exit 1
fi

# Ask for source path
[ $SOURCE_DIR ] && DEFAULT_SOURCE_DIR=$SOURCE_DIR
if [ "${NON_INTERACTIVE_CONFIGURE}" == "" ]
then
  echo -n "Source path ? [$DEFAULT_SOURCE_DIR] "
  read SOURCE_DIR
  [ $SOURCE_DIR ] || SOURCE_DIR="$DEFAULT_SOURCE_DIR"
else
  SOURCE_DIR="$DEFAULT_SOURCE_DIR"
fi

# guess build and install suffix if not provided
[ $PRODUCT_SUFFIX ] || PRODUCT_SUFFIX=`(cd $SOURCE_DIR;$GIT branch)|grep "^* "|cut -b3-`
[ $DEFAULT_RELEASE_DIR ] || DEFAULT_RELEASE_DIR="/usr/local/$PRODUCT_NAME"

# detect if this is a release branch
PRODUCT_IS_RELEASE=`echo $PRODUCT_SUFFIX|cut -b5-11`
if [ "$PRODUCT_IS_RELEASE" == "release" ]
then
    PRODUCT_RELEASE_VERSION=`echo $PRODUCT_SUFFIX|cut -b13-`
    echo "Release v$PRODUCT_RELEASE_VERSION detected"
    DEFAULT_INSTALL_DIR="$DEFAULT_RELEASE_DIR/${PRODUCT_RELEASE_VERSION}/${TARGET_OS}${PARTITION_PATH_SUFFIX}"
else
  if [ "$PRODUCT_SUFFIX" == "master" ]
  then
	PRODUCT_SUFFIX=""
  fi
fi
[ "$PRODUCT_SUFFIX" != "" ] && PRODUCT_SUFFIX="-$PRODUCT_SUFFIX"

# build default paths for build and install
[ $DEFAULT_BUILD_DIR ] || DEFAULT_BUILD_DIR="$TEMP_DIR/build/$PRODUCT_NAME${PRODUCT_SUFFIX}/${TARGET_OS}${PARTITION_PATH_SUFFIX}${EXASTAMP_INSTALL_SUFFIX}"
[ $DEFAULT_INSTALL_DIR ] || DEFAULT_INSTALL_DIR="$TEMP_DIR/local/$PRODUCT_NAME${PRODUCT_SUFFIX}/${TARGET_OS}${PARTITION_PATH_SUFFIX}${EXASTAMP_INSTALL_SUFFIX}"

# ask for build path
[ $BUILD_DIR ] && DEFAULT_BUILD_DIR=${BUILD_DIR}
if [ "${NON_INTERACTIVE_CONFIGURE}" == "" ]
then
  echo -n "Build dir ? [$DEFAULT_BUILD_DIR] "
  read BUILD_DIR
  [ $BUILD_DIR ] || BUILD_DIR="$DEFAULT_BUILD_DIR"
else
  BUILD_DIR="$DEFAULT_BUILD_DIR"
fi

# test if build dir already exists
if [ "${NON_INTERACTIVE_CONFIGURE}" == "" ]
then
  [ -d ${BUILD_DIR} ] && EXISTING_BUILD_DIR=true
  if [ ${EXISTING_BUILD_DIR} ]
  then
    echo -n "Build dir already exists, erase it (y/n) ? [n] "
    read ERASE_BUILD_DIR
    [ "${ERASE_BUILD_DIR}" == "y" ] && echo -n "Erase ${BUILD_DIR} ... " && rm -fr ${BUILD_DIR} && EXISTING_BUILD_DIR="" && echo " ok"
  fi
else
  rm -fr $BUILD_DIR
fi

# prepare build dir and change dir to it
echo "Create directory ${BUILD_DIR} ..."
mkdir -p ${BUILD_DIR}
if [ ! -d ${BUILD_DIR} ]
then
	echo "Error: Directory ${BUILD_DIR} could not be created, aborting"
	exit 1
fi

# ask for install path
[ $INSTALL_DIR ] && DEFAULT_INSTALL_DIR=$INSTALL_DIR
if [ "${NON_INTERACTIVE_CONFIGURE}" == "" ]
then
  echo -n "Install prefix ? [$DEFAULT_INSTALL_DIR] "
  read INSTALL_DIR
  [ $INSTALL_DIR ] || INSTALL_DIR="$DEFAULT_INSTALL_DIR"
else
  INSTALL_DIR="$DEFAULT_INSTALL_DIR"
fi

# Preparation des options pour cmake a partir des variables definies dans ${CMAKE_CONFIG_VARS}
CMAKE_CONFIG_OPTS=""
for cmd in `tr " \t\n" "@@ "<<<"${CMAKE_CONFIG_VARS}${CMAKE_VARS_ADDON}"`
do
    CMAKE_VAR_NAME=`cut -d'=' -f1 <<<"$cmd" | sed s/^@*//g | sed s/@*\$//g`
    CMAKE_VAR_VALUE=`cut -d'=' -f2- <<<"$cmd" | sed s/^@*//g | tr "@" " "`
    CMAKE_CONFIG_OPTS="$CMAKE_CONFIG_OPTS -D${CMAKE_VAR_NAME}='${CMAKE_VAR_VALUE}'"
done

PROJECT_SETUP_ENV_COMMANDS=""
# Excution des commandes dans ${ENV_SETUP_COMMANDS}
for cmd in `tr " \t\n" "@@ "<<<"${ENV_SETUP_COMMANDS}"`
do
    END_CMD=`sed s/^@*//g <<<"$cmd" | tr "@" " "`
    [ "${PROJECT_SETUP_ENV_COMMANDS}" != "" ] && PROJECT_SETUP_ENV_COMMANDS="$PROJECT_SETUP_ENV_COMMANDS ;"
    PROJECT_SETUP_ENV_COMMANDS="${PROJECT_SETUP_ENV_COMMANDS}${END_CMD}"
    if [ $PRINT_ENV_COMMANDS ]
    then
        echo "+ $END_CMD"
    fi
    $END_CMD
done

# add environment setup commands in a CMake variable so cmake can generate a launch script
if [ "${PROJECT_SETUP_ENV_COMMANDS}" != "" ]
then
    CMAKE_CONFIG_OPTS="$CMAKE_CONFIG_OPTS -DPROJECT_SETUP_ENV_COMMANDS='${PROJECT_SETUP_ENV_COMMANDS}'"
fi
PROJECT_SETUP_ENV_SCRIPT="${BUILD_DIR}/setup-env.sh"

# project build type
[ $PROJECT_BUILD_TYPE ] || PROJECT_BUILD_TYPE=Release

# default is no verbose
[ $CMAKE_VERBOSE_MAKEFILE ] || CMAKE_VERBOSE_MAKEFILE=OFF

# jump to build directory
cd ${BUILD_DIR}

# configuration avec CMake (ou ccmake is pas de pre-configuration fournie)
echo "Configuring project with CMake ..."
[ "${CMAKE_CONFIG_OPTS}" == "" ] && FORCE_INTERACTIVE=1
CONFIG_TOOL=${CMAKE}
[ $FORCE_INTERACTIVE ] && CONFIG_TOOL=${CCMAKE}
if [ ! ${EXISTING_BUILD_DIR} ]
then
  set -x
  eval ${CONFIG_TOOL} \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
  ${CMAKE_CONFIG_OPTS} \
  -DCMAKE_BUILD_TYPE=${PROJECT_BUILD_TYPE} \
  -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE} \
  -DHOST_ALWAYS_USE_MPIRUN=${HOST_ALWAYS_USE_MPIRUN} \
  -DHOST_HW_THREADS=${NTHREADS} \
  -DHOST_HW_CORES=${NCORES} \
  -DHOST_HW_PARTITION=${HOST_PARTITION} \
  -DXNB_PRODUCT_VARIANT=${PRODUCT_VARIANT} \
  ${SOURCE_DIR}
  set +x
else
  set -x
  ${CONFIG_TOOL} .
  set +x
fi

# print generated help if available
[ -f README.txt ] && cat README.txt

