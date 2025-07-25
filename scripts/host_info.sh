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

FUNCTION_SCRIPT="$(dirname $0)/functions.sh"
if [ ! -f ${FUNCTION_SCRIPT} ]
then
	echo "Cannot find script file '${FUNCTION_SCRIPT}'"
	exit 1
fi

source ${FUNCTION_SCRIPT}

# host machine type
#set -x
CCC_OS_CMD=`which ccc_os 2>/dev/null`
CEA_OS_CMD=`which cea_os 2>/dev/null`
LSB_RELEASE_CMD=`which lsb_release 2>/dev/null`
#__nvcc_device_query

[ $CCC_OS_CMD ] && HOST_OS=`$CCC_OS_CMD`
[ $CEA_OS_CMD ] && [ -z "$HOST_OS" ] && HOST_OS=`$CEA_OS_CMD`
[ $LSB_RELEASE_CMD ] && [ -z "$HOST_OS" ] && HOST_OS=`$LSB_RELEASE_CMD -d|sed 's/Description:[ \t]*//g'|cut -d'.' -f1,2|tr " " "-"|cut -d'-' -f1,2`
[ $HOST_OS ] || HOST_OS=`uname -m`

HOST_PARTITION="$SLURM_JOB_PARTITION"
[ $HOST_PARTITION ] || HOST_PARTITION="none"

if [ `which hwloc-ls` ]
then
  HWLOC_INFO=`hwloc-ls --of console --whole-system -p`
  NCORES=`echo "$HWLOC_INFO"|grep "Core P#"|wc -l`
  NTHREADS=`echo "$HWLOC_INFO"|grep "PU P#"|wc -l`
else
  NCORES=4
  NTHREADS=4
fi

echo "${HOST_OS} ${HOST_PARTITION} ${NCORES} ${NTHREADS}"

