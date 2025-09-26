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

# usefull functions
function version { echo "$@" | gawk -F. '{ printf("%03d%03d%03d\n", $1,$2,$3); }'; }

function mpi_command
{
   NUMPROCS=$1
   NUMTHREADS=$2
   shift ; shift
   CMD=`eval echo $MPI_LAUNCHER_CMD`
   echo $CMD
}

function mpi_preallocate_command
{
   NUMPROCS=$1
   NUMTHREADS=$2
   shift ; shift
   CMD=`eval echo $MPI_PREALLOCATE_CMD`
   echo $CMD
}

function run_mpi_job
{
  CMD=`mpi_command $1 $2`
  shift ; shift
  if [ "$1" == "DEBUG" ]
  then
    shift
    echo $CMD $*
  else
    $CMD $*
  fi
}


