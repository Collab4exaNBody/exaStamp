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


