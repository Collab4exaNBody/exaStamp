#!/bin/bash

[ $GDBSERVER ] || GDBSERVER=cuda-gdbserver

[ $MYRANK ] || [ $PMI_RANK ] && MYRANK=$PMI_RANK
[ $MYRANK ] || [ $PMIX_RANK ] && MYRANK=$PMIX_RANK
[ $MYRANK ] || [ $SLURM_PROCID ] && MYRANK=$SLURM_PROCID

if [ ! $MYRANK ]
then
    echo "Error: could not determine MPI rank"
    env | grep -i mpi
    exit 1
fi

let HOSTPORT=28000+${MYRANK}
CLIENTHOST=$1
shift
SERVERHOST=`hostname`

touch dbgtargets/${MYRANK}.${SERVERHOST}.${HOSTPORT}

set -x
${GDBSERVER} ${CLIENTHOST}:${HOSTPORT} $*

