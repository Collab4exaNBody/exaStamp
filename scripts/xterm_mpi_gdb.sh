#!/bin/bash

[ $GDB ] || GDB=gdb

[ $POFFSETX ] || POFFSETX=0
[ $POFFSETY ] || POFFSETY=32
[ $PSIZEX ] || PSIZEX=302
[ $PSIZEY ] || PSIZEY=100
[ $TSIZEX ] || TSIZEX=80
[ $TSIZEY ] || TSIZEY=16
[ $NX ] || NX=11
[ $XTERMOPT ] || XTERMOPT="-fn -misc-fixed-medium-r-normal--8-80-75-75-c-50-iso10646-1 -b 0 +sb"

[ $MYRANK ] || [ $PMI_RANK ] && MYRANK=$PMI_RANK
[ $MYRANK ] || [ $PMIX_RANK ] && MYRANK=$PMIX_RANK
[ $MYRANK ] || [ $SLURM_PROCID ] && MYRANK=$SLURM_PROCID

if [ ! $MYRANK ]
then
    echo "Error: could not determine MPI rank"
    env | grep -i mpi
    exit 1
fi

let PX=${MYRANK}%${NX}
let PY=${MYRANK}/${NX}
let GEOMPOSX=${POFFSETX}+${PX}*${PSIZEX}
let GEOMPOSY=${POFFSETY}+${PY}*${PSIZEY}

TITLE="MASTER"
if [ $MYRANK -gt 0 ]
then
  TITLE="P$MYRANK"
fi

set -x
xterm -geom ${TSIZEX}x${TSIZEY}+${GEOMPOSX}+${GEOMPOSY} -title "${TITLE}" ${XTERMOPT} -e ${GDB} -ex 'set pagination off' -ex run --args $*

