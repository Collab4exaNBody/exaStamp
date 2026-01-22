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

[ $GDB ] || GDB=gdb

[ $POFFSETX ] || POFFSETX=0
[ $POFFSETY ] || POFFSETY=32
[ $PSIZEX ] || PSIZEX=302
[ $PSIZEY ] || PSIZEY=100
[ $TSIZEX ] || TSIZEX=80
[ $TSIZEY ] || TSIZEY=16
[ $NX ] || NX=4
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
xterm -geom ${TSIZEX}x${TSIZEY}+${GEOMPOSX}+${GEOMPOSY} -title "${TITLE}" ${XTERMOPT} -e ${GDB} -ex 'set pagination off' -ex "set confirm off" -ex run --args $*

