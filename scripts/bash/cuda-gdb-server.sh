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

