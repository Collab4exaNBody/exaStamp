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

#[ "$OPNAME" ] || OPNAME="meam_force"
[ "$OPNAME" ] || OPNAME="meam_lj_force;chunk_neighbors;particle_displ_over;check_and_update_particles;loop"
[ "$PREFIX" ] || PREFIX="result_"

#CUSTOMCOUNTERS="chk_ratio;amr_res;amr_dens;nb_iter"
CUSTOMCOUNTERS="nb_iter"
chk_ratio="grep '<d / chunk ratio       = ' \${LOGFILE}|cut -d'=' -f2|cut -d',' -f1|tr '.' ','|tr -d ' '"
amr_res="grep 'AMR res        =' \${LOGFILE}|cut -d'=' -f2|tr '.' ','|tr -d ' '"
amr_dens="grep 'AMR density    =' \${LOGFILE}|cut -d'=' -f2|tr '.' ','|tr -d ' '"
chk_mem="grep 'chunk_neighbors.*GridChunkNeighbors' \${LOGFILE}|sed s/chunk_neighbors[^0-9]*//g |cut -d' ' -f1"
nb_iter="grep '    scheme ........................................' \${LOGFILE}|tr -s ' '|cut -d' ' -f7"
#dispover="grep '        particle_displ_over .......................' \${LOGFILE}|tr -s ' '|cut -d' ' -f4"

FIRSTLOG=`ls ${PREFIX}*|head -1`
#echo "FIRSTLOG=${FIRSTLOG}"
LISTVAR=`echo "${FIRSTLOG}" | sed -e "s/^${PREFIX}//g" -e "s/.txt\$//g"  | tr "_" " " | sed "s/=[0-9]*//g"`
#echo "LISTVAR=${LISTVAR}"

# print CSV header row
for V in $LISTVAR; do echo -n "$V;" ; done
echo -n "${OPNAME}"
[ "${CUSTOMCOUNTERS}" ] && echo -n ";${CUSTOMCOUNTERS}"
echo ""

# populate CSV rows parsing logs
for LOGFILE in ${PREFIX}*
do
#  LOGFILE=`ls -t $R/log*.txt|head -1`
#  echo "REP=$R , LOG=$LOGFILE"
  for V in $LISTVAR; do eval $V="N/A" ; done
  eval `echo $LOGFILE|sed -e "s/^${PREFIX}//g" -e "s/.txt\$//g" |sed "s/_/;/g"`
  for V in $LISTVAR; do eval "echo -n \$$V\;" ; done
  SEP=""
  for OP in `echo $OPNAME|tr ';' ' '`
  do
    echo -n "$SEP\""
    COUNT=`grep -c "${OP} \." ${LOGFILE}`
    if [ $COUNT -ge 1 ]
    then
      grep "${OP} \." ${LOGFILE} |tail -n1|tr -s " "|cut -d' ' -f4|tr "." ","| tr -d " \n"
    else
      echo -n "N/A"
    fi
    echo -n '"'
    SEP=";"
  done
  for CCNT in `echo $CUSTOMCOUNTERS|tr ';' ' '`
  do
    echo -n "$SEP\""
    eval "${!CCNT}" | tr -d " \n"
    echo -n '"'
    SEP=";"
  done
  echo ""
done
