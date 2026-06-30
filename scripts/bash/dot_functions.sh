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
mdot()
{
  OFMT=$1
  shift
  for filename in $*
  do
    EXT="${filename##*.}"
    TOOL=`echo ${EXT}|cut -d "-" -f1`
    # echo "TOOL $TOOL"
    ARGS=`echo "${EXT##${TOOL}}"|sed "s/-/ -/g"`
    # echo "ARGS $ARGS"
    BASE=`basename ${filename} .${EXT}`
    OUTPUTFILEOPT="-o ${BASE}.${OFMT}"
    [ ${OFMT} == "xlib" ] && OUTPUTFILEOPT=""
    if [ ${TOOL} == "dot" ]
    then
      CMD="${TOOL} ${ARGS} ${filename} | gvpack -u -- | neato -T${OFMT} -n2 ${OUTPUTFILEOPT}"
    else
      CMD="${TOOL} -T${OFMT} ${ARGS} ${filename} ${OUTPUTFILEOPT}"
    fi
    echo "${CMD}"
    eval ${CMD}
  done
}

