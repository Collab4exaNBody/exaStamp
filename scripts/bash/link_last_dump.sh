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
PROTNUM=`ls -t *.MpiIO | head -1 | sed "s/\(.*\).MpiIO/\1/g"`
for PROT in `ls ${PROTNUM}*`
do
	PART=`echo ${PROT}|sed "s/.*\.MpiIO//g"`
	echo "${PROT} -> lastLegacyDump${PART}"
	ln -sf ${PROT} lastLegacyDump${PART}
done

