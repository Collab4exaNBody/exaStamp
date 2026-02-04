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

DIRNAME=$1

FILES=`ls $DIRNAME|sort`
MD5DB=""

for FC in $FILES
do
    F=`basename $FC .xz`
    [ "$F" != "$FC" ] && [ -f "$DIRNAME/$F" ] && echo "remove $FC ($F exists)" && rm "$DIRNAME/$FC"
done

for FC in $FILES
do
  if [ "$FC" != "data_dir.md5" ]
  then
    F=`basename $FC .xz`
    FDATE=`stat -c "%Y" $DIRNAME/$FC`
    COMPUTE_HASH="YES"
    if [ -f $DIRNAME/data_dir.md5 ]
    then
	DB=`cat $DIRNAME/data_dir.md5 | grep "^$F " | cut -d' ' -f2,3`
        DBDATE=0
        DBHASH=0

        if [ "$DB" != "" ]
	then
	  DBDATE=`echo $DB |cut -d' ' -f2`
        fi

	if [ $DBDATE -ge $FDATE ]
        then
	  printf "%-30s %-20s" "$F" "(cached)"
	  HASH=`echo $DB |cut -d' ' -f1`
          COMPUTE_HASH=""
	fi
    fi

    if [ $COMPUTE_HASH ]
    then
      if [ "$F" != "$FC" ]
      then
        printf "%-30s %-20s" "$F" "(compressed)"
        HASH=`xz -T 0 -d -c $DIRNAME/$FC | md5sum | cut -d' ' -f1`
      else
        printf "%-30s %-20s" "$F" ""
        HASH=`md5sum $DIRNAME/$FC | cut -d' ' -f1`
      fi
    fi

    printf " %-33s%-15s\n" "${HASH}" "${FDATE}"
    MD5DB="${MD5DB}${F}:${HASH}:${FDATE};"
  fi
done

echo -n ${MD5DB} | tr ":;" " \n" | sort > $DIRNAME/data_dir.md5

