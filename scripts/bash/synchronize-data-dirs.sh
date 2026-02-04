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

DBG=""
if [ "$1" == "--fake" ]
then
  DBG="echo"
  shift
fi

SRCDIR=$1
DSTDIR=$2

SCRIPT_DIR=`dirname $0`
echo ""
echo "================== $SRCDIR ================="
$SCRIPT_DIR/hash-data-dir.sh $SRCDIR
echo ""
echo "================== $DSTDIR ================"
$SCRIPT_DIR/hash-data-dir.sh $DSTDIR
echo ""

echo "================== SYNCHRONIZE ================"
for SF in `cat $SRCDIR/data_dir.md5 $DSTDIR/data_dir.md5|cut -d' ' -f1|sort -u`
do
  SH=0
  SD=0
  SB=`cat $SRCDIR/data_dir.md5|grep "^$SF "|cut -d' ' -f2,3`
  if [ "$SB" != "" ]
  then
   SH=`echo $SB|cut -d' ' -f1`
   SD=`echo $SB|cut -d' ' -f2`
  fi
  DB=`cat $DSTDIR/data_dir.md5|grep "^$SF "|cut -d' ' -f2,3`
  DD=0
  DH=0
  if [ "$DB" != "" ]
  then
   DH=`echo $DB|cut -d' ' -f1`
   DD=`echo $DB|cut -d' ' -f2`
  fi
  if [ "x$SH" != "x$DH" ]
  then
    if [ $SD -ge $DD ]
    then
	printf "%-30s SRC -> DST\n" "$SF"
	rm -f ${DSTDIR}/${SF} ${DSTDIR}/${SF}.xz
	[ -f ${SRCDIR}/${SF} ] && cp ${SRCDIR}/${SF} ${DSTDIR}/
	[ -f ${SRCDIR}/${SF}.xz ] && cp ${SRCDIR}/${SF}.xz ${DSTDIR}/
    else
	printf "%-30s DST -> SRC\n" "$SF"
	rm -f ${SRCDIR}/${SF} ${SRCDIR}/${SF}.xz
	[ -f ${DSTDIR}/${SF} ] && cp ${DSTDIR}/${SF} ${SRCDIR}/
	[ -f ${DSTDIR}/${SF}.xz ] && cp ${DSTDIR}/${SF}.xz ${SRCDIR}/
    fi
  else
    printf "%-30s synchronized\n" "$SF"
  fi
done

