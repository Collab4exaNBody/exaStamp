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

