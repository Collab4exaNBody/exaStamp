#!/bin/sh

DIRNAME=$1
FILES=`ls $DIRNAME|sort`

for FC in $FILES
do
  if [ "$FC" != "data_dir.md5" ]
  then
    F=`basename $FC .xz`
    if [ "$F" == "$FC" ]
    then
      echo -n "compress $F ..."
      xz -T 0 $DIRNAME/$F
      echo " done"
    fi
  fi
done

