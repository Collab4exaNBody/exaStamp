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

