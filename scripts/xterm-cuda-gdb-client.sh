#!/bin/bash


[ $GDB ] || GDB=cuda-gdb

[ $POFFSETX ] || POFFSETX=0
[ $POFFSETY ] || POFFSETY=32
[ $PSIZEX ] || PSIZEX=302
[ $PSIZEY ] || PSIZEY=100
[ $TSIZEX ] || TSIZEX=80
[ $TSIZEY ] || TSIZEY=25
[ $NX ] || NX=11
[ $XTERMOPT ] || XTERMOPT="-fn -misc-fixed-medium-r-normal--8-80-75-75-c-50-iso10646-1 -b 0 +sb"

for TARGET in `ls dbgtargets`
do
	MYRANK=`echo $TARGET| cut -d'.' -f1`
	HOST=`echo $TARGET| cut -d'.' -f2`
	PORT=`echo $TARGET| cut -d'.' -f3`

	let PX=${MYRANK}%${NX}
	let PY=${MYRANK}/${NX}
	let GEOMPOSX=${POFFSETX}+${PX}*${PSIZEX}
	let GEOMPOSY=${POFFSETY}+${PY}*${PSIZEY}

	TITLE="MASTER"
	if [ $MYRANK -gt 0 ]
	then
	  TITLE="P$MYRANK"
	fi

	echo "RANK ${MYRANK} -> ${HOST}:${PORT}"
	xterm -geom ${TSIZEX}x${TSIZEY}+${GEOMPOSX}+${GEOMPOSY} -title "${TITLE}" ${XTERMOPT} -e ${GDB} -ex 'set pagination off' -ex "target remote ${HOST}:${PORT}" -ex continue &
done

echo "Waiting for debug session to end ..."
wait

