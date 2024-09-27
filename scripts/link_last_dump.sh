PROTNUM=`ls -t $1*.MpiIO | head -1 | sed "s/\(.*\).MpiIO/\1/g"`
for PROT in `ls ${PROTNUM}*`
do
	PART=`echo ${PROT}|sed "s/.*\.MpiIO//g"`
	echo "${PROT} -> lastLegacyDump${PART}"
	ln -sf ${PROT} lastLegacyDump${PART}
done

