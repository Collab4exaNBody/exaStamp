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

