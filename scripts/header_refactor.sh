#!/usr/bin/bash

REALLYDOIT=0

for HDR in `find src -name "*.h"`
do
  FNAME=`basename ${HDR}`

  if [ `echo "$HDR"|grep -c "/deprecated/"` -gt 0 ]
  then
	  echo "SKIP $HDR"
  else

  # find include PATH
  INCPATH=`echo $HDR | sed 's/^.*include\///g'`

  # does it use a specific namespace ?
  NAMESPACE=`cat $HDR | grep "namespace" | head -n 1| sed 's/^namespace[ \t]*\([a-zA-Z0-9_]*\).*/\1/g'`
  HASNAMESPACE=0
  [ "$NAMESPACE" != "" ] && HASNAMESPACE=1

  # what software package is it in ?
  PACKAGE=`echo $HDR | cut -d'/' -f2`

  MESG=""
  MOVE=0
  COMPTREE=0

  if [ "$INCPATH" == "$HDR" ]
  then

#    echo "SKIP internal header $FNAME in" `dirname $HDR`
    MESG="$INCPATH (INTERNAL) : NAMESPACE=$NAMESPACE : PKG=$PACKAGE : COMP=$COMPONENT"

  else

    # what software component is it in ?
    if [ "$PACKAGE" == "onika" ]
    then
      COMPONENT=`echo ${INCPATH} | cut -d'/' -f2`
      HASCOMP=1
      COMPTREE=0
      [ "$COMPONENT" == "$FNAME" ] && HASCOMP=0 && COMPONENT=""
    else
      COMPDIR=`echo $HDR|sed 's/\/include\/.*//g'`
      COMPONENT=`echo $COMPDIR | cut -d'/' -f3-`
      HASCOMP=1
      COMPTREE=1
      [ -d src/$PACKAGE/$COMPONENT ] || HASCOMP=0
      [ "$COMPONENT" == "include" ] && HASCOMP=0
      [ $HASCOMP -eq 1 ] || COMPONENT=""
      [ "$COMPONENT" == "" ] && HASCOMP=0
    fi

    MESG="$INCPATH : NAMESPACE=$NAMESPACE : PKG=$PACKAGE : COMP=${COMPONENT} (HC=${HASCOMP})"

    # does include respect include path convention ?
    if [ $HASCOMP -eq 1 ]
    then
      INCDST=${PACKAGE}/${COMPONENT}/${FNAME}
      if [ "${INCPATH}" != "${PACKAGE}/${COMPONENT}/${FNAME}" ]
      then
	      [ "$MESG" != "" ] && echo "$MESG" && MESG=""
	      MOVE=1
	      [ $COMPTREE -eq 1 ] && MOVEDST=src/${PACKAGE}/${COMPONENT}/include/${PACKAGE}/${COMPONENT}/${FNAME}
	      [ $COMPTREE -eq 0 ] && MOVEDST=src/${PACKAGE}/include/${PACKAGE}/${COMPONENT}/${FNAME}
      fi
    else
      INCDST=${PACKAGE}/${FNAME}
      if [ "${INCPATH}" != "${PACKAGE}/${FNAME}" ]
      then
	      [ "$MESG" != "" ] && echo "$MESG" && MESG=""
	      MOVE=1
	      MOVEDST=src/${PACKAGE}/include/${PACKAGE}/${FNAME}
      fi
    fi

    # 2. Check for remaining bad include paths
    NBOCC=`find src -name "${FNAME}" -print|wc -l`
    if [ $NBOCC -eq 1 ]
    then
      NBINCFORMS=`grep -r "^#include.*[\"</]${FNAME}[\">]" src|grep -v ".cmake:"|cut -d':' -f2|grep -v "^#include <${INCDST//\//\\\/}>"|sed "s/#include[ \t]*[\"<]\(.*${FNAME}\)[\">].*/\1/g" |sort -u|wc -l`
      if [ $NBINCFORMS -gt 0 ]
      then
	[ "$MESG" != "" ] && echo "$MESG" && MESG=""
        echo "file $FNAME is unique, but included with $NBINCFORMS different bad forms :"
	grep -r "^#include.*[\"</]${FNAME}[\">]" src|grep -v ".cmake:"|cut -d':' -f2|sed "s/#include[ \t]*[\"<]\(.*${FNAME}\)[\">].*/\1/g" |sort -u
        echo "need to patch following files"
	grep -r "^#include.*[\"</]${FNAME}[\">]" src|grep -v ".cmake:"|cut -d':' -f1
	echo -n "fix include path (y/n) ? [n] " ; read DOIT
	if [ "$DOIT" == "y" ]
	then
          grep -r "^#include.*[\"</]${FNAME}[\">]" src|grep -v ".cmake:"|cut -d':' -f1 | xargs sed -i "s/^#include[ \t]*[\"<]\(.*\/\)*\(${FNAME}\)[\">]/#include <${INCDST//\//\\\/}>/g"
	fi
      fi
    fi
  fi

  if [ $MOVE -eq 1 ]
  then
	  if [ -f $MOVEDST ]
	  then
		  echo "Can't move $HDR to $MOVEDST, file already exists"
		  exit 1
	  fi
	  NBPLACES=`grep -h -r "^#include.*[/\"\<]${FNAME}[\"\>]" src apps | sed 's/\([">]\)[ \t]*\/\/.*/\1/g' |  tr "<>" "\"\"" | sort -u|wc -l`
	  if [ $NBPLACES -eq 0 ]
	  then
		  echo "Error: $FNAME not included anywhere"
		  exit 1
	  fi
	  if [ $NBPLACES -gt 1 ]
	  then
		  echo "Error: $FNAME included with $NBPLACES different paths"
		  grep -h -r "^#include.*[/\"\<]${FNAME}[\"\>]" src apps | sed 's/\([">]\)[ \t]*\/\/.*/\1/g' | sort -u
		  exit 1
	  fi
          
	  echo "    MOVE FILE : $HDR => $MOVEDST"
	  if [ $REALLYDOIT -eq 1 ]
	  then
	    set -x
	    mkdir -p `dirname ${MOVEDST}`
	    git mv $HDR $MOVEDST
	    set +x
	  fi

	  echo "    MOVE INC  : $INCPATH => $INCDST :"
	  for TOPATCH in `grep -r -l "^#include.*[\"<]${INCPATH}[\">]" src`
	  do
	    echo "        patch $TOPATCH"
	    if [ $REALLYDOIT -eq 1 ]
	    then
	      set -x
	      sed -i "s/^#include[ ]*.${INCPATH//\//\\\/}./#include \<${INCDST//\//\\\/}\>/g" $TOPATCH
	      set +x
	    fi
	  done
  fi

  if [ ${HASNAMESPACE} -eq 1 ] && [ "${NAMESPACE}" != "${PACKAGE}" ]
  then
    [ "$MESG" != "" ] && echo "$MESG" && MESG=""
    echo "    namespace : ${NAMESPACE} => ${PACKAGE}"
  fi

  fi # if deprecated

done

