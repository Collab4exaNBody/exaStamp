#!/bin/sh

for f in `find . -not -user $USER`
do
   [ -d $f ] && echo "dir $f : && mv $f ${f}.old && mkdir $f && mv ${f}.old/* $f/ && rm -r ${f}.old"
   [ -f $f ] && echo "file $f"
done

