#!/bin/sh

for OBJ in `find . -name "*.o"`
do
    ( nm $OBJ > /dev/null ) || ( echo "delete $OBJ" && rm $OBJ )
done

