#!/bin/bash -eu

# param covers 2 examples:
# 1. naive example (all parameters are ones)
# 2. /mnt/cesimat/exastamp/SNAP/CodeSaclaySNAP example

for param in "1. 2. 3. 5.0 1" "-2.6630 -2.6630 -2.6630 5.0000 7"
do
  echo ""
  echo "param: $param"
  IFS=" " read -r x y z rcut jmax <<< "$param"
  echo "x: $x"
  echo "y: $y"
  echo "z: $z"
  echo "rcut: $rcut"
  echo "jmax: $jmax"

  # Run fortran reference code => produce gsh.ref

  echo ""
  ./snap-compute-gsh-f90 "$x" "$y" "$z" "$rcut" "$jmax"

  # Run cpp code (+ diff with gsh.ref)

  echo ""
  ./snap-compute-gsh-cpp "$x" "$y" "$z" "$rcut" "$jmax"

  echo ""
  echo "============================================="
done
