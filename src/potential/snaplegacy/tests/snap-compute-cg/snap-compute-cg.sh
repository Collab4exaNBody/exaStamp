#!/bin/bash -eu

# param covers 2 examples:
# 1. naive example (all parameters are ones)
# 2. /mnt/cesimat/exastamp/SNAP/CodeSaclaySNAP example
#    "3.5 2" (CodeSaclaySNAP example) <=> param "7. 1" (same thing but with ntype = 1)

for param in "1.0 1" "3.5 2"
do
  echo ""
  echo "param: $param"
  IFS=" " read -r jmax ntype <<< "$param"
  echo "jmax: $jmax"
  echo "ntype: $ntype"

  # Run fortran reference code => produce cg.ref

  echo ""
  ./snap-compute-cg-f90 "$jmax" "$ntype"

  # Run cpp code (+ diff with cg.ref)

  echo ""
  ./snap-compute-cg-cpp "$jmax" "$ntype"

  echo ""
  echo "============================================="



done
