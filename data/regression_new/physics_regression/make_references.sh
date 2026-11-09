#!/bin/bash
# (Re)generates the committed references of physics regression cases.
# usage: make_references.sh <exaStamp launcher> <case ...|all>
#   e.g. make_references.sh ~/local/exaStamp/bin/exaStamp pair_lj_single
# Each case runs once (1 MPI process, 1 thread) in a temporary copy; on success ref_values.dat and ref_thermo.yaml
# are copied back into the case folder. Only regenerate references for an intentional change of the physics.
here=$(cd "$(dirname "$0")" && pwd)
exe=$1 ; shift
[ -n "$exe" ] && [ $# -gt 0 ] || { echo "usage: $0 <exaStamp launcher> <case ...|all>" >&2 ; exit 1 ; }
cases="$*"
[ "$cases" = "all" ] && cases=$(cd "$here" && ls -d */ | tr -d '/')
status=0
for c in $cases ; do
  [ -f "$here/$c/run.msp" ] || { echo "$c: no run.msp" >&2 ; status=1 ; continue ; }
  tmp=$(mktemp -d)
  cp -rL "$here/$c/." "$tmp/" && rm -f "$tmp/ref_values.dat" "$tmp/ref_thermo.yaml"
  ln -s "$here/../../potentials" "$tmp/data"
  if ( cd "$tmp" && OMP_NUM_THREADS=1 mpirun -np 1 "$exe" "$here/$c/run.msp" "$here/create_references.msp" > log.txt 2>&1 ) ; then
    cp "$tmp/ref_values.dat" "$tmp/ref_thermo.yaml" "$here/$c/" && rm -rf "$tmp" && echo "$c: references written"
  else
    echo "$c: FAILED, see $tmp/log.txt" ; status=1
  fi
done
exit $status
