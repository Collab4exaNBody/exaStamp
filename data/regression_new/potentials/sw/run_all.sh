#!/bin/bash
# Run every exaStamp Stillinger-Weber case and diff it against its LAMMPS log.lammps.
# usage: run_all.sh [case ...]   env: NP (default 1), EXASTAMP, extra exaStamp args in EXA_ARGS (e.g. "--nogpu true")
here=$(cd "$(dirname "$0")" && pwd)
EXASTAMP=${EXASTAMP:-$(command -v exaStamp)}
cases=${*:-Si GaN}
for c in $cases; do
  cd "$here/$c" || exit 1
  mpirun -np ${NP:-1} bash "$EXASTAMP" ${c}_exaStamp.msp $EXA_ARGS > exa.log 2>&1 || { echo "$c: exaStamp FAILED (see $c/exa.log)"; continue; }
  python3 "$here/../../mliaps/snap/compare.py" "$here/$c"
done
