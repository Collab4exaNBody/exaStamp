#!/bin/bash
# det_shift then noise: LAMMPS, lmp2xyz, exaStamp (CPU, 1 thread), with timing
H=$(cd "$(dirname "$0")" && pwd); source ~/myv2env/bin/activate
for v in det_shift noise; do
  cd $H/$v
  t0=$(date +%s); OMP_NUM_THREADS=1 nice /local_home/lafourcadep/CODES/ATOMISTIC/LAMMPS/lammps/build/lmp -in in.lammps > lmp.out 2>&1; echo "$v lmp rc=$? $(( $(date +%s)-t0 )) s, $(grep 'Neighbor list builds' log.lammps)"
  python3 ../../../snap/lmp2xyz.py init.dump init.xyz Au
  t0=$(date +%s); OMP_NUM_THREADS=1 timeout 1800 nice mpirun -np 1 bash ~/local/exaStamp/bin/exaStamp run.msp --nogpu true > exa.log 2>&1; rc=$?
  echo "$v exa rc=$rc wall $(( $(date +%s)-t0 )) s; last csv write $(( $(stat -c %Y thermodynamic_state.csv)-t0 )) s after start"
  python3 ../../compare.py . || true
done
