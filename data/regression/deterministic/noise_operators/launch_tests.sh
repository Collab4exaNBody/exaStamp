#!/bin/bash

source ~/local/exaStamp/bin/exaStamp
execcmd=/home/lafourcadep/local/onika/bin/onika-exec

for name in uniform_r uniform_v uniform_f gaussian_v gaussian_v gaussian_f;
do
    for OMP in 1 4;
    do
        for MPI in 1 2;
        do
            echo "Testing with ${OMP} threads and ${MPI} MPI processes"
            OMP_NUM_THREADS=${OMP} mpirun -np ${MPI} ${execcmd} ${name}.msp
        done
    done
done


