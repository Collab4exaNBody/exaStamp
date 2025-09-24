#!/bin/bash

execcmd=/home/lafourcadep/local/onika/bin/onika-exec

for name in random_gaussian_r random_gaussian_v random_uniform_r random_uniform_v random_langevin;
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


