# OpenMP

Tutorial the OpenMP version

## USE

The openmp version corresponds, for the most part, to a rewriting of include/parallel/thread/tbb.hxx . 
The reduction functions have been rewritten. 

WARNING 

Mutex read and written are not checked.

advice:
   - KMP_AFFINITY=scatter,1 (to bind one thread per core)
   - OPT += -D__use_lib_omp -qopenmp (KNL->) -qoffload-arch=mic-avx512 

## MUTEX

One of the current defects of openmp is the additional cost of mutex. Indeed on small cases, using the mutex provided by the openMP library, we can be almost 50% behind the TBB version.

4 types are therefore used: 

  - mutex openmp
  - mutex from the standard library
  - mutex from the TBB library
  - mutex handwritten.

## NOTE

 - On small cases (256 000/4M) of atoms similar (or less efficient) results are obtained on sandy and broadwell. On a KNL the performances are currently less interesting (20%).
 - The values of the grains are not necessarily adapted and when the code uses a parallel_for without grain, if an element (vector type) is defined before the loop, it will be defined n times. In tbb this is defined once per block.  




