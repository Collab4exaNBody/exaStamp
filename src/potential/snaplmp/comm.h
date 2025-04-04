#ifndef LMP_COMM_H
#define LMP_COMM_H

namespace LAMMPS_NS
{

  struct CommunicatorInfo
  {
    int me = 0;
    int nprocs = 1;
  };

}

#endif

