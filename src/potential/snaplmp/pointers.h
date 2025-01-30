#ifndef LMP_POINTERS_H
#define LMP_POINTERS_H

#include <memory>
#include <cstdint>
#include <onika/cuda/cuda_math.h>

#include "error.h"
#include "comm.h"

#define FLERR __FILE__,__LINE__
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

namespace LAMMPS_NS
{
  using bigint = int64_t;

  class Memory;

  struct LAMMPS
  {
    ErrorLogWrapper* error = nullptr;    
    Memory* memory = nullptr;
    CommunicatorInfo* comm = nullptr;
  };

	class Pointers
	{
	public:
	  inline Pointers(LAMMPS* ptr)
	    : lmp(ptr)
	    , error(ptr->error)
	    , memory(ptr->memory)
	    , comm(ptr->comm)
	  {}
	    
    virtual ~Pointers() = default;
    Pointers() = delete;
    Pointers(const Pointers &) = default;
    Pointers(Pointers &&) = delete;
    Pointers & operator=(const Pointers&) = delete;
    Pointers & operator=(Pointers&&) = delete;

   protected:
    LAMMPS *lmp = nullptr;
    ErrorLogWrapper* error = nullptr;
    Memory* memory = nullptr;
    CommunicatorInfo* comm = nullptr;
	};
	
}

#endif

