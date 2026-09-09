/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

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

