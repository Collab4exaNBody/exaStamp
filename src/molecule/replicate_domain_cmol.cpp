#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/geometry.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/log.h>
#include <exanb/core/thread.h>
#include <exanb/grid_cell_particles/replicate_domain.h>

#include <vector>
#include <iostream>
#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  struct ShiftParticleIdAndCMolFunctor
  {
    template<class FieldTupleT>
    inline void operator () ( FieldTupleT& tp , int64_t id_offset ) const
    {
      if constexpr ( tp.has_field( field::id ) )
      {
        tp[field::id] += id_offset;
      }
      if constexpr ( tp.has_field( field::cmol ) )
      {
        auto cmol = tp[field::cmol];
        for(size_t i=0;i<cmol.size();i++)
        {
          if( cmol[i] != std::numeric_limits<uint64_t>::max() ) cmol[i] += id_offset;
        }
        tp[field::cmol] = cmol;
      }
    }
  };

  template<class GridT> using ReplicateDomainCMol = ReplicateDomain<GridT,ShiftParticleIdAndCMolFunctor>;

   // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "replicate_domain_cmol", make_grid_variant_operator< ReplicateDomainCMol > );
  }

}

