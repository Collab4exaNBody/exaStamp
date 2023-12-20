#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <exanb/core/basic_types_stream.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <vector>

namespace exaStamp
{

  using namespace exanb;

  template<typename GridT, class = AssertGridHasFields< GridT, field::_type, field::_charge> >
  struct CopyChargeFromSpecyToParticleNode : public OperatorNode
  {
    ADD_SLOT( ParticleSpecies , species , INPUT , REQUIRED );
    ADD_SLOT( GridT           , grid    , INPUT_OUTPUT );

    inline void execute () override final
    {
      const ParticleSpecie* __restrict__ sp = species->data();

      auto cells = grid->cells();
      IJK dims = grid->dimension();
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,i,loc)
        {
          const uint8_t * __restrict__ ptype = cells[i][field::type];
          double* __restrict__ charge = cells[i][field::charge];
          size_t n = cells[i].size();
#         pragma omp simd
          for(size_t j=0;j<n;j++)
          {
            assert( ptype[j]>=0 && ptype[j]<species->size() );
            charge[j] = sp[ ptype[j] ].m_charge;
          }
        }
        GRID_OMP_FOR_END
      }
    }

  };
  
  namespace tmplhelper { template<class GridT> using CopyChargeFromSpecyToParticleNode = ::exaStamp::CopyChargeFromSpecyToParticleNode<GridT>; }

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "copy_charge_specy_to_particle", make_grid_variant_operator< tmplhelper::CopyChargeFromSpecyToParticleNode > );
  }

}

