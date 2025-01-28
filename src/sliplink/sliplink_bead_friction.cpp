#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <onika/memory/allocator.h>
#include <onika/parallel/random.h>
#include <exaStamp/sliplink/sliplink.h>

#include <memory>


namespace exaStamp
{

  template< class GridT >
  struct SliplinkForceOverdamped : public OperatorNode
  {      
    ADD_SLOT( GridT , grid ,INPUT_OUTPUT);
    ADD_SLOT(SlipLinkParameters , sliplink_config , INPUT, REQUIRED );

    inline void execute () override final
    {
      GridT& grid = *(this->grid);
      
      auto cells = grid.cells();
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();

      const double sigma1 = sliplink_config->sigma1;

#     pragma omp parallel
      {
        auto& re = rand::random_engine();
        std::normal_distribution<double> gaussian(0.0,sigma1);

        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc)
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          size_t n = cells[i].size();
          
          auto* __restrict__ rx = cells[i][ field::rx ]; ONIKA_ASSUME_ALIGNED(rx);
          auto* __restrict__ ry = cells[i][ field::ry ]; ONIKA_ASSUME_ALIGNED(ry);
          auto* __restrict__ rz = cells[i][ field::rz ]; ONIKA_ASSUME_ALIGNED(rz);

          for(size_t k=0;k<n;k++)
          {
			rx[k] += gaussian(re);
			ry[k] += gaussian(re);
			rz[k] += gaussian(re);
          }
        }
        GRID_OMP_FOR_END
      }
    }

  };

  template<class GridT> using SliplinkForceOverdampedTmpl = SliplinkForceOverdamped<GridT>;
  
 // === register factories ===  
  ONIKA_AUTORUN_INIT(sliplink_bead_friction)
  {
   OperatorNodeFactory::instance()->register_factory( "sliplink_bead_friction", make_grid_variant_operator< SliplinkForceOverdampedTmpl > );
  }

}

