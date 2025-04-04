#include <memory>

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/physics/units.h>
#include <onika/physics/units.h>
#include <onika/memory/allocator.h>

#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types_stream.h>
#include <exaStamp/npt/npt.h>
#include <exanb/core/domain.h>
#include <onika/physics/constants.h>
#include <exaStamp/compute/thermodynamic_state.h>

#include <iostream>
#include <string>

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_ax, field::_ay, field::_az, field::_vx, field::_vy, field::_vz >
    >
  
  struct NHVPressNode : public OperatorNode
  {
    using PointerTuple = onika::soatl::FieldPointerTuple<
      GridT::CellParticles::Alignment , GridT::CellParticles::ChunkSize , 
      field::_rx, field::_ry, field::_rz,
      field::_vx, field::_vy, field::_vz,
      field::_ax, field::_ay, field::_az >;
    
    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;

    using has_id_field_t = typename GridT::CellParticles::template HasField < field::_id > ;
    static constexpr bool has_id_field = has_id_field_t::value;
    
    ADD_SLOT( GridT          , grid    , INPUT_OUTPUT);
    ADD_SLOT( NPTContext              , npt_ctx    , INPUT_OUTPUT );
    ADD_SLOT( Vec3d                   , vscale     , OUTPUT );
    ADD_SLOT( Vec3d                   , vadd       , OUTPUT );

    inline void execute () override final
    {

      if( grid->number_of_cells() == 0 ) return;
      
      GridT& grid = *(this->grid);
      auto cells = grid.cells();
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();      
      Vec3d factor;
      factor.x = exp(-npt_ctx->dt4*(npt_ctx->omega_dot[0]+npt_ctx->mtk_term2));
      factor.y = exp(-npt_ctx->dt4*(npt_ctx->omega_dot[1]+npt_ctx->mtk_term2));
      factor.z = exp(-npt_ctx->dt4*(npt_ctx->omega_dot[2]+npt_ctx->mtk_term2));

#     pragma omp parallel
      {
        PointerTuple ptrs;   	
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
          {
            size_t i = grid_ijk_to_index( dims , loc + gl );
            const unsigned int n = cells[i].size();
            cells[i].capture_pointers( ptrs );
            
            auto* __restrict__ vx = cells[i][field::vx];
            auto* __restrict__ vy = cells[i][field::vy];
            auto* __restrict__ vz = cells[i][field::vz];
            
            for(unsigned int j=0;j<n;j++)
              {
                vx[j] *= factor.x;
                vy[j] *= factor.y;
                vz[j] *= factor.z;
                if (npt_ctx->pstyle == "TRICLINIC") {
                  vx[j] -= npt_ctx->dthalf*(vy[j]*npt_ctx->omega_dot[5] + vz[j]*npt_ctx->omega_dot[4]);
                  vy[j] -= npt_ctx->dthalf*vz[j]*npt_ctx->omega_dot[3];
                }
                vx[j] *= factor.x;
                vy[j] *= factor.y;
                vz[j] *= factor.z;		
              }
          }
        GRID_OMP_FOR_END
          }
      
    }
  };

  template<class GridT> using NHVPressNodeTmpl = NHVPressNode<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(nh_v_press)
  {
   OperatorNodeFactory::instance()->register_factory( "nh_v_press", make_grid_variant_operator< NHVPressNodeTmpl > );
  }

}
