#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/quantity_yaml.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/source_term.h>

#include <memory>

namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class InitTTM : public OperatorNode
  {
    ADD_SLOT( GridT       , grid         , INPUT , REQUIRED );
    ADD_SLOT( Domain         , domain       , INPUT , REQUIRED );
    ADD_SLOT( double         , physical_time, INPUT , REQUIRED );

    ADD_SLOT( ScalarSourceTermInstance , te_source  , INPUT_OUTPUT );
    ADD_SLOT( ScalarSourceTermInstance , ti_source  , INPUT_OUTPUT );
    
    ADD_SLOT( long           , grid_subdiv      , INPUT , 3 );
    ADD_SLOT( GridCellValues , grid_cell_values , INPUT_OUTPUT );

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {           
      // retreive field data accessor. create data field if needed
      const int subdiv = *grid_subdiv;
      if( ! grid_cell_values->has_field("te") )
      {
        grid_cell_values->add_field("te",subdiv,1);
      }
      assert( size_t(subdiv) == grid_cell_values->field("te").m_subdiv );
      assert( size_t(subdiv * subdiv * subdiv) == grid_cell_values->field("te").m_components );
      auto cell_te_data = grid_cell_values->field_data("te");
      
      const Mat3d xform = domain->xform();
      const double subcell_size = domain->cell_size() / subdiv;
      const IJK dims = grid->dimension();
      
      const auto& te_source_func = * (*te_source);
      
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,cell_i,cell_loc, schedule(static) )
        {
	        const Vec3d cell_origin = grid->cell_position( cell_loc );

          for(int ck=0;ck<subdiv;ck++)
          for(int cj=0;cj<subdiv;cj++)
          for(int ci=0;ci<subdiv;ci++)
          {
            IJK sc { ci, cj, ck };
	          Vec3d scr = { ci+0.5, cj+0.5, ck+0.5 };
            const size_t j = cell_i*cell_te_data.m_stride +  grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , sc );
            const Vec3d center = xform * ( cell_origin + scr * subcell_size );
            cell_te_data.m_data_ptr[j] = te_source_func ( center, *physical_time );
          }
        }
        GRID_OMP_FOR_END
      }

    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
Initializes cell_te grid values and source terms for ttm model
)EOF";
    }

  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory("init_ttm", make_grid_variant_operator< InitTTM > );
  }

}
