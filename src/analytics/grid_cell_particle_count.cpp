#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/domain.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/i64_as_double.h>

namespace exaStamp
{

  template<class GridT>
  class GridCellParticleCount : public OperatorNode
  {
    ADD_SLOT( GridT          , grid            , INPUT , REQUIRED );
    ADD_SLOT( Domain         , domain          , INPUT , REQUIRED );
    ADD_SLOT( GridCellValues , grid_cell_values , INPUT_OUTPUT );
    ADD_SLOT( std::string    , field_name     , INPUT , "pc" ); // particle count field name
    ADD_SLOT( long           , grid_subdiv     , INPUT , 1 );

  public:
    inline void execute () override final
    {      
      //const IJK dom_dims = domain->grid_dimension();
      const IJK dims = grid->dimension();
      //const IJK grid_offset = grid->offset();
      const auto cells = grid->cells();
      //const ssize_t n_cells = grid->number_of_cells();

      const double cell_size = domain->cell_size();
      const int subdiv = *grid_subdiv;
      const double subcell_size = cell_size / subdiv;
      //const int n_subcells = subdiv * subdiv * subdiv;

//      const int n_subcells = subdiv * subdiv * subdiv;
      if( ! grid_cell_values->has_field(*field_name) )
      {
        ldbg << "create grid cell field "<< *field_name << std::endl;
        grid_cell_values->add_field(*field_name,subdiv,1);
      }
      assert( subdiv == static_cast<ssize_t>(grid_cell_values->field(*field_name).m_subdiv) );
      assert( (subdiv * subdiv * subdiv) == static_cast<ssize_t>(grid_cell_values->field(*field_name).m_components) );

      auto field_data = grid_cell_values->field_data(*field_name);
      auto* __restrict__ pc_ptr = field_data.m_data_ptr;
      const size_t pc_stride = field_data.m_stride;

#     ifndef NDEBUG
      ssize_t n_cells = grid->number_of_cells();
      ssize_t n_subcells = subdiv * subdiv * subdiv;
#     endif

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,i,cell_loc, schedule(dynamic) )
        {
	  const Vec3d cell_origin = grid->cell_position( cell_loc );
          const size_t n = cells[i].size();
          const double* __restrict__ rx = cells[i][field::rx]; ONIKA_ASSUME_ALIGNED(rx);
          const double* __restrict__ ry = cells[i][field::ry]; ONIKA_ASSUME_ALIGNED(ry);
          const double* __restrict__ rz = cells[i][field::rz]; ONIKA_ASSUME_ALIGNED(rz);
          for(size_t j=0;j<n;j++)          
          {
            const Vec3d r { rx[j] , ry[j] , rz[j] };

            IJK center_cell_loc;
            IJK center_subcell_loc;
            Vec3d rco = r - cell_origin;
            localize_subcell( rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc );
            center_cell_loc += cell_loc;

            ssize_t center_cell_i = grid_ijk_to_index( dims , center_cell_loc );
            ssize_t center_subcell_i = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , center_subcell_loc );
            assert( center_cell_i>=0 && center_cell_i<n_cells );
            assert( center_subcell_i>=0 && center_subcell_i<n_subcells );
            size_t scindex = center_cell_i * pc_stride + center_subcell_i;

#           pragma omp atomic update
            pc_ptr[ scindex ] += 1.0; // FIXME: data race here
          }
        }
        GRID_OMP_FOR_END
      }

    }

  private:
    static inline void localize_subcell( const Vec3d& r, double cell_size, double sub_cellsize, ssize_t subdiv, IJK& cell_loc, IJK& subcell_loc )
    {
      cell_loc = make_ijk( r / cell_size );
      Vec3d ro = r - (cell_loc*cell_size);
      subcell_loc = vclamp( make_ijk(ro / sub_cellsize) , 0 , subdiv-1 );
    }

  };

  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "grid_cell_particle_count", make_grid_variant_operator<GridCellParticleCount> );
  }

}

