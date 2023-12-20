#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/domain.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/grid_algorithm.h>

#include <mpi.h>
#include <onika/silent_use.h>

namespace exaStamp
{
  using namespace exanb;

  class GridConnectedComponents : public OperatorNode
  {
    ADD_SLOT( MPI_Comm       , mpi                    , INPUT , MPI_COMM_WORLD , DocString{"MPI communicator"} );
    ADD_SLOT( Domain         , domain                 , INPUT , REQUIRED );
    ADD_SLOT( double         , density_threshold      , INPUT , 1. , DocString{"density threshold. below this threshold is considered empty"} );
    ADD_SLOT( GridCellValues , grid_cell_values       , INPUT_OUTPUT );
    
  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      MPI_Comm comm = *mpi; ONIKA_SILENT_USE(comm);
    
      // cell size (alle cells are cubical)
      const double cell_size = domain->cell_size(); 
      
      // local processor's grid dimensions, including ghost cells
      const IJK grid_dims = grid_cell_values->grid_dims(); 
      
      // local processor's grid position in simulation's grid
      const IJK grid_offset = grid_cell_values->grid_offset();
      
      // simulation's grid dimensions
      const IJK domain_dims = domain->grid_dimension(); 
      
       // ghost layers
      const ssize_t gl = grid_cell_values->ghost_layers();
      
      // number of subdivisions, in each directions, applied on cells
      const ssize_t subdiv = grid_cell_values->field("density").m_subdiv;
      
      // number of sub-cells per cell
      const ssize_t n_subcells = subdiv*subdiv*subdiv; ONIKA_SILENT_USE(n_subcells);
      
      // side size of a sub-cell
      const double subcell_size = cell_size / subdiv;
      
      // dimension of the subdivided simulation's grid
      const IJK domain_subdiv_dims = domain_dims * subdiv; 

      // some debug information
      ldbg << "ghost_layers="<<gl<<", cell_size="<<cell_size<<", subdiv="<<subdiv<<", subcell_size="
           <<subcell_size<<", grid_dims="<<grid_dims<<", grid_offset="<<grid_offset<<", domain_dims="<<domain_dims<<", domain_subdiv_dims="<<domain_subdiv_dims<<std::endl;

      // note: if we are to add new data fields, they must be added BEFORE we retreive access information

      // create additional data field for connected component label
      if( ! grid_cell_values->has_field("cc_label") )
      {
        grid_cell_values->add_field("cc_label",subdiv,1);
      }
      
      // retreive cc_label field data accessor.
      const auto cc_label_accessor = grid_cell_values->field_data("cc_label");
      double * __restrict__ cc_label_ptr = cc_label_accessor.m_data_ptr;
      const size_t cc_label_stride = cc_label_accessor.m_stride; ONIKA_SILENT_USE(cc_label_stride);

      // retreive density field data accessor.
      const auto density_accessor = grid_cell_values->field_data("density");
      const double  * __restrict__ density_ptr = density_accessor.m_data_ptr;
      const size_t density_stride = density_accessor.m_stride;
      
      // sanity check
      assert( cc_label_stride == density_stride );
      const size_t stride = density_stride; // for simplicity

      // triple loop to enumerate grid cells, excluding ghost layers  
      for( ssize_t k=0 ; k < grid_dims.k ; k++)
      for( ssize_t j=0 ; j < grid_dims.j ; j++)
      for( ssize_t i=0 ; i < grid_dims.i ; i++)
      {
        // informations that might be usefull later on ...
        // tells if cell is in the ghost layers ()
        // bool is_ghost = i<gl || i>=(grid_dims.i-gl) || j<gl || j>=(grid_dims.j-gl) || k<gl || k>=(grid_dims.k-gl);
        
        // position of the cell in the simulation grid, which size is 'domain_dims'
        IJK cell_location = IJK{i,j,k} + grid_offset;

        // triple loop to enumerate sub cells inside a cell
        for( ssize_t sk=0 ; sk<subdiv ; sk++)
        for( ssize_t sj=0 ; sj<subdiv ; sj++)
        for( ssize_t si=0 ; si<subdiv ; si++)
        {
          // location of the subcell in simulation's subcell grid, which dimension is 'domain_subdiv_dims'
          IJK subcell_global_location = cell_location * subdiv + IJK{si,sj,sk};
          ssize_t subcell_global_index = grid_ijk_to_index( domain_subdiv_dims , subcell_global_location ); if(subcell_global_index<0){ lerr<<"Bad index\n"; }

          // computation of subcell index
          ssize_t cell_index = grid_ijk_to_index( grid_dims , IJK{i,j,k} );
          ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , IJK{si,sj,sk} );
          // value index, is the index of current subcell for local processor's grid
          ssize_t value_index = cell_index * stride + subcell_index; if(value_index<0){ lerr<<"Bad index\n"; }

          // compute a label to assign to subcell, if it is dense enough.
          // TODO: change this so it is invariant from domain decomposition
          double label = double(value_index);
          
          // assign a label to cells where density is above given threshold, otherwise assign -1
          cc_label_ptr[ value_index ] = density_ptr[ value_index ] > (*density_threshold) ? label : -1.0;
        }
 
      }
  
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
Computes cell connected components information
)EOF";
    }

  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory("grid_connected_components", make_simple_operator< GridConnectedComponents > );
  }

}
