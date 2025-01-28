#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/grid.h>
#include <onika/math/basic_types_stream.h>
#include <onika/log.h>
#include <exanb/core/domain.h>
#include <exanb/core/quantity.h>
#include <exanb/core/string_utils.h>
#include <exanb/core/grid_fields.h>

#include <exanb/defbox/deformation.h>
#include <exanb/defbox/deformation_stream.h>
#include <exanb/defbox/deformation_yaml.h>
#include <exanb/defbox/deformation_math.h>

#include <exaStamp/mechanical/cell_particles_local_metrics.h>

#include <onika/soatl/packed_field_arrays.h>
#include <memory>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <mpi.h>
#include <string>
#include <iomanip>
#include <experimental/filesystem>
#include <math.h>
#include <algorithm>

namespace exaStamp
{
  using namespace exanb;
  
  template< class GridT
           >
  class ComputeLocalMetricsOperator : public OperatorNode
  {
    //field id
    using has_id_field_t = typename GridT::CellParticles::template HasField < field::_id > ;
    static constexpr bool has_id_field = has_id_field_t::value;
    //field type
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;
    //field charge
    using has_charge_field_t = typename GridT::CellParticles::template HasField < field::_charge > ;
    static constexpr bool has_charge_field = has_charge_field_t::value;
    //field virial
    using has_virial_field_t = typename GridT::CellParticles::template HasField < field::_virial > ;
    static constexpr bool has_virial_field = has_virial_field_t::value;

    // fore field
    using has_forcex_field_t = typename GridT::CellParticles::template HasField < field::_fx > ;
    using has_forcey_field_t = typename GridT::CellParticles::template HasField < field::_fy > ;
    using has_forcez_field_t = typename GridT::CellParticles::template HasField < field::_fz > ;
    static constexpr bool has_force_field = has_forcex_field_t::value && has_forcey_field_t::value && has_forcez_field_t::value;

    // fore field
    using has_vx_field_t = typename GridT::CellParticles::template HasField < field::_vx > ;
    using has_vy_field_t = typename GridT::CellParticles::template HasField < field::_vy > ;
    using has_vz_field_t = typename GridT::CellParticles::template HasField < field::_vz > ;
    static constexpr bool has_velocity_field = has_vx_field_t::value && has_vy_field_t::value && has_vz_field_t::value;

    //field ep
    using has_ep_field_t = typename GridT::CellParticles::template HasField < field::_ep > ;
    static constexpr bool has_ep_field = has_ep_field_t::value;

    //field rxf
    using has_rxf_field_t = typename GridT::CellParticles::template HasField < field::_rxf > ;
    static constexpr bool has_rxf_field = has_rxf_field_t::value;
    
    //field ryf
    using has_ryf_field_t = typename GridT::CellParticles::template HasField < field::_ryf > ;
    static constexpr bool has_ryf_field = has_ryf_field_t::value;
    
    //field rzf
    using has_rzf_field_t = typename GridT::CellParticles::template HasField < field::_rzf > ;
    static constexpr bool has_rzf_field = has_rzf_field_t::value;
    
    using VariablesVec = std::vector<std::string>;
    
    ADD_SLOT( MPI_Comm                , mpi                 , INPUT );
    ADD_SLOT( GridT                   , grid                , INPUT );
    ADD_SLOT( Domain                  , domain              , INPUT );
    ADD_SLOT( Deformation             , defbox              , INPUT );    
    ADD_SLOT( long                    , timestep            , INPUT );
    ADD_SLOT( bool                    , is_ghosts           , INPUT , false);
    ADD_SLOT( GridParticleLocalMetrics, local_data          , OUTPUT );
    ADD_SLOT( VariablesVec            , per_atom_data       , INPUT );    
  public:
  
    inline void execute () override final
    {      
      GridT& grid = *(this->grid);
      bool is_ghosts = *(this->is_ghosts);
      
      size_t n_cells = grid.number_of_cells();
      auto cells = grid.cells();      

      GridParticleLocalMetrics& local_data = *(this->local_data);
      local_data.resize( n_cells );

      VariablesVec per_atom_data = *(this->per_atom_data);
      std::sort(per_atom_data.begin(), per_atom_data.end());
      
      // ****************************************************************************** //
      // output of per-atom potential energy
      // ****************************************************************************** //      
      if constexpr (has_ep_field) if (std::binary_search(per_atom_data.begin(), per_atom_data.end(), "pe"))
	    {
	      lout << "\t- Gathering per-atom potential energy" << std::endl;
	      for(size_t c=0; c<n_cells;++c)
	        {
	          local_data[c].pe.resize( cells[c].size() );
	          int np = 0;
	          if( !grid.is_ghost_cell(c) || is_ghosts )
		        {
		          const auto * __restrict__ ep = cells[c][field::ep];		  
		          np = cells[c].size();
		          for(int p=0;p<np;++p) local_data[c].pe[p] = ep[p];
		        }
	        }
	    }

      // ****************************************************************************** //
      // output of per-atom force vector
      // ****************************************************************************** //      
      if constexpr (has_force_field) if (std::binary_search(per_atom_data.begin(), per_atom_data.end(), "f"))
	    {
	      lout << "\t- Gathering per-atom force vector" << std::endl;
	      for(size_t c=0; c<n_cells;++c)
	        {
	          local_data[c].f.resize( cells[c].size() );
	          int np = 0;
	          if( !grid.is_ghost_cell(c) || is_ghosts )
		        {
		          const auto * __restrict__ fx = cells[c][field::fx];
		          const auto * __restrict__ fy = cells[c][field::fy];
		          const auto * __restrict__ fz = cells[c][field::fz];		  		  
		          np = cells[c].size();
		          for(int p=0;p<np;++p) local_data[c].f[p] = Vec3d{fx[p], fy[p], fz[p]};
		        }
	        }
	    }

      // ****************************************************************************** //
      // output of per-atom velocity vector
      // ****************************************************************************** //      
      if constexpr (has_velocity_field) if (std::binary_search(per_atom_data.begin(), per_atom_data.end(), "v"))
	    {
	      lout << "\t- Gathering per-atom velocity vector" << std::endl;
	      for(size_t c=0; c<n_cells;++c)
	        {
	          local_data[c].v.resize( cells[c].size() );
	          int np = 0;
	          if( !grid.is_ghost_cell(c) || is_ghosts )
		        {
		          const auto * __restrict__ vx = cells[c][field::vx];
		          const auto * __restrict__ vy = cells[c][field::vy];
		          const auto * __restrict__ vz = cells[c][field::vz];		  		  
		          np = cells[c].size();
		          for(int p=0;p<np;++p) local_data[c].v[p] = Vec3d{vx[p], vy[p], vz[p]};
		        }
	        }
	    }

      // ****************************************************************************** //
      // output of per-atom charge
      // ****************************************************************************** //      
      if (std::binary_search(per_atom_data.begin(), per_atom_data.end(), "q"))
	    {
	      lout << "\t- Gathering per-atom charge" << std::endl;
	      for(size_t c=0; c<n_cells;++c)
	        {
	          local_data[c].q.resize( cells[c].size() );
	          int np = 0;
	          if( !grid.is_ghost_cell(c) || is_ghosts )
		        {
		          const auto * __restrict__ q = cells[c][field::charge];
		          np = cells[c].size();
		          for(int p=0;p<np;++p) local_data[c].q[p] = q[p];
		        }
	        }
	    }

      // ****************************************************************************** //
      // output of per-atom virial
      // ****************************************************************************** //      
      if constexpr (has_virial_field) if (std::binary_search(per_atom_data.begin(), per_atom_data.end(), "virial"))
	{
	  lout << "\t- Gathering per-atom virial" << std::endl;
	  for(size_t c=0; c<n_cells;++c)
	    {
	      local_data[c].virial.resize( cells[c].size() );
	      int np = 0;
	      if( !grid.is_ghost_cell(c) || is_ghosts )
		{
		  const auto * __restrict__ virial = cells[c][field::virial];
		  np = cells[c].size();
		  for(int p=0;p<np;++p) local_data[c].virial[p] = virial[p];
		}
	    }
	}

      // ****************************************************************************** //
      // output of per-atom atom id
      // ****************************************************************************** //      
      if constexpr (has_id_field) if (std::binary_search(per_atom_data.begin(), per_atom_data.end(), "id"))
	{
	  lout << "\t- Gathering per-atom id" << std::endl;
	  for(size_t c=0; c<n_cells;++c)
	    {
	      local_data[c].id.resize( cells[c].size() );
	      int np = 0;
	      if( !grid.is_ghost_cell(c) || is_ghosts )
		{
		  const auto * __restrict__ id = cells[c][field::id];
		  np = cells[c].size();
		  for(int p=0;p<np;++p) local_data[c].id[p] = id[p];
		}
	    }
	}

      // ****************************************************************************** //
      // output of per-atom type
      // ****************************************************************************** //      
      if constexpr (has_type_field) if (std::binary_search(per_atom_data.begin(), per_atom_data.end(), "type"))
	{
	  lout << "\t- Gathering per-atom type" << std::endl;
	  for(size_t c=0; c<n_cells;++c)
	    {
	      local_data[c].type.resize( cells[c].size() );
	      int np = 0;
	      if( !grid.is_ghost_cell(c) || is_ghosts )
		{
		  const auto * __restrict__ type = cells[c][field::type];
		  np = cells[c].size();
		  for(int p=0;p<np;++p) local_data[c].type[p] = type[p];
		}
	    }
	}
      
    }
};

  template<class GridT> using ComputeLocalMetricsOperatorTmpl = ComputeLocalMetricsOperator<GridT>;
  
  // === register factories ===  
  ONIKA_AUTORUN_INIT(compute_local_metrics)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_local_metrics", make_grid_variant_operator< ComputeLocalMetricsOperatorTmpl > );
  }

}
