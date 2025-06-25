
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_stream.h>
#include <exanb/grid_cell_particles/particle_region.h>
#include <exaStamp/particle_species/particle_specie_yaml.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <iostream>
#include <string>
#include <algorithm>
#include <mpi.h>

namespace exaStamp
{

  using namespace exanb;

  /*
  WARNING: only one region can be tracked by now. algothim works and guarantees contiguous particle ids all over the domain,
  BUT it cannot conserve original id order or values, nor other tracked region id intervals
  */
  template<class GridT , class = AssertGridHasFields<GridT,field::_id,field::_type> >
  class TrackSpeciesAsRegions : public OperatorNode
  {  
    ADD_SLOT( MPI_Comm          , mpi    , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( GridT             , grid   , INPUT_OUTPUT);

    ADD_SLOT( ParticleSpecies   , species , INPUT , REQUIRED );
    ADD_SLOT( ParticleRegions   , particle_regions , INPUT_OUTPUT , ParticleRegions{} );

  public:
    inline void execute() override final
    {    
      int nprocs = 1;
      int rank = 0;
      MPI_Comm_rank(*mpi,&rank);
      MPI_Comm_size(*mpi,&nprocs);

      size_t n_cells = grid->number_of_cells();
      auto cells = grid->cells();

      std::vector<unsigned long long> local_type_count( species->size() , 0 );
      std::vector<unsigned long long> global_type_count( species->size() , 0 );
      std::vector<unsigned long long> global_type_offset( species->size() , 0 );

      for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( ! grid->is_ghost_cell(cell_i) )
      {
        size_t n_particles = cells[cell_i].size();
        for(size_t p=0;p<n_particles;p++)
        {
          const unsigned int type = cells[cell_i][field::type][p];
          auto & tcount = local_type_count[type];
#         pragma omp atomic update
          ++ tcount;
        }
      }

      MPI_Allreduce(local_type_count.data(), global_type_count.data()  , species->size(), MPI_UNSIGNED_LONG_LONG, MPI_SUM, *mpi);      
      MPI_Exscan(   local_type_count.data(), global_type_offset.data() , species->size(), MPI_UNSIGNED_LONG_LONG, MPI_SUM, *mpi);

      size_t all_type_count = 0;
      for(size_t i=0;i<species->size();i++)
      {
        global_type_offset[i] += all_type_count;
        all_type_count += global_type_count[i];
        ldbg << species->at(i).name()<<" count = " << global_type_count[i] << std::endl;
        ldbg << species->at(i).name()<<" offset = " << global_type_offset[i] << std::endl;
      }

      ldbg << "All types count = " << all_type_count << std::endl;

#     pragma omp parallel for schedule(dynamic)
      for(size_t cell_i=0;cell_i<n_cells;cell_i++) if( ! grid->is_ghost_cell(cell_i) )
      {
        size_t n_particles = cells[cell_i].size();
        for(size_t p=0;p<n_particles;p++)
        {
          const unsigned int type = cells[cell_i][field::type][p];
          auto & tcount = global_type_offset[type];
          cells[cell_i][field::id][p] = tcount;
#         pragma omp atomic update
          ++ tcount;
        }
      }

      all_type_count = 0;
      for(size_t i=0;i<species->size();i++)
      {
        ParticleRegion r = {};
        r.set_name( species->at(i).name() );
        r.m_id_start = all_type_count;
        all_type_count += global_type_count[i];
        r.m_id_end = all_type_count;
        r.m_id_range_flag = true;
        ldbg << "add tracking region "<<species->at(i).name()<<" with "<< (r.m_id_end - r.m_id_start) <<" particles starting at id="<<r.m_id_start<<std::endl;
        particle_regions->push_back( r );
      }
    }

  };

  template<class GridT> using TrackSpeciesAsRegionsTmpl = TrackSpeciesAsRegions< GridT >;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(track_species_regions)
  {
    OperatorNodeFactory::instance()->register_factory( "track_species_regions", make_grid_variant_operator<TrackSpeciesAsRegionsTmpl> );
  }

}

