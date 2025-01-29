#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <onika/memory/allocator.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/particle_neighbors/chunk_neighbors_iterator.h>
#include <exanb/particle_neighbors/chunk_neighbors_stream_check.h>

#include <onika/string_utils.h>
#include <exanb/core/particle_type_pair.h>
#include <exaStamp/particle_species/particle_specie_yaml.h>
#include <exanb/core/particle_type_id.h>
#include <exaStamp/particle_species/distance_map.h>

#include <iomanip>
#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  template<typename GridT , class = AssertGridHasFields< GridT, field::_type > >
  class AtomNeighborsStats : public OperatorNode
  {  
    ADD_SLOT( MPI_Comm           , mpi             , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( GridT              , grid            , INPUT , REQUIRED );
    ADD_SLOT( GridChunkNeighbors , chunk_neighbors , INPUT , GridChunkNeighbors{} );
    ADD_SLOT( double             , nbh_dist        , INPUT , 0.0 );
    ADD_SLOT( bool               , ghost           , INPUT , false );
    ADD_SLOT( ParticleSpecies    , species         , INPUT , REQUIRED );
    ADD_SLOT( ParticleTypeMap    , particle_type_map  , INPUT , REQUIRED);
    
    ADD_SLOT( DistanceMap        , pair_distances  , INPUT , OPTIONAL , DocString{"optional set of distances. if set, more statistics about neighbor distances will be accounted, for each distance threshold in 'distances'"} );

  public:
    inline void execute () override final
    {  
      if( grid->number_of_cells() == 0 ) return;
      assert( grid->number_of_cells() == chunk_neighbors->number_of_cells() );
    
      ldbg << "atoms = "<<grid->number_of_particles()<<", cells = "<<grid->number_of_cells()<<std::endl;

      size_t cs = chunk_neighbors->m_chunk_size;
      size_t cs_log2 = 0;
      while( cs > 1 )
      {
        assert( (cs&1)==0 );
        cs = cs >> 1;
        ++ cs_log2;
      }
      cs = chunk_neighbors->m_chunk_size;
      ldbg << "cs="<<cs<<", log2(cs)="<<cs_log2<<std::endl;

      // pair distance stats
      const unsigned int ntypes = species->size();
      const unsigned int npairs = unique_pair_count( ntypes );
      const double nbh_d2 = (*nbh_dist) * (*nbh_dist) ;
      const IJK dims = grid->dimension();
      unsigned long long total_particles = grid->number_of_particles();

      std::vector<unsigned long long> pair_count( npairs , 0 );
      std::vector<unsigned long long> pair_dist_count_min( npairs , total_particles );
      std::vector<unsigned long long> pair_dist_count_max( npairs , 0 );
      
      std::vector<double> pairdist( npairs , 0.0 );
      std::vector<double> pair_min( npairs , std::numeric_limits<double>::max() );
      std::vector<double> pair_max( npairs , 0.0 );
      std::vector<double> pair_sum( npairs , 0.0 );
      ldbg << "atom types = "<<ntypes<<", npairs = "<<npairs<<std::endl;

      if( pair_distances.has_value() )
      {
        ldbg << "found " << pair_distances->size() << " pair distances"<<std::endl;
        for(const auto& pdist : *pair_distances)
        {
          std::string ts1 = pdist.first.substr(0,pdist.first.find(','));
          std::string ts2 = pdist.first.substr(pdist.first.find(',')+1);
          const unsigned int t1 = particle_type_map->at( ts1 );
          const unsigned int t2 = particle_type_map->at( ts2 );
          const unsigned int pair_id = unique_pair_id(t1,t2);
          ldbg << ts1 << " / " << ts2 << " -> "<<pdist.second<<std::endl;
          pairdist[pair_id] = pdist.second * pdist.second;
        }
      }
      else
      {
        pairdist.assign( npairs , nbh_d2 );
      }

      unsigned long long total_nbh = 0;
      unsigned long long total_nbh_d2 = 0;
      //size_t total_nbh_chunk = 0;
      unsigned long long total_nbh_cells = 0;
      
      unsigned long long nbh_d2_min = total_particles;
      unsigned long long nbh_d2_max = 0;

      unsigned long long nbh_min = total_particles;
      unsigned long long nbh_max = 0;
      total_particles = 0;

      auto cells = grid->cells();
      using CellT = std::remove_cv_t< std::remove_reference_t< decltype(cells[0]) > >;
      ChunkParticleNeighborsIterator<CellT> chunk_nbh_it_in = { grid->cells() , chunk_neighbors->data() , dims , chunk_neighbors->m_chunk_size };

#     pragma omp parallel
      {
        std::vector<unsigned long long> local_pair_count( npairs , 0 );
        std::vector<unsigned long long> local_pair_dist_count_min( npairs , total_particles );
        std::vector<unsigned long long> local_pair_dist_count_max( npairs , 0 );
        std::vector<double> local_pair_min( npairs , std::numeric_limits<double>::max() );
        std::vector<double> local_pair_max( npairs , 0.0 );
        std::vector<double> local_pair_sum( npairs , 0.0 );
        std::vector<unsigned long long> p_a_nbh_d2_pair( npairs , 0 );

#       ifndef NDEBUG
        std::vector<unsigned int> cell_a_particle_chunk_count;
#       endif

        auto chunk_nbh_it = chunk_nbh_it_in;

        GRID_OMP_FOR_BEGIN(dims,cell_a,loc_a, schedule(dynamic) \
                                              reduction(+:total_particles,total_nbh,total_nbh_d2,total_nbh_cells) \
                                              reduction(min:nbh_d2_min,nbh_min) \
                                              reduction(max:nbh_d2_max,nbh_max) )
        {
          // std::cout<<"dims="<<dims<<" cell_a="<<cell_a<<" loc_a="<<loc_a<<std::endl;
          size_t n_particles_a = cells[cell_a].size();

#         ifndef NDEBUG
          cell_a_particle_chunk_count.clear();
          chunk_neighbors_stream_check( *grid, chunk_neighbors->m_chunk_size, cell_a, chunk_neighbors->cell_stream(cell_a), chunk_neighbors->cell_stream_size(cell_a), cell_a_particle_chunk_count );
          assert( cell_a_particle_chunk_count.size() == n_particles_a );
#         endif

          const double* __restrict__ rx_a = cells[cell_a][field::rx];
          const double* __restrict__ ry_a = cells[cell_a][field::ry];
          const double* __restrict__ rz_a = cells[cell_a][field::rz];
          const uint8_t* __restrict__ type_a = cells[cell_a][field::type];

          const double* __restrict__ rx_b = nullptr; 
          const double* __restrict__ ry_b = nullptr;
          const double* __restrict__ rz_b = nullptr;
          const uint8_t* __restrict__ type_b = nullptr;

          // decode compacted chunks
          chunk_nbh_it.start_cell( cell_a , n_particles_a );
          for(size_t p_a=0;p_a<n_particles_a;p_a++)
          {
            if( (*ghost) ||  ! grid->is_ghost_cell(loc_a) )
            {
              ++ total_particles;
            }
            //size_t total_nbh_before = total_nbh;
            unsigned long long p_a_nbh_d2 = 0;
            unsigned long long p_a_nbh = 0;
            p_a_nbh_d2_pair.assign( npairs , 0 );

            chunk_nbh_it.start_particle( p_a );
            size_t last_cell = std::numeric_limits<size_t>::max();
            while( ! chunk_nbh_it.end() )
            {
              size_t cell_b=0, p_b=0;
              chunk_nbh_it.get_nbh( cell_b , p_b );
              //std::cout<<"C"<<cell_a<<"P"<<p_a<<" -> C"<<cell_b<<"P"<<p_b<<std::endl;
              if( cell_b != last_cell )
              {
                rx_b = cells[cell_b][field::rx];
                ry_b = cells[cell_b][field::ry];
                rz_b = cells[cell_b][field::rz];
                type_b = cells[cell_b][field::type];
                last_cell = cell_b;
                if( (*ghost) || ! grid->is_ghost_cell(loc_a) )
                {
                  ++ total_nbh_cells;
                }
              }
              
              const double dx = rx_b[p_b] - rx_a[p_a];
              const double dy = ry_b[p_b] - ry_a[p_a];
              const double dz = rz_b[p_b] - rz_a[p_a];
              const double d2 = dx*dx+dy*dy+dz*dz;
              const double d = sqrt(d2);                            
              const unsigned int pair_id = unique_pair_id(type_a[p_a],type_b[p_b]);
              assert( pair_id < npairs );
              
              if( d2 <= nbh_d2 )
      	      {	      
                ++ local_pair_count[pair_id];
                local_pair_sum[pair_id] += d;
                local_pair_min[pair_id] = std::min( local_pair_min[pair_id] , d );
                local_pair_max[pair_id] = std::max( local_pair_min[pair_id] , d );
      	      }
              
              ++ p_a_nbh;
              if( d2 <= nbh_d2 )
              {
                ++ p_a_nbh_d2;
              }
              if( d2 <= pairdist[pair_id])
              {
                ++ p_a_nbh_d2_pair[pair_id];
              }
              
              chunk_nbh_it.next();
            }
                        
            if( (*ghost) ||  ! grid->is_ghost_cell(loc_a) )
            {
              total_nbh_d2 += p_a_nbh_d2;
              total_nbh += p_a_nbh;

              nbh_d2_min = std::min( nbh_d2_min , p_a_nbh_d2 );
              nbh_d2_max = std::max( nbh_d2_max , p_a_nbh_d2 );

              nbh_min = std::min( nbh_min , p_a_nbh );
              nbh_max = std::max( nbh_min , p_a_nbh );
              
              for(unsigned int i=0;i<npairs;i++)
              {
                local_pair_dist_count_min[i] = std::min( local_pair_dist_count_min[i] , p_a_nbh_d2_pair[i] );
                local_pair_dist_count_max[i] = std::max( local_pair_dist_count_max[i] , p_a_nbh_d2_pair[i] );
              }
            }
          }
                    
        }
        GRID_OMP_FOR_END

#       pragma omp critical(reduce_pair_distance_stats)
        {
          for(unsigned int i=0;i<npairs;i++)
          {
            pair_count[i] += local_pair_count[i];
            pair_dist_count_min[i] = std::min( pair_dist_count_min[i] , local_pair_dist_count_min[i] );
            pair_dist_count_max[i] = std::max( pair_dist_count_max[i] , local_pair_dist_count_max[i] );
            pair_sum[i] += local_pair_sum[i];
            pair_min[i] = std::min( pair_min[i] , local_pair_min[i] );
            pair_max[i] = std::max( pair_max[i] , local_pair_max[i] );
          }
        }

      } // end of parallel section

      MPI_Allreduce(MPI_IN_PLACE,&total_particles,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&total_nbh,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&total_nbh_d2,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&total_nbh_cells,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,*mpi);

      MPI_Allreduce(MPI_IN_PLACE,&nbh_d2_min,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&nbh_min,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,*mpi);

      MPI_Allreduce(MPI_IN_PLACE,&nbh_d2_max,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&nbh_max,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,*mpi);

      MPI_Allreduce(MPI_IN_PLACE,pair_count.data(),npairs,MPI_UNSIGNED_LONG_LONG,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,pair_dist_count_min.data(),npairs,MPI_UNSIGNED_LONG_LONG,MPI_MIN,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,pair_dist_count_max.data(),npairs,MPI_UNSIGNED_LONG_LONG,MPI_MAX,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,pair_sum.data(),npairs,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,pair_min.data(),npairs,MPI_DOUBLE,MPI_MIN,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,pair_max.data(),npairs,MPI_DOUBLE,MPI_MAX,*mpi);

      lout << "===== atom neighbors stats =====" << std::fixed << std::setprecision(2) << std::endl;
	    lout << "Chunk size             = "<<cs <<std::endl;
      lout << "Particles              = "<<total_particles<<std::endl;
      if(total_particles>0) lout << "Nbh cells (tot./avg)   = "<< total_nbh_cells <<" / "<< (total_nbh_cells*1.0/total_particles) <<std::endl;
      lout << "Neighbors (chunk/<d)   = "<<total_nbh <<" / "<<total_nbh_d2 << std::endl;
	    if(total_nbh>0) lout << "<d / chunk ratio       = " << (total_nbh_d2*100/total_nbh)*0.01 << " , storage eff. = "<< (total_nbh_d2*1.0/total_nbh)*cs <<std::endl;
      if(total_particles>0) lout << "Avg nbh (chunk/<d)     = "<< (total_nbh*1.0/total_particles) <<" / "<< (total_nbh_d2*1.0/total_particles) <<std::endl;
      lout << "min [chunk;<d] / Max [chunk;<d] = ["<< nbh_min<<";"<<nbh_d2_min <<"] / ["<< nbh_max <<";"<<nbh_d2_max<<"]"  <<std::endl;

      for(unsigned int ta=0;ta<ntypes;ta++)
      for(unsigned int tb=ta;tb<ntypes;tb++)
      {
        const unsigned int pair_id = unique_pair_id(ta,tb);
        if( pair_min[pair_id] != std::numeric_limits<double>::max() )
        {
          lout << species->at(ta).name() << " / " << species->at(tb).name()
               << " : count = " << large_integer_to_string(pair_count[pair_id])
               << " , min = " << pair_min[pair_id]
               << " , max = " << pair_max[pair_id]
               << " , count<"<<onika::format_string("% .3e",sqrt(pairdist[pair_id])) <<" = [ " << large_integer_to_string(pair_dist_count_min[pair_id])<<" ; "<< large_integer_to_string(pair_dist_count_max[pair_id])<<" ]"
               << " , avg = " << (pair_sum[pair_id] / pair_count[pair_id]) << std::endl;
        }
      }

      lout << "=================================" << std::defaultfloat << std::endl;
    }
  
  };

  template<class GridT> using AtomNeighborsStatsTmpl = AtomNeighborsStats<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(atom_neighbors_stats)
  {
   OperatorNodeFactory::instance()->register_factory("atom_neighbors_stats", make_grid_variant_operator< AtomNeighborsStatsTmpl > );
  }

}

