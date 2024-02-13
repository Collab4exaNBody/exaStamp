#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/domain.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/grid_algorithm.h>
#include <exanb/core/simple_block_rcb.h>

#include <exanb/core/grid_connected_components.h>

#include <unordered_map>
#include <unordered_set>
#include <vector> 
#include <utility>
#include <cstdint>
#include <mpi.h>
#include <omp.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <time.h>
//#include <sys/sysinfo.h>

namespace exaStamp
{
  using namespace exanb;

  class GridConnectedComponents : public OperatorNode
  {
    // using CellOwner  = std::vector<int>;
    // ADD_SLOT( CellOwner           , cell_owner               , INPUT , REQUIRED );
    
    ADD_SLOT( MPI_Comm              , mpi               , INPUT , MPI_COMM_WORLD , DocString{"MPI communicator"} );
    ADD_SLOT( Domain                , domain            , INPUT , REQUIRED , DocString{"The global domain"} );
    ADD_SLOT( GridCellValues        , grid_cell_values  , INPUT_OUTPUT , DocString{"the local grid"} );
    
    ADD_SLOT( double                , density_threshold , INPUT , 1. , DocString{"density threshold. below this threshold is considered empty"} );
    ADD_SLOT( bool                  , enable_debug      , INPUT , false );
    
    ADD_SLOT( ConnectedComponentMap , loc_label_map     , OUTPUT , DocString{"local map between the labels in the grid and their memory id"}); // needed by the analysis module.
    ADD_SLOT( ssize_t               , nb_cc             , OUTPUT , DocString{"Total number of CC"}); // for now it is only an output of the MPI root to avoid some computation on the analysis module.  

    ADD_SLOT( std::string           , filter            , INPUT_OUTPUT , "agregate" );
    
    ADD_SLOT( std::string    , contribution_field_name  , INPUT , REQUIRED );
    ADD_SLOT( std::string    , mass_field_name          , INPUT , REQUIRED );
    
    ADD_SLOT( std::string    , cc_label_field_name      , INPUT_OUTPUT , "cc_label" );
    
  public:

    // -----------------------------------------------
    // -----------------------------------------------

    inline void execute ()  override final
    {
      const ssize_t subdiv = grid_cell_values->field(*contribution_field_name).m_subdiv;

      // cell size (all cells are cubical)
      const double cell_size = domain->cell_size(); 
      // side size of a sub-cell
      const double subcell_size = cell_size / subdiv;
      // sub-cell's volume
      const double subcell_vol  = subcell_size * subcell_size * subcell_size;
      
      // -- Creation :
      // additional data field for connected component label
      if( ! grid_cell_values->has_field(*cc_label_field_name) )
        {
          grid_cell_values->add_field(*cc_label_field_name,subdiv,1);
        }
      // additional data field for thread id
      if( *enable_debug && ! grid_cell_values->has_field("thread_id") )
        {
          grid_cell_values->add_field("thread_id",subdiv,1);
        }
      
      // mass field data accessor.
      const auto mass_accessor = grid_cell_values->field_data(*mass_field_name);
      const double  * __restrict__ mass_ptr = mass_accessor.m_data_ptr;
      // const size_t mass_stride = mass_accessor.m_stride;

      // filter function, identifies "plain" subcells
      auto filter_func_agregate = [mass_ptr, threshold_value=*density_threshold * subcell_vol] ( ssize_t subcell_id ) -> bool { return (mass_ptr[subcell_id] > threshold_value); } ;
      
      // filter function, identifies "empty" subcells
      auto filter_func_void = [mass_ptr, threshold_value=*density_threshold * subcell_vol] ( ssize_t subcell_id ) -> bool { return (mass_ptr[subcell_id] < threshold_value); } ;

      if( (*filter) == "agregate" )
        {
          this->execute_with_filter_func( filter_func_agregate );
        }
      else if ( (*filter) == "void" )
        {
          this->execute_with_filter_func( filter_func_void );
        }
      else
        {
          std::cout << "\n\n WARNING !!!\n\n Unknown filter option for grid_connected_components module ! Available filters are \"agregate\" and \"void\"\n\n";
          abort();
        }
    }

    template<class FilterFuncT>
    inline void execute_with_filter_func( FilterFuncT filter_func )
    // inline void execute ()  override final
    {
      static_assert( sizeof(ssize_t) == sizeof(double) ); // the module works only in this case
      static_assert( sizeof(long long) == sizeof(ssize_t) ); // ensures we can send ssize_t with MPI_LONG_LONG
    
      MPI_Comm comm = *mpi;

      // Constants definition :

      // -- Domain :
      // simulation's grid dimensions
      const IJK domain_dims = domain->grid_dimension();
      // cell size (all cells are cubical)
      const double cell_size = domain->cell_size(); 
      
      // -- Grid : 
      // local processor's grid dimensions, including ghost cells
      const IJK grid_dims = grid_cell_values->grid_dims(); 
      // local processor's grid position in simulation's grid
      const IJK grid_offset = grid_cell_values->grid_offset();
      // ghost layers
      const ssize_t gl = grid_cell_values->ghost_layers();

      const ssize_t n_cells = grid_cell_values->number_of_cells();
      
      if( n_cells == 0 )
        {
          return;
        }
      
      // number of subdivisions, in each directions, applied on cells
      const ssize_t subdiv = grid_cell_values->field(*contribution_field_name).m_subdiv;
      // local processor's total grid dimension 
      // const IJK total_dims = grid_dims * subdiv;

      // -- Usefull constants construction :
      // local grid boundaries
      // const GridBlock local_space = {grid_offset, grid_offset + grid_dims};
      // number of sub-cells per cell
      const size_t n_subcells = subdiv*subdiv*subdiv;
      // dimension of the subdivided simulation's grid
      const IJK domain_subdiv_dims = domain_dims * subdiv; 
      // side size of a sub-cell
      const double subcell_size = cell_size / subdiv;
      // sub-cell's volume
      // const double subcell_vol  = subcell_size * subcell_size * subcell_size;

      // double threshold = *density_threshold * subcell_vol;
      
      // check that the integer maximum precision with double will be enough to avoid collision with the hashmap of our label.
      assert( ((unsigned long long int) 1 << 52) > (domain_subdiv_dims.i * domain_subdiv_dims.j * domain_subdiv_dims.k + 1));
      
      // -----------------------------------------------
      // some debug information
      ldbg << "ghost_layers="<<gl<<", cell_size="<<cell_size<<", subdiv="<<subdiv<<", subcell_size="
           <<subcell_size<<", grid_dims="<<grid_dims<<", grid_offset="<<grid_offset<<", domain_dims="<<domain_dims
           <<", domain_subdiv_dims="<<domain_subdiv_dims<<std::endl;
      // -----------------------------------------------
      
      // Grid data fields :
      // note: if we are to add new data fields, they must be added BEFORE we retreive access ANY information

      // // -- Creation :
      // // additional data field for connected component label
      // if( ! grid_cell_values->has_field(*cc_label_field_name) )
      //   {
      //     grid_cell_values->add_field(*cc_label_field_name,subdiv,1);
      //   }
      // // additional data field for thread id
      // if( *enable_debug && ! grid_cell_values->has_field("thread_id") )
      //   {
      //     grid_cell_values->add_field("thread_id",subdiv,1);
      //   }
      
      // -- Retrieve grid data field accessor :
      // cc_label field data accessor.
      const auto cc_label_accessor = grid_cell_values->field_data(*cc_label_field_name);
      // uint64_t * __restrict__ cc_label_ptr = (uint64_t *) cc_label_accessor.m_data_ptr;
      double * __restrict__ cc_label_ptr = cc_label_accessor.m_data_ptr;
      const size_t cc_label_stride = cc_label_accessor.m_stride;
      // set to 0
      for(ssize_t i=0;i<n_cells;i++)
        {
          for(int j=0;j<n_subcells;j++) { cc_label_ptr[ i * cc_label_stride + j ] = 0.0; }
        }
      
      // thread_id field data accessor.
      double * __restrict__ thread_id_ptr = nullptr;
      if (*enable_debug)
        {
          const auto thread_id_accessor = grid_cell_values->field_data("thread_id");
          thread_id_ptr = thread_id_accessor.m_data_ptr;
          // const size_t thread_id_stride = thread_id_accessor.m_stride;
        }

      // mass field data accessor.
      const auto mass_accessor = grid_cell_values->field_data(*mass_field_name);
      // const double  * __restrict__ mass_ptr = mass_accessor.m_data_ptr;
      const size_t mass_stride = mass_accessor.m_stride;

      // -- Sanity check :
      assert( cc_label_accessor.m_stride == mass_stride );
      const size_t stride = mass_stride; // for simplicity
      
      // ===============================================================================
      //timespec t_start, t_end, t_delta;
      //timespec t_start_in, t_end_in, t_delta_in;
      //clock_gettime(CLOCK_MONOTONIC, &t_start); 

      int rank=0, np = 0;
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &np);

      int OMP_NUM_THREADS = omp_get_max_threads();

      GridBlock work_block[OMP_NUM_THREADS];

      std::unordered_map<double, double> eq_table[OMP_NUM_THREADS];
      std::unordered_map<double, double> join_table; // used to connect in shared memory
      std::unordered_map<ssize_t, duet> interface_table; // used to connect in distributed memory, will be send to the root.
      interface_table.reserve((grid_dims.i * grid_dims.j * 2 + grid_dims.j * grid_dims.k * 2 + grid_dims.i * grid_dims.k * 2) ); // / 8 ?
      
      std::unordered_set<double> threads_root_set[OMP_NUM_THREADS]; // used to find the number of connected component in each thread space.
      std::unordered_set<double> domain_root_set; // used to find the number of connected component in the local domain.

      // std::unordered_map<double, ssize_t> loc_relab_map;
      std::unordered_map<double, ssize_t> global_relab_map; // used to relab all the labels

      std::unique_ptr<std::vector<ssize_t>> local_buf = std::make_unique<std::vector<ssize_t>>(); // new (ssize_t *) malloc(sizeof(ssize_t) * 3 * nb_triplet + nb_root);

      std::unordered_map<double, double> trad_map; // used to transmit the local relabelisation
      
      ssize_t *master_buf = nullptr;
      //size_t msg_size;
#     pragma omp parallel
      {
        int nthreads = omp_get_num_threads(); 
        int thread_id = omp_get_thread_num();
        
        assert(nthreads == OMP_NUM_THREADS); // check that there is the number of threads that we expect

        // Each thread retrieve its work block
        work_block[thread_id] = simple_block_rcb(GridBlock{IJK{gl,gl,gl}, grid_dims-gl}, nthreads, thread_id);
        
        eq_table[thread_id].reserve( GridBlock_size(work_block[thread_id]) * n_subcells / 8); // need to perform test on the parity of each direction. 
        
        // First pass : 
        for( ssize_t k=work_block[thread_id].start.k ; k < work_block[thread_id].end.k ; k++)
          for( ssize_t j=work_block[thread_id].start.j ; j < work_block[thread_id].end.j ; j++)
            for( ssize_t i=work_block[thread_id].start.i ; i < work_block[thread_id].end.i ; i++)
              {
                // informations that might be usefull later on ...
                // tells if cell is in the ghost layers ()
                // bool is_ghost = i<gl || i>=(grid_dims.i-gl) || j<gl || j>=(grid_dims.j-gl) || k<gl || k>=(grid_dims.k-gl);
                
                // position of the cell in the simulation grid, which size is 'domain_dims'
                //IJK cell_location = IJK{i,j,k} + grid_offset;

                // computation of cell index
                IJK cell_loc = {i,j,k};
                ssize_t cell_index = grid_ijk_to_index( grid_dims , cell_loc );
                
                // triple loop to enumerate sub cells inside a cell
                for( ssize_t sk=0 ; sk<subdiv ; sk++)
                  for( ssize_t sj=0 ; sj<subdiv ; sj++)
                    for( ssize_t si=0 ; si<subdiv ; si++)
                      {
                        // location of the subcell in simulation's subcell grid, which dimension is 'domain_subdiv_dims'
                        //IJK subcell_global_location = cell_location * subdiv + IJK{si,sj,sk};
                        //ssize_t subcell_global_index = grid_ijk_to_index( domain_subdiv_dims , subcell_global_location );

                        // computation of subcell index
                        IJK subcell_loc = {si,sj,sk};
                        ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , subcell_loc );
                        
                        // value index, is the index of current subcell for local processor's grid
                        ssize_t value_index = cell_index * stride + subcell_index;

                        if( filter_func( value_index ) )
                        // if(mass_ptr[ value_index ] < threshold)
                          {
                            find_par(domain_dims, grid_dims, grid_offset, work_block[thread_id], cell_loc, subcell_loc,
                                     subdiv, stride, cc_label_ptr, eq_table[thread_id]);
                          }
                        
                        if(*enable_debug)
                          thread_id_ptr[value_index] = rank*nthreads + thread_id+1; 
                      }
              }

        // One half pass :
        union_loc(eq_table[thread_id]);

        // join tables
#       pragma omp barrier
        if(nthreads > 1 && thread_id == 0)
          {
            size_t ncell_to_check = 0; // roughly
            for(int t = 0; t < nthreads; t++)
              {
                if (work_block[t].end.i < (grid_dims.i - gl))
                  ncell_to_check += (work_block[t].end.j - work_block[t].start.j) * (work_block[t].end.k - work_block[t].start.k);
                
                if (work_block[t].end.j < (grid_dims.j - gl))
                  ncell_to_check += (work_block[t].end.i - work_block[t].start.i) * (work_block[t].end.k - work_block[t].start.k);
                
                if (work_block[t].end.k < (grid_dims.k - gl))
                  ncell_to_check += (work_block[t].end.j - work_block[t].start.j) * (work_block[t].end.i - work_block[t].start.i);
              }
            
            join_table.reserve(ncell_to_check * subdiv * subdiv / 8);
              
            for(int t = 0; t < nthreads; t++)
              { 
                if (work_block[t].end.i < grid_dims.i - gl)
                  for(ssize_t k = work_block[t].start.k; k < work_block[t].end.k; k++)
                    for(ssize_t j = work_block[t].start.j; j < work_block[t].end.j; j++)
                      {
                        IJK cell_loc = {work_block[t].end.i-1,j,k};
                        ssize_t cell_index = grid_ijk_to_index( grid_dims , cell_loc );
                        for( ssize_t sk=0 ; sk<subdiv ; sk++)
                          for( ssize_t sj=0 ; sj<subdiv ; sj++)
                            {
                              IJK subcell_loc = IJK{subdiv-1,sj,sk};
                              ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , subcell_loc );
                              // value index, is the index of current subcell for local processor's grid
                              ssize_t value_index = cell_index * stride + subcell_index;
                              if(cc_label_ptr[ value_index ] > 0.0)
                                {
                                  double root_t = eq_table[t].at(cc_label_ptr[value_index]);
                                  for(int nk = -1; nk < 2; nk++)
                                    for(int nj = -1; nj < 2; nj++)
                                      {
                                        IJK nbh_cell_loc;
                                        IJK nbh_subcell_loc;
                                        subcell_neighbor( cell_loc, subcell_loc, subdiv, IJK{1,nj,nk}, nbh_cell_loc, nbh_subcell_loc );
                                        size_t nbh_id = grid_ijk_to_index(grid_dims, nbh_cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                                        if(is_inside(GridBlock{IJK{gl,gl,gl}, grid_dims-gl}, nbh_cell_loc) && cc_label_ptr[ nbh_id ] > 0.0)
                                          {
                                            int nbh_t = 0;
                                            while (!is_inside(work_block[nbh_t], nbh_cell_loc))
                                              nbh_t++;
                                            
                                            double root_nbh = eq_table[nbh_t].at(cc_label_ptr[nbh_id]);
                                            
                                            if (root_nbh < root_t)
                                              {
                                                double tmp = root_t;
                                                root_t = root_nbh;
                                                root_nbh = tmp;
                                              }
                                            
                                            insert(join_table, root_t, root_nbh);
                                          }
                                      }   
                                }
                            }
                      }
                
                if (work_block[t].end.j < grid_dims.j - gl)
                  for(ssize_t k = work_block[t].start.k; k < work_block[t].end.k; k++)
                    for(ssize_t i = work_block[t].start.i; i < work_block[t].end.i; i++)
                      {
                        IJK cell_loc = {i,work_block[t].end.j-1,k};
                        ssize_t cell_index = grid_ijk_to_index( grid_dims , cell_loc );
                        for( ssize_t sk=0 ; sk<subdiv ; sk++)
                          for( ssize_t si=0 ; si<subdiv ; si++)
                            {
                              IJK subcell_loc = IJK{si,subdiv-1,sk};
                              ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , subcell_loc );
                              // value index, is the index of current subcell for local processor's grid
                              ssize_t value_index = cell_index * stride + subcell_index;
                              if(cc_label_ptr[ value_index ] > 0.0)
                                {
                                  double root_t = eq_table[t].at(cc_label_ptr[value_index]);
                                  for(int nk = -1; nk < 2; nk++)
                                    for(int ni = -1; ni < 2; ni++)
                                      {
                                        IJK nbh_cell_loc;
                                        IJK nbh_subcell_loc;
                                        subcell_neighbor( cell_loc, subcell_loc, subdiv, IJK{ni,1,nk}, nbh_cell_loc, nbh_subcell_loc );
                                        size_t nbh_id = grid_ijk_to_index(grid_dims, nbh_cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                                        if(is_inside(GridBlock{IJK{gl,gl,gl}, grid_dims-gl}, nbh_cell_loc) && cc_label_ptr[ nbh_id ] > 0.0)
                                          {
                                            int nbh_t = 0;
                                            while (!is_inside(work_block[nbh_t], nbh_cell_loc))
                                              nbh_t++;
                                              
                                            double root_nbh = eq_table[nbh_t].at(cc_label_ptr[nbh_id]);
                                            
                                            if (root_nbh < root_t)
                                              {
                                                double tmp = root_t;
                                                root_t = root_nbh;
                                                root_nbh = tmp;
                                              }
                                            
                                            insert(join_table, root_t, root_nbh);
                                          }
                                      } 
                                }
                            }
                      }
                
                if (work_block[t].end.k < grid_dims.k - gl)
                  for(ssize_t j = work_block[t].start.j; j < work_block[t].end.j; j++)
                    for(ssize_t i = work_block[t].start.i; i < work_block[t].end.i; i++)
                      {
                        IJK cell_loc = {i,j,work_block[t].end.k-1};
                        ssize_t cell_index = grid_ijk_to_index( grid_dims , cell_loc );
                        for( ssize_t sj=0 ; sj<subdiv ; sj++)
                          for( ssize_t si=0 ; si<subdiv ; si++)
                            {
                              IJK subcell_loc = IJK{si,sj,subdiv-1};
                              ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , subcell_loc );
                              // value index, is the index of current subcell for local processor's grid
                              ssize_t value_index = cell_index * stride + subcell_index;
                              if(cc_label_ptr[ value_index ] > 0.0)
                                {
                                  double root_t = eq_table[t].at(cc_label_ptr[value_index]);
                                  for(int nj = -1; nj < 2; nj++)
                                    for(int ni = -1; ni < 2; ni++)
                                      {
                                        IJK nbh_cell_loc;
                                        IJK nbh_subcell_loc;
                                        subcell_neighbor( cell_loc, subcell_loc, subdiv, IJK{ni,nj,1}, nbh_cell_loc, nbh_subcell_loc );
                                        size_t nbh_id = grid_ijk_to_index(grid_dims, nbh_cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                                        if(is_inside(GridBlock{IJK{gl,gl,gl}, grid_dims-gl}, nbh_cell_loc) && cc_label_ptr[ nbh_id ] > 0.0)
                                          {
                                            int nbh_t = 0;
                                            while (!is_inside(work_block[nbh_t], nbh_cell_loc))
                                              nbh_t++;
                                            
                                            double root_nbh = eq_table[nbh_t].at(cc_label_ptr[nbh_id]);
                                            
                                            if (root_nbh < root_t)
                                              {
                                                double tmp = root_t;
                                                root_t = root_nbh;
                                                root_nbh = tmp;
                                              }
                                            
                                            insert(join_table, root_t, root_nbh);
                                          }
                                      } 
                                }
                            }
                      }
                union_loc(join_table);
              }     
          }
#       pragma omp barrier
        
        join_and_extract_local_thread_root(join_table, eq_table[thread_id], threads_root_set[thread_id]);

#       pragma omp barrier

#       pragma omp single        
        if(nthreads > 1)
          {
            size_t sum = 0;
            for(int t = 0; t < nthreads; t++)
              sum += threads_root_set[t].size();
            
            domain_root_set.reserve(sum);
            for(int t = 0; t < nthreads; t++)
              {
                domain_root_set.insert(threads_root_set[t].begin(), threads_root_set[t].end());
              }
          }
        else // should not happen
          {
            domain_root_set = threads_root_set[0];
          }
        

        if(np > 1 && thread_id == 0)
          {   
            //clock_gettime(CLOCK_MONOTONIC, &t_start_in);

            for(int t = 0; t < nthreads; t++)
              (const std::unordered_map<double, double>) eq_table[t];
            
            // check "i-faces" of the cuboid 
            for( ssize_t k = gl ; k < grid_dims.k-gl ; k++)
              for( ssize_t j = gl ; j < grid_dims.j-gl ; j++)
                {             
                  for( ssize_t sk=0 ; sk<subdiv ; sk++)
                    for( ssize_t sj=0 ; sj<subdiv ; sj++)
                      {
                        IJK cell_loc = {gl,j,k};
                        IJK subcell_loc = IJK{0,sj,sk};
                        ssize_t cell_index = grid_ijk_to_index( grid_dims , cell_loc );
                        ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , subcell_loc );
                        // value index, is the index of current subcell for local processor's grid
                        ssize_t value_index = cell_index * stride + subcell_index;
                        if(cc_label_ptr[ value_index ] > 0.0)
                          {
                            for(int nk = -1; nk < 2; nk++)
                              for(int nj = -1; nj < 2; nj++)
                                {
                                  IJK nbh_cell_loc;
                                  IJK nbh_subcell_loc;
                                  subcell_neighbor( cell_loc, subcell_loc, subdiv, IJK{-1,nj,nk}, nbh_cell_loc, nbh_subcell_loc );
                                  size_t nbh_id = grid_ijk_to_index(grid_dims, nbh_cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                                  if(filter_func(nbh_id))
                                    {
                                      // compute global location of the sub-cells
                                      ssize_t subcell_global_index = grid_ijk_to_index(domain_subdiv_dims, (cell_loc + grid_offset) * subdiv + subcell_loc);
                                      ssize_t nbh_subcell_global_index = grid_ijk_to_index(domain_subdiv_dims, (nbh_cell_loc + grid_offset) * subdiv + nbh_subcell_loc);
                                        
                                      // find in which table is store the root label of the subcell
                                      int t = 0;
                                      while ( !is_inside(work_block[t], cell_loc) )
                                        t++;
                                        
                                      // write the triplet in the interface_table
                                      interface_table.try_emplace(subcell_global_index, duet{nbh_subcell_global_index, eq_table[t].at(cc_label_ptr[value_index])});
                                    }
                                }
                          }
                      }
                }
            for( ssize_t k = gl ; k < grid_dims.k-gl ; k++)
              for( ssize_t j = gl ; j < grid_dims.j-gl ; j++)
                {             
                  for( ssize_t sk=0 ; sk<subdiv ; sk++)
                    for( ssize_t sj=0 ; sj<subdiv ; sj++)
                      {
                        IJK cell_loc = {grid_dims.i-gl-1,j,k};
                        IJK subcell_loc = IJK{subdiv-1,sj,sk};
                        ssize_t cell_index = grid_ijk_to_index( grid_dims , cell_loc );
                        ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , subcell_loc );
                        // value index, is the index of current subcell for local processor's grid
                        ssize_t value_index = cell_index * stride + subcell_index;
                        if(cc_label_ptr[ value_index ] > 0.0)
                          {
                            for(int nk = -1; nk < 2; nk++)
                              for(int nj = -1; nj < 2; nj++)
                                {
                                  IJK nbh_cell_loc;
                                  IJK nbh_subcell_loc;
                                  subcell_neighbor( cell_loc, subcell_loc, subdiv, IJK{1,nj,nk}, nbh_cell_loc, nbh_subcell_loc );
                                  size_t nbh_id = grid_ijk_to_index(grid_dims, nbh_cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                                  if(filter_func(nbh_id))
                                    {
                                      // compute global location of the sub-cells
                                      ssize_t subcell_global_index = grid_ijk_to_index(domain_subdiv_dims, (cell_loc + grid_offset) * subdiv + subcell_loc);
                                      ssize_t nbh_subcell_global_index = grid_ijk_to_index(domain_subdiv_dims, (nbh_cell_loc + grid_offset) * subdiv + nbh_subcell_loc);
                                        
                                      // find in which table is store the root label of the subcell
                                      int t = 0;
                                      while ( !is_inside(work_block[t], cell_loc) )
                                        t++;
                                        
                                      // write the triplet in the interface_table
                                      interface_table.try_emplace(subcell_global_index, duet{nbh_subcell_global_index, eq_table[t].at(cc_label_ptr[value_index])});
                                    }
                                }
                          }
                      }
                }

            // check "j-faces" of the cuboid 
            for( ssize_t k = gl ; k < grid_dims.k-gl ; k++)
              for( ssize_t i = gl ; i < grid_dims.i-gl ; i++)
                {             
                  for( ssize_t sk=0 ; sk<subdiv ; sk++)
                    for( ssize_t si=0 ; si<subdiv ; si++)
                      {
                        IJK cell_loc = {i,gl,k};
                        IJK subcell_loc = IJK{si,0,sk};
                        ssize_t cell_index = grid_ijk_to_index( grid_dims , cell_loc );
                        ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , subcell_loc );
                        // value index, is the index of current subcell for local processor's grid
                        ssize_t value_index = cell_index * stride + subcell_index;
                        if(cc_label_ptr[ value_index ] > 0.0)
                          {
                            for(int nk = -1; nk < 2; nk++)
                              for(int ni = -1; ni < 2; ni++)
                                {
                                  IJK nbh_cell_loc;
                                  IJK nbh_subcell_loc;
                                  subcell_neighbor( cell_loc, subcell_loc, subdiv, IJK{ni,-1,nk}, nbh_cell_loc, nbh_subcell_loc );
                                  size_t nbh_id = grid_ijk_to_index(grid_dims, nbh_cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                                  if(filter_func(nbh_id))
                                    {
                                      // compute global location of the sub-cells
                                      ssize_t subcell_global_index = grid_ijk_to_index(domain_subdiv_dims, (cell_loc + grid_offset) * subdiv + subcell_loc);
                                      ssize_t nbh_subcell_global_index = grid_ijk_to_index(domain_subdiv_dims, (nbh_cell_loc + grid_offset) * subdiv + nbh_subcell_loc);
                                        
                                      // find in which table is store the root label of the subcell
                                      int t = 0;
                                      while ( !is_inside(work_block[t], cell_loc) )
                                        t++;
                                        
                                      // write the triplet in the interface_table
                                      interface_table.try_emplace(subcell_global_index, duet{nbh_subcell_global_index, eq_table[t].at(cc_label_ptr[value_index])});
                                    }
                                }
                          }
                      }
                }
            for( ssize_t k = gl ; k < grid_dims.k-gl ; k++)
              for( ssize_t i = gl ; i < grid_dims.i-gl ; i++)
                {             
                  for( ssize_t sk=0 ; sk<subdiv ; sk++)
                    for( ssize_t si=0 ; si<subdiv ; si++)
                      {
                        IJK cell_loc = {i,grid_dims.j-gl-1,k};
                        IJK subcell_loc = IJK{si,subdiv-1,sk};
                        ssize_t cell_index = grid_ijk_to_index( grid_dims , cell_loc );
                        ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , subcell_loc );
                        // value index, is the index of current subcell for local processor's grid
                        ssize_t value_index = cell_index * stride + subcell_index;
                        if(cc_label_ptr[ value_index ] > 0.0)
                          {
                            for(int nk = -1; nk < 2; nk++)
                              for(int ni = -1; ni < 2; ni++)
                                {
                                  IJK nbh_cell_loc;
                                  IJK nbh_subcell_loc;
                                  subcell_neighbor( cell_loc, subcell_loc, subdiv, IJK{ni,1,nk}, nbh_cell_loc, nbh_subcell_loc );
                                  size_t nbh_id = grid_ijk_to_index(grid_dims, nbh_cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                                  if(filter_func(nbh_id))
                                    {
                                      // compute global location of the sub-cells
                                      ssize_t subcell_global_index = grid_ijk_to_index(domain_subdiv_dims, (cell_loc + grid_offset) * subdiv + subcell_loc);
                                      ssize_t nbh_subcell_global_index = grid_ijk_to_index(domain_subdiv_dims, (nbh_cell_loc + grid_offset) * subdiv + nbh_subcell_loc);
                                        
                                      // find in which table is store the root label of the subcell
                                      int t = 0;
                                      while ( !is_inside(work_block[t], cell_loc) )
                                        t++;
                                        
                                      // write the triplet in the interface_table
                                      interface_table.try_emplace(subcell_global_index, duet{nbh_subcell_global_index, eq_table[t].at(cc_label_ptr[value_index])});
                                    }
                                }
                          }
                      }
                }

            // check "k-faces" of the cuboid 
            for( ssize_t j = gl ; j < grid_dims.j-gl ; j++)
              for( ssize_t i = gl ; i < grid_dims.i-gl ; i++)
                {             
                  for( ssize_t sj=0 ; sj<subdiv ; sj++)
                    for( ssize_t si=0 ; si<subdiv ; si++)
                      {
                        IJK cell_loc = {i,j,gl};
                        IJK subcell_loc = IJK{si,sj,0};
                        ssize_t cell_index = grid_ijk_to_index( grid_dims , cell_loc );
                        ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , subcell_loc );
                        // value index, is the index of current subcell for local processor's grid
                        ssize_t value_index = cell_index * stride + subcell_index;
                        if(cc_label_ptr[ value_index ] > 0.0)
                          {
                            for(int nj = -1; nj < 2; nj++)
                              for(int ni = -1; ni < 2; ni++)
                                {
                                  IJK nbh_cell_loc;
                                  IJK nbh_subcell_loc;
                                  subcell_neighbor( cell_loc, subcell_loc, subdiv, IJK{ni,nj,-1}, nbh_cell_loc, nbh_subcell_loc );
                                  size_t nbh_id = grid_ijk_to_index(grid_dims, nbh_cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                                  if(filter_func(nbh_id))
                                    {
                                      // compute global location of the sub-cells
                                      ssize_t subcell_global_index = grid_ijk_to_index(domain_subdiv_dims, (cell_loc + grid_offset) * subdiv + subcell_loc);
                                      ssize_t nbh_subcell_global_index = grid_ijk_to_index(domain_subdiv_dims, (nbh_cell_loc + grid_offset) * subdiv + nbh_subcell_loc);
                                        
                                      // find in which table is store the root label of the subcell
                                      int t = 0;
                                      while ( !is_inside(work_block[t], cell_loc) )
                                        t++;
                                        
                                      // write the triplet in the interface_table
                                      interface_table.try_emplace(subcell_global_index, duet{nbh_subcell_global_index, eq_table[t].at(cc_label_ptr[value_index])});
                                    }
                                }
                          }
                      }
                }
            for( ssize_t j = gl ; j < grid_dims.j-gl ; j++)
              for( ssize_t i = gl ; i < grid_dims.i-gl ; i++)
                {             
                  for( ssize_t sj=0 ; sj<subdiv ; sj++)
                    for( ssize_t si=0 ; si<subdiv ; si++)
                      {
                        IJK cell_loc = {i,j,grid_dims.k-gl-1};
                        IJK subcell_loc = IJK{si,sj,subdiv-1};
                        ssize_t cell_index = grid_ijk_to_index( grid_dims , cell_loc );
                        ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , subcell_loc );
                        // value index, is the index of current subcell for local processor's grid
                        ssize_t value_index = cell_index * stride + subcell_index;
                        if(cc_label_ptr[ value_index ] > 0.0)
                          {
                            for(int nj = -1; nj < 2; nj++)
                              for(int ni = -1; ni < 2; ni++)
                                {
                                  IJK nbh_cell_loc;
                                  IJK nbh_subcell_loc;
                                  subcell_neighbor( cell_loc, subcell_loc, subdiv, IJK{ni,nj,1}, nbh_cell_loc, nbh_subcell_loc );
                                  size_t nbh_id = grid_ijk_to_index(grid_dims, nbh_cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                                  if(filter_func(nbh_id))
                                    {
                                      // compute global location of the sub-cells
                                      ssize_t subcell_global_index = grid_ijk_to_index(domain_subdiv_dims, (cell_loc + grid_offset) * subdiv + subcell_loc);
                                      ssize_t nbh_subcell_global_index = grid_ijk_to_index(domain_subdiv_dims, (nbh_cell_loc + grid_offset) * subdiv + nbh_subcell_loc);
                                        
                                      // find in which table is store the root label of the subcell
                                      int t = 0;
                                      while ( !is_inside(work_block[t], cell_loc) )
                                        t++;
                                        
                                      // write the triplet in the interface_table
                                      interface_table.try_emplace(subcell_global_index, duet{nbh_subcell_global_index, eq_table[t].at(cc_label_ptr[value_index])});
                                    }
                                }
                          }
                      }
                }
            
            int nb_triplet = interface_table.size();
            //int *nb_triplet_buf = nullptr;
            
            int nb_root = domain_root_set.size();
            int message_elements[2] = {nb_triplet, nb_root};
            int *message_buf = nullptr;

            //ssize_t *local_buf = nullptr;
            // std::unique_ptr<ssize_t[]> local_buf = std::make_unique<ssize_t[]>( 3 * nb_triplet); // new (ssize_t *) malloc(sizeof(ssize_t) * 3 * nb_triplet);
            // std::unique_ptr<ssize_t[]> local_buf = std::make_unique<ssize_t[]>( 3 * nb_triplet + nb_root); // new (ssize_t *) malloc(sizeof(ssize_t) * 3 * nb_triplet);

            (*local_buf).reserve(3 * nb_triplet + nb_root);
            
            int *displs = nullptr;
            int *vcounts = nullptr;
            MPI_Request requests[2];
            MPI_Status status[2];
            
            if (rank != 0)
              {
                MPI_Igather(&message_elements, 2, MPI_INT, message_buf, 2, MPI_INT, 0, comm, &requests[0]);
                
                // Write data in the sendbuf that is local_buf
                for(auto& cur:interface_table) // We use push_back instead of emplace_back to avoid unnecessary work from the compilator, if it produce warning we will have to use direct memory access
                  {
                    (*local_buf).push_back(cur.first);
                    (*local_buf).push_back(cur.second.first);
                    (*local_buf).push_back(convert(cur.second.second)); // use convert to serialize different type in MPI buffer and avoid a compilator warning made by *(double*)&  
                  }
                for(auto& cur:domain_root_set)
                  {
                    (*local_buf).push_back(convert(cur)); 
                  }

                MPI_Wait(&requests[0], &status[0]);

                MPI_Gatherv((*local_buf).data(), 3 * nb_triplet + nb_root, MPI_LONG_LONG, master_buf, vcounts, displs, MPI_LONG_LONG, 0, comm);

                // re-use of the send buffer as a recive buffer
                (*local_buf).resize(nb_root);
                (*local_buf).shrink_to_fit();
                
                MPI_Scatterv(master_buf, vcounts, displs, MPI_LONG_LONG, (*local_buf).data(), nb_root, MPI_LONG_LONG, 0, comm);
              }
            else // the root
              {
                message_buf = new int[np*2];
                MPI_Igather(&message_elements, 2, MPI_INT, message_buf, 2, MPI_INT, 0, comm, &requests[0]);
                
                // Write data in the sendbuf that is local_buf
                for(auto& cur:interface_table) // We use push_back instead of emplace_back to avoid unnecessary work from the compilator, if it produce warning we will have to use direct memory access
                  {
                    (*local_buf).push_back(cur.first);
                    (*local_buf).push_back(cur.second.first);
                    (*local_buf).push_back(convert(cur.second.second)); // use convert to serialize different type in MPI buffer and avoid a compilator warning made by *(double*)&  
                  }
                for(auto& cur:domain_root_set)
                  {
                    (*local_buf).push_back(convert(cur)); 
                  }
                
                // displs = (int *) malloc(sizeof(int) * np);
                displs = new int[np];
                // vcounts = (int *) malloc(sizeof(int) * np);
                vcounts = new int[np];
                
                MPI_Wait(&requests[0], &status[0]);
                
                displs[0] = 0;
                vcounts[0] = 3 * message_buf[0] + message_buf[1];
                for(int i = 1; i < np; i++)
                  {
                    displs[i] = displs[i-1] + vcounts[i-1];
                    vcounts[i] = 3 * message_buf[2*i] + message_buf[2*i+1];
                  }
                
                // master_buf = (ssize_t *) malloc(sizeof(ssize_t) * 3 * cnt);
                master_buf = new ssize_t[displs[np-1] + vcounts[np-1]];
                
                MPI_Gatherv((*local_buf).data(), 3 * nb_triplet + nb_root, MPI_LONG_LONG, master_buf, vcounts, displs, MPI_LONG_LONG, 0, comm);
                
                // local_buf.reset(nullptr); // free(local_buf);
                
                std::unordered_map<double, double> link_table;
                std::unordered_map<ssize_t, double> linker;
                
                for(int i = 0; i < np; i++)
                  {
                    for(int j = 0; j < (vcounts[i] - message_buf[(2*i) + 1]); j+=3) // vcounts[i] - message_buf[(2*i) + 1] is equivalent to 3 * message_buf[2*i]
                      {
                        auto search = linker.find(master_buf[displs[i] + j + 1]);
                        if(search == linker.end())
                          {
                            if(linker.try_emplace(master_buf[displs[i] + j], *(double*)& master_buf[displs[i] + j + 2]).second == true )
                              link_table.try_emplace(*(double*)& master_buf[displs[i] + j + 2], *(double*)& master_buf[displs[i] + j + 2]);
                          }
                        else
                          {
                            double root2 = search->second;
                            while(root2 != link_table[root2])
                              root2 = link_table[root2];

                            double root1 = *(double*)& master_buf[displs[i] + j + 2];
                            if(link_table.find(root1) == link_table.end())
                              {
                                if(root1 < root2)
                                  {
                                    link_table[root1] = root1;
                                    link_table[root2] = root1;
                                  }
                                else
                                  {
                                    link_table[root1] = root2;
                                  }
                              }
                            else
                              {
                                while(root1 != link_table[root1])
                                  root1 = link_table[root1];
                                
                                if(root1 < root2)
                                  link_table[root2] = root1;
                                else
                                  link_table[root1] = root2;
                              }
                          }
                      }
                  }
                union_loc(link_table);

                {
                  ssize_t label_id = 1;
                  for(int i = 0; i < np; i++)
                    {
                      for(int j = 0; j < message_buf[(2*i) + 1]; j++)
                        {
                          double label = *(double*)& master_buf[displs[i] + 3 * message_buf[2*i] + j];
                          auto test =  link_table.find(label);
                          if( (test == link_table.end() || test->first == test->second) ) // true mean it is a root (allow to pop out the final table leaves)
                            if( global_relab_map.try_emplace(label, label_id).second == true)
                              {
                                label_id++;
                              } 
                        }
                    }
                }

                // Set the number of connected components in the domain for the root to avoid some computations in the analysis module (and in the future for a better allocation in this module)   
                (*nb_cc) = global_relab_map.size();

                free_container(linker);

                // adapt infos on the message to send. We re-use the master buffer but only the parts with the domains roots sets of each procs.
                for(int i = 0; i < np; i++ )
                  {
                    displs[i] += 3 * message_buf[2*i];
                    vcounts[i] = message_buf[2*i + 1];
                  }

                // adapt our local_buf 
                (*local_buf).resize(nb_root);
                (*local_buf).shrink_to_fit();

                for(int i = 0; i < np; i++)
                  {
                    for(int j = 0; j < vcounts[i]; j++)
                      {
                        double cur_label = *(double *)& master_buf[displs[i] + j];
                        master_buf[displs[i] + j] = ( link_table.find(cur_label) == link_table.end() ) ? global_relab_map[cur_label] : global_relab_map[link_table[cur_label]];
                      }
                  }
                
                MPI_Scatterv(master_buf, vcounts, displs, MPI_LONG_LONG, (*local_buf).data(), nb_root, MPI_LONG_LONG, 0, comm);
              }
            //clock_gettime(CLOCK_MONOTONIC, &t_end_in);
          }

#       pragma omp barrier
        
#       pragma omp single nowait
        {
          free_container(join_table);
        }

#       pragma omp single
        {
          ConnectedComponentMap new_map;
          *loc_label_map = new_map;
          trad_map.reserve(domain_root_set.size());
          if(np > 1)
            {
              int i = 0;
              ssize_t j = 0;
              for(auto &cur:domain_root_set)
                {
                  trad_map[cur] = (double) (*local_buf)[i];
                  i++;
                  
                  if((*loc_label_map).try_emplace(trad_map[cur] ,j).second == true)
                    {
                      j++;
                    }
                }
              // std::cout << "loc_label_map out of " << rank <<": " << (*loc_label_map).size() << "\n";
            }
          else
            {
              double i = 1.0;
              ssize_t j = 0;
              for(auto &cur:domain_root_set)
                {
                  trad_map[cur] = i;
                  i += 1.0;

                  (*loc_label_map)[cur] = j;
                  j++;
                }
            }
        } // implicit barrier that is needed  
        
        // Second pass :
        for( ssize_t k=work_block[thread_id].start.k ; k < work_block[thread_id].end.k ; k++)
          for( ssize_t j=work_block[thread_id].start.j ; j < work_block[thread_id].end.j ; j++)
            for( ssize_t i=work_block[thread_id].start.i ; i < work_block[thread_id].end.i ; i++)
              {
                // IJK cell_location = IJK{i,j,k} + grid_offset;

                for( ssize_t sk=0 ; sk<subdiv ; sk++)
                  for( ssize_t sj=0 ; sj<subdiv ; sj++)
                    for( ssize_t si=0 ; si<subdiv ; si++)
                      {
                        ssize_t cell_index = grid_ijk_to_index( grid_dims , IJK{i,j,k} );
                        ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , IJK{si,sj,sk} );

                        ssize_t value_index = cell_index * stride + subcell_index;
                      
                        if(cc_label_ptr[ value_index ] > 0.0)
                          {
                            cc_label_ptr[value_index] = trad_map[eq_table[thread_id][cc_label_ptr[value_index]]];
                          }
                      }
              }
      }// end parallel region
      
      //clock_gettime(CLOCK_MONOTONIC, &t_end);
      /* 
         t_delta.tv_nsec = t_end.tv_nsec - t_start.tv_nsec;
         t_delta.tv_sec  = t_end.tv_sec - t_start.tv_sec;
         if (t_delta.tv_sec > 0 && t_delta.tv_nsec < 0)
         {
         t_delta.tv_nsec += 1000000000;
         t_delta.tv_sec--;
         }

         t_delta_in.tv_nsec = t_end_in.tv_nsec - t_start_in.tv_nsec;
         t_delta_in.tv_sec  = t_end_in.tv_sec - t_start_in.tv_sec;
         if (t_delta_in.tv_sec > 0 && t_delta_in.tv_nsec < 0)
         {
         t_delta_in.tv_nsec += 1000000000;
         t_delta_in.tv_sec--;
         }
      */
      //fprintf(stdout, "\nAlgorithm execution time : %d.%.9ld\t pid : %d\n", (int)t_delta.tv_sec, t_delta.tv_nsec, rank);
      //fprintf(stdout, "\nAlgorithm execution time in distributed code: %d.%.9ld\t pid : %d\n", (int)t_delta_in.tv_sec, t_delta_in.tv_nsec, rank);
    }

    // inline ssize_t subcell_value_index(const IJK& total_dims, ssize_t subdiv, size_t stride, const IJK& coord)
    // {
    //   ssize_t si = coord.i % subdiv;
    //   ssize_t sj = coord.j % subdiv;
    //   ssize_t sk = coord.k % subdiv;
    //   ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , IJK{si,sj,sk} );
    //   IJK coarse_coord = coord / subdiv;
    //   IJK coarse_dims = total_dims / subdiv;
    //   ssize_t coarse_cell_index = grid_ijk_to_index( coarse_dims , coarse_coord );
    //   return coarse_cell_index * stride + subcell_index;
    // }

    static inline size_t GridBlock_size(GridBlock &block)
    {
      return (block.end.i-block.start.i) * (block.end.j-block.start.j) * (block.end.k-block.start.k);
    }

    static inline void subcell_neighbor( const IJK& cell_loc, const IJK& subcell_loc, ssize_t subdiv, IJK ninc, IJK& nbh_cell_loc, IJK& nbh_subcell_loc )
    {
      nbh_cell_loc = cell_loc;
      nbh_subcell_loc = subcell_loc + ninc;
      if(nbh_subcell_loc.i<0) { -- nbh_cell_loc.i; } else if(nbh_subcell_loc.i>=subdiv) { ++ nbh_cell_loc.i; }
      if(nbh_subcell_loc.j<0) { -- nbh_cell_loc.j; } else if(nbh_subcell_loc.j>=subdiv) { ++ nbh_cell_loc.j; }
      if(nbh_subcell_loc.k<0) { -- nbh_cell_loc.k; } else if(nbh_subcell_loc.k>=subdiv) { ++ nbh_cell_loc.k; }
      nbh_subcell_loc.i = ( nbh_subcell_loc.i + subdiv ) % subdiv;
      nbh_subcell_loc.j = ( nbh_subcell_loc.j + subdiv ) % subdiv;
      nbh_subcell_loc.k = ( nbh_subcell_loc.k + subdiv ) % subdiv;      
    }

    static inline bool is_inside(const GridBlock& block, const IJK& cell_loc)
    {
      return (cell_loc.i >= block.start.i) && (cell_loc.j >= block.start.j) && (cell_loc.k >= block.start.k) &&
        (cell_loc.i < block.end.i) && (cell_loc.j < block.end.j) && (cell_loc.k < block.end.k);
    }
    
    static inline void find_par(const IJK& domain_dims,
                                const IJK& grid_dims,
                                const IJK& offset,
                                const GridBlock& work_block,
                                const IJK& cell_loc, 
                                const IJK& subcell_loc, 
                                const ssize_t subdiv, 
                                const size_t stride,
                                double * __restrict__ cc_label_ptr,
                                std::unordered_map<double, double>& eq_table
                                )
    {
      size_t id = grid_ijk_to_index(grid_dims, cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, subcell_loc);
      size_t label_boundary = grid_ijk_to_index(domain_dims * subdiv, domain_dims * subdiv) + 2;
      double prov = label_boundary;
      for( ssize_t nk=-1 ; nk < 2 ; nk++)
        for( ssize_t nj=-1 ; nj < 2; nj++)
          for( ssize_t ni=-1 ; ni < 2; ((nk==0 && nj==0)? ni+=2 :ni++))
            {
              IJK nbh_cell_loc;
              IJK nbh_subcell_loc;
              subcell_neighbor( cell_loc, subcell_loc, subdiv, IJK{ni,nj,nk}, nbh_cell_loc, nbh_subcell_loc );
              size_t nbh_id = grid_ijk_to_index(grid_dims, nbh_cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
              if(is_inside(work_block, nbh_cell_loc) && cc_label_ptr[nbh_id] > 0.0)
                {
                  double cursor = cc_label_ptr[nbh_id];
                  while(eq_table[cursor] != cursor)
                    {
                      cursor = eq_table[cursor];
                    }

                  prov = (cursor < prov)? cursor : prov;
                }
            }
      
      if (prov == label_boundary)
        {
          cc_label_ptr[id] = (double) grid_ijk_to_index(domain_dims * subdiv, (cell_loc + offset) * subdiv + subcell_loc) + 1.0;
          eq_table[cc_label_ptr[id]] = cc_label_ptr[id];
        }
      else
        {
          cc_label_ptr[id] = prov;
          
          for( ssize_t nk=-1 ; nk < 2 ; nk++)
            for( ssize_t nj=-1 ; nj < 2; nj++)
              for( ssize_t ni=-1 ; ni < 2; ((nk==0 && nj==0)? ni+=2 :ni++))
                {
                  IJK nbh_cell_loc;
                  IJK nbh_subcell_loc;
                  subcell_neighbor( cell_loc, subcell_loc, subdiv, IJK{ni,nj,nk}, nbh_cell_loc, nbh_subcell_loc );
                  size_t nbh_id = grid_ijk_to_index(grid_dims, nbh_cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
                  if(is_inside(work_block, nbh_cell_loc) && cc_label_ptr[nbh_id] > 0.0)
                    {
                      double cursor = cc_label_ptr[nbh_id];
                      while(eq_table[cursor] != cursor)
                        {
                          cursor = eq_table[cursor];
                        }
                      eq_table[cursor] = prov;
                    }
                }
        }
    }

    static inline void union_loc(std::unordered_map<double, double>& eq_table)
    {
      for(auto& cur : eq_table)
        {
          double next = eq_table[cur.first]; // <=> next = cur.second
          while(eq_table[next] != next)
            {
              next = eq_table[next];
            }
          eq_table[cur.first] = next; // <=> cur.second = next
        }
    }
    
    template<typename T=double>
    static inline void union_loc_and_root_find(std::unordered_map<T, T>& eq_table, std::unordered_set<T>& root_set)
    {
      for(auto& cur : eq_table)
        {
          if(cur.first == cur.second)
            {
              root_set.insert(cur.first);
            }
          else
            {
              T next = eq_table[cur.first]; // <=> next = cur.second
              while(eq_table[next] != next)
                {
                  next = eq_table[next];
                }
              eq_table[cur.first] = next; // <=> cur.second = next
            }
        }
    }

    static inline void insert(std::unordered_map<double,double>& table, double& a, double& b)
    {
      double root1;
      if(table.try_emplace(a, a).second == false)
        {
          root1 = table[a];
          while(table[root1] != root1)
            root1 = table[root1];
        }
      else
        root1 = a;
                                                                                                                                                                           
      if(table.try_emplace(b, root1).second == false)
        {
          double root2 = table[b];
          while(table[root2] != root2)
            root2 = table[root2];

          if(root2 < root1)
            table[root1] = root2;
          else
            table[root2] = root1;
        }
    }

    static inline void join(const std::unordered_map<double,double>& join_table, std::unordered_map<double,double>& loc_eq_table)
    {
      for(auto& cur:loc_eq_table)
        {
          if(cur.first == cur.second && join_table.find(cur.first) != join_table.end())
            {
              loc_eq_table[cur.first] = join_table.at(cur.first);
              loc_eq_table[join_table.at(cur.first)] = join_table.at(cur.first);
            }
        }
      union_loc(loc_eq_table);
    }

    template<typename T=double>
    static inline void join_and_extract_local_thread_root(const std::unordered_map<T,T>& join_table, std::unordered_map<T,T>& loc_eq_table, std::unordered_set<T>& root_set)
    {
      for(auto& cur:loc_eq_table)
        {
          if(cur.first == cur.second && join_table.find(cur.first) != join_table.end())
            {
              loc_eq_table[cur.first] = join_table.at(cur.first);
              loc_eq_table[join_table.at(cur.first)] = join_table.at(cur.first);
            }
        }
      union_loc_and_root_find(loc_eq_table, root_set);
    }
    
    static inline size_t relab(std::unordered_map<double, double>& eq_table)
    {
      std::unordered_map<double, double> red_eq_table;

      size_t i = 0;
      for(auto& cur : eq_table)
        if (cur.first == cur.second) 
          {
            i++;
            red_eq_table[cur.first] = i;
          }

      for(auto& cur : eq_table)
        eq_table[cur.first] = red_eq_table[eq_table[cur.first]];

      return red_eq_table.size();
    }
    
    template<typename T>
    static inline void free_container(T& p_container)
    {
      T empty;
      std::swap(p_container, empty);
    }
    
    struct duet
    {
      ssize_t first;
      double  second;

      // bool operator==(const duet_key &other) const
      // {
      //   return (first == other.first && second == other.second)
      // }
    };
    
    static inline double convert(ssize_t in)
    {
      static_assert(sizeof(ssize_t) == sizeof(double), "double and ssize_t has different size :'(");
      double result;
      memcpy(&result, &in, sizeof(result));
      return result;
    }
    
    static inline ssize_t convert(double in)
    {
      static_assert(sizeof(ssize_t) == sizeof(double), "double and ssize_t has different size :'(");
      ssize_t result;
      memcpy(&result, &in, sizeof(result));
      return result;
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
