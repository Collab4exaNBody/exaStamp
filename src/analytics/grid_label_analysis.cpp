#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/domain.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/grid_algorithm.h>
#include <exanb/core/simple_block_rcb.h>
#include <exanb/core/physics_constants.h>

#include <exanb/core/grid_connected_components.h>

#include <string>
#include <unordered_map>
#include <utility>
#include <cstdint>
#include <mpi.h>
#include <omp.h>

namespace exaStamp
{
  using namespace exanb;
      
  class ConnectedGridAnalysis : public OperatorNode
  { 
    ADD_SLOT( MPI_Comm              , mpi                            , INPUT , MPI_COMM_WORLD , DocString{"MPI communicator"} );
    ADD_SLOT( Domain                , domain                         , INPUT , REQUIRED );
    ADD_SLOT( GridCellValues        , grid_cell_values               , INPUT_OUTPUT , REQUIRED );

    ADD_SLOT( ConnectedComponentMap , loc_label_map                  , INPUT , REQUIRED );
    ADD_SLOT( ssize_t               , nb_cc                          , INPUT , REQUIRED );
    
    ADD_SLOT( double                , nominal_density                , INPUT , 1. , DocString{"Standard density of your material in your used case. Will be used to evaluate voidness in a cell."} );

    ADD_SLOT( std::string           , filename                       , INPUT , REQUIRED );
    ADD_SLOT( std::string           , data_separator                 , INPUT , "," );
    ADD_SLOT( bool                  , enable_volume                  , INPUT , true );
    // ADD_SLOT( bool                  , enable_mass                    , INPUT , true );
    ADD_SLOT( bool                  , enable_mass_center             , INPUT , true );
    ADD_SLOT( bool                  , enable_geometric_center        , INPUT , true );
    ADD_SLOT( bool                  , enable_velocity                , INPUT , true );
    ADD_SLOT( bool                  , enable_density                 , INPUT , true );
    ADD_SLOT( bool                  , enable_temperature             , INPUT , true );
    ADD_SLOT( bool                  , enable_pressure                , INPUT , true );
    
    ADD_SLOT( bool                  , enable_debug                   , INPUT , false );
    ADD_SLOT( bool                  , grid_info                      , INPUT , true );

    ADD_SLOT( std::string           , filter                         , INPUT , REQUIRED );
    
    ADD_SLOT( std::string           , cc_label_field_name            , INPUT , REQUIRED );
    ADD_SLOT( std::string           , mass_field_name                , INPUT , REQUIRED );
    ADD_SLOT( std::string           , contribution_field_name        , INPUT , REQUIRED );
    ADD_SLOT( std::string           , momentum_vector_field_name     , INPUT , REQUIRED );
    ADD_SLOT( std::string           , m_v2_field_name                , INPUT , REQUIRED );
    ADD_SLOT( std::string           , tr_virial_field_name           , INPUT , REQUIRED );
    
    ADD_SLOT( std::string           , cc_volume_field_name           , INPUT_OUTPUT , "cc_volume (A^3)" );
    ADD_SLOT( std::string           , cc_density_field_name          , INPUT_OUTPUT , "cc_density (g / cm^3)" );
    ADD_SLOT( std::string           , cc_mass_center_field_name      , INPUT_OUTPUT , "cc_mass_center (A)" );
    ADD_SLOT( std::string           , cc_geometric_center_field_name , INPUT_OUTPUT , "cc_geometric_center (A)" );
    ADD_SLOT( std::string           , cc_velocity_field_name         , INPUT_OUTPUT , "cc_velocity (A/ps)" );
    ADD_SLOT( std::string           , cc_temperature_field_name      , INPUT_OUTPUT , "cc_temperature (K)" );
    ADD_SLOT( std::string           , cc_pressure_field_name         , INPUT_OUTPUT , "cc_pressure (GPa)" );
    
    inline bool is_sink() const override final { return true; } // not a suppressable operator
    
  public:

    bool first_time = true;
    
    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      static_assert( sizeof(ssize_t) == sizeof(double) ); // the module works only in this case
      static_assert( sizeof(long long) == sizeof(ssize_t) ); // ensures we can send ssize_t with MPI_LONG_LONG
      
      // ---- checks on the capacity of the module to execute it self as it has been ask by the user. (Users can be dumb. Especially me...)
      if( ! grid_cell_values->has_field(*cc_label_field_name) )
        {
          if(first_time == true)
            {
              lout << std::endl <<"WARNING !" << std::endl << std::endl << "You are trying to analyse a connected grid that hasn't been connected. This module will not execute itself."
                   << std::endl << std::endl;
              
              
              first_time = false;
            }
          return;
        }

      if( *enable_mass_center && ! grid_cell_values->has_field(*contribution_field_name) )
        {
          if(first_time == true)
            {
              lout << std::endl <<"WARNING !" << std::endl << std::endl << "You are trying to compute the center of mass from regions connected trough connected components on a grid that hasn't mass "
                   << "infomation. The module will skip this computation." << std::endl << std::endl;
            }
          *enable_mass_center = false;
        }

      if( *enable_geometric_center && *enable_volume == false )
        {
          if(first_time == true)
            {
              lout << std::endl <<"WARNING !" << std::endl << std::endl << "You are trying to compute the geometric center from regions connected trough connected components without having enable volume computation. The volume computation is needed. "
                   << "The module will skip this computation." << std::endl << std::endl;
            }
          *enable_geometric_center = false;
        }
      
      if( *enable_velocity &&  ! grid_cell_values->has_field(*momentum_vector_field_name) )
        {
          if(first_time == true)
            {
              lout << std::endl <<"WARNING !" << std::endl << std::endl << "You are trying to compute the velocity of regions connected trough connected components on a grid that hasn't momentum vector "
                   << "infomation. The module will skip this computation." << std::endl << std::endl;
            }
          *enable_velocity = false;
        }
      
      if( *enable_density &&  ! grid_cell_values->has_field(*mass_field_name) )
        {
          if(first_time == true)
            {
              lout << std::endl <<"WARNING !" << std::endl << std::endl << "You are trying to compute the density of regions connected trough connected components on a grid that hasn't mass "
                   << "infomation. The module will skip this computation." << std::endl << std::endl;
            }
          *enable_density = false;
        }

      if(*enable_volume == false && *enable_density == false && *enable_mass_center == false && *enable_velocity == false && *enable_geometric_center == false)
        {
          if(first_time == true)
            {
              lout << std::endl <<"WARNING !" << std::endl << std::endl << "There is no computation to do or that can be done in this module (connected_grid_analysis)!"
                   << std::endl << std::endl;
            }
          
          first_time = false;
          return;
        }
      first_time = false;
      // ---- checks' end
      
      MPI_Comm comm = *mpi;

      // Constants definition :

      // -- Domain :
      // simulation's grid dimensions
      // const IJK domain_dims = domain->grid_dimension();
      // cell size (all cells are cubical)
      const double cell_size = domain->cell_size(); 
      
      // -- Grid : 
      // local processor's grid dimensions, including ghost cells
      const IJK grid_dims = grid_cell_values->grid_dims(); 
      // local processor's grid position in simulation's grid
      const IJK grid_offset = grid_cell_values->grid_offset();
      // ghost layers
      const ssize_t gl = grid_cell_values->ghost_layers();
      // number of subdivisions, in each directions, applied on cells
      const ssize_t subdiv = grid_cell_values->field(*contribution_field_name).m_subdiv;
      // local processor's total grid dimension 
      // const IJK total_dims = grid_dims * subdiv;

      // -- Usefull constants construction :
      // local grid boundaries
      // const GridBlock local_space = {grid_offset, grid_offset + grid_dims};
      // number of sub-cells per cell
      //const size_t n_subcells = subdiv*subdiv*subdiv;
      // dimension of the subdivided simulation's grid
      //const IJK domain_subdiv_dims = domain_dims * subdiv; 
      // side size of a sub-cell
      const double subcell_size = cell_size / subdiv;
      // sub-cell's volume
      const double subcell_vol  = subcell_size * subcell_size * subcell_size;

      double tmp_inv_nominal_density = 0.0;
      if( (*filter) == "void" )
        {
          tmp_inv_nominal_density = 1.0 / (*nominal_density);
        }
      const double inv_nominal_density = tmp_inv_nominal_density; 
      
      // // -----------------------------------------------
      // // some debug information
      // ldbg << "ghost_layers="<<gl<<", cell_size="<<cell_size<<", subdiv="<<subdiv<<", subcell_size="
      //      <<subcell_size<<", grid_dims="<<grid_dims<<", grid_offset="<<grid_offset<<", domain_dims="<<domain_dims
      //      <<", domain_subdiv_dims="<<domain_subdiv_dims<<std::endl;
      // // -----------------------------------------------

      // Grid data fields :
      // note: if we are to add new data fields, they must be added BEFORE we retreive access ANY information

      // -- Creation :
      // additional data field for connected component label
      if(*grid_info == true)
        {
          if( *enable_volume && ! grid_cell_values->has_field(*cc_volume_field_name) )
            {
              grid_cell_values->add_field(*cc_volume_field_name,subdiv,1);
            }
          if( *enable_density && ! grid_cell_values->has_field(*cc_density_field_name) )
            {
              grid_cell_values->add_field(*cc_density_field_name,subdiv,1);
            }
          if( *enable_mass_center && ! grid_cell_values->has_field(*cc_mass_center_field_name) )
            {
              grid_cell_values->add_field(*cc_mass_center_field_name,subdiv,3);
            }
          if( *enable_geometric_center && ! grid_cell_values->has_field(*cc_geometric_center_field_name) )
            {
              grid_cell_values->add_field(*cc_geometric_center_field_name,subdiv,3);
            }
          if( *enable_velocity && ! grid_cell_values->has_field(*cc_velocity_field_name) )
            {
              grid_cell_values->add_field(*cc_velocity_field_name,subdiv,3);
            }
          if( *enable_temperature && ! grid_cell_values->has_field(*cc_temperature_field_name) )
            {
              grid_cell_values->add_field(*cc_temperature_field_name,subdiv,1);
            }
          if( *enable_pressure && ! grid_cell_values->has_field(*cc_pressure_field_name) )
            {
              grid_cell_values->add_field(*cc_pressure_field_name,subdiv,1);
            }
        }
      
      int nb_field = 0;
      
      // -- Retrieve grid data field accessor :
      // cc_label field data accessor. <--- Is needed in every case.
      const auto cc_label_accessor = grid_cell_values->field_data(*cc_label_field_name);
      // uint64_t * __restrict__ cc_label_ptr = (uint64_t *) cc_label_accessor.m_data_ptr;
      const double * __restrict__ cc_label_ptr = cc_label_accessor.m_data_ptr;
      const size_t stride = cc_label_accessor.m_stride; // for simplicity

      nb_field += 1;
      
      // contribution field data accessor.
      double * __restrict__ tmp_ctb_ptr = nullptr;
      if (*enable_velocity || *enable_mass_center || *enable_density) // enable density ?
        {
          tmp_ctb_ptr = grid_cell_values->field_data(*contribution_field_name).m_data_ptr;
          nb_field += 1;
        }
      const double * __restrict__ ctb_ptr = tmp_ctb_ptr;

      // momentum field data accessor.
      double * __restrict__ tmp_momentum_ptr = nullptr;
      if (*enable_velocity)
        {
          tmp_momentum_ptr = grid_cell_values->field_data(*momentum_vector_field_name).m_data_ptr;
          nb_field += 3;
        }
      const double * __restrict__ momentum_ptr = tmp_momentum_ptr;

      // mass field data accessor.
      double * __restrict__ tmp_mass_ptr = nullptr;
      if (*enable_density || *enable_volume)
        {
          tmp_mass_ptr = grid_cell_values->field_data(*mass_field_name).m_data_ptr;
          nb_field += 1;
        }
      const double * __restrict__ mass_ptr = tmp_mass_ptr;

      if(*enable_volume)
        nb_field += 1;

      if(*enable_mass_center)
        nb_field += 3;

      if(*enable_geometric_center)
        nb_field += 3;

      double * __restrict__ tmp_m_v2_ptr = nullptr;
      if(*enable_temperature)
        {
          tmp_m_v2_ptr = grid_cell_values->field_data(*m_v2_field_name).m_data_ptr;
          nb_field += 1;
        }
      const double * __restrict__ m_v2_ptr = tmp_m_v2_ptr;

      double * __restrict__ tmp_tr_virial_ptr = nullptr;
      if(*enable_pressure)
        {
          tmp_tr_virial_ptr = grid_cell_values->field_data(*tr_virial_field_name).m_data_ptr;
          nb_field += 1;
        }
      const double * __restrict__ tr_virial_ptr = tmp_tr_virial_ptr;
      
      // -----------------------------------------------
      
      ssize_t loc_nb_label = (*loc_label_map).size();
      
      double *local_connected_components = (double *)calloc(loc_nb_label * nb_field, sizeof(double));

      double *contribution       = nullptr;
      double *volume             = nullptr;
      double *mass               = nullptr;
      double *mass_center_x      = nullptr;
      double *mass_center_y      = nullptr;
      double *mass_center_z      = nullptr;
      double *geometric_center_x = nullptr;
      double *geometric_center_y = nullptr;
      double *geometric_center_z = nullptr;
      double *mvx                = nullptr;
      double *mvy                = nullptr;
      double *mvz                = nullptr;
      double *temperature        = nullptr;
      double *pressure           = nullptr;

      {
        int i = 1;
        if(*enable_velocity || *enable_mass_center || *enable_density)
          {
            contribution       = &local_connected_components[loc_nb_label *  i   ];
            i += 1;
          }
        if(*enable_volume)
          {
            volume             = &local_connected_components[loc_nb_label *  i   ];
            i += 1;
          }
        if(*enable_density)
          {
            mass               = &local_connected_components[loc_nb_label *  i   ];
            i += 1;
          }
        if(*enable_mass_center)
          {
            mass_center_x      = &local_connected_components[loc_nb_label *  i   ];
            mass_center_y      = &local_connected_components[loc_nb_label * (i+1)];
            mass_center_z      = &local_connected_components[loc_nb_label * (i+2)];
            i += 3;
          }
        if(*enable_geometric_center)
          {
            geometric_center_x = &local_connected_components[loc_nb_label *  i   ];
            geometric_center_y = &local_connected_components[loc_nb_label * (i+1)];
            geometric_center_z = &local_connected_components[loc_nb_label * (i+2)];
            i += 3;
          }
        if(*enable_velocity)
          {
            mvx                = &local_connected_components[loc_nb_label *  i   ];
            mvy                = &local_connected_components[loc_nb_label * (i+1)];
            mvz                = &local_connected_components[loc_nb_label * (i+2)];
            i += 3;
          }
        if(*enable_temperature)
          {
            temperature        = &local_connected_components[loc_nb_label *  i   ];
            i += 1;
          }
        if(*enable_pressure)
          {
            pressure           = &local_connected_components[loc_nb_label *  i   ];
            i += 1;
          }
        assert( i == nb_field);
      }
      
      double *cc                    = nullptr;
      
      double *cc_volume             = nullptr;
      double *cc_contribution       = nullptr;
      double *cc_density            = nullptr;
      double *cc_mass_center_x      = nullptr;
      double *cc_mass_center_y      = nullptr;
      double *cc_mass_center_z      = nullptr;
      double *cc_geometric_center_x = nullptr;
      double *cc_geometric_center_y = nullptr;
      double *cc_geometric_center_z = nullptr;
      double *cc_vx                 = nullptr;
      double *cc_vy                 = nullptr;
      double *cc_vz                 = nullptr;
      double *cc_temperature        = nullptr;
      double *cc_pressure           = nullptr;
      
      int rank=0, np = 0;
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &np); 

      // int OMP_NUM_THREADS = omp_get_max_threads(); 
      
      // #     pragma omp parallel
      //       {
      //         int nthreads = omp_get_num_threads(); 
      //         int thread_id = omp_get_thread_num();
        
      //        assert(nthreads == OMP_NUM_THREADS); // check that there is the number of threads that we expect

      //         // Each thread retrieve its working block
      //         GridBlock working_block = simple_block_rcb(GridBlock{IJK{gl,gl,gl}, grid_dims-gl}, nthreads, thread_id);

      GridBlock working_block = GridBlock{IJK{gl,gl,gl}, grid_dims-gl};
      for( ssize_t k=working_block.start.k ; k < working_block.end.k ; k++)
        for( ssize_t j=working_block.start.j ; j < working_block.end.j ; j++)
          for( ssize_t i=working_block.start.i ; i < working_block.end.i ; i++)
            {
              ssize_t cell_index = grid_ijk_to_index( grid_dims , IJK{i,j,k} );

              for( ssize_t sk=0 ; sk<subdiv ; sk++)
                for( ssize_t sj=0 ; sj<subdiv ; sj++)
                  for( ssize_t si=0 ; si<subdiv ; si++)
                    {
                      ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , IJK{si,sj,sk} );

                      ssize_t value_index = cell_index * stride + subcell_index;

                      if(cc_label_ptr[value_index] > 0.0)
                        {
                          ssize_t scindex_vec_0 = cell_index * stride + subcell_index * 3 + 0;
                          ssize_t scindex_vec_1 = cell_index * stride + subcell_index * 3 + 1;
                          ssize_t scindex_vec_2 = cell_index * stride + subcell_index * 3 + 2;

                          ssize_t loc_id = (*loc_label_map).at(cc_label_ptr[value_index]); 

                          if(*enable_velocity || *enable_mass_center || *enable_density)
                            {
                              contribution       [loc_id] += ctb_ptr[value_index];
                            }
                          if(*enable_volume)
                            {
                              volume             [loc_id] += subcell_vol - mass_ptr[value_index] * inv_nominal_density;
                            }
                          if(*enable_density)
                            {
                              mass               [loc_id] += mass_ptr[value_index];
                            }
                          if(*enable_mass_center)
                            {
                              mass_center_x      [loc_id] += ((i + grid_offset.i) * subdiv + si) * subcell_size * mass_ptr[value_index];
                              mass_center_y      [loc_id] += ((j + grid_offset.j) * subdiv + sj) * subcell_size * mass_ptr[value_index];
                              mass_center_z      [loc_id] += ((k + grid_offset.k) * subdiv + sk) * subcell_size * mass_ptr[value_index];
                            }
                          if(*enable_geometric_center)
                            {
                              geometric_center_x [loc_id] += ((i + grid_offset.i) * subdiv + si) * subcell_size;
                              geometric_center_y [loc_id] += ((j + grid_offset.j) * subdiv + sj) * subcell_size;
                              geometric_center_z [loc_id] += ((k + grid_offset.k) * subdiv + sk) * subcell_size;
                            }
                          if(*enable_velocity)
                            {
                              mvx                [loc_id] += momentum_ptr[scindex_vec_0];
                              mvy                [loc_id] += momentum_ptr[scindex_vec_1];
                              mvz                [loc_id] += momentum_ptr[scindex_vec_2];
                            }
                          if(*enable_temperature)
                            {
                              temperature        [loc_id] += m_v2_ptr[value_index];
                            }
                          if(*enable_pressure)
                            {
                              pressure           [loc_id] += tr_virial_ptr[value_index];
                            }
                        }
                    }
            }
        
        
      // #       pragma omp single
      //        {
      int msg_elements = loc_nb_label;
      int *msgs_size = nullptr;
          
      MPI_Request requests[1];
      MPI_Status status[1];
          
      int *displs = nullptr, *vcounts = nullptr;
      double *global_buffer = nullptr;

      if(np > 1 && rank != 0)
        {
          MPI_Igather(&msg_elements, 1, MPI_INT, msgs_size, 1, MPI_INT, 0, comm, &requests[0]);
          {
            for(auto& cur:(*loc_label_map))
              {
                local_connected_components[cur.second] = cur.first;
              }
          }
          MPI_Wait(&requests[0], &status[0]);
          MPI_Gatherv(local_connected_components, nb_field * loc_nb_label, MPI_DOUBLE, global_buffer, vcounts, displs, MPI_DOUBLE, 0, comm);

          if(*grid_info == true)
            {
              assert(&local_connected_components[loc_nb_label] == memset(&local_connected_components[loc_nb_label], 0, sizeof(double)*loc_nb_label*(nb_field-1)));
              MPI_Scatterv(global_buffer, vcounts, displs, MPI_DOUBLE, local_connected_components, nb_field * loc_nb_label, MPI_DOUBLE, 0, comm);
            }
        }
      if(np > 1 && rank == 0)
        {
          msgs_size = new int[np];
          MPI_Igather(&msg_elements, 1, MPI_INT, msgs_size, 1, MPI_INT, 0, comm, &requests[0]);
          {
            for(auto& cur:(*loc_label_map))
              {
                local_connected_components[cur.second] = cur.first;
              }
          }
          displs = new int[np];
          vcounts = new int[np];

          MPI_Wait(&requests[0], &status[0]);

          displs[0] = 0;
          vcounts[0] = nb_field * msgs_size[0];
          for(int i = 1; i < np; i++)
            {
              displs[i] = displs[i-1] + vcounts[i-1];
              vcounts[i] = nb_field * msgs_size[i];
            }
            
          global_buffer = new double[displs[np-1] + vcounts[np-1]];

          MPI_Gatherv(local_connected_components, nb_field * loc_nb_label, MPI_DOUBLE, global_buffer, vcounts, displs, MPI_DOUBLE, 0, comm);

          // TO DO : move the next part in a thread parallel region.

          cc = (double *) calloc((nb_field-1) * (*nb_cc) , sizeof(double));

          {
            int i = 0;
            if(*enable_velocity || *enable_mass_center || *enable_density)
              {
                cc_contribution  = &cc[(*nb_cc) *  i   ];
                i += 1;
              }
            if(*enable_volume)
              {
                cc_volume        = &cc[(*nb_cc) *  i   ];
                i += 1;
              }
            if(*enable_density)
              {
                cc_density       = &cc[(*nb_cc) *  i   ];
                i += 1;
              }
            if(*enable_mass_center)
              {
                cc_mass_center_x = &cc[(*nb_cc) *  i   ];
                cc_mass_center_y = &cc[(*nb_cc) * (i+1)];
                cc_mass_center_z = &cc[(*nb_cc) * (i+2)];
                i += 3;
              }
            if(*enable_geometric_center)
              {
                cc_geometric_center_x = &cc[(*nb_cc) *  i   ];
                cc_geometric_center_y = &cc[(*nb_cc) * (i+1)];
                cc_geometric_center_z = &cc[(*nb_cc) * (i+2)];
                i += 3;
              }
            if(*enable_velocity)
              {
                cc_vx            = &cc[(*nb_cc) *  i   ];
                cc_vy            = &cc[(*nb_cc) * (i+1)];
                cc_vz            = &cc[(*nb_cc) * (i+2)];
                i += 3;
              }
            if(*enable_temperature)
              {
                cc_temperature   = &cc[(*nb_cc) *  i   ];
                i += 1;
              }
            if(*enable_pressure)
              {
                cc_pressure      = &cc[(*nb_cc) *  i   ];
                i += 1;
              }
            assert((nb_field-1) == i);
          }

          double **g_label = new double*[np];
          
          double **g_contribution       = nullptr;
          double **g_volume             = nullptr;
          double **g_density            = nullptr;
          double **g_mass_center_x      = nullptr;
          double **g_mass_center_y      = nullptr;
          double **g_mass_center_z      = nullptr;
          double **g_geometric_center_x = nullptr;
          double **g_geometric_center_y = nullptr;
          double **g_geometric_center_z = nullptr;
          double **g_vx                 = nullptr;
          double **g_vy                 = nullptr;
          double **g_vz                 = nullptr;
          double **g_temperature        = nullptr;
          double **g_pressure           = nullptr;

          if(*enable_velocity || *enable_mass_center || *enable_density)
            g_contribution       = new double*[np];
          if(*enable_volume)
            g_volume             = new double*[np];
          if(*enable_density)
            g_density            = new double*[np];
          if(*enable_mass_center)
            {
              g_mass_center_x      = new double*[np];
              g_mass_center_y      = new double*[np];
              g_mass_center_z      = new double*[np];
            }
          if(*enable_geometric_center)
            {
              g_geometric_center_x = new double*[np];
              g_geometric_center_y = new double*[np];
              g_geometric_center_z = new double*[np];
            }
          if(*enable_velocity)
            {
              g_vx                 = new double*[np];
              g_vy                 = new double*[np];
              g_vz                 = new double*[np];
            }
          if(*enable_temperature)
            g_temperature          = new double*[np];
          if(*enable_pressure)
            g_pressure             = new double*[np];
          
          for(int i = 0; i < np; i++)
            {
              g_label[i] = &global_buffer[displs[i]];

              int j = 1;
              if(*enable_velocity || *enable_mass_center || *enable_density)
                {
                  g_contribution[i]       = g_label[i] +  j    * msgs_size[i];
                  j += 1;
                }
              if(*enable_volume)
                {
                  g_volume[i]             = g_label[i] +  j    * msgs_size[i];
                  j += 1;
                }
              if(*enable_density)
                {
                  g_density[i]            = g_label[i] +  j    * msgs_size[i];
                  j += 1;
                }
              if(*enable_mass_center)
                {
                  g_mass_center_x[i]      = g_label[i] +  j    * msgs_size[i];
                  g_mass_center_y[i]      = g_label[i] + (j+1) * msgs_size[i];
                  g_mass_center_z[i]      = g_label[i] + (j+2) * msgs_size[i];
                  j += 3;
                }
              if(*enable_geometric_center)
                {
                  g_geometric_center_x[i] = g_label[i] +  j    * msgs_size[i];
                  g_geometric_center_y[i] = g_label[i] + (j+1) * msgs_size[i];
                  g_geometric_center_z[i] = g_label[i] + (j+2) * msgs_size[i];
                  j += 3;
                }
              if(*enable_velocity)
                {
                  g_vx[i]                 = g_label[i] +  j    * msgs_size[i];
                  g_vy[i]                 = g_label[i] + (j+1) * msgs_size[i];
                  g_vz[i]                 = g_label[i] + (j+2) * msgs_size[i];
                  j += 3;
                }
              if(*enable_temperature)
                {
                  g_temperature[i]        = g_label[i] +  j    * msgs_size[i];
                  j += 1;
                }
              if(*enable_pressure)
                {
                  g_pressure[i]           = g_label[i] +  j    * msgs_size[i];
                  j += 1;
                }
              assert(j == nb_field);
            }
          
          for(int i = 0; i < np; i++)
            {
              for(int j = 0; j < msgs_size[i]; j++)
                {
                  double lab = g_label[i][j];
                  ssize_t lab_key = ((ssize_t) lab) - 1;
                  assert(lab_key >= 0 && lab_key < *nb_cc );
                  
                  if(*enable_velocity || *enable_mass_center || *enable_density)
                    {
                      cc_contribution       [lab_key] += g_contribution[i][j];
                    }
                  if(*enable_volume)
                    {
                      cc_volume             [lab_key] += g_volume[i][j];
                    }
                  if(*enable_density)
                    {
                      cc_density            [lab_key] += g_density[i][j]; // the mass
                    }
                  if(*enable_mass_center)
                    {
                      cc_mass_center_x      [lab_key] += g_mass_center_x[i][j];
                      cc_mass_center_y      [lab_key] += g_mass_center_y[i][j];
                      cc_mass_center_z      [lab_key] += g_mass_center_z[i][j];
                    }
                  if(*enable_geometric_center)
                    {
                      cc_geometric_center_x [lab_key] += g_geometric_center_x[i][j];
                      cc_geometric_center_y [lab_key] += g_geometric_center_y[i][j];
                      cc_geometric_center_z [lab_key] += g_geometric_center_z[i][j];
                    }
                  if(*enable_velocity)
                    {
                      cc_vx                 [lab_key] += g_vx[i][j];
                      cc_vy                 [lab_key] += g_vy[i][j];
                      cc_vz                 [lab_key] += g_vz[i][j];
                    }
                  if(*enable_temperature)
                    {
                      cc_temperature        [lab_key] += g_temperature[i][j];
                    }
                  if(*enable_pressure)
                    {
                      cc_pressure           [lab_key] += g_pressure[i][j];
                    }
                }
            }

          const double density_cst = 1.660539066;
          const double temp_cst = legacy_constant::atomicMass * 10000.0 / (3.0*legacy_constant::boltzmann);
          const double pressure_cst = 1.660539066 * 0.01;
          for(size_t i = 0; i < (*nb_cc); i++)
            {
              if(*enable_mass_center)
                {
                  cc_mass_center_x      [i] /= cc_density[i];
                  cc_mass_center_y      [i] /= cc_density[i];
                  cc_mass_center_z      [i] /= cc_density[i];
                }
              if(*enable_geometric_center)
                {
                  cc_geometric_center_x [i] /= (cc_volume[i] / subcell_vol); // (cc_volume[i] / subcell_vol) = nb_cell
                  cc_geometric_center_y [i] /= (cc_volume[i] / subcell_vol);
                  cc_geometric_center_z [i] /= (cc_volume[i] / subcell_vol);
                }
              if(*enable_velocity)
                {
                  cc_vx                 [i] /= cc_density[i];
                  cc_vy                 [i] /= cc_density[i];
                  cc_vz                 [i] /= cc_density[i];
                }
              if(*enable_temperature)
                {
                  cc_temperature        [i] = (cc_temperature[i] - cc_density[i] * (cc_vx[i] * cc_vx[i] + cc_vy[i] * cc_vy[i] + cc_vz[i] * cc_vz[i]));
                }
              if(*enable_pressure)
                {
                  cc_pressure           [i] = (cc_pressure[i] + cc_temperature[i]) / (cc_volume[i] * 3.0);
                }

              if(*enable_density)
                {
                  cc_density            [i] /= cc_volume[i]; // mass to density
                  cc_density            [i] *= density_cst; // in g / cm^3
                }
              if(*enable_temperature)
                {
                  cc_temperature        [i] *= temp_cst; // in K
                }
              if(*enable_pressure)
                {
                  cc_pressure           [i] *= pressure_cst; // in GPa
                }
            }
          
          if(*grid_info == true)
            {
              for(int i = 0; i < np; i++)
                {
                  for(int j = 0; j < msgs_size[i]; j++)
                    {
                      double lab = g_label[i][j];
                      ssize_t lab_key = ((ssize_t) lab) - 1;
                      assert(lab_key >= 0 && lab_key < *nb_cc );
                  
                      if(*enable_velocity || *enable_mass_center || *enable_density)
                        {
                          g_contribution[i][j] = cc_contribution  [lab_key];
                        }
                      if(*enable_volume)
                        {
                          g_volume[i][j] = cc_volume[lab_key];
                        }
                      if(*enable_density)
                        {
                          g_density[i][j] = cc_density[lab_key];
                        }
                      if(*enable_mass_center)
                        {
                          g_mass_center_x[i][j] = cc_mass_center_x [lab_key];
                          g_mass_center_y[i][j] = cc_mass_center_y [lab_key];
                          g_mass_center_z[i][j] = cc_mass_center_z [lab_key];
                        }
                      if(*enable_geometric_center)
                        {
                          g_geometric_center_x[i][j] = cc_geometric_center_x [lab_key];
                          g_geometric_center_y[i][j] = cc_geometric_center_y [lab_key];
                          g_geometric_center_z[i][j] = cc_geometric_center_z [lab_key];
                        }
                      if(*enable_velocity)
                        {
                          g_vx[i][j] = cc_vx [lab_key];
                          g_vy[i][j] = cc_vy [lab_key];
                          g_vz[i][j] = cc_vz [lab_key];
                        }
                      if(*enable_temperature)
                        {
                          g_temperature[i][j] = cc_temperature[lab_key];
                        }
                      if(*enable_pressure)
                        {
                          g_pressure[i][j] = cc_pressure[lab_key];
                        }
                    }
                }
              MPI_Scatterv(global_buffer, vcounts, displs, MPI_DOUBLE, local_connected_components, nb_field * loc_nb_label, MPI_DOUBLE, 0, comm);
            }
              
          delete msgs_size;
          delete displs;
          delete vcounts;
          
          FILE *output = fopen((*filename).c_str(), "w");
          fprintf(output,"id");
          if(*enable_volume)
            {
              fprintf(output,(*data_separator).c_str());
              fprintf(output,"volume(A^3)");
            }
          if(*enable_density)
            {
              fprintf(output,(*data_separator).c_str());
              fprintf(output,"density(g/cm^3)");
            }
          if(*enable_mass_center)
            {
              fprintf(output,(*data_separator).c_str());
              fprintf(output,"mass_center_x(A)");
              fprintf(output,(*data_separator).c_str());
              fprintf(output,"mass_center_y(A)");
              fprintf(output,(*data_separator).c_str());
              fprintf(output,"mass_center_z(A)");
            }
          if(*enable_geometric_center)
            {
              fprintf(output,(*data_separator).c_str());
              fprintf(output,"geometric_center_x(A)");
              fprintf(output,(*data_separator).c_str());
              fprintf(output,"geometric_center_y(A)");
              fprintf(output,(*data_separator).c_str());
              fprintf(output,"geometric_center_z(A)");
            }
          if(*enable_velocity)
            {
              fprintf(output,(*data_separator).c_str());
              fprintf(output,"vx(A/ps)");
              fprintf(output,(*data_separator).c_str());
              fprintf(output,"vy(A/ps)");
              fprintf(output,(*data_separator).c_str());
              fprintf(output,"vz(A/ps)");
            }
          if(*enable_temperature)
            {
              fprintf(output,(*data_separator).c_str());
              fprintf(output,"temperature(K)");
            }
          if(*enable_pressure)
            {
              fprintf(output,(*data_separator).c_str());
              fprintf(output,"pressure(GPa)");
            }
          fprintf(output,"\n");
          
          for(size_t i = 0; i < size_t(*nb_cc); i++)
            {
              fprintf(output,"%ld",i+1);
              if(*enable_volume)
                {
                  fprintf(output,(*data_separator).c_str());
                  fprintf(output,"%.16lf",cc_volume[i]);
                }
              if(*enable_density)
                {
                  fprintf(output,(*data_separator).c_str());
                  fprintf(output,"%.16lf",cc_density[i]);
                }
              if(*enable_mass_center)
                {
                  fprintf(output,(*data_separator).c_str());
                  fprintf(output,"%.16lf",cc_mass_center_x[i]);
                  fprintf(output,(*data_separator).c_str());
                  fprintf(output,"%.16lf",cc_mass_center_y[i]);
                  fprintf(output,(*data_separator).c_str());
                  fprintf(output,"%.16lf",cc_mass_center_z[i]);
                }
              if(*enable_geometric_center)
                {
                  fprintf(output,(*data_separator).c_str());
                  fprintf(output,"%.16lf",cc_geometric_center_x[i]);
                  fprintf(output,(*data_separator).c_str());
                  fprintf(output,"%.16lf",cc_geometric_center_y[i]);
                  fprintf(output,(*data_separator).c_str());
                  fprintf(output,"%.16lf",cc_geometric_center_z[i]);
                }
              if(*enable_velocity)
                {
                  fprintf(output,(*data_separator).c_str());
                  fprintf(output,"%.16lf",cc_vx[i]);
                  fprintf(output,(*data_separator).c_str());
                  fprintf(output,"%.16lf",cc_vy[i]);
                  fprintf(output,(*data_separator).c_str());
                  fprintf(output,"%.16lf",cc_vz[i]);
                }
              if(*enable_temperature)
                {
                  fprintf(output,(*data_separator).c_str());
                  fprintf(output,"%.16lf", cc_temperature[i]);
                }
              if(*enable_pressure)
                {
                  fprintf(output,(*data_separator).c_str());
                  fprintf(output,"%.16lf", cc_pressure[i]);
                }
              fprintf(output,"\n");
            }
          fclose(output);
        }

      if(*grid_info == true)
        {
          // retrieve field data accessor for grid info.
          double * __restrict__ cc_volume_ptr = nullptr;
          if (*enable_volume)
            {
              cc_volume_ptr = grid_cell_values->field_data(*cc_volume_field_name).m_data_ptr;
            }
          double * __restrict__ cc_density_ptr = nullptr;
          if (*enable_density)
            {
              cc_density_ptr = grid_cell_values->field_data(*cc_density_field_name).m_data_ptr;
            }
          double * __restrict__ cc_mass_center_ptr = nullptr;
          if (*enable_mass_center)
            {
              cc_mass_center_ptr = grid_cell_values->field_data(*cc_mass_center_field_name).m_data_ptr;
            }
          double * __restrict__ cc_geometric_center_ptr = nullptr;
          if (*enable_geometric_center)
            {
              cc_geometric_center_ptr = grid_cell_values->field_data(*cc_geometric_center_field_name).m_data_ptr;
            }
          double * __restrict__ cc_velocity_ptr = nullptr;
          if (*enable_velocity)
            {
              cc_velocity_ptr = grid_cell_values->field_data(*cc_velocity_field_name).m_data_ptr;
            }
          double * __restrict__ cc_temperature_ptr = nullptr;
          if (*enable_temperature)
            {
              cc_temperature_ptr = grid_cell_values->field_data(*cc_temperature_field_name).m_data_ptr;
            }
          double * __restrict__ cc_pressure_ptr = nullptr;
          if (*enable_pressure)
            {
              cc_pressure_ptr = grid_cell_values->field_data(*cc_pressure_field_name).m_data_ptr;
            }

          for( ssize_t k=working_block.start.k ; k < working_block.end.k ; k++)
            for( ssize_t j=working_block.start.j ; j < working_block.end.j ; j++)
              for( ssize_t i=working_block.start.i ; i < working_block.end.i ; i++)
                {
                  // IJK cell_location = IJK{i,j,k} + grid_offset;

                  for( ssize_t sk=0 ; sk<subdiv ; sk++)
                    for( ssize_t sj=0 ; sj<subdiv ; sj++)
                      for( ssize_t si=0 ; si<subdiv ; si++)
                        {
                          ssize_t cell_index = grid_ijk_to_index( grid_dims , IJK{i,j,k} );
                          ssize_t subcell_index = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , IJK{si,sj,sk} );

                          ssize_t value_index = cell_index * stride + subcell_index;

                          ssize_t scindex_vec_0 = cell_index * stride + subcell_index * 3 + 0;
                          ssize_t scindex_vec_1 = cell_index * stride + subcell_index * 3 + 1;
                          ssize_t scindex_vec_2 = cell_index * stride + subcell_index * 3 + 2;

                          if(cc_label_ptr[value_index] > 0.0)
                            {
                              ssize_t loc_id = (*loc_label_map).at(cc_label_ptr[value_index]);

                              // assert(loc_id == (ssize_t)local_connected_components[loc_id]);
                              
                              if( *enable_volume )
                                {
                                  cc_volume_ptr[value_index] = volume[loc_id];
                                }
                              if( *enable_density )
                                {
                                  cc_density_ptr[value_index] = mass[loc_id];
                                }
                              if( *enable_mass_center )
                                {
                                  cc_mass_center_ptr[scindex_vec_0] = mass_center_x[loc_id];
                                  cc_mass_center_ptr[scindex_vec_1] = mass_center_y[loc_id];
                                  cc_mass_center_ptr[scindex_vec_2] = mass_center_z[loc_id];
                                }
                              if( *enable_geometric_center )
                                {
                                  cc_geometric_center_ptr[scindex_vec_0] = geometric_center_x[loc_id];
                                  cc_geometric_center_ptr[scindex_vec_1] = geometric_center_y[loc_id];
                                  cc_geometric_center_ptr[scindex_vec_2] = geometric_center_z[loc_id];
                                }
                              if( *enable_velocity )
                                {
                                  cc_velocity_ptr[scindex_vec_0] = mvx[loc_id];
                                  cc_velocity_ptr[scindex_vec_1] = mvy[loc_id];
                                  cc_velocity_ptr[scindex_vec_2] = mvz[loc_id];
                                }
                              if( *enable_temperature )
                                {
                                  cc_temperature_ptr[value_index] = temperature[loc_id];
                                }
                              if( *enable_pressure )
                                {
                                  cc_pressure_ptr[value_index] = pressure[loc_id];
                                }
                            }
                          else
                            {
                              if( *enable_volume )
                                {
                                  cc_volume_ptr[value_index] = 0.0;
                                }
                              if( *enable_density )
                                {
                                  cc_density_ptr[value_index] = 0.0;
                                }
                              if( *enable_mass_center )
                                {
                                  cc_mass_center_ptr[scindex_vec_0] = 0.0;
                                  cc_mass_center_ptr[scindex_vec_1] = 0.0;
                                  cc_mass_center_ptr[scindex_vec_2] = 0.0;
                                }
                              if( *enable_geometric_center )
                                {
                                  cc_geometric_center_ptr[scindex_vec_0] = 0.0;
                                  cc_geometric_center_ptr[scindex_vec_1] = 0.0;
                                  cc_geometric_center_ptr[scindex_vec_2] = 0.0;
                                }
                              if( *enable_velocity )
                                {
                                  cc_velocity_ptr[scindex_vec_0] = 0.0;
                                  cc_velocity_ptr[scindex_vec_1] = 0.0;
                                  cc_velocity_ptr[scindex_vec_2] = 0.0;
                                }
                              if( *enable_temperature )
                                {
                                  cc_temperature_ptr[value_index] = 0.0;
                                }
                              if( *enable_pressure )
                                {
                                  cc_pressure_ptr[value_index] = 0.0;
                                }
                            }
                        }
                }
        }
      // }// end pragma omp single
        
      // } // end omp parallel 
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
    
    // static inline void find_par(const IJK& domain_dims,
    //                             const IJK& grid_dims,
    //                             const IJK& offset,
    //                             const GridBlock& working_block,
    //                             const IJK& cell_loc, 
    //                             const IJK& subcell_loc, 
    //                             const ssize_t subdiv, 
    //                             const size_t stride,
    //                             double * __restrict__ cc_label_ptr,
    //                             const double  * __restrict__ density_ptr,
    //                             double density_threshold,
    //                             std::unordered_map<double, double>& eq_table
    //                             )
    // {
    //   size_t id = grid_ijk_to_index(grid_dims, cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, subcell_loc);
    //   size_t label_boundary = grid_ijk_to_index(domain_dims * subdiv, domain_dims * subdiv) + 2;
    //   double prov = label_boundary;
    //   for( ssize_t nk=-1 ; nk < 2 ; nk++)
    //     for( ssize_t nj=-1 ; nj < 2; nj++)
    //       for( ssize_t ni=-1 ; ni < 2; ((nk==0 && nj==0)? ni+=2 :ni++))
    //         {
    //           IJK nbh_cell_loc;
    //           IJK nbh_subcell_loc;
    //           subcell_neighbor( cell_loc, subcell_loc, subdiv, IJK{ni,nj,nk}, nbh_cell_loc, nbh_subcell_loc );
    //           size_t nbh_id = grid_ijk_to_index(grid_dims, nbh_cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
    //           if(is_inside(working_block, nbh_cell_loc) && density_ptr[nbh_id] > density_threshold && cc_label_ptr[nbh_id] != 0.0)
    //             {
    //               double cursor = cc_label_ptr[nbh_id];
    //               while(eq_table[cursor] != cursor)
    //                 {
    //                   cursor = eq_table[cursor];
    //                 }

    //               prov = (cursor < prov)? cursor : prov;
    //             }
    //         }
      
    //   if (prov == label_boundary)
    //     {
    //       cc_label_ptr[id] = (double) grid_ijk_to_index(domain_dims * subdiv, (cell_loc + offset) * subdiv + subcell_loc) + 1.0;
    //       eq_table[cc_label_ptr[id]] = cc_label_ptr[id];
    //     }
    //   else
    //     {
    //       cc_label_ptr[id] = prov;
          
    //       for( ssize_t nk=-1 ; nk < 2 ; nk++)
    //         for( ssize_t nj=-1 ; nj < 2; nj++)
    //           for( ssize_t ni=-1 ; ni < 2; ((nk==0 && nj==0)? ni+=2 :ni++))
    //             {
    //               IJK nbh_cell_loc;
    //               IJK nbh_subcell_loc;
    //               subcell_neighbor( cell_loc, subcell_loc, subdiv, IJK{ni,nj,nk}, nbh_cell_loc, nbh_subcell_loc );
    //               size_t nbh_id = grid_ijk_to_index(grid_dims, nbh_cell_loc) * stride + grid_ijk_to_index(IJK{subdiv, subdiv, subdiv}, nbh_subcell_loc);
    //               if(is_inside(working_block, nbh_cell_loc) && density_ptr[nbh_id] > density_threshold && cc_label_ptr[nbh_id] != 0.0)
    //                 {
    //                   double cursor = cc_label_ptr[nbh_id];
    //                   while(eq_table[cursor] != cursor)
    //                     {
    //                       cursor = eq_table[cursor];
    //                     }
    //                   eq_table[cursor] = prov;
    //                 }
    //             }
    //     }
    // }

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
Computes analytical informations on a connected grid
)EOF";
    }

  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory("connected_grid_analysis", make_simple_operator< ConnectedGridAnalysis > );
  }

}
