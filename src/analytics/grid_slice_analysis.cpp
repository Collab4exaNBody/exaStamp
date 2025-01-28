#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/grid.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/grid_algorithm.h>
#include <exanb/core/simple_block_rcb.h>
#include <exanb/core/physics_constants.h>

#include <exanb/core/make_grid_variant_operator.h>

#include <exanb/core/grid_connected_components.h>

#include <string>
#include <unordered_map>
#include <utility>
#include <cstdint>
#include <cstdlib>
#include <mpi.h>
#include <omp.h>

namespace exaStamp
{
  // template<
  //   class GridT ,
  //   class = AssertGridHasFields< GridT, /*field::_id, field::_rx, field::_ry, field::_rz,*/ field::_vx, field::_vy, field::_vz >
  //   >
  template<class GridT>
  class GridSliceAnalysis : public OperatorNode
  {
    ADD_SLOT( MPI_Comm       , mpi                , INPUT        , MPI_COMM_WORLD , DocString{"MPI communicator"} );
    ADD_SLOT( Domain         , domain             , INPUT        , REQUIRED );
    ADD_SLOT( GridCellValues , grid_cell_values   , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( GridT          , grid               , INPUT        , REQUIRED );

    ADD_SLOT( ssize_t        , slice_thickness    , INPUT        , 0            , DocString{"Thickness of slices for the analysis. Thickness in subcells. Default value 0 <=> will be equal to subdiv."} );
    ADD_SLOT( IJK            , slice_direction    , INPUT        , IJK{1, 0, 0} , DocString{"Direction in which we will slice our simulation. For now it can only be {1, 0, 0}, {0, 1, 0} or {0, 0, 1}"} );
    
    ADD_SLOT( std::string    , filename           , INPUT , REQUIRED );
    ADD_SLOT( std::string    , data_separator     , INPUT , "," );
    ADD_SLOT( bool           , slice_file         , INPUT , true );
    
    // ADD_SLOT( bool                  , enable_mass             , INPUT , true );
    // ADD_SLOT( bool                  , enable_mass_center      , INPUT , false);
    ADD_SLOT( bool           , enable_velocity      , INPUT , true  );
    ADD_SLOT( bool           , enable_density       , INPUT , true  );
    ADD_SLOT( bool           , enable_temperature   , INPUT , true  );
    ADD_SLOT( bool           , enable_pressure      , INPUT , true  );
    ADD_SLOT( bool           , enable_stress_tensor , INPUT , false );
    ADD_SLOT( bool           , enable_porosity      , INPUT , false );
    
    ADD_SLOT( bool           , enable_debug         , INPUT , false );
    ADD_SLOT( bool           , grid_info            , INPUT , false );

    ADD_SLOT( bool           , fields_time_diagrams , INPUT , false );
    ADD_SLOT( long           , timestep             , INPUT , REQUIRED );
    ADD_SLOT( std::string    , simulation_name      , INPUT , "Nobody" ); // Its name is Nobody (obscur ref) if you don't give it a name.  
    ADD_SLOT( std::string    , file_extension       , INPUT , "csv" );

    ADD_SLOT( std::string    , contribution_field_name        , INPUT , REQUIRED );
    // ADD_SLOT( std::string    , density_field_name             , INPUT , REQUIRED );
    ADD_SLOT( std::string    , mass_field_name                , INPUT , REQUIRED );
    ADD_SLOT( std::string    , momentum_vector_field_name     , INPUT , REQUIRED );
    ADD_SLOT( std::string    , virial_tensor_field_name       , INPUT , REQUIRED );
    ADD_SLOT( std::string    , m_v2_field_name                , INPUT , REQUIRED );
    ADD_SLOT( std::string    , m_v2_tensor_field_name         , INPUT , REQUIRED );
    ADD_SLOT( std::string    , tr_virial_field_name           , INPUT , REQUIRED );
    
    ADD_SLOT( std::string    , sl_density_field_name          , INPUT_OUTPUT , "sl_density (g / cm^3)" );
    ADD_SLOT( std::string    , sl_velocity_field_name         , INPUT_OUTPUT , "sl_velocity (m/s)" );
    ADD_SLOT( std::string    , sl_temperature_field_name      , INPUT_OUTPUT , "sl_temperature (K)" );
    ADD_SLOT( std::string    , sl_pressure_field_name         , INPUT_OUTPUT , "sl_pressure (GPa)" );
    ADD_SLOT( std::string    , sl_stress_tensor_field_name    , INPUT_OUTPUT , "sl_stress_tensor (GPa)" );
    ADD_SLOT( std::string    , sl_porosity_field_name         , INPUT_OUTPUT , "sl_porosity (-)" );
    
    inline bool is_sink() const override final { return true; } // not a suppressable operator
    
  public:

    bool first_time = true;

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      // compile time constant indicating if grid has type field
      //using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
      //using has_vx_field_t = typename GridT::CellParticles::template HasField < field::_vx > ;
      //using has_vy_field_t = typename GridT::CellParticles::template HasField < field::_vy > ;
      //using has_vz_field_t = typename GridT::CellParticles::template HasField < field::_vz > ;
      //using has_virial_field_t = typename GridT::CellParticles::template HasField < field::_virial > ;
      // static_assert( sizeof(ssize_t) == sizeof(double) ); // the module works only in this case
      // static_assert( sizeof(long long) == sizeof(ssize_t) ); // ensures we can send ssize_t with MPI_LONG_LONG
      
      // // ---- checks on the capacity of the module to execute it self as it has been ask by the user. (Users can be dumb. Especially me...)
      // if( ! grid_cell_values->has_field("cc_label") )
      //   {
      //     if(first_time == true)
      //       {
      //         lout << std::endl <<"WARNING !" << std::endl << std::endl << "You are trying to analyse a connected grid that hasn't been connected. This module will not execute itself."
      //              << std::endl << std::endl;
              
              
      //         first_time = false;
      //       }
      //     return;
      //   }

      // if( *enable_mass_center && ! grid_cell_values->has_field(*contribution_field_name) )
      //   {
      //     if(first_time == true)
      //       {
      //         lout << std::endl <<"WARNING !" << std::endl << std::endl << "You are trying to compute the center of mass from regions connected trough connected components on a grid that hasn't mass "
      //              << "infomation. The module will skip this computation." << std::endl << std::endl;
      //       }
      //     *enable_mass_center = false;
      //   }

      // if( *enable_geometric_center && *enable_volume == false )
      //   {
      //     if(first_time == true)
      //       {
      //         lout << std::endl <<"WARNING !" << std::endl << std::endl << "You are trying to compute the geometric center from regions connected trough connected components without having enable volume computation. The volume computation is needed. "
      //              << "The module will skip this computation." << std::endl << std::endl;
      //       }
      //     *enable_geometric_center = false;
      //   }
      
      // if( *enable_velocity &&  ! grid_cell_values->has_field(*velocity_vector_field_name) )
      //   {
      //     if(first_time == true)
      //       {
      //         lout << std::endl <<"WARNING !" << std::endl << std::endl << "You are trying to compute the velocity of regions connected trough connected components on a grid that hasn't velocity vector "
      //              << "infomation. The module will skip this computation." << std::endl << std::endl;
      //       }
      //     *enable_velocity = false;
      //   }
      
      // if( *enable_density &&  ! grid_cell_values->has_field(*density_field_name) )
      //   {
      //     if(first_time == true)
      //       {
      //         lout << std::endl <<"WARNING !" << std::endl << std::endl << "You are trying to compute the density of regions connected trough connected components on a grid that hasn't mass "
      //              << "infomation. The module will skip this computation." << std::endl << std::endl;
      //       }
      //     *enable_density = false;
      //   }

      // if(*enable_volume == false && *enable_density == false && *enable_mass_center == false && *enable_velocity == false && *enable_geometric_center == false)
      //   {
      //     if(first_time == true)
      //       {
      //         lout << std::endl <<"WARNING !" << std::endl << std::endl << "There is no computation to do or that can be done in this module (connected_grid_analysis)!"
      //              << std::endl << std::endl;
      //       }
          
      //     first_time = false;
      //     return;
      //   }
      // first_time = false;
      // // ---- checks' end
      
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
      const double subcell_vol  = subcell_size * subcell_size * subcell_size;
      
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
      if(*grid_info == true && first_time == true)
        {
          if( *enable_density && ! grid_cell_values->has_field(*sl_density_field_name) )
            {
              grid_cell_values->add_field(*sl_density_field_name,subdiv,1);
            }
          if( *enable_temperature && ! grid_cell_values->has_field(*sl_temperature_field_name) )
            {
              grid_cell_values->add_field(*sl_temperature_field_name,subdiv,1);
            }
          if( *enable_velocity && ! grid_cell_values->has_field(*sl_velocity_field_name) )
            {
              grid_cell_values->add_field(*sl_velocity_field_name,subdiv,3);
            }
          if( *enable_pressure && ! grid_cell_values->has_field(*sl_pressure_field_name) )
            {
              grid_cell_values->add_field(*sl_pressure_field_name,subdiv,1);
            }
          if( *enable_stress_tensor && ! grid_cell_values->has_field(*sl_stress_tensor_field_name) )
            {
              grid_cell_values->add_field(*sl_stress_tensor_field_name,subdiv,9);
            }
          if( *enable_porosity && ! grid_cell_values->has_field(*sl_porosity_field_name) )
            {
              grid_cell_values->add_field(*sl_porosity_field_name,subdiv,1);
            }
          first_time = false;
        }
      
      int nb_field = 0;
      
      // -- Retrieve grid data field accessor :
      const auto ctb_accessor = grid_cell_values->field_data(*contribution_field_name);
      // contribution field data accessor.
      double * __restrict__ tmp_ctb_ptr = nullptr;
      // if (*enable_velocity || *enable_temperature /*|| *enable_density*/) // enable density ? <-- always needed now
        {
          tmp_ctb_ptr = grid_cell_values->field_data(*contribution_field_name).m_data_ptr;
          nb_field += 1;
        }
      const double * __restrict__ ctb_ptr = tmp_ctb_ptr;
      const size_t stride = ctb_accessor.m_stride; // for simplicity

      // // density field data accessor.
      // double * __restrict__ tmp_density_ptr = nullptr;
      // if (*enable_density)
      //   {
      //     tmp_density_ptr = grid_cell_values->field_data(*density_field_name).m_data_ptr;
      //     nb_field += 1;
      //   }
      // const double * __restrict__ density_ptr = tmp_density_ptr;

      // mass field data accessor.
      double * __restrict__ tmp_mass_ptr = nullptr;
      if (*enable_density || *enable_porosity || *enable_temperature || *enable_pressure || *enable_stress_tensor)
        {
          tmp_mass_ptr = grid_cell_values->field_data(*mass_field_name).m_data_ptr;
          nb_field += 1;
        }
      const double * __restrict__ mass_ptr = tmp_mass_ptr;
      
      // momentum field data accessor.
      double * __restrict__ tmp_momentum_ptr = nullptr;
      if (*enable_velocity || *enable_temperature || *enable_pressure || *enable_stress_tensor)
        {
          tmp_momentum_ptr = grid_cell_values->field_data(*momentum_vector_field_name).m_data_ptr;
          nb_field += 3;
        }
      const double * __restrict__ momentum_ptr = tmp_momentum_ptr;

      // m_v2 field data accessor.
      double * __restrict__ tmp_m_v2_ptr = nullptr;
      if (*enable_temperature || *enable_pressure)
        {
          tmp_m_v2_ptr = grid_cell_values->field_data(*m_v2_field_name).m_data_ptr;
          nb_field += 1;
        }
      const double * __restrict__ m_v2_ptr = tmp_m_v2_ptr;

      // tr_virial field data accessor.
      double * __restrict__ tmp_tr_virial_ptr = nullptr;
      if (*enable_pressure)
        {
          tmp_tr_virial_ptr = grid_cell_values->field_data(*tr_virial_field_name).m_data_ptr;
          nb_field += 1;
        }
      const double * __restrict__ tr_virial_ptr = tmp_tr_virial_ptr;

      // virial field data accessor.
      double * __restrict__ tmp_virial_ptr = nullptr;
      if (*enable_stress_tensor)
        {
          tmp_virial_ptr = grid_cell_values->field_data(*virial_tensor_field_name).m_data_ptr;
          nb_field += 9;
        }
      const double * __restrict__ virial_ptr = tmp_virial_ptr;

      // m_v2_tensor field data accessor.
      double * __restrict__ tmp_m_v2_tensor_ptr = nullptr;
      if (*enable_stress_tensor)
        {
          tmp_m_v2_tensor_ptr = grid_cell_values->field_data(*m_v2_tensor_field_name).m_data_ptr;
          nb_field += 9;
        }
      const double * __restrict__ m_v2_tensor_ptr = tmp_m_v2_tensor_ptr;

      if(*enable_porosity)
        nb_field += 1;
      
      // -----------------------------------------------

      if (*slice_thickness == 0)
        *slice_thickness = subdiv;
      
      // Length in subcells in the slice direction. Works only with {1, 0, 0} or {0, 1, 0} or {0, 0, 1}
      const ssize_t subcell_in_sl_dir = ((slice_direction->i * grid_dims.i + slice_direction->j * grid_dims.j + slice_direction->k * grid_dims.k) - 2*gl) * subdiv ;
      // std::cout <<"subcell_in_sl_dir : " << subcell_in_sl_dir << std::endl;

      const ssize_t local_orig = (slice_direction->i * grid_offset.i + slice_direction->j * grid_offset.j + slice_direction->k * grid_offset.k + gl) * subdiv;
      // std::cout <<"local_orig : " << local_orig << std::endl;
      const std::lldiv_t dv = lldiv( local_orig , *slice_thickness);
      // std::cout <<"dv : " << dv.quot << " " << dv.rem << std::endl;
      // Global id of the first local slice. 
      const ssize_t first_sl_id = dv.quot;
      // Thickness of the first slice.
      const ssize_t first_sl_thickness = *slice_thickness - dv.rem;

      const std::lldiv_t dv2 = lldiv(subcell_in_sl_dir - first_sl_thickness, *slice_thickness);
      // std::cout <<"dv2 : " << dv2.quot << " " << dv2.rem << std::endl;
      const ssize_t last_sl_thickness = dv2.rem;
      const ssize_t nb_mid_slices = (last_sl_thickness != 0)? dv2.quot : dv2.quot - 1;
      const ssize_t loc_nb_slices = (first_sl_id != (local_orig + subcell_in_sl_dir) / *slice_thickness) ? 2 + nb_mid_slices : 1;

      const IJK mask = IJK{1,1,1} - *slice_direction;
      const IJK subcell_in_mask = IJK{mask.i * domain_dims.i * subdiv, mask.j * domain_dims.j * subdiv, mask.k * domain_dims.k * subdiv};
 
      const double slice_vol = (subcell_in_mask.i + *slice_thickness * slice_direction->i) * (subcell_in_mask.j + *slice_thickness * slice_direction->j) * (subcell_in_mask.k + *slice_thickness * slice_direction->k) * subcell_vol;
      
      double *sl              = nullptr; 
      
      double *sl_contribution = nullptr;
      double *sl_density      = nullptr;
      double *sl_vx           = nullptr;
      double *sl_vy           = nullptr;
      double *sl_vz           = nullptr;
      double *sl_temperature  = nullptr;
      double *sl_pressure     = nullptr;
      double *sl_sigma_11     = nullptr;
      double *sl_sigma_12     = nullptr;
      double *sl_sigma_13     = nullptr;
      double *sl_sigma_21     = nullptr;
      double *sl_sigma_22     = nullptr;
      double *sl_sigma_23     = nullptr;
      double *sl_sigma_31     = nullptr;
      double *sl_sigma_32     = nullptr;
      double *sl_sigma_33     = nullptr;
      double *sl_m_v1_v1      = nullptr;
      double *sl_m_v1_v2      = nullptr;
      double *sl_m_v1_v3      = nullptr;
      double *sl_m_v2_v1      = nullptr;
      double *sl_m_v2_v2      = nullptr;
      double *sl_m_v2_v3      = nullptr;
      double *sl_m_v3_v1      = nullptr;
      double *sl_m_v3_v2      = nullptr;
      double *sl_m_v3_v3      = nullptr;
      double *sl_porosity     = nullptr;
      
      int rank=0, np = 0;
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &np); 

      int OMP_NUM_THREADS = omp_get_max_threads(); 

      // memory pool
      double *local_slices = (double *)calloc(OMP_NUM_THREADS * loc_nb_slices * nb_field, sizeof(double));

      double *contribution = nullptr;
      double *mass         = nullptr;
      double *mvx          = nullptr;
      double *mvy          = nullptr;
      double *mvz          = nullptr;
      double *m_v2         = nullptr;
      double *tr_virial    = nullptr;
      double *virial_11    = nullptr;
      double *virial_12    = nullptr;
      double *virial_13    = nullptr;
      double *virial_21    = nullptr;
      double *virial_22    = nullptr;
      double *virial_23    = nullptr;
      double *virial_31    = nullptr;
      double *virial_32    = nullptr;
      double *virial_33    = nullptr;
      double *m_v1_v1      = nullptr;
      double *m_v1_v2      = nullptr;
      double *m_v1_v3      = nullptr;
      double *m_v2_v1      = nullptr;
      double *m_v2_v2      = nullptr;
      double *m_v2_v3      = nullptr;
      double *m_v3_v1      = nullptr;
      double *m_v3_v2      = nullptr;
      double *m_v3_v3      = nullptr;
      double *void_volume  = nullptr;
      
      {
        int i = 0;
        // if(*enable_velocity || *enable_temperature || *enable_density) // <-- old test, now contribution is always needed
          {
            contribution = &local_slices[loc_nb_slices *  i   ];
            i += 1;
          }
        if(*enable_density || *enable_temperature || *enable_pressure || *enable_stress_tensor)
          {
            mass         = &local_slices[loc_nb_slices *  i   ];
            i += 1;
          }
        if(*enable_velocity || *enable_temperature || *enable_pressure || *enable_stress_tensor)
          {
            mvx          = &local_slices[loc_nb_slices *  i   ];
            mvy          = &local_slices[loc_nb_slices * (i+1)];
            mvz          = &local_slices[loc_nb_slices * (i+2)];
            i += 3;
          }
        if(*enable_temperature || *enable_pressure)
          {
            m_v2         = &local_slices[loc_nb_slices *  i   ];
            i += 1;
          }
        if(*enable_pressure)
          {
            tr_virial    = &local_slices[loc_nb_slices *  i   ];
            i += 1;
          }
        if(*enable_stress_tensor)
          {
            virial_11    = &local_slices[loc_nb_slices *  i   ];
            virial_12    = &local_slices[loc_nb_slices * (i+1)];
            virial_13    = &local_slices[loc_nb_slices * (i+2)];
            virial_21    = &local_slices[loc_nb_slices * (i+3)];
            virial_22    = &local_slices[loc_nb_slices * (i+4)];
            virial_23    = &local_slices[loc_nb_slices * (i+5)];
            virial_31    = &local_slices[loc_nb_slices * (i+6)];
            virial_32    = &local_slices[loc_nb_slices * (i+7)];
            virial_33    = &local_slices[loc_nb_slices * (i+8)];
            i += 9;
          }
        if(*enable_stress_tensor)
          {
            m_v1_v1      = &local_slices[loc_nb_slices *  i   ];
            m_v1_v2      = &local_slices[loc_nb_slices * (i+1)];
            m_v1_v3      = &local_slices[loc_nb_slices * (i+2)];
            m_v2_v1      = &local_slices[loc_nb_slices * (i+3)];
            m_v2_v2      = &local_slices[loc_nb_slices * (i+4)];
            m_v2_v3      = &local_slices[loc_nb_slices * (i+5)];
            m_v3_v1      = &local_slices[loc_nb_slices * (i+6)];
            m_v3_v2      = &local_slices[loc_nb_slices * (i+7)];
            m_v3_v3      = &local_slices[loc_nb_slices * (i+8)];
            i += 9;
          }
        if(*enable_porosity)
          {
            void_volume  = &local_slices[loc_nb_slices *  i   ];
            i += 1;
          }
        assert( i == nb_field);
      }
      
#     pragma omp parallel
      {
        int nthreads = omp_get_num_threads(); 
        int thread_id = omp_get_thread_num();
        
        assert(nthreads == OMP_NUM_THREADS); // check that there is the number of threads that we expect

        // Each thread retrieve its working block
        GridBlock working_block = simple_block_rcb(GridBlock{IJK{gl,gl,gl}, grid_dims-gl}, nthreads, thread_id);

        // GridBlock working_block = GridBlock{IJK{gl,gl,gl}, grid_dims-gl};
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

                        ssize_t pos = (slice_direction->i * ((i*subdiv)+si) + slice_direction->j * ((j*subdiv)+sj) + slice_direction->k * ((k*subdiv)+sk)) - gl*subdiv; 
                        ssize_t sl_local_index = (pos >= first_sl_thickness)? (pos - first_sl_thickness) / *slice_thickness + 1 : 0;
                        // ssize_t sl_local_index = (pos > first_sl_thickness)? (pos - first_sl_thickness) / *slice_thickness : 0;

                        ssize_t loc_id = thread_id * nb_field * loc_nb_slices + sl_local_index;
                        
                        ssize_t scindex_vec_0 = cell_index * stride + subcell_index * 3 + 0;
                        ssize_t scindex_vec_1 = cell_index * stride + subcell_index * 3 + 1;
                        ssize_t scindex_vec_2 = cell_index * stride + subcell_index * 3 + 2;
                        
                        if(*enable_velocity || *enable_temperature || *enable_density)
                          {
                            contribution [loc_id] += ctb_ptr[value_index];
                          }
                        if(*enable_density || *enable_temperature || *enable_pressure || *enable_stress_tensor)
                          {
                            // density      [loc_id] += density_ptr[value_index] * ctb_ptr[value_index];
                            // density      [loc_id] += density_ptr[value_index] * subcell_vol; // here density store the mass
                            mass         [loc_id] += mass_ptr[value_index];
                          }
                        if(*enable_velocity || *enable_temperature || *enable_pressure || *enable_stress_tensor)
                          {
                            mvx          [loc_id] += momentum_ptr[scindex_vec_0];
                            mvy          [loc_id] += momentum_ptr[scindex_vec_1];
                            mvz          [loc_id] += momentum_ptr[scindex_vec_2];
                          }
                        if(*enable_temperature || *enable_pressure)
                          {
                            m_v2         [loc_id] += m_v2_ptr[value_index];
                          }
                        if(*enable_pressure)
                          {
                            tr_virial    [loc_id] += tr_virial_ptr[value_index];
                          }
                        if(*enable_stress_tensor)
                          {
                            ssize_t scindex_tens_0 = cell_index * stride + subcell_index * 9 + 0;
                            ssize_t scindex_tens_1 = cell_index * stride + subcell_index * 9 + 1;
                            ssize_t scindex_tens_2 = cell_index * stride + subcell_index * 9 + 2;
                            ssize_t scindex_tens_3 = cell_index * stride + subcell_index * 9 + 3;
                            ssize_t scindex_tens_4 = cell_index * stride + subcell_index * 9 + 4;
                            ssize_t scindex_tens_5 = cell_index * stride + subcell_index * 9 + 5;
                            ssize_t scindex_tens_6 = cell_index * stride + subcell_index * 9 + 6;
                            ssize_t scindex_tens_7 = cell_index * stride + subcell_index * 9 + 7;
                            ssize_t scindex_tens_8 = cell_index * stride + subcell_index * 9 + 8;
                        
                            virial_11    [loc_id] += virial_ptr[scindex_tens_0];
                            virial_12    [loc_id] += virial_ptr[scindex_tens_1];
                            virial_13    [loc_id] += virial_ptr[scindex_tens_2];
                            virial_21    [loc_id] += virial_ptr[scindex_tens_3];
                            virial_22    [loc_id] += virial_ptr[scindex_tens_4];
                            virial_23    [loc_id] += virial_ptr[scindex_tens_5];
                            virial_31    [loc_id] += virial_ptr[scindex_tens_6];
                            virial_32    [loc_id] += virial_ptr[scindex_tens_7];
                            virial_33    [loc_id] += virial_ptr[scindex_tens_8];

                            m_v1_v1      [loc_id] += m_v2_tensor_ptr[scindex_tens_0];
                            m_v1_v2      [loc_id] += m_v2_tensor_ptr[scindex_tens_1];
                            m_v1_v3      [loc_id] += m_v2_tensor_ptr[scindex_tens_2];
                            m_v2_v1      [loc_id] += m_v2_tensor_ptr[scindex_tens_3];
                            m_v2_v2      [loc_id] += m_v2_tensor_ptr[scindex_tens_4];
                            m_v2_v3      [loc_id] += m_v2_tensor_ptr[scindex_tens_5];
                            m_v3_v1      [loc_id] += m_v2_tensor_ptr[scindex_tens_6];
                            m_v3_v2      [loc_id] += m_v2_tensor_ptr[scindex_tens_7];
                            m_v3_v3      [loc_id] += m_v2_tensor_ptr[scindex_tens_8];
                          }
                        if(*enable_porosity && mass_ptr[value_index] < 1e-2) // arbitrary mass threshold in Da
                          {
                            void_volume  [loc_id] += subcell_vol;
                          }
                      }
              }
#       pragma omp barrier

        // homemade reduction, we use threads to handle each fields (a field can be handled only by a unique thread)
        {
          int my_work = thread_id;
          while(my_work < nb_field)
            {
              for(ssize_t t = 1; t < nthreads; t++)
                for(ssize_t i = 0; i < loc_nb_slices; i++)
                  {
                    //          [   field         +  slice_id]                [       field              +  loop over threads + slice_id]
                    local_slices[my_work * loc_nb_slices +i] += local_slices[nb_field * loc_nb_slices * t + my_work * loc_nb_slices + i] ;
                  }
              my_work += nthreads;
            }
        }
#       pragma omp barrier
        
#       pragma omp single
        {
          int msg_elements = loc_nb_slices;
          int *msgs_size = nullptr;
          
          MPI_Request requests[1];
          MPI_Status status[1];
          
          int *displs = nullptr, *vcounts = nullptr;
          double *global_buffer = nullptr;

          if(np > 1 && rank != 0)
            {
              MPI_Igather(&msg_elements, 1, MPI_INT, msgs_size, 1, MPI_INT, 0, comm, &requests[0]);
              {
                local_slices[loc_nb_slices * nb_field] = convert(first_sl_id);
              }
              MPI_Wait(&requests[0], &status[0]);
              
              MPI_Gatherv(local_slices, nb_field * loc_nb_slices + 1, MPI_DOUBLE, global_buffer, vcounts, displs, MPI_DOUBLE, 0, comm);

              if(*grid_info == true)
                {
                  MPI_Scatterv(global_buffer, vcounts, displs, MPI_DOUBLE, local_slices, nb_field * loc_nb_slices + 1, MPI_DOUBLE, 0, comm);
                }
            }
          if(np > 1 && rank == 0)
            {
              msgs_size = new int[np];
              MPI_Igather(&msg_elements, 1, MPI_INT, msgs_size, 1, MPI_INT, 0, comm, &requests[0]);
              {
                local_slices[loc_nb_slices * nb_field] = convert(first_sl_id);
                displs = new int[np];
                vcounts = new int[np];
              }
              MPI_Wait(&requests[0], &status[0]);

              displs[0] = 0;
              vcounts[0] = nb_field * msgs_size[0] + 1;
              for(int i = 1; i < np; i++)
                {
                  displs[i] = displs[i-1] + vcounts[i-1];
                  vcounts[i] = nb_field * msgs_size[i] + 1;
                }
            
              global_buffer = new double[displs[np-1] + vcounts[np-1]];

              MPI_Gatherv(local_slices, nb_field * loc_nb_slices + 1, MPI_DOUBLE, global_buffer, vcounts, displs, MPI_DOUBLE, 0, comm);

              // TO DO : move the next part in a thread parallel region.

              const std::lldiv_t dv3 = lldiv( (slice_direction->i * domain_dims.i + slice_direction->j * domain_dims.j + slice_direction->k * domain_dims.k) * subdiv , *slice_thickness);
              
              const size_t g_nb_slices = (dv3.rem != 0 )? dv3.quot + 1 : dv3.quot;
              
              sl = (double *) calloc(nb_field * g_nb_slices, sizeof(double));

              {
                int i = 0;
                if(*enable_velocity || *enable_temperature || *enable_density)
                  {
                    sl_contribution = &sl[g_nb_slices *  i   ];
                    i += 1;
                  }
                if(*enable_density || *enable_temperature || *enable_pressure || *enable_stress_tensor)
                  {
                    sl_density      = &sl[g_nb_slices *  i   ];
                    i += 1;
                  }
                if(*enable_velocity || *enable_temperature || *enable_pressure || *enable_stress_tensor)
                  {
                    sl_vx           = &sl[g_nb_slices *  i   ];
                    sl_vy           = &sl[g_nb_slices * (i+1)];
                    sl_vz           = &sl[g_nb_slices * (i+2)];
                    i += 3;
                  }
                if(*enable_temperature || *enable_pressure)
                  {
                    sl_temperature  = &sl[g_nb_slices *  i   ];
                    i += 1;
                  }
                if(*enable_pressure)
                  {
                    sl_pressure     = &sl[g_nb_slices *  i   ];
                    i += 1;
                  }
                if(*enable_stress_tensor)
                  {
                    sl_sigma_11     = &sl[g_nb_slices *  i   ];
                    sl_sigma_12     = &sl[g_nb_slices * (i+1)];
                    sl_sigma_13     = &sl[g_nb_slices * (i+2)];
                    sl_sigma_21     = &sl[g_nb_slices * (i+3)];
                    sl_sigma_22     = &sl[g_nb_slices * (i+4)];
                    sl_sigma_23     = &sl[g_nb_slices * (i+5)];
                    sl_sigma_31     = &sl[g_nb_slices * (i+6)];
                    sl_sigma_32     = &sl[g_nb_slices * (i+7)];
                    sl_sigma_33     = &sl[g_nb_slices * (i+8)];
                    i += 9;
                  }
                if(*enable_stress_tensor)
                  {
                    sl_m_v1_v1      = &sl[g_nb_slices *  i   ];
                    sl_m_v1_v2      = &sl[g_nb_slices * (i+1)];
                    sl_m_v1_v3      = &sl[g_nb_slices * (i+2)];
                    sl_m_v2_v1      = &sl[g_nb_slices * (i+3)];
                    sl_m_v2_v2      = &sl[g_nb_slices * (i+4)];
                    sl_m_v2_v3      = &sl[g_nb_slices * (i+5)];
                    sl_m_v3_v1      = &sl[g_nb_slices * (i+6)];
                    sl_m_v3_v2      = &sl[g_nb_slices * (i+7)];
                    sl_m_v3_v3      = &sl[g_nb_slices * (i+8)];
                    i += 9;
                  }
                if(*enable_porosity)
                  {
                    sl_porosity     = &sl[g_nb_slices *  i   ];
                  }
                assert(nb_field == i);
              }

              for(int i = 0; i < np; i++)
                {
                  ssize_t slice_offset_msg = convert(global_buffer[displs[i] + msgs_size[i] * nb_field]);
                  for(int j = 0; j < msgs_size[i]; j++)
                    {
                      int k = 0; 
                      if(*enable_velocity || *enable_temperature || *enable_density)
                        {
                          sl_contribution[slice_offset_msg + j] += global_buffer[displs[i] + msgs_size[i] *  k    + j];
                          k += 1;
                        }
                      if(*enable_density || *enable_temperature || *enable_pressure || *enable_stress_tensor)
                        {
                          sl_density[slice_offset_msg + j]      += global_buffer[displs[i] + msgs_size[i] *  k    + j];
                          k += 1;
                        }
                      if(*enable_velocity || *enable_temperature || *enable_pressure || *enable_stress_tensor)
                        {
                          sl_vx[slice_offset_msg + j]           += global_buffer[displs[i] + msgs_size[i] *  k    + j];
                          sl_vy[slice_offset_msg + j]           += global_buffer[displs[i] + msgs_size[i] * (k+1) + j];
                          sl_vz[slice_offset_msg + j]           += global_buffer[displs[i] + msgs_size[i] * (k+2) + j];
                          k += 3;
                        }
                      if(*enable_temperature || *enable_pressure)
                        {
                          sl_temperature[slice_offset_msg + j]  += global_buffer[displs[i] + msgs_size[i] *  k    + j];
                          k += 1;
                        }
                      if(*enable_pressure)
                        {
                          sl_pressure[slice_offset_msg + j]     += global_buffer[displs[i] + msgs_size[i] *  k    + j];
                          k += 1;
                        }
                      if(*enable_stress_tensor)
                        {
                          sl_sigma_11[slice_offset_msg + j]     += global_buffer[displs[i] + msgs_size[i] *  k    + j];
                          sl_sigma_12[slice_offset_msg + j]     += global_buffer[displs[i] + msgs_size[i] * (k+1) + j];
                          sl_sigma_13[slice_offset_msg + j]     += global_buffer[displs[i] + msgs_size[i] * (k+2) + j];
                          sl_sigma_21[slice_offset_msg + j]     += global_buffer[displs[i] + msgs_size[i] * (k+3) + j];
                          sl_sigma_22[slice_offset_msg + j]     += global_buffer[displs[i] + msgs_size[i] * (k+4) + j];
                          sl_sigma_23[slice_offset_msg + j]     += global_buffer[displs[i] + msgs_size[i] * (k+5) + j];
                          sl_sigma_31[slice_offset_msg + j]     += global_buffer[displs[i] + msgs_size[i] * (k+6) + j];
                          sl_sigma_32[slice_offset_msg + j]     += global_buffer[displs[i] + msgs_size[i] * (k+7) + j];
                          sl_sigma_33[slice_offset_msg + j]     += global_buffer[displs[i] + msgs_size[i] * (k+8) + j];
                          k += 9;
                        }
                      if(*enable_stress_tensor)
                        {
                          sl_m_v1_v1[slice_offset_msg + j]      += global_buffer[displs[i] + msgs_size[i] *  k    + j];
                          sl_m_v1_v2[slice_offset_msg + j]      += global_buffer[displs[i] + msgs_size[i] * (k+1) + j];
                          sl_m_v1_v3[slice_offset_msg + j]      += global_buffer[displs[i] + msgs_size[i] * (k+2) + j];
                          sl_m_v2_v1[slice_offset_msg + j]      += global_buffer[displs[i] + msgs_size[i] * (k+3) + j];
                          sl_m_v2_v2[slice_offset_msg + j]      += global_buffer[displs[i] + msgs_size[i] * (k+4) + j];
                          sl_m_v2_v3[slice_offset_msg + j]      += global_buffer[displs[i] + msgs_size[i] * (k+5) + j];
                          sl_m_v3_v1[slice_offset_msg + j]      += global_buffer[displs[i] + msgs_size[i] * (k+6) + j];
                          sl_m_v3_v2[slice_offset_msg + j]      += global_buffer[displs[i] + msgs_size[i] * (k+7) + j];
                          sl_m_v3_v3[slice_offset_msg + j]      += global_buffer[displs[i] + msgs_size[i] * (k+8) + j];
                          k += 9;
                        }
                      if(*enable_porosity)
                        {
                          sl_porosity[slice_offset_msg + j]     += global_buffer[displs[i] + msgs_size[i] *  k    + j];
                          k += 1;
                        }
                      assert(k == nb_field);
                    }
                }

              // const double temp_cst = legacy_constant::atomicMass * 10000.0 / (3.0*legacy_constant::boltzmann);
              const double temp_cst = legacy_constant::atomicMass * 10000.0 / (3.0*legacy_constant::boltzmann);
              const double density_cst = 1.660539066;
              const double velocity_cst = 100.0;
              const double pressure_cst = 1.660539066 * 0.01;
              for(ssize_t i = 0; i < g_nb_slices-1; i++) // all exept the last slice
                {
                  if((*enable_velocity || *enable_temperature || *enable_pressure || *enable_stress_tensor) && sl_density[i] > 0.0)
                    {
                      sl_vx[i]          /= sl_density[i];
                      sl_vy[i]          /= sl_density[i];
                      sl_vz[i]          /= sl_density[i];
                    }
                  if(*enable_temperature || *enable_pressure)
                    {
                      sl_temperature[i] = (sl_temperature[i] - sl_density[i] * (sl_vx[i] * sl_vx[i] + sl_vy[i] * sl_vy[i] + sl_vz[i] * sl_vz[i])); // sum m_v2 - sum_i mi vi2
                    }
                  if(*enable_pressure)
                    {
                      sl_pressure[i] = (sl_pressure[i] + sl_temperature[i]) / (slice_vol * 3.0); // in ?
                    }
                  if(*enable_stress_tensor)
                    {
                      sl_sigma_11[i] += (sl_m_v1_v1[i] - sl_density[i] * sl_vx[i] * sl_vx[i]);
                      sl_sigma_12[i] += (sl_m_v1_v2[i] - sl_density[i] * sl_vx[i] * sl_vy[i]);
                      sl_sigma_13[i] += (sl_m_v1_v3[i] - sl_density[i] * sl_vx[i] * sl_vz[i]);
                      sl_sigma_21[i] += (sl_m_v2_v1[i] - sl_density[i] * sl_vy[i] * sl_vx[i]);
                      sl_sigma_22[i] += (sl_m_v2_v2[i] - sl_density[i] * sl_vy[i] * sl_vy[i]);
                      sl_sigma_23[i] += (sl_m_v2_v3[i] - sl_density[i] * sl_vy[i] * sl_vz[i]);
                      sl_sigma_31[i] += (sl_m_v3_v1[i] - sl_density[i] * sl_vz[i] * sl_vx[i]);
                      sl_sigma_32[i] += (sl_m_v3_v2[i] - sl_density[i] * sl_vz[i] * sl_vy[i]);
                      sl_sigma_33[i] += (sl_m_v3_v3[i] - sl_density[i] * sl_vz[i] * sl_vz[i]);
                    }
                  if(*enable_density)
                    {
                      sl_density[i]     /= slice_vol;
                    }
                  
                  
                  if(*enable_density)
                    {
                      sl_density[i]     *= density_cst; // switch Da/A^3 to g/cm^3
                    }
                  if(*enable_velocity)
                    {
                      sl_vx[i]          *= velocity_cst; // switch A/ps to m/s
                      sl_vy[i]          *= velocity_cst;
                      sl_vz[i]          *= velocity_cst;
                    }
                  if(*enable_temperature && sl_contribution[i] > 0.0)
                    {
                      sl_temperature[i] *= (temp_cst / sl_contribution[i]); // switch ? to K
                    }
                  
                  if(*enable_pressure)
                    {
                      sl_pressure[i] *= pressure_cst; // switch ? to GPa
                    }
                  if(*enable_stress_tensor)
                    {
                      sl_sigma_11[i] *= -pressure_cst / slice_vol; // switch ? to GPa
                      sl_sigma_12[i] *= -pressure_cst / slice_vol; 
                      sl_sigma_13[i] *= -pressure_cst / slice_vol; 
                      sl_sigma_21[i] *= -pressure_cst / slice_vol;
                      sl_sigma_22[i] *= -pressure_cst / slice_vol;
                      sl_sigma_23[i] *= -pressure_cst / slice_vol;
                      sl_sigma_31[i] *= -pressure_cst / slice_vol;
                      sl_sigma_32[i] *= -pressure_cst / slice_vol;
                      sl_sigma_33[i] *= -pressure_cst / slice_vol;
                    }
                  if(*enable_porosity)
                    {
                      sl_porosity[i] /= slice_vol; // no unit
                    }
                }
              { // do the last slice
                const ssize_t g_last_thickness = (domain_dims.i * subdiv * slice_direction->i + domain_dims.j * subdiv * slice_direction->j + domain_dims.k * subdiv * slice_direction->k) % *slice_thickness;
                const double g_last_slice_vol = (g_last_thickness != 0) ? (subcell_in_mask.i + g_last_thickness * slice_direction->i) * (subcell_in_mask.j + g_last_thickness * slice_direction->j) * (subcell_in_mask.k + g_last_thickness * slice_direction->k) * subcell_vol : slice_vol;

                if((*enable_velocity || *enable_temperature || *enable_pressure || *enable_stress_tensor) && sl_density[g_nb_slices-1] > 0.0)
                  {
                    sl_vx[g_nb_slices-1]          /= sl_density[g_nb_slices-1];
                    sl_vy[g_nb_slices-1]          /= sl_density[g_nb_slices-1];
                    sl_vz[g_nb_slices-1]          /= sl_density[g_nb_slices-1];
                  }
                if(*enable_temperature || *enable_pressure)
                  {
                    sl_temperature[g_nb_slices-1] = (sl_temperature[g_nb_slices-1] - sl_density[g_nb_slices-1] * (sl_vx[g_nb_slices-1] * sl_vx[g_nb_slices-1] + sl_vy[g_nb_slices-1] * sl_vy[g_nb_slices-1] + sl_vz[g_nb_slices-1] * sl_vz[g_nb_slices-1])); // sum m_v2 - sum mi vi2
                  }
                if(*enable_pressure)
                  {
                    sl_pressure[g_nb_slices-1] = (sl_pressure[g_nb_slices-1] + sl_temperature[g_nb_slices-1]) / (g_last_slice_vol * 3.0); // in ?
                  }
                if(*enable_density)
                  {
                    sl_density[g_nb_slices-1]     /= g_last_slice_vol;
                  }
                if(*enable_stress_tensor)
                  {
                    sl_sigma_11[g_nb_slices-1] += sl_m_v1_v1[g_nb_slices-1];
                    sl_sigma_12[g_nb_slices-1] += sl_m_v1_v2[g_nb_slices-1];
                    sl_sigma_13[g_nb_slices-1] += sl_m_v1_v3[g_nb_slices-1];
                    sl_sigma_21[g_nb_slices-1] += sl_m_v2_v1[g_nb_slices-1];
                    sl_sigma_22[g_nb_slices-1] += sl_m_v2_v2[g_nb_slices-1];
                    sl_sigma_23[g_nb_slices-1] += sl_m_v2_v3[g_nb_slices-1];
                    sl_sigma_31[g_nb_slices-1] += sl_m_v3_v1[g_nb_slices-1];
                    sl_sigma_32[g_nb_slices-1] += sl_m_v3_v2[g_nb_slices-1];
                    sl_sigma_33[g_nb_slices-1] += sl_m_v3_v3[g_nb_slices-1];
                  }
                  
                if(*enable_density)
                  {
                    sl_density[g_nb_slices-1]     *= density_cst; // switch Da/A^3 to g/cm^3
                  }
                if(*enable_velocity)
                  {
                    sl_vx[g_nb_slices-1]          *= velocity_cst; // switch A/ps to m/s
                    sl_vy[g_nb_slices-1]          *= velocity_cst;
                    sl_vz[g_nb_slices-1]          *= velocity_cst;
                  }
                if(*enable_temperature && sl_contribution[g_nb_slices-1] > 0.0)
                  {
                    sl_temperature[g_nb_slices-1] *= (temp_cst / sl_contribution[g_nb_slices-1]); // switch ? to K
                  }
                if(*enable_pressure)
                  {
                    sl_pressure[g_nb_slices-1] *= pressure_cst; // switch ? to GPa
                  }
                if(*enable_stress_tensor)
                    {
                      sl_sigma_11[g_nb_slices-1] *= -pressure_cst / g_last_slice_vol; // switch ? to GPa
                      sl_sigma_12[g_nb_slices-1] *= -pressure_cst / g_last_slice_vol; 
                      sl_sigma_13[g_nb_slices-1] *= -pressure_cst / g_last_slice_vol; 
                      sl_sigma_21[g_nb_slices-1] *= -pressure_cst / g_last_slice_vol;
                      sl_sigma_22[g_nb_slices-1] *= -pressure_cst / g_last_slice_vol;
                      sl_sigma_23[g_nb_slices-1] *= -pressure_cst / g_last_slice_vol;
                      sl_sigma_31[g_nb_slices-1] *= -pressure_cst / g_last_slice_vol;
                      sl_sigma_32[g_nb_slices-1] *= -pressure_cst / g_last_slice_vol;
                      sl_sigma_33[g_nb_slices-1] *= -pressure_cst / g_last_slice_vol;
                    }
                if(*enable_porosity)
                    {
                      sl_porosity[g_nb_slices-1] /= g_last_slice_vol; // no unit
                    }
              }
              
              if(*grid_info == true)
                {
                  for(int i = 0; i < np; i++)
                    {
                      ssize_t slice_offset_msg = convert(global_buffer[displs[i] + msgs_size[i] * nb_field]);
                      for(int j = 0; j < msgs_size[i]; j++)
                        {
                          int k = 0; 
                          if(*enable_velocity || *enable_temperature || *enable_density)
                            {
                              global_buffer[displs[i] + msgs_size[i] *  k    + j] = sl_contribution[slice_offset_msg + j];
                              k += 1;
                            }
                          if(*enable_density)
                            {
                              global_buffer[displs[i] + msgs_size[i] *  k    + j] = sl_density[slice_offset_msg + j];
                              k += 1;
                            }
                          if(*enable_velocity)
                            {
                              global_buffer[displs[i] + msgs_size[i] *  k    + j] = sl_vx[slice_offset_msg + j];
                              global_buffer[displs[i] + msgs_size[i] * (k+1) + j] = sl_vy[slice_offset_msg + j];
                              global_buffer[displs[i] + msgs_size[i] * (k+2) + j] = sl_vz[slice_offset_msg + j];
                              k += 3;
                            }
                          if(*enable_temperature)
                            {
                              global_buffer[displs[i] + msgs_size[i] *  k    + j] = sl_temperature[slice_offset_msg + j];
                              k += 1;
                            }
                          if(*enable_pressure)
                            {
                              global_buffer[displs[i] + msgs_size[i] *  k    + j] = sl_pressure[slice_offset_msg + j];
                              k += 1;
                            }
                          if(*enable_stress_tensor)
                            {
                              global_buffer[displs[i] + msgs_size[i] *  k    + j] = sl_sigma_11[slice_offset_msg + j];
                              global_buffer[displs[i] + msgs_size[i] * (k+1) + j] = sl_sigma_12[slice_offset_msg + j];
                              global_buffer[displs[i] + msgs_size[i] * (k+2) + j] = sl_sigma_13[slice_offset_msg + j];
                              global_buffer[displs[i] + msgs_size[i] * (k+3) + j] = sl_sigma_21[slice_offset_msg + j];
                              global_buffer[displs[i] + msgs_size[i] * (k+4) + j] = sl_sigma_22[slice_offset_msg + j];
                              global_buffer[displs[i] + msgs_size[i] * (k+5) + j] = sl_sigma_23[slice_offset_msg + j];
                              global_buffer[displs[i] + msgs_size[i] * (k+6) + j] = sl_sigma_31[slice_offset_msg + j];
                              global_buffer[displs[i] + msgs_size[i] * (k+7) + j] = sl_sigma_32[slice_offset_msg + j];
                              global_buffer[displs[i] + msgs_size[i] * (k+8) + j] = sl_sigma_33[slice_offset_msg + j];
                              k += 9;
                            }
                          if(*enable_porosity)
                            {
                              global_buffer[displs[i] + msgs_size[i] *  k    + j] = sl_porosity[slice_offset_msg + j];
                              k += 1;
                            }
                        }
                    }
                  MPI_Scatterv(global_buffer, vcounts, displs, MPI_DOUBLE, local_slices, nb_field * loc_nb_slices + 1, MPI_DOUBLE, 0, comm);
                }
              
              delete msgs_size;
              delete displs;
              delete vcounts;

              if(*slice_file)
                {
                  FILE *output = fopen((*filename).c_str(), "w");
                  fprintf(output,"pos(A)");
                  if(*enable_density)
                    {
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"density(g/cm3)");
                    }
                  if(*enable_velocity)
                    {
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"vx(m/s)");
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"vy(m/s)");
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"vz(m/s)");
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
                  if(*enable_stress_tensor)
                    {
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"sigma_xx(GPa)");
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"sigma_xy(GPa)");
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"sigma_xz(GPa)");
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"sigma_yx(GPa)");
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"sigma_yy(GPa)");
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"sigma_yz(GPa)");
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"sigma_zx(GPa)");
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"sigma_zy(GPa)");
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"sigma_zz(GPa)");
                    }
                  if(*enable_porosity)
                    {
                      fprintf(output,(*data_separator).c_str());
                      fprintf(output,"porosity(-)");
                    }
                  fprintf(output,"\n");
          
                  for(size_t i = 0; i < g_nb_slices; i++)
                    {
                      fprintf(output,"%.16lf", i * (*slice_thickness) * subcell_size);
                      if(*enable_density)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_density[i]);
                        }
                      if(*enable_velocity)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_vx[i]);
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_vy[i]);
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_vz[i]);
                        }
                      if(*enable_temperature)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_temperature[i]);
                        }
                      if(*enable_pressure)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_pressure[i]);
                        }
                      if(*enable_stress_tensor)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_11[i]);
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_12[i]);
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_13[i]);
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_21[i]);
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_22[i]);
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_23[i]);
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_31[i]);
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_32[i]);
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_33[i]);
                        }
                      if(*enable_porosity)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_porosity[i]);
                        }
                      fprintf(output,"\n");
                    }
                  fclose(output);
                }
    
              if(*fields_time_diagrams == true)
                {
                  if(*enable_density)
                    {
                      std::string field_filename = *simulation_name + ".density-time_diagram." + *file_extension;
                      FILE *output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_density[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);
                    }
                  if(*enable_velocity)
                    {
                      std::string field_filename = *simulation_name + ".x_velocity-time_diagram." + *file_extension;
                      FILE *output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_vx[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);
                  
                      field_filename = *simulation_name + ".y_velocity-time_diagram." + *file_extension;
                      output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_vy[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);

                      field_filename = *simulation_name + ".z_velocity-time_diagram." + *file_extension;
                      output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_vz[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);
                    }
                  if(*enable_temperature)
                    {
                      std::string field_filename = *simulation_name + ".temperature-time_diagram." + *file_extension;
                      FILE *output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_temperature[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);
                    }
                  if(*enable_pressure)
                    {
                      std::string field_filename = *simulation_name + ".pressure-time_diagram." + *file_extension;
                      FILE *output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_pressure[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);
                    }
                  if(*enable_stress_tensor)
                    {
                      std::string field_filename = *simulation_name + ".sigma_xx-time_diagram." + *file_extension;
                      FILE *output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_11[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);
                  
                      field_filename = *simulation_name + ".sigma_xy-time_diagram." + *file_extension;
                      output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_12[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);

                      field_filename = *simulation_name + ".sigma_xz-time_diagram." + *file_extension;
                      output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_13[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);

                      field_filename = *simulation_name + ".sigma_yx-time_diagram." + *file_extension;
                      output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_21[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);
                  
                      field_filename = *simulation_name + ".sigma_yy-time_diagram." + *file_extension;
                      output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_22[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);

                      field_filename = *simulation_name + ".sigma_yz-time_diagram." + *file_extension;
                      output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_23[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);

                      field_filename = *simulation_name + ".sigma_zx-time_diagram." + *file_extension;
                      output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_31[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);
                  
                      field_filename = *simulation_name + ".sigma_zy-time_diagram." + *file_extension;
                      output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_32[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);

                      field_filename = *simulation_name + ".sigma_zz-time_diagram." + *file_extension;
                      output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_sigma_33[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);
                    }
                  if(*enable_porosity)
                    {
                      std::string field_filename = *simulation_name + ".porosity-time_diagram." + *file_extension;
                      FILE *output = fopen(field_filename.c_str(), "a");
                      fprintf(output,"%ld",*timestep);
                      for(int i=0; i < g_nb_slices; i++)
                        {
                          fprintf(output,(*data_separator).c_str());
                          fprintf(output,"%.16lf",sl_porosity[i]);
                        }
                      fprintf(output,"\n");
                      fclose(output);
                    }
                }
            }
          free(sl);
        }// end pragma omp single

        if(*grid_info == true)
          {
            // retrieve field data accessor for grid info.
            double * __restrict__ sl_density_ptr = nullptr;
            if (*enable_density)
              {
                sl_density_ptr = grid_cell_values->field_data(*sl_density_field_name).m_data_ptr;
              }
            double * __restrict__ sl_velocity_ptr = nullptr;
            double *vx, *vy, *vz = nullptr;
            if (*enable_velocity)
              {
                sl_velocity_ptr = grid_cell_values->field_data(*sl_velocity_field_name).m_data_ptr;
                vx = mvx; // same as mv adress but know it is v so we change the name. 
                vy = mvy;
                vz = mvz;
              }
            double * __restrict__ sl_temperature_ptr = nullptr;
            double *temperature = nullptr;
            if (*enable_temperature)
              {
                sl_temperature_ptr = grid_cell_values->field_data(*sl_temperature_field_name).m_data_ptr;
                temperature = m_v2; // same as m_v2 adress but know it is the temperature so we change the name. 
              }
            double * __restrict__ sl_pressure_ptr = nullptr;
            double *pressure = nullptr;
            if (*enable_pressure)
              {
                sl_pressure_ptr = grid_cell_values->field_data(*sl_pressure_field_name).m_data_ptr;
                pressure = tr_virial; // same as tr_virial adress but know it is the pressure so we change the name. 
              }
            double * __restrict__ sl_stress_tensor_ptr = nullptr;
            double *sigma_xx, *sigma_xy, *sigma_xz, *sigma_yx, *sigma_yy, *sigma_yz, *sigma_zx, *sigma_zy, *sigma_zz = nullptr;
            if (*enable_stress_tensor)
              {
                sl_stress_tensor_ptr = grid_cell_values->field_data(*sl_stress_tensor_field_name).m_data_ptr;
                sigma_xx = virial_11; // same as virial adress but know it is the stress tensor so we change the name. 
                sigma_xy = virial_12;
                sigma_xz = virial_13;
                sigma_yx = virial_21;
                sigma_yy = virial_22;
                sigma_yz = virial_23;
                sigma_zx = virial_31;
                sigma_zy = virial_32;
                sigma_zz = virial_33;
              }
            double * __restrict__ sl_porosity_ptr = nullptr;
            double *porosity = nullptr;
            if (*enable_porosity)
              {
                sl_porosity_ptr = grid_cell_values->field_data(*sl_porosity_field_name).m_data_ptr;
                porosity = void_volume; // same as void_volume adress but know it is the porosity so we change the name. 
              }
            
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

                            ssize_t pos = (slice_direction->i * ((i*subdiv)+si) + slice_direction->j * ((j*subdiv)+sj) + slice_direction->k * ((k*subdiv)+sk)) - gl*subdiv; 
                            ssize_t sl_local_index = (pos >= first_sl_thickness)? (pos - first_sl_thickness) / *slice_thickness + 1 : 0;
                            
                            ssize_t scindex_vec_0 = cell_index * stride + subcell_index * 3 + 0;
                            ssize_t scindex_vec_1 = cell_index * stride + subcell_index * 3 + 1;
                            ssize_t scindex_vec_2 = cell_index * stride + subcell_index * 3 + 2;
                            
                            if(*enable_density)
                              {
                                sl_density_ptr[value_index] = mass[sl_local_index];
                              }
                            if(*enable_velocity)
                              {
                                sl_velocity_ptr[scindex_vec_0] = vx[sl_local_index];
                                sl_velocity_ptr[scindex_vec_1] = vy[sl_local_index];
                                sl_velocity_ptr[scindex_vec_2] = vz[sl_local_index];
                              }
                            if(*enable_temperature)
                              {
                                sl_temperature_ptr[value_index] = temperature[sl_local_index];
                              }
                            if(*enable_pressure)
                              {
                                sl_pressure_ptr[value_index] = pressure[sl_local_index];
                              }
                            if(*enable_stress_tensor)
                              {
                                ssize_t scindex_tens_0 = cell_index * stride + subcell_index * 9 + 0;
                                ssize_t scindex_tens_1 = cell_index * stride + subcell_index * 9 + 1;
                                ssize_t scindex_tens_2 = cell_index * stride + subcell_index * 9 + 2;
                                ssize_t scindex_tens_3 = cell_index * stride + subcell_index * 9 + 3;
                                ssize_t scindex_tens_4 = cell_index * stride + subcell_index * 9 + 4;
                                ssize_t scindex_tens_5 = cell_index * stride + subcell_index * 9 + 5;
                                ssize_t scindex_tens_6 = cell_index * stride + subcell_index * 9 + 6;
                                ssize_t scindex_tens_7 = cell_index * stride + subcell_index * 9 + 7;
                                ssize_t scindex_tens_8 = cell_index * stride + subcell_index * 9 + 8;
                        
                                sl_stress_tensor_ptr[scindex_tens_0] = sigma_xx[sl_local_index];
                                sl_stress_tensor_ptr[scindex_tens_1] = sigma_xy[sl_local_index];
                                sl_stress_tensor_ptr[scindex_tens_2] = sigma_xz[sl_local_index];
                                sl_stress_tensor_ptr[scindex_tens_3] = sigma_yx[sl_local_index];
                                sl_stress_tensor_ptr[scindex_tens_4] = sigma_yy[sl_local_index];
                                sl_stress_tensor_ptr[scindex_tens_5] = sigma_yz[sl_local_index];
                                sl_stress_tensor_ptr[scindex_tens_6] = sigma_zx[sl_local_index];
                                sl_stress_tensor_ptr[scindex_tens_7] = sigma_zy[sl_local_index];
                                sl_stress_tensor_ptr[scindex_tens_8] = sigma_zz[sl_local_index];
                              }
                            if(*enable_porosity)
                              {
                                sl_porosity_ptr[value_index] = porosity[sl_local_index];
                              }
                          }
                  }
          }
        
      } // end omp parallel
      
      free(local_slices);
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
  ONIKA_AUTORUN_INIT(grid_slice_analysis)
  {
    // OperatorNodeFactory::instance()->register_factory("grid_slice_analysis", make_simple_operator< GridSliceAnalysis > );
    OperatorNodeFactory::instance()->register_factory("grid_slice_analysis", make_grid_variant_operator< GridSliceAnalysis > );
  }

}
