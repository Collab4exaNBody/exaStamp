#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/fields.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/memory/allocator.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/physics_constants.h>

#include <exaStamp/mechanical/cell_particles_local_mechanical_metrics.h>

#include <mpi.h>
#include <iomanip>
#include <memory>

namespace exaStamp
{
  using namespace exanb; 
  
  // get particle virial tensor. assume the virial is null if particle hasn't virial field
  template<bool> static Mat3d get_particle_virial(const Mat3d* __restrict__, size_t);
  template<> inline Mat3d get_particle_virial<false>(const Mat3d* __restrict__ virials, size_t p_i) { return Mat3d(); }
  template<> inline Mat3d get_particle_virial<true>(const Mat3d* __restrict__ virials, size_t p_i) { return virials[p_i]; }
  
  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_ep, field::_ax, field::_ay, field::_az, field::_vx, field::_vy, field::_vz >
    >
  class GridReduction : public OperatorNode
  {
    // // compile time constant indicating if grid has type field
    // using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    // static constexpr bool has_type_field = has_type_field_t::value;

    // compile time constant indicating if grid has type virial
    using has_virial_field_t = typename GridT::CellParticles::template HasField < field::_virial > ;
    static constexpr bool has_virial_field = has_virial_field_t::value;
    
    ADD_SLOT( MPI_Comm       , mpi                        , INPUT , MPI_COMM_WORLD );
    
    ADD_SLOT( GridT          , grid                       , INPUT , REQUIRED );
    ADD_SLOT( Domain         , domain                     , INPUT , REQUIRED );
    // ADD_SLOT( ParticleSpecies, species                    , INPUT , REQUIRED );
    // ADD_SLOT( double         , splat_size                 , INPUT , REQUIRED );
    ADD_SLOT( long           , grid_subdiv                , INPUT , REQUIRED );
    ADD_SLOT( GridCellValues , grid_cell_values           , INPUT_OUTPUT );

    ADD_SLOT( bool           , enable_density             , INPUT , true );
    ADD_SLOT( std::string    , density_field_name         , INPUT_OUTPUT , "density (g/cm^3)" );
    ADD_SLOT( std::string    , mass_field_name            , INPUT , "mass (Da)" );
    
    ADD_SLOT( bool           , enable_velocity_x          , INPUT , false );
    ADD_SLOT( std::string    , velocity_x_field_name      , INPUT_OUTPUT , "velocity_x (m/s)" );
    ADD_SLOT( std::string    , momentum_x_field_name      , INPUT , "momentum_x (A/ps)" );

    ADD_SLOT( bool           , enable_velocity_y          , INPUT , false );
    ADD_SLOT( std::string    , velocity_y_field_name      , INPUT_OUTPUT , "velocity_y (m/s)" );
    ADD_SLOT( std::string    , momentum_y_field_name      , INPUT , "momentum_y (A/ps)" );    

    ADD_SLOT( bool           , enable_velocity_z          , INPUT , false );
    ADD_SLOT( std::string    , velocity_z_field_name      , INPUT_OUTPUT , "velocity_z (m/s)" );
    ADD_SLOT( std::string    , momentum_z_field_name      , INPUT , "momentum_z (A/ps)" );    

    ADD_SLOT( bool           , enable_velocity_vector     , INPUT , true );
    ADD_SLOT( std::string    , velocity_vector_field_name , INPUT_OUTPUT , "velocity_vector (m/s)" );
    ADD_SLOT( std::string    , momentum_vector_field_name , INPUT , "momentum_vector (A/ps)" );

    ADD_SLOT( bool           , enable_virial_tensor       , INPUT , false );
    ADD_SLOT( std::string    , virial_tensor_field_name   , INPUT , "virial_tensor (Da * A^2 / ps^2)?" );

    ADD_SLOT( std::string    , contribution_field_name    , INPUT , "contribution" );

    // -----------------------------------------------
    // ADD_SLOT( GridParticleLocalMechanicalMetrics, local_mechanical_data , INPUT );

    ADD_SLOT( bool           , enable_deformation_gradient_tensor     , INPUT , false);
    ADD_SLOT( std::string    , deformation_gradient_tensor_field_name , INPUT , "deformation_gradient_tensor (?)" );
    
  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      // compile time constant indicating if grid has type field
      using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
      using has_vx_field_t = typename GridT::CellParticles::template HasField < field::_vx > ;
      using has_vy_field_t = typename GridT::CellParticles::template HasField < field::_vy > ;
      using has_vz_field_t = typename GridT::CellParticles::template HasField < field::_vz > ;

      static constexpr has_type_field_t has_type_field{};
      static constexpr std::integral_constant<bool, has_vx_field_t::value && has_vy_field_t::value && has_vz_field_t::value> has_velocity{};
      
      // GridParticleLocalMechanicalMetrics& local_mechanical_data = *(this->local_mechanical_data);

      const ssize_t n_cells = grid->number_of_cells();
      if( n_cells == 0 )
        {
          return;
        }

      const double cell_size = domain->cell_size();
      const ssize_t subdiv = *grid_subdiv;
      const ssize_t n_subcells = subdiv * subdiv * subdiv;
      const double subcell_size = cell_size / subdiv;
      const double subcell_volume = subcell_size * subcell_size * subcell_size;

      // -------------------------------------------------------------------------------------------------------------

      // create cell data fields if needed
      if( *enable_density && ! grid_cell_values->has_field(*density_field_name))
        {
          grid_cell_values->add_field(*density_field_name,subdiv,1);
        }
      if( *enable_velocity_x && ! grid_cell_values->has_field(*velocity_x_field_name) )
        {
          grid_cell_values->add_field(*velocity_x_field_name,subdiv,1);
        }
      if( *enable_velocity_y && ! grid_cell_values->has_field(*velocity_y_field_name) )
        {
          grid_cell_values->add_field(*velocity_y_field_name,subdiv,1);
        }
      if( *enable_velocity_z && ! grid_cell_values->has_field(*velocity_z_field_name) )
        {
          grid_cell_values->add_field(*velocity_z_field_name,subdiv,1);
        }
      if( *enable_velocity_vector && ! grid_cell_values->has_field(*velocity_vector_field_name) )
        {
          grid_cell_values->add_field(*velocity_vector_field_name,subdiv,3);
        } 
      
      // -------------------------------------------------------------------------------------------------------------

      const bool need_contrib = (*enable_density || *enable_virial_tensor || *enable_deformation_gradient_tensor );
      const bool look_for_velocity = (*enable_velocity_x || *enable_velocity_y || *enable_velocity_z || *enable_velocity_vector);
      const bool need_mass = (*enable_density || look_for_velocity );
      
      // -------------------------------------------------------------------------------------------------------------
      
      // look for the grid fields
      const bool has_contrib = grid_cell_values->has_field(*contribution_field_name);
      const bool has_mass = grid_cell_values->has_field(*mass_field_name);
      const bool has_momentum_x = grid_cell_values->has_field(*momentum_x_field_name);
      const bool has_momentum_y = grid_cell_values->has_field(*momentum_y_field_name);
      const bool has_momentum_z = grid_cell_values->has_field(*momentum_z_field_name);
      const bool has_momentum_vector = grid_cell_values->has_field(*momentum_vector_field_name);
      const bool has_virial_tensor = grid_cell_values->has_field(*virial_tensor_field_name) ;
      const bool has_deformation_gradient_tensor = grid_cell_values->has_field(*deformation_gradient_tensor_field_name);
      
      
      // retreive contribution field data accessor.
      double * __restrict__ ctb_ptr = nullptr;
      size_t ctb_stride = 0;
      if( has_contrib && need_contrib)
        {
          auto cell_contribution_accessor = grid_cell_values->field_data(*contribution_field_name);
          ctb_ptr = cell_contribution_accessor.m_data_ptr;
          ctb_stride = cell_contribution_accessor.m_stride;
        }
      else if( !(has_contrib) && need_contrib )
        {
          return;
        }

      // retreive mass field data accessor.
      double * __restrict__ mass_ptr = nullptr;
      size_t mass_stride = ctb_stride;
      if( has_mass && need_mass)
        {
          auto cell_mass_accessor = grid_cell_values->field_data(*mass_field_name);
          mass_ptr = cell_mass_accessor.m_data_ptr;
        }
      else if( !(has_mass) && need_mass)
        {
          return;
        }

      double * __restrict__ density_ptr = nullptr;
      if( *enable_density )
        {
          auto cell_density_accessor = grid_cell_values->field_data(*density_field_name);
          density_ptr = cell_density_accessor.m_data_ptr;
        }
      
      // retreive momentum_x field data accessor.
      double * __restrict__ momentum_x_ptr = nullptr;
      size_t momentum_x_stride = ctb_stride;
      if( has_momentum_x && *enable_velocity_x)
        {
          auto cell_momentum_x_accessor = grid_cell_values->field_data(*momentum_x_field_name);
          momentum_x_ptr = cell_momentum_x_accessor.m_data_ptr;
        }
      else if( !(has_momentum_x) && *enable_velocity_x)
        {
          return;
        }

      double * __restrict__ velocity_x_ptr = nullptr;
      if( *enable_velocity_x )
        {
          auto cell_velocity_x_accessor = grid_cell_values->field_data(*velocity_x_field_name);
          velocity_x_ptr = cell_velocity_x_accessor.m_data_ptr;
        }
      
      // retreive momentum_y field data accessor.
      double * __restrict__ momentum_y_ptr = nullptr;
      size_t momentum_y_stride = ctb_stride;
      if( has_momentum_y && *enable_velocity_y)
        {
          auto cell_momentum_y_accessor = grid_cell_values->field_data(*momentum_y_field_name);
          momentum_y_ptr = cell_momentum_y_accessor.m_data_ptr;
        }
      else if( !(has_momentum_y) && *enable_velocity_y)
        {
          return;
        }

      double * __restrict__ velocity_y_ptr = nullptr;
      if( *enable_velocity_y )
        {
          auto cell_velocity_y_accessor = grid_cell_values->field_data(*velocity_y_field_name);
          velocity_y_ptr = cell_velocity_y_accessor.m_data_ptr;
        }
      
      // retreive momentum_z field data accessor.
      double * __restrict__ momentum_z_ptr = nullptr;
      size_t momentum_z_stride = ctb_stride;
      if( has_momentum_z && *enable_velocity_z)
        {
          auto cell_momentum_z_accessor = grid_cell_values->field_data(*momentum_z_field_name);
          momentum_z_ptr = cell_momentum_z_accessor.m_data_ptr;
        }
      else if( !(has_momentum_z) && *enable_velocity_z)
        {
          return;
        }

      double * __restrict__ velocity_z_ptr = nullptr;
      if( *enable_velocity_z )
        {
          auto cell_velocity_z_accessor = grid_cell_values->field_data(*velocity_z_field_name);
          velocity_z_ptr = cell_velocity_z_accessor.m_data_ptr;
        }
      
      // retreive momentum vector field data accessor. create data field if needed
      double * __restrict__ momentum_vector_ptr = nullptr;
      size_t momentum_vector_stride = ctb_stride;
      if( has_momentum_vector && *enable_velocity_vector)
        {
          auto cell_momentum_vector_accessor = grid_cell_values->field_data(*momentum_vector_field_name);
          momentum_vector_ptr = cell_momentum_vector_accessor.m_data_ptr;
        }
      else if( !(has_momentum_vector) && *enable_velocity_vector)
        {
          return;
        }

      double * __restrict__ velocity_vector_ptr = nullptr;
      if( *enable_velocity_vector )
        {
          auto cell_velocity_vector_accessor = grid_cell_values->field_data(*velocity_vector_field_name);
          velocity_vector_ptr = cell_velocity_vector_accessor.m_data_ptr;
        }

      // retreive virial tensor field data accessor. create data field if needed
      double * __restrict__ virial_tensor_ptr = nullptr;
      size_t virial_tensor_stride = ctb_stride;
      if( has_virial_tensor && *enable_virial_tensor )
        {
          auto cell_virial_tensor_accessor = grid_cell_values->field_data(*virial_tensor_field_name);
          virial_tensor_ptr = cell_virial_tensor_accessor.m_data_ptr;
        }
      else if( !(has_virial_tensor) && *enable_virial_tensor )
        {
          return;
        }
      
      // retreive deformation gradient tensor field data accessor. create data field if needed
      double * __restrict__ deformation_gradient_tensor_ptr = nullptr;
      size_t deformation_gradient_tensor_stride = ctb_stride;
      if( has_deformation_gradient_tensor && *enable_deformation_gradient_tensor )
        {
          auto cell_deformation_gradient_tensor_accessor = grid_cell_values->field_data(*deformation_gradient_tensor_field_name);
          deformation_gradient_tensor_ptr = cell_deformation_gradient_tensor_accessor.m_data_ptr;
        }
      else if( has_deformation_gradient_tensor && *enable_deformation_gradient_tensor )
        {
          return;
        }
      
      // // particle splatting size
      // const double sp_size = *splat_size;
      // if( sp_size/2 > subcell_size )
      //   {
      //     lerr << "in " << pathname() << std::endl
      //          << "splat_size = "<<sp_size << " is 2 times larger than subcell size = "<<subcell_size<< std::endl
      //          << "Choose smaller splat_size or coarser grid_subdiv value" << std::endl;
      //     std::abort();
      //   }
  
      // store specie masses in temporary array
      // size_t nSpecies = species->size();
      // if( nSpecies != 1 && !has_type_field )
      //   {
      //     lerr << pathname() << std::endl;
      //     lerr << "no type information, can't retreive masses" << std::endl;
      //     std::abort();
      //   }
      // double masses[nSpecies];
      // for(size_t i=0;i<nSpecies;i++) { masses[i] = species->at(i).m_mass; }

      auto cells = grid->cells();
      const IJK dims = grid->dimension();
      const ssize_t gl = grid->ghost_layers();      
      
      // ldbg << "ghost_layers="<<gl<<", cell_size="<<cell_size<<", splat_size="<<sp_size<<", subdiv="<<subdiv<<", subcell_size="<<subcell_size<<std::endl;

      const double temp_cst = legacy_constant::atomicMass * 10000.0 / (3.0*legacy_constant::boltzmann);
      const double density_cst = 1.660539066;
      const double velocity_cst = 100.0;
      const double pressure_cst = 1.660539066 * 0.01;
      
      // computes per cell Mass, per cell Mass*Velocity^2, per cell x-component of Momentum and per cell x-component or Virial
#     pragma omp parallel
      {        
        _Pragma(USTAMP_STR(omp for collapse(2) schedule(dynamic)))
          for(ssize_t i = 0; i < dims.i*dims.j*dims.k; ++i)
            for(ssize_t j = 0; j < subdiv*subdiv*subdiv; ++j)
              {
                size_t index = i * ctb_stride + j;

                size_t scindex_vec_0 = i * ctb_stride + j * 3 + 0;
                size_t scindex_vec_1 = i * ctb_stride + j * 3 + 1;
                size_t scindex_vec_2 = i * ctb_stride + j * 3 + 2;
                                
                size_t scindex_tns_0 = i * ctb_stride + j * 9 + 0;
                size_t scindex_tns_1 = i * ctb_stride + j * 9 + 1;
                size_t scindex_tns_2 = i * ctb_stride + j * 9 + 2;
                size_t scindex_tns_3 = i * ctb_stride + j * 9 + 3;
                size_t scindex_tns_4 = i * ctb_stride + j * 9 + 4;
                size_t scindex_tns_5 = i * ctb_stride + j * 9 + 5;
                size_t scindex_tns_6 = i * ctb_stride + j * 9 + 6;
                size_t scindex_tns_7 = i * ctb_stride + j * 9 + 7;
                size_t scindex_tns_8 = i * ctb_stride + j * 9 + 8;
                
                if(ctb_ptr[index] > 0.0)
                  {
                    if(*enable_density)
                      {
                        density_ptr[index] = density_cst * mass_ptr[index] / subcell_volume;
                      }
                    if(*enable_velocity_x)
                      {
                        velocity_x_ptr[index] = velocity_cst * momentum_x_ptr[index] / mass_ptr[index];
                      }
                    if(*enable_velocity_y)
                      {
                        velocity_y_ptr[index] = velocity_cst * momentum_y_ptr[index] / mass_ptr[index];
                      }
                    if(*enable_velocity_z)
                      {
                        velocity_z_ptr[index] = velocity_cst * momentum_z_ptr[index] / mass_ptr[index];
                      }
                    if(*enable_velocity_vector)
                      {
                        velocity_vector_ptr[scindex_vec_0] = velocity_cst * momentum_vector_ptr[scindex_vec_0] / mass_ptr[index];
                        velocity_vector_ptr[scindex_vec_1] = velocity_cst * momentum_vector_ptr[scindex_vec_1] / mass_ptr[index];
                        velocity_vector_ptr[scindex_vec_2] = velocity_cst * momentum_vector_ptr[scindex_vec_2] / mass_ptr[index];
                      }
                    if(*enable_virial_tensor)
                      {
                        virial_tensor_ptr[scindex_tns_0] /= ctb_ptr[index];
                        virial_tensor_ptr[scindex_tns_1] /= ctb_ptr[index];
                        virial_tensor_ptr[scindex_tns_2] /= ctb_ptr[index];
                        virial_tensor_ptr[scindex_tns_3] /= ctb_ptr[index];
                        virial_tensor_ptr[scindex_tns_4] /= ctb_ptr[index];
                        virial_tensor_ptr[scindex_tns_5] /= ctb_ptr[index];
                        virial_tensor_ptr[scindex_tns_6] /= ctb_ptr[index];
                        virial_tensor_ptr[scindex_tns_7] /= ctb_ptr[index];
                        virial_tensor_ptr[scindex_tns_8] /= ctb_ptr[index];
                      }
                    if(*enable_deformation_gradient_tensor)
                      {
                        deformation_gradient_tensor_ptr[scindex_tns_0] /= ctb_ptr[index];
                        deformation_gradient_tensor_ptr[scindex_tns_1] /= ctb_ptr[index];
                        deformation_gradient_tensor_ptr[scindex_tns_2] /= ctb_ptr[index];
                        deformation_gradient_tensor_ptr[scindex_tns_3] /= ctb_ptr[index];
                        deformation_gradient_tensor_ptr[scindex_tns_4] /= ctb_ptr[index];
                        deformation_gradient_tensor_ptr[scindex_tns_5] /= ctb_ptr[index];
                        deformation_gradient_tensor_ptr[scindex_tns_6] /= ctb_ptr[index];
                        deformation_gradient_tensor_ptr[scindex_tns_7] /= ctb_ptr[index];
                        deformation_gradient_tensor_ptr[scindex_tns_8] /= ctb_ptr[index];
                      }
                  }
              }
      }
      
    }
    
    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(Computes particle density)EOF";
    }


  private:
    static inline double get_mass( unsigned int j, const uint8_t* type_ptr, const double* masses, std::true_type )
    {
      return masses[type_ptr[j]];
    }

    static inline double get_mass( unsigned int, const uint8_t*, const double* masses, std::false_type )
    {
      return masses[0];
    }

    static inline double get_square_velocity( unsigned int j, const double * __restrict__ vx, const double * __restrict__ vy, const double * __restrict__ vz, std::true_type )
    {
      return vx[j]*vx[j] + vy[j]*vy[j] + vz[j]*vz[j];
    }
    
    static inline constexpr double get_square_velocity( unsigned int, const double * __restrict__, const double * __restrict__, const double * __restrict__, std::false_type )
    {
      return 1.0;
    }

    static inline double get_velocity_x( unsigned int j, const double * __restrict__ vx, std::true_type )
    {
      return vx[j];
    }
    
    static inline constexpr double get_velocity_x( unsigned int, const double * __restrict__, std::false_type )
    {
      return 1.0;
    }

    static inline double get_velocity_y( unsigned int j, const double * __restrict__ vy, std::true_type )
    {
      return vy[j];
    }
    
    static inline constexpr double get_velocity_y( unsigned int, const double * __restrict__, std::false_type )
    {
      return 1.0;
    }

    static inline double get_velocity_z( unsigned int j, const double * __restrict__ vz, std::true_type )
    {
      return vz[j];
    }
    
    static inline constexpr double get_velocity_z( unsigned int, const double * __restrict__, std::false_type )
    {
      return 1.0;
    }    
    
    static inline void localize_subcell( const Vec3d& r, double cell_size, double sub_cellsize, ssize_t subdiv, IJK& cell_loc, IJK& subcell_loc )
    {
      cell_loc = make_ijk( r / cell_size );
      Vec3d ro = r - (cell_loc*cell_size);
      subcell_loc = vclamp( make_ijk(ro / sub_cellsize) , 0 , subdiv-1 );
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

    // @return how much of this particle contributes to region cell_box.
    // sum of contributions for all disjoint cell_box paving the domain is guaranteed to be 1.0
    static inline double particle_smoothing(const Vec3d& r, double sp_size, const AABB& cell_box)
    {
      AABB contrib_box = { r - sp_size*0.5 , r + sp_size*0.5 };
      AABB sub_contrib_box = intersection( contrib_box , cell_box );
      double w = 0.0;
      if( ! is_nil(sub_contrib_box) ) { w = bounds_volume(sub_contrib_box) / (sp_size*sp_size*sp_size); }
      assert( w>=0. && w<=(1.0+1.e-11) );
      return w;
    }

  };

  template<class GridT> using GridReductionTmpl = GridReduction<GridT>;
  
  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory("grid_reduction", make_grid_variant_operator< GridReductionTmpl > );
  }

}
