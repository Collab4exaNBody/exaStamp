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
  class GridCellParticleSplatting : public OperatorNode
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
    ADD_SLOT( ParticleSpecies, species                    , INPUT , REQUIRED );
    ADD_SLOT( double         , splat_size                 , INPUT , REQUIRED );
    ADD_SLOT( long           , grid_subdiv                , INPUT_OUTPUT , 1 );
    ADD_SLOT( GridCellValues , grid_cell_values           , INPUT_OUTPUT );

    ADD_SLOT( bool           , enable_mass                , INPUT , true );
    ADD_SLOT( std::string    , mass_field_name            , INPUT_OUTPUT , "mass (Da)" );
    
    ADD_SLOT( bool           , enable_momentum_x          , INPUT , false );
    ADD_SLOT( std::string    , momentum_x_field_name      , INPUT_OUTPUT , "momentum_x (A/ps)" );    

    ADD_SLOT( bool           , enable_momentum_y          , INPUT , false );
    ADD_SLOT( std::string    , momentum_y_field_name      , INPUT_OUTPUT , "momentum_y (A/ps)" );    

    ADD_SLOT( bool           , enable_momentum_z          , INPUT , false );
    ADD_SLOT( std::string    , momentum_z_field_name      , INPUT_OUTPUT , "momentum_z (A/ps)" );    

    ADD_SLOT( bool           , enable_momentum_vector     , INPUT , true );
    ADD_SLOT( std::string    , momentum_vector_field_name , INPUT_OUTPUT , "momentum_vector (A/ps)" );

    ADD_SLOT( bool           , enable_m_v2                , INPUT , true );
    ADD_SLOT( std::string    , m_v2_field_name            , INPUT_OUTPUT , "m * v^2 (Da * (A/ps)^2)" );

    ADD_SLOT( bool           , enable_m_v2_tensor         , INPUT , false );
    ADD_SLOT( std::string    , m_v2_tensor_field_name     , INPUT_OUTPUT , "m * v^2 tensor (Da * (A/ps)^2)" );

    ADD_SLOT( bool           , enable_tr_virial           , INPUT , true );
    ADD_SLOT( std::string    , tr_virial_field_name       , INPUT_OUTPUT , "Tr(virial) (Da * A^2 / ps^2)?" );

    ADD_SLOT( bool           , enable_virial_tensor       , INPUT , false );
    ADD_SLOT( std::string    , virial_tensor_field_name   , INPUT_OUTPUT , "virial_tensor (Da * A^2 / ps^2)?" );

    ADD_SLOT( std::string    , contribution_field_name    , INPUT_OUTPUT , "contribution" );

    // -----------------------------------------------
    ADD_SLOT( GridParticleLocalMechanicalMetrics, local_mechanical_data , INPUT, OPTIONAL );

    ADD_SLOT( bool           , enable_deformation_gradient_tensor     , INPUT , false);
    ADD_SLOT( std::string    , deformation_gradient_tensor_field_name , INPUT_OUTPUT , "deformation_gradient_tensor (?)" );
    
  public:

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

      const bool output_contrib = *enable_mass || *enable_m_v2 || *enable_momentum_x || *enable_momentum_y || *enable_momentum_z || *enable_momentum_vector || *enable_m_v2_tensor || *enable_tr_virial || *enable_virial_tensor || *enable_deformation_gradient_tensor;
      const bool output_mass = *enable_mass;
      const bool output_momentum_x = has_velocity && *enable_momentum_x;
      const bool output_momentum_y = has_velocity && *enable_momentum_y;
      const bool output_momentum_z = has_velocity && *enable_momentum_z;
      const bool output_momentum_vector = has_velocity && *enable_momentum_vector;
      const bool output_virial_tensor = /*has_virial &&*/ *enable_virial_tensor;
      const bool output_m_v2 = has_velocity && *enable_m_v2;
      const bool output_m_v2_tensor = has_velocity && *enable_m_v2_tensor;
      const bool output_tr_virial = /*has_virial &&*/ *enable_tr_virial;
      const bool output_deformation_gradient_tensor = *enable_deformation_gradient_tensor;
      
      // create cell data fields if needed
      if(!output_contrib)
	{
	  return;
	}
      if(output_contrib && ! grid_cell_values->has_field(*contribution_field_name))
	{
	  grid_cell_values->add_field(*contribution_field_name, subdiv,1);
	}
      if( output_mass && ! grid_cell_values->has_field(*mass_field_name))
	{
	  grid_cell_values->add_field(*mass_field_name,subdiv,1);
	}
      if( output_momentum_x && ! grid_cell_values->has_field(*momentum_x_field_name) )
	{
	  grid_cell_values->add_field(*momentum_x_field_name,subdiv,1);
	}
      if( output_momentum_y && ! grid_cell_values->has_field(*momentum_y_field_name) )
	{
	  grid_cell_values->add_field(*momentum_y_field_name,subdiv,1);
	}
      if( output_momentum_z && ! grid_cell_values->has_field(*momentum_z_field_name) )
	{
	  grid_cell_values->add_field(*momentum_z_field_name,subdiv,1);
	}
      if( output_momentum_vector && ! grid_cell_values->has_field(*momentum_vector_field_name) )
	{
	  grid_cell_values->add_field(*momentum_vector_field_name,subdiv,3);
	}
      if( output_virial_tensor && ! grid_cell_values->has_field(*virial_tensor_field_name) )
	{
	  grid_cell_values->add_field(*virial_tensor_field_name,subdiv,9);
	} 
      if( output_m_v2 && ! grid_cell_values->has_field(*m_v2_field_name) )
	{
	  grid_cell_values->add_field(*m_v2_field_name,subdiv,1);
	}
      if( output_m_v2_tensor && ! grid_cell_values->has_field(*m_v2_tensor_field_name) )
	{
	  grid_cell_values->add_field(*m_v2_tensor_field_name,subdiv,9);
	}
      if( output_tr_virial && ! grid_cell_values->has_field(*tr_virial_field_name) )
	{
	  grid_cell_values->add_field(*tr_virial_field_name,subdiv,1);
	}
      if( output_deformation_gradient_tensor && ! grid_cell_values->has_field(*deformation_gradient_tensor_field_name) )
	{
	  grid_cell_values->add_field(*deformation_gradient_tensor_field_name,subdiv,9);
	}

      // -------------------------------------------------------------------------------------------------------------
      
      // retreive contribution field data accessor. create data field if needed
      double * __restrict__ ctb_ptr = nullptr;
      size_t ctb_stride = 0;
      if( output_contrib )
        {
          assert( size_t(subdiv) == grid_cell_values->field(*contribution_field_name).m_subdiv );
          assert( size_t(subdiv * subdiv * subdiv) == grid_cell_values->field(*contribution_field_name).m_components );
          auto cell_contribution_accessor = grid_cell_values->field_data(*contribution_field_name);
          ctb_ptr = cell_contribution_accessor.m_data_ptr;
          ctb_stride = cell_contribution_accessor.m_stride;
          // field is set to 0 in the parallel region
        }

      // retreive mass field data accessor. create data field if needed
      double * __restrict__ mass_ptr = nullptr;
      size_t mass_stride = ctb_stride;
      if( output_mass )
        {
          assert( size_t(subdiv) == grid_cell_values->field(*mass_field_name).m_subdiv );
          assert( size_t(subdiv * subdiv * subdiv) == grid_cell_values->field(*mass_field_name).m_components );
          auto cell_mass_accessor = grid_cell_values->field_data(*mass_field_name);
          mass_ptr = cell_mass_accessor.m_data_ptr;
          mass_stride = cell_mass_accessor.m_stride;
          assert( ctb_stride == mass_stride );
          // field is set to 0 in the parallel region
        }

      // retreive momentum_x field data accessor. create data field if needed
      double * __restrict__ momentum_x_ptr = nullptr;
      size_t momentum_x_stride = ctb_stride;
      if( output_momentum_x )
        {
          assert( subdiv == grid_cell_values->field(*momentum_x_field_name).m_subdiv );
          assert( (subdiv * subdiv * subdiv) == grid_cell_values->field(*momentum_x_field_name).m_components );
          auto cell_momentum_x_accessor = grid_cell_values->field_data(*momentum_x_field_name);
          momentum_x_ptr = cell_momentum_x_accessor.m_data_ptr;
          momentum_x_stride = cell_momentum_x_accessor.m_stride;
          assert( momentum_x_stride == ctb_stride );
          // field is set to 0 in the parallel region
        }

      // retreive momentum_y field data accessor. create data field if needed
      double * __restrict__ momentum_y_ptr = nullptr;
      size_t momentum_y_stride = ctb_stride;
      if( output_momentum_y )
        {
          assert( subdiv == grid_cell_values->field(*momentum_y_field_name).m_subdiv );
          assert( (subdiv * subdiv * subdiv) == grid_cell_values->field(*momentum_y_field_name).m_components );
          auto cell_momentum_y_accessor = grid_cell_values->field_data(*momentum_y_field_name);
          momentum_y_ptr = cell_momentum_y_accessor.m_data_ptr;
          momentum_y_stride = cell_momentum_y_accessor.m_stride;
          assert( momentum_y_stride == ctb_stride );
          // field is set to 0 in the parallel region
        }

      // retreive momentum_z field data accessor. create data field if needed
      double * __restrict__ momentum_z_ptr = nullptr;
      size_t momentum_z_stride = ctb_stride;
      if( output_momentum_z )
        {
          assert( subdiv == grid_cell_values->field(*momentum_z_field_name).m_subdiv );
          assert( (subdiv * subdiv * subdiv) == grid_cell_values->field(*momentum_z_field_name).m_components );
          auto cell_momentum_z_accessor = grid_cell_values->field_data(*momentum_z_field_name);
          momentum_z_ptr = cell_momentum_z_accessor.m_data_ptr;
          momentum_z_stride = cell_momentum_z_accessor.m_stride;
          assert( momentum_z_stride == ctb_stride );
          // field is set to 0 in the parallel region
        }
      
      // retreive momentum vector field data accessor. create data field if needed
      double * __restrict__ momentum_vector_ptr = nullptr;
      size_t momentum_vector_stride = 0;
      if( output_momentum_vector )
        {
          assert( subdiv == grid_cell_values->field(*momentum_vector_field_name).m_subdiv );
          assert( (subdiv * subdiv * subdiv * 3) == grid_cell_values->field(*momentum_vector_field_name).m_components );
          assert( (subdiv * subdiv * subdiv) == n_subcells );
          auto cell_momentum_vector_accessor = grid_cell_values->field_data(*momentum_vector_field_name);
          momentum_vector_ptr = cell_momentum_vector_accessor.m_data_ptr;
          momentum_vector_stride = cell_momentum_vector_accessor.m_stride;
          assert( momentum_vector_stride == ctb_stride );
          // field is set to 0 in the parallel region
        }

      // retreive virial tensor field data accessor. create data field if needed
      double * __restrict__ virial_tensor_ptr = nullptr;
      size_t virial_tensor_stride = ctb_stride;
      if( output_virial_tensor )
        {
          assert( subdiv == grid_cell_values->field(*virial_tensor_field_name).m_subdiv );
          assert( (subdiv * subdiv * subdiv * 9) == grid_cell_values->field(*virial_tensor_field_name).m_components );
          assert( (subdiv * subdiv * subdiv) == n_subcells );
          auto cell_virial_tensor_accessor = grid_cell_values->field_data(*virial_tensor_field_name);
          virial_tensor_ptr = cell_virial_tensor_accessor.m_data_ptr;
          virial_tensor_stride = cell_virial_tensor_accessor.m_stride;
          assert( virial_tensor_stride == ctb_stride );
          // field is set to 0 in the parallel region
          int n_comps_per_subcell = grid_cell_values->field(*virial_tensor_field_name).m_components / (subdiv*subdiv*subdiv);
          assert( n_comps_per_subcell == 9 );
          // lout << "n comps per subcell : " << n_comps_per_subcell << std::endl;
          for(ssize_t i=0;i<n_cells;i++)
            {
              for(int j=0;j<n_subcells;j++)
                {
                  for(int k=0;k<n_comps_per_subcell;k++)
                    {
                      virial_tensor_ptr[ i * virial_tensor_stride + j * n_comps_per_subcell + k] = 0.0;
                    }
                }
            }
        }

      // retreive mass * v^2 field data accessor. create data field if needed
      double * __restrict__ m_v2_ptr = nullptr;
      size_t m_v2_stride = ctb_stride;
      if( output_m_v2 )
        {
          assert( subdiv == grid_cell_values->field(*m_v2_field_name).m_subdiv );
          assert( (subdiv * subdiv * subdiv) == grid_cell_values->field(*m_v2_field_name).m_components );
          auto cell_m_v2_accessor = grid_cell_values->field_data(*m_v2_field_name);
          m_v2_ptr = cell_m_v2_accessor.m_data_ptr;
          m_v2_stride = cell_m_v2_accessor.m_stride;
          assert( m_v2_stride == ctb_stride );
          // field is set to 0 in the parallel region
        }

      // retreive m_v2 tensor field data accessor. create data field if needed
      double * __restrict__ m_v2_tensor_ptr = nullptr;
      size_t m_v2_tensor_stride = ctb_stride;
      if( output_m_v2_tensor )
        {
          assert( subdiv == grid_cell_values->field(*m_v2_tensor_field_name).m_subdiv );
          assert( (subdiv * subdiv * subdiv * 9) == grid_cell_values->field(*m_v2_tensor_field_name).m_components );
          assert( (subdiv * subdiv * subdiv) == n_subcells );
          auto cell_m_v2_tensor_accessor = grid_cell_values->field_data(*m_v2_tensor_field_name);
          m_v2_tensor_ptr = cell_m_v2_tensor_accessor.m_data_ptr;
          m_v2_tensor_stride = cell_m_v2_tensor_accessor.m_stride;
          assert( m_v2_tensor_stride == ctb_stride );
          // field is set to 0 in the parallel region
          int n_comps_per_subcell = grid_cell_values->field(*m_v2_tensor_field_name).m_components / (subdiv*subdiv*subdiv);
          assert( n_comps_per_subcell == 9 );
          // lout << "n comps per subcell : " << n_comps_per_subcell << std::endl;
          for(ssize_t i=0;i<n_cells;i++)
            {
              for(int j=0;j<n_subcells;j++)
                {
                  for(int k=0;k<n_comps_per_subcell;k++)
                    {
                      m_v2_tensor_ptr[ i * m_v2_tensor_stride + j * n_comps_per_subcell + k] = 0.0;
                    }
                }
            }
        }
      
      // retreive tr_virial field data accessor. create data field if needed
      double * __restrict__ tr_virial_ptr = nullptr;
      size_t tr_virial_stride = ctb_stride;
      if( output_tr_virial )
        {
          assert( subdiv == grid_cell_values->field(*tr_virial_field_name).m_subdiv );
          assert( (subdiv * subdiv * subdiv) == grid_cell_values->field(*tr_virial_field_name).m_components );
          auto cell_tr_virial_accessor = grid_cell_values->field_data(*tr_virial_field_name);
          tr_virial_ptr = cell_tr_virial_accessor.m_data_ptr;
          tr_virial_stride = cell_tr_virial_accessor.m_stride;
          assert( tr_virial_stride == ctb_stride );
          // field is set to 0 in the parallel region
        }

      // retreive deformation gradient tensor field data accessor. create data field if needed
      double * __restrict__ deformation_gradient_tensor_ptr = nullptr;
      size_t deformation_gradient_tensor_stride = ctb_stride;
      if( output_deformation_gradient_tensor )
        {
          assert( subdiv == grid_cell_values->field(*deformation_gradient_tensor_field_name).m_subdiv );
          assert( (subdiv * subdiv * subdiv * 9) == grid_cell_values->field(*deformation_gradient_tensor_field_name).m_components );
          assert( (subdiv * subdiv * subdiv) == n_subcells );
          auto cell_deformation_gradient_tensor_accessor = grid_cell_values->field_data(*deformation_gradient_tensor_field_name);
          deformation_gradient_tensor_ptr = cell_deformation_gradient_tensor_accessor.m_data_ptr;
          deformation_gradient_tensor_stride = cell_deformation_gradient_tensor_accessor.m_stride;
          assert( deformation_gradient_tensor_stride == ctb_stride );
          // field is set to 0 in the parallel region
        }
      
      // particle splatting size
      const double sp_size = *splat_size;
      if( sp_size/2 > subcell_size )
        {
          lerr << "in " << pathname() << std::endl
               << "splat_size = "<<sp_size << " is 2 times larger than subcell size = "<<subcell_size<< std::endl
               << "Choose smaller splat_size or coarser grid_subdiv value" << std::endl;
          std::abort();
        }
  
      // store specie masses in temporary array
      size_t nSpecies = species->size();
      if( nSpecies != 1 && !has_type_field )
        {
          lerr << pathname() << std::endl;
          lerr << "no type information, can't retreive masses" << std::endl;
          std::abort();
        }
      double masses[nSpecies];
      for(size_t i=0;i<nSpecies;i++) { masses[i] = species->at(i).m_mass; }

      auto cells = grid->cells();
      const IJK dims = grid->dimension();
      const ssize_t gl = grid->ghost_layers();      
      
      ldbg << "ghost_layers="<<gl<<", cell_size="<<cell_size<<", splat_size="<<sp_size<<", subdiv="<<subdiv<<", subcell_size="<<subcell_size<<std::endl;
      
      // computes per cell Mass, per cell Mass*Velocity^2, per cell x-component of Momentum and per cell x-component or Virial
#     pragma omp parallel
      {
	// we need to set our fields to 0 
	_Pragma(USTAMP_STR(omp for collapse(2) schedule(dynamic)))
	  for(ssize_t i = 0; i < dims.i*dims.j*dims.k; ++i)
	    for(ssize_t j = 0; j < subdiv*subdiv*subdiv; ++j)
	      {
		size_t index = i * ctb_stride + j;
		if( output_contrib )
		  {
		    ctb_ptr[ index ] = 0.0;
		  }
		if( output_mass )
		  {
		    mass_ptr[ index ] = 0.0;
		  }
		if( output_momentum_x )
		  {
		    momentum_x_ptr[ index ] = 0.0;
		  }
		if( output_momentum_y )
		  {
		    momentum_y_ptr[ index ] = 0.0;
		  }
		if( output_momentum_z )
		  {
		    momentum_z_ptr[ index ] = 0.0;
		  }
		if( output_momentum_vector )
		  {
		    size_t index_vec_0 = i * ctb_stride + j * 3 + 0;
		    size_t index_vec_1 = i * ctb_stride + j * 3 + 1;
		    size_t index_vec_2 = i * ctb_stride + j * 3 + 2;
		    momentum_vector_ptr[ index_vec_0 ] = 0.0;
		    momentum_vector_ptr[ index_vec_1 ] = 0.0;
		    momentum_vector_ptr[ index_vec_2 ] = 0.0;
		  }
		if( output_virial_tensor )
		  {
		    size_t index_tens_0 = i * ctb_stride + j * 9 + 0;
		    size_t index_tens_1 = i * ctb_stride + j * 9 + 1;
		    size_t index_tens_2 = i * ctb_stride + j * 9 + 2;
		    size_t index_tens_3 = i * ctb_stride + j * 9 + 3;
		    size_t index_tens_4 = i * ctb_stride + j * 9 + 4;
		    size_t index_tens_5 = i * ctb_stride + j * 9 + 5;
		    size_t index_tens_6 = i * ctb_stride + j * 9 + 6;
		    size_t index_tens_7 = i * ctb_stride + j * 9 + 7;
		    size_t index_tens_8 = i * ctb_stride + j * 9 + 8;
		    virial_tensor_ptr[ index_tens_0 ] = 0.0;
		    virial_tensor_ptr[ index_tens_1 ] = 0.0;
		    virial_tensor_ptr[ index_tens_2 ] = 0.0;
		    virial_tensor_ptr[ index_tens_3 ] = 0.0;
		    virial_tensor_ptr[ index_tens_4 ] = 0.0;
		    virial_tensor_ptr[ index_tens_5 ] = 0.0;
		    virial_tensor_ptr[ index_tens_6 ] = 0.0;
		    virial_tensor_ptr[ index_tens_7 ] = 0.0;
		    virial_tensor_ptr[ index_tens_8 ] = 0.0;
		  }
		if( output_m_v2 )
		  {
		    m_v2_ptr[ index ] = 0.0;
		  }
		if( output_m_v2_tensor )
		  {
		    size_t index_tens_0 = i * ctb_stride + j * 9 + 0;
		    size_t index_tens_1 = i * ctb_stride + j * 9 + 1;
		    size_t index_tens_2 = i * ctb_stride + j * 9 + 2;
		    size_t index_tens_3 = i * ctb_stride + j * 9 + 3;
		    size_t index_tens_4 = i * ctb_stride + j * 9 + 4;
		    size_t index_tens_5 = i * ctb_stride + j * 9 + 5;
		    size_t index_tens_6 = i * ctb_stride + j * 9 + 6;
		    size_t index_tens_7 = i * ctb_stride + j * 9 + 7;
		    size_t index_tens_8 = i * ctb_stride + j * 9 + 8;
		    m_v2_tensor_ptr[ index_tens_0 ] = 0.0;
		    m_v2_tensor_ptr[ index_tens_1 ] = 0.0;
		    m_v2_tensor_ptr[ index_tens_2 ] = 0.0;
		    m_v2_tensor_ptr[ index_tens_3 ] = 0.0;
		    m_v2_tensor_ptr[ index_tens_4 ] = 0.0;
		    m_v2_tensor_ptr[ index_tens_5 ] = 0.0;
		    m_v2_tensor_ptr[ index_tens_6 ] = 0.0;
		    m_v2_tensor_ptr[ index_tens_7 ] = 0.0;
		    m_v2_tensor_ptr[ index_tens_8 ] = 0.0;
		  }
		if( output_tr_virial )
		  {
		    tr_virial_ptr[ index ] = 0.0;
		  }
		if( output_deformation_gradient_tensor )
		  {
		    size_t index_tens_0 = i * ctb_stride + j * 9 + 0;
		    size_t index_tens_1 = i * ctb_stride + j * 9 + 1;
		    size_t index_tens_2 = i * ctb_stride + j * 9 + 2;
		    size_t index_tens_3 = i * ctb_stride + j * 9 + 3;
		    size_t index_tens_4 = i * ctb_stride + j * 9 + 4;
		    size_t index_tens_5 = i * ctb_stride + j * 9 + 5;
		    size_t index_tens_6 = i * ctb_stride + j * 9 + 6;
		    size_t index_tens_7 = i * ctb_stride + j * 9 + 7;
		    size_t index_tens_8 = i * ctb_stride + j * 9 + 8;
		    deformation_gradient_tensor_ptr[ index_tens_0 ] = 0.0;
		    deformation_gradient_tensor_ptr[ index_tens_1 ] = 0.0;
		    deformation_gradient_tensor_ptr[ index_tens_2 ] = 0.0;
		    deformation_gradient_tensor_ptr[ index_tens_3 ] = 0.0;
		    deformation_gradient_tensor_ptr[ index_tens_4 ] = 0.0;
		    deformation_gradient_tensor_ptr[ index_tens_5 ] = 0.0;
		    deformation_gradient_tensor_ptr[ index_tens_6 ] = 0.0;
		    deformation_gradient_tensor_ptr[ index_tens_7 ] = 0.0;
		    deformation_gradient_tensor_ptr[ index_tens_8 ] = 0.0;
		  }
	      }
	
        GRID_OMP_FOR_BEGIN(dims,i,cell_loc, schedule(dynamic) )
          {
            const Vec3d cell_origin = grid->cell_position( cell_loc );
            const auto* __restrict__ rx = cells[i][field::rx];
            const auto* __restrict__ ry = cells[i][field::ry];
            const auto* __restrict__ rz = cells[i][field::rz];
            const auto* __restrict__ atom_type = cells[i].field_pointer_or_null(field::type);
            const auto* __restrict__ vx = cells[i].field_pointer_or_null(field::vx);
            const auto* __restrict__ vy = cells[i].field_pointer_or_null(field::vy);
            const auto* __restrict__ vz = cells[i].field_pointer_or_null(field::vz);
            const auto* __restrict__ vir = cells[i].field_pointer_or_null(field::virial);
            const unsigned int n = cells[i].size();
            for(unsigned int j=0;j<n;j++)
              {
                Vec3d r { rx[j] , ry[j] , rz[j] };
                const double mass = get_mass( j, atom_type, masses, has_type_field );
                const double v2   = get_square_velocity( j, vx,vy,vz, has_velocity );
                const double velx = get_velocity_x( j, vx, has_velocity );
                const double vely = get_velocity_y( j, vy, has_velocity );
                const double velz = get_velocity_z( j, vz, has_velocity );
                
                const Mat3d  pvir = get_particle_virial<has_virial_field>( vir, j );
                const double tr_vir = pvir.m11 + pvir.m22 + pvir.m33;

		Mat3d deformation_gradient_tensor_tmp;
		if (*enable_deformation_gradient_tensor)
		  {
		    deformation_gradient_tensor_tmp = (*local_mechanical_data)[i].F[j];
		  }
		const Mat3d deformation_gradient_tensor = deformation_gradient_tensor_tmp;

                IJK center_cell_loc;
                IJK center_subcell_loc;
                Vec3d rco = r - cell_origin;
                localize_subcell( rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc );
                center_cell_loc += cell_loc;

                for(int ck=-1;ck<=1;ck++)
                  for(int cj=-1;cj<=1;cj++)
                    for(int ci=-1;ci<=1;ci++)
                      {
                        IJK nbh_cell_loc;
                        IJK nbh_subcell_loc;
                        subcell_neighbor( center_cell_loc, center_subcell_loc, subdiv, IJK{ci,cj,ck}, nbh_cell_loc, nbh_subcell_loc );
                        if( grid->contains(nbh_cell_loc) )
                          {
                            ssize_t nbh_cell_i = grid_ijk_to_index( dims , nbh_cell_loc );
                            ssize_t nbh_subcell_i = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , nbh_subcell_loc );
                            assert( nbh_cell_i>=0 && nbh_cell_i<n_cells );
                            assert( nbh_subcell_i>=0 && nbh_subcell_i<n_subcells );

                            // compute weighted contribution of particle to sub cell
                            Vec3d nbh_cell_origin = grid->cell_position(nbh_cell_loc);
                            AABB subcell_box = { nbh_cell_origin + nbh_subcell_loc*subcell_size , nbh_cell_origin + (nbh_subcell_loc+1)*subcell_size };
                            const double w = particle_smoothing(r, sp_size, subcell_box);
                            double mass_contrib = mass * w;
                            double momentum_x_contrib = velx * mass_contrib;
                            double momentum_y_contrib = vely * mass_contrib;
                            double momentum_z_contrib = velz * mass_contrib;
                            double virial_11_contrib = pvir.m11 * w;
                            double virial_12_contrib = pvir.m12 * w;
                            double virial_13_contrib = pvir.m13 * w;
                            double virial_21_contrib = pvir.m21 * w;
                            double virial_22_contrib = pvir.m22 * w;
                            double virial_23_contrib = pvir.m23 * w;
                            double virial_31_contrib = pvir.m31 * w;
                            double virial_32_contrib = pvir.m32 * w;
                            double virial_33_contrib = pvir.m33 * w;
                            double m_v2_contrib = v2 * mass_contrib;
                            double m_v1_v1_contrib = velx * velx * mass_contrib;
                            double m_v1_v2_contrib = velx * vely * mass_contrib;
                            double m_v1_v3_contrib = velx * velz * mass_contrib;
                            double m_v2_v1_contrib = vely * velx * mass_contrib;
                            double m_v2_v2_contrib = vely * vely * mass_contrib;
                            double m_v2_v3_contrib = vely * velz * mass_contrib;
                            double m_v3_v1_contrib = velz * velx * mass_contrib;
                            double m_v3_v2_contrib = velz * vely * mass_contrib;
                            double m_v3_v3_contrib = velz * velz * mass_contrib;
                            double tr_virial_contrib = tr_vir * w;

			    double F_11_contrib = deformation_gradient_tensor.m11 * w; // w is choosen for now but may be the mass contrib is more judicious
			    double F_12_contrib = deformation_gradient_tensor.m12 * w;
			    double F_13_contrib = deformation_gradient_tensor.m13 * w;
			    double F_21_contrib = deformation_gradient_tensor.m21 * w;
			    double F_22_contrib = deformation_gradient_tensor.m22 * w;
			    double F_23_contrib = deformation_gradient_tensor.m23 * w;
			    double F_31_contrib = deformation_gradient_tensor.m31 * w;
			    double F_32_contrib = deformation_gradient_tensor.m32 * w;
			    double F_33_contrib = deformation_gradient_tensor.m33 * w;
			    
                            size_t scindex = nbh_cell_i * ctb_stride + nbh_subcell_i;
                            
                            if( output_contrib )
                              {
#                               pragma omp atomic update
                                ctb_ptr[ scindex ] += w;
                              }

                            if( output_mass )
                              {
#                               pragma omp atomic update
                                mass_ptr[ scindex ] += mass_contrib;
                              }
                            
                            if( output_momentum_x )
                              {
#                               pragma omp atomic update
                                momentum_x_ptr[ scindex ] += momentum_x_contrib;
                              }

                            if( output_momentum_y )
                              {
#                               pragma omp atomic update
                                momentum_y_ptr[ scindex ] += momentum_y_contrib;
                              }

                            if( output_momentum_z )
                              {
#                               pragma omp atomic update
                                momentum_z_ptr[ scindex ] += momentum_z_contrib;
                              }                         

                            if( output_momentum_vector )
                              {
                                size_t scindex_vec_0 = nbh_cell_i * ctb_stride + nbh_subcell_i * 3 + 0;
                                size_t scindex_vec_1 = nbh_cell_i * ctb_stride + nbh_subcell_i * 3 + 1;
                                size_t scindex_vec_2 = nbh_cell_i * ctb_stride + nbh_subcell_i * 3 + 2;
#                               pragma omp atomic update
                                momentum_vector_ptr[ scindex_vec_0 ] += momentum_x_contrib;
#                               pragma omp atomic update
                                momentum_vector_ptr[ scindex_vec_1 ] += momentum_y_contrib;
#                               pragma omp atomic update
                                momentum_vector_ptr[ scindex_vec_2 ] += momentum_z_contrib;
                              }

                            if( output_virial_tensor )
                              {
                                size_t scindex_vec_0 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 0;
                                size_t scindex_vec_1 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 1;
                                size_t scindex_vec_2 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 2;
                                size_t scindex_vec_3 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 3;
                                size_t scindex_vec_4 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 4;
                                size_t scindex_vec_5 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 5;
                                size_t scindex_vec_6 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 6;
                                size_t scindex_vec_7 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 7;
                                size_t scindex_vec_8 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 8;
#                               pragma omp atomic update
                                virial_tensor_ptr[ scindex_vec_0 ] += virial_11_contrib;
#                               pragma omp atomic update
                                virial_tensor_ptr[ scindex_vec_1 ] += virial_12_contrib;
#                               pragma omp atomic update
                                virial_tensor_ptr[ scindex_vec_2 ] += virial_13_contrib;
#                               pragma omp atomic update
                                virial_tensor_ptr[ scindex_vec_3 ] += virial_21_contrib;
#                               pragma omp atomic update
                                virial_tensor_ptr[ scindex_vec_4 ] += virial_22_contrib;
#                               pragma omp atomic update
                                virial_tensor_ptr[ scindex_vec_5 ] += virial_23_contrib;
#                               pragma omp atomic update
                                virial_tensor_ptr[ scindex_vec_6 ] += virial_31_contrib;
#                               pragma omp atomic update
                                virial_tensor_ptr[ scindex_vec_7 ] += virial_32_contrib;
#                               pragma omp atomic update
                                virial_tensor_ptr[ scindex_vec_8 ] += virial_33_contrib;
                              }

                            if( output_m_v2 && m_v2_contrib>0.0 )
                              {
#                               pragma omp atomic update
                                m_v2_ptr[ scindex ] += m_v2_contrib; 
                              }

                            if( output_m_v2_tensor )
                              {
                                size_t scindex_vec_0 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 0;
                                size_t scindex_vec_1 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 1;
                                size_t scindex_vec_2 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 2;
                                size_t scindex_vec_3 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 3;
                                size_t scindex_vec_4 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 4;
                                size_t scindex_vec_5 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 5;
                                size_t scindex_vec_6 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 6;
                                size_t scindex_vec_7 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 7;
                                size_t scindex_vec_8 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 8;
#                               pragma omp atomic update
                                m_v2_tensor_ptr[ scindex_vec_0 ] += m_v1_v1_contrib;
#                               pragma omp atomic update
                                m_v2_tensor_ptr[ scindex_vec_1 ] += m_v1_v2_contrib;
#                               pragma omp atomic update
                                m_v2_tensor_ptr[ scindex_vec_2 ] += m_v1_v3_contrib;
#                               pragma omp atomic update
                                m_v2_tensor_ptr[ scindex_vec_3 ] += m_v2_v1_contrib;
#                               pragma omp atomic update
                                m_v2_tensor_ptr[ scindex_vec_4 ] += m_v2_v2_contrib;
#                               pragma omp atomic update
                                m_v2_tensor_ptr[ scindex_vec_5 ] += m_v2_v3_contrib;
#                               pragma omp atomic update
                                m_v2_tensor_ptr[ scindex_vec_6 ] += m_v3_v1_contrib;
#                               pragma omp atomic update
                                m_v2_tensor_ptr[ scindex_vec_7 ] += m_v3_v2_contrib;
#                               pragma omp atomic update
                                m_v2_tensor_ptr[ scindex_vec_8 ] += m_v3_v3_contrib;
                              }
                            
                            if( output_tr_virial )
                              {
#                               pragma omp atomic update
                                tr_virial_ptr[ scindex ] += tr_virial_contrib; 
                              }
			    
			    if( output_deformation_gradient_tensor )
                              {
                                size_t scindex_vec_0 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 0;
                                size_t scindex_vec_1 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 1;
                                size_t scindex_vec_2 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 2;
                                size_t scindex_vec_3 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 3;
                                size_t scindex_vec_4 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 4;
                                size_t scindex_vec_5 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 5;
                                size_t scindex_vec_6 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 6;
                                size_t scindex_vec_7 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 7;
                                size_t scindex_vec_8 = nbh_cell_i * ctb_stride + nbh_subcell_i * 9 + 8;
#                               pragma omp atomic update
                                deformation_gradient_tensor_ptr[ scindex_vec_0 ] += F_11_contrib;
#                               pragma omp atomic update
                                deformation_gradient_tensor_ptr[ scindex_vec_1 ] += F_12_contrib;
#                               pragma omp atomic update
                                deformation_gradient_tensor_ptr[ scindex_vec_2 ] += F_13_contrib;
#                               pragma omp atomic update
                                deformation_gradient_tensor_ptr[ scindex_vec_3 ] += F_21_contrib;
#                               pragma omp atomic update
                                deformation_gradient_tensor_ptr[ scindex_vec_4 ] += F_22_contrib;
#                               pragma omp atomic update
                                deformation_gradient_tensor_ptr[ scindex_vec_5 ] += F_23_contrib;
#                               pragma omp atomic update
                                deformation_gradient_tensor_ptr[ scindex_vec_6 ] += F_31_contrib;
#                               pragma omp atomic update
                                deformation_gradient_tensor_ptr[ scindex_vec_7 ] += F_32_contrib;
#                               pragma omp atomic update
                                deformation_gradient_tensor_ptr[ scindex_vec_8 ] += F_33_contrib;
                              }
                          }
                      }
              }
          }
        GRID_OMP_FOR_END;
        
        // _Pragma(USTAMP_STR(omp for collapse(2) schedule(dynamic)))
        //   for(ssize_t i = 0; i < dims.i*dims.j*dims.k; ++i)
        //     for(ssize_t j = 0; j < subdiv*subdiv*subdiv; ++j)
        //       {
        //         size_t index = i * ctb_stride + j;
        //         if(ctb_ptr[index] > 0.0);
        //       }
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

  template<class GridT> using GridCellParticleSplattingTmpl = GridCellParticleSplatting<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory("grid_cell_particle_splatting", make_grid_variant_operator< GridCellParticleSplattingTmpl > );
  }

}
