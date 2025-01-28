#include <chrono>
#include <ctime>
#include <mpi.h>
#include <string>
#include <numeric>
#include <memory>
#include <random>
#include <cmath>

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/domain.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/simple_block_rcb.h>
//#include "exanb/container_utils.h"

#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types_stream.h>
#include <onika/log.h>
//#include "exanb/vector_utils.h"
#include <exanb/core/file_utils.h>
#include <exanb/core/check_particles_inside_cell.h>
#include <onika/physics/constants.h>
#include <onika/parallel/random.h>

#include <exanb/core/thread.h>

namespace exaStamp
{
  using namespace exanb;

  // Generate Orthorhombic lattice, i.e. a, b, c may have different lengths but alpha, beta and gamma equal 90 degrees.
  // Requirements :
  // structure : SC, BCC, FCC, HCP, DIAMOND, ROCKSALT, FLUORITE, C15, PEROVSKITE, ST, BCT, FCT, 2BCT
  // types : number of types depends on the structure and types must be consistent with the species defined beforehand
  // size : vector containing a, b, c lengths of the lattice
  // repeats : number of repetition along three laboratory directions

  template<typename GridT>
  class OrthorhombicLatticeOperator : public OperatorNode
  {
    using LatticeAtoms = std::vector<std::string>;
    //using LatticePositions = std::vector<Vec3d>;

    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------    
    ADD_SLOT( MPI_Comm        , mpi          , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( ReadBoundsSelectionMode, bounds_mode   , INPUT , ReadBoundsSelectionMode::FILE_BOUNDS );
    ADD_SLOT( Domain          , domain       , INPUT_OUTPUT );
    ADD_SLOT( GridT           , grid         , INPUT_OUTPUT );
    ADD_SLOT( double          , enlarge_bounds, INPUT , 0.0 );
    ADD_SLOT( bool            , pbc_adjust_xform , INPUT , false );
    ADD_SLOT( ParticleSpecies , species      , INPUT ); // optional. if no species given, type ids are allocated automatically

    // Variables related to the crystal structure
    ADD_SLOT( std::string      , structure , INPUT , REQUIRED );
    ADD_SLOT( LatticeAtoms     , types     , INPUT , REQUIRED );    
    ADD_SLOT( IJK              , repeats   , INPUT , REQUIRED );
    ADD_SLOT( Vec3d            , size      , INPUT , REQUIRED );    
    ADD_SLOT( bool             , gaussian_noise        , INPUT , false);    
    ADD_SLOT( double           , gaussian_noise_sigma  , INPUT , 0.0);
    ADD_SLOT( Vec3d            , position_shift  , INPUT , Vec3d{0., 0., 0.} );

    // Variables related to the special geometry, here a plane above/below which we keep/remove the particles. WARNING : be careful with the PBC
    ADD_SLOT( bool             , geo_plane    , INPUT , false);
    ADD_SLOT( Vec3d            , plane_normal , INPUT , Vec3d{0., 0., 1.});
    ADD_SLOT( double           , plane_offset , INPUT , 0.);

    // Variables related to the special geometry, here a cylinder inside/outside which we keep/remove the particles. WARNING : be careful with the PBC    
    ADD_SLOT( bool             , geo_cylinder    , INPUT , false);
    ADD_SLOT( bool             , cylinder_inside , INPUT , true);    
    ADD_SLOT( Vec3d            , cylinder_axis   , INPUT , Vec3d{0., 0., 1.});
    ADD_SLOT( double           , cylinder_radius , INPUT , 0.);

    // Variables related to the special geometry, here a cylinder inside/outside which we keep/remove the particles. WARNING : be careful with the PBC    
    ADD_SLOT( bool             , geo_void           , INPUT , false);
    ADD_SLOT( std::string      , void_mode          , INPUT , "simple"); //simple is the one void mode, possiblities: simple, porosity
    ADD_SLOT( Vec3d            , void_center        , INPUT , Vec3d{0., 0., 0.});
    ADD_SLOT( double           , void_radius        , INPUT , 0.);
    ADD_SLOT( double           , void_porosity      , INPUT , 0.);
    ADD_SLOT( double           , void_mean_diameter , INPUT , 0.);        
    
  public:
    inline void execute () override final
    {

      Domain& domain = *(this->domain);
      GridT& grid = *(this->grid);
      std::string structure = *(this->structure);

      bool is_noise = *(this->gaussian_noise);            
      double sigma_noise = *(this->gaussian_noise_sigma);

      if( *pbc_adjust_xform )
      {
        if( ! domain.xform_is_identity() )
        {
          lout << "pbc_adjust_xform needs initial XForm to be identity, resetting XForm"<<std::endl;
          domain.set_xform( make_identity_matrix() );
        }
      }
      
      using ParticleTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_id, field::_type>;
      using ParticleTuple = decltype( grid.cells()[0][0] );

      assert( species->size() > 0 );
      assert( grid.number_of_particles() == 0 );

      // MPI Initialization
      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      // uint64_t n_particles = 0;

      // Verification du nombre de types d'atomes en fonction du type de maille cubique voulu
      // ***choix**  **nb_types**
      //         SC   1  
      //        BCC   2
      //        FCC   4
      //        HCP   4
      //    DIAMOND   8
      //   ROCKSALT   8
      //   FLUORITE  12
      //        C15  24
      // PEROVSKITE   5
      //         ST   1
      //        BCT   2
      //        FCT   4
      //       2BCT   4
      Vec3d size = *(this->size);
      uint64_t n_particles_cell = 0.;
      std::vector<Vec3d> positions;
      if (structure == "SC" ) {
        n_particles_cell = 1;
        if( types->size() != n_particles_cell )
          {
            lerr << "Parameter types for simple cubic lattice must contain exactly 1 name" << std::endl;
            std::abort();
          }
        positions.resize(n_particles_cell);
        positions = {
          {.25, .25, .25} };
      } else if (structure == "BCC" ) {
        n_particles_cell = 2;
        if( types->size() != n_particles_cell )
          {
            lerr << "Parameter types for body centered cubic lattice must contain exactly 2 names" << std::endl;
            std::abort();
          }
        positions.resize(n_particles_cell);
        positions = {
          {.25, .25, .25} ,
          {.75, .75, .75} };          
      } else if (structure == "FCC" ) {
        n_particles_cell = 4;
        if( types->size() != n_particles_cell )
          {
            lerr << "Parameter types for face centered cubic lattice must contain exactly 4 names" << std::endl;
            std::abort();
          }
        positions.resize(n_particles_cell);
        positions = {
          {.25, .25, .25} ,
          {.25, .75, .75} ,
          {.75, .25, .75} ,
          {.75, .75, .25} };
      } else if (structure == "HCP" ) {
        n_particles_cell = 4;
        if( types->size() != n_particles_cell )
          {
            lerr << "Parameter types for hexagonal compact lattice must contain exactly 4 names" << std::endl;
            std::abort();
          }
        positions.resize(n_particles_cell);
        positions = {
	  {0.25000000,    0.25000000,    0.25000000} ,
	  {0.75000000,    0.75000000,    0.25000000} ,
	  {0.25000000,    0.58333333,    0.75000000} ,
	  {0.75000000,    0.08333333,    0.75000000} };

	size.y = 2. * size.y * sin(120. * M_PI / 180.);
	
	lout << "hcp cell = " << size << std::endl;
      } else if (structure == "DIAMOND" ) {
        n_particles_cell = 8;
        if( types->size() != n_particles_cell )
          {
            lerr << "Parameter types for diamond lattice must contain exactly 8 names" << std::endl;
            std::abort();
          }
        positions.resize(n_particles_cell);
        positions = {
          {.00, .00, .00} ,
          {.50, .50, .00} ,
          {.00, .50, .50} ,
          {.50, .00, .50} ,
          {.25, .25, .25} ,
          {.75, .75, .25} ,          
          {.75, .25, .75} ,
          {.25, .75, .75} };
      } else if (structure == "ROCKSALT" ) {
        n_particles_cell = 8;
        if( types->size() != n_particles_cell )
          {
            lerr << "Parameter types for rocksalt lattice must contain exactly 8 names" << std::endl;
            std::abort();
          }
        positions.resize(n_particles_cell);
        positions = {
          {.00, .00, .00} ,
          {.50, .50, .00} ,
          {.00, .50, .50} ,
          {.50, .00, .50} ,
          {.50, .00, .00} ,
          {.00, .50, .00} ,          
          {.00, .00, .50} ,
          {.50, .50, .50} };
      } else if (structure == "FLUORITE" ) {
        n_particles_cell = 12;
        if( types->size() != n_particles_cell )
          {
            lerr << "Parameter types for fluorite lattice must contain exactly 12 names" << std::endl;
            std::abort();
          }
        positions.resize(n_particles_cell);
        positions = {
          {.00, .00, .00} ,
          {.00, .50, .50} ,
          {.50, .00, .50} ,
          {.50, .50, .00} ,
          {.25, .25, .25} ,
          {.75, .25, .25} ,
          {.25, .75, .25} ,
          {.75, .75, .25} ,
          {.25, .25, .75} ,
          {.75, .25, .75} ,
          {.25, .75, .75} ,
          {.75, .75, .75} };
      } else if (structure == "C15" ) {
        n_particles_cell = 24;
        if( types->size() != n_particles_cell )
          {
            lerr << "Parameter types for C15 (lava) lattice must contain exactly 24 names" << std::endl;
            std::abort();
          }
        positions.resize(n_particles_cell);
        positions = {
          // We consider here the Cu2Mg structure
          // Cu atoms
          {.50, .50, .50} ,
          {.50, .75, .75} ,
          {.75, .50, .75} ,
          {.75, .75, .50} ,
          {.25, .50, .25} ,
          {.50, .25, .25} ,
          {.25, .25, .50} ,
          {.25, .75, .00} ,
          {.50, .00, .00} ,
          {.25, .00, .75} ,
          {.00, .50, .00} ,
          {.75, .25, .00} ,
          {.00, .25, .75} ,
          {.00, .75, .25} ,
          {.75, .00, .25} ,
          {.00, .00, .50} ,
          // Mg atoms are at positions of the diamond lattice
          {.125, .125, .125} ,
          {.875, .875, .875} ,
          {.875, .375, .375} ,
          {.375, .875, .375} ,
          {.625, .125, .625} ,
          {.375, .375, .875} ,
          {.125, .625, .625} ,          
          {.625, .625, .125} };
      } else if (structure == "PEROVSKITE" ) {
        n_particles_cell = 5;
        if( types->size() != n_particles_cell )
          {
            lerr << "Parameter types for perovskite lattice must contain exactly 5 names" << std::endl;
            std::abort();
          }
        positions.resize(n_particles_cell);
        positions = {
          {.50, .50, .50} ,
          {.00, .00, .00} ,
          {.50, .00, .00} ,
          {.00, .50, .00} ,
          {.00, .00, .50} };
      } else if (structure == "ST" ) {
        n_particles_cell = 1;
        if( types->size() != n_particles_cell )
          {
            lerr << "Parameter types for simple tetragonal lattice must contain exactly 1 name" << std::endl;
            std::abort();
          }
        positions.resize(n_particles_cell);
        positions = {
          {.25, .25, .25} };
      } else if (structure == "BCT" ) {
        n_particles_cell = 2;
        if( types->size() != n_particles_cell )
          {
            lerr << "Parameter types for body centered tetragonal  lattice must contain exactly 2 names" << std::endl;
            std::abort();
          }
        positions.resize(n_particles_cell);
        positions = {
        {.25, .25, .25},
        {.75, .75, .75} };          
      } else if (structure == "FCT" ) {
        n_particles_cell = 4;
        if( types->size() != n_particles_cell )
          {
            lerr << "Parameter types for face centered tetragonal lattice must contain exactly 4 names" << std::endl;
            std::abort();
          }
        positions.resize(n_particles_cell);
        positions = {
          {.25, .25, .25} ,
          {.25, .75, .75} ,
          {.75, .25, .75} ,
          {.75, .75, .25} };	
      }  else if (structure == "2BCT" ) {
        n_particles_cell = 4;
        if( types->size() != n_particles_cell )
          {
            lerr << "Parameter types for double body centered tetragonal lattice must contain exactly 4 names" << std::endl;
            std::abort();
          }
        positions.resize(n_particles_cell);
        positions = {
          {.00, .00, .00} ,
          {.00, .50, .25} ,
          {.50, .50, .50} ,
          {.50, .00, .75} };	
      }
      
      // optional shifting
      for(auto & r : positions) { r += *position_shift ; }

      std::map<std::string,unsigned int> typeMap;
      if( species.has_value() )
      {
        for(size_t i=0;i<species->size();i++)
        {
          typeMap[species->at(i).m_name]=i;
        }
      }

      // Get max and min positions
      // Need to define the size of the box
      // NOTE : only one processor need to do that      
      double box_size_x = repeats->i * size.x;
      double box_size_y = repeats->j * size.y;
      double box_size_z = repeats->k * size.z;
      Vec3d domain_size { box_size_x, box_size_y, box_size_z };
      
      if(rank==0) {
	
	AABB lattice_bounds  = { { 0., 0., 0. } , {box_size_x,box_size_y,box_size_z} };
	ldbg << "Lattice bounds      = "<<lattice_bounds<<std::endl;

	AABB computed_bounds = lattice_bounds;
	ldbg << "Computed_bounds  = " << computed_bounds << std::endl;

        compute_domain_bounds(domain,*bounds_mode,*enlarge_bounds,lattice_bounds,computed_bounds, *pbc_adjust_xform );
      }

      ldbg << "xform = : " << domain.xform() << std::endl;
      
      //send bounds and size_box values to all cores
      MPI_Bcast(&domain, sizeof(Domain), MPI_CHARACTER, 0, *mpi);
      assert( check_domain(domain) );

      // compute local processor's grid size and location so that cells are evenly distributed
      GridBlock in_block = { IJK{0,0,0} , domain.grid_dimension() };
      GridBlock out_block = simple_block_rcb( in_block, np, rank );
      ldbg<<"Domain = "<<domain<<std::endl;
      ldbg<<"In  block = "<< in_block << std::endl;
      ldbg<<"Out block = "<< out_block << std::endl;

      // initializes local processor's grid
      grid.set_offset( out_block.start );
      grid.set_origin( domain.bounds().bmin );
      grid.set_cell_size( domain.cell_size() );
      IJK local_grid_dim = out_block.end - out_block.start;
      grid.set_dimension( local_grid_dim );      

      // local processor's bounds
      AABB local_bounds = grid.grid_bounds();
      ldbg<<"local_bounds = "<< local_bounds << std::endl;
      ldbg<<"size = "<< size << std::endl;      
      Vec3d size_corr = domain.inv_xform() * size;
      ldbg<<"size_corr = "<< size_corr << std::endl;            

      ldbg << "local_bounds.bmin.x " << local_bounds.bmin.x << std::endl;
      ldbg << "local_bounds.bmax.x " << local_bounds.bmax.x << std::endl << std::endl;
      ldbg << "local_bounds.bmin.y " << local_bounds.bmin.y << std::endl;
      ldbg << "local_bounds.bmax.y " << local_bounds.bmax.y << std::endl << std::endl;
      ldbg << "local_bounds.bmin.z " << local_bounds.bmin.z << std::endl;
      ldbg << "local_bounds.bmax.z " << local_bounds.bmax.z << std::endl;

      // compute local processor's lattice portion
      ssize_t repeat_i_start = static_cast<ssize_t>( std::floor( local_bounds.bmin.x / size_corr.x ) ) ;
      ssize_t repeat_i_end   = static_cast<ssize_t>( std::ceil ( local_bounds.bmax.x / size_corr.x ) ) ;
      ssize_t repeat_j_start = static_cast<ssize_t>( std::floor( local_bounds.bmin.y / size_corr.y ) ) ;
      ssize_t repeat_j_end   = static_cast<ssize_t>( std::ceil ( local_bounds.bmax.y / size_corr.y ) ) ;
      ssize_t repeat_k_start = static_cast<ssize_t>( std::floor( local_bounds.bmin.z / size_corr.z ) ) ;
      ssize_t repeat_k_end   = static_cast<ssize_t>( std::ceil ( local_bounds.bmax.z / size_corr.z ) ) ;
      ldbg<<"lattice range = "<< IJK{repeat_i_start,repeat_j_start,repeat_k_start} << " - " << IJK{repeat_i_end,repeat_j_end,repeat_k_end} << std::endl;

      // =============== Section that concerns the porosity mode ========================

      std::random_device rd{};
      std::mt19937 gen{rd()};
      std::default_random_engine generator;
      std::normal_distribution<double> distribution_radius(*void_mean_diameter, *void_mean_diameter/2.);
      std::uniform_real_distribution<double> distribution_center_x(((*plane_offset) + 2.*(*void_mean_diameter)) * dot(*plane_normal,{1.,0.,0.}), domain_size.x);
      std::uniform_real_distribution<double> distribution_center_y(((*plane_offset) + 2.*(*void_mean_diameter)) * dot(*plane_normal,{0.,1.,0.}), domain_size.y);
      std::uniform_real_distribution<double> distribution_center_z(((*plane_offset) + 2.*(*void_mean_diameter)) * dot(*plane_normal,{0.,0.,1.}), domain_size.z);      
      
      std::vector<Vec3d> void_centers;
      std::vector<double> void_radiuses;
      
      if ((*geo_void == true) and ( *void_mode == "porosity" ))
      {
	      double calc_porosity = 0.;
	      double void_volume = 0.;
	      double volume = 1.0;
	      const auto a = column1( domain.xform() ); // m11, domain.xform().m21, domain.xform().m31 };
	      const auto b = column2( domain.xform() ); // { domain.xform().m12, domain.xform().m22, domain.xform().m32 };
	      const auto c = column3( domain.xform() ); // { domain.xform().m13, domain.xform().m23, domain.xform().m33 };
	      volume = dot( cross(a,b) , c );
	      volume *= bounds_volume( domain.bounds() );          
	      while ( calc_porosity < *void_porosity)
	      {
	        double x = distribution_center_x(generator);
	        double y = distribution_center_y(generator);
	        double z = distribution_center_z(generator);            
	        void_centers.push_back( {x, y, z} );
	        void_radiuses.push_back( distribution_radius(generator) );
	        void_volume += 4.*M_PI*pow(distribution_radius(generator),3.)/3.;
	        calc_porosity = void_volume/volume;
	      }
      }

      // =============== Generate particles into grid ========================

      std::vector<unsigned int> particle_type_id(n_particles_cell);
      for(size_t i=0;i<n_particles_cell;i++) { particle_type_id[i] = typeMap.at( types->at(i) ); }
      spin_mutex_array cell_locks;
      cell_locks.resize( grid.number_of_cells() );
      auto cells = grid.cells();

      ldbg << "total particles = "<< repeats->i * repeats->j * repeats->k * n_particles_cell<<std::endl;

#     pragma omp parallel for collapse(3)
      for (int i=repeat_i_start; i<repeat_i_end; i++)
	{
	  for (int j=repeat_j_start; j<repeat_j_end; j++)
	    {
	      for (int k=repeat_k_start; k<repeat_k_end; k++)
		{
		  size_t global_lattice_cell_index = i * repeats->j * repeats->k + j * repeats->k + k;
		  size_t global_particle_index = global_lattice_cell_index * n_particles_cell;
		  for (size_t l=0; l<n_particles_cell;l++)
		    {
		      Vec3d reduced_pos = {
			(i + positions[l].x ) ,
			(j + positions[l].y ) ,
			(k + positions[l].z ) };

		      Vec3d lattice = { size.x , size.y , size.z };

		      Vec3d pbc_adjust_invxform_lattice = domain.inv_xform() * lattice;		      
		      
		      Vec3d pos = {
			reduced_pos.x * pbc_adjust_invxform_lattice.x ,
			reduced_pos.y * pbc_adjust_invxform_lattice.y ,
			reduced_pos.z * pbc_adjust_invxform_lattice.z };			
		      
		      if (is_noise) {

			auto& re = onika::parallel::random_engine();
			std::normal_distribution<double> f_rand(0.,sigma_noise);
			pos.x += f_rand(re);
			pos.y += f_rand(re);
			pos.z += f_rand(re);
			
		      }
		      
		      /* Special geometry, compute flags to know if we keep or not the particle */
		      bool keep_particle=true;
		      if (*geo_plane == true or *geo_cylinder == true or *geo_void == true) {
		      	if (*geo_void == true and *void_mode == "porosity") {
		      	  live_or_die_void_porosity(domain.xform() * pos, domain_size, void_centers, void_radiuses, keep_particle);
		      	}
		      	if (*geo_void == true and *void_mode == "simple") {
		      	  live_or_die_void(domain.xform() * pos, domain_size, *void_center, *void_radius, keep_particle);
		      	}
		      	if (*geo_cylinder == true) {
		      	  live_or_die_cylinder(domain.xform() * pos, domain_size, *cylinder_inside, *cylinder_axis, *cylinder_radius, keep_particle);
		      	}						
		      	if (*geo_plane == true) {
			  if (*geo_void == true) {
			    live_or_die_plane_with_voids(domain.xform() * pos, *plane_normal, *plane_offset, keep_particle);
			  } else {
			    live_or_die_plane_without_voids(domain.xform() * pos, *plane_normal, *plane_offset, keep_particle);
			  }
		      	}
		      }
		      
		      IJK loc = grid.locate_cell(pos); //domain_periodic_location( domain , pos );
		      if( loc.i>=0 && loc.i<local_grid_dim.i &&
			  loc.j>=0 && loc.j<local_grid_dim.j &&
			  loc.k>=0 && loc.k<local_grid_dim.k &&
			  pos.x>=0.0 && pos.x<domain_size.x && 
			  pos.y>=0.0 && pos.y<domain_size.y && 
			  pos.z>=0.0 && pos.z<domain_size.z )
			{
			  if (*geo_plane == true or *geo_cylinder == true or *geo_void == true) {
			    if( keep_particle == true )
			      {
			  	uint64_t id = global_particle_index + l;
			  	unsigned int type_id = particle_type_id[l];
			  	assert( grid.contains(loc) );
			  	assert( min_distance2_between( pos, grid.cell_bounds(loc) ) < grid.epsilon_cell_size2() );
			  	ParticleTuple t = ParticleTupleIO(pos.x,pos.y,pos.z,id,type_id);
			  	size_t cell_i = grid_ijk_to_index(local_grid_dim, loc);
                      
			  	cell_locks[cell_i].lock();
			  	cells[cell_i].push_back( t );
			  	cell_locks[cell_i].unlock();
			      }
			  } else {
			    uint64_t id = global_particle_index + l;
			    unsigned int type_id = particle_type_id[l];
			    assert( grid.contains(loc) );
			    assert( min_distance2_between( pos, grid.cell_bounds(loc) ) < grid.epsilon_cell_size2() );
			    ParticleTuple t = ParticleTupleIO(pos.x,pos.y,pos.z,id,type_id);
			    size_t cell_i = grid_ijk_to_index(local_grid_dim, loc);
			  
			    cell_locks[cell_i].lock();
			    cells[cell_i].push_back( t );
			    cell_locks[cell_i].unlock();
			  }
			}
		    }
		}
	    }
	}
      
      grid.rebuild_particle_offsets();
      
    }

  private:

    static inline void live_or_die_plane_without_voids(const Vec3d pos, const Vec3d N, const double offset, bool& keep)
    {
      double dist = (dot(pos, N) - offset)/norm(N);
      if (dist > 0.) keep = false;
    }

    static inline void live_or_die_plane_with_voids(const Vec3d pos, const Vec3d N, const double offset, bool& keep)
    {
      double dist = (dot(pos, N) - offset)/norm(N);
      if (dist < 0.) keep = true;
    }
    
    static inline void live_or_die_cylinder(const Vec3d pos, const Vec3d domain_size, const bool inside, const Vec3d N, const double cylinder_radius, bool& keep)
    {
      double dist = norm(cross(pos-domain_size/2., N))/norm(N);
      if ((inside == true) && (dist > cylinder_radius)) {
	keep = false;
      } else if ((inside == false) && (dist < cylinder_radius)) {
	keep = false;
      }
    }

    static inline void live_or_die_void(const Vec3d pos, const Vec3d domain_size, const Vec3d void_center, const double void_radius, bool& keep)
    {
      Vec3d dist_vec = pos-void_center;
      if (dist_vec.x > domain_size.x/2.) dist_vec.x -= domain_size.x;
      if (dist_vec.y > domain_size.y/2.) dist_vec.y -= domain_size.y;
      if (dist_vec.z > domain_size.z/2.) dist_vec.z -= domain_size.z;
      double dist = norm(dist_vec);      
      if (dist < void_radius) keep = false;
    }

    static inline void live_or_die_void_porosity(const Vec3d pos, const Vec3d domain_size, const std::vector<Vec3d> void_center, const std::vector<double> void_radius, bool& keep)
    {
      for (size_t i=0; i < void_center.size(); i++) {
        Vec3d dist_vec = pos-void_center[i];
        if (dist_vec.x > domain_size.x/2.) dist_vec.x -= domain_size.x;
        if (dist_vec.y > domain_size.y/2.) dist_vec.y -= domain_size.y;
        if (dist_vec.z > domain_size.z/2.) dist_vec.z -= domain_size.z;
        double dist = norm(dist_vec);
        if (dist < void_radius[i]) {
          keep = false;
          break;
        }
      }
    }
    
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(orthorhombic_lattice)
  {
    OperatorNodeFactory::instance()->register_factory("orthorhombic_lattice", make_grid_variant_operator< OrthorhombicLatticeOperator >);
  }

}
