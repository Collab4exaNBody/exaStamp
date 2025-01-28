#include <chrono>
#include <ctime>
#include <mpi.h>
#include <string>
#include <numeric>

#include <onika/math/basic_types_yaml.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <onika/math/basic_types_stream.h>
#include <onika/log.h>
//#include "exanb/vector_utils.h"
#include <exanb/core/file_utils.h>
#include <exanb/core/domain.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/check_particles_inside_cell.h>
#include <onika/parallel/random.h>

namespace exaStamp
{
  using namespace exanb;

  // Read XYZ files.
  // This files must be ASCII files
  // first line : number of particles
  // second line : a comment, not read by the code
  // next lines : type X Y Z
  // example : CHNO2 5.3 9.23 -3.5
  // Note : the type have to match a specie name

  template<typename GridT>
  class ReadXYZwXFormNode : public OperatorNode
  {
    ADD_SLOT( MPI_Comm        , mpi          , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( std::string     , file         , INPUT , REQUIRED );
    ADD_SLOT( Domain          , domain       , INPUT_OUTPUT );
    ADD_SLOT( GridT           , grid         , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies , species      , INPUT ); // optional. if no species given, type ids are allocated automatically
    ADD_SLOT( ReadBoundsSelectionMode, bounds_mode   , INPUT , ReadBoundsSelectionMode::FILE_BOUNDS );
    ADD_SLOT( bool             , gaussian_noise        , INPUT , false);    
    ADD_SLOT( double           , gaussian_noise_sigma  , INPUT , 0.0);

  public:
    inline void execute () override final
    {
      // static const double coord_conv = UnityConverterHelper::convert(1.0, "ang"); // conversion factor between input unity to exaStamp internal unity

      //-------------------------------------------------------------------------------------------
      // Reading datas from YAML or previous input
      std::string file_name = data_file_path( *file );
      Domain& domain = *(this->domain);
      GridT& grid = *(this->grid);

      bool is_noise = *(this->gaussian_noise);            
      double sigma_noise = *(this->gaussian_noise_sigma);
            
      //-------------------------------------------------------------------------------------------
      std::string basename;
      std::string::size_type p = file_name.rfind("/");
      if( p != std::string::npos ) basename = file_name.substr(p+1);
      else basename=file_name;      
      lout << "======== " << basename << " ========" << std::endl;
      //-------------------------------------------------------------------------------------------

      using ParticleTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_id, field::_type>;
      using ParticleTuple = decltype( grid.cells()[0][0] );

      assert( grid.number_of_particles() == 0 );

      // MPI Initialization
      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      //useful later
      std::string line;
      uint64_t n_particles = 0;

      std::map<std::string,unsigned int> typeMap;
      unsigned int nextTypeId = 0;
      if( species.has_value() )
      {
        for(size_t i=0;i<species->size();i++)
        {
          typeMap[species->at(i).m_name]=i;
        }
        nextTypeId = species->size();
      }
           
      std::vector<ParticleTupleIO> particle_data;

      // Get max and min positions
      // Need to define the size of the box
      // NOTE : only one processor need to do that
      bool uniform_scale = false;
      Mat3d domain_xform = make_identity_matrix();
      if(rank==0)
      {
        //Open xyz file
        std::ifstream file;
        file.open(file_name, std::ifstream::in);
        if(!file.is_open())
        {
          lerr << "Error in reading xyz : file "<< file_name << " not found !" << std::endl;
          std::abort();
        }

        // Read number of atoms
        ssize_t n_atoms = -1;
        getline(file,line);
        std::stringstream(line) >> n_atoms;

        double box_size_x=-1.0;
        double box_size_y=-1.0;
        double box_size_z=-1.0;
        
        Mat3d H = make_identity_matrix();	
        std::getline(file,line);

	      unsigned quote;
	      quote = line.find("\"");
	      quote += 1;
	      std::string secondPart = line.substr(quote);
	      quote = secondPart.find("\"");
	      std::string lineclean = secondPart.substr(0,quote);

        std::stringstream(line) >> H.m11 >> H.m12 >> H.m13 >> H.m21 >> H.m22 >> H.m23 >> H.m31 >> H.m32 >> H.m33;

        lout << "H = " << H << std::endl;

        Vec3d a = Vec3d{H.m11, H.m12, H.m13};
        Vec3d b = Vec3d{H.m21, H.m22, H.m23};
        Vec3d c = Vec3d{H.m31, H.m32, H.m33};

        box_size_x = norm(a);
        box_size_y = norm(b);
        box_size_z = norm(c);

        lout << "a = " << box_size_x << std::endl;        
        lout << "b = " << box_size_y << std::endl;        
        lout << "c = " << box_size_z << std::endl;        

        Mat3d Ht = transpose(H);

        // read one line at a time
        while (std::getline(file,line))
        {
          std::string type;
          double x=0.0, y=0.0, z=0.0;

          //first value not necessary here
          std::stringstream(line) >> type >> x >> y >> z;

	        Vec3d r{x, y, z};

	        r = inverse(Ht) * r;

	        if (r.x < 0.) r.x += 1.;
	        if (r.y < 0.) r.y += 1.;
	        if (r.z < 0.) r.z += 1.;

	        if (r.x >= 1.) r.x -= 1.;
	        if (r.y >= 1.) r.y -= 1.;
	        if (r.z >= 1.) r.z -= 1.;

                x = box_size_x * r.x;
                y = box_size_y * r.y;
                z = box_size_z * r.z;

	        if (is_noise) {

	          auto& re = rand::random_engine();
	          std::normal_distribution<double> f_rand(0.,sigma_noise);
	          x += f_rand(re);
	          y += f_rand(re);
	          z += f_rand(re);
			      
	        }
	  
          if( typeMap.find(type) == typeMap.end() )
          {
            typeMap[type] = nextTypeId;
            ++ nextTypeId;
          }
          
          particle_data.push_back( ParticleTupleIO(x,y,z,n_particles++,typeMap[type]) );
        }
        
        if( ! domain.xform_is_identity() )
        {
          lerr << "needs initial XForm, resetting XForm"<<std::endl;
          domain.set_xform( make_identity_matrix() );
        }

        AABB file_bounds  = { { 0., 0., 0. } , {box_size_x,box_size_y,box_size_z} };
        compute_domain_bounds(domain,*bounds_mode,0.0, file_bounds ,file_bounds, true );
        
        if( !domain.xform_is_identity() )
        {
          Mat3d inv_xform = domain.inv_xform();
          for( auto& p : particle_data )
          {
            Vec3d r = inv_xform * Vec3d{ p[field::rx] , p[field::ry] , p[field::rz] };
            p[field::rx] = r.x;
            p[field::ry] = r.y;
            p[field::rz] = r.z;            
          }
        }

        const Mat3d D = diag_matrix(Vec3d{box_size_x, box_size_y, box_size_z});
        const Mat3d Hbis = Ht * inverse(D);
        domain_xform = Hbis * domain.xform();
        
        uniform_scale = is_uniform_scale( domain_xform );
        lout << "Uniform scale    = "<<std::boolalpha<<uniform_scale << std::endl;
        if( uniform_scale )
        {
          domain.set_xform( make_identity_matrix() );
          domain.set_cell_size( domain.cell_size() * domain_xform.m11 );
          domain.set_bounds( { domain.origin() * domain_xform.m11 , domain.extent() * domain_xform.m11 } );
        }
        else
        {
          domain.set_xform( domain_xform );
        }
        
        lout << "Particles        = "<<particle_data.size()<<std::endl;
        lout << "Domain XForm     = "<<domain.xform()<<std::endl;
        lout << "Domain bounds    = "<<domain.bounds()<<std::endl;
        lout << "Domain size      = "<<bounds_size(domain.bounds()) <<std::endl;
        lout << "Real size        = "<<bounds_size(domain.bounds()) * Vec3d{domain.xform().m11,domain.xform().m22,domain.xform().m33} <<std::endl;
        lout << "Cell size        = "<<domain.cell_size()<<std::endl;
        lout << "Grid dimensions  = "<<domain.grid_dimension()<<" ("<<grid_cell_count(domain.grid_dimension())<<" cells)"<< std::endl;
      }

      //send bounds and size_box values to all cores
      MPI_Bcast(&domain, sizeof(Domain), MPI_CHARACTER, 0, *mpi);
      assert( check_domain(domain) );

      grid.set_offset( IJK{0,0,0} );
      grid.set_origin( domain.bounds().bmin );
      grid.set_cell_size( domain.cell_size() );
      grid.set_dimension( domain.grid_dimension() );

      if(rank==0)
      {
        for( auto p : particle_data )
        {
          Vec3d r = { p[field::rx] , p[field::ry] , p[field::rz] };
          if( uniform_scale ) r = r * domain_xform.m11;
          IJK loc = domain_periodic_location( domain , r ); 
          assert( grid.contains(loc) );
          assert( min_distance2_between( r, grid.cell_bounds(loc) ) < grid.epsilon_cell_size2() );
          p[field::rx] = r.x;
          p[field::ry] = r.y;
          p[field::rz] = r.z;
          ParticleTuple t = p;
          grid.cell( loc ).push_back( t );
        }
      }

      lout << "============================" << std::endl;
     
      grid.rebuild_particle_offsets();

#     ifndef NDEBUG
      bool particles_inside_cell = check_particles_inside_cell(grid);
      assert( particles_inside_cell );
#     endif
    }
    
  };

  // === register factories ===
  __attribute__((constructor)) static void register_factories()
  {
    OperatorNodeFactory::instance()->register_factory("read_xyz_file_with_xform", make_grid_variant_operator< ReadXYZwXFormNode >);
  }

}
