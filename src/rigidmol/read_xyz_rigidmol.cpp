#include <chrono>
#include <ctime>
#include <mpi.h>
#include <string>
#include <numeric>

#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/log.h>
//#include "exanb/vector_utils.h"
#include <exanb/core/file_utils.h>
#include <exanb/core/domain.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/check_particles_inside_cell.h>

#include <exanb/core/quaternion_operators.h>

namespace exaStamp
{

  // Read XYZ files.
  // This files must be ASCII files
  // first line : number of particles
  // second line : a comment, not read by the code
  // next lines : type X Y Z
  // example : CHNO2 5.3 9.23 -3.5
  // Note : the type have to match a specie name

  template<class ParticleTupleIO>
  static inline AABB read_xyz_particles(
    const std::string& dir_name,
    const std::string& file_name,
    std::unordered_map<std::string,unsigned int>& typeMap,
    unsigned int& nextTypeId,
    ParticleSpecies& species,
    std::vector<ParticleTupleIO>& particle_data,
    size_t& n_particles,
    double& min_x, double& min_y, double& min_z, 
    double& max_x, double& max_y, double& max_z )
  {
    static const double coord_conv = UnityConverterHelper::convert(1.0, "ang"); // conversion factor between input unity to exaStamp internal unity
    using MoleculeTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_type>;
    using LocalParticleTuple = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_id, field::_type, field::_orient >;

    ldbg << "read "<<file_name <<std::endl;
    
    std::ifstream file;
    file.open( data_file_path(file_name), std::ifstream::in);
    if(!file.is_open())
    {
      lerr << "Error in reading xyz : file "<< file_name << " not found !" << std::endl;
      std::abort();
    }

    // line buffer
    std::string line;

    // Read number of atoms
    getline(file,line);
    ssize_t n_atoms = -1;
    std::stringstream(line) >> n_atoms;

    std::getline(file,line);
    ldbg << "File comment = "<<line<<std::endl;
    ldbg << "N Atoms = "<<n_atoms<<std::endl;

    // read one line at a time
    while ( std::getline(file,line) )
    {
      std::string type;
      double x=0.0, y=0.0, z=0.0;
      Quaternion Q = { 0., 0., 0., 0. };

      //first value not necessary here
      std::stringstream line_stream(line);
      line_stream >> type >> x >> y >> z;

      if( ! type.empty() )
      {
        x *= coord_conv; // XYZ files use angstrom
        min_x = std::min( min_x , x );
        max_x = std::max( max_x , x );

        y *= coord_conv;
        min_y = std::min( min_y , y );
        max_y = std::max( max_y , y );

        z *= coord_conv;
        min_z = std::min( min_z , z );
        max_z = std::max( max_z , z );
        
        auto type_it = typeMap.find(type);
        if( type_it == typeMap.end() )
        {
          line_stream >> Q.w >> Q.x >> Q.y >> Q.z ;
          std::string sub_file = dir_name + type + ".xyz";
          ldbg <<"'"<< type << "' not found, implicitly a rigid molecule name, read " << sub_file << std::endl;
          std::vector<MoleculeTupleIO> molecule_data;
          double mol_min_x = std::numeric_limits<double>::max();
          double mol_min_y = std::numeric_limits<double>::max();
          double mol_min_z = std::numeric_limits<double>::max();
          double mol_max_x = std::numeric_limits<double>::lowest();
          double mol_max_y = std::numeric_limits<double>::lowest();
          double mol_max_z = std::numeric_limits<double>::lowest();
          size_t mol_atoms = 0;
          read_xyz_particles( dir_name, sub_file, typeMap, nextTypeId, species, molecule_data, mol_atoms, mol_min_x,mol_min_y,mol_min_z, mol_max_x,mol_max_y,mol_max_z );
          ParticleSpecie mol_specy;
          mol_specy.set_name( type );
          mol_specy.m_mass = 0.0;
          mol_specy.m_minert = {0.0,0.0,0.0};
          mol_specy.m_charge = 0.0;
          mol_specy.m_z = 0;
          mol_specy.m_rigid_atom_count = molecule_data.size();
          ldbg << "molecule "<<mol_specy.m_name<<" has "<<mol_specy.m_rigid_atom_count<<" atoms"<<std::endl;
          assert( mol_specy.m_rigid_atom_count <= MAX_RIGID_MOLECULE_ATOMS );
          //mol_specy.m_rigid_atom_names.resize(mol_specy.m_rigid_atom_count);
          for(unsigned int i=0;i<mol_specy.m_rigid_atom_count;i++)
          {
            mol_specy.m_rigid_atoms[i].m_pos = { molecule_data[i][field::rx] , molecule_data[i][field::ry] , molecule_data[i][field::rz] };
            mol_specy.m_rigid_atoms[i].m_atom_type = molecule_data[i][field::type];
            if(mol_specy.m_rigid_atoms[i].m_atom_type<0 || size_t(mol_specy.m_rigid_atoms[i].m_atom_type)>=species.size() )
            {
              lerr << "bad rigid atom type found ("<<mol_specy.m_rigid_atoms[i].m_atom_type<<")"<<std::endl;
              std::abort();
            }
            mol_specy.set_rigid_atom_name( i , species.at(mol_specy.m_rigid_atoms[i].m_atom_type).name() );
            mol_specy.m_mass += species.at(mol_specy.m_rigid_atoms[i].m_atom_type).m_mass;
          }
          species.push_back( mol_specy );
          typeMap[type] = nextTypeId;
          ++ nextTypeId;
          assert( species.size() == nextTypeId );
        }
        else if( species.at(type_it->second).m_rigid_atom_count > 1 )
        {
          ldbg << "read orientation for molecule "<<type<<std::endl;
          line_stream >> Q.w >> Q.x >> Q.y >> Q.z ;      
        }
        
        
//fprintf(stdout," RIGIDMOL_normalize_quat:     OR_0(%14.10lf %14.10lf %14.10lf %14.10lf)",Q.w,Q.x,Q.y,Q.z);//NP_TEST
Q = normalize(Q);//NP_TEST
//fprintf(stdout," => OR_1(%14.10lf %14.10lf %14.10lf %14.10lf)\n",Q.w,Q.x,Q.y,Q.z);//NP_TEST
        LocalParticleTuple t(x,y,z,n_particles++,typeMap[type],Q);
        particle_data.push_back( t );
      }      
    }
    return { { min_x, min_y, min_z } , { max_x, max_y, max_z } };
  }

  template<typename GridT>
  class ReadXYZRigidMol : public OperatorNode
  {
    ADD_SLOT( MPI_Comm        , mpi          , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( std::string     , file         , INPUT , REQUIRED );
    ADD_SLOT( ReadBoundsSelectionMode, bounds_mode   , INPUT , ReadBoundsSelectionMode::FILE_BOUNDS );
    ADD_SLOT( Domain          , domain       , INPUT_OUTPUT );
    ADD_SLOT( GridT           , grid         , INPUT_OUTPUT );
    ADD_SLOT( double          , enlarge_bounds, INPUT , 0.0 );
    ADD_SLOT( bool            , pbc_adjust_xform , INPUT , false );
    ADD_SLOT( ParticleSpecies , species      , INPUT_OUTPUT ); // optional. if no species given, type ids are allocated automatically

  public:
    inline void execute () override final
    {

      //-------------------------------------------------------------------------------------------
      // Reading datas from YAML or previous input
      std::string file_name = *file;
      Domain& domain = *(this->domain);
      GridT& grid = *(this->grid);
      
      if( *pbc_adjust_xform )
      {
        if( ! domain.xform_is_identity() )
        {
          lerr << "pbc_adjust_xform needs initial XForm to be identity, resetting XForm"<<std::endl;
          domain.set_xform( make_identity_matrix() );
        }
      }

//      grid = GridT();
      //-------------------------------------------------------------------------------------------
      std::string basename;
      std::string dirname = "";
      std::string::size_type p = file_name.rfind("/");
      if( p != std::string::npos )
      {
        dirname = file_name.substr(0,p+1);
        basename = file_name.substr(p+1);
      }
      else
      {
        basename = file_name;
      }
      lout << "======== " << basename << " ========" << std::endl;
      //-------------------------------------------------------------------------------------------

      using ParticleTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_id, field::_type, field::_orient >;
      using ParticleTuple = decltype( grid.cells()[0][0] );

      assert( grid.number_of_particles() == 0 );

      // MPI Initialization
      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      uint64_t n_particles = 0;

      std::unordered_map<std::string,unsigned int> typeMap;
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
      if(rank==0)
      {
        // We need to define limits of the domain from the .xyz file
        // and the position of particles
        // For that, we check the max and min position from all particles
        double min_x = std::numeric_limits<double>::max();
        double min_y = std::numeric_limits<double>::max();
        double min_z = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        double max_y = std::numeric_limits<double>::lowest();
        double max_z = std::numeric_limits<double>::lowest();

        //Open xyz file
        AABB file_bounds = read_xyz_particles(dirname, file_name, typeMap, nextTypeId, *species, particle_data, n_particles, min_x, min_y, min_z, max_x, max_y, max_z);

        ldbg << "min position xyz file : " << min_x << " " << min_y << " " << min_z << std::endl;
        ldbg << "max position xyz file : " << max_x << " " << max_y << " " << max_z << std::endl;

        //DOMAIN
        AABB computed_bounds = { {min_x, min_y, min_z} , {max_x, max_y, max_z} };
        ldbg << "computed_bounds  = " << computed_bounds << std::endl;
 
        //domain.m_bounds = bounds;
        compute_domain_bounds(domain,*bounds_mode,*enlarge_bounds,file_bounds,computed_bounds, *pbc_adjust_xform );
        if( *pbc_adjust_xform && !domain.xform_is_identity() )
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
          Vec3d r{ p[field::rx] , p[field::ry] , p[field::rz] };
          IJK loc = domain_periodic_location( domain , r ); //grid.locate_cell(r);
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
    OperatorNodeFactory::instance()->register_factory("read_xyz_rigidmol", make_grid_variant_operator< ReadXYZRigidMol >);
  }

}
