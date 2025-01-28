

#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/grid_fields.h>
#include <onika/math/basic_types_stream.h>
#include <onika/log.h>
#include <exanb/core/unityConverterHelper.h>

#include <exanb/io/mpi_file_io.h>
#include <exaStamp/molecule/stampv4_io.h>
#include <exanb/mpi/all_value_equal.h>
#include <onika/oarray.h>

#include <onika/math/basic_types_yaml.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>

//#include "exanb/vector_utils.h"
#include <exanb/core/file_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/molecule/molecule_species.h>

#include <exanb/core/check_particles_inside_cell.h>
#include <onika/parallel/random.h>
#include <exanb/core/math_utils.h>
#include <onika/math/basic_types_operators.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <string>
#include <limits>
#include <memory>
#include <map>

#include <ctime>
#include <mpi.h>
#include <string>
#include <numeric>

namespace exaStamp
{
  using namespace exanb;

  // Read XYZ files.
  // This files must be ASCII files
  // first line : number of particles
  // second line : xform of simulation domain
  // next lines : typeFF X Y Z molid moltype connectivity[5]
  // example : Cno2 5.3 9.23 -3.5 1 0 6 2 7 -1 -1
  // Note : the typeFF have to match a specie name

  template<typename GridT>
  class ReadXYZMolecules : public OperatorNode
  {
    ADD_SLOT( MPI_Comm        , mpi          , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( std::string     , filename     , INPUT , REQUIRED );
    ADD_SLOT( Domain          , domain       , INPUT_OUTPUT );
    ADD_SLOT( GridT           , grid         , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies , species      , INPUT ); // optional. if no species given, type ids are allocated automatically
    ADD_SLOT( MoleculeSpeciesVector , molecules , INPUT_OUTPUT, MoleculeSpeciesVector{} , DocString{"Molecule descriptions"} );

  public:
    inline void execute () override final
    {
      //-------------------------------------------------------------------------------------------
      std::string file_name = data_file_path( *filename );
      std::string basename;
      std::string::size_type p = file_name.rfind("/");
      if( p != std::string::npos ) basename = file_name.substr(p+1);
      else basename = file_name;      
      lout << "======== " << basename << " ========" << std::endl;
      //-------------------------------------------------------------------------------------------

      using MoleculeTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_vx, field::_vy, field::_vz, field::_id, field::_type, field::_idmol, field::_cmol, field::_charge>;
      //using MoleculeTuple = decltype( grid.cells()[0][0] );      
      assert( grid->number_of_particles() == 0 );

      // MPI Initialization
      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);
      
      // converts type name to type index
      auto type_name_to_index = [&]( const std::string& typeAtome ) -> int 
      {
        int nSpecies = species->size();
        int at_type = 0;
        for(at_type = 0; at_type < nSpecies ; ++at_type)
        {
          if( species->at(at_type).name() == typeAtome ) return at_type;
        }
        lerr << "Atom of type " << typeAtome << " is unknown, the "<<nSpecies<<" available types are :" << std::endl;
        for(at_type = 0; at_type < nSpecies ; ++at_type) lerr << "\t" << species->at(at_type).name()<<std::endl << std::flush;
        std::abort();
        return -1;
      };

      // 1. Reads header, setup domain, and broadcast domain information to others MPI procs
      uint64_t count = 0;      
      if(rank==0)
      {
        //Open xyz file
        std::ifstream file(file_name);
        if( ! file )
        {
          fatal_error() << "Error in reading xyz : file "<< file_name << " not found !" << std::endl;
        }

        // Read number of atoms
        file >> count;
        ldbg << "count = "<<count<<std::endl;

        Mat3d H = make_identity_matrix();         
        file >> H.m11 >> H.m12 >> H.m13 >> H.m21 >> H.m22 >> H.m23 >> H.m31 >> H.m32 >> H.m33;
        ldbg << "H = " << H << std::endl;
        double smin=1.0, smax=1.0;
        matrix_scale_min_max(H,smin,smax);
        ldbg << "scale = "<<smin<<" / "<<smax<<std::endl;

        const Vec3d a = Vec3d{H.m11, H.m12, H.m13};
        const Vec3d b = Vec3d{H.m21, H.m22, H.m23};
        const Vec3d c = Vec3d{H.m31, H.m32, H.m33};
        const double Lx = norm(a);
        const double Ly = norm(b);
        const double Lz = norm(c);
        ldbg << "a = " << a << " norm=" << Lx << std::endl;        
        ldbg << "b = " << b << " norm=" << Ly << std::endl;        
        ldbg << "c = " << c << " norm=" << Lz << std::endl;
        const double angle_a_b = (180/M_PI) * acos( dot(a,b) / ( Lx*Ly ) );
        const double angle_a_c = (180/M_PI) * acos( dot(a,c) / ( Lx*Lz ) );
        const double angle_b_c = (180/M_PI) * acos( dot(b,c) / ( Ly*Lz ) );
        ldbg << "angle a.b = " << angle_a_b << std::endl;        
        ldbg << "angle a.c = " << angle_a_c << std::endl;        
        ldbg << "angle b.c = " << angle_b_c << std::endl;        
        
        if( is_diagonal(H) && Lx==Ly && Ly==Lz )
        {
          ldbg << "cubic & diagonal H matrix" << std::endl;
        }
        
        const Mat3d inv_H = inverse(H);
        ldbg << "H^-1 = " << inv_H << std::endl;

        // read one line at a time
        std::map< std::string , int > molecule_name_map;
        std::string typeMol;
        std::string typeAtom;
        std::vector<MoleculeTupleIO> atom_data;
        atom_data.reserve( count );

        AABB phys_bbox = { {0.0,0.0,0.0} , {0.0,0.0,0.0} };
        
        for(uint64_t cnt=0;cnt<count;cnt++)
        {
          int64_t at_id = -1;
          int64_t moleculeid = -1;
          double x=0.0, y=0.0, z=0.0;
          int64_t c0=-1 , c1=-1 , c2=-1 , c3=-1 , c4=-1;
          double charge = 0.0;

          assert( ! file.eof() && file.good() );
          file >> typeMol >> typeAtom >> moleculeid >> at_id >> x >> y >> z >> c0 >> c1 >> c2 >> c3 >> c4 >> charge;
          //lout <<"eof="<<file.eof()<<" good="<<file.good()<<" typeMol="<<typeMol<<" typeAtom="<<typeAtom<<" moleculeid="<<moleculeid<<" at_id="<<at_id<<" r="<<x<<","<<y<<","<<z<<" con="<<c0<<","<<c1<<","<<c2<<","<<c3 << std::endl;

          if( c4 != -1 )
          {
            fatal_error() << "atom connectivity with more than 4 bonds is not supported" << std::endl;
          }

          // convert position to grid's base
          Vec3d r{x, y, z};
          if( cnt == 0 )
          {
            phys_bbox.bmin = phys_bbox.bmax = r;
          }
          else
          {
            phys_bbox.bmin = min( phys_bbox.bmin , r );
            phys_bbox.bmax = max( phys_bbox.bmax , r );
          }
          r = inv_H * Vec3d{ x, y, z };

          auto it = molecule_name_map.find( typeMol );
          int itypeMol = -1;
          if( it != molecule_name_map.end() )
          {
            itypeMol = it->second;
          }
          else
          {
            itypeMol = molecule_name_map.size();
            molecule_name_map.insert( { typeMol , itypeMol } );
          }
          const int itypeAtom = type_name_to_index(typeAtom);

          const Vec3d v = { 0., 0., 0. }; // not read from file by now
          MoleculeTupleIO tp( r.x, r.y, r.z , v.x, v.y, v.z, at_id, itypeAtom, 0, std::array<uint64_t,4>{uint64_t(c0),uint64_t(c1),uint64_t(c2),uint64_t(c3)} , charge );
          atom_data.push_back( tp );     
        }
        file.close();

        ldbg << "phys_bbox = " << phys_bbox << std::endl;    

        AABB unit_bbox = { {0.0,0.0,0.0} , {0.0,0.0,0.0} };
        if( ! atom_data.empty() ) unit_bbox.bmin = unit_bbox.bmax = Vec3d{ atom_data[0][field::rx] , atom_data[0][field::ry] , atom_data[0][field::rz] };
        for(const auto& tp : atom_data)
        {
          Vec3d r = { tp[field::rx] , tp[field::ry] , tp[field::rz] };
          unit_bbox.bmin = min( unit_bbox.bmin , r );
          unit_bbox.bmax = max( unit_bbox.bmax , r );
        }
        
        const Vec3d unit_size = unit_bbox.bmax - unit_bbox.bmin;
        ldbg << "unit_bbox=" << unit_bbox << " , unit_size="<<unit_size<<std::endl;
        const double cell_size = domain->cell_size() > 0.0 ? domain->cell_size() : 8.0 ;
        domain->set_cell_size( cell_size );
        const IJK grid_dims = { static_cast<ssize_t>(ceil(Lx/cell_size)) , static_cast<ssize_t>(ceil(Ly/cell_size)) , static_cast<ssize_t>(ceil(Lz/cell_size)) };
        domain->set_grid_dimension( grid_dims );
        ldbg << "domain grid = "<<grid_dims<< std::endl;
        Vec3d dom_size = grid_dims * cell_size;
        ldbg << "domain size = "<< dom_size << std::endl;
        
        domain->set_bounds( { {0.0,0.0,0.0} , dom_size } );
        const Mat3d dom_H = { dom_size.x , 0.0 , 0.0
                            , 0.0 , dom_size.y , 0.0
                            , 0.0 , 0.0 , dom_size.z };
        domain->set_xform( H * inverse( dom_H ) );
        ldbg << "domain = " << *domain << std::endl;

        const Mat3d dom_box = domain->xform() * dom_H;

        const Vec3d dom_a = { dom_box.m11, dom_box.m12, dom_box.m13 };
        const Vec3d dom_b = { dom_box.m21, dom_box.m22, dom_box.m23 };
        const Vec3d dom_c = { dom_box.m31, dom_box.m32, dom_box.m33 };
        const double dom_Lx = norm(dom_a);
        const double dom_Ly = norm(dom_b);
        const double dom_Lz = norm(dom_c);
        ldbg << "dom_a = " << dom_a << " , norm = " << dom_Lx << std::endl;        
        ldbg << "dom_b = " << dom_b << " , norm = " << dom_Ly << std::endl;        
        ldbg << "dom_c = " << dom_c << " , norm = " << dom_Lz << std::endl;        
        ldbg << "dom angle a.b = " << (180/M_PI) * acos( dot(dom_a,dom_b) / ( dom_Lx*dom_Ly ) ) << std::endl;        
        ldbg << "dom angle a.c = " << (180/M_PI) * acos( dot(dom_a,dom_c) / ( dom_Lx*dom_Lz ) ) << std::endl;        
        ldbg << "dom angle b.c = " << (180/M_PI) * acos( dot(dom_b,dom_c) / ( dom_Ly*dom_Lz ) ) << std::endl;        
                            
        const Mat3d inv_xform = domain->inv_xform();
        ldbg << "file to grid inv_xform = " << inv_xform << std::endl;        
        
        grid->set_offset( IJK{0,0,0} );
        grid->set_origin( domain->bounds().bmin );
        grid->set_cell_size( domain->cell_size() );
        grid->set_dimension( domain->grid_dimension() );

        AABB dom_bbox = { {0.0,0.0,0.0} , {0.0,0.0,0.0} }; bool dom_bbox_init=true;
        for(auto& tp : atom_data)
        {
          Vec3d r = { tp[field::rx] , tp[field::ry] , tp[field::rz] };
          r = r * dom_size;
          IJK loc = domain_periodic_location( *domain, r ); // grid->locate_cell( r );
          tp[field::rx] = r.x; tp[field::ry] = r.y; tp[field::rz] = r.z;
          if( dom_bbox_init )
          {
            dom_bbox_init = false;
            dom_bbox.bmin = dom_bbox.bmax = r;
          }
          else
          {
            dom_bbox.bmin = min( dom_bbox.bmin , r );
            dom_bbox.bmax = max( dom_bbox.bmax , r );
          }
          const uint64_t at_id = tp[field::id];
          if( ! grid->contains(loc) )
          {
            fatal_error() <<"particle #"<<at_id<<", ro="<<r<<", r="<<r<<"<< in cell "<<loc<<" not in grid : offset="<<grid->offset()<< " dims="<<grid->dimension() << std::endl;
          }
          if( ! is_inside(grid->cell_bounds(loc),r) )
          {
            fatal_error()<<"particle #"<<at_id<<", ro="<<r<<", r="<<r<<"<< in cell "<<loc<<" not inside cell bounds ="<<grid->cell_bounds(loc)<<std::endl;
          }
        }
        ldbg << "dom_bbox=" << dom_bbox <<std::endl;

        for(const auto& tp : atom_data)
        {
          Vec3d r = { tp[field::rx] , tp[field::ry] , tp[field::rz] };
          IJK loc = grid->locate_cell( r ); // domain_periodic_location( *domain, r );
          grid->cell(loc).push_back( tp );
        }

        //fatal_error() << "not fully implemented" << std::endl;
        
        molecules->m_molecules.clear();
        for(const auto& p : molecule_name_map)
        {
          MoleculeSpecies mol = {};
          mol.set_name( p.first );
          molecules->m_molecules.push_back( mol );
        }
      }

      // broadcast molecule specuies names
      int nmol = molecules->m_molecules.size();
      MPI_Bcast(&nmol, 1, MPI_INT, 0, *mpi);
      molecules->m_bridge_molecules.clear();
      molecules->m_molecules.resize( nmol );
      MPI_Bcast( molecules->m_molecules.data() , sizeof(MoleculeSpecies)*nmol , MPI_CHARACTER , 0 , *mpi );

      // Broadcast domain parameters computed on MPI rank 0
      MPI_Bcast( domain.get_pointer() , sizeof(Domain), MPI_CHARACTER, 0, *mpi);
      assert( check_domain(*domain) );

      grid->set_offset( IJK{0,0,0} );
      grid->set_origin( domain->bounds().bmin );
      grid->set_cell_size( domain->cell_size() );
      grid->set_dimension( domain->grid_dimension() );
      grid->rebuild_particle_offsets();

      lout << "Particles        = "<<count<<std::endl;
      lout << "Molecules types  = "<<nmol<<std::endl;
      lout << "Domain XForm     = "<<domain->xform()<<std::endl;
      lout << "Domain bounds    = "<<domain->bounds()<<std::endl;
      lout << "Domain size      = "<<bounds_size(domain->bounds()) <<std::endl;
      lout << "Real size        = "<<bounds_size(domain->bounds()) * Vec3d{domain->xform().m11,domain->xform().m22,domain->xform().m33} <<std::endl;
      lout << "Cell size        = "<<domain->cell_size()<<std::endl;
      lout << "Grid dimensions  = "<<domain->grid_dimension()<<" ("<<grid_cell_count(domain->grid_dimension())<<" cells)"<< std::endl;
      lout << "================================="<< std::endl;
#     ifndef NDEBUG
      bool particles_inside_cell = check_particles_inside_cell(*grid);
      assert( particles_inside_cell );
#     endif
    }

  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(read_xyz_molecules)
  {
    OperatorNodeFactory::instance()->register_factory("read_xyz_molecules", make_grid_variant_operator< ReadXYZMolecules >);
  }

}
