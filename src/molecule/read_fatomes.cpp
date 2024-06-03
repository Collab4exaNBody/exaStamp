

#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/fields.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/log.h>
#include <exanb/core/unityConverterHelper.h>

#include <exanb/io/mpi_file_io.h>
#include <exaStamp/molecule/stampv4_io.h>
#include <exanb/mpi/all_value_equal.h>
#include <onika/oarray.h>

#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>

//#include "exanb/vector_utils.h"
#include <exanb/core/file_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/molecule/molecule_species.h>
#include <exaStamp/molecule/bonds_potentials_parameters.h>
#include <exaStamp/molecule/bends_potentials_parameters.h>
#include <exaStamp/molecule/torsions_potentials_parameters.h>
#include <exaStamp/molecule/impropers_potentials_parameters.h>

#include <exanb/core/check_particles_inside_cell.h>
#include <exanb/core/parallel_random.h>
#include <exanb/core/math_utils.h>
#include <exanb/core/basic_types_operators.h>
#include <exanb/defbox/deformation.h>
#include <exanb/defbox/deformation_math.h>

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
  class ReadFAtomesMolecules : public OperatorNode
  {
    ADD_SLOT( MPI_Comm        , mpi          , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( std::string     , filename     , INPUT , REQUIRED );
    ADD_SLOT( Domain          , domain       , INPUT_OUTPUT );
    ADD_SLOT( GridT           , grid         , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies , species      , INPUT ); // optional. if no species given, type ids are allocated automatically

    ADD_SLOT( MoleculeSpeciesVector , molecules , INPUT_OUTPUT, MoleculeSpeciesVector{} , DocString{"Molecule descriptions"} );

    ADD_SLOT( BondsPotentialParameters     , potentials_for_bonds     , INPUT_OUTPUT, BondsPotentialParameters{} );
    ADD_SLOT( BendsPotentialParameters     , potentials_for_angles    , INPUT_OUTPUT, BendsPotentialParameters{} );
    ADD_SLOT( TorsionsPotentialParameters  , potentials_for_torsions  , INPUT_OUTPUT, TorsionsPotentialParameters{} );
    ADD_SLOT( ImpropersPotentialParameters , potentials_for_impropers , INPUT_OUTPUT, ImpropersPotentialParameters{} );

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
      int nb_atom_types = 0;
      std::string atom_name;
      double atom_mass = std::numeric_limits<double>::quiet_NaN();
      double atom_charge = std::numeric_limits<double>::quiet_NaN();
      double Lx=0.0, Ly=0.0, Lz=0.0;
      double angle_a_b=90.0, angle_a_c=90.0, angle_b_c=90.0;
      std::string atom_speed_unit;
      std::string atom_charge_unit;
      std::string atom_pos_unit;
      bool force_field_description = false;
      bool atom_bonds_description = false;
      
      if(rank==0)
      {
        //Open xyz file
        std::ifstream file(file_name);
        if( ! file )
        {
          fatal_error() << "FAtomes file "<< file_name << " not found !" << std::endl;
        }

        int64_t lineno = -1;
        auto next_line = [&file,&lineno]() -> std::istringstream
        {
          std::string line;
          do
          {
            char buffer[1024];
            file.getline(buffer,1023,'\n');
            ++ lineno;
            line = buffer;
            auto com = line.find('*');
            if( com != std::string::npos ) line = line.substr(0,com);
            while( ! line.empty() && ( line.back()==' ' ||  line.back()=='\t' || line.back()=='\n') ) line.pop_back();
          } while( line.empty() && ! file.eof() );
          return std::istringstream(line);
        };

        auto read_unit = [this] ( std::string u ) -> std::string
        {
          const std::map<std::string,std::string> replace = { {"ang2","ang^2"} , {"deg","degree"} };
          for(const auto& r : replace)
          {
            auto p = u.find(r.first);
            while( p != std::string::npos )
            {
              //ldbg << "found "<<r.first<<" @"<<p<<std::endl;
              if( r.second.find(r.first)!=0 || p != u.find(r.second,p) ) u.replace( p , r.second.length() , r.second );
              p = u.find( r.first , p + r.second.length() );
            }
          }
          return u;
        };

        auto process_kw_line = [&] ( std::istringstream && iss ) 
        {
          std::string kw;
          std::string tmp;
          std::string unit;
          iss >> kw;
          if( kw == "NbTypesAtomes" ) { iss >> nb_atom_types; ldbg << "NbTypesAtomes = "<<nb_atom_types<<std::endl; }
          else if( kw == "structure" ) { iss >> tmp; ldbg << "structure = "<<tmp<<std::endl; }
          else if( kw == "maille_long" ) { iss >> Lx >> Ly >> Lz; ldbg << "maille_long = "<<Lx<<" , " <<Ly << " , "<<Lz <<std::endl; }
          else if( kw == "maille_angle" ) { iss >> angle_a_b >> angle_a_c >> angle_b_c; ldbg << "maille_angle = "<<angle_a_b<<" , " <<angle_a_c<< " , "<<angle_b_c <<std::endl; }
          else if( kw == "maille_orient" ) { iss >> tmp; ldbg << "maille_orient = "<<tmp<<std::endl; }
          else if( kw == "maille_ref" ) { iss >> tmp; ldbg << "maille_ref = "<<tmp<<std::endl; }
          else if( kw == "nom" ) { iss >> atom_name; ldbg << "nom = "<<atom_name<<std::endl; }
          else if( kw == "nomXYZ" ) { iss >> tmp; ldbg << "nomFF = "<<tmp<<std::endl; }
          else if( kw == "nomFF" ) { iss >> tmp; ldbg << "nomFF = "<<tmp<<std::endl; }
          else if( kw == "type" ) { iss >> tmp; ldbg << "type = "<<tmp<<std::endl; }
          else if( kw == "masse" )
          {
            iss >> atom_mass >> unit;
            unit = read_unit(unit);
            atom_mass = exanb::make_quantity( atom_mass , unit ).convert();
            ldbg << "masse = " << atom_mass << " (unit="<<unit<<")"<< std::endl;
          }
          else if( kw == "charge" )
          {
            iss >> atom_charge >> unit;
            unit = read_unit(unit);
            atom_charge = exanb::make_quantity( atom_charge , unit ).convert();
            ldbg << "charge = " << atom_charge << " (unit="<<unit<<")"<< std::endl;
          }
          else if( kw == "Potentiel" )
          {
            int ta=0, tb=0;
            iss >> ta >> tb >> kw;
            ldbg << "Skip pair potential "<<kw<<" for pair "<<ta<<" , "<<tb<<std::endl;
          }
          else if( kw == "Regle_melange" ) { iss >> tmp; ldbg << "Regle_melange = "<<tmp<<std::endl; }
          else if( kw == "PositionDesAtomesCart" ) { iss >> atom_pos_unit; ldbg << "PositionDesAtomesCart = "<<atom_pos_unit<<std::endl; }
          else if( kw == "VitesseDesAtomes" ) { iss >> atom_speed_unit; ldbg << "VitesseDesAtomes = "<<atom_speed_unit<<std::endl; }
          else if( kw == "ModificationChargeDesAtomes" ) { iss >> atom_charge_unit; ldbg << "ModificationChargeDesAtomes = "<<atom_charge_unit<<std::endl; }
          else if( kw == "ChampDeForces" ) { force_field_description = true; }
          else if( kw == "Zmatrice" ) { atom_bonds_description = true; }
          else
          {
            fatal_error() << "unexpected token "<<kw << " at line "<<lineno <<std::endl;
          }
        };        

        while( nb_atom_types == 0 )
        {
          process_kw_line( next_line() );
        }
        ldbg << "nb_atom_types="<<nb_atom_types<<std::endl;
        species->clear();
        species->resize(nb_atom_types);

        int atom_type_count = 0;
        while( atom_type_count < nb_atom_types )
        {
          //ldbg << "read atom types : atom_type_count="<<atom_type_count<<std::endl;
          process_kw_line( next_line() );
          if( !atom_name.empty() && !std::isnan(atom_mass) && !std::isnan(atom_charge) )
          {
            ldbg << "Add atom type "<<atom_name<<" : mass="<<atom_mass<<" , charge="<<atom_charge<<std::endl;
            species->at(atom_type_count).set_name( atom_name );
            species->at(atom_type_count).m_charge = atom_charge;
            species->at(atom_type_count).m_mass = atom_mass;
            species->at(atom_type_count).m_rigid_atoms[ 0 ] = { Vec3d{0.,0.,0.} , -1 };
            species->at(atom_type_count).set_rigid_atom_name( 0 , atom_name );
            species->at(atom_type_count).m_rigid_atom_count = 1;
            atom_name.clear();
            atom_mass = std::numeric_limits<double>::quiet_NaN();
            atom_charge = std::numeric_limits<double>::quiet_NaN();
            ++ atom_type_count;
          }
        }

        while( atom_pos_unit.empty() )
        {
          process_kw_line( next_line() );
        }
        
        next_line() >> count;
        ldbg << "Number of atoms = "<<count<<std::endl;        

        Mat3d H = make_identity_matrix();        
        if( angle_a_b!=90 && angle_a_c!=90 && angle_b_c!=90 )
        {
          H = deformation_to_matrix( Deformation{ { angle_a_b*M_PI/180, angle_a_c*M_PI/180, angle_b_c*M_PI/180 } , { Lx , Ly , Lz } } );
        }
        else
        {
          H = make_diagonal_matrix( { Lx , Ly , Lz } );
        }
                
        ldbg << "H = " << H << std::endl;
        double smin=1.0, smax=1.0;
        matrix_scale_min_max(H,smin,smax);
        ldbg << "scale = "<<smin<<" / "<<smax<<std::endl;

        if( is_diagonal(H) && Lx==Ly && Ly==Lz )
        {
          ldbg << "cubic & diagonal H matrix" << std::endl;
        }
        
        const Mat3d inv_H = inverse(H);
        ldbg << "H^-1 = " << inv_H << std::endl;

        // read one line at a time
        std::string typeMol;

        AABB file_bbox = { {0.0,0.0,0.0} , {0.0,0.0,0.0} };

        std::map< std::string , int > molecule_name_map;
        std::vector<MoleculeTupleIO> atom_data;
        atom_data.reserve( count );
        std::string typeAtom;
      
        for(uint64_t cnt=0;cnt<count;cnt++)
        {
          //int64_t moleculeid = -1;
          double x=0.0, y=0.0, z=0.0;
          //int64_t c0=-1 , c1=-1 , c2=-1 , c3=-1 , c4=-1;
          double charge = 0.0;

          assert( ! file.eof() && file.good() );
          next_line() >> typeAtom >> x >> y >> z;
          int64_t at_id = cnt;

          // convert position to grid's base
          Vec3d r{x, y, z};
          if( cnt == 0 )
          {
            file_bbox.bmin = file_bbox.bmax = r;
          }
          else
          {
            file_bbox.bmin = min( file_bbox.bmin , r );
            file_bbox.bmax = max( file_bbox.bmax , r );
          }
          r = inv_H * r; // reduced coordinates in [0;1]

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
          if( itypeAtom < 0 || size_t(itypeAtom) >= species->size() )
          {
            fatal_error() << "bad atom type index "<<itypeAtom<<std::endl;
          }
          charge = species->at(itypeAtom).m_charge;

          MoleculeTupleIO tp( r.x, r.y, r.z , 0.0, 0.0, 0.0, at_id, itypeAtom, 0, std::array<uint64_t,4>{uint64_t(-1),uint64_t(-1),uint64_t(-1),uint64_t(-1)} , charge );
          atom_data.push_back( tp );     
        }
        ldbg << "file bbox = " << file_bbox << " , size = " << ( file_bbox.bmax - file_bbox.bmin ) << std::endl;       
       
        // read atom speeds
        if( atom_speed_unit.empty() ) process_kw_line( next_line() );
        if( atom_speed_unit.empty() )
        {
          fatal_error() << "Keyword VitesseDesAtomes not found as expected, at line "<<lineno<<std::endl;
        }
        atom_speed_unit = read_unit(atom_speed_unit);
        for(uint64_t cnt=0;cnt<count;cnt++)
        {
          double vx=0.0, vy=0.0, vz=0.0;
          next_line() >> vx >> vy >> vz;
          atom_data[cnt][field::vx] = exanb::make_quantity( vx , atom_speed_unit ).convert();
          atom_data[cnt][field::vy] = exanb::make_quantity( vy , atom_speed_unit ).convert();
          atom_data[cnt][field::vz] = exanb::make_quantity( vz , atom_speed_unit ).convert();
        }
        
        process_kw_line( next_line() );
        if( ! atom_charge_unit.empty() )
        {
          atom_charge_unit = read_unit(atom_charge_unit); 
          uint64_t nb_charge_modif = 0;
          next_line() >> nb_charge_modif;
          ldbg << "Modifying individual charges of "<<nb_charge_modif<<" atoms"<<std::endl;
          for(uint64_t cnt=0;cnt<nb_charge_modif;cnt++)
          {
            int64_t at_id = -1;
            double charge = 0.0;
            next_line() >> at_id >> charge;
            if( at_id < 0 || size_t(at_id) >= count )
            {
              fatal_error() << "invalid atom id "<<at_id<<" found at line "<<lineno<<std::endl;
            }
            atom_data[at_id][field::charge] = exanb::make_quantity( charge , atom_charge_unit ).convert();
          }
          process_kw_line( next_line() );
        }
        
        if( force_field_description )
        {
          uint64_t nb_force_fields = 0;
          next_line() >> nb_force_fields;
          ldbg << "Reading "<<nb_force_fields<<" force fields"<<std::endl;
          for(uint64_t cnt=0;cnt<nb_force_fields;cnt++)
          {
            auto iss = next_line();
            std::string fftype;
            iss >> fftype;
            if( fftype == "bond_opls" )
            {
              double p1, p2, p3, p4;
              std::string t1, t2;
              std::string u1, u2, u3, u4;
              iss >> t1 >> t2 >> p1 >> u1 >> p2 >> u2 >> p3 >> u3 >> p4 >> u4;
              double r0 = exanb::make_quantity( p1 , read_unit(u1) ).convert(); if( std::isnan(r0) ) r0=0.0;
              double k  = exanb::make_quantity( p2 , read_unit(u2) ).convert(); if( std::isnan(k)  ) k=0.0;
              // p3 = exanb::make_quantity( p3 , read_unit(u3) ).convert(); if( std::isnan(p3) ) p3=0.0; // 2nd and 3rd parameter ignored
              // p4 = exanb::make_quantity( p4 , read_unit(u4) ).convert(); if( std::isnan(p4) ) p4=0.0;
              ldbg << "Intramolecular "<<fftype<<" for "<<t1<<","<<t2<<" : r0="<<r0<<" , k="<<k << std::endl;
              potentials_for_bonds->m_bond_desc.push_back( BondPotential{ fftype , {t1,t2} , std::make_shared<IntraMolecularBondOPLSFunctional>(k,r0) } );
            }
            else if( fftype == "angle_opls" )
            {
              double p1, p2, p3, p4;
              std::string t1, t2, t3;
              std::string u1, u2, u3, u4;
              iss >> t1 >> t2 >> t3 >> p1 >> u1 >> p2 >> u2 >> p3 >> u3 >> p4 >> u4;
              double phi0 = exanb::make_quantity( p1 , read_unit(u1) ).convert(); if( std::isnan(phi0) ) phi0=0.0;
              double k    = exanb::make_quantity( p2 , read_unit(u2) ).convert(); if( std::isnan(k)    ) k=0.0;
              // p3 = exanb::make_quantity( p3 , read_unit(u3) ).convert(); if( std::isnan(p3) ) p3=0.0;  // 2nd and 3rd parameter ignored
              // p4 = exanb::make_quantity( p4 , read_unit(u4) ).convert(); if( std::isnan(p4) ) p4=0.0;
              ldbg << "Intramolecular "<<fftype<<" for "<<t1<<","<<t2<<" : "<<t3<<" : phi0="<<phi0<<" , k="<<k << std::endl;
              potentials_for_angles->m_potentials.push_back( BendPotential{ fftype , {t1,t2,t3} , std::make_shared<IntraMolecularBondOPLSFunctional>(k,phi0) } );
            }
            else if( fftype == "torsion_opls" )
            {
              double p1, p2, p3;
              std::string t1, t2, t3, t4;
              std::string u1, u2, u3;
              iss >> t1 >> t2 >> t3 >> t4 >> p1 >> u1 >> p2 >> u2 >> p3 >> u3 ;
              p1 = exanb::make_quantity( p1 , read_unit(u1) ).convert(); if( std::isnan(p1) ) p1=0.0;
              p2 = exanb::make_quantity( p2 , read_unit(u2) ).convert(); if( std::isnan(p2) ) p2=0.0;
              p3 = exanb::make_quantity( p3 , read_unit(u3) ).convert(); if( std::isnan(p3) ) p3=0.0;
              ldbg << "Intramolecular force fireld "<<fftype<<" for "<<t1<<","<<t2<<","<<t3<<","<<t4<<" : p1="<<p1<<" , p2="<<p2<<" , p3="<<p3<< std::endl;
              potentials_for_torsions->m_potentials.push_back( TorsionPotential{ fftype , {t1,t2,t3,t4} , std::make_shared<IntraMolecularCosOPLSFunctional>(p1,p2,p3) } );
            }
          }
          process_kw_line( next_line() );
        }
        
        if( atom_bonds_description )
        {
          uint64_t nb_con = 0;
          next_line() >> nb_con;
          ldbg << "Reading connectivity for "<<nb_con<<" atoms"<<std::endl;
          for(uint64_t cnt=0;cnt<nb_con;cnt++)
          {
            int64_t at_id=-1, c1=-1, c2=-1, c3=-1, c4=-1;
            next_line() >> at_id >> c1 >> c2 >> c3 >> c4;
            if( at_id < 0 || size_t(at_id) >= count )
            {
              fatal_error() << "invalid atom id "<<at_id<<" found at line "<<lineno<<std::endl;
            }
            atom_data[at_id][field::cmol] = { uint64_t(c1), uint64_t(c2), uint64_t(c3), uint64_t(c4) };
          }
        }
        
        file.close();

        
        double cell_size = domain->cell_size() > 0.0 ? domain->cell_size() : 8.0 ;
        ldbg << "cell_size = "<<cell_size << std::endl;
        const IJK grid_dims = { static_cast<ssize_t>(ceil(Lx/cell_size)) , static_cast<ssize_t>(ceil(Ly/cell_size)) , static_cast<ssize_t>(ceil(Lz/cell_size)) };
        if( grid_dims.i==grid_dims.j && grid_dims.j==grid_dims.k )
        {
          double cell_size_ratio = ( Lx / grid_dims.i ) / cell_size;
          if( cell_size_ratio > 0.5 && cell_size_ratio < 1.5 )
          {
            cell_size = Lx / grid_dims.i;
            ldbg << "modify cell_size to "<<cell_size<<" to avoid XForm"<<std::endl;
          }
        }
        domain->set_cell_size( cell_size );

        domain->set_grid_dimension( grid_dims );
        ldbg << "domain grid = "<<grid_dims<< std::endl;
        Vec3d dom_size = grid_dims * cell_size;
        ldbg << "domain size = "<< dom_size << std::endl;
        domain->set_bounds( { {0.0,0.0,0.0} , dom_size } );
        const Mat3d dom_H = make_diagonal_matrix( dom_size );
        domain->set_xform( H * inverse( dom_H ) );
        ldbg << "domain = " << *domain << std::endl;
        const Mat3d dom_xform = domain->xform() ;

        const Vec3d xform_a = { dom_xform.m11, dom_xform.m12, dom_xform.m13 };
        const Vec3d xform_b = { dom_xform.m21, dom_xform.m22, dom_xform.m23 };
        const Vec3d xform_c = { dom_xform.m31, dom_xform.m32, dom_xform.m33 };
        const double xform_Lx = norm(xform_a);
        const double xform_Ly = norm(xform_b);
        const double xform_Lz = norm(xform_c);
        ldbg << "xform_a = " << xform_a << " , norm = " << xform_Lx << std::endl;        
        ldbg << "xform_b = " << xform_b << " , norm = " << xform_Ly << std::endl;        
        ldbg << "xform_x = " << xform_c << " , norm = " << xform_Lz << std::endl;        
        ldbg << "xform angle a.b = " << (180/M_PI) * acos( dot(xform_a,xform_b) / ( xform_Lx*xform_Ly ) ) << std::endl;        
        ldbg << "xform angle a.c = " << (180/M_PI) * acos( dot(xform_a,xform_c) / ( xform_Lx*xform_Lz ) ) << std::endl;        
        ldbg << "xform angle b.c = " << (180/M_PI) * acos( dot(xform_b,xform_c) / ( xform_Ly*xform_Lz ) ) << std::endl;        
                                    
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
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory("read_fatomes_mol", make_grid_variant_operator< ReadFAtomesMolecules >);
  }

}
