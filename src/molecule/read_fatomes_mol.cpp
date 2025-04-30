#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/grid_fields.h>
#include <onika/math/basic_types_stream.h>
#include <onika/log.h>
#include <onika/physics/units.h>

#include <exanb/io/mpi_file_io.h>
#include <exaStamp/molecule/stampv4_io.h>
#include <onika/mpi/all_value_equal.h>
#include <onika/oarray.h>

#include <onika/math/basic_types_yaml.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>

//#include "exanb/vector_utils.h"
#include <onika/file_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/molecule/molecule_species.h>
#include <exaStamp/molecule/bonds_potentials_parameters.h>
#include <exaStamp/molecule/bends_potentials_parameters.h>
#include <exaStamp/molecule/torsions_potentials_parameters.h>
#include <exaStamp/molecule/impropers_potentials_parameters.h>
#include <exaStamp/molecule/intramolecular_pair_weight.h>
#include <exaStamp/molecule/molecule_optional_header_io.h>
#include <exaStamp/molecule/iFFname.h>

#include <exaStamp/potential/ljexp6rf/ljexp6rf.h>

#include <exanb/core/check_particles_inside_cell.h>
#include <onika/parallel/random.h>
#include <onika/math/math_utils.h>
#include <onika/math/basic_types_operators.h>
#include <onika/string_utils.h>
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

  // Read FAtomes files.
  // This files must be ASCII files
  // ATTENTION : we make many assumptions (not made in Stamp but valid for our test case so far) in this  function.
  // We assume that the order of the sections is
  // 1. atom types (also true in Stamp)
  // 2. "Potentiel"
  // 3. PositionDesAtomes
  // 4. Vitesses
  // 5. ModificationDesCharges
  // 6. ChampDeForce (intra)
  // 7. Zmatrice
  // 8. ContribDispRepIntra

  template<typename GridT>
  class ReadFAtomesMolecules : public OperatorNode
  {
    ADD_SLOT( MPI_Comm        , mpi          , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( std::string     , filename     , INPUT , REQUIRED );
    ADD_SLOT( std::string     , stampinput   , INPUT , OPTIONAL );
    ADD_SLOT( Domain          , domain       , INPUT_OUTPUT );
    ADD_SLOT( GridT           , grid         , INPUT_OUTPUT );

    ADD_SLOT( ParticleSpecies , species      , INPUT_OUTPUT ); 
    ADD_SLOT( MoleculeSpeciesVector , molecules , INPUT_OUTPUT, MoleculeSpeciesVector{} , DocString{"Molecule descriptions"} );

    ADD_SLOT( BondsPotentialParameters         , potentials_for_bonds     , INPUT_OUTPUT, BondsPotentialParameters{} );
    ADD_SLOT( BendsPotentialParameters         , potentials_for_angles    , INPUT_OUTPUT, BendsPotentialParameters{} );
    ADD_SLOT( TorsionsPotentialParameters      , potentials_for_torsions  , INPUT_OUTPUT, TorsionsPotentialParameters{} );
    ADD_SLOT( ImpropersPotentialParameters     , potentials_for_impropers , INPUT_OUTPUT, ImpropersPotentialParameters{} );    
    ADD_SLOT( LJExp6RFMultiParms               , potentials_for_pairs     , INPUT_OUTPUT , LJExp6RFMultiParms{} );
    ADD_SLOT( IntramolecularPairWeighting      , mol_pair_weights         , INPUT_OUTPUT , IntramolecularPairWeighting{} );
    ADD_SLOT( IntraFFtypes                     , ffnameVector             , INPUT_OUTPUT);

    ADD_SLOT( double                           , rcut_max                 , INPUT_OUTPUT , 0.0 );

  public:
    inline void execute () override final
    {
      //-------------------------------------------------------------------------------------------
      std::string file_name = onika::data_file_path( *filename );
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
      std::string atom_FFname;
      double atom_mass = std::numeric_limits<double>::quiet_NaN();
      double atom_charge = std::numeric_limits<double>::quiet_NaN();
      double Lx=0.0, Ly=0.0, Lz=0.0;
      double angle_a_b=90.0, angle_a_c=90.0, angle_b_c=90.0;
      std::string atom_speed_unit;
      std::string atom_charge_unit;
      std::string atom_pos_unit;
      std::string mixingRule = "0";
      bool force_field_description = false;
      bool atom_bonds_description = false;
      bool intramol_weight_description = false;
      
      if(rank==0)
      {
        //Open xyz file
        std::ifstream file(file_name);
        if( ! file )
        {
          fatal_error() << "FAtomes file "<< file_name << " not found !" << std::endl;
        }

        // lambda function to get the next line of FAtomes file
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

        // lambda function to read units
        auto read_unit = [this] ( std::string u ) -> std::string
        {
          //const std::map<std::string,std::string> replace = { {"J/mol.m6","J*m^6   "} , {"ang2","ang^2"} , {"deg","degree"} , {"m6","m^6"} , {"m-1","m^-1"} , {"kcal/mol", "kcal"}};
          const std::map<std::string,std::string> replace = { {"J.m6","J*m^6   "} , {"ang2","ang^2"} , {"deg","degree"} , {"m6","m^6"} , {"m-1","m^-1"} , {"metre", "m"} };

          for(const auto& r : replace)
          {
            auto p = u.find(r.first);
            while( p != std::string::npos )
            {
              //ldbg << "found "<<r.first<<" @"<<p<<std::endl;
              //if( r.second.find(r.first)!=0 || p != u.find(r.second,p) ) u.replace( p , r.second.length() , r.second );
	      u.replace( p , r.first.length() , r.second );
              p = u.find( r.first , p + r.second.length() );
            }
          }
          return u;
        };

        // lambda function to read all keywords of FAtomes file, save right away what can be saved 
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
          else if( kw == "nomFF" ) { iss >> atom_FFname; ldbg << "nomFF = "<< atom_FFname <<std::endl; }
          else if( kw == "type" ) { iss >> tmp; ldbg << "type = "<<tmp<<std::endl; }
          else if( kw == "masse" )
          {
            iss >> atom_mass >> unit;
            unit = read_unit(unit);
            atom_mass = onika::physics::make_quantity( atom_mass , unit ).convert();
            ldbg << "masse = " << atom_mass << " (unit="<<unit<<")"<< std::endl;
          }
          else if( kw == "charge" )
          {
            iss >> atom_charge >> unit;
            unit = read_unit(unit);
            atom_charge = onika::physics::make_quantity( atom_charge , unit ).convert();
            ldbg << "charge = " << atom_charge << " (unit="<<unit<<")"<< std::endl;
          }
          else if( kw == "Potentiel" )
          {
            int ita=0, itb=0;
	    std::string ta, tb;
            std::string typepot;
            iss >> ta >> tb >> typepot;

            // if ta and tb are atom type ids,
            if (isdigit(ta[0])) {
              // we convert the string read as an integer (the atom type id)
              ita = stoi(ta);
	      // and the type id into the correct string corresponding to the type
	      ta = species->at(ita).name();
	    }
            if (isdigit(tb[0])) {
              itb = stoi(tb);
	      tb = species->at(itb).name();
            }

            //std::string ta = species->at(ita).name();
            //std::string tb = species->at(itb).name();
            ldbg << "pair potential "<<typepot<<" for pair "<<ta<<" , "<<tb<<std::endl;
            std::map<std::string,double> params;
            LJExp6RFMultiParmsPair pot = { ta, tb }; // declare 'mixed' potential
            if( typepot == "LJ" )
            {
              std::string p1,u1, p2, u2, p3, u3;
              double v1, v2, v3;
              iss >> p1 >> v1 >> u1 >> p2 >> v2 >> u2 >> p3 >> v3 >> u3;
              params[p1] = onika::physics::make_quantity( v1 , read_unit(u1) ).convert();
              params[p2] = onika::physics::make_quantity( v2 , read_unit(u2) ).convert();
              params[p3] = onika::physics::make_quantity( v3 , read_unit(u3) ).convert();
              double sigma = params["sigma"];
              double epsilon = params["epsilon"];
              double rc = params["rc"];
              ldbg<< "LJ params : sigma="<<sigma<<" , epsilon="<<epsilon<<" , rc="<<rc<<std::endl;
              std::string ext;
              iss >> ext;
              if( ext == "CONV_EXP6v1" )
              {
                double alpha=0.0, D=0.0;
                iss >> p1 >> alpha >> p2 >> D >> u2;
                D = onika::physics::make_quantity( D , read_unit(u2) ).convert();
                double rmin = sigma * pow( 2. , 1./6. );
                double A = 6. * epsilon * exp(alpha) / ( alpha - 6. );
                double B = alpha / rmin;
                double C = epsilon * alpha / ( alpha - 6. ) * pow(rmin,6);
                ldbg << "Exp6 conversion : alpha="<<alpha<<" , D="<<D<<" , rmin="<<rmin<<" , A="<<A<<" , B="<<B<<" , C="<<C << std::endl;
                pot.m_params.set_exp6_parameters( A, B, C, D, rc ); // set exp6 parameters
              }
              else
              {
                if( !ext.empty() ) { fatal_error()<<"unexpected LJ modifier '"<<ext<<"'"<<std::endl; }
                pot.m_params.set_lj_parameters( sigma, epsilon, rc ); // set LJ parameters
              }
            }
            else if( typepot == "Exp6v1" )
            {
              std::string p1,u1, p2, u2, p3, u3, p4, u4, p5, u5;
              double v1, v2, v3, v4, v5;
              iss >> p1 >> v1 >> u1 >> p2 >> v2 >> u2 >> p3 >> v3 >> u3 >> p4 >> v4 >> u4 >> p5 >> v5 >> u5;
              u1 = read_unit(u1);
              u2 = read_unit(u2);
              u3 = read_unit(u3);
              u4 = read_unit(u4);
              u5 = read_unit(u5);
              ldbg << p1<<"=" << v1 <<" "<< u1<<" , "<< p2<<"=" << v2<<" " << u2<<" , " << p3<<"=" << v3<<" " << u3<<" , " << p4<<"=" << v4<<" " << u4<<" , " << p5<<"=" << v5<<" " << u5<<std::endl;
              params[p1] = onika::physics::make_quantity( v1 , u1 ).convert();
              params[p2] = onika::physics::make_quantity( v2 , u2 ).convert();
              params[p3] = onika::physics::make_quantity( v3 , u3 ).convert();
              params[p4] = onika::physics::make_quantity( v4 , u4 ).convert();
              params[p5] = onika::physics::make_quantity( v5 , u5 ).convert();
              const double A=params["a"], B=params["b"], C=params["c"], D=params["d"], rc=params["rc"];
              ldbg << "Exp6v1 : A="<<A<<" , B="<<B<<" , C="<<C<<" , D="<<D<<" , rc="<<rc<<std::endl;
              pot.m_params.set_exp6_parameters( A, B, C, D, rc ); // set exp6 parameters
            }
            else
            {
              fatal_error() << "Potential "<<typepot<<" is not supported"<<std::endl;
            }
            pot.m_params.update_ecut();
            potentials_for_pairs->m_potentials.push_back(pot); // Important: add pair potential parameters
          }
          else if( kw == "Regle_melange" )
          {
            iss >> mixingRule; ldbg << "Regle_melange = "<< mixingRule <<std::endl;
          }
          else if( kw == "PositionDesAtomesCart" ) { iss >> atom_pos_unit; ldbg << "PositionDesAtomesCart = "<<atom_pos_unit<<std::endl; }
          else if( kw == "VitesseDesAtomes" ) { iss >> atom_speed_unit; ldbg << "VitesseDesAtomes = "<<atom_speed_unit<<std::endl; }
          else if( kw == "ModificationChargeDesAtomes" ) { iss >> atom_charge_unit; ldbg << "ModificationChargeDesAtomes = "<<atom_charge_unit<<std::endl; }
          else if( kw == "ChampDeForces" ) { force_field_description = true; }
          else if( kw == "Zmatrice" ) { atom_bonds_description = true; }
          else if( kw == "ContribDispRepIntra" ) { intramol_weight_description = true; }
          else
          {
            fatal_error() << "unexpected token "<<kw << " at line "<<lineno <<std::endl;
          }
        };        

        // to make sure that the number of atom types is read
        while( nb_atom_types == 0 )
        {
          process_kw_line( next_line() );
        }
        ldbg << "nb_atom_types="<<nb_atom_types<<std::endl;
        species->clear();
        species->resize(nb_atom_types);

        // to save the info (charge, mass, name, type,...) of each atom type
        int atom_type_count = 0;
        std::map<std::string,size_t> atom_type_index;
        while( atom_type_count < nb_atom_types )
        {
          //ldbg << "read atom types : atom_type_count="<<atom_type_count<<std::endl;
          process_kw_line( next_line() );
	  /* Attention: this line assumes that either nom or masse or charge is the last useful feature read for an atom type */
	  /* if a feature is located after that, it is not saved. The next "while" goes directly to the next line without saving the piece of info */
          if( !atom_name.empty() && !std::isnan(atom_mass) && !std::isnan(atom_charge) )
          {
            ldbg << "Add atom type "<<atom_name<<" : mass="<<atom_mass<<" , charge="<<atom_charge<<std::endl;
            species->at(atom_type_count).set_name( atom_name );
	    if (atom_FFname.empty()) species->at(atom_type_count).set_FFname( atom_name );
	    else {
              species->at(atom_type_count).set_FFname( atom_FFname );

	    }
            species->at(atom_type_count).m_charge = atom_charge;
            species->at(atom_type_count).m_mass = atom_mass;
            species->at(atom_type_count).m_rigid_atoms[ 0 ] = { Vec3d{0.,0.,0.} , -1 };
            species->at(atom_type_count).set_rigid_atom_name( 0 , atom_name );
            species->at(atom_type_count).m_rigid_atom_count = 1;
            atom_type_index[atom_name] = atom_type_count;
            atom_name.clear();
            atom_FFname.clear();
            atom_mass = std::numeric_limits<double>::quiet_NaN();
            atom_charge = std::numeric_limits<double>::quiet_NaN();
            ++ atom_type_count;
          }
        }

	//section to fill in ffnameVector (such that ffnameVector[FFtypeId] = FFname)
	// The fact that it is a vector is useful to ease the future MPI_broadcast
	// Locally, we use an equivalent map, ffTypeNameToFFTypeId, encoding the same info
	std::string FFname_loc;

	// local map
	std::map< std::string , int > ffTypeNameToFFTypeId;

	// Loop on all (full) atom types
	for (atom_type_count = 0; atom_type_count < nb_atom_types; atom_type_count++) {
	  // string correspondig to FFname of the current atomic type
          FFname_loc = species->at(atom_type_count).FFname();

	  //default value for m_FFtypeId in case FFname is not provided
	  int ffTypeId = atom_type_count;

          // find or allocate an id for the FFType name string
	  auto it = ffTypeNameToFFTypeId.find( FFname_loc );
	  if( it == ffTypeNameToFFTypeId.end() )
	  {
	    ffTypeId = ffTypeNameToFFTypeId.size();
	    ffTypeNameToFFTypeId[FFname_loc] = ffTypeId; // create the entry
	  }
	  else
	  {
	    ffTypeId = ffTypeNameToFFTypeId[FFname_loc];
	  }

	  // fill in species with the corresponding ffTypeId
	  species->at(atom_type_count).m_FFtypeId = ffTypeId;

	  ldbg << "In read_fatomes_mol.cpp full atom type id " <<  atom_type_count << " (" << species->at(atom_type_count).name() << ") corresponds to FFtype name " << species->at(atom_type_count).FFname() << " (ff id = "<< ffTypeId << ")" << std::endl;
	}
	
	// initialize ffnameVector with ffTypeNameToFFTypeId.size() empty strings
	ffnameVector->clear();
	ffnameVector->resize( ffTypeNameToFFTypeId.size() ); 
	for( const auto & ffpair : ffTypeNameToFFTypeId )
	{
	  ffnameVector->at( ffpair.second ) = ffpair.first;
	}
	for(const auto & ffname : *ffnameVector)
	{
	  assert( ! ffname.empty() );
	}
	

        // read file until section "PositionDesAtomes" and fill all variables corresponding to the sections "TypeAtomes" and "Potentiel" read before
        while( atom_pos_unit.empty() )
        {
          process_kw_line( next_line() );
        }

        // Take into account the mixing rule
        int ita, itb;
        // Loop on all "cross" pairs
        for (ita = 0; ita < atom_type_count - 1; ita++) {
          for (itb = ita+1 ; itb < atom_type_count; itb++) {
             LJExp6RFParms* params_for_pair = potentials_for_pairs->params_for_pair( species->at(ita).name(), species->at(itb).name() );
             // if this pair is not yet defined
             if( params_for_pair == nullptr )
             {
               // initialize the potential for the cross term
               LJExp6RFMultiParmsPair pot = { species->at(ita).name(), species->at(itb).name() };
               //LJExp6RFMultiParmsPair pot = { species->at(itb).name(), species->at(ita).name() };

               // get the parameters of the interactions (ita, ita) and (itb, itb)
               LJExp6RFParms* params_ita = potentials_for_pairs->params_for_pair( species->at(ita).name(), species->at(ita).name() );
               LJExp6RFParms* params_itb = potentials_for_pairs->params_for_pair( species->at(itb).name(), species->at(itb).name() );

               // check that ita and itb have the same potential type (LJ ou exp6v1)
               if (params_ita->m_B_ISEXP6 == 0. && params_itb->m_B_ISEXP6 == 0. ) {
                 // LJ potential

                 double epsilon_a, epsilon_b, epsilon_ab;
                 double sigma_a, sigma_b, sigma_ab;
                 double rcut_a, rcut_b, rcut_ab;

                 epsilon_a = params_ita->m_C_EPSILON;
                 sigma_a   = params_ita->m_D_SIGMA;
                 rcut_a    = params_ita->m_rcut;

                 epsilon_b = params_itb->m_C_EPSILON;
                 sigma_b   = params_itb->m_D_SIGMA;
                 rcut_b    = params_itb->m_rcut;

                 epsilon_ab = 0.; sigma_ab = 0.; rcut_ab = 0.;

                 // definition of a default mixing rule if undefined
                 if (mixingRule == "0") mixingRule = "Lorentz-Berthelot";

                 if (mixingRule == "geometrique") {
                   epsilon_ab = sqrt(epsilon_a*epsilon_b);
                   sigma_ab   = sqrt(sigma_a*sigma_b);
                   rcut_ab    = sqrt(rcut_a*rcut_b);

                 } else if (mixingRule == "Lorentz-Berthelot" || mixingRule == "LB") {
                   epsilon_ab = sqrt(epsilon_a*epsilon_b);
                   sigma_ab   = 0.5*(sigma_a+sigma_b);
                   rcut_ab    = 0.5*(rcut_a*rcut_b);

                 } else {
                   fatal_error() << "Unexpected mixing rule (" << mixingRule << ") for LJ potential for cross term (" << species->at(ita).name() << ", " << species->at(itb).name() << ")" << std::endl;
                 }

                 pot.m_params.set_lj_parameters( sigma_ab, epsilon_ab, rcut_ab );

               } else if (params_ita->m_B_ISEXP6 != 0. && params_ita->m_B_ISEXP6 != 0.) {
                 // exp6v1 potential
                 double A_a, A_b, A_ab;
                 double B_a, B_b, B_ab;
                 double C_a, C_b, C_ab;
                 double D_a, D_b, D_ab;
                 double rcut_a, rcut_b, rcut_ab;

                 A_a    = params_ita->m_A;
                 B_a    = params_ita->m_B_ISEXP6;
                 C_a    = params_ita->m_C_EPSILON;
                 D_a    = params_ita->m_D_SIGMA;
                 rcut_a = params_ita->m_rcut;

                 A_b    = params_itb->m_A;
                 B_b    = params_itb->m_B_ISEXP6;
                 C_b    = params_itb->m_C_EPSILON;
                 D_b    = params_itb->m_D_SIGMA;
                 rcut_b = params_itb->m_rcut;

                 A_ab = 0.; B_ab = 0.; C_ab = 0.; D_ab = 0.; rcut_ab= 0.;

                 // definition of a default mixing rule if undefined
                 if (mixingRule == "0") mixingRule = "Waldman-Hagler";

                 if (mixingRule == "geometrique") {
                   A_ab = sqrt(A_a*A_b);
                   B_ab = sqrt(B_a*B_b);
                   C_ab = sqrt(C_a*C_b);
                   D_ab = sqrt(D_a*D_b);
                   rcut_ab    = sqrt(rcut_a*rcut_b);

                 } else if (mixingRule == "Waldman-Hagler") {

                   // B parameter
                   B_ab = pow(B_a, -6) + pow(B_b, -6);
                   B_ab = 2. / B_ab;
                   B_ab = pow(B_ab, 1./6.);

                   // A parameter
                   A_ab  = sqrt(A_a*A_b);
                   A_ab *= pow(B_ab, 6);
                   A_ab /= pow(B_a, 3); 
                   A_ab /= pow(B_b, 3); 

                   // C, D and rcut parameters
                   C_ab = sqrt(C_a*C_b);
                   D_ab = sqrt(D_a*D_b);
                   rcut_ab    = sqrt(rcut_a*rcut_b);

                 } else {
                   fatal_error() << "Unexpected mixing rule (" << mixingRule << ") for Exp-6 potential for cross term (" << species->at(ita).name() << ", " << species->at(itb).name() << ")" << std::endl;
                 }

                 pot.m_params.set_exp6_parameters( A_ab, B_ab, C_ab, D_ab, rcut_ab );

               } else {
                 fatal_error() << "Atom type " << species->at(ita).name() << " and atom type " << species->at(itb).name() << " do not have the same potential. The cross term cannot be defined." << std::endl;
               }
               pot.m_params.update_ecut();
               potentials_for_pairs->m_potentials.push_back(pot);
             }
            
          }
        }
        
        next_line() >> count;
        ldbg << "Number of atoms = "<<count<<std::endl;        

        // Save box shape and dimensions
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
      
        // Read atom positions
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

          if( typeMol.empty() ) typeMol="MOL";
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
       
        // read atom velocity
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
          atom_data[cnt][field::vx] = onika::physics::make_quantity( vx , atom_speed_unit ).convert();
          atom_data[cnt][field::vy] = onika::physics::make_quantity( vy , atom_speed_unit ).convert();
          atom_data[cnt][field::vz] = onika::physics::make_quantity( vz , atom_speed_unit ).convert();
        }
        
        // Read charges
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
            atom_data[at_id][field::charge] = onika::physics::make_quantity( charge , atom_charge_unit ).convert();
          }
          process_kw_line( next_line() );
        }
        
        // Read intramolecular force field
        if( force_field_description )
        {
          uint64_t nb_force_fields = 0;
          next_line() >> nb_force_fields;
          ldbg << "Reading "<<nb_force_fields<<" force fields"<<std::endl;
          for(uint64_t cnt=0;cnt<nb_force_fields;cnt++)
          {
            auto iss = next_line();
            std::string ff_pot_kind;
            iss >> ff_pot_kind;
            if( ff_pot_kind == "bond_opls" )
            {
              double p1, p2, p3, p4;
              std::string t1, t2;
              std::string u1, u2, u3, u4;
              iss >> t1 >> t2 >> p1 >> u1 >> p2 >> u2 >> p3 >> u3 >> p4 >> u4;
              double r0 = onika::physics::make_quantity( p1 , read_unit(u1) ).convert(); if( std::isnan(r0) ) r0=0.0;
              double k  = onika::physics::make_quantity( p2 , read_unit(u2) ).convert(); if( std::isnan(k)  ) k=0.0;
              // p3 = onika::physics::make_quantity( p3 , read_unit(u3) ).convert(); if( std::isnan(p3) ) p3=0.0; // 2nd and 3rd parameter ignored
              // p4 = onika::physics::make_quantity( p4 , read_unit(u4) ).convert(); if( std::isnan(p4) ) p4=0.0;
              ldbg << "Intramolecular "<<ff_pot_kind<<" for "<<t1<<","<<t2<<" : r0="<<r0<<" , k="<<k << std::endl;
              potentials_for_bonds->m_bond_desc.push_back( BondPotential{ ff_pot_kind , {t1,t2} , std::make_shared<IntraMolecularBondOPLSFunctional>(k,r0) } );
            }
            else if( ff_pot_kind == "angle_opls" )
            {
              double p1, p2, p3, p4;
              std::string t1, t2, t3;
              std::string u1, u2, u3, u4;
              iss >> t1 >> t2 >> t3 >> p1 >> u1 >> p2 >> u2 >> p3 >> u3 >> p4 >> u4;
              double phi0 = onika::physics::make_quantity( p1 , read_unit(u1) ).convert(); if( std::isnan(phi0) ) phi0=0.0;
              double k    = onika::physics::make_quantity( p2 , read_unit(u2) ).convert(); if( std::isnan(k)    ) k=0.0;
              // p3 = onika::physics::make_quantity( p3 , read_unit(u3) ).convert(); if( std::isnan(p3) ) p3=0.0;  // 2nd and 3rd parameter ignored
              // p4 = onika::physics::make_quantity( p4 , read_unit(u4) ).convert(); if( std::isnan(p4) ) p4=0.0;
              ldbg << "Intramolecular "<<ff_pot_kind<<" for "<<t1<<","<<t2<<" : "<<t3<<" : phi0="<<phi0<<" , k="<<k << std::endl;
              potentials_for_angles->m_potentials.push_back( BendPotential{ ff_pot_kind , {t1,t2,t3} , std::make_shared<IntraMolecularBondOPLSFunctional>(k,phi0) } );
            }
            else if( ff_pot_kind == "torsion_opls" )
            {
              double p1, p2, p3;
              std::string t1, t2, t3, t4;
              std::string u1, u2, u3;
              iss >> t1 >> t2 >> t3 >> t4 >> p1 >> u1 >> p2 >> u2 >> p3 >> u3 ;
              p1 = onika::physics::make_quantity( p1 , read_unit(u1) ).convert(); if( std::isnan(p1) ) p1=0.0;
              p2 = onika::physics::make_quantity( p2 , read_unit(u2) ).convert(); if( std::isnan(p2) ) p2=0.0;
              p3 = onika::physics::make_quantity( p3 , read_unit(u3) ).convert(); if( std::isnan(p3) ) p3=0.0;
              ldbg << "Intramolecular "<<ff_pot_kind<<" for "<<t1<<","<<t2<<","<<t3<<","<<t4<<" : p1="<<p1<<" , p2="<<p2<<" , p3="<<p3<< std::endl;
              potentials_for_torsions->m_potentials.push_back( TorsionPotential{ ff_pot_kind , {t1,t2,t3,t4} , std::make_shared<IntraMolecularCosOPLSFunctional>(p1,p2,p3) } );
            }
            else
            {
              fatal_error() << "unknown force field type "<<ff_pot_kind<<std::endl;
            }
          }
          process_kw_line( next_line() );
        }
        
        // Read intramolecular connectivity
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
             
        // read intramolecular weight for non-bonding potentials
        while( !intramol_weight_description && !file.eof() )
        {
          process_kw_line( next_line() );
        }
        if( intramol_weight_description )
        {
          int n = 0;
          next_line() >> n;
          ldbg << "Reading "<<n<<" intramolecular pair weights"<<std::endl;
          mol_pair_weights->m_molecule_weight.clear();
          for(int i=0;i<n;i++)
          {
            std::string ta, tb;
            MolecularPairWeight w;
            next_line() >> ta >> tb >> w.m_bond_weight >> w.m_bend_weight >> w.m_torsion_weight >> w.m_rf_bond_weight >> w.m_rf_bend_weight >> w.m_rf_torsion_weight;
            ldbg << ta << " - " << tb << " weights : pair bond="<<w.m_bond_weight<<", angle="<<w.m_bend_weight<<", torsion="<<w.m_torsion_weight<<" , rf bond="<<w.m_rf_bond_weight<<", angle="<<w.m_rf_bend_weight<<", torsion="<<w.m_rf_torsion_weight<<std::endl;
            for(const auto& p : molecule_name_map)
            {
              const auto it = mol_pair_weights->m_molecule_weight.find(p.first);
              if( it != mol_pair_weights->m_molecule_weight.end() )
              {
                if( it->second.m_bond_weight != w.m_bond_weight || it->second.m_bend_weight != w.m_bend_weight || it->second.m_torsion_weight != w.m_torsion_weight
                  ||it->second.m_rf_bond_weight != w.m_rf_bond_weight || it->second.m_rf_bend_weight != w.m_rf_bend_weight || it->second.m_rf_torsion_weight != w.m_rf_torsion_weight )
                {
                  fatal_error() << "Different pair weights for distinct particle type pairs is not supported yet" << std::endl;
                }
              }
              else
              {
                ldbg << "Add intramolecular weights for molecule '"<<p.first<<"'"<<std::endl;
                mol_pair_weights->m_molecule_weight[ p.first ] = w;
              }
            }
          }
        }

        // close main file
        file.close();

        // if available, open and read Stamp init file to get reaction field parameters
        if( stampinput.has_value() )
        {
          //ldbg << "Init file = "<< *stampinput << std::endl;
          file.open( onika::data_file_path( *stampinput ) );
        }
        if( file.good() )
        {
          ldbg << "Reading additional init file "<< *stampinput << std::endl;
          lineno = -1;
          std::map< std::string , double > init;
          while( !file.eof() )
          {
            std::string k;
            double v = 0.0;
            next_line() >> k >> v;
            //ldbg << "READ "<<k<<" = "<<v<<std::endl;
            init[k] = v;
          }
          if( init["ReactionField"] == 1 )
          {
            double rf_epsilon = init["ReactionField_epsilon"];
            double rf_rc = init["ReactionField_rc"];
            rf_rc = EXASTAMP_QUANTITY( rf_rc * m ); // StampV4 has implicit meter unit for distances
            ldbg << "found RF parameters : epsilon="<<rf_epsilon<<" , rcut="<<rf_rc<<std::endl;
            std::set< std::pair<std::string,std::string> > existing_pair_pots;
            for( auto & pot : potentials_for_pairs->m_potentials )
            {
              existing_pair_pots.insert( { pot.m_type_a , pot.m_type_b } );
              existing_pair_pots.insert( { pot.m_type_b , pot.m_type_a } );
              pot.m_params.set_reaction_field_parameters( rf_rc , rf_epsilon );
              ldbg << "update RF parameters for pair "<<pot.m_type_a<<" / "<<pot.m_type_b<<", rc="<<rf_rc<<", epsilon="<<rf_epsilon<< std::endl;
            }
            unsigned int n_species = species->size();
            for(unsigned int i=0;i<n_species;i++)
            {
              for(unsigned int j=i;j<n_species;j++)
              {
                if( existing_pair_pots.find( { species->at(i).name() , species->at(j).name() } ) == existing_pair_pots.end() )
                {
                  LJExp6RFMultiParmsPair pot = { species->at(i).name() , species->at(j).name() , {} };
                  pot.m_params.set_reaction_field_parameters( rf_rc , rf_epsilon );
                  potentials_for_pairs->m_potentials.push_back( pot );
                  ldbg<<"add RF parameters for pair "<<pot.m_type_a<<" / "<<pot.m_type_b<<", rc="<<rf_rc<<", epsilon="<<rf_epsilon<< std::endl;
                }
              }
            }
          }
        }
        
        for(const auto& pot : potentials_for_pairs->m_potentials )
        {
          ldbg<<"potentials_for_pairs: "<<pot.m_type_a<<" / "<<pot.m_type_b<<" => " << pot.m_params << std::endl;
          *rcut_max = std::max( *rcut_max , std::max( pot.m_params.m_rcut , pot.m_params.m_rf.rc ) );
        }
        ldbg << "global pair potential rcut = "<< *rcut_max << std::endl;
        
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

        ldbg <<std::endl << "Atomes maille elementaire : "<<atom_data.size()<<std::endl;
        for(const auto& tp : atom_data)
        {
          const Vec3d r = { tp[field::rx] , tp[field::ry] , tp[field::rz] };
          IJK loc = grid->locate_cell( r ); // domain_periodic_location( *domain, r );
//          const Vec3d v = { tp[field::vx] , tp[field::vy] , tp[field::vz] };
//          ldbg << onika::format_string("\t%d : %s , Pos=(% .5e,% .5e,% .5e) , V=(% .5e,% .5e,% .5e)\n",tp[field::id], species->at(tp[field::type]).name(), r.x,r.y,r.z, v.x,v.y,v.z );
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

      // broadcast molecule parameters
      using MolIOExt = MoleculeOptionalHeaderIO<decltype(ldbg)>;
      double bond_max_dist = 0.0;
      double bond_max_stretch = 0.0;
      MolIOExt molecule_io = {
        bond_max_dist ,
        bond_max_stretch ,
        molecules.get_pointer() ,
        ldbg,
        potentials_for_bonds.get_pointer() ,
        potentials_for_angles.get_pointer() ,
        potentials_for_torsions.get_pointer() ,
        potentials_for_impropers.get_pointer() ,
        potentials_for_pairs.get_pointer() ,
        mol_pair_weights.get_pointer() };
      std::vector<uint8_t> serialize_buffer;
      unsigned long molecule_header_size = 0;
      if( rank == 0 )
      {
        molecule_header_size = molecule_io.serialize_molecule_header( serialize_buffer );
        if( molecule_header_size != serialize_buffer.size() ) fatal_error() << "serialize buffer size error" << std::endl;
      }
      MPI_Bcast(&molecule_header_size, 1, MPI_UNSIGNED_LONG, 0, *mpi);
      serialize_buffer.resize( molecule_header_size , 0 );
      MPI_Bcast(serialize_buffer.data(), molecule_header_size, MPI_CHARACTER, 0, *mpi);
      molecule_header_size = molecule_io.deserialize_molecule_header( serialize_buffer );
      if( molecule_header_size != serialize_buffer.size() ) fatal_error() << "deserialize buffer size error" << std::endl;

      int nmol = molecules->m_molecules.size();
//      MPI_Bcast(&nmol, 1, MPI_INT, 0, *mpi);
//      molecules->m_bridge_molecules.clear();
//      molecules->m_molecules.resize( nmol );
//      MPI_Bcast( molecules->m_molecules.data() , sizeof(MoleculeSpecies)*nmol , MPI_CHARACTER , 0 , *mpi );

      // Broadcast of IntraFFtypes
      int nFFtypes = -1;
      if( rank == 0) nFFtypes = ffnameVector->size(); // only meaningfull for rank==0
      MPI_Bcast(&nFFtypes, 1, MPI_INT, 0, *mpi);
      assert( nFFtypes >= 0 );
      static constexpr size_t FF_TYPE_NAME_MAX_LEN = 32;
      std::vector<char> ff_names( nFFtypes * FF_TYPE_NAME_MAX_LEN , '\0' );
      if( rank == 0)
      {
        for(size_t i=0;i<ffnameVector->size();i++)
        {
          std::strncpy( & ff_names[i*FF_TYPE_NAME_MAX_LEN] , ffnameVector->at(i).c_str() , FF_TYPE_NAME_MAX_LEN );
          ff_names[ (i+1)*FF_TYPE_NAME_MAX_LEN -1 ] = '\0';
        }
      }
      MPI_Bcast(ff_names.data(), ff_names.size(), MPI_CHAR, 0, *mpi);
      ffnameVector->clear();
      ffnameVector->resize(nFFtypes);
      for(int i=0;i<nFFtypes;i++)
      {
        ffnameVector->at(i) = & ff_names[i*FF_TYPE_NAME_MAX_LEN];
      }


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
  ONIKA_AUTORUN_INIT(read_fatomes_mol)
  {
    OperatorNodeFactory::instance()->register_factory("read_fatomes_mol", make_grid_variant_operator< ReadFAtomesMolecules >);
  }

}
