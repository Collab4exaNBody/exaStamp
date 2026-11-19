/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <onika/math/basic_types.h>
#include <exanb/core/particle_id_codec.h>
#include <onika/log.h>
#include <exanb/core/particle_type_id.h>
#include <exanb/core/domain.h>

#include <exaStamp/particle_species/particle_specie.h>

#include <exaStamp/molecule/mol_connectivity.h>
#include <exaStamp/molecule/molecule_species.h>
#include <exaStamp/molecule/molecule_compute_param.h>

#include <exaStamp/molecule/bonds_potentials_parameters.h>
#include <exaStamp/molecule/bends_potentials_parameters.h>
#include <exaStamp/molecule/torsions_potentials_parameters.h>
#include <exaStamp/molecule/impropers_potentials_parameters.h>
#include <exaStamp/molecule/id_map.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/molecule/iFFname.h>

#include <exaStamp/molecule/intramolecular_pair_weight.h>
#include <exaStamp/potential/ljexp6rf/ljexp6rf.h>

#include <onika/oarray_stream.h>

#include <unordered_set>
#include <vector>
#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;
  
  template< class GridT >
  class IntramolecularSetup : public OperatorNode
  {
    ADD_SLOT( MPI_Comm                     , mpi               , INPUT , MPI_COMM_WORLD);

    ADD_SLOT( Domain                       , domain            , INPUT , REQUIRED );
    ADD_SLOT( ParticleSpecies              , species           , INPUT , REQUIRED );
    ADD_SLOT( ParticleTypeMap              , particle_type_map , INPUT , REQUIRED );
    
    ADD_SLOT( BondsPotentialParameters     , potentials_for_bonds     , INPUT, OPTIONAL , DocString{"Bond potentials of the force field"} );
    ADD_SLOT( BendsPotentialParameters     , potentials_for_angles    , INPUT, OPTIONAL , DocString{"Bending angle potentials of the force field"} );
    ADD_SLOT( TorsionsPotentialParameters  , potentials_for_torsions  , INPUT, OPTIONAL , DocString{"Torsion potentials of the force field"} );
    ADD_SLOT( ImpropersPotentialParameters , potentials_for_impropers , INPUT, OPTIONAL , DocString{"Improper torsion potentials of the force field"} );
    ADD_SLOT( LJExp6RFMultiParms           , potentials_for_pairs     , INPUT_OUTPUT, LJExp6RFMultiParms{} , DocString{"Pair potentials (LJ or Exp6, plus reaction field) for each pair of species. With long_range_correction, their energy shift at the cutoff is removed"} );
    ADD_SLOT( IntraFFtypes                 , ffnameVector             , INPUT, OPTIONAL , DocString{"Force field type names, indexed by force field type id"} );
    
    ADD_SLOT( IntramolecularPairWeighting  , mol_pair_weights  , INPUT , IntramolecularPairWeighting{} , DocString{"Intramolecular pair weights for each pair of species"} );

    ADD_SLOT( GridT                        , grid              , INPUT_OUTPUT);
    ADD_SLOT( MoleculeSpeciesVector        , molecules         , INPUT_OUTPUT , REQUIRED , DocString{"Molecule descriptions"} );
    
    ADD_SLOT( bool                         , long_range_correction , INPUT_OUTPUT, true , DocString{"Compute long range corrections if set to true"} );    
    ADD_SLOT( MoleculeComputeParameterSet  , molecule_compute_parameters , INPUT_OUTPUT, DocString{"Intramolecular functionals' parameters"} );

  public:
    inline bool is_sink() const override final { return true; }
    
    inline std::string documentation() const override final
    {
      return R"EOF(
Prepares the parameters of intramolecular terms (molecule_compute_parameters) : parameter sets of bonds, angles, torsions
and impropers indexed by force field types, and pair parameters with their bond, bend and torsion weights (mol_pair_weights).
With long_range_correction, computes the long range correction of the energy and virial per atom type (for LJ or Exp6
pair potentials beyond their cutoff), applied with the current box volume by zero_force_energy_with_correction, and
removes the energy shift of pair potentials at their cutoff.
)EOF";
    }

    inline void execute ()  override final
    {
      const unsigned int nmol = molecules->m_molecules.size();
      ldbg << "Number of molecule species : "<<nmol<<std::endl;
      if( nmol == 1 && molecules->m_molecules.front().name().empty() ) molecules->m_molecules.front().set_name("MOL");
      for(unsigned int m=0;m<nmol;m++)
      {
        ldbg << "molecule #"<<m<<" : '"<<molecules->m_molecules.at(m).name()<<"'" << std::endl;
        if( ! molecules->m_molecules[m].has_connectivity() )
        {
          ldbg << "Update connectivity for molecule #"<<m<< std::endl;
          molecules->m_molecules[m].update_connectivity();
          molecules->m_molecules[m].print( ldbg , *species );
        }
      }

      const auto & tmap = *particle_type_map;
      auto str2type = [&tmap]( const std::string& s ) -> int
        {
          auto it=tmap.find(s);
          if(it!=tmap.end()) return it->second;
          else return -1;
        };


      /********* global energy and virial correction ***************/
      molecule_compute_parameters->m_energy_correction.assign( species->size() , 0.0 );
      molecule_compute_parameters->m_virial_correction.assign( species->size() , Mat3d{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0} );
      const size_t n_cells = grid->number_of_cells();
      // count number of atoms of each type
      std::vector<long> type_count( species->size() , 0 );
      for(size_t ci=0;ci<n_cells;ci++)
      {
        if( ! grid->is_ghost_cell(ci) )
        {
          const auto& cell = grid->cell(ci);
          size_t n_particles = cell.size();
          for(size_t pi=0;pi<n_particles;pi++)
          {
            int t = cell[field::type][pi];
            ++ type_count[t];
          }
        }
      }
      MPI_Allreduce(MPI_IN_PLACE,type_count.data(), type_count.size(), MPI_LONG, MPI_SUM, *mpi );
      for(size_t t=0;t<species->size();t++) ldbg << "Type "<<species->at(t).name()<<" has "<<type_count[t]<<" atoms"<<std::endl;

      // the long range correction needs the same functional form (LJ or Exp6) for all pairs with a pair term
      bool all_lj = true;
      bool all_exp6 = true;
      for(const auto& potelem : potentials_for_pairs->m_potentials)
      {
        const auto & pot = potelem.m_params;
        ldbg << "pot "<<potelem.m_type_a<<"/"<<potelem.m_type_b<<" : PAIR : A="<<pot.m_A<<", B="<<pot.m_B_ISEXP6<<", C="<<pot.m_C_EPSILON<<", D="<<pot.m_D_SIGMA<<", ecut="<<pot.m_ecut
             << ", RF : RF0="<<pot.m_rf.RF0<<", RF1="<<pot.m_rf.RF1<<", RF2="<<pot.m_rf.RF2<<", ecut="<<pot.m_rf.ecut<<std::endl;
        if( pot.has_pair_term() )
        {
          all_lj = all_lj && pot.is_lj();
          all_exp6 = all_exp6 && pot.is_exp6();
        }
      }

      double torsion_pair_weight=1.0, torsion_rf_weight=1.0, angle_pair_weight=1.0, angle_rf_weight=1.0, bond_pair_weight=1.0, bond_rf_weight=1.0;

      const unsigned int n_type_pairs = unique_pair_count( species->size() );
      molecule_compute_parameters->m_pair_params.assign( n_type_pairs * 3 , IntramolecularPairParams{ LJExp6RFParms{}, 0.0, 0.0 } );
      for(auto& potelem : potentials_for_pairs->m_potentials)
      {
        auto & pot = potelem.m_params;

        if( str2type(potelem.m_type_a)==-1 ) { fatal_error()<<"unknown type "<<potelem.m_type_a<<" in potential description"<<std::endl; }
        if( str2type(potelem.m_type_b)==-1 ) { fatal_error()<<"unknown type "<<potelem.m_type_b<<" in potential description"<<std::endl; }
        const unsigned int ta = str2type(potelem.m_type_a);
        const unsigned int tb = str2type(potelem.m_type_b); 
        const unsigned int pair_id = unique_pair_id(ta,tb);

        auto it = mol_pair_weights->m_molecule_weight.find( potelem.m_type_a + "_" + potelem.m_type_b );
        if( it == mol_pair_weights->m_molecule_weight.end() ) it = mol_pair_weights->m_molecule_weight.find( potelem.m_type_b + "_" + potelem.m_type_a );
        if( it == mol_pair_weights->m_molecule_weight.end() )
        {
          fatal_error() << "intramolecular weight for pair (" << potelem.m_type_a << ", " << potelem.m_type_b << ") not found !" << std::endl;
        }
        bond_pair_weight = it->second.m_bond_weight;
        bond_rf_weight = it->second.m_rf_bond_weight;
        angle_pair_weight = it->second.m_bend_weight;
        angle_rf_weight = it->second.m_rf_bend_weight;
        torsion_pair_weight = it->second.m_torsion_weight;
        torsion_rf_weight = it->second.m_rf_torsion_weight;

        // with the long range correction, the pair potential (not the reaction field) is not shifted at its cutoff
        if( *long_range_correction )
        {
          pot.m_ecut = 0.0;
        }
        
        molecule_compute_parameters->m_pair_params[ n_type_pairs * 0 + pair_id ] = IntramolecularPairParams{ pot , bond_pair_weight, bond_rf_weight };
        molecule_compute_parameters->m_pair_params[ n_type_pairs * 1 + pair_id ] = IntramolecularPairParams{ pot , angle_pair_weight, angle_rf_weight };
        molecule_compute_parameters->m_pair_params[ n_type_pairs * 2 + pair_id ] = IntramolecularPairParams{ pot , torsion_pair_weight, torsion_rf_weight };

        // no correction for a pair without LJ or Exp6 term (rf-only pair)
        if( ! *long_range_correction || ! pot.has_pair_term() ) continue;

        // per atom of type ta, from the Nb atoms of type tb beyond rc, times the box volume V (divided by the current
        // volume in zero_force_energy_with_correction):
        // ecorr = 2 pi Nb int_rc^inf r^2 U(r) dr , vcorr (each diagonal term) = - 2/3 pi Nb int_rc^inf r^3 dU/dr dr
        double ecorr = 0.;
        double vcorr = 0.;
        if( all_lj )
        {
          const double rcut = pot.m_rcut;
          const double epsilon = pot.m_C_EPSILON;
          const double sigma = pot.m_D_SIGMA;
          // U(r) = 4 epsilon ( (sigma/r)^12 - (sigma/r)^6 )
          ecorr = 4.*epsilon * (1./9. *pow(sigma, 12)/pow(rcut, 9) - 1./3.* pow(sigma, 6)/pow(rcut, 3));
          vcorr = 4.*epsilon * (4./3.*pow(sigma, 12)/pow(rcut, 9) - 2.*pow(sigma, 6)/pow(rcut, 3));
        }
        else if( all_exp6 )
        {
          const double rc = pot.m_rcut;
          const double A = pot.m_A;
          const double B = pot.m_B_ISEXP6;
          const double C = pot.m_C_EPSILON;
          const double D = pot.m_D_SIGMA;
          // U(r) = A exp(-B r) - C/r^6 + D (12/(B r))^12
          ecorr = A/B * exp(-B*rc)* (pow(rc, 2) + 2.* rc/B + 2./pow(B, 2)) -C/3. *pow(rc, -3) + D*pow(12./B, 12) *1./9. *pow(rc, -9);
          vcorr = - (A* exp(-B*rc) * (pow(rc, 3) + 3.*pow(rc, 2)/B + 6.*rc/pow(B, 2) + 6./pow(B,3))
                     + 2.*C / pow(rc, 3)
                     - 4./3. *D *pow(12./B, 12) * pow(rc,-9));
        }
        else
        {
          fatal_error() << "pair potentials must be all LJ or all Exp6" << std::endl;
        }
        ecorr *= 4. * M_PI * 1./2. ;
        vcorr *= 4. * M_PI * 1./6. ;

        // the couple (ta, tb) is listed once: correct atoms of type ta, and atoms of type tb if different
        const double nb = type_count[tb];
        const double na = type_count[ta];
        molecule_compute_parameters->m_energy_correction[ta] += ecorr * nb;
        molecule_compute_parameters->m_virial_correction[ta] += Mat3d{ vcorr*nb, 0.0, 0.0, 0.0, vcorr*nb, 0.0, 0.0, 0.0, vcorr*nb };
        if( ta != tb )
        {
          molecule_compute_parameters->m_energy_correction[tb] += ecorr * na;
          molecule_compute_parameters->m_virial_correction[tb] += Mat3d{ vcorr*na, 0.0, 0.0, 0.0, vcorr*na, 0.0, 0.0, 0.0, vcorr*na };
        }
      }

      for(size_t t=0;t<species->size();t++)
      {
        ldbg << "Correction for type "<<species->at(t).name()<<" (times volume) : energy="<<molecule_compute_parameters->m_energy_correction[t]<<" , virial="<<molecule_compute_parameters->m_virial_correction[t]<<std::endl;
      }
      /************************************************************/        

     
      auto& intramol_param_map = molecule_compute_parameters->m_intramol_param_map; // types to functional parameters index
      intramol_param_map.clear();
      molecule_compute_parameters->m_null_param_keys.clear();
      molecule_compute_parameters->m_has_improper_potentials = potentials_for_impropers.has_value() && ! potentials_for_impropers->m_potentials.empty();      
      std::map< MoleculeGenericFuncParam , int > intramol_param_id_map ; // functional parameters to its index

      unsigned int parameter_id = 0;

      // bond (angle, torsion and improper).species contains FFname, it is converted to FFtypeId by ffnameVector (constructed such that ffnameVector[FFtypeId] = FFname).
      // To access the info easily, we construct locally the reverse map, ffnameToFFtypeId
      std::map< std::string , int64_t > ffnameToFFtypeId;
      if (ffnameVector.has_value()) {
        for (size_t iFF = 0; iFF < ffnameVector->size(); iFF++) {
          ffnameToFFtypeId[ ffnameVector->at(iFF) ] = iFF;
        }
      }
      auto ff_type_id = [&ffnameToFFtypeId]( const std::string& ffname ) -> int
      {
        const auto it = ffnameToFFtypeId.find( ffname );
        if( it == ffnameToFFtypeId.end() ) { fatal_error() << "unknown force field type '"<<ffname<<"' in intramolecular potential description" << std::endl; }
        return it->second;
      };


      
      if( potentials_for_bonds.has_value() )
      for(const auto& bond : potentials_for_bonds->m_bond_desc)
      {
        int a=0, b=0;
        // bond.species contains FFname, it is converted to FFtypeId by ffnameToFFtypeId.
        a = ff_type_id( bond.species[0] );
        b = ff_type_id( bond.species[1] );

        const auto types = bond_key(a,b);
        ldbg << "read bond " << bond.species[0] << ","<< bond.species[1] <<" -> "<<(const void*)types<<std::endl;
        auto param = bond.potential->generic_parameters();
        if( ! param.is_null() )
        {
          if( intramol_param_id_map.find(param)==intramol_param_id_map.end() )
          {
            ldbg<<"bond: add parameter pack "<<param<<" -> "<<parameter_id<<" , for key "<<(const void*)types<<std::endl;
            intramol_param_id_map[param] = parameter_id++;
          }
          intramol_param_map[ types ] = intramol_param_id_map[param];
        }
        else
        {
          molecule_compute_parameters->m_null_param_keys.insert( types ); // explicitly null : term is intentionally ignored
        }
      }
      
      if( potentials_for_angles.has_value() )
      for(const auto& angle : potentials_for_angles->m_potentials)
      {
        int a=0, b=0, c=0;

        // angle.species contains FFname, it is converted to FFtypeId by ffnameToFFtypeId.
        a = ff_type_id( angle.species[0] );
        b = ff_type_id( angle.species[1] );
        c = ff_type_id( angle.species[2] );

        const auto types = angle_key(a,b,c);
        ldbg << "read angle "<<angle.species[0]<<","<<angle.species[1]<<","<<angle.species[2]<<" -> "<<(const void*)types<<std::endl;
        auto param = angle.m_potential_function->generic_parameters();
        if( ! param.is_null() )
        {
          if( intramol_param_id_map.find(param)==intramol_param_id_map.end() )
          {
            ldbg<<"angle: add parameter pack "<<param<<" -> "<<parameter_id<<" , for key "<<(const void*)types<<std::endl;
            intramol_param_id_map[param] = parameter_id++;
          }
          intramol_param_map[ types ] = intramol_param_id_map[param];
        }
        else
        {
          molecule_compute_parameters->m_null_param_keys.insert( types ); // explicitly null : term is intentionally ignored
        }
      }
      
      if( potentials_for_torsions.has_value() )
      for(const auto& torsion : potentials_for_torsions->m_potentials)
      {
        int a=0, b=0, c=0, d=0;

        // torsion.species contains FFname, it is converted to FFtypeId by ffnameToFFtypeId.
        a = ff_type_id( torsion.species[0] );
        b = ff_type_id( torsion.species[1] );
        c = ff_type_id( torsion.species[2] );
        d = ff_type_id( torsion.species[3] );

        const auto types = torsion_key(a,b,c,d);
        ldbg << "read torsion "<<torsion.species[0]<<","<<torsion.species[1]<<","<<torsion.species[2]<<","<<torsion.species[3]<<" -> "<<(const void*)types<<std::endl;
        auto param = torsion.m_potential_function->generic_parameters();
        if( ! param.is_null() )
        {
          if( intramol_param_id_map.find(param)==intramol_param_id_map.end() )
          {
            ldbg<<"torsion: add parameter pack "<<param<<" -> "<<parameter_id<<" , for key "<<(const void*)types<<std::endl;
            intramol_param_id_map[param] = parameter_id++;
          }
          intramol_param_map[ types ] = intramol_param_id_map[param];
        }
        else
        {
          molecule_compute_parameters->m_null_param_keys.insert( types ); // explicitly null : term is intentionally ignored
        }
      }
      
      if( potentials_for_impropers.has_value() )
      for(const auto& improper : potentials_for_impropers->m_potentials)
      {
        int a=0, b=0, c=0, d=0;

        // torsion.species contains FFname, it is converted to FFtypeId by ffnameToFFtypeId.
        a = ff_type_id( improper.species[0] );
        b = ff_type_id( improper.species[1] );
        c = ff_type_id( improper.species[2] );
        d = ff_type_id( improper.species[3] );

        const auto types = improper_key(a,b,c,d);
        ldbg << "read improper "<<improper.species[0]<<","<<improper.species[1]<<","<<improper.species[2]<<","<<improper.species[3]<<" -> "<<(const void*)types<<std::endl;
        auto param = improper.m_potential_function->generic_parameters();
        if( ! param.is_null() )
        {
          if( intramol_param_id_map.find(param)==intramol_param_id_map.end() )
          {
            ldbg<<"improper: add parameter pack "<<param<<" -> "<<parameter_id<<" , for key "<<(const void*)types<<std::endl;
            intramol_param_id_map[param] = parameter_id++;
          }
          intramol_param_map[ types ] = intramol_param_id_map[param];
        }
        else
        {
          molecule_compute_parameters->m_null_param_keys.insert( types ); // explicitly null : term is intentionally ignored
        }
      }

      unsigned int nbparams = intramol_param_id_map.size();
      ldbg << nbparams << " different parameter sets"<<std::endl;
      molecule_compute_parameters->m_func_params.resize( nbparams );
      for(const auto &p : intramol_param_id_map)
      {
        unsigned int pidx = p.second;
        assert( pidx < nbparams );
        molecule_compute_parameters->m_func_params[ pidx ] = p.first;
      }
            
      // print the intramolecular pair weights
      if( mol_pair_weights.has_value() )
      {
        ldbg << "Intramolecular pair weighting map :"<<std::endl;      
        for( const auto & p : mol_pair_weights->m_molecule_weight )
        {
          ldbg << p.first << " : bond="<<p.second.m_bond_weight<<" , bond_rf="<<p.second.m_rf_bond_weight
                          <<" , bend="<<p.second.m_bend_weight<<" , bend_rf="<<p.second.m_rf_bend_weight
                          <<" , torsion="<<p.second.m_torsion_weight<<" , torsion_rf="<<p.second.m_rf_torsion_weight<<std::endl;
        }
      }

      ldbg << molecule_compute_parameters->m_pair_params.size()<<" pair parameters and "<<molecule_compute_parameters->m_func_params.size()<<" intra parameters"<<std::endl;     
    }

  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(intramolecular_setup)
  {
    OperatorNodeFactory::instance()->register_factory( "intramolecular_setup", make_grid_variant_operator< IntramolecularSetup > );
  }

}

