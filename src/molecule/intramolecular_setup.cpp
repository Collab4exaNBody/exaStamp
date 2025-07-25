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
    
    ADD_SLOT( BondsPotentialParameters     , potentials_for_bonds     , INPUT, OPTIONAL );
    ADD_SLOT( BendsPotentialParameters     , potentials_for_angles    , INPUT, OPTIONAL );
    ADD_SLOT( TorsionsPotentialParameters  , potentials_for_torsions  , INPUT, OPTIONAL );
    ADD_SLOT( ImpropersPotentialParameters , potentials_for_impropers , INPUT, OPTIONAL );
    ADD_SLOT( LJExp6RFMultiParms           , potentials_for_pairs     , INPUT, LJExp6RFMultiParms{} );
    
    ADD_SLOT( IntramolecularPairWeighting  , mol_pair_weights  , INPUT , IntramolecularPairWeighting{} );

    ADD_SLOT( GridT                        , grid              , INPUT_OUTPUT);
    ADD_SLOT( MoleculeSpeciesVector        , molecules         , INPUT_OUTPUT , REQUIRED , DocString{"Molecule descriptions"} );
    
    ADD_SLOT( bool                         , long_range_correction , INPUT_OUTPUT, true , DocString{"Compute long range corrections if set to true"} );    
    ADD_SLOT( MoleculeComputeParameterSet  , molecule_compute_parameters , INPUT_OUTPUT, DocString{"Intramolecular functionals' parameters"} );

  public:
    inline bool is_sink() const override final { return true; }
    
    inline void execute ()  override final
    {
      //static constexpr MoleculeGenericFuncParam null_param = {0.0,0.0,0.0,0.0};

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

//      molecule_compute_parameters->m_molecules.assign( nmol , MoleculeComputeParams{} );   
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

      // compute global energy correction term
      bool all_lj = true;
      bool all_exp6 = true;
      for(const auto& potelem : potentials_for_pairs->m_potentials)
      {
        const auto & pot = potelem.m_params;
        ldbg << "pot "<<potelem.m_type_a<<"/"<<potelem.m_type_b<<" : PAIR : A="<<pot.m_A<<", B="<<pot.m_B_ISEXP6<<", C="<<pot.m_C_EPSILON<<", D="<<pot.m_D_SIGMA<<", ecut="<<pot.m_ecut
             << ", RF : RF0="<<pot.m_rf.RF0<<", RF1="<<pot.m_rf.RF1<<", RF2="<<pot.m_rf.RF2<<", ecut="<<pot.m_rf.ecut<<std::endl;
        if( ! pot.pair_is_null() )
        {
          ldbg << "pot "<<potelem.m_type_a<<"/"<<potelem.m_type_b<<" is "<< ( pot.is_lj()?"LJ":"Exp6" ) << std::endl;
          all_lj = all_lj && pot.is_lj();
          all_exp6 = all_exp6 && pot.is_exp6();
        }
        else
        {
          ldbg << "pot "<<potelem.m_type_a<<"/"<<potelem.m_type_b<<" is empty"<<std::endl;
        }
      }

      float torsion_pair_weight=1.0f, torsion_rf_weight=1.0f, angle_pair_weight=1.0f, angle_rf_weight=1.0f, bond_pair_weight=1.0f, bond_rf_weight=1.0f;
      for( const auto & p : mol_pair_weights->m_molecule_weight )
      {
        ldbg << p.first << " : bond="<<p.second.m_bond_weight<<" , bond_rf="<<p.second.m_rf_bond_weight
                        <<" , bend="<<p.second.m_bend_weight<<" , bend_rf="<<p.second.m_rf_bend_weight
                        <<" , torsion="<<p.second.m_torsion_weight<<" , torsion_rf="<<p.second.m_rf_torsion_weight<<std::endl;
        bond_pair_weight = p.second.m_bond_weight;
        bond_rf_weight = p.second.m_rf_bond_weight;
        angle_pair_weight = p.second.m_bend_weight;
        angle_rf_weight = p.second.m_rf_bend_weight;
        torsion_pair_weight = p.second.m_torsion_weight;
        torsion_rf_weight = p.second.m_rf_torsion_weight;
      }

      const unsigned int n_type_pairs = unique_pair_count( species->size() );
      molecule_compute_parameters->m_pair_params.assign( n_type_pairs * 3 , IntramolecularPairParams{ LJExp6RFParms{}, 0.0f, 0.0f } );
      for(auto& potelem : potentials_for_pairs->m_potentials)
      {
        auto & pot = potelem.m_params;
        if( str2type(potelem.m_type_a)==-1 ) { fatal_error()<<"unknown type "<<potelem.m_type_a<<" in potential description"<<std::endl; }
        if( str2type(potelem.m_type_b)==-1 ) { fatal_error()<<"unknown type "<<potelem.m_type_b<<" in potential description"<<std::endl; }
        const unsigned int ta = str2type(potelem.m_type_a);
        const unsigned int tb = str2type(potelem.m_type_b); 
        const unsigned int pair_id = unique_pair_id(ta,tb);
        
        // if long_range_correction is enabled, we must set non-RF pair potenital's ecut to 0
        if( *long_range_correction )
        {
          pot.m_ecut = 0.0;
        }
        
        molecule_compute_parameters->m_pair_params[ n_type_pairs * 0 + pair_id ] = IntramolecularPairParams{ pot , bond_pair_weight, bond_rf_weight };
        molecule_compute_parameters->m_pair_params[ n_type_pairs * 1 + pair_id ] = IntramolecularPairParams{ pot , angle_pair_weight, angle_rf_weight };
        molecule_compute_parameters->m_pair_params[ n_type_pairs * 2 + pair_id ] = IntramolecularPairParams{ pot , torsion_pair_weight, torsion_rf_weight };
        
        const size_t na = type_count[ta]; // number of atom of type A
        const size_t nb = type_count[tb]; // number of atom of type B
        double Vtot = 1. ; // total volume of the simulation box
        double ecorr = 0.; // basis for the energy correction per atom 
	      // global formula for the long range energy correction due to the couple (ta, tb) for EACH atom of type ta:
	      // ecorr = int_rc^(+inf) r^2 * U(r) 4 pi/Vtot 1/2 Nb
	      double vcorr = 0.; //virial correction
        // global formula for the long range correction to EACH term of the diagonal component of the virial tensor due to the couple (ta, tb) for EACH atom of type ta:
        // vcorr = - int_rc^(+inf) r^3 * dU/dr 4 pi/Vtot 1/6 Nb

        // total volume calculation
	      if( ! domain->xform_is_identity() )
        {
          Mat3d mat = domain->xform();
          Vec3d a { mat.m11, mat.m21, mat.m31 };
          Vec3d b { mat.m12, mat.m22, mat.m32 };
          Vec3d c { mat.m13, mat.m23, mat.m33 };
          Vtot = dot( cross(a,b) , c );
        }
        Vtot *= bounds_volume( domain->bounds() );


	      ldbg << "Compute correction for pair "<<potelem.m_type_a<<" , "<<potelem.m_type_b<<std::endl;
        if( all_lj )
        {
          const double rcut = pot.m_rcut;
          const double epsilon = pot.m_C_EPSILON;
          const double sigma = pot.m_D_SIGMA;
          //molecule_compute_parameters->m_energy_correction[ta] += 0.0; // so something here ...
	        // U(r) = 4*epsilon ( (sigma/r)^12 - (sigma/r)^6)
	        // int_rc^(+inf) r^2 * U(r)  = 4 * epsilon (1/9 sigma^12/rc^9 - 1/3 sigma^6/rc^3)
	        ecorr = 4.*epsilon * (1./9. *pow(sigma, 12)/pow(rcut, 9) - 1./3.* pow(sigma, 6)/pow(rcut, 3));
          //molecule_compute_parameters->m_virial_correction[ta] += Mat3d{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // so something here ...
      	  // - int_rc^(+inf) r^3 * dU/dr = 4 * epsilon (4/3 sigma^12/rc^9 -2 sigma^6/rc^3)
          vcorr =  4.*epsilon * (4./3.*pow(sigma, 12)/pow(rcut, 9) - 2.*pow(sigma, 6)/pow(rcut, 3));
        }
        else if( all_exp6 )
        {
          const double rc = pot.m_rf.rc;
          const double A = pot.m_A;
          const double B = pot.m_B_ISEXP6;
          const double C = pot.m_C_EPSILON;
          const double D = pot.m_D_SIGMA;
          // U(r) = A exp(-B*r) - C/r^6 + D (12/(B*r))^12
          // int_rc^(+inf) r^2 * U(r) = A/B  exp(-B*rc) (rc^2 + 2 rc/B + 2/B^2) 
          //                           -C/3 1/rc^3
	        //                           +D (12/B)^12 1/9 1/rc^9
	        ecorr = A/B * exp(-B*rc)* (pow(rc, 2) + 2.* rc/B + 2./pow(B, 2)) -C/3. *pow(rc, -3) + D*pow(12./B, 12) *1./9. *pow(rc, -9);
          // - int_rc^(+inf) r^3 * dU/dr = - (A exp(-B*rc) * (rc^3 + 3 rc^2/B + 6 rc/B^2 + 6/B^3)
          //                                  + 2 C / rc^3
	        //                                  - 4/3 D (12/B)^12 1/rc^9)
          vcorr = - (A* exp(-B*rc) * (pow(rc, 3) + 3.*pow(rc, 2)/B + 6.*rc/pow(B, 2) + 6./pow(B,3))
		     + 2.*C / pow(rc, 3)
		     - 4./3. *D *pow(12./B, 12) * pow(rc,-9));
          //molecule_compute_parameters->m_energy_correction[ta] += 0.0; // so something here ...
          //molecule_compute_parameters->m_virial_correction[ta] += Mat3d{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // so something here ...
        }
        else
        {
          fatal_error() << "pair potentials must be all LJ or all Exp6" << std::endl;
        }

        //multiplication of the energy and virial corrections by factors common to the different potentials
	      ecorr *= 4. * M_PI / Vtot * 1./2. ;
	      vcorr *= 4. * M_PI / Vtot * 1./6. ;

	      if (*long_range_correction) {
	        // addition of the corrections due to the couple [ta, tb] to the corrections of atoms of type ta
	        molecule_compute_parameters->m_energy_correction[ta] += ecorr *nb;
	        molecule_compute_parameters->m_virial_correction[ta] += Mat3d{ vcorr*nb, 0.0, 0.0, 0.0, vcorr*nb, 0.0, 0.0, 0.0, vcorr*nb };
	      } else {
		      molecule_compute_parameters->m_energy_correction[ta] += 0.;
		      molecule_compute_parameters->m_virial_correction[ta] += Mat3d{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	      }

	      //In the case where ta is different from tb, addition of the corrections to atoms of type tb
	      // the couple (ta, tb) is encountered only once
        if (ta != tb) {
          if( *long_range_correction )
          {
	    molecule_compute_parameters->m_energy_correction[tb] += ecorr *na;
            molecule_compute_parameters->m_virial_correction[tb] += Mat3d{ vcorr*na, 0.0, 0.0, 0.0, vcorr*na, 0.0, 0.0, 0.0, vcorr*na };
          } else {
            molecule_compute_parameters->m_energy_correction[tb] += 0.;
	    molecule_compute_parameters->m_virial_correction[tb] += Mat3d{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	  }
        }


        for(size_t t=0;t<species->size();t++)
        {
          ldbg << "Correction for type "<<species->at(t).name()<<" : energy="<<molecule_compute_parameters->m_energy_correction[t]<<" , virial="<<molecule_compute_parameters->m_virial_correction[t]<<std::endl;
        }
      }
      /************************************************************/        

     
      auto& intramol_param_map = molecule_compute_parameters->m_intramol_param_map; // types to functional parameters index
      intramol_param_map.clear();      
      std::map< MoleculeGenericFuncParam , int > intramol_param_id_map ; // functional parameters to its index

      unsigned int parameter_id = 0;
      
      if( potentials_for_bonds.has_value() )
      for(const auto& bond : potentials_for_bonds->m_bond_desc)
      {
        int a = str2type( bond.species[0] );
        int b = str2type( bond.species[1] );
        const auto types = bond_key(a,b);
        ldbg << "read bond "<<bond.species[0]<<","<<bond.species[1]<<" -> "<<(const void*)types<<std::endl;
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
      }
      
      if( potentials_for_angles.has_value() )
      for(const auto& angle : potentials_for_angles->m_potentials)
      {
        int a = str2type( angle.species[0] );
        int b = str2type( angle.species[1] );
        int c = str2type( angle.species[2] );
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
      }
      
      if( potentials_for_torsions.has_value() )
      for(const auto& torsion : potentials_for_torsions->m_potentials)
      {
        int a = str2type( torsion.species[0] );
        int b = str2type( torsion.species[1] );
        int c = str2type( torsion.species[2] );
        int d = str2type( torsion.species[3] );
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
      }
      
      if( potentials_for_impropers.has_value() )
      for(const auto& improper : potentials_for_impropers->m_potentials)
      {
        int a = str2type( improper.species[0] );
        int b = str2type( improper.species[1] );
        int c = str2type( improper.species[2] );
        int d = str2type( improper.species[3] );
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
            
      // look for and initialize pair potentials
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

      for(unsigned int m=0;m<nmol;m++)
      {
        if( mol_pair_weights.has_value() )
        {
          if( mol_pair_weights->m_molecule_weight.find( molecules->m_molecules.at(m).name() ) == mol_pair_weights->m_molecule_weight.end() )
          {
            fatal_error() << "mol_pair_weights has no entry for molecule name '"<<molecules->m_molecules.at(m).name()<<"'" << std::endl;
          }
        }
      }
      
      // count number of different pair potential parameters
      std::map< std::pair<int,int> , LJExp6RFMultiParmsPair > pair_param_map;
      for( const auto & pp : potentials_for_pairs->m_potentials )
      {
        ldbg << "molecule pair : " << pp.m_type_a <<" / "<<pp.m_type_b<<std::endl;
        int ta = str2type(pp.m_type_a);
        int tb = str2type(pp.m_type_b);
        if( ta > tb ) std::swap(ta,tb);
        pair_param_map[ {ta,tb} ] = pp;
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

