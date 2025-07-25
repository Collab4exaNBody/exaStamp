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

#pragma once

#include <exanb/compute/compute_pair_function.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/grid_fields.h>
#include <onika/log.h>
#include <exanb/core/particle_type_pair.h>
#include <onika/flat_tuple.h>

#include <map>

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <functional>

namespace exaStamp
{
  using namespace exanb;

  struct PairPotentialResult
  {
    double E;
    double dE;
  };

  struct PairPotentialParameters
  {
    ParticleSpecie m_atom_a;
    ParticleSpecie m_atom_b;
  };

  struct PairPotentialMinimalParameters
  {
    MinimalParticleSpecies m_atom_a;
    MinimalParticleSpecies m_atom_b;
    PairPotentialMinimalParameters() = default;
    PairPotentialMinimalParameters(const PairPotentialMinimalParameters&) = default;
    PairPotentialMinimalParameters(PairPotentialMinimalParameters&&) = default;
    inline PairPotentialMinimalParameters(const PairPotentialParameters& sp) : m_atom_a(sp.m_atom_a), m_atom_b(sp.m_atom_b) {}
    PairPotentialMinimalParameters& operator = (const PairPotentialMinimalParameters&) = default;
    PairPotentialMinimalParameters& operator = (PairPotentialMinimalParameters&&) = default;
    inline PairPotentialMinimalParameters& operator = (const PairPotentialParameters& sp) { m_atom_a=sp.m_atom_a; m_atom_b=sp.m_atom_b; return *this; }    
  };

  using PairPotentialFieldSet = FieldSet< field::_ep, field::_fx, field::_fy, field::_fz >;
  using PairPotentialVirialFieldSet = FieldSet< field::_ep, field::_fx, field::_fy, field::_fz, field::_virial >;

  class PairPotentialComputeOperator : public ComputePairOperator<PairPotentialFieldSet>
  {
  public:
      virtual void set_parameters(const PairPotentialParameters& params) =0;
      virtual void set_rcut( double rcut ) =0;
      virtual uint64_t signature() const =0;
      virtual const std::string& name() const =0;
  };

  class PairPotentialComputeVirialOperator : public ComputePairOperator<PairPotentialVirialFieldSet>
  {
  public:
      virtual void set_parameters(const PairPotentialParameters& params) =0;
      virtual void set_rcut( double rcut ) =0;
      virtual uint64_t signature() const =0;
      virtual const std::string& name() const =0;
  };

  class PairPotential
  {
  public:
    virtual ~PairPotential() = default;
    virtual std::shared_ptr<PairPotentialComputeOperator> force_op();
    virtual std::shared_ptr<PairPotentialComputeVirialOperator> force_virial_op();
  };

  struct PairPotentialUserParameters
  {
    std::string m_type_a;
    std::string m_type_b;
    std::shared_ptr<PairPotential> m_pair_potential;
    double m_rcut = 0.0;
  };

  template<bool has_virial_field>
  struct CompiledMultiPairPotentials
  {
    using WorkFieldFieldSet = std::conditional_t< has_virial_field , PairPotentialVirialFieldSet , PairPotentialFieldSet >;  
    using BaseForceOp = ComputePairOperator<WorkFieldFieldSet>;
    using BaseForceOps = std::vector< std::shared_ptr<BaseForceOp> >;

    std::vector<double> m_rcuts;
    BaseForceOps m_force_ops;
    std::vector<int> m_pair_id_map;
  };

  struct UserMultiPairPotentials
  {
    std::vector<PairPotentialUserParameters> m_user_potentials;
    CompiledMultiPairPotentials<false> m_compiled_potentials;
    CompiledMultiPairPotentials<true> m_compiled_potentials_virial;
    struct PotentialCompileInfo { const double rcut_max; const unsigned int NTypes; };
    
    
    inline CompiledMultiPairPotentials<false>& compiled_potentials( std::false_type ) { return m_compiled_potentials; }
    inline CompiledMultiPairPotentials<true>& compiled_potentials( std::true_type ) { return m_compiled_potentials_virial; }
        
    template<class ZeroPot, bool yn>
    inline PotentialCompileInfo
    compile_potentials_for_species( const ParticleSpecies& species, ZeroPot zero_pot_vir, std::integral_constant<bool,yn> has_virial_field )
    {
      //using PairForceOps = std::conditional_t<has_virial_field, std::vector< std::shared_ptr<PairPotentialComputeVirialOperator> >, std::vector< std::shared_ptr<PairPotentialComputeOperator> > >;    
      //auto& ldbg = std::cout; // debug, override debug stream
      
      double rcut_max = 0.0;
      auto& user_potentials = m_user_potentials;
      auto& compiled_potentials = this->compiled_potentials( has_virial_field );
      
      unsigned int single_atom_types = 0;
      for(;single_atom_types<species.size() && species.at(single_atom_types).m_rigid_atom_count==1 ;single_atom_types++);
      ldbg << single_atom_types << " single atoms among "<<species.size()<<" species"<<std::endl;
      for(unsigned int i=single_atom_types;i<species.size() ;i++)
      {
        if( species.at(i).m_rigid_atom_count <= 1 )
        {
          lerr << "Bad ordering of species, "<<species.at(i).m_name<<" (idx "<<i<<") after single atoms has "<<species.at(i).m_rigid_atom_count<<" sub atoms"<<std::endl;
          std::abort();
        }
      }

      const unsigned int NTypes = single_atom_types;
      const unsigned int NTypePairs = unique_pair_count(NTypes);

      // update operator parameters from atom types and user potnetials' parameters
      if( compiled_potentials.m_force_ops.size() != NTypePairs
       || compiled_potentials.m_force_ops.size() != compiled_potentials.m_rcuts.size() )
      {
        compiled_potentials.m_rcuts.assign( NTypePairs, 0.0 );
        
        compiled_potentials.m_force_ops.clear();
        compiled_potentials.m_force_ops.resize( NTypePairs , nullptr );
        auto& force_ops = compiled_potentials.m_force_ops;
        
        std::map<std::string,unsigned int> typeIdMap;
        for(unsigned int i=0;i<NTypes;i++)
        {
          typeIdMap[ species.at(i).m_name ] = i;
        }

        const unsigned int n_input_pairs = user_potentials.size();
        for(unsigned int i=0;i<n_input_pairs;i++)
        {
          if( typeIdMap.find(user_potentials.at(i).m_type_a) == typeIdMap.end() )
          {
            lerr << "cannot find particle type '"<<user_potentials.at(i).m_type_a<<"'"<<std::endl;
            std::abort();
          }
          if( typeIdMap.find(user_potentials.at(i).m_type_b) == typeIdMap.end() )
          {
            lerr << "cannot find particle type '"<<user_potentials.at(i).m_type_b<<"'"<<std::endl;
            std::abort();
          }
          unsigned int ta = typeIdMap[ user_potentials.at(i).m_type_a ];
          unsigned int tb = typeIdMap[ user_potentials.at(i).m_type_b ];
          assert( ta < NTypes );
          assert( tb < NTypes );
          unsigned int pair_id = unique_pair_id(ta,tb);
          auto force_op = get_potential_force_op( * user_potentials.at(i).m_pair_potential , has_virial_field );          
          force_op->set_parameters( PairPotentialParameters{ species.at(ta) , species.at(tb) } );
          compiled_potentials.m_rcuts[pair_id] = user_potentials.at(i).m_rcut;
          force_op->set_rcut( compiled_potentials.m_rcuts[pair_id] );          
          rcut_max = std::max( compiled_potentials.m_rcuts[pair_id] , rcut_max );
          if( force_ops[pair_id] != nullptr )
          {
            lerr << "Fatal: potential for types "<<user_potentials.at(i).m_type_a<<" and "<<user_potentials.at(i).m_type_b<<" duplicated"<<std::endl;
            std::abort();
          }
          force_ops[pair_id] = force_op;
        }

        for(unsigned int i=0;i<NTypePairs;i++)
        {
          if( force_ops[i] == nullptr )
          {
            unsigned int ta=0,tb=0;
            pair_id_to_type_pair(i,ta,tb);
            ldbg << "automatic zero potential for pair "<<species.at(ta).m_name<<" / "<<species.at(tb).m_name<<std::endl;
            force_ops[i] = zero_pot_vir;
            compiled_potentials.m_rcuts[i] = 0.0;
          }
        }
              
        // ---- here, all potentials are initialized ----
        compiled_potentials.m_pair_id_map.clear();        
        assert( compiled_potentials.m_force_ops.size() == compiled_potentials.m_rcuts.size() );

        ldbg << "rcut_max = "<< rcut_max << std::endl;
      }

      return { rcut_max , NTypes };
    }
  };


  template<class PotentialParamT , class PairPotExtractT >
  inline std::vector<char> serialize_pair_potential(double rcut, const PotentialParamT& pot_params, const PairPotExtractT& pair_params)
  {
    static_assert(sizeof(char)==1,"Assumption that sizeof(char) is 1 failed");
    auto pp = onika::make_flat_tuple( rcut , pair_params );

    std::vector<char> signature_buffer;
    signature_buffer.reserve( sizeof(pot_params) + sizeof(pp) );
    
    // add data from the specific pair potential parameter structure
    const char* sptr = reinterpret_cast<const char*>( &pot_params );
    signature_buffer.assign( sptr, sptr+sizeof(pot_params) );
    
    // add data from the atom type pair parameters
    sptr = reinterpret_cast<const char*>( &pp );
    signature_buffer.insert( signature_buffer.end() , sptr , sptr+sizeof(pp) );
        
    return signature_buffer;
  }

  template<class PotentialParamT , class PairPotExtractT >
  inline uint64_t hash_pair_potential(double rcut, const PotentialParamT& pot_params, const PairPotExtractT& pair_params)
  {
    auto serial = serialize_pair_potential(rcut,pot_params,pair_params);
    return std::hash<std::string_view>{}( std::string_view( serial.data() , serial.size() ) );
  }

}


