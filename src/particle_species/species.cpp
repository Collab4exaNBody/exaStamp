#include <exaStamp/particle_species/particle_specie_yaml.h>
#include <onika/math/basic_types_stream.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/log.h>
#include <onika/yaml/yaml_utils.h>
#include <exanb/core/particle_type_id.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <algorithm>

namespace exaStamp
{
  using namespace exanb;

  struct ParticleSpeciesNode : public OperatorNode
  {
    ADD_SLOT( ParticleSpecies , species , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( ParticleTypeMap , particle_type_map , INPUT_OUTPUT );
    ADD_SLOT( ParticleTypeProperties , particle_type_properties , INPUT_OUTPUT , ParticleTypeProperties{} );
    ADD_SLOT( bool , verbose , INPUT , true );
    ADD_SLOT( bool , fail_if_empty , INPUT , false );

    inline void execute () override final
    {
    
      if( *fail_if_empty && species->empty() )
      {
        lerr<<"No species found, aborting"<<std::endl;
        std::abort();
      }
    
      if(*verbose) lout<<"======== Atom species ========"<<std::endl;

      if( species->empty() )
      {
        if(*verbose) lout<<"No species defined"<<std::endl;
      }

      std::sort( species->begin() , species->end() , [](const ParticleSpecie& a, const ParticleSpecie& b) -> bool { return a.m_rigid_atom_count < b.m_rigid_atom_count; } );

      std::unordered_map<std::string,unsigned int> type_map;
      for(unsigned int i=0;i<species->size();i++)
      {
        type_map[species->at(i).m_name] = i;
      }
      
      int first_rigid_molecule = -1;
      for(unsigned int a=0;a<species->size();a++)
      {
        ParticleSpecie& atom = species->at(a);
      	//lout << "Type #"<<a<<" = "<<atom.name()<<std::endl;
        if( atom.m_rigid_atom_count == 1 && first_rigid_molecule != -1 )
        {
          lerr << "single atom encountered after the first rigid molecule. rigid molecules must be at end of the list"<<std::endl;
          std::abort();
        }
        if( atom.m_rigid_atom_count < 1 )
        {
          lerr << "rigid atom count must be >= 1"<<std::endl;
          std::abort();
        }
        /*if( atom.m_rigid_atom_count != atom.m_rigid_atom_names.size() )
        {
          lerr << "insconsistent rigid atom names count. count="<<atom.m_rigid_atom_count<<", names=";
          for(const auto& s:atom.m_rigid_atom_names) lerr<<s<<" ";
          lerr<<std::endl;
          std::abort();          
        }*/
        for(unsigned int i=0;i<atom.m_rigid_atom_count;i++)
        {
          unsigned int j = type_map[ atom.rigid_atom_name(i) ];
          if( j >= species->size() )
          {
            lerr << "bad rigid atom name '"<<atom.rigid_atom_name(i)<<"'"<<std::endl;
            std::abort();
          }
          if( species->at(j).m_rigid_atom_count > 1 )
          {
            lerr << "rigid atom "<<atom.rigid_atom_name(i)<<"' is not a single atom"<<std::endl;
            std::abort();
          }
          atom.m_rigid_atoms[i].m_atom_type = j;
        }
        if( atom.m_rigid_atom_count==1 && atom.m_rigid_atoms[0].m_atom_type!=int(a) )
        {
          lerr<<"single atom "<<atom.name()<<" points to "<< atom.rigid_atom_name(0)<<" as its only rigid atom, should refer to itself instead"<<std::endl;
          std::abort();
        }
        if( atom.m_rigid_atom_count > 1 )
        {
          ldbg << "atome type "<< atom.name() <<" ("<<a<<") is rigid molecule, with "<<atom.m_rigid_atom_count<<" atoms"<<std::endl;
          atom.m_mass = 0.0;
          atom.m_z = 0;
          atom.m_charge = 0.0;
          for(unsigned int i=0;i<atom.m_rigid_atom_count;i++)
          {
            atom.m_mass += species->at(atom.m_rigid_atoms[i].m_atom_type).m_mass;
            atom.m_charge += species->at(atom.m_rigid_atoms[i].m_atom_type).m_charge;
          }
        }
        
        if( atom.m_rigid_atom_count > 1 && first_rigid_molecule==-1 ) { first_rigid_molecule = a; }
      }
      
      if(*verbose)
      {
        for(ParticleSpecie& atom : *species)
        {
          lout << std::setprecision(8);
          if( atom.m_rigid_atom_count <= 1 )
          {
            lout<<"- "<<atom.name()<<": { mass: "<<atom.m_mass<<" , z: "<<atom.m_z<<" , charge: "<<atom.m_charge;
            if(!atom.molecule_name().empty()) { lout<<" , molecule:"<<atom.molecule_name(); }
            lout<<" }"<<std::endl;
          }
          else
          {
            double max_radius = 0.0;
            for(unsigned int i=0;i<atom.m_rigid_atom_count;i++) max_radius = std::max( max_radius , norm2( atom.m_rigid_atoms[i].m_pos ) );

            lout<<"- "<<atom.name()<<":"<<std::endl
                <<"    mass: "<<atom.m_mass<<std::endl
                <<"    z: "<<atom.m_z<<std::endl
                <<"    charge: "<<atom.m_charge<<std::endl
                <<"    radius: "<<sqrt(max_radius)<<std::endl;
            if(!atom.molecule_name().empty()) { lout<<"    molecule:"<<atom.molecule_name()<<std::endl; }
            lout<<"    rigid_molecule:"<<std::endl;

            for(unsigned int i=0;i<atom.m_rigid_atom_count;i++)
            {
              lout<<"      - "<<atom.rigid_atom_name(i)<<": ["<<atom.m_rigid_atoms[i].m_pos<<"]"<<std::endl;
            }
          }
        }
        lout<<"=============================="<<std::endl<<std::endl;
      }
      
      // update generic (exaNBody level) particle type name map and particle type properties
      particle_type_map->clear();
      particle_type_properties->m_name_map.clear();
      particle_type_properties->m_scalars.clear();
      particle_type_properties->m_vectors.clear();
      particle_type_properties->m_names.assign( species->size() , "" );
      for(unsigned int a=0;a<species->size();a++)
      {
        const auto & spec = species->at(a);
        const std::string type_name = spec.name();
        (*particle_type_map) [ type_name ] = a;
        particle_type_properties->m_names[a] = type_name;
        
        particle_type_properties->m_name_map[ type_name ].m_scalars["mass"] = spec.m_mass;
        particle_type_properties->m_name_map[ type_name ].m_scalars["z"] = spec.m_z;
        particle_type_properties->m_name_map[ type_name ].m_scalars["charge"] = spec.m_charge;

        particle_type_properties->scalar_property("mass")[a] = spec.m_mass;
        particle_type_properties->scalar_property("z")[a] = spec.m_z;
        particle_type_properties->scalar_property("charge")[a] = spec.m_charge;
      }
            
    }

    inline void yaml_initialize(const YAML::Node& node) override final
    {
      YAML::Node tmp;
      if( node.IsSequence() )
      {
        tmp["species"] = node;
      }
      else { tmp = node; }
      this->OperatorNode::yaml_initialize(tmp);
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(species)
  {
    OperatorNodeFactory::instance()->register_factory( "species", make_simple_operator<ParticleSpeciesNode> );
  }

}

