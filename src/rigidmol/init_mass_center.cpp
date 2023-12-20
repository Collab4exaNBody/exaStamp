#include <exaStamp/particle_species/particle_specie_yaml.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/basic_types_operators.h>

#include <exanb/core/operator.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/log.h>
#include <exanb/core/yaml_utils.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <unordered_map>

namespace exaStamp
{
  using namespace exanb;

  struct InitMassCenter : public OperatorNode
  {
    ADD_SLOT( ParticleSpecies , species , INPUT_OUTPUT , REQUIRED );

    inline void execute () override final
    {
      for(unsigned int a=0;a<species->size();a++)
      {
        ParticleSpecie& mol = species->at(a);
        if( mol.m_rigid_atom_count < 1 )
        {
          lerr << "rigid atom count must be >= 1"<<std::endl;
          std::abort();
        }
        /*if( mol.m_rigid_atom_count != mol.m_rigid_atom_names.size() )
        {
          lerr << "insconsistent rigid atom names count. count="<<mol.m_rigid_atom_count<<", names=";
          for(const auto& s:mol.m_rigid_atom_names) lerr<<s<<" ";
          lerr<<std::endl;
          std::abort();          
        }*/
        
        Vec3d masscenter_pos { 0.,0.,0. };
        double totalmass = 0.;
        for(unsigned int i=0;i<mol.m_rigid_atom_count;i++)
        {
          size_t site_atom_type = mol.m_rigid_atoms[i].m_atom_type;
          if( site_atom_type >= species->size() )
          {
            lerr << "bad rigid molecule site atom type "<<site_atom_type <<std::endl;
            std::abort();
          }
          double site_mass = species->at( site_atom_type ).m_mass;
          masscenter_pos += mol.m_rigid_atoms[i].m_pos * site_mass;
          totalmass += site_mass;
        }
        masscenter_pos /= totalmass;
        for(unsigned int i=0;i<mol.m_rigid_atom_count;i++)
        {
          //int site_atom_type = mol.m_rigid_atoms[i].m_atom_type;
          mol.m_rigid_atoms[i].m_pos -= masscenter_pos;
        }
      }
    }

  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "init_mass_center", make_simple_operator< InitMassCenter > );
  }

}

