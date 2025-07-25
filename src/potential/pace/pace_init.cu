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

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exaStamp/particle_species/particle_specie.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <onika/file_utils.h>
#include <exanb/core/particle_type_id.h>

#include <iostream>

#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"
#include "ace-evaluator/ace_recursive.h"
#include "ace-evaluator/ace_version.h"
#include "ace/ace_b_basis.h"

#include "pace_params.h"
#include "pace_config.h"

namespace exaStamp
{

  bool hasExtension_bis(const std::string& filename, const std::string& extension) {
    if (filename.length() >= extension.length()) {
      return std::equal(extension.rbegin(), extension.rend(), filename.rbegin());
    }
    return false;
  }

  static char const *const elements_pace_bis[] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
    "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};
  static constexpr int elements_num_pace = sizeof(elements_pace_bis) / sizeof(const char *);
  
  static int AtomicNumberByName_pace_bis(const char *elname)
  {
    for (int i = 1; i < elements_num_pace; i++)
      if (strcmp(elname, elements_pace_bis[i]) == 0) return i;
    return -1;
  }
  
  using namespace exanb;
  
  class PaceInit : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( PaceParams            , parameters        , INPUT , REQUIRED );
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( ParticleSpecies, species, INPUT, REQUIRED );
    ADD_SLOT( PaceContext           , pace_ctx          , OUTPUT );
        
  public:
    
    // Operator execution
    inline void execute () override final
    {

      std::cout << "Entering ACE potential initialization" << std::endl;
      bool cTildeBasis = false;
      (*pace_ctx).recursive = (*parameters).recursive;
      (*pace_ctx).aceimpl = new ACEImpl;
      (*pace_ctx).aceimpl->basis_set = new ACECTildeBasisSet;
      (*pace_ctx).aceimpl->ace = new ACERecursiveEvaluator();

      ACEBBasisSet bBasisSet;
      ACECTildeBasisSet cTildeBasisSet;

      auto potential_file_name = (*parameters).pace_coef;
      const auto& sp = *species;
      const int nspecies = sp.size();
      (*pace_ctx).nspecies = nspecies;
      
      if (hasExtension_bis(potential_file_name, ".yaml")) {
        bBasisSet = ACEBBasisSet(potential_file_name);
        cTildeBasisSet = bBasisSet.to_ACECTildeBasisSet();
        *(*pace_ctx).aceimpl->basis_set = cTildeBasisSet;
        cTildeBasis = true;
      } else {
        cTildeBasis = false;          
        *(*pace_ctx).aceimpl->basis_set = ACECTildeBasisSet(potential_file_name);
      }
      (*pace_ctx).aceimpl->ace->set_recursive((*pace_ctx).recursive);
      (*pace_ctx).aceimpl->ace->element_type_mapping.init((*pace_ctx).nspecies + 1);
      
      for (int i = 1; i <= nspecies; i++) {
        const char *elemname = sp[i-1].m_name;
        int atomic_number = AtomicNumberByName_pace_bis(elemname);
        if (atomic_number == -1) std::cout << elemname << "is not a valid element" << std::endl;
        SPECIES_TYPE mu = (*pace_ctx).aceimpl->basis_set->get_species_index_by_name(elemname);
        if (mu != -1) {
          std::cout << "Mapping LAMMPS atom type #"<< i << "("<<elemname<<") -> ACE species type #"<< mu << std::endl;
          (*pace_ctx).aceimpl->ace->element_type_mapping(i) = mu;
        } else {
          std::cout << "Element "<< elemname << " is not supported by ACE-potential from file " << potential_file_name << std::endl;
        }
           }
      (*pace_ctx).aceimpl->ace->set_basis(*(*pace_ctx).aceimpl->basis_set, 1);
      
      double cutoff=0.;
      for (int i = 0; i < nspecies; i++) {
        for (int j = 0; j < nspecies; j++) {
          *rcut_max = std::max( cutoff, (*pace_ctx).aceimpl->basis_set->radial_functions->cut(i,j) );
        }
      }
      
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(pace_init)
  {
    OperatorNodeFactory::instance()->register_factory( "pace_init" ,make_simple_operator< PaceInit > );
  }

}


