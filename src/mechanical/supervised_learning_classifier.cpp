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

#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/grid.h>
#include <onika/math/basic_types_stream.h>
#include <onika/log.h>
#include <exanb/core/domain.h>
#include <onika/physics/units.h>
#include <onika/string_utils.h>
#include <exaStamp/compute/thermodynamic_state.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/grid_fields.h>

#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/defbox/deformation.h>
#include <exanb/defbox/deformation_stream.h>
#include <exanb/defbox/deformation_yaml.h>
#include <exanb/defbox/deformation_math.h>

#include <exaStamp/mechanical/cell_particles_local_structural_metrics.h>

#include <exaStamp/potential/snaplegacy/SnapLegacyBS.h>

#include <onika/soatl/packed_field_arrays.h>
#include <memory>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <mpi.h>
#include <string>
#include <iomanip>
#include <experimental/filesystem>
#include <math.h>

namespace exaStamp
{
  using namespace exanb;
  
  template< class GridT
           >
  class SupervisedLearningClassifierOperator : public OperatorNode
  {

    using VariablesVec = std::vector<std::string>;
    using LDAVec = std::vector<double>;
    using VecVec = std::vector<Vec3d>;
        
    ADD_SLOT( MPI_Comm           , mpi                 , INPUT );
    ADD_SLOT( GridT              , grid                , INPUT );
    ADD_SLOT( Domain             , domain              , INPUT );
    ADD_SLOT( Deformation        , defbox              , INPUT );    
    ADD_SLOT( long               , timestep            , INPUT);
    ADD_SLOT( ParticleSpecies    , species             , INPUT , REQUIRED );

    ADD_SLOT( GridParticleLocalStructuralMetrics, local_structural_data , INPUT_OUTPUT, OPTIONAL);        
    ADD_SLOT( VariablesVec                      , per_atom_data         , INPUT, OPTIONAL);

    ADD_SLOT( VecVec             , lda_scalings , INPUT  , REQUIRED );
    ADD_SLOT( LDAVec             , overall_mean , INPUT  , REQUIRED );    
    ADD_SLOT( VecVec             , decision     , INPUT  , REQUIRED );        
    ADD_SLOT( LDAVec             , biais        , INPUT  , REQUIRED );        
    ADD_SLOT( Vec3d              , mean_bcc     , INPUT  , REQUIRED );
    ADD_SLOT( Vec3d              , mean_fcc     , INPUT  , REQUIRED );
    ADD_SLOT( Vec3d              , mean_hcp     , INPUT  , REQUIRED );
    ADD_SLOT( Vec3d              , mean_sc      , INPUT  , REQUIRED );            
    ADD_SLOT( Mat3d              , cova_bcc     , INPUT  , REQUIRED );
    ADD_SLOT( Mat3d              , cova_fcc     , INPUT  , REQUIRED );
    ADD_SLOT( Mat3d              , cova_hcp     , INPUT  , REQUIRED );
    ADD_SLOT( Mat3d              , cova_sc      , INPUT  , REQUIRED );                
    ADD_SLOT( double             , distance     , INPUT  , REQUIRED );                
    
  public:
  
    inline void execute () override final
    {

      static constexpr bool is_ghosts = false;
      
      GridT& grid = *(this->grid);
      size_t n_cells = grid.number_of_cells();
      auto cells = grid.cells();

      GridParticleLocalStructuralMetrics& local_structural_data = *(this->local_structural_data);
      VariablesVec per_atom_variables;

      int jmax=3;
      int ncompbs = SnapLegacyBS::n_idx_bs( 2*jmax );		      

      for(size_t c=0; c<n_cells;++c)
      {
        int np = 0;
        if( !grid.is_ghost_cell(c) || is_ghosts )
        {
          np = cells[c].size();
          local_structural_data[c].crystal_structure.resize(np);
          for(int pos=0;pos<np;++pos)
          {

	    // First part : perform LDA reduction dimension
	    Vec3d reddata = {0.,0.,0.};

	    for (int i=0;i<ncompbs;i++) {
	      reddata.x += (*lda_scalings)[i].x * (local_structural_data[c].bispectrum[pos][i] - (*overall_mean)[i]);
	      reddata.y += (*lda_scalings)[i].y * (local_structural_data[c].bispectrum[pos][i] - (*overall_mean)[i]);
	      reddata.z += (*lda_scalings)[i].z * (local_structural_data[c].bispectrum[pos][i] - (*overall_mean)[i]);
	    }

	    vector<double> scores = {0.,0.,0.,0.};
	    for (int i=0;i<4;i++) {
	      scores[i] = (*biais)[i] + (*decision)[i].x * reddata.x + (*decision)[i].y * reddata.y + (*decision)[i].z * reddata.z;
	    }

	    vector<double> expscores = {0.,0.,0.,0.};
	    double sumexpscores = 0.;
	    for (int i=0;i<4;i++) {
	      expscores[i] = exp(scores[i]);
	      sumexpscores += expscores[i];
	    }
	    
	    vector<double> coefs = {0.,0.,0.,0.};
	    for (int i=0;i<4;i++) {
	      coefs[i] = expscores[i]/sumexpscores;
	    }

	    vector<int> types = {0,1,2,3,4};	    
	    int crystal_structure = types[4];
	    if (coefs[0]>0.5) crystal_structure = types[0];
	    if (coefs[1]>0.5) crystal_structure = types[1];
	    if (coefs[2]>0.5) crystal_structure = types[2];
	    if (coefs[3]>0.5) crystal_structure = types[3];		

	    // Second part : compute mahalanobis distance
	    double maha_bcc = 0.;
	    double maha_fcc = 0.;
	    double maha_hcp = 0.;
	    double maha_sc = 0.;

	    if (crystal_structure == 0) {
	      Vec3d xmu_bcc = reddata-(*mean_bcc);
	      Mat3d inv_cova_bcc = inverse(*cova_bcc);
	      Vec3d right_bcc = inv_cova_bcc * xmu_bcc;		
	      maha_bcc = sqrt(dot(xmu_bcc, right_bcc));
	      if (maha_bcc > *distance) crystal_structure=4;
	    } else if (crystal_structure == 1) {
	      Vec3d xmu_fcc = reddata-(*mean_fcc);
	      Mat3d inv_cova_fcc = inverse(*cova_fcc);
	      Vec3d right_fcc = inv_cova_fcc * xmu_fcc;		
	      maha_fcc = sqrt(dot(xmu_fcc, right_fcc));
	      if (maha_fcc > *distance) crystal_structure=4;		  
	    } else if (crystal_structure == 2) {
	      Vec3d xmu_hcp = reddata-(*mean_hcp);
	      Mat3d inv_cova_hcp = inverse(*cova_hcp);
	      Vec3d right_hcp = inv_cova_hcp * xmu_hcp;		
	      maha_hcp = sqrt(dot(xmu_hcp, right_hcp));
	      if (maha_hcp > *distance) crystal_structure=4;		  
	    } else if (crystal_structure == 3) {
	      Vec3d xmu_sc = reddata-(*mean_sc);
	      Mat3d inv_cova_sc = inverse(*cova_sc);
	      Vec3d right_sc = inv_cova_sc * xmu_sc;		
	      maha_sc = sqrt(dot(xmu_sc, right_sc));
	      if (maha_sc > *distance) crystal_structure=4;		  
	    }	    
	    local_structural_data[c].crystal_structure[pos] = crystal_structure;
	  }
	}
      }
      
      
    }
    
};

  template<class GridT> using SupervisedLearningClassifierOperatorTmpl = SupervisedLearningClassifierOperator<GridT>;
  
  // === register factories ===  
  ONIKA_AUTORUN_INIT(supervised_learning_classifier)
  {
    OperatorNodeFactory::instance()->register_factory( "supervised_learning_classifier", make_grid_variant_operator< SupervisedLearningClassifierOperatorTmpl > );
  }

}
