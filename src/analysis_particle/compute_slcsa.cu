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

#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_yaml.h>

#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/compute_cell_particles.h>

#include <cmath>

// GPU-compatible supervised learning crystal structure analysis (SLCSA): classifies each
// particle as BCC/FCC/HCP/SC/other from an already-computed per-particle bispectrum (see
// compute_bispectrum) via a pretrained LDA projection + softmax + Mahalanobis-distance
// rejection. Rewrite of supervised_learning_classifier.cpp: pointwise (no neighbor list),
// reads the bispectrum flat buffer directly (via cell_particle_offset, same convention as
// compute_bispectrum/snap_check_bispectrum) instead of the old ragged
// GridParticleLocalStructuralMetrics storage, and inverts the 4 covariance matrices once per
// call instead of once per particle.
namespace exaStamp
{
  using namespace exanb;

  struct SLCSAFunctor
  {
    const size_t * const __restrict__ cell_particle_offset = nullptr;
    const double * const __restrict__ bispectrum = nullptr;   // [ncoeff*(cell_particle_offset[cell]+part)+i]
    const long ncoeff = 0;
    const Vec3d * const __restrict__ lda_scalings = nullptr;  // [ncoeff], LDA projection matrix rows
    const double * const __restrict__ overall_mean = nullptr; // [ncoeff], per-component centering

    // fixed 4-class (BCC,FCC,HCP,SC) model, embedded by value -- small enough to just copy
    const Vec3d decision[4] = {};
    const double biais[4] = {};
    const Vec3d mean[4] = {};
    const Mat3d inv_cova[4] = {};
    const double max_distance = 0.0;

    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( size_t cell, unsigned int part, double& crystal_structure ) const
    {
      const long off = ncoeff * static_cast<long>( cell_particle_offset[cell] + part );

      Vec3d reddata = { 0., 0., 0. };
      for( long i=0; i<ncoeff; i++ )
      {
        reddata += lda_scalings[i] * ( bispectrum[off+i] - overall_mean[i] );
      }

      double scores[4];
      for( int k=0; k<4; k++ ) { scores[k] = biais[k] + dot( decision[k], reddata ); }

      double sumexp = 0.0;
      double expscores[4];
      for( int k=0; k<4; k++ ) { expscores[k] = exp( scores[k] ); sumexp += expscores[k]; }

      int cs = 4; // default: unclassified
      for( int k=0; k<4; k++ ) { if( ( expscores[k] / sumexp ) > 0.5 ) { cs = k; } }

      if( cs < 4 )
      {
        const Vec3d xmu = reddata - mean[cs];
        const double maha = sqrt( dot( xmu, inv_cova[cs] * xmu ) );
        if( maha > max_distance ) { cs = 4; }
      }

      crystal_structure = static_cast<double>( cs );
    }
  };
}

namespace exanb
{
  template<> struct ComputeCellParticlesTraits<exaStamp::SLCSAFunctor> { static inline constexpr bool CudaCompatible = true; };
}

namespace exaStamp
{
  template<class GridT>
  class ComputeSLCSA : public OperatorNode
  {
    ADD_SLOT( GridT   , grid          , INPUT_OUTPUT );
    ADD_SLOT( onika::memory::CudaMMVector<double> , bispectrum , INPUT , REQUIRED , DocString{"Flat per-particle bispectrum buffer (see compute_bispectrum)"} );
    ADD_SLOT( long     , ncoeff       , INPUT , REQUIRED , DocString{"Number of bispectrum components per particle (see compute_bispectrum)"} );

    ADD_SLOT( onika::memory::CudaMMVector<Vec3d>  , lda_scalings , INPUT , REQUIRED , DocString{"LDA projection matrix, ncoeff rows of 3 (one Vec3d per bispectrum component)"} );
    ADD_SLOT( onika::memory::CudaMMVector<double> , overall_mean , INPUT , REQUIRED , DocString{"Per-component bispectrum mean used to center data before LDA projection, ncoeff values"} );
    ADD_SLOT( std::vector<Vec3d> , decision , INPUT , REQUIRED , DocString{"Softmax decision vectors in LDA space, one Vec3d per class [BCC,FCC,HCP,SC]"} );
    ADD_SLOT( std::vector<double>, biais    , INPUT , REQUIRED , DocString{"Softmax bias, one value per class [BCC,FCC,HCP,SC]"} );
    ADD_SLOT( Vec3d   , mean_bcc , INPUT , REQUIRED , DocString{"BCC reference mean in LDA space"} );
    ADD_SLOT( Vec3d   , mean_fcc , INPUT , REQUIRED , DocString{"FCC reference mean in LDA space"} );
    ADD_SLOT( Vec3d   , mean_hcp , INPUT , REQUIRED , DocString{"HCP reference mean in LDA space"} );
    ADD_SLOT( Vec3d   , mean_sc  , INPUT , REQUIRED , DocString{"SC reference mean in LDA space"} );
    ADD_SLOT( Mat3d   , cova_bcc , INPUT , REQUIRED , DocString{"BCC reference covariance in LDA space"} );
    ADD_SLOT( Mat3d   , cova_fcc , INPUT , REQUIRED , DocString{"FCC reference covariance in LDA space"} );
    ADD_SLOT( Mat3d   , cova_hcp , INPUT , REQUIRED , DocString{"HCP reference covariance in LDA space"} );
    ADD_SLOT( Mat3d   , cova_sc  , INPUT , REQUIRED , DocString{"SC reference covariance in LDA space"} );
    ADD_SLOT( double  , distance , INPUT , REQUIRED , DocString{"Mahalanobis distance threshold beyond which a particle is rejected back to 'other' (class 4)"} );

    ADD_SLOT( std::string , crystal_structure_field , INPUT , std::string("crystal_structure") , DocString{"Name of the resulting classification scalar field (0=BCC,1=FCC,2=HCP,3=SC,4=other)"} );

  public:
    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 ) return;

      if( lda_scalings->size() != static_cast<size_t>(*ncoeff) || overall_mean->size() != static_cast<size_t>(*ncoeff) )
      {
        fatal_error() << "compute_slcsa: lda_scalings/overall_mean size ("<<lda_scalings->size()<<"/"<<overall_mean->size()<<") does not match ncoeff ("<<*ncoeff<<")" << std::endl;
      }
      if( decision->size() != 4 || biais->size() != 4 )
      {
        fatal_error() << "compute_slcsa: decision/biais must each have exactly 4 entries (BCC,FCC,HCP,SC)" << std::endl;
      }

      // covariances inverted once here, not per particle (the old code inverted them inside
      // the per-particle loop -- a real perf bug this rewrite fixes along the way)
      SLCSAFunctor func {
        grid->cell_particle_offset_data(), bispectrum->data(), *ncoeff, lda_scalings->data(), overall_mean->data(),
        { decision->at(0), decision->at(1), decision->at(2), decision->at(3) },
        { biais->at(0), biais->at(1), biais->at(2), biais->at(3) },
        { *mean_bcc, *mean_fcc, *mean_hcp, *mean_sc },
        { inverse(*cova_bcc), inverse(*cova_fcc), inverse(*cova_hcp), inverse(*cova_sc) },
        *distance
      };

      auto cs_acc = grid->field_accessor( field::mk_generic_real( *crystal_structure_field ) );
      compute_cell_particles( *grid, false, func, onika::make_flat_tuple( cs_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Supervised learning crystal structure analysis (SLCSA): classifies each particle as
BCC(0)/FCC(1)/HCP(2)/SC(3)/other(4) from its per-particle SNAP bispectrum (see
compute_bispectrum), via a pretrained LDA dimensionality reduction to 3D, a 4-class
softmax decision, and a final Mahalanobis-distance rejection against each class's
reference distribution. Pointwise (no neighbor list): reads the bispectrum flat buffer
directly via cell_particle_offset, writes a single scalar field. GPU-compatible.

Usage example:

compute_bispectrum:
  parameters: { param: "W.snapparam", coef: "W.snapcoeff" }
compute_slcsa:
  lda_scalings: [ ... ]   # ncoeff Vec3d rows
  overall_mean: [ ... ]   # ncoeff values
  decision: [ ... ]       # 4 Vec3d
  biais: [ ... ]          # 4 values
  mean_bcc: [ ... ]
  mean_fcc: [ ... ]
  mean_hcp: [ ... ]
  mean_sc:  [ ... ]
  cova_bcc: [ ... ]
  cova_fcc: [ ... ]
  cova_hcp: [ ... ]
  cova_sc:  [ ... ]
  distance: 3.0

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_slcsa)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_slcsa", make_grid_variant_operator< ComputeSLCSA > );
  }

}
