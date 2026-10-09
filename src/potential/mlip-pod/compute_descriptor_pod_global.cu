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
#include <exanb/core/grid_fields.h>
#include <exanb/core/domain.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>

#include <cstdint>
#include <string>
#include <vector>
#include <mpi.h>

#include "pod_params.h"
#include "pod_config.h"

// Global linear-fitting design matrix for POD, built from the output of compute_descriptor_pod
// (pod_descriptors + compute_derivative pda_* aggregate fields), without recomputing them.
//   rows             : 1 + 3*natoms + 6 (natoms = total atom count, all MPI ranks)
//   columns          : nCoeffPerElement * nelements, one block per element
//   row 0            : summed descriptor, including the nl1 one-body atom-count term
//   rows 1..3*natoms : force-signed gradient w.r.t. atom id's x/y/z (row = 1 + 3*id + xyz),
//                      F = +coeff . row
//   last 6 rows      : sum over slots of r . gradient, Voigt order [xx,yy,zz,yz,xz,xy]
// The pda_* aggregate holds one slot per CENTRAL atom species, so the gradient rows loop over every
// species. Each rank accumulates its owned and ghost slots (ghosts weighted by their periodic-image
// position) before the final Allreduce: this operator must run BEFORE update_opt_from_ghost folds
// the ghost aggregates into their owners.
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class ComputeDescriptorPodGlobal : public OperatorNode
  {
    ADD_SLOT( MPI_Comm , mpi  , INPUT , REQUIRED );
    ADD_SLOT( GridT    , grid , INPUT , REQUIRED );
    ADD_SLOT( Domain   , domain , INPUT , REQUIRED );
    ADD_SLOT( PodContext , pod_ctx , INPUT , REQUIRED , DocString{"POD context built by pod_init"} );
    ADD_SLOT( onika::memory::CudaMMVector<double> , pod_descriptors , INPUT , OPTIONAL , DocString{"see compute_descriptor_pod; required"} );
    ADD_SLOT( long     , ncoeff , INPUT , OPTIONAL , DocString{"see compute_descriptor_pod; required (= Mdesc*nClusters)"} );
    ADD_SLOT( std::string , deriv_agg_field_prefix , INPUT , std::string("pda_")
            , DocString{"Must match compute_descriptor_pod's deriv_agg_field_prefix. Must run before update_opt_from_ghost on these fields."} );

    ADD_SLOT( onika::memory::CudaMMVector<double> , pod_global , OUTPUT
            , DocString{"Row-major (1+3*natoms+6) x ncoeff_all global design matrix (natoms = total atom count, all MPI ranks). Row 0 = summed descriptor; rows 1..3*natoms = force-signed gradient (row=1+3*id+xyz); last 6 rows = virial, Voigt order [xx,yy,zz,yz,xz,xy]."} );
    ADD_SLOT( long , ncoeff_all , OUTPUT , DocString{"Number of columns = nCoeffPerElement*nelements"} );

  public:
    inline void execute() override final
    {
      if( ! pod_descriptors.has_value() || ! ncoeff.has_value() )
      {
        fatal_error() << "compute_descriptor_pod_global: pod_descriptors/ncoeff unavailable -- run compute_descriptor_pod first" << std::endl;
      }

      auto & eapod0 = *pod_ctx->m_eapod[0];

      const long Mdesc = eapod0.Mdesc;
      const long nClusters = eapod0.nClusters;
      const long nl1 = eapod0.nl1;
      const long nelements = eapod0.nelements;
      const long ncols = static_cast<long>(eapod0.nCoeffPerElement) * nelements;
      *ncoeff_all = ncols;
      const long nc = *ncoeff;   // == Mdesc*nClusters
      const long nc3 = nc * 3 * nelements;   // widened by nelements, matches compute_descriptor_pod.cu

      std::vector<const double*> agg_ptr( static_cast<size_t>(nc3), nullptr );
      for( long k=0; k<nc3; k++ )
      {
        agg_ptr[k] = grid->flat_array_data_nocreate( field::mk_generic_real( *deriv_agg_field_prefix + std::to_string(k) ) );
        if( agg_ptr[k] == nullptr )
        {
          fatal_error() << "compute_descriptor_pod_global: field '"<<*deriv_agg_field_prefix<<k<<"' not found -- run compute_descriptor_pod with compute_derivative: true first" << std::endl;
        }
      }

      const auto * cell_particle_offset = grid->cell_particle_offset_data();
      const size_t n_cells = grid->number_of_cells();

      // local owned (non-ghost) atom count -> global total via one small Allreduce, so every rank
      // allocates the same full-size array before the main per-atom Allreduce below
      long local_owned = 0;
      for( size_t ci=0; ci<n_cells; ci++ )
      {
        if( grid->is_ghost_cell(ci) ) continue;
        local_owned += static_cast<long>( grid->cell(ci).size() );
      }
      long natoms_global = 0;
      MPI_Allreduce( &local_owned, &natoms_global, 1, MPI_LONG, MPI_SUM, *mpi );

      const long grad_rows = 1 + 3*natoms_global;
      const long virial_row0 = grad_rows;
      const long rows = grad_rows + 6;

      pod_global->clear();
      pod_global->resize( static_cast<size_t>(rows) * static_cast<size_t>(ncols), 0.0 );
      double * const arr = pod_global->data();
      const Mat3d xform = domain->xform();

      for( size_t ci=0; ci<n_cells; ci++ )
      {
        const bool is_ghost = grid->is_ghost_cell(ci);
        const auto & cell = grid->cell(ci);
        const size_t np = cell.size();
        for( size_t pi=0; pi<np; pi++ )
        {
          const size_t p = cell_particle_offset[ci] + pi;
          const uint64_t id = cell[field::id][pi];

          if( ! is_ghost )
          {
            // row 0 (one-body + descriptor sum) only needs this atom's own species
            const long ti0 = pod_ctx->type_map[ cell[field::type][pi] ] - 1;
            if( nl1 > 0 ) arr[ eapod0.nCoeffPerElement*ti0 ] += 1.0; // one-body atom-count term, row 0
            const double * const src = pod_descriptors->data() + static_cast<size_t>(nc) * p;
            for( long m=0; m<Mdesc; m++ )
            for( long j=0; j<nClusters; j++ )
            {
              const long col = eapod0.nCoeffPerElement*ti0 + nl1 + m + j*Mdesc;
              arr[ col ] += src[ m + Mdesc*j ];
            }
          }

          // real-frame position of this slot (owned atom or ghost image), same frame as the
          // pda_* gradients, which were computed on xform-applied pair vectors
          const Vec3d r = xform * Vec3d{ cell[field::rx][pi], cell[field::ry][pi], cell[field::rz][pi] };

          // Gradient rows: the pda_* aggregate has one slot per CENTRAL-atom species, loop over all.
          // Virial rows: accumulated per slot, owned AND ghost, before the id-collapse. A ghost slot
          // holds the neighbor-side terms of the pairs that reached that periodic image, so it is
          // weighted by the image position: collapsing by id first would drop the box-vector term
          // of every boundary-crossing pair.
          const long grad_row0 = 1 + 3*static_cast<long>(id);
          for( long ti0=0; ti0<nelements; ti0++ )
          for( long m=0; m<Mdesc; m++ )
          for( long j=0; j<nClusters; j++ )
          {
            const long col = eapod0.nCoeffPerElement*ti0 + nl1 + m + j*Mdesc;
            const long comp = (m + Mdesc*j + Mdesc*nClusters*ti0)*3;
            const double dx = agg_ptr[comp+0][p];
            const double dy = agg_ptr[comp+1][p];
            const double dz = agg_ptr[comp+2][p];
            arr[ (grad_row0+0)*ncols + col ] += dx;
            arr[ (grad_row0+1)*ncols + col ] += dy;
            arr[ (grad_row0+2)*ncols + col ] += dz;

            arr[ (virial_row0+0)*ncols + col ] += dx*r.x; // xx
            arr[ (virial_row0+1)*ncols + col ] += dy*r.y; // yy
            arr[ (virial_row0+2)*ncols + col ] += dz*r.z; // zz
            arr[ (virial_row0+3)*ncols + col ] += dz*r.y; // yz
            arr[ (virial_row0+4)*ncols + col ] += dz*r.x; // xz
            arr[ (virial_row0+5)*ncols + col ] += dy*r.x; // xy
          }
        }
      }

      // Combine every rank's local partial contribution to all rows (descriptor, gradient, virial).
      MPI_Allreduce( MPI_IN_PLACE, arr, static_cast<int>(rows*ncols), MPI_DOUBLE, MPI_SUM, *mpi );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Global linear-fitting design matrix for POD. Row-major (1+3*natoms+6) x ncoeff_all array,
ncoeff_all = nCoeffPerElement*nelements (one column block per element):

  row 0             -- summed descriptor, including the one-body atom-count term.
                       coeff . row = total energy.
  rows 1..3*natoms  -- force-signed gradient w.r.t. atom m's x/y/z, at row 1+3*m+xyz
                       (m = particle id, 0-indexed). coeff . row = force component.
  rows 3N+1..3N+6   -- virial, Voigt order [xx,yy,zz,yz,xz,xy]. coeff . row = virial component.

The matrix has no reference-label column. Must run after compute_descriptor_pod
(compute_derivative: true) and before any update_opt_from_ghost on its aggregate fields.

Usage example:

init_parameters:
  - species
  - pod_init: { parameters: { pod_file: "Ta_param.pod", coeff_file: "Ta_coefficients.pod" } }

compute_descriptor_pod: { compute_derivative: true }
compute_descriptor_pod_global
write_descriptor_pod_global: { filename: "pod_global.txt" }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(compute_descriptor_pod_global)
  {
    OperatorNodeFactory::instance()->register_factory("compute_descriptor_pod_global", make_grid_variant_operator<ComputeDescriptorPodGlobal>);
  }

}
