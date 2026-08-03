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

#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/compute_cell_particles.h>

#include <exaStamp/delaunay/dxa_lattice_correspondence.h>

#include <string>

// Materializes compute_dxa_lattice_clusters' own flat DXALatticeClusters::atom_cluster buffer into
// a named per-particle grid field (field::mk_generic_real), so any consumer expecting a real grid
// field (write_xyz, write_delaunay_vtk's color_fields, ...) can read it directly -- same bridge
// cna_fields/ptm_fields already provide for their own flat output buffers (src/cna/compute_cna.cu,
// src/ptm/compute_ptm.cu), previously absent for DXALatticeClusters. Purely pointwise, no neighbor
// search -- same exanb::compute_cell_particles pattern as those two.
namespace exaStamp
{
  using namespace exanb;

  struct DXALatticeClusterFieldsFunctor
  {
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    const int32_t * const __restrict__ m_atom_cluster = nullptr;

    inline void operator () ( size_t cell, unsigned int part, double& cluster_id_out ) const
    {
      cluster_id_out = static_cast<double>( m_atom_cluster[ m_cell_particle_offset[cell] + part ] );
    }
  };

  template<class GridT>
  class DXALatticeClusterFields : public OperatorNode
  {
    ADD_SLOT( GridT                , grid                 , INPUT_OUTPUT );
    ADD_SLOT( DXALatticeClusters   , dxa_lattice_clusters , INPUT , REQUIRED );
    ADD_SLOT( std::string          , cluster_field        , INPUT , std::string("dxa_cluster_id") , DocString{"Name of the resulting per-particle cluster-id grid field"} );

  public:
    inline void execute () override final
    {
      auto cluster_acc = grid->field_accessor( field::mk_generic_real( *cluster_field ) );
      DXALatticeClusterFieldsFunctor func = { grid->cell_particle_offset_data(), dxa_lattice_clusters->atom_cluster.data() };
      compute_cell_particles( *grid, false, func, onika::make_flat_tuple( cluster_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Copies compute_dxa_lattice_clusters' flat per-particle output buffer (DXALatticeClusters::
atom_cluster, 0 = unresolved) into a named per-particle grid field, so any consumer expecting a real
grid field (write_xyz, write_delaunay_vtk's color_fields, ...) can read it directly. Purely pointwise
(no neighbor search) -- run any time after compute_dxa_lattice_clusters. Only touches owned
particles -- ghost copies of this field are not synchronized across MPI ranks.

Usage example:

compute_dxa_lattice_correspondence: { rcut: 5.0 ang }
compute_dxa_lattice_clusters: {}
dxa_lattice_cluster_fields: { cluster_field: dxa_cluster_id }

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(dxa_lattice_cluster_fields)
  {
    OperatorNodeFactory::instance()->register_factory( "dxa_lattice_cluster_fields", make_grid_variant_operator< DXALatticeClusterFields > );
  }

}
