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

#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/compute_cell_particles.h>

#include <exaStamp/delaunay/delaunay_tessellation.h>
#include <exaStamp/delaunay/dxa_core_ownership.h>

#include <string>
#include <vector>

// User-requested operator: marks every atom with the id of the dislocation it belongs to (-1 for
// an atom not part of any dislocation's own core). Final step of the core-atom marking pipeline --
// see compute_dxa_core_atoms.cpp's own header comment for the actual per-tet marking mechanism
// (a combinatorial analog of OVITO's own markCoreAtoms) and compute_dxa_mpi_stitch_lines.cpp's own
// header comment on DXALocalToFinalDislocationId for why a per-rank-local id needs remapping here.
//
// Purely a translation step: DXACoreTetOwnership::tet_dislocation_id is already exactly right,
// just in this rank's own LOCAL (pre-MPI-stitch) numbering -- resolve each tet's own local id
// through DXALocalToFinalDislocationId, then assign the FINAL id to that tet's own 4 vertex atoms
// (mesh.vertex_particle_index). An atom shared between two tets whose own final ids happen to
// differ (a genuine boundary atom between two dislocations' own claimed core regions) gets
// whichever tet's assignment is processed last -- not a meaningful ambiguity in practice (both
// ids are equally "correct" for that shared atom), not worth breaking a tie over.
namespace exaStamp
{
  using namespace exanb;

  // Namespace-scope (not a local struct inside execute()) -- same convention
  // dxa_lattice_cluster_fields.cpp's own functor already uses, matching this codebase's own
  // established rule for anything handed to compute_cell_particles (a local/function-scope type
  // caused this exact field to silently never materialize for write_xyz, even though the operator
  // itself ran and computed real values -- found by direct A/B diff against that known-working file).
  struct DXAMarkCoreAtomsFunctor
  {
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    const int32_t * const __restrict__ m_atom_dislocation_id = nullptr;
    inline void operator () ( size_t cell, unsigned int part, double& out ) const
    {
      out = static_cast<double>( m_atom_dislocation_id[ m_cell_particle_offset[cell] + part ] );
    }
  };

  template<class GridT>
  class DXAMarkCoreAtoms : public OperatorNode
  {
    ADD_SLOT( GridT                          , grid                               , INPUT_OUTPUT );
    ADD_SLOT( DelaunayTessellation            , delaunay_tessellation              , INPUT , REQUIRED );
    ADD_SLOT( DXACoreTetOwnership             , dxa_core_tet_ownership             , INPUT , REQUIRED );
    ADD_SLOT( DXALocalToFinalDislocationId    , dxa_local_to_final_dislocation_id  , INPUT , REQUIRED );
    ADD_SLOT( std::string                     , dislocation_id_field               , INPUT , std::string("dxa_disloc_id") , DocString{"Name of the resulting per-particle dislocation-id grid field. IMPORTANT: exanb's own dynamic grid-field name storage (grid_fields.h's own XNB_DECLARE_DYNAMIC_FIELD, NAME_MAX_LEN=16) silently truncates any name over 15 characters with NO error or warning -- a longer name here will create a field whose REAL stored name differs from what you asked for, so a downstream consumer (write_xyz, etc.) requesting the FULL name you specified will never find it. Keep this at 15 characters or fewer."} );
    ADD_SLOT( long                            , n_core_atoms                       , OUTPUT , DocString{"Number of atoms marked with a real (non -1) dislocation id"} );

  public:
    inline void execute () override final
    {
      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const DXACoreTetOwnership& core = *dxa_core_tet_ownership;
      const std::vector<int32_t>& remap = dxa_local_to_final_dislocation_id->final_id;

      std::vector<int32_t> atom_dislocation_id( grid->number_of_particles(), -1 );
      for(size_t t=0; t<core.tet_dislocation_id.size(); t++)
      {
        const int32_t local_id = core.tet_dislocation_id[t];
        if( local_id < 0 ) { continue; }
        if( static_cast<size_t>(local_id) >= remap.size() ) { continue; } // shouldn't happen, guard anyway
        const int32_t final_id = remap[local_id];
        if( final_id < 0 ) { continue; } // this local dislocation didn't survive stitching (e.g. dropped as a redundant duplicate with no valid resolution)
        for( uint32_t v : mesh.tetrahedra[t] )
        {
          atom_dislocation_id[ mesh.vertex_particle_index[v] ] = final_id;
        }
      }

      long n_marked = 0;
      for( int32_t d : atom_dislocation_id ) { if( d >= 0 ) { ++n_marked; } }
      *n_core_atoms = n_marked;
      lout << "dxa_mark_core_atoms: " << n_marked << " / " << atom_dislocation_id.size() << " atoms marked with a real dislocation id" << std::endl;

      auto field_acc = grid->field_accessor( field::mk_generic_real( *dislocation_id_field ) );
      DXAMarkCoreAtomsFunctor func = { grid->cell_particle_offset_data(), atom_dislocation_id.data() };
      compute_cell_particles( *grid, false, func, onika::make_flat_tuple( field_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Marks every atom with the id of the dislocation it belongs to (-1 if not part of any dislocation's
own core) -- final step of the core-atom marking pipeline, translating compute_dxa_core_atoms' own
per-tet, per-rank-LOCAL dislocation ownership into the FINAL (post-MPI-stitch) numbering via
compute_dxa_mpi_stitch_lines' own DXALocalToFinalDislocationId remap, then writing it out as a named
per-particle grid field. Only touches owned particles -- ghost copies of this field are not
synchronized across MPI ranks (same convention as dxa_lattice_cluster_fields).

IMPORTANT, dislocation_id_field's own name length: keep it to 15 characters or fewer -- see that
slot's own DocString for why (a real, silent-truncation gotcha in exanb's own dynamic grid-field
name storage, found the hard way: a longer name creates a field under a truncated name with zero
error, so any consumer requesting the FULL name you asked for never finds it).

Usage example:

compute_dxa_circuit_sweep: {}
compute_dxa_core_atoms: {}
compute_dxa_mpi_stitch_lines: {}
dxa_mark_core_atoms: { dislocation_id_field: dxa_dislocation_id }

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(dxa_mark_core_atoms)
  {
    OperatorNodeFactory::instance()->register_factory( "dxa_mark_core_atoms", make_grid_variant_operator< DXAMarkCoreAtoms > );
  }

}
