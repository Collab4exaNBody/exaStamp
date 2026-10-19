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
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/mpi/grid_update_ghosts.h>
#include <exanb/mpi/update_ghost_functors.h>
#include <exaStamp/molecule/id_map.h>
#include <exaStamp/molecule/molecule_species.h>
#include <mpi.h>
#include <limits>

namespace exaStamp
{

  template<typename GridT
    , class = AssertGridHasFields< GridT, field::_id, field::_idmol, field::_cmol >
    >
  class MoleculeCheckNumbering : public OperatorNode
  {
    ADD_SLOT( MPI_Comm                 , mpi                 , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( Domain                   , domain              , INPUT , REQUIRED );
    ADD_SLOT( GridT                    , grid                , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( IdMap                    , id_map              , INPUT , REQUIRED );
    ADD_SLOT( IdMapGhosts              , id_map_ghosts       , INPUT , REQUIRED );
    ADD_SLOT( UpdateGhostConfig        , update_ghost_config , INPUT_OUTPUT , UpdateGhostConfig{} );
    ADD_SLOT( GhostCommunicationScheme , ghost_comm_scheme   , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( MoleculeSpeciesVector    , molecules           , INPUT , REQUIRED , DocString{"Molecule descriptions and type of each molecule"} );
    ADD_SLOT( long                     , expected_molecules  , INPUT , -1 , DocString{"if >= 0, expected number of molecules"} );

  public:
    // a check has no consumer : without this, the operator would be removed from the graph
    inline bool is_sink() const override final { return true; }

    inline std::string documentation() const override final
    {
      return R"EOF(
Checks the molecule numbering (idmol field) and stops with an error if it is wrong :
every idmol is in [0, number of molecules), atoms connected by a bond have the same idmol (also across MPI
processes), every molecule has at least one atom, molecules are numbered in increasing order of their smallest
atom id, the type of every molecule is a valid molecule species and, if given, the number of molecules is
expected_molecules. Requires id_map and up to date ghost particles.
)EOF";
    }

    inline void execute () override final
    {
      auto cells = grid->cells_accessor();
      const size_t n_cells = grid->number_of_cells();
      const auto idmol_field = grid->field_accessor( field::idmol );
      const auto cmol_field = grid->field_accessor( field::cmol );
      static constexpr uint64_t null_loc = std::numeric_limits<uint64_t>::max();
      static constexpr uint64_t no_id = std::numeric_limits<uint64_t>::max();
      static constexpr size_t null_index = std::numeric_limits<size_t>::max();
      static constexpr GridCellValues * null_cell_values = nullptr;

      const size_t n_mol = molecules->number_of_molecule_instances();
      if( n_mol == 0 ) fatal_error() << "molecule_check_numbering: no molecule numbering (empty molecule type table)" << std::endl;

      // ghost copies need the idmol of their owner
      auto pecfunc = [self=this](auto ... args) { return self->parallel_execution_context(args ...); };
      auto peqfunc = [self=this]() -> onika::parallel::ParallelExecutionQueue& { return self->parallel_execution_queue(); };
      const auto ghost_update_fields = onika::make_flat_tuple( idmol_field );
      auto ghost_comm_buffers = UpdateGhostsUtils::make_forward_ghost_update_manager( cells, ghost_update_fields );
      grid_update_ghosts( ldbg, *mpi, *ghost_comm_scheme, grid.get_pointer(), *domain, null_cell_values,
                          ghost_comm_buffers, pecfunc, peqfunc, ghost_update_fields, *update_ghost_config );

      unsigned long n_out_of_range = 0, n_bond_mismatch = 0, n_bonds = 0, n_bad_type = 0, n_unordered = 0, n_empty = 0;
      std::vector<uint64_t> min_atom_id( n_mol , no_id );
      for(size_t c=0;c<n_cells;c++) if( ! grid->is_ghost_cell(c) )
      {
        const auto * __restrict__ ids = cells[c][field::id];
        const auto * __restrict__ cmol = cells[c][cmol_field];
        const auto * __restrict__ mol_ids = cells[c][idmol_field];
        const size_t n = cells[c].size();
        for(size_t p=0;p<n;p++)
        {
          const uint64_t idmol = mol_ids[p];
          if( idmol >= n_mol ) { ++ n_out_of_range; continue; }
          min_atom_id[idmol] = std::min( min_atom_id[idmol] , uint64_t(ids[p]) );
          for(size_t j=0;j<cmol[p].size();j++)
          {
            const auto conloc = atom_from_idmap_if_found( cmol[p][j] , *id_map , *id_map_ghosts , null_loc );
            if( conloc == null_loc ) continue;
            size_t cell=null_index, pos=null_index;
            decode_cell_particle( conloc , cell, pos );
            ++ n_bonds;
            if( cells[cell][idmol_field][pos] != idmol ) ++ n_bond_mismatch;
          }
        }
      }
      MPI_Allreduce( MPI_IN_PLACE, min_atom_id.data(), n_mol, MPI_UINT64_T, MPI_MIN, *mpi );
      for(size_t m=0;m<n_mol;m++)
      {
        if( min_atom_id[m] == no_id ) ++ n_empty;
        else if( m > 0 && min_atom_id[m-1] != no_id && min_atom_id[m] <= min_atom_id[m-1] ) ++ n_unordered;
        if( molecules->m_molecule_type[m] >= molecules->m_molecules.size() ) ++ n_bad_type;
      }
      unsigned long counts[3] = { n_out_of_range , n_bond_mismatch , n_bonds };
      MPI_Allreduce( MPI_IN_PLACE, counts, 3, MPI_UNSIGNED_LONG, MPI_SUM, *mpi );
      n_out_of_range = counts[0]; n_bond_mismatch = counts[1]; n_bonds = counts[2];

      lout << "molecule_check_numbering: "<<n_mol<<" molecules, "<<n_bonds<<" bonds checked"<<std::endl;
      bool ok = true;
      if( n_out_of_range > 0 ) { lout << "ERROR: "<<n_out_of_range<<" atoms with idmol out of range"<<std::endl; ok = false; }
      if( n_bond_mismatch > 0 ) { lout << "ERROR: "<<n_bond_mismatch<<" bonds between atoms with different idmol"<<std::endl; ok = false; }
      if( n_empty > 0 ) { lout << "ERROR: "<<n_empty<<" molecules without atoms"<<std::endl; ok = false; }
      if( n_unordered > 0 ) { lout << "ERROR: "<<n_unordered<<" molecules not in increasing order of smallest atom id"<<std::endl; ok = false; }
      if( n_bad_type > 0 ) { lout << "ERROR: "<<n_bad_type<<" molecules with an invalid type"<<std::endl; ok = false; }
      if( *expected_molecules >= 0 && n_mol != size_t(*expected_molecules) )
      {
        lout << "ERROR: "<<n_mol<<" molecules, expected "<<*expected_molecules<<std::endl;
        ok = false;
      }
      if( ! ok ) fatal_error() << "molecule_check_numbering: wrong molecule numbering" << std::endl;
      lout << "molecule_check_numbering: OK" << std::endl;
    }
  };

  template<class GridT> using MoleculeCheckNumberingTmpl = MoleculeCheckNumbering<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(molecule_check_numbering)
  {
   OperatorNodeFactory::instance()->register_factory( "molecule_check_numbering", make_grid_variant_operator<MoleculeCheckNumberingTmpl> );
  }

}
