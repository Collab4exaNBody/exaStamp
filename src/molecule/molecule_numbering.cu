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
#include <exaStamp/molecule/id_map.h>
#include <exaStamp/molecule/molecule_species.h>
#include <exaStamp/molecule/molecule_min_atom_id.h>
#include <exaStamp/molecule/molecule_numbering.h>
#include <mpi.h>

namespace exaStamp
{

  template<typename GridT
    , class = AssertGridHasFields< GridT, field::_id, field::_idmol, field::_cmol >
    >
  class MoleculeNumbering : public OperatorNode
  {
    ADD_SLOT( MPI_Comm                 , mpi                 , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( Domain                   , domain              , INPUT , REQUIRED );
    ADD_SLOT( GridT                    , grid                , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( IdMap                    , id_map              , INPUT , REQUIRED );
    ADD_SLOT( IdMapGhosts              , id_map_ghosts       , INPUT , REQUIRED );
    ADD_SLOT( UpdateGhostConfig        , update_ghost_config , INPUT_OUTPUT , UpdateGhostConfig{} );
    ADD_SLOT( GhostCommunicationScheme , ghost_comm_scheme   , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( MoleculeSpeciesVector    , molecules           , INPUT_OUTPUT , REQUIRED , DocString{"Molecule descriptions, receives the type of each molecule"} );

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Numbers molecules : the idmol field of each particle is set to the index of its molecule, 0 <= idmol < number of
molecules, molecules being numbered in increasing order of their smallest atom id. Molecules are the groups of atoms
connected by bonds (cmol field), across MPI processes. The type of each molecule is kept (from the molecule type
table and the previous idmol values, type 0 if unknown).
Use it after replicate_domain_cmol (replicas share the idmol of the original molecules) or after reading a dump
written before molecule ids were dense. Requires id_map and up to date ghost particles.
)EOF";
    }

    inline void execute () override final
    {
      auto cells = grid->cells_accessor();
      const size_t n_cells = grid->number_of_cells();
      const auto idmol_field = grid->field_accessor( field::idmol );

      // previous molecule type of each owned particle
      const auto & old_types = molecules->m_molecule_type;
      std::vector<size_t> cell_offset( n_cells , 0 );
      std::vector<uint8_t> particle_type;
      size_t n_unknown_type = 0;
      for(size_t c=0;c<n_cells;c++) if( ! grid->is_ghost_cell(c) )
      {
        cell_offset[c] = particle_type.size();
        const size_t n = cells[c].size();
        for(size_t p=0;p<n;p++)
        {
          const uint64_t idmol = cells[c][idmol_field][p];
          uint8_t t = 0;
          if( ! old_types.empty() )
          {
            if( idmol < old_types.size() ) t = old_types[idmol];
            else ++ n_unknown_type;
          }
          particle_type.push_back( t );
        }
      }
      MPI_Allreduce( MPI_IN_PLACE, &n_unknown_type, 1, MPI_UNSIGNED_LONG, MPI_SUM, *mpi );
      if( n_unknown_type > 0 )
      {
        lout << "WARNING: molecule_numbering: "<<n_unknown_type<<" atoms with unknown molecule type, type 0 assumed"<<std::endl;
      }

      auto pecfunc = [self=this](auto ... args) { return self->parallel_execution_context(args ...); };
      auto peqfunc = [self=this]() -> onika::parallel::ParallelExecutionQueue& { return self->parallel_execution_queue(); };
      propagate_molecule_min_atom_id( ldbg, *mpi, *ghost_comm_scheme, grid.get_pointer(), *domain, *update_ghost_config,
                                      pecfunc, peqfunc, *id_map, *id_map_ghosts );

      const uint8_t * ptype = particle_type.data();
      const size_t * poffset = cell_offset.data();
      molecules->m_molecule_type = densify_molecule_ids( *mpi , *grid , [cells,idmol_field,ptype,poffset](size_t c, size_t p) -> std::pair<uint64_t,uint64_t>
        {
          return { cells[c][idmol_field][p] , ptype[ poffset[c] + p ] };
        } );
      ldbg << "molecule_numbering: "<<molecules->m_molecule_type.size()<<" molecules"<<std::endl;
    }
  };

  template<class GridT> using MoleculeNumberingTmpl = MoleculeNumbering<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(molecule_numbering)
  {
   OperatorNodeFactory::instance()->register_factory( "molecule_numbering", make_grid_variant_operator<MoleculeNumberingTmpl> );
  }

}
