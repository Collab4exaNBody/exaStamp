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

#pragma once

#include <vector>
#include <array>
#include <cstdint>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <algorithm>
#include <map>
#include <string>
#include <onika/log.h>
#include <exaStamp/fields.h>
#include <mpi.h>

namespace exaStamp
{

  // Dense molecule numbering of a set of atoms held by a single process (readers) :
  // molecules are the connected components of the bond graph (cmol), numbered 0 .. n-1 in increasing order of their
  // smallest atom id. Bonds to atoms that are not in the set are ignored.
  // returns the number of molecules, mol_index[i] is the molecule index of atom i
  inline size_t number_molecules_by_connectivity( const std::vector<uint64_t>& ids,
                                                  const std::vector< std::array<uint64_t,4> >& cmol,
                                                  std::vector<uint64_t>& mol_index )
  {
    const size_t n = ids.size();
    std::unordered_map<uint64_t,size_t> id_to_index;
    id_to_index.reserve( n );
    for(size_t i=0;i<n;i++) id_to_index[ ids[i] ] = i;

    std::vector<size_t> parent( n );
    std::iota( parent.begin(), parent.end(), size_t(0) );
    auto find_root = [&parent]( size_t i ) -> size_t
    {
      while( parent[i] != i ) { parent[i] = parent[parent[i]]; i = parent[i]; }
      return i;
    };
    for(size_t i=0;i<n;i++)
    {
      for(const auto c : cmol[i])
      {
        if( c == std::numeric_limits<uint64_t>::max() ) continue;
        const auto it = id_to_index.find( c );
        if( it == id_to_index.end() ) continue;
        const size_t a = find_root(i), b = find_root(it->second);
        if( a != b ) parent[ std::max(a,b) ] = std::min(a,b);
      }
    }

    // smallest atom id of each component, then rank components by it
    std::vector<uint64_t> min_id( n , std::numeric_limits<uint64_t>::max() );
    for(size_t i=0;i<n;i++) { const size_t r = find_root(i); min_id[r] = std::min( min_id[r] , ids[i] ); }
    std::vector<size_t> roots;
    for(size_t i=0;i<n;i++) if( find_root(i) == i ) roots.push_back( i );
    std::sort( roots.begin(), roots.end(), [&min_id](size_t a, size_t b){ return min_id[a] < min_id[b]; } );
    std::vector<uint64_t> root_rank( n , 0 );
    for(size_t k=0;k<roots.size();k++) root_rank[ roots[k] ] = k;

    mol_index.assign( n , 0 );
    for(size_t i=0;i<n;i++) mol_index[i] = root_rank[ find_root(i) ];
    return roots.size();
  }


  // Readers : sets the idmol field of atom tuples (which have id and cmol fields) to dense molecule indices and returns
  // the molecule type of each molecule, as the index of its molecule name in the name-ordered molecule_name_map
  // (the order in which readers build MoleculeSpeciesVector::m_molecules).
  template<class AtomTupleVectorT>
  inline std::vector<uint8_t> number_molecules_of_atoms( AtomTupleVectorT& atom_data,
                                                         const std::vector<std::string>& atom_molecule_name,
                                                         const std::map<std::string,int>& molecule_name_map )
  {
    using namespace exanb;
    const size_t n_atoms = atom_data.size();
    std::vector<uint64_t> ids( n_atoms );
    std::vector< std::array<uint64_t,4> > con( n_atoms );
    for(size_t i=0;i<n_atoms;i++)
    {
      ids[i] = atom_data[i][field::id];
      const auto c = atom_data[i][field::cmol];
      con[i] = { c[0] , c[1] , c[2] , c[3] };
    }
    std::vector<uint64_t> mol_index;
    const size_t n_molecules = number_molecules_by_connectivity( ids , con , mol_index );
    std::vector<int> mol_type( n_molecules , -1 );
    for(size_t i=0;i<n_atoms;i++)
    {
      const auto it = molecule_name_map.find( atom_molecule_name[i] );
      if( it == molecule_name_map.end() ) fatal_error() << "unknown molecule name "<<atom_molecule_name[i]<<std::endl;
      const int t = std::distance( molecule_name_map.begin() , it );
      if( mol_type[ mol_index[i] ] != -1 && mol_type[ mol_index[i] ] != t )
      {
        fatal_error() << "atom #"<<ids[i]<<" ("<<atom_molecule_name[i]<<") is bonded to atoms of another molecule type"<<std::endl;
      }
      mol_type[ mol_index[i] ] = t;
      atom_data[i][field::idmol] = mol_index[i];
    }
    return std::vector<uint8_t>( mol_type.begin() , mol_type.end() );
  }


  // Distributed dense numbering : key_type(cell,p) gives, for each owned particle, a key identifying its molecule
  // (unique per molecule over all processes) and the molecule type. Molecules are numbered 0 .. n-1 in increasing key
  // order, the idmol field is set to this index and the type of each molecule is returned (same on all processes).
  template<class GridT, class KeyTypeFuncT>
  inline std::vector<uint8_t> densify_molecule_ids( MPI_Comm comm, GridT& grid, KeyTypeFuncT key_type )
  {
    using namespace exanb;
    auto cells = grid.cells_accessor();
    const auto idmol_field = grid.field_accessor( field::idmol );
    const size_t n_cells = grid.number_of_cells();

    // distinct (key,type) pairs held by this process
    std::map<uint64_t,uint64_t> local;
    for(size_t c=0;c<n_cells;c++) if( ! grid.is_ghost_cell(c) )
    {
      const size_t n = cells[c].size();
      for(size_t p=0;p<n;p++)
      {
        const auto [key,type] = key_type(c,p);
        const auto it = local.find(key);
        if( it != local.end() && it->second != type ) fatal_error() << "molecule key "<<key<<" has atoms of different molecule types"<<std::endl;
        local[key] = type;
      }
    }
    std::vector<uint64_t> send; send.reserve( 2*local.size() );
    for(const auto& kt : local) { send.push_back(kt.first); send.push_back(kt.second); }

    int nprocs=1; MPI_Comm_size(comm,&nprocs);
    int nsend = send.size();
    std::vector<int> counts(nprocs), displs(nprocs,0);
    MPI_Allgather( &nsend, 1, MPI_INT, counts.data(), 1, MPI_INT, comm );
    for(int i=1;i<nprocs;i++) displs[i] = displs[i-1] + counts[i-1];
    std::vector<uint64_t> all( displs[nprocs-1] + counts[nprocs-1] );
    MPI_Allgatherv( send.data(), nsend, MPI_UINT64_T, all.data(), counts.data(), displs.data(), MPI_UINT64_T, comm );

    std::map<uint64_t,uint64_t> global;
    for(size_t i=0;i<all.size();i+=2)
    {
      const auto it = global.find(all[i]);
      if( it != global.end() && it->second != all[i+1] ) fatal_error() << "molecule key "<<all[i]<<" has atoms of different molecule types"<<std::endl;
      global[ all[i] ] = all[i+1];
    }
    std::vector<uint64_t> keys; keys.reserve( global.size() );
    std::vector<uint8_t> types; types.reserve( global.size() );
    for(const auto& kt : global)
    {
      if( kt.second > 255 ) fatal_error() << "molecule type "<<kt.second<<" out of range"<<std::endl;
      keys.push_back(kt.first); types.push_back( kt.second );
    }

    for(size_t c=0;c<n_cells;c++) if( ! grid.is_ghost_cell(c) )
    {
      const size_t n = cells[c].size();
      for(size_t p=0;p<n;p++)
      {
        const uint64_t key = key_type(c,p).first;
        cells[c][idmol_field][p] = std::lower_bound( keys.begin(), keys.end(), key ) - keys.begin();
      }
    }
    return types;
  }

}
