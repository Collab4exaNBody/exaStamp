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

// Flat atoms and CSR bonds for the staged (gather -> kernels -> scatter) machine-learning potentials.
//
// Flat atom index = grid cell particle offset + index in cell, over all cells (ghost cells included). The central
// (non-ghost) atoms get the bonds to their neighbours within the model cutoff, r_ij = r_j - r_i in the real frame
// (domain xform applied); ghosts get none. The kernels write forces on all atoms (ghosts included), per-atom energies
// and per-atom virials (the whole bond virial on the central atom); mlip_scatter_outputs adds them into the grid
// fields, so that the force operator needs, in the simulation:
//
//   compute_force_prolog: [ zero_force_energy: { ghost: true } ]
//   compute_force_epilog: [ update_force_energy_from_ghost, force_to_accel ]

#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/grid_fields.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/particle_neighbors/chunk_neighbors_apply.h>
#include <exaStamp/unit_system.h>

#include <onika/log.h>
#include <onika/memory/allocator.h>

#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

namespace exaStamp
{
  namespace mlip_utils
  {
    // CSR bonds of the central atoms, read by the kernels (device addressable when the arrays are)
    struct BondsView
    {
      int inum = 0;                     // number of central atoms
      const int* ilist = nullptr;       // [inum] flat index of the central atoms
      int n_nodes = 0;                  // number of flat atoms (ghosts included)
      const int* type = nullptr;        // [n_nodes] model type (0 .. ntypes-1), -1: not in the model
      const int64_t* offset = nullptr;  // [n_nodes+1] bonds of atom i: offset[i] .. offset[i+1]
      const int* j = nullptr;           // [nbonds] flat index of the neighbour
      const double* rij = nullptr;      // [nbonds][3] r_j - r_i
    };

    // outputs, in eV and eV/ang (added into, zeroed by the gather)
    struct OutputsView
    {
      double* f = nullptr;              // [n_nodes][3]
      double* eatom = nullptr;          // [n_nodes]
      double* vatom = nullptr;          // [n_nodes][6] xx yy zz xy xz yz, or nullptr
    };

    // storage of the flat arrays, kept from step to step
    struct FlatBonds
    {
      onika::memory::CudaMMVector<int> type;
      onika::memory::CudaMMVector<int> ilist;
      onika::memory::CudaMMVector<int64_t> offset;
      onika::memory::CudaMMVector<int> nbh;
      onika::memory::CudaMMVector<double> rij;
      onika::memory::CudaMMVector<double> f;
      onika::memory::CudaMMVector<double> eatom;
      onika::memory::CudaMMVector<double> vatom;
    };

    // Gathers the flat atoms and the bonds of the central atoms within rcut. species_type: exaStamp species index ->
    // model type (-1: not in the model, such atoms are neither centres nor neighbours).
    template<class GridT>
    static inline void gather_bonds( GridT& grid, const exanb::GridChunkNeighbors& chunk_neighbors, const exanb::Domain& domain,
                                     const std::vector<int>& species_type, double rcut, bool want_vatom,
                                     FlatBonds& fb, BondsView& bv, OutputsView& out, int nth )
    {
      using namespace exanb;
      grid.rebuild_particle_offsets();
      const size_t ncells = grid.number_of_cells();
      const size_t* poff = grid.cell_particle_offset_data();
      const size_t n_nodes = grid.number_of_particles();
      if( n_nodes >= size_t(std::numeric_limits<int>::max()) ) fatal_error() << "MLIP flat bonds: too many atoms for int indices" << std::endl;
      auto* cells = grid.cells();
      const IJK dims = grid.dimension();
      const Mat3d xform = domain.xform();
      const bool identity = domain.xform_is_identity();
      const double rc2 = rcut * rcut;

      fb.type.resize( n_nodes );
      fb.offset.assign( n_nodes + 1 , 0 );
      int* type = fb.type.data();
      int64_t* offset = fb.offset.data();
      const int nspecies = species_type.size();
      const int* sp_type = species_type.data();

      // pass 1: types, bond counts
      size_t ncentral = 0;
#     pragma omp parallel for schedule(static) num_threads(nth)
      for(size_t c=0;c<ncells;c++)
      {
        const size_t np = cells[c].size();
        const auto* __restrict__ ty = cells[c][field::type];
        for(size_t p=0;p<np;p++)
        {
          const int s = ty[p];
          type[poff[c]+p] = ( s >= 0 && s < nspecies ) ? sp_type[s] : -1;
        }
      }
#     pragma omp parallel for schedule(dynamic) reduction(+:ncentral) num_threads(nth)
      for(size_t c=0;c<ncells;c++)
      {
        if( grid.is_ghost_cell(c) ) continue;
        const size_t np = cells[c].size();
        ncentral += np;
        const auto* __restrict__ rx = cells[c][field::rx];
        const auto* __restrict__ ry = cells[c][field::ry];
        const auto* __restrict__ rz = cells[c][field::rz];
        apply_cell_particle_neighbors( grid, chunk_neighbors, c, grid_index_to_ijk(dims,c), std::false_type{},
          [&](unsigned int p, size_t cb, unsigned int pb, size_t)
          {
            if( type[poff[c]+p] < 0 || type[poff[cb]+pb] < 0 ) return;
            Vec3d d = { cells[cb][field::rx][pb] - rx[p], cells[cb][field::ry][pb] - ry[p], cells[cb][field::rz][pb] - rz[p] };
            if( !identity ) d = xform * d;
            if( norm2(d) < rc2 ) ++ offset[ poff[c] + p + 1 ];
          } );
      }

      // central atoms, CSR offsets
      fb.ilist.resize( ncentral );
      {
        size_t k = 0;
        for(size_t c=0;c<ncells;c++)
        {
          if( grid.is_ghost_cell(c) ) continue;
          const size_t np = cells[c].size();
          for(size_t p=0;p<np;p++) fb.ilist[k++] = int( poff[c] + p );
        }
      }
      for(size_t i=0;i<n_nodes;i++) offset[i+1] += offset[i];
      const size_t nbonds = offset[n_nodes];
      fb.nbh.resize( nbonds );
      fb.rij.resize( 3 * nbonds );
      int* nbh = fb.nbh.data();
      double* rij = fb.rij.data();

      // pass 2: bonds
#     pragma omp parallel for schedule(dynamic) num_threads(nth)
      for(size_t c=0;c<ncells;c++)
      {
        if( grid.is_ghost_cell(c) ) continue;
        const auto* __restrict__ rx = cells[c][field::rx];
        const auto* __restrict__ ry = cells[c][field::ry];
        const auto* __restrict__ rz = cells[c][field::rz];
        int64_t q = 0;
        unsigned int plast = ~0u;
        apply_cell_particle_neighbors( grid, chunk_neighbors, c, grid_index_to_ijk(dims,c), std::false_type{},
          [&](unsigned int p, size_t cb, unsigned int pb, size_t)
          {
            if( p != plast ) { q = offset[ poff[c] + p ]; plast = p; }
            if( type[poff[c]+p] < 0 || type[poff[cb]+pb] < 0 ) return;
            Vec3d d = { cells[cb][field::rx][pb] - rx[p], cells[cb][field::ry][pb] - ry[p], cells[cb][field::rz][pb] - rz[p] };
            if( !identity ) d = xform * d;
            if( norm2(d) < rc2 )
            {
              nbh[q] = int( poff[cb] + pb );
              rij[3*q] = d.x; rij[3*q+1] = d.y; rij[3*q+2] = d.z;
              ++ q;
            }
          } );
      }

      fb.f.assign( 3 * n_nodes , 0.0 );
      fb.eatom.assign( n_nodes , 0.0 );
      if( want_vatom ) fb.vatom.assign( 6 * n_nodes , 0.0 );
      else fb.vatom.clear();

      bv = BondsView{};
      bv.inum = int( ncentral );
      bv.ilist = fb.ilist.data();
      bv.n_nodes = int( n_nodes );
      bv.type = type;
      bv.offset = offset;
      bv.j = nbh;
      bv.rij = rij;
      out = OutputsView{};
      out.f = fb.f.data();
      out.eatom = fb.eatom.data();
      out.vatom = want_vatom ? fb.vatom.data() : nullptr;
    }

    // model outputs (eV, eV/ang) added into the grid fields, ghosts included (folded back by update_*_from_ghost)
    template<class GridT>
    static inline void scatter_outputs( GridT& grid, const FlatBonds& fb, int nth )
    {
      using namespace exanb;
      static constexpr bool has_virial_field = GridHasField<GridT, field::_virial>::value;
      static constexpr double conv_energy = EXASTAMP_CONST_QUANTITY( 1. * eV );
      const size_t ncells = grid.number_of_cells();
      const size_t* poff = grid.cell_particle_offset_data();
      auto* cells = grid.cells();
      const double* f = fb.f.data();
      const double* ea = fb.eatom.data();
      const double* va = fb.vatom.empty() ? nullptr : fb.vatom.data();
#     pragma omp parallel for schedule(static) num_threads(nth)
      for(size_t c=0;c<ncells;c++)
      {
        const size_t np = cells[c].size();
        auto* __restrict__ fx = cells[c][field::fx];
        auto* __restrict__ fy = cells[c][field::fy];
        auto* __restrict__ fz = cells[c][field::fz];
        auto* __restrict__ ep = cells[c][field::ep];
        for(size_t p=0;p<np;p++)
        {
          const size_t i = poff[c] + p;
          fx[p] += f[3*i] * conv_energy;
          fy[p] += f[3*i+1] * conv_energy;
          fz[p] += f[3*i+2] * conv_energy;
          ep[p] += ea[i] * conv_energy;
          if constexpr ( has_virial_field )
          {
            if( va != nullptr )
            {
              const double* v = va + 6*i;    // xx yy zz xy xz yz
              cells[c][field::virial][p] += Mat3d{ v[0], v[3], v[4], v[3], v[1], v[5], v[4], v[5], v[2] } * conv_energy;
            }
          }
        }
      }
    }

  }
}
