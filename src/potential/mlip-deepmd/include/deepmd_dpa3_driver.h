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

// Phases of a DPA3 evaluation (see deepmd_dpa3.h) over the atoms of a grid: CSR bonds of the central atoms (flat atom
// indices, ghosts included), node features of every layer with ghost updates in between, fitting network, reverse
// phases with the node gradients folded from the ghosts, then forces, energies and per-atom virials into the grid.

#include "deepmd_dpa3.h"
#include "deepmd_force_op.h"

#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  struct DPA3Buffers
  {
    std::vector<int> type, ilist, nbh;
    std::vector<int64_t> offset;
    std::vector<double> rij, gb, f, eatom, vatom;
    std::vector< std::vector<double> > N, G;        // node features / gradients per layer [n_nodes][n_dim]
    std::vector<unsigned char> locks;
    // store mode: edge / angle inputs of every layer and the gradient of the current one, per central atom
    std::vector<int64_t> eoff, aoff;                // [inum + 1], in values (edges x e_dim, angle pairs x a_dim)
    std::vector< std::vector<double> > SE, SA;      // [layer][..]
    std::vector<double> dSE, dSA;
    bool stored = false;
  };

  struct DPA3RowLock
  {
    unsigned char* l;
    inline explicit DPA3RowLock(unsigned char* p) : l(p) { while( __atomic_test_and_set( l, __ATOMIC_ACQUIRE ) ) {} }
    inline ~DPA3RowLock() { __atomic_clear( l, __ATOMIC_RELEASE ); }
  };

  template<class GridT>
  inline void dpa3_compute( GridT& grid, const GridChunkNeighbors& chunk_neighbors, const Domain& domain, const DPA3Model& m,
                            const DeepmdModelView& fit, const GhostCommunicationScheme* scheme, MPI_Comm comm, DPA3Buffers& B,
                            double store_gb = 4.0 )
  {
    static constexpr bool has_virial_field = GridHasField<GridT, field::_virial>::value;
    static constexpr double conv_energy = EXASTAMP_CONST_QUANTITY( 1. * eV );
    grid.rebuild_particle_offsets();
    const size_t ncells = grid.number_of_cells();
    const size_t* poff = grid.cell_particle_offset_data();
    const size_t n_nodes = grid.number_of_particles();
    if( n_nodes >= size_t( std::numeric_limits<int>::max() ) ) fatal_error() << "DPA3: too many atoms for int indices" << std::endl;
    auto* cells = grid.cells();
    const IJK dims = grid.dimension();
    const Mat3d xform = domain.xform();
    const bool identity = domain.xform_is_identity();
    const double rc2 = m.rc * m.rc;
    const int nd = m.n_dim, NL = m.layers.size();

    // ---- CSR bonds of the central atoms within rcut, r_ij = r_j - r_i (real frame)
    B.type.resize( n_nodes );
    B.offset.assign( n_nodes + 1 , 0 );
    size_t ncentral = 0;
#   pragma omp parallel for schedule(dynamic) reduction(+:ncentral)
    for(size_t c = 0; c < ncells; c++)
    {
      const size_t np = cells[c].size();
      const auto* __restrict__ ty = cells[c][field::type];
      for(size_t p = 0; p < np; p++) B.type[poff[c] + p] = fit.species_type[ ty[p] ];
      if( grid.is_ghost_cell(c) ) continue;
      ncentral += np;
      const auto* __restrict__ rx = cells[c][field::rx];
      const auto* __restrict__ ry = cells[c][field::ry];
      const auto* __restrict__ rz = cells[c][field::rz];
      apply_cell_particle_neighbors( grid, chunk_neighbors, c, grid_index_to_ijk(dims, c), std::false_type{},
        [&](unsigned int p, size_t cb, unsigned int pb, size_t)
        {
          Vec3d d = { cells[cb][field::rx][pb] - rx[p], cells[cb][field::ry][pb] - ry[p], cells[cb][field::rz][pb] - rz[p] };
          if( !identity ) d = xform * d;
          if( norm2(d) < rc2 ) ++ B.offset[ poff[c] + p + 1 ];
        } );
    }
    B.ilist.resize( ncentral );
    {
      size_t k = 0;
      for(size_t c = 0; c < ncells; c++)
      {
        if( grid.is_ghost_cell(c) ) continue;
        for(size_t p = 0; p < cells[c].size(); p++) B.ilist[k++] = int( poff[c] + p );
      }
    }
    for(size_t i = 0; i < n_nodes; i++) B.offset[i+1] += B.offset[i];
    const size_t nbonds = B.offset[n_nodes];
    B.nbh.resize( nbonds );
    B.rij.resize( 3 * nbonds );
#   pragma omp parallel for schedule(dynamic)
    for(size_t c = 0; c < ncells; c++)
    {
      if( grid.is_ghost_cell(c) ) continue;
      const auto* __restrict__ rx = cells[c][field::rx];
      const auto* __restrict__ ry = cells[c][field::ry];
      const auto* __restrict__ rz = cells[c][field::rz];
      int64_t q = 0;
      unsigned int plast = ~0u;
      apply_cell_particle_neighbors( grid, chunk_neighbors, c, grid_index_to_ijk(dims, c), std::false_type{},
        [&](unsigned int p, size_t cb, unsigned int pb, size_t)
        {
          if( p != plast ) { q = B.offset[ poff[c] + p ]; plast = p; }
          Vec3d d = { cells[cb][field::rx][pb] - rx[p], cells[cb][field::ry][pb] - ry[p], cells[cb][field::rz][pb] - rz[p] };
          if( !identity ) d = xform * d;
          if( norm2(d) < rc2 )
          {
            B.nbh[q] = int( poff[cb] + pb );
            B.rij[3*q] = d.x; B.rij[3*q+1] = d.y; B.rij[3*q+2] = d.z;
            ++ q;
          }
        } );
    }
    const int inum = ncentral;

    // ---- node features of every layer, gradients
    B.N.resize( NL + 1 );
    B.G.resize( NL + 1 );
    for(auto& v : B.N) v.assign( n_nodes * nd , 0.0 );
    for(size_t a = 0; a < n_nodes; a++)
      for(int c = 0; c < nd; c++) B.N[0][a*nd + c] = dpa3_act( m, m.tebd[ size_t( B.type[a] ) * nd + c ] );
    auto exchange = [&](std::vector<double>& buf, bool to_ghost)
    {
      if( scheme == nullptr ) fatal_error() << "DPA3: no ghost communication scheme (a ghost update must run before deepmd_force)" << std::endl;
      GhostFlatBuffer fb;
      fb.m_data = buf.data(); fb.m_row_size = nd; fb.m_row_stride = nd; fb.m_nrows = n_nodes;
      fb.m_type = GhostFlatBuffer::ValueType::Double;
      if( to_ghost ) flat_buffer_ghost_exchange<true>( comm, *scheme, poff, fb, 313 );
      else flat_buffer_ghost_exchange<false>( comm, *scheme, poff, fb, 314 );
    };
    const int64_t* offset = B.offset.data();
    const int* nbh = B.nbh.data();
    const double* rij = B.rij.data();

    // ---- store mode: per-atom edge / angle counts, offsets, if the store fits the budget
    B.stored = false;
    if( store_gb > 0.0 )
    {
      std::vector<int> cn( inum ), ca( inum );
#     pragma omp parallel
      {
        static thread_local DPA3AtomScratch S;
#       pragma omp for schedule(dynamic, 16)
        for(int q = 0; q < inum; q++) { dpa3_geometry( m, B.ilist[q], B.type[B.ilist[q]], offset, nbh, rij, S ); cn[q] = S.g.n; ca[q] = S.g.na; }
      }
      B.eoff.assign( inum + 1, 0 ); B.aoff.assign( inum + 1, 0 );
      for(int q = 0; q < inum; q++) { B.eoff[q+1] = B.eoff[q] + int64_t(cn[q]) * m.e_dim; B.aoff[q+1] = B.aoff[q] + int64_t(ca[q]) * ca[q] * m.a_dim; }
      const double bytes = 8.0 * double( B.eoff[inum] + B.aoff[inum] ) * ( NL + 1 );
      if( bytes <= store_gb * 1e9 )
      {
        B.stored = true;
        B.SE.resize( NL ); B.SA.resize( NL );
        for(int l = 0; l < NL; l++) { B.SE[l].resize( B.eoff[inum] ); B.SA[l].resize( B.aoff[inum] ); }
        B.dSE.assign( B.eoff[inum], 0.0 ); B.dSA.assign( B.aoff[inum], 0.0 );
      }
    }
    const bool stored = B.stored;
    // inputs of layer l of atom q (index in ilist) from the store, or recomputed (chain of layers 0 .. l-1)
    auto layer_inputs = [&](int q, int l, DPA3AtomScratch& S, int recompute_mode = 0)
    {
      const int i = B.ilist[q];
      if( int(S.L.size()) < NL + 1 ) S.L.resize( NL + 1 );
      DPA3Timer T;
      dpa3_geometry( m, i, B.type[i], offset, nbh, rij, S );
      T.mark(0);
      if( !stored || l == 0 )
      {
        dpa3_init_features( m, B.N[0].data() + size_t(i) * nd, S );
        if( stored ) { std::copy( S.L[0].e.begin(), S.L[0].e.begin() + ( B.eoff[q+1] - B.eoff[q] ), B.SE[0].begin() + B.eoff[q] );
                       std::copy( S.L[0].a.begin(), S.L[0].a.begin() + ( B.aoff[q+1] - B.aoff[q] ), B.SA[0].begin() + B.aoff[q] ); }
        else for(int k = 0; k < l; k++) dpa3_layer_forward( m, k, B.N[k].data(), S, recompute_mode );
        return;
      }
      DPA3Work& w = S.L[l];
      w.n.assign( B.N[l].begin() + size_t(i) * nd, B.N[l].begin() + size_t(i+1) * nd );
      w.e.assign( B.SE[l].begin() + B.eoff[q], B.SE[l].begin() + B.eoff[q+1] );
      w.a.assign( B.SA[l].begin() + B.aoff[q], B.SA[l].begin() + B.aoff[q+1] );
    };

    // ---- forward layers
    for(int l = 0; l < NL; l++)
    {
#     pragma omp parallel
      {
        static thread_local DPA3AtomScratch S;
#       pragma omp for schedule(dynamic, 1)
        for(int q = 0; q < inum; q++)
        {
          const int i = B.ilist[q];
          layer_inputs( q, l, S );
          dpa3_layer_forward( m, l, B.N[l].data(), S );
          for(int c = 0; c < nd; c++) B.N[l+1][size_t(i)*nd + c] = S.L[l+1].n[c];
          if( stored && l + 1 < NL )
          {
            std::copy( S.L[l+1].e.begin(), S.L[l+1].e.begin() + ( B.eoff[q+1] - B.eoff[q] ), B.SE[l+1].begin() + B.eoff[q] );
            std::copy( S.L[l+1].a.begin(), S.L[l+1].a.begin() + ( B.aoff[q+1] - B.aoff[q] ), B.SA[l+1].begin() + B.aoff[q] );
          }
        }
      }
      exchange( B.N[l+1], true );
    }

    // ---- fitting network: energies, dE/dn of the last layer
    B.eatom.assign( n_nodes, 0.0 );
    B.G[NL].assign( n_nodes * nd, 0.0 );
    const int din = nd + ( m.concat_tebd ? nd : 0 );
#   pragma omp parallel
    {
      std::vector<double> x( din ), dx( din ), g0( std::max( din, DP_MAX_FIT_WIDTH ) ), g1( std::max( din, DP_MAX_FIT_WIDTH ) );
      auto wk = std::make_unique< DPMlpWork<DP_MAX_FIT_LAYERS, DP_MAX_FIT_WIDTH> >();
#     pragma omp for schedule(static)
      for(int q = 0; q < inum; q++)
      {
        const int i = B.ilist[q], ti = B.type[i];
        for(int c = 0; c < nd; c++) x[c] = B.N[NL][size_t(i)*nd + c];
        if( m.concat_tebd ) for(int c = 0; c < nd; c++) x[nd + c] = m.tebd[ size_t(ti) * nd + c ];
        const DPNet& net = fit.fitting[ fit.fitting_mixed_types ? 0 : ti ];
        B.eatom[i] = dp_mlp_forward( net, fit.layers, fit.weights, x.data(), *wk )[0] + fit.energy_bias[ti];
        const double one = 1.0;
        dp_mlp_backward( net, fit.layers, fit.weights, x.data(), *wk, &one, dx.data(), g0.data(), g1.data() );
        for(int c = 0; c < nd; c++) B.G[NL][size_t(i)*nd + c] = dx[c];
      }
    }

    // ---- reverse phases
    B.gb.assign( 3 * nbonds, 0.0 );
    if( B.locks.size() != n_nodes ) B.locks.assign( n_nodes, 0 );
    for(int L = NL - 1; L >= 0; L--)
    {
      if( L > 0 ) B.G[L].assign( n_nodes * nd, 0.0 );
#     pragma omp parallel
      {
        static thread_local DPA3AtomScratch S;
#       pragma omp for schedule(dynamic, 1)
        for(int q = 0; q < inum; q++)
        {
          const int i = B.ilist[q];
          // store mode: layer L only, the gradients of its outputs from the store; else the whole chain, the products
          // from the top layer down to L
          const int ltop = stored ? L : NL - 1;
          layer_inputs( q, stored ? L : NL - 1, S, 2 );
          if( stored ) dpa3_layer_forward( m, L, B.N[L].data(), S, 1 );
          else dpa3_layer_forward( m, NL - 1, B.N[NL-1].data(), S, 2 );
          const int n = S.g.n, na = S.g.na;
          S.gg.reset( n, na );
          dpa3_grow( S.de, size_t(n) * m.e_dim ); dpa3_grow( S.da, size_t(na) * na * m.a_dim ); dpa3_grow( S.dn, nd );
          if( stored && L < NL - 1 )
          {
            std::copy( B.dSE.begin() + B.eoff[q], B.dSE.begin() + B.eoff[q+1], S.de.begin() );
            std::copy( B.dSA.begin() + B.aoff[q], B.dSA.begin() + B.aoff[q+1], S.da.begin() );
          }
          else
          {
            std::fill( S.de.begin(), S.de.begin() + size_t(n) * m.e_dim, 0.0 );
            std::fill( S.da.begin(), S.da.begin() + size_t(na) * na * m.a_dim, 0.0 );
          }
          for(int l = ltop; l >= L; l--)
          {
            const double* Gt = B.G[l+1].data() + size_t(i) * nd;
            for(int c = 0; c < nd; c++) S.dn[c] = Gt[c];
            dpa3_layer_vjp( m, l, B.N[l].data(), S, l == L );
            if( l > L )
            {
              std::copy( S.de_in.begin(), S.de_in.begin() + size_t(n) * m.e_dim, S.de.begin() );
              std::copy( S.da_in.begin(), S.da_in.begin() + size_t(na) * na * m.a_dim, S.da.begin() );
            }
          }
          if( L > 0 )
          {
            {
              DPA3RowLock lk( &B.locks[i] );
              double* gi = B.G[L].data() + size_t(i) * nd;
              for(int c = 0; c < nd; c++) gi[c] += S.dn_in[c];
            }
            for(int k = 0; k < n; k++)
            {
              const int j = S.g.j[k];
              DPA3RowLock lk( &B.locks[j] );
              double* gj = B.G[L].data() + size_t(j) * nd;
              const double* dj = S.dnj.data() + size_t(k) * nd;
              for(int c = 0; c < nd; c++) gj[c] += dj[c];
            }
          }
          DPA3Timer T2;
          if( L == 0 ) dpa3_init_vjp( m, S.de_in.data(), S.da_in.data(), S );
          if( stored && L > 0 )
          {
            std::copy( S.de_in.begin(), S.de_in.begin() + ( B.eoff[q+1] - B.eoff[q] ), B.dSE.begin() + B.eoff[q] );
            std::copy( S.da_in.begin(), S.da_in.begin() + ( B.aoff[q+1] - B.aoff[q] ), B.dSA.begin() + B.aoff[q] );
          }
          dpa3_geom_to_bonds( m, S, B.gb.data() );
          T2.mark(13);
        }
      }
      if( L > 0 ) exchange( B.G[L], false );
    }

    // ---- forces (dE/dd per bond, d = r_j - r_i), per-atom virials on the central atom, into the grid
    B.f.assign( 3 * n_nodes, 0.0 );
    B.vatom.assign( 9 * n_nodes, 0.0 );
    std::vector<unsigned char>& lk = B.locks;
#   pragma omp parallel for schedule(dynamic, 16)
    for(int q = 0; q < inum; q++)
    {
      const int i = B.ilist[q];
      double fi[3] = { 0, 0, 0 }, vi[9] = { 0,0,0, 0,0,0, 0,0,0 };
      for(int64_t b = offset[i]; b < offset[i+1]; b++)
      {
        const double* g = B.gb.data() + 3 * b;
        if( g[0] == 0.0 && g[1] == 0.0 && g[2] == 0.0 ) continue;
        const double* d = rij + 3 * b;
        for(int x = 0; x < 3; x++) fi[x] += g[x];
        for(int x = 0; x < 3; x++) for(int y = 0; y < 3; y++) vi[3*x + y] += -g[x] * d[y];
        const int j = nbh[b];
        DPA3RowLock l( &lk[j] );
        for(int x = 0; x < 3; x++) B.f[3*size_t(j) + x] -= g[x];
      }
      DPA3RowLock l( &lk[i] );
      for(int x = 0; x < 3; x++) B.f[3*size_t(i) + x] += fi[x];
      for(int x = 0; x < 9; x++) B.vatom[9*size_t(i) + x] += vi[x];
    }
#   pragma omp parallel for schedule(static)
    for(size_t c = 0; c < ncells; c++)
    {
      const size_t np = cells[c].size();
      auto* __restrict__ fx = cells[c][field::fx];
      auto* __restrict__ fy = cells[c][field::fy];
      auto* __restrict__ fz = cells[c][field::fz];
      auto* __restrict__ ep = cells[c][field::ep];
      for(size_t p = 0; p < np; p++)
      {
        const size_t i = poff[c] + p;
        fx[p] += B.f[3*i] * conv_energy;
        fy[p] += B.f[3*i+1] * conv_energy;
        fz[p] += B.f[3*i+2] * conv_energy;
        ep[p] += B.eatom[i] * conv_energy;
        if constexpr ( has_virial_field )
        {
          const double* v = B.vatom.data() + 9*i;
          cells[c][field::virial][p] += Mat3d{ v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8] } * conv_energy;
        }
      }
    }
  }

}
