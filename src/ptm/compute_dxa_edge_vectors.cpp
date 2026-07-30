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

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <exaStamp/delaunay/delaunay_tessellation.h>
#include <exaStamp/delaunay/dxa_edge_vectors.h>

// Reference/ideal lattice templates only -- refdata_t's const data tables (structure_bcc.points
// etc), no ptm_index()/ptm_initialize_* calls, so this doesn't link exaStampPTM's compiled
// library any differently than compute_ptm.cu already does, just reuses its headers (this file
// lives in src/ptm/, not src/delaunay/, specifically so ptm_constants.h/ptm_initialize_data.h --
// a flat, non-namespaced include dir wired only into this plugin -- are reachable without adding
// a cross-plugin include path just for this).
#include <ptm_constants.h>
#include <ptm_initialize_data.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// DXA pipeline step (iii) (Stukowski, Bulatov, Arsenlis, Model. Simul. Mater. Sci. Eng. 20 (2012)
// 085007): assign an ideal lattice vector to each Delaunay tessellation edge, w.r.t. a single
// user-chosen target_structure (e.g. "BCC") -- NOT whatever PTM happened to locally match, so the
// same target can be probed even where other structures/defects coexist. For each edge (u,v):
//   - both endpoints must have PTM_MATCH_<target_structure> as their own local structure type
//     (ptm_type field, from ptm_fields) -- otherwise the edge is left unresolved.
//   - rotate the actual bond vector into the LOWER-INDEXED endpoint's own "ideal lattice frame"
//     (ptm_orientation rotation tensor, R maps ideal->actual, so actual->ideal is R^T for a
//     rotation), and snap to the target structure's own nearest ideal neighbor direction (PTM's
//     own reference templates, ptm_initialize_data.h -- exactly what PTM matched atoms against,
//     not re-derived). angle_tolerance is now just a sanity bound on THIS one-sided snap (reject
//     only if even the nearest ideal direction is implausibly far off), not a requirement that
//     the OTHER endpoint independently agrees.
//
//   Originally this DID require both endpoints to independently agree on the same direction
//   within angle_tolerance (15 deg) -- checked against two independent reference DXA
//   implementations (LAMMPS fix_disloc's 2014 CNA-triangulation method, and Stukowski's own
//   2010 CNA-based DXA predecessor, DXA1.3.6) after finding the resulting "bad" (defective)
//   region on a real dislocation quadrupole test case came out ~3x wider (by both triangle count
//   and surface area) than OVITO's own DXA output on the identical configuration. Neither
//   reference algorithm has any per-edge cross-atom orientation-agreement gate -- both classify
//   each atom's structure independently, propagate the ideal-lattice mapping one-sidedly (stamp
//   the matching atom's own orientation onto its own edges, no rejection at assignment time), and
//   defer ALL defect detection to a later global circuit-closure test. A per-edge two-sided
//   agreement gate is far more sensitive to a dislocation's genuine long-range elastic
//   lattice-rotation gradient than either reference design, which is exactly why it was
//   overclassifying strained-but-not-actually-defective atoms as "unresolved" well outside the
//   true disordered core. Fixed to match: resolution is now per-edge/one-sided; defect detection
//   is deferred to compute_dxa_burgers_circuits' global spanning-tree circuit closure, which
//   already implements exactly this "defer to the end" principle.
namespace exaStamp
{
  using namespace exanb;

  static const ptm::refdata_t* dxa_target_structure_refdata( const std::string& name )
  {
    if( name == "FCC" ) return &ptm::structure_fcc;
    if( name == "HCP" ) return &ptm::structure_hcp;
    if( name == "BCC" ) return &ptm::structure_bcc;
    if( name == "ICO" ) return &ptm::structure_ico;
    if( name == "SC"  ) return &ptm::structure_sc;
    return nullptr;
  }

  template<class GridT>
  class ComputeDXAEdgeVectors : public OperatorNode
  {
    ADD_SLOT( GridT                , grid                 , INPUT , REQUIRED );
    ADD_SLOT( Domain               , domain               , INPUT , REQUIRED );
    ADD_SLOT( DelaunayTessellation , delaunay_tessellation, INPUT , REQUIRED );
    ADD_SLOT( DXAEdgeVectors       , dxa_edge_vectors     , OUTPUT );
    ADD_SLOT( std::string          , target_structure     , INPUT , std::string("BCC") , DocString{"Crystal structure edges are classified against: FCC, HCP, BCC, ICO or SC. Only edges whose both endpoints locally match this structure (ptm_type) can be resolved."} );
    ADD_SLOT( double               , angle_tolerance      , INPUT , 40.0 , DocString{"Sanity bound (degrees): max angular deviation between an edge's bond direction (rotated into its lower-indexed endpoint's own ideal lattice frame) and the nearest ideal neighbor direction, for the edge to be considered resolved. NOT a cross-atom agreement check -- the other endpoint's own orientation doesn't need to agree (see file header comment); this only rejects a snap so far off it's ambiguous which ideal direction it's even closest to."} );
    ADD_SLOT( std::string          , struct_field         , INPUT , std::string("ptm_type")        , DocString{"Name of the per-particle structure-type grid field (written by ptm_fields)"} );
    ADD_SLOT( std::string          , orient_field         , INPUT , std::string("ptm_orientation") , DocString{"Name of the per-particle lattice-orientation rotation-tensor grid field (written by ptm_fields)"} );
    ADD_SLOT( long                 , n_edges_resolved     , OUTPUT , DocString{"Number of tessellation edges resolved against target_structure (out of dxa_edge_vectors->edges.size())"} );

  public:
    inline void execute () override final
    {
      const ptm::refdata_t* refdata = dxa_target_structure_refdata( *target_structure );
      if( refdata == nullptr )
      {
        fatal_error() << "compute_dxa_edge_vectors: unknown target_structure '" << *target_structure
                      << "', expected one of FCC, HCP, BCC, ICO, SC" << std::endl;
      }
      const double target_type = static_cast<double>( refdata->type );
      const int num_nbrs = refdata->num_nbrs;

      // ideal neighbor directions for the target structure, unit-normalized once (points[0][0]
      // is the template's own central atom, at the origin -- neighbors start at index 1)
      std::vector<Vec3d> ideal_dir( num_nbrs );
      for(int k=0;k<num_nbrs;k++)
      {
        const auto& p = refdata->points[0][k+1];
        Vec3d v { p[0], p[1], p[2] };
        ideal_dir[k] = v / norm(v);
      }

      // flatten ptm_type/ptm_orientation into per-particle arrays, same convention as
      // write_delaunay_vtk's color_fields (cell_particle_offset[cell]+particle)
      auto cells = grid->cells_accessor();
      auto struct_acc = grid->field_const_accessor( field::mk_generic_real( *struct_field ) );
      auto orient_acc = grid->field_const_accessor( field::mk_generic_mat3( *orient_field ) );
      const size_t * const cell_particle_offset = grid->cell_particle_offset_data();
      const size_t n_cells_grid = grid->number_of_cells();
      const size_t n_particles = grid->number_of_particles();

      std::vector<double> flat_struct_type( n_particles, 0.0 );
      std::vector<Mat3d> flat_orient( n_particles, make_identity_matrix() );
      for(size_t c=0;c<n_cells_grid;c++)
      {
        const size_t np = cells[c].size();
        for(size_t p=0;p<np;p++)
        {
          const size_t i = cell_particle_offset[c] + p;
          flat_struct_type[i] = cells[c][struct_acc][p];
          flat_orient[i] = cells[c][orient_acc][p];
        }
      }

      const DelaunayTessellation& mesh = *delaunay_tessellation;
      DXAEdgeVectors& result = *dxa_edge_vectors;
      result.edges.clear();
      result.ideal_vector.clear();
      result.resolved.clear();
      result.edge_index.clear();

      // deduplicate edges: every tetrahedron contributes its 6 vertex pairs, sort+unique
      static constexpr int edge_lv[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };
      std::vector<std::array<uint32_t,2>> all_edges;
      all_edges.reserve( mesh.tetrahedra.size() * 6 );
      for(const auto& tet : mesh.tetrahedra)
      {
        for(int e=0;e<6;e++)
        {
          uint32_t a = tet[ edge_lv[e][0] ];
          uint32_t b = tet[ edge_lv[e][1] ];
          all_edges.push_back( a<b ? std::array<uint32_t,2>{a,b} : std::array<uint32_t,2>{b,a} );
        }
      }
      std::sort( all_edges.begin(), all_edges.end() );
      all_edges.erase( std::unique( all_edges.begin(), all_edges.end() ), all_edges.end() );

      result.edges = std::move( all_edges );
      result.ideal_vector.resize( result.edges.size() );
      result.resolved.resize( result.edges.size() );

      result.vertex_matches_target.resize( mesh.vertices.size() );
      for(size_t v=0; v<mesh.vertices.size(); v++)
      {
        result.vertex_matches_target[v] = ( flat_struct_type[ mesh.vertex_particle_index[v] ] == target_type ) ? 1 : 0;
      }

      size_t n_resolved = 0;
      const double angle_tol = *angle_tolerance;
      for(size_t e=0;e<result.edges.size();e++)
      {
        const uint32_t vu = result.edges[e][0];
        const uint32_t vv = result.edges[e][1];
        result.edge_index[ DXAEdgeVectors::key(vu,vv) ] = static_cast<uint32_t>(e);

        const uint32_t pu = mesh.vertex_particle_index[vu];
        const uint32_t pv = mesh.vertex_particle_index[vv];

        result.resolved[e] = false;
        result.ideal_vector[e] = Vec3d{0.,0.,0.};

        if( flat_struct_type[pu] != target_type || flat_struct_type[pv] != target_type ) { continue; }

        const Vec3d d = domain->xform() * ( mesh.vertices[vv] - mesh.vertices[vu] );
        const double dlen = norm(d);
        if( dlen <= 0.0 ) { continue; }
        const Vec3d d_hat = d / dlen;

        // one-sided: only the lower-indexed endpoint's own frame decides this edge's ideal
        // vector -- vv's own orientation is not consulted here at all (see file header comment).
        const Vec3d d_ideal_u = transpose( flat_orient[pu] ) * d_hat;

        int best_k = -1;
        double best_dot = -2.0;
        for(int k=0;k<num_nbrs;k++)
        {
          const double dot = d_ideal_u.x*ideal_dir[k].x + d_ideal_u.y*ideal_dir[k].y + d_ideal_u.z*ideal_dir[k].z;
          if( dot > best_dot ) { best_dot = dot; best_k = k; }
        }
        const double angle_deg = std::acos( std::clamp(best_dot,-1.0,1.0) ) * (180.0/M_PI);

        if( angle_deg <= angle_tol )
        {
          result.resolved[e] = true;
          const auto& p = refdata->points[0][best_k+1];
          result.ideal_vector[e] = Vec3d{ p[0], p[1], p[2] };
          ++n_resolved;
        }
      }

      *n_edges_resolved = static_cast<long>( n_resolved );
      lout << "compute_dxa_edge_vectors: " << n_resolved << " / " << result.edges.size()
           << " edges resolved against " << *target_structure << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

DXA pipeline step (iii): assigns an ideal lattice vector to each Delaunay tessellation edge
(compute_delaunay), w.r.t. a single user-chosen target_structure (FCC/HCP/BCC/ICO/SC) -- not
whatever PTM happened to locally match, so a specific structure can be probed even in a system
with mixed/defective regions. An edge is resolved if both endpoints locally match target_structure
(ptm_type, from ptm_fields) and its bond direction, rotated into the LOWER-INDEXED endpoint's own
ideal lattice frame (ptm_orientation) only, snaps unambiguously (within angle_tolerance) to one of
PTM's own reference-template directions -- the other endpoint's own orientation is not consulted,
by design (see file header comment: this used to require both endpoints to independently agree,
which turned out far more sensitive to a dislocation's ordinary long-range elastic
lattice-rotation field than either reference DXA implementation checked against). Defect detection
itself is deferred entirely to compute_dxa_burgers_circuits' global spanning-tree circuit closure.

Usage example:

compute_ptm: { rcut: 3.6 ang }
ptm_fields: {}
compute_delaunay: {}
compute_dxa_edge_vectors: { target_structure: BCC, angle_tolerance: 40.0 }

)EOF";
    }
  };

  // DXA pipeline step (iv): a tetrahedron is "good" (part of the undistorted target-structure
  // lattice) iff all 4 of its vertices individually match target_structure
  // (DXAEdgeVectors::vertex_matches_target, from step iii) -- "bad" (a defect: dislocation core,
  // grain boundary, stacking fault, second phase, ...) otherwise. No grid access needed at all
  // (purely a DelaunayTessellation + DXAEdgeVectors post-process), so this is a plain
  // (non-grid-variant) operator.
  //
  // First attempt (superseded, see git history) additionally required all six edges to resolve
  // AND their six ideal vectors to close consistently around the tetrahedron (an edge-vector
  // closure test). Checked against two independent reference DXA implementations
  // (LAMMPS fix_disloc's 2014 CNA-triangulation method, and Stukowski's own 2010 CNA-based DXA
  // predecessor, DXA1.3.6) after finding the resulting "bad" region on a real dislocation
  // quadrupole test case came out much wider than OVITO's own DXA output on the identical
  // configuration. Neither reference defines "good/bad" via anything like an edge-vector-closure
  // test -- both use exactly this simpler per-atom criterion (their atomic-structure-typing
  // result alone) and defer ALL consistency/defect detection to a later, global circuit-closure
  // test -- compute_dxa_burgers_circuits' job here.
  class ComputeDXATetClassification : public OperatorNode
  {
    ADD_SLOT( DelaunayTessellation , delaunay_tessellation , INPUT , REQUIRED );
    ADD_SLOT( DXAEdgeVectors       , dxa_edge_vectors       , INPUT , REQUIRED );
    ADD_SLOT( DXATetClassification , dxa_tet_classification , OUTPUT );
    ADD_SLOT( long                 , n_tets_good            , OUTPUT , DocString{"Number of tetrahedra classified good (out of delaunay_tessellation->tetrahedra.size())"} );

  public:
    inline void execute () override final
    {
      const DelaunayTessellation& mesh = *delaunay_tessellation;
      const DXAEdgeVectors& ev = *dxa_edge_vectors;
      const size_t n_tets = mesh.tetrahedra.size();

      DXATetClassification& result = *dxa_tet_classification;
      result.good.assign( n_tets, 0.0 );

      size_t n_good = 0;
      for(size_t t=0;t<n_tets;t++)
      {
        const auto& tet = mesh.tetrahedra[t];
        const bool all_match = ev.vertex_matches_target[tet[0]] && ev.vertex_matches_target[tet[1]]
                             && ev.vertex_matches_target[tet[2]] && ev.vertex_matches_target[tet[3]];
        if( all_match ) { result.good[t] = 1.0; ++n_good; }
      }

      *n_tets_good = static_cast<long>( n_good );
      lout << "compute_dxa_tet_classification: " << n_good << " / " << n_tets << " tetrahedra classified good" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

DXA pipeline step (iv): classifies each Delaunay tetrahedron (compute_delaunay) as "good" (part of
the undistorted target-structure lattice) or "bad" (part of a defect -- dislocation core, grain
boundary, stacking fault, second phase, ...): good means all 4 of its vertices individually match
target_structure (DXAEdgeVectors::vertex_matches_target, from compute_dxa_edge_vectors) -- not an
edge-vector-closure test (that overclassifies a dislocation's ordinary elastic strain field as
defective, see this file's own comment). Output is a plain 0/1 array parallel to the tessellation's
own tetrahedra, directly usable as write_delaunay_vtk's dxa_tet_classification input for
visualization.

Usage example:

compute_ptm: { rcut: 3.6 ang }
ptm_fields: {}
compute_delaunay: {}
compute_dxa_edge_vectors: { target_structure: BCC, angle_tolerance: 40.0 }
compute_dxa_tet_classification: {}
write_delaunay_vtk: { filename: "paraview/delaunay" }   # dxa_tet_classification auto-wires in

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_dxa_edge_vectors)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_dxa_edge_vectors", make_grid_variant_operator< ComputeDXAEdgeVectors > );
    OperatorNodeFactory::instance()->register_factory( "compute_dxa_tet_classification", make_simple_operator< ComputeDXATetClassification > );
  }

}
