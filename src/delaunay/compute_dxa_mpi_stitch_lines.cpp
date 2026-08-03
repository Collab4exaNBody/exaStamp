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

#include <exanb/core/domain.h>
#include <exaStamp/delaunay/dxa_dislocation_lines.h>

#include <mpi.h>
#include <algorithm>
#include <vector>
#include <array>
#include <utility>
#include <cmath>
#include <unordered_map>
#include <map>
#include <cstdint>
#include <functional>

// Stitches dislocation lines that compute_dxa_circuit_sweep left fragmented across an MPI domain-
// decomposition boundary back into single, continuous lines. compute_dxa_circuit_sweep runs
// entirely per-rank, on that rank's own local interface mesh only -- a real dislocation whose path
// crosses into a neighboring rank's own owned cells simply stops there (StopReason::OpenEdge,
// tagged as DXADislocationLines::open_boundary_front/back), reported as if it were a genuine dead
// end. This produced a real, user-visible symptom: dislocation/segment counts measurably grew with
// MPI rank count on the identical input (8/11/18 physical dislocations at 1/2/4 ranks on the
// quadrupole test case) purely from this fragmentation, not from any real change in the underlying
// physics.
//
// No OVITO mechanism to port here -- OVITO's own DXA implementation isn't domain-decomposed, so this
// is a new problem specific to running it inside exaStamp's MPI-parallel framework.
//
// Mechanism: gather every rank's own lines (points, core_size, Burgers vector, open_boundary flags,
// and each open-boundary end's own boundary-loop atom ids/positions) to rank 0 via MPI_Gatherv (same
// pattern write_dxa_ca_file/write_ovito_interface_mesh already use -- no existing multi-piece
// convention fits this either, since stitching is fundamentally a cross-rank operation).
//
// Matching is PRIMARILY by shared ghost atom id, not fuzzy position: two ranks' own independently-
// grown circuits can stop with genuinely different final loop shapes/centroids even at the exact
// same real crossing (measured: ~2.6-2.8 Å apart on the real quadrupole case) -- close enough to look
// like a match, but nowhere near the near-zero gap a real, physically identical crossing point should
// have. compute_dxa_circuit_sweep now records, for every OpenEdge-stopped end, the DelaunayTessellation
// ::vertex_global_id (a persistent, cross-rank-identical particle id, since exaStamp's own ghost
// duplication preserves id verbatim) of every mesh vertex in that end's own final loop. Two ends on
// DIFFERENT ranks that share even one such atom id are, by construction, looking at the exact same
// real physical location -- an EXACT match, no distance tolerance needed for correctness. The
// resulting "snap position" (see below) is the (very tightly clustered, since it's real ghost-shared
// coordinates) average of every atom both ends have in common. Only when an end shares NO atom id
// with anything (shouldn't normally happen, but not assumed impossible) does this fall back to the
// old fuzzy position+Burgers-vector proximity match (within `match_tolerance`, `burgers_tolerance`)
// as a strictly lower-confidence safety net. Matching is opportunistic either way, never required:
// an end that finds no partner at all (e.g. a genuine single-rank physical boundary, not an MPI cut)
// is simply left as a normal standalone end, unchanged from today's behavior.
//
// Closing the gap for real, not just picking a looser tolerance: once matched, the snap position is
// APPENDED as one new point between the two fragments (not used to overwrite either fragment's own
// last recorded point, same "don't discard real data" principle already established for the same-rank
// junction-endpoint reconciliation fix) -- so the assembled line passes exactly through the one
// verified-real crossing coordinate, with only ordinary point-to-point spacing on either side of it,
// rather than ever placing two supposedly-identical points measurably apart.
//
// Since each line has only 2 ends and each end matches at most 1 partner, the resulting match graph
// (nodes = lines, edges = matches) has max degree 2 per node -- a union of simple paths and simple
// cycles, nothing more complex is possible. Paths are walked from their one free end inward,
// concatenating each fragment's own points in the direction that keeps the assembled line
// continuous (reversing a fragment's own point order when arriving via its "back" end), inserting
// the snap point at each transition; a fully closed multi-rank loop is detected when a walk's next
// partner is the line it started from.
//
// Burgers vector for an assembled multi-fragment line is simply whichever fragment the walk started
// from -- no sign reconciliation performed, matching the same simplification compute_dxa_circuit_
// sweep's own do_merge already accepts for same-rank merges (its own `burgers` field is never
// revisited after a merge either).
//
// Result lives entirely on rank 0 afterward: every other rank's own DXADislocationLines is cleared
// to empty. This composes for free with every existing consumer (smooth_dxa_dislocation_lines,
// write_dxa_ca_file, write_ovito_interface_mesh's sibling write_dxa_dislocation_lines) since they
// already either operate per-line (a no-op on an empty rank) or already gather-and-concatenate every
// rank's own contribution themselves -- rank 0 contributing everything and every other rank
// contributing nothing produces the correct total with zero changes needed on their end. Run this
// operator right after compute_dxa_circuit_sweep and before smooth_dxa_dislocation_lines, so
// coarsening/smoothing sees each real dislocation's full, stitched length at once rather than
// smoothing two disjoint fragments independently.
namespace exaStamp
{
  using namespace exanb;

  class ComputeDXAMPIStitchLines : public OperatorNode
  {
    ADD_SLOT( MPI_Comm             , mpi                  , INPUT , REQUIRED );
    ADD_SLOT( Domain                , domain               , INPUT , REQUIRED , DocString{"Only used for periodic-image-aware wrapping when averaging a shared ghost atom's own position (the same real atom can be recorded at a different periodic image on each rank)"} );
    ADD_SLOT( DXADislocationLines  , dxa_dislocation_lines , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( double               , match_tolerance       , INPUT , 5.0 , DocString{"Max distance (same units as line_positions) between two ranks' own open-boundary endpoints for them to be considered the same real crossing point -- both ranks trace the identical dislocation core through the same ghost-shared atoms there, so a real match should coincide much more tightly than this; kept generous to tolerate ordinary circuit-centroid noise without needing exact agreement"} );
    ADD_SLOT( double               , burgers_tolerance     , INPUT , 0.1 , DocString{"Max |difference| (or |sum|, to allow an opposite-sign convention from the other rank's own sweep direction) between two candidate ends' Burgers vectors for them to be considered compatible"} );
    ADD_SLOT( long                 , n_mpi_stitches        , OUTPUT , DocString{"Number of cross-rank end-to-end matches actually performed"} );

  public:
    inline void execute () override final
    {
      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      if( np <= 1 )
      {
        *n_mpi_stitches = 0;
        return; // nothing to stitch -- single rank already has everything it will ever have
      }

      DXADislocationLines& dl = *dxa_dislocation_lines;
      const size_t n_local = dl.line_positions.size();

      // Serialize this rank's own lines into one flat, self-delimiting double buffer: per line,
      // burgers.x/y/z, open_front, open_back, is_loop, npoints, then npoints*4 doubles (x,y,z,
      // core_size) -- same style as write_dxa_ca_file's own serialization, extended with the 2 new
      // open_boundary flags and is_loop (needed here, unlike write_dxa_ca_file, to correctly carry a
      // standalone line's own pre-existing self-closure through unchanged) -- then, for each end,
      // n_atoms followed by n_atoms*4 doubles (atom_id as a double -- exact for any realistic particle
      // count, well within a double's 53-bit mantissa -- then x,y,z).
      std::vector<double> local_buf;
      for(size_t i=0;i<n_local;i++)
      {
        const Vec3d& b = dl.burgers_vector[i];
        local_buf.push_back(b.x); local_buf.push_back(b.y); local_buf.push_back(b.z);
        local_buf.push_back( dl.open_boundary_front.empty() ? 0. : static_cast<double>( dl.open_boundary_front[i] ) );
        local_buf.push_back( dl.open_boundary_back.empty()  ? 0. : static_cast<double>( dl.open_boundary_back[i] ) );
        local_buf.push_back( dl.is_loop.empty() ? 0. : static_cast<double>( dl.is_loop[i] ) );
        const auto& pts = dl.line_positions[i];
        local_buf.push_back( static_cast<double>( pts.size() ) );
        const bool has_core = !dl.core_size.empty() && dl.core_size[i].size() == pts.size();
        for(size_t k=0;k<pts.size();k++)
        {
          local_buf.push_back(pts[k].x); local_buf.push_back(pts[k].y); local_buf.push_back(pts[k].z);
          local_buf.push_back( has_core ? static_cast<double>( dl.core_size[i][k] ) : 0. );
        }
        auto write_boundary_atoms = [&]( const std::vector<std::vector<uint64_t>>& id_lists, const std::vector<std::vector<Vec3d>>& pos_lists )
        {
          const bool have = !id_lists.empty() && !pos_lists.empty() && pos_lists[i].size() == id_lists[i].size();
          const size_t n_atoms = have ? id_lists[i].size() : 0;
          local_buf.push_back( static_cast<double>( n_atoms ) );
          for(size_t k=0;k<n_atoms;k++)
          {
            local_buf.push_back( static_cast<double>( id_lists[i][k] ) );
            local_buf.push_back( pos_lists[i][k].x ); local_buf.push_back( pos_lists[i][k].y ); local_buf.push_back( pos_lists[i][k].z );
          }
        };
        write_boundary_atoms( dl.boundary_loop_atom_id_front, dl.boundary_loop_atom_pos_front );
        write_boundary_atoms( dl.boundary_loop_atom_id_back, dl.boundary_loop_atom_pos_back );
      }
      const int local_count = static_cast<int>( local_buf.size() );

      std::vector<int> recvcounts, displs;
      if( rank == 0 ) { recvcounts.resize(np); }
      MPI_Gather( &local_count, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, *mpi );

      std::vector<double> all_buf;
      if( rank == 0 )
      {
        displs.resize(np);
        int off = 0;
        for(int r=0;r<np;r++) { displs[r] = off; off += recvcounts[r]; }
        all_buf.resize(off);
      }
      MPI_Gatherv( local_buf.data(), local_count, MPI_DOUBLE,
                   rank==0 ? all_buf.data() : nullptr,
                   rank==0 ? recvcounts.data() : nullptr,
                   rank==0 ? displs.data() : nullptr,
                   MPI_DOUBLE, 0, *mpi );

      // Every rank's own DXADislocationLines is empty afterward except rank 0's -- every existing
      // consumer either operates per-line (no-op on empty) or already gathers+concatenates every
      // rank's own contribution itself, so this composes for free (see this file's own header).
      if( rank != 0 )
      {
        dl.lines.clear(); dl.line_positions.clear(); dl.core_size.clear(); dl.is_loop.clear();
        dl.junction_vertices.clear(); dl.burgers_vector.clear(); dl.dislocation_id.clear();
        dl.open_boundary_front.clear(); dl.open_boundary_back.clear();
        dl.boundary_loop_atom_id_front.clear(); dl.boundary_loop_atom_id_back.clear();
        dl.boundary_loop_atom_pos_front.clear(); dl.boundary_loop_atom_pos_back.clear();
        *n_mpi_stitches = 0;
        return;
      }

      struct BoundaryAtom { uint64_t id; Vec3d pos; };
      struct DecodedLine { Vec3d burgers; std::vector<Vec3d> pts; std::vector<int32_t> core; bool open_front=false, open_back=false, is_loop=false; std::vector<BoundaryAtom> front_atoms, back_atoms; };
      std::vector<DecodedLine> all_lines;
      for(int r=0;r<np;r++)
      {
        size_t pos = static_cast<size_t>( displs[r] );
        const size_t end = pos + static_cast<size_t>( recvcounts[r] );
        while( pos < end )
        {
          DecodedLine line;
          line.burgers = Vec3d{ all_buf[pos], all_buf[pos+1], all_buf[pos+2] };
          line.open_front = all_buf[pos+3] != 0.;
          line.open_back  = all_buf[pos+4] != 0.;
          line.is_loop    = all_buf[pos+5] != 0.;
          pos += 6;
          const size_t npoints = static_cast<size_t>( all_buf[pos] );
          pos += 1;
          line.pts.reserve(npoints); line.core.reserve(npoints);
          for(size_t k=0;k<npoints;k++)
          {
            line.pts.push_back( Vec3d{ all_buf[pos], all_buf[pos+1], all_buf[pos+2] } );
            line.core.push_back( static_cast<int32_t>( all_buf[pos+3] ) );
            pos += 4;
          }
          for( std::vector<BoundaryAtom>* atoms : { &line.front_atoms, &line.back_atoms } )
          {
            const size_t n_atoms = static_cast<size_t>( all_buf[pos] );
            pos += 1;
            atoms->reserve(n_atoms);
            for(size_t k=0;k<n_atoms;k++)
            {
              atoms->push_back( BoundaryAtom{ static_cast<uint64_t>( all_buf[pos] ), Vec3d{ all_buf[pos+1], all_buf[pos+2], all_buf[pos+3] } } );
              pos += 4;
            }
          }
          all_lines.push_back( std::move(line) );
        }
      }
      const int N = static_cast<int>( all_lines.size() );

      // The SAME real (ghost-shared) atom can be recorded at a different periodic image on each
      // rank's own side of a match (e.g. one rank's own ghost halo wraps it at z=110, the other's at
      // z=0 -- the same real location under periodic boundary conditions, ~110 Å apart in raw
      // Cartesian coordinates if never wrapped). Averaging two such positions naively would land the
      // snap point at a nonsensical mid-box position instead of the correct shared location. Wrap `p`
      // to whichever periodic image sits closest to `ref` (minimum-image convention, per axis) before
      // ever averaging two positions for the same atom id.
      // domain->bounds_size() is in the SAME raw/reduced (pre-xform) frame as field::rx/ry/rz --
      // NOT real space (see compute_delaunay.cpp's own header comment on this, a real, previously-
      // unapplied xform bug found and fixed the same day this file's own wrap logic was checked).
      // mesh.vertices (hence every Vec3d position this operator handles) IS real space after that
      // fix, so the per-axis wrap size here must be the REAL edge length too -- same xform-applied
      // edge-vector convention write_dxa_ca_file's own SIMULATION_CELL_MATRIX already uses. This
      // correctly handles a diagonal (scaled) xform; a genuinely sheared (non-diagonal) box would
      // need full lattice-vector minimum-image search, not just a per-axis wrap -- not implemented,
      // since every real xform in this project so far is identity, but flagged honestly rather than
      // silently assumed away.
      const Mat3d xform = domain->xform();
      const Vec3d reduced_size = domain->bounds_size();
      const Vec3d box_size { norm( xform * Vec3d{reduced_size.x,0.,0.} ),
                             norm( xform * Vec3d{0.,reduced_size.y,0.} ),
                             norm( xform * Vec3d{0.,0.,reduced_size.z} ) };
      const bool periodic_x = domain->periodic_boundary_x();
      const bool periodic_y = domain->periodic_boundary_y();
      const bool periodic_z = domain->periodic_boundary_z();
      auto wrap_axis = [&]( double d, double box, bool periodic ) -> double
      {
        if( !periodic || box <= 0. ) { return d; }
        while( d >  0.5*box ) { d -= box; }
        while( d < -0.5*box ) { d += box; }
        return d;
      };
      auto wrap_to_reference = [&]( const Vec3d& p, const Vec3d& ref ) -> Vec3d
      {
        const Vec3d d{ wrap_axis( p.x-ref.x, box_size.x, periodic_x ),
                       wrap_axis( p.y-ref.y, box_size.y, periodic_y ),
                       wrap_axis( p.z-ref.z, box_size.z, periodic_z ) };
        return ref + d;
      };

      // Match candidates: every open_boundary-flagged end. end_partner[li][0/1] = the matched
      // (other_li, other_end, snap_pos) once matched, (-1,-1,{}) otherwise.
      const double pos_tol = *match_tolerance;
      const double b_tol = *burgers_tolerance;
      struct EndPartner { int li=-1; int end=-1; Vec3d snap_pos{0.,0.,0.}; };
      std::vector<std::array<EndPartner,2>> end_partner( N );

      struct Candidate { int li; int end; }; // end: 0=front, 1=back
      std::vector<Candidate> candidates;
      for(int li=0; li<N; li++)
      {
        if( all_lines[li].open_front ) { candidates.push_back({li,0}); }
        if( all_lines[li].open_back )  { candidates.push_back({li,1}); }
      }
      auto end_atoms = [&]( int li, int end ) -> const std::vector<BoundaryAtom>& { return end==0 ? all_lines[li].front_atoms : all_lines[li].back_atoms; };
      auto end_pos   = [&]( int li, int end ) -> const Vec3d& { return end==0 ? all_lines[li].pts.front() : all_lines[li].pts.back(); };
      {
        size_t n_with_atoms = 0;
        for( const auto& c : candidates ) { if( !end_atoms(c.li,c.end).empty() ) { ++n_with_atoms; } }
        lout << "compute_dxa_mpi_stitch_lines: " << candidates.size() << " open-boundary end candidates total ("
             << n_with_atoms << " carry boundary-loop atom data)" << std::endl;
      }

      // Primary pass: EXACT match via a shared ghost atom id -- but a real 3+-way junction can sit
      // exactly at (or straddle) an MPI domain boundary just as easily as a clean 2-way pass-through
      // can, and every arm meeting there independently gets its own OpenEdge stop, potentially all
      // touching some of the SAME boundary atoms. A first version of this matching greedily paired
      // each candidate with whichever OTHER candidate shared the most atoms with it -- with 3+ arms
      // present, this silently spliced together two arbitrary (whichever the greedy pass happened to
      // pick, itself sensitive to mesh non-determinism) arms of what is really a junction, as if they
      // were one continuous line -- the actual cause of the user-reported "Burgers vector changes"
      // and "still doesn't work at MPI boundaries" bugs: the reported vector, and the resulting
      // topology, depended on an arbitrary/wrong pairing among genuinely different dislocations, not
      // on a real defect in the seam. Fixed by grouping candidates that share atoms into connected
      // components (union-find, mirroring the same real-vs-junction distinction compute_dxa_
      // circuit_sweep's own do_merge already makes for same-rank stops) BEFORE deciding anything:
      // a component of exactly 2 members is a clean pass-through (stitch them); a component of 3+ is
      // a real junction straddling the boundary (leave every member standalone, exactly like a
      // same-rank junction stop already does -- no attempt at reconstructing shared junction-node
      // geometry across ranks here, out of scope for this fix).
      //
      // In a dense network (the quadrupole test, unlike the screw dipole's 2 well-separated lines),
      // uniting on ANY single shared atom is too permissive: two genuinely DIFFERENT dislocations
      // passing close together near the same boundary region can each have a boundary loop that
      // happens to touch one of the same nearby atoms, without being the same crossing at all --
      // measured as several truly-independent candidates all getting lumped into one fake "3+-way
      // junction" and left unstitched, when the real structure is 2-3 separate clean crossings sitting
      // near each other. A real single crossing's two independently-recorded loops share SEVERAL
      // atoms (3, in every verified screw-dipole match) since both sides trace the same core
      // cross-section; require at least MIN_SHARED_ATOMS_FOR_UNION before treating two candidates as
      // connected at all, so incidental single-atom proximity between unrelated dislocations no longer
      // creates a spurious junction.
      static constexpr int MIN_SHARED_ATOMS_FOR_UNION = 2;
      std::unordered_map<uint64_t, std::vector<int>> atom_to_candidates; // atom id -> candidate indices
      for(size_t c=0;c<candidates.size();c++)
      {
        for( const auto& a : end_atoms( candidates[c].li, candidates[c].end ) ) { atom_to_candidates[a.id].push_back( static_cast<int>(c) ); }
      }
      std::map<std::pair<int,int>,int> pair_shared_count;
      for( const auto& kv : atom_to_candidates )
      {
        const auto& cs = kv.second;
        for(size_t i=0;i<cs.size();i++) { for(size_t j=i+1;j<cs.size();j++)
        {
          int a = cs[i], b = cs[j]; if( a > b ) { std::swap(a,b); }
          ++pair_shared_count[{a,b}];
        } }
      }
      std::vector<int> uf( candidates.size() );
      for(size_t c=0;c<candidates.size();c++) { uf[c] = static_cast<int>(c); }
      std::function<int(int)> find_root = [&]( int x ) -> int { while( uf[x] != x ) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; };
      auto unite = [&]( int a, int b ) { a = find_root(a); b = find_root(b); if( a != b ) { uf[a] = b; } };
      for( const auto& [pr, cnt] : pair_shared_count )
      {
        if( cnt >= MIN_SHARED_ATOMS_FOR_UNION ) { unite( pr.first, pr.second ); }
      }
      std::unordered_map<int,std::vector<int>> components;
      for(size_t c=0;c<candidates.size();c++) { components[ find_root(static_cast<int>(c)) ].push_back( static_cast<int>(c) ); }

      long n_stitches = 0, n_exact = 0, n_fallback = 0, n_junctions_skipped = 0;
      std::vector<bool> candidate_used( candidates.size(), false );

      auto pair_shared = [&]( int a, int b ) -> int
      {
        int x=a, y=b; if( x>y ) { std::swap(x,y); }
        const auto it = pair_shared_count.find({x,y});
        return it != pair_shared_count.end() ? it->second : 0;
      };

      // Tries to stitch candidates a/b: Burgers-compatibility check (direction-agnostic -- see the
      // comment below for why), snap-position averaging, end_partner bookkeeping. Returns false
      // (nothing done, no candidate marked used) if the pair fails the compatibility check -- the
      // caller decides what "not stitched" means in its own context (leave standalone, or fall
      // through to try a different pairing).
      auto try_stitch_pair = [&]( int a, int b ) -> bool
      {
        const int li_a = candidates[a].li, end_a = candidates[a].end;
        const int li_b = candidates[b].li, end_b = candidates[b].end;
        if( li_a == li_b ) { return false; } // a line's own two ends sharing an atom with only each other -- not a real cross-rank stitch

        // Shared atoms are necessary but NOT sufficient evidence this is a genuine pass-through: two
        // physically different dislocation lines can pass close enough together (e.g. near an
        // otherwise-undetected junction, or just densely-packed nearby defects) that their own
        // boundary loops happen to touch a couple of the same atoms without actually being the same
        // line. Cross-check with the Burgers vector -- but DIRECTION-AGNOSTIC (either equal or
        // exactly opposite passes), not end-topology-aware. An earlier version of this check tried
        // to PREDICT which relation should hold from end topology (front/back), on the theory that
        // "front-back" is a natural non-reversed continuation (raw vectors equal) while "front-front/
        // back-back" needs one side reversed (raw vectors opposite) -- measured, directly, to be
        // wrong: on the screw-dipole test, the SAME physical pair of fragments showed up in one run
        // as "same-sense" with opposite vectors (rejected, "expected opposite" -- so this one actually
        // matched by luck) and in another run as "opposite-sense" with the SAME (still opposite-to-
        // each-other) vectors (rejected, "expected equal" -- wrongly, this time). Root cause: a
        // fragment's own "front"/"back" label is just which end of its own points array happens to be
        // index 0 vs last -- an arbitrary artifact of which direction its own independent circuit
        // trace happened to grow, with no fixed relationship to its own stored Burgers vector's sign
        // (see burgers_of_loop's own doc comment in compute_dxa_circuit_sweep.cpp: the loop's own
        // internal array order determines the sign, and that order is itself arbitrary per trace).
        // Two different <111>/2 families are neither equal nor exact negatives of each other, so
        // accepting either relation still correctly rejects a genuine junction/incompatible pair --
        // it just stops trying to guess a fixed relation from a label that carries no such guarantee.
        const Vec3d& b_a = all_lines[li_a].burgers;
        const Vec3d& b_b = all_lines[li_b].burgers;
        const bool compatible = ( norm( b_a - b_b ) < b_tol ) || ( norm( b_a + b_b ) < b_tol );
        if( !compatible )
        {
          lout << "compute_dxa_mpi_stitch_lines: REJECTED line " << li_a << " end " << end_a << " <-> line "
               << li_b << " end " << end_b << ": shares ghost atoms but Burgers vectors ("
               << b_a.x << " " << b_a.y << " " << b_a.z << " vs " << b_b.x << " " << b_b.y << " " << b_b.z
               << ") aren't physically compatible (neither equal nor opposite) -- likely 2 distinct "
               << "nearby dislocations, not one continuous line; leaving both standalone" << std::endl;
          return false;
        }

        // Snap position: average every atom the two ends actually have in common (not just the ones
        // that happened to trigger the union -- a real crossing usually shares several), wrapping
        // side b's own recorded position to whichever periodic image sits closest to side a's before
        // averaging -- the same real (ghost-shared) atom can be recorded at a different periodic
        // image on each rank's own side (see this section's own header comment above).
        std::unordered_map<uint64_t,Vec3d> atoms_a; for( const auto& at : end_atoms(li_a,end_a) ) { atoms_a[at.id] = at.pos; }
        Vec3d snap{0.,0.,0.}; int n_shared = 0;
        for( const auto& at : end_atoms(li_b,end_b) )
        {
          const auto it = atoms_a.find( at.id );
          if( it == atoms_a.end() ) { continue; }
          snap = snap + it->second + wrap_to_reference( at.pos, it->second );
          ++n_shared;
        }
        if( n_shared == 0 ) { return false; } // shouldn't happen (union implies at least one shared atom), guard anyway
        snap = snap / static_cast<double>( 2 * n_shared );

        candidate_used[static_cast<size_t>(a)] = true; candidate_used[static_cast<size_t>(b)] = true;
        end_partner[li_a][end_a] = { li_b, end_b, snap };
        end_partner[li_b][end_b] = { li_a, end_a, snap };
        ++n_stitches; ++n_exact;
        lout << "compute_dxa_mpi_stitch_lines: EXACT match line " << li_a << " end " << end_a
             << " <-> line " << li_b << " end " << end_b << " (" << n_shared << " shared ghost atoms), "
             << "pre-snap centroid offset " << ( norm( end_pos(li_a,end_a)-snap ) + norm( end_pos(li_b,end_b)-snap ) )
             << " Ang (each side's own raw circuit-centroid distance from the shared atom -- expected, "
             << "not the final seam gap: both sides get overwritten to the exact same shared coordinate)" << std::endl;
        return true;
      };

      for( const auto& kv : components )
      {
        std::vector<int> remaining = kv.second;
        if( remaining.size() < 2 ) { continue; } // no other candidate anywhere shares an atom with this one

        // A domain-decomposition corner/edge (more than 2 MPI ranks sharing a boundary there) can
        // make a single real 2-way crossing show up as a 3+-member group: every rank near that corner
        // independently records an open-boundary end, and several of them end up sharing >=2 atoms
        // purely from being spatially close together, not because they're distinct dislocation arms.
        // Found via the rectangular-loop test (provably ONE physical dislocation, so any 3+ group
        // there is unambiguously this artifact, not a real junction): 3 members, but one pair shared
        // 3 atoms (and landed at bit-identical snap positions) while the other two pairings each only
        // shared 2 -- a clear, dominant real match with a redundant third recording nearby, not 3
        // genuinely competing arms. Before declaring "real junction, leave everyone unstitched", peel
        // off a pair ONLY if it's unambiguously dominant (strictly more shared atoms than any OTHER
        // pairing touching either of its own two members -- a real junction's arms don't have this
        // property, since no single pairing there is privileged over the others) AND Burgers-
        // compatible; repeat until <=2 members remain. A genuine 3+-way junction (no dominant pairing
        // to peel, e.g. the quadrupole's real cases) falls straight through unchanged.
        while( remaining.size() > 2 )
        {
          int best_a=-1, best_b=-1, best_cnt=0;
          for(size_t i=0;i<remaining.size();i++) { for(size_t j=i+1;j<remaining.size();j++)
          {
            const int cnt = pair_shared( remaining[i], remaining[j] );
            if( cnt > best_cnt ) { best_cnt = cnt; best_a = remaining[i]; best_b = remaining[j]; }
          } }
          if( best_a < 0 ) { break; }

          bool dominant = true;
          for( int m : remaining )
          {
            if( m == best_a || m == best_b ) { continue; }
            if( pair_shared(best_a,m) >= best_cnt || pair_shared(best_b,m) >= best_cnt ) { dominant = false; break; }
          }
          if( !dominant ) { break; }

          lout << "compute_dxa_mpi_stitch_lines: dominant pair (" << best_cnt << " shared atoms) found within a "
               << remaining.size() << "-member group at an MPI boundary -- treating as one real crossing plus "
               << (remaining.size()-2) << " redundant nearby recording(s), not a junction" << std::endl;
          if( !try_stitch_pair( best_a, best_b ) ) { break; } // dominant pairing itself isn't Burgers-compatible -- give up peeling, treat as a real junction below

          remaining.erase( std::remove( remaining.begin(), remaining.end(), best_a ), remaining.end() );
          remaining.erase( std::remove( remaining.begin(), remaining.end(), best_b ), remaining.end() );
        }

        if( remaining.size() > 2 )
        {
          // real 3+-way junction straddling the boundary -- leave every arm standalone, don't guess
          // which two (of 3+) actually continue one another.
          ++n_junctions_skipped;
          lout << "compute_dxa_mpi_stitch_lines: 3+-way junction detected at an MPI boundary ("
               << remaining.size() << " arms sharing ghost atoms there) -- leaving all "
               << remaining.size() << " arms standalone, not stitching" << std::endl;
          for( int c : remaining ) { candidate_used[static_cast<size_t>(c)] = true; } // exclude from the fallback pass too
          continue;
        }
        if( remaining.size() == 2 ) { try_stitch_pair( remaining[0], remaining[1] ); }
        // remaining.size() == 1: nothing left to pair here -- falls through to the fallback pass below.
      }

      // Fallback pass: fuzzy position+Burgers proximity, only for candidates that found no exact
      // atom-id match above (shouldn't normally trigger, kept as a safety net -- see this file's own
      // header comment).
      for(size_t a=0; a<candidates.size(); a++)
      {
        if( candidate_used[a] ) { continue; }
        const int li_a = candidates[a].li, end_a = candidates[a].end;
        const Vec3d& pos_a = end_pos(li_a,end_a);
        const Vec3d& b_a = all_lines[li_a].burgers;

        int best_b = -1; double best_d = pos_tol;
        for(size_t b=a+1; b<candidates.size(); b++)
        {
          if( candidate_used[b] || candidates[b].li == li_a ) { continue; }
          const int li_b = candidates[b].li, end_b = candidates[b].end;
          const Vec3d& pos_b = end_pos(li_b,end_b);
          const double d = norm( pos_a - pos_b );
          if( d >= best_d ) { continue; }
          const Vec3d& b_b = all_lines[li_b].burgers;
          // Direction-agnostic, same as the exact-match pass above -- see that pass's own comment
          // for why end-topology (front/back) can't reliably predict which relation should hold.
          const bool burgers_ok = ( norm(b_a - b_b) < b_tol ) || ( norm(b_a + b_b) < b_tol );
          if( !burgers_ok ) { continue; }
          best_d = d; best_b = static_cast<int>(b);
        }

        if( best_b >= 0 )
        {
          candidate_used[a] = true; candidate_used[static_cast<size_t>(best_b)] = true;
          const int li_b = candidates[best_b].li, end_b = candidates[best_b].end;
          const Vec3d snap = ( pos_a + end_pos(li_b,end_b) ) / 2.;
          end_partner[li_a][end_a] = { li_b, end_b, snap };
          end_partner[li_b][end_b] = { li_a, end_a, snap };
          ++n_stitches; ++n_fallback;
          lout << "compute_dxa_mpi_stitch_lines: FALLBACK (no shared ghost atom) match line " << li_a << " end " << end_a
               << " <-> line " << li_b << " end " << end_b << ", gap " << best_d << " Ang" << std::endl;
        }
      }
      // Walk the match graph (max degree 2 per line -> simple paths + simple cycles) and assemble
      // the final, stitched line set. Paths first (lines with at most 1 matched end -- a genuine
      // chain endpoint, or a standalone line untouched by any match), then whatever remains must be
      // pure cycles.
      //
      // Seam handling: OVERWRITE (not append) each fragment's own connecting point with the matched
      // `snap_pos`, and DROP every non-first fragment's own leading point in the walk (it's the exact
      // same shared location the previous fragment's own trailing point was just overwritten to) --
      // this is safe HERE, unlike the earlier same-rank junction-endpoint bug (where overwriting used
      // a merge-chain-resolved, unverified position): the ghost atom id match is a real, direct
      // physical identity check, not a guess, so replacing either side's own approximate loop-centroid
      // with the actual shared coordinate is strictly more correct, not less. Net effect: the two
      // fragments share the literal same coordinate at the join (bit-identical after the overwrite),
      // not two independently-computed points a few Å apart.
      auto ordered_points = [&]( const DecodedLine& L, bool forward ) -> std::pair<std::vector<Vec3d>,std::vector<int32_t>>
      {
        if( forward ) { return { L.pts, L.core }; }
        return { std::vector<Vec3d>( L.pts.rbegin(), L.pts.rend() ), std::vector<int32_t>( L.core.rbegin(), L.core.rend() ) };
      };

      std::vector<bool> visited( N, false );
      struct AssembledLine { Vec3d burgers; std::vector<Vec3d> pts; std::vector<int32_t> core; bool is_loop; };
      std::vector<AssembledLine> assembled;

      for(int li=0; li<N; li++)
      {
        if( visited[li] ) { continue; }
        const bool front_matched = end_partner[li][0].li != -1;
        const bool back_matched  = end_partner[li][1].li != -1;
        if( front_matched && back_matched ) { continue; } // part of a cycle -- second pass below

        AssembledLine chain; chain.is_loop = all_lines[li].is_loop;
        int cur = li;
        int enter_end = front_matched ? 1 : 0; // the line's own free end (0=front free, 1=back free)
        bool first_fragment = true;
        bool have_chain_burgers = false;
        for(;;)
        {
          visited[cur] = true;
          const bool forward = ( enter_end == 0 );
          auto ord = ordered_points( all_lines[cur], forward );
          const int exit_end = forward ? 1 : 0;
          const auto& partner = end_partner[cur][exit_end];
          if( partner.li != -1 && !ord.first.empty() ) { ord.first.back() = partner.snap_pos; }
          const size_t skip = first_fragment ? 0 : 1;
          for(size_t k=skip;k<ord.first.size();k++) { chain.pts.push_back(ord.first[k]); chain.core.push_back(ord.second[k]); }

          // A fragment's own recorded Burgers vector is tied to ITS OWN original tangent direction
          // (front-to-back); if this walk emits it reversed (back-to-front), the vector must be
          // negated to stay physically consistent with the assembled chain's own overall tangent --
          // otherwise the reported vector for a stitched line depends arbitrarily on which raw
          // fragment happened to become the walk's starting point (itself a function of unrelated
          // mesh/seed non-determinism), which is exactly the "Burgers vector changes" bug the user
          // caught: it wasn't OVITO changing anything, this operator was reporting an inconsistent,
          // occasionally sign-flipped value for the very same physical dislocation from run to run.
          const Vec3d& raw_b = all_lines[cur].burgers;
          const Vec3d eff_b = forward ? raw_b : Vec3d{ -raw_b.x, -raw_b.y, -raw_b.z };
          if( !have_chain_burgers ) { chain.burgers = eff_b; have_chain_burgers = true; }
          else if( norm( eff_b - chain.burgers ) > 0.2 && norm( eff_b + chain.burgers ) > 0.2 )
          {
            lout << "compute_dxa_mpi_stitch_lines: WARNING fragment " << cur << " own (direction-corrected) "
                 << "Burgers vector " << eff_b.x << " " << eff_b.y << " " << eff_b.z << " disagrees with chain's "
                 << chain.burgers.x << " " << chain.burgers.y << " " << chain.burgers.z
                 << " -- possible bad stitch, reporting the chain's own starting value regardless" << std::endl;
          }

          if( partner.li == -1 ) { break; }
          cur = partner.li; enter_end = partner.end; first_fragment = false;
        }
        if( chain.pts.size() > all_lines[li].pts.size() ) { chain.is_loop = false; } // a real multi-fragment path, not a pre-existing self-closed loop
        assembled.push_back( std::move(chain) );
      }
      for(int li=0; li<N; li++)
      {
        if( visited[li] ) { continue; }
        // pure cycle -- arbitrary cut at this line's own front, walk until the ring closes back here
        AssembledLine chain; chain.is_loop = true;
        int cur = li; int enter_end = 0; bool closed = false; bool first_fragment = true;
        bool have_chain_burgers = false;
        for(;;)
        {
          visited[cur] = true;
          const bool forward = ( enter_end == 0 );
          auto ord = ordered_points( all_lines[cur], forward );
          const int exit_end = forward ? 1 : 0;
          const auto& partner = end_partner[cur][exit_end];
          if( partner.li != -1 && !ord.first.empty() ) { ord.first.back() = partner.snap_pos; }
          const size_t skip = first_fragment ? 0 : 1;
          for(size_t k=skip;k<ord.first.size();k++) { chain.pts.push_back(ord.first[k]); chain.core.push_back(ord.second[k]); }

          // Same Burgers-vector direction correction as the path-walk above -- see its own comment.
          const Vec3d& raw_b = all_lines[cur].burgers;
          const Vec3d eff_b = forward ? raw_b : Vec3d{ -raw_b.x, -raw_b.y, -raw_b.z };
          if( !have_chain_burgers ) { chain.burgers = eff_b; have_chain_burgers = true; }
          else if( norm( eff_b - chain.burgers ) > 0.2 && norm( eff_b + chain.burgers ) > 0.2 )
          {
            lout << "compute_dxa_mpi_stitch_lines: WARNING fragment " << cur << " own (direction-corrected) "
                 << "Burgers vector " << eff_b.x << " " << eff_b.y << " " << eff_b.z << " disagrees with chain's "
                 << chain.burgers.x << " " << chain.burgers.y << " " << chain.burgers.z
                 << " -- possible bad stitch, reporting the chain's own starting value regardless" << std::endl;
          }

          if( partner.li == -1 ) { break; } // shouldn't happen in a genuine cycle, guard anyway
          if( partner.li == li ) { closed = true; break; }
          cur = partner.li; enter_end = partner.end; first_fragment = false;
        }
        chain.is_loop = closed;
        assembled.push_back( std::move(chain) );
      }

      // Drop short lines whose ENTIRE path is redundant with (retraces the same physical stretch
      // as) a longer line -- found via the rectangular-loop test (a single physical closed
      // dislocation by construction): with keep_ghost_tets, a rank whose own local view legitimately
      // OWNS a small sliver of territory near a domain-decomposition corner can still independently
      // grow its own short segment duplicating a stretch a NEIGHBORING rank's own longer fragment
      // already covers -- both seeds are genuinely owned territory, so the seed-ownership
      // restriction in compute_dxa_circuit_sweep (try_seed_from only from an owned vertex) doesn't
      // prevent this specific case. Surfaces as a visually "duplicated"/loose-end short line right
      // where two ranks' own coverage overlaps. Treat a short line as redundant, not a real second
      // dislocation, if essentially every one of its own points sits within ordinary in-line point
      // spacing of SOME point on a strictly longer line -- a genuinely separate, real dislocation
      // nearby would not have this property (its own points would trace a physically distinct path).
      {
        static constexpr double REDUNDANT_POINT_TOL = 5.0; // Ang -- same scale as match_tolerance's own default
        static constexpr double REDUNDANT_COVERAGE_FRACTION = 0.9;
        std::vector<bool> drop( assembled.size(), false );
        for(size_t i=0;i<assembled.size();i++)
        {
          if( assembled[i].pts.empty() ) { continue; }
          for(size_t j=0;j<assembled.size();j++)
          {
            if( i==j || drop[j] || assembled[i].pts.size() >= assembled[j].pts.size() ) { continue; }
            size_t n_covered = 0;
            for( const auto& p : assembled[i].pts )
            {
              for( const auto& q : assembled[j].pts ) { if( norm(p-q) < REDUNDANT_POINT_TOL ) { ++n_covered; break; } }
            }
            if( static_cast<double>(n_covered) / static_cast<double>(assembled[i].pts.size()) >= REDUNDANT_COVERAGE_FRACTION )
            {
              drop[i] = true;
              lout << "compute_dxa_mpi_stitch_lines: dropping line " << i << " (" << assembled[i].pts.size()
                   << " points) as redundant with line " << j << " (" << assembled[j].pts.size() << " points) -- "
                   << n_covered << "/" << assembled[i].pts.size() << " of its own points sit within "
                   << REDUNDANT_POINT_TOL << " Ang of the longer line's own path, likely a duplicate ghost-overlap "
                   << "recording near an MPI boundary, not a separate physical dislocation" << std::endl;
              break;
            }
          }
        }
        if( std::any_of( drop.begin(), drop.end(), []( bool b ){ return b; } ) )
        {
          std::vector<AssembledLine> kept;
          for(size_t i=0;i<assembled.size();i++) { if( !drop[i] ) { kept.push_back( std::move(assembled[i]) ); } }
          assembled = std::move(kept);
        }
      }

      dl.lines.clear();
      dl.line_positions.clear();
      dl.core_size.clear();
      dl.is_loop.clear();
      dl.junction_vertices.clear();
      dl.burgers_vector.clear();
      dl.dislocation_id.clear();
      dl.open_boundary_front.clear();
      dl.open_boundary_back.clear();
      dl.boundary_loop_atom_id_front.clear();
      dl.boundary_loop_atom_id_back.clear();
      dl.boundary_loop_atom_pos_front.clear();
      dl.boundary_loop_atom_pos_back.clear();
      for(size_t i=0;i<assembled.size();i++)
      {
        dl.lines.push_back( {} );
        dl.line_positions.push_back( std::move( assembled[i].pts ) );
        dl.core_size.push_back( std::move( assembled[i].core ) );
        dl.is_loop.push_back( assembled[i].is_loop ? 1 : 0 );
        dl.burgers_vector.push_back( assembled[i].burgers );
        dl.dislocation_id.push_back( static_cast<int32_t>(i) );
        dl.open_boundary_front.push_back(0);
        dl.open_boundary_back.push_back(0);
        dl.boundary_loop_atom_id_front.push_back( {} );
        dl.boundary_loop_atom_id_back.push_back( {} );
        dl.boundary_loop_atom_pos_front.push_back( {} );
        dl.boundary_loop_atom_pos_back.push_back( {} );
      }

      *n_mpi_stitches = n_stitches;
      lout << "compute_dxa_mpi_stitch_lines: " << N << " raw fragments across " << np << " ranks, "
           << n_stitches << " cross-rank end-to-end matches (" << n_exact << " exact ghost-atom, "
           << n_fallback << " fuzzy fallback), " << n_junctions_skipped << " boundary junctions left "
           << "unstitched, " << assembled.size() << " lines after stitching" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Stitches dislocation lines that compute_dxa_circuit_sweep left fragmented at an MPI domain-
decomposition boundary back into single, continuous lines -- see this file's own header comment for
the full mechanism. Run right after compute_dxa_circuit_sweep, before smooth_dxa_dislocation_lines.
A no-op when running on a single MPI rank.

Usage example:

compute_dxa_circuit_sweep: { min_burgers_norm: 0.3 }
compute_dxa_mpi_stitch_lines: {}
smooth_dxa_dislocation_lines: { target_point_interval: 2.5, target_smoothing_level: 1 }

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(compute_dxa_mpi_stitch_lines)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_dxa_mpi_stitch_lines", make_simple_operator< ComputeDXAMPIStitchLines > );
  }

}
