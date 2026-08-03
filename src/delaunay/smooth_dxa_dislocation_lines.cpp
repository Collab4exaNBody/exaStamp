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

#include <exaStamp/delaunay/dxa_dislocation_lines.h>

#include <vector>

// DXA line post-processing: reduces the tortuosity of raw swept dislocation lines
// (compute_dxa_circuit_sweep, one point per elementary sweep move -- typically ~1 interatomic
// spacing apart) via the exact same two-stage mechanism OVITO's real DislocationAnalysisModifier
// applies before display/export, ported from ovito/src/ovito/crystalanalysis/objects/
// DislocationNetwork.cpp's own coarsenDislocationLine()/smoothDislocationLine():
//
//   1. Coarsening (target_point_interval, OVITO's own "linePointInterval", default 2.5): merges
//      consecutive raw points into single averaged points. The merge group size adapts to
//      DXADislocationLines::core_size (the sweeping circuit's own vertex count when each point was
//      recorded) -- a point recorded while the circuit was wide/stretched (typically near a
//      junction, noisier) gets merged more aggressively than one recorded while the circuit was at
//      its normal, narrow width. The two endpoints of an open segment are always kept fixed (so
//      junction connectivity between segments isn't disturbed); a closed loop instead gets one
//      "seam" point built from a combined start+end averaging window.
//   2. Smoothing (target_smoothing_level, OVITO's own "lineSmoothingLevel", default 1): a 2D Taubin
//      "signal processing" mesh-smoothing pass (Taubin, SIGGRAPH 95) on the now-coarsened points --
//      alternates a small forward (lambda) and backward (mu, |mu|>|lambda|) Laplacian relaxation
//      step, which smooths high-frequency zig-zag without the net shrinkage a naive single-direction
//      Laplacian smooth would introduce. Endpoints of an open segment are pinned (zero Laplacian);
//      a closed loop wraps around instead.
//
// Consumes and overwrites DXADislocationLines::line_positions/core_size in place (matching OVITO's
// own segment->line = std::move(line); segment->coreSize.clear();) -- run once, after
// compute_dxa_circuit_sweep, before any consumer (write_dxa_dislocation_lines, write_dxa_ca_file)
// that cares about the final line shape. Only touches lines with a non-empty core_size entry (i.e.
// compute_dxa_circuit_sweep's own output) -- the other, vertex-index-based extractors don't
// populate line_positions/core_size at all and are left untouched.
namespace exaStamp
{
  using namespace exanb;

  namespace
  {
    // Ported from DislocationNetwork::coarsenDislocationLine(). See this file's own header
    // comment for what `core_size` (OVITO's own coreSize) means and why it weights the merge.
    void coarsen_dislocation_line( double target_point_interval,
                                    const std::vector<Vec3d>& input, const std::vector<int32_t>& core_size,
                                    std::vector<Vec3d>& output, std::vector<int32_t>& output_core_size,
                                    bool is_closed_loop )
    {
      if( target_point_interval <= 0.0 || input.size() < 4 )
      {
        output = input;
        output_core_size = core_size;
        return;
      }

      output.clear();
      output_core_size.clear();

      // Always keep the endpoints of an open segment fixed so junction connectivity isn't
      // disturbed; a closed loop's own "seam" point is built below instead.
      if( !is_closed_loop )
      {
        output.push_back( input.front() );
        output_core_size.push_back( core_size.front() );
      }

      const int min_num_points = is_closed_loop ? 4 : 2;
      const int n = static_cast<int>( input.size() );

      size_t fwd_ptr = 0, fwd_core_ptr = 0;
      int sum = 0, count = 0;
      Vec3d com{0.0,0.0,0.0};

      // Average over a half interval, starting from the beginning of the segment.
      do
      {
        sum += core_size[fwd_core_ptr];
        com = com + ( input[fwd_ptr] - input.front() );
        count++;
        fwd_ptr++; fwd_core_ptr++;
      }
      while( 2*count*count < static_cast<int>(target_point_interval*sum) && count+1 < n/min_num_points/2 );

      // Average over a half interval, starting from the end of the segment.
      size_t bwd_ptr = input.size()-1, bwd_core_ptr = core_size.size()-1;
      while( count*count < static_cast<int>(target_point_interval*sum) && count < n/min_num_points )
      {
        sum += core_size[bwd_core_ptr];
        com = com + ( input[bwd_ptr] - input.back() );
        count++;
        bwd_ptr--; bwd_core_ptr--;
      }

      if( is_closed_loop )
      {
        output.push_back( input.front() + com / static_cast<double>(count) );
        output_core_size.push_back( sum/count );
      }

      // Middle of the segment, in full (not half) intervals.
      while( fwd_ptr < bwd_ptr )
      {
        int gsum = 0, gcount = 0;
        Vec3d gcom{0.0,0.0,0.0};
        do
        {
          gsum += core_size[fwd_core_ptr];
          gcom = gcom + input[fwd_ptr];
          gcount++;
          fwd_ptr++; fwd_core_ptr++;
        }
        while( gcount*gcount < static_cast<int>(target_point_interval*gsum) && gcount+1 < n/min_num_points && fwd_ptr != bwd_ptr );
        output.push_back( gcom / static_cast<double>(gcount) );
        output_core_size.push_back( gsum/gcount );
      }

      if( !is_closed_loop )
      {
        output.push_back( input.back() );
        output_core_size.push_back( core_size.back() );
      }
      else
      {
        // Deliberately reuses the OUTER sum/count/com from the two boundary half-interval passes
        // above (not the middle loop's own last group) -- this is what OVITO's own C++ scoping
        // does (the middle loop declares its OWN local sum/count/com, shadowing these), and it's
        // the mechanism that makes a loop's closing point land near its own opening point (both
        // built from the same combined start+end averaging window).
        output.push_back( input.back() + com / static_cast<double>(count) );
        output_core_size.push_back( sum/count );
      }
    }

    // Ported from DislocationNetwork::smoothDislocationLine() -- 2D Taubin mesh smoothing
    // (Taubin, "A Signal Processing Approach To Fair Surface Design", SIGGRAPH 95).
    void smooth_dislocation_line( int smoothing_level, std::vector<Vec3d>& line, bool is_loop )
    {
      if( smoothing_level <= 0 || line.size() <= 2 ) { return; }
      if( is_loop && line.size() <= 4 ) { return; } // don't smooth loops with too few points to mean anything

      const double k_PB = 0.1;
      const double lambda = 0.5;
      const double mu = 1.0 / (k_PB - 1.0/lambda);
      const double prefactors[2] = { lambda, mu };

      const size_t n = line.size();
      std::vector<Vec3d> laplacians(n);

      for(int iter=0; iter<smoothing_level; iter++)
      {
        for(int pass=0; pass<=1; pass++)
        {
          size_t li = 0;
          if( !is_loop ) { laplacians[li++] = Vec3d{0.0,0.0,0.0}; }
          else { laplacians[li++] = ( (line[n-2]-line[n-3]) + (line[1]-line[0]) ) * 0.5; }

          size_t p1 = 0, p2 = 1;
          for(;;)
          {
            const size_t p0 = p1;
            ++p1; ++p2;
            if( p2 == n ) { break; }
            laplacians[li++] = ( (line[p0]-line[p1]) + (line[p2]-line[p1]) ) * 0.5;
          }
          laplacians[li++] = laplacians[0];

          for(size_t k=0;k<n;k++) { line[k] = line[k] + laplacians[k]*prefactors[pass]; }
        }
      }
    }
  }

  class SmoothDXADislocationLines : public OperatorNode
  {
    ADD_SLOT( DXADislocationLines , dxa_dislocation_lines   , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( double              , target_point_interval   , INPUT , 2.5 , DocString{"Coarsening target -- OVITO's own linePointInterval. Larger values merge more raw sweep points into each output point (shorter, straighter final lines); 0 or negative disables coarsening entirely."} );
    ADD_SLOT( long                , target_smoothing_level  , INPUT , 1 , DocString{"Number of Taubin smoothing iterations applied after coarsening -- OVITO's own lineSmoothingLevel. 0 disables smoothing."} );

  public:
    inline void execute () override final
    {
      DXADislocationLines& result = *dxa_dislocation_lines;
      const double interval = *target_point_interval;
      const int smoothing = static_cast<int>( *target_smoothing_level );

      auto line_length = []( const std::vector<Vec3d>& pts ) -> double
      {
        double l = 0.0;
        for(size_t i=1;i<pts.size();i++) { l += norm( pts[i] - pts[i-1] ); }
        return l;
      };

      long n_lines_processed = 0;
      long points_before = 0, points_after = 0;
      double length_before = 0.0, length_after = 0.0;

      for(size_t i=0; i<result.line_positions.size(); i++)
      {
        if( i >= result.core_size.size() || result.core_size[i].size() != result.line_positions[i].size() ) { continue; }
        if( result.line_positions[i].size() < 2 ) { continue; }

        points_before += static_cast<long>( result.line_positions[i].size() );
        length_before += line_length( result.line_positions[i] );

        std::vector<Vec3d> coarsened;
        std::vector<int32_t> coarsened_core;
        coarsen_dislocation_line( interval, result.line_positions[i], result.core_size[i], coarsened, coarsened_core, result.is_loop[i] != 0 );
        smooth_dislocation_line( smoothing, coarsened, result.is_loop[i] != 0 );

        result.line_positions[i] = std::move(coarsened);
        result.core_size[i].clear(); // no longer meaningful once smoothed, matches OVITO's own segment->coreSize.clear()
        points_after += static_cast<long>( result.line_positions[i].size() );
        length_after += line_length( result.line_positions[i] );
        ++n_lines_processed;
      }

      lout << "smooth_dxa_dislocation_lines: " << n_lines_processed << " lines coarsened+smoothed, "
           << points_before << " -> " << points_after << " total points, total length "
           << length_before << " -> " << length_after << " Ang" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Reduces the tortuosity of raw swept dislocation lines (compute_dxa_circuit_sweep) via OVITO's real
two-stage line post-processing: adaptive coarsening (merge raw sweep points, weighted by local
circuit width) followed by Taubin smoothing. See this file's own header comment for the full
mechanism -- port of DislocationNetwork::coarsenDislocationLine()/smoothDislocationLine().

Usage example:

compute_dxa_circuit_sweep: {}
smooth_dxa_dislocation_lines: { target_point_interval: 2.5, target_smoothing_level: 1 }
write_dxa_dislocation_lines: { filename: "paraview/dxa_lines" }

)EOF";
    }
  };

  // === register factory ===
  ONIKA_AUTORUN_INIT(smooth_dxa_dislocation_lines)
  {
    OperatorNodeFactory::instance()->register_factory( "smooth_dxa_dislocation_lines", make_simple_operator< SmoothDXADislocationLines > );
  }

}
