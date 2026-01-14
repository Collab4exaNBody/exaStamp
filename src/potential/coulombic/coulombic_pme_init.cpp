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

/*
================================================================================
PME Initialization Operator
================================================================================

This operator initializes the Particle Mesh Ewald (PME) parameters including:
- Grid dimensions for the charge mesh
- B-spline interpolation order
- Ewald splitting parameter (alpha)
- Precomputed influence function (Green's function with B-spline correction)
- B-spline moduli for each dimension

The operator is called once at simulation startup and whenever domain
parameters change.

================================================================================
*/

#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_stream.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/domain.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <exaStamp/potential/coulombic/pme.h>
#include <exaStamp/potential/coulombic/pme_bspline.h>
#include <exaStamp/compute/thermodynamic_state.h>

namespace exaStamp
{
  using namespace exanb;

  class PMEInitOperator : public OperatorNode
  {
    // ========= I/O slots =======================
    
    // Input parameters
    ADD_SLOT( double     , accuracy_relative , INPUT , REQUIRED );
    ADD_SLOT( double     , alpha             , INPUT , REQUIRED );  // Ewald splitting parameter (0 = auto)
    ADD_SLOT( double     , rcut              , INPUT , REQUIRED );  // Real-space cutoff
    ADD_SLOT( long       , spline_order      , INPUT , 4 );         // B-spline order (default: 4 = cubic)
    ADD_SLOT( long       , mesh_nx           , INPUT , 0 );         // Mesh size X (0 = auto)
    ADD_SLOT( long       , mesh_ny           , INPUT , 0 );         // Mesh size Y (0 = auto)
    ADD_SLOT( long       , mesh_nz           , INPUT , 0 );         // Mesh size Z (0 = auto)
    ADD_SLOT( long       , mesh_size         , INPUT , 0 );         // Uniform mesh size (overrides nx/ny/nz if set)
    
    // Domain and system information
    ADD_SLOT( Domain     , domain            , INPUT , OPTIONAL );
    ADD_SLOT( double     , sum_square_charge , INPUT );
    ADD_SLOT( double     , sum_charge        , INPUT );
    ADD_SLOT( uint64_t   , natoms            , INPUT , OPTIONAL );
    
    // Output configuration
    ADD_SLOT( PMEParameters , pme_config     , INPUT_OUTPUT );
    ADD_SLOT( double     , rcut_out          , OUTPUT );
    ADD_SLOT( double     , rcut_max          , INPUT_OUTPUT , 0.0 );

  public:

    // -----------------------------------------------
    // Execute: Initialize PME parameters
    // -----------------------------------------------
    inline void execute() override final
    {
      using pme_constants::fpe0;
      using pme_constants::epsilonZero;
      using pme_constants::ONE_4PI_EPS0;

      if (!domain.has_value())
      {
        ldbg << "PME Init: domain not available, skipping initialization" << std::endl;
        *rcut_out = *rcut;
        *rcut_max = std::max(*rcut_max, *rcut);
        return;
      }

      auto& p = *pme_config;
      
      // Check if reinitialization is needed
      bool needs_init = false;
      
      // Check if alpha changed
      if (*alpha > 0.0 && *alpha != p.alpha) needs_init = true;
      
      // Check if volume is uninitialized
      if (p.volume == 0.0) needs_init = true;
      
      // Check if cutoff changed
      if (*rcut != p.rcut) needs_init = true;
      
      // Check if accuracy changed
      if (*accuracy_relative != p.accuracy_relative) needs_init = true;
      
      // Check if spline order changed
      if (*spline_order != p.spline_order) needs_init = true;
      
      // Check if mesh size changed
      int target_nx = (*mesh_size > 0) ? *mesh_size : ((*mesh_nx > 0) ? *mesh_nx : 0);
      int target_ny = (*mesh_size > 0) ? *mesh_size : ((*mesh_ny > 0) ? *mesh_ny : 0);
      int target_nz = (*mesh_size > 0) ? *mesh_size : ((*mesh_nz > 0) ? *mesh_nz : 0);
      
      if (target_nx > 0 && target_nx != p.grid_size.nx) needs_init = true;
      if (target_ny > 0 && target_ny != p.grid_size.ny) needs_init = true;
      if (target_nz > 0 && target_nz != p.grid_size.nz) needs_init = true;

      if (!needs_init)
      {
        ldbg << "PME Init: parameters unchanged, skipping reinitialization" << std::endl;
        *rcut_out = *rcut;
        *rcut_max = std::max(*rcut_max, *rcut);
        return;
      }

      // Get domain size
      auto domainSize = domain->bounds_size();
      auto xform = domain->xform();
      
      if (!is_diagonal(xform))
      {
        lerr << "PME Error: Domain XForm is not diagonal, cannot compute domain box size" << std::endl;
        std::abort();
      }
      
      domainSize = xform * domainSize;
      
      // Check periodic boundaries
      if (!(domain->periodic_boundary_x() && domain->periodic_boundary_y() && domain->periodic_boundary_z()))
      {
        lerr << "PME Error: Domain must be entirely periodic for PME" << std::endl;
        std::abort();
      }

      // Determine grid size
      PMEGridSize grid_size;
      
      if (*mesh_size > 0)
      {
        // Uniform mesh size specified
        grid_size.nx = grid_size.ny = grid_size.nz = *mesh_size;
      }
      else if (*mesh_nx > 0 && *mesh_ny > 0 && *mesh_nz > 0)
      {
        // Individual mesh sizes specified
        grid_size.nx = *mesh_nx;
        grid_size.ny = *mesh_ny;
        grid_size.nz = *mesh_nz;
      }
      else
      {
        // Auto-compute mesh size based on domain size and accuracy
        // Rule of thumb: ~1 grid point per Angstrom, rounded to FFT-friendly size
        auto compute_mesh_dim = [](double length, double accuracy) -> int
        {
          // Base estimate: 1 grid point per unit length
          int n = static_cast<int>(std::ceil(length));
          
          // Adjust for accuracy (finer mesh for higher accuracy)
          if (accuracy < 1e-6) n = static_cast<int>(n * 1.5);
          else if (accuracy < 1e-5) n = static_cast<int>(n * 1.2);
          
          // Round to FFT-friendly size (multiple of 2, 3, or 5 is good for cuFFT)
          // Find nearest power of 2 that's >= n
          int power2 = 1;
          while (power2 < n) power2 *= 2;
          
          // Also consider powers that include factors of 3 and 5
          int alt1 = (n / 6 + 1) * 6;   // Multiple of 6
          int alt2 = (n / 10 + 1) * 10; // Multiple of 10
          
          // Choose the smallest that's >= n
          int result = power2;
          if (alt1 >= n && alt1 < result) result = alt1;
          if (alt2 >= n && alt2 < result) result = alt2;
          
          // Ensure minimum size
          return std::max(result, 16);
        };
        
        grid_size.nx = compute_mesh_dim(domainSize.x, *accuracy_relative);
        grid_size.ny = compute_mesh_dim(domainSize.y, *accuracy_relative);
        grid_size.nz = compute_mesh_dim(domainSize.z, *accuracy_relative);
      }

      // Validate spline order
      int order = *spline_order;
      if (order < 3 || order > 8)
      {
        lerr << "PME Error: spline_order must be between 3 and 8, got " << order << std::endl;
        std::abort();
      }

      // Initialize PME parameters
      pme_init_parameters(
        domainSize,
        grid_size,
        order,
        *alpha,
        *rcut,
        *accuracy_relative,
        natoms.has_value() ? *natoms : 1,
        *sum_square_charge,
        p,
        ldbg << ""
      );

      // Log configuration
      if (p.volume > 0.0)
      {
        lout << "======== PME Configuration ========" << std::endl;
        lout << "Domain size     = " << domainSize << std::endl;
        lout << "Periodic        = " << std::boolalpha << domain->periodic_boundary_x()
             << " , " << std::boolalpha << domain->periodic_boundary_y()
             << " , " << std::boolalpha << domain->periodic_boundary_z() << std::endl;
        lout << "Grid size       = " << p.grid_size.nx << " x " << p.grid_size.ny << " x " << p.grid_size.nz << std::endl;
        lout << "Grid spacing    = " << p.grid_spacing.x << " x " << p.grid_spacing.y << " x " << p.grid_spacing.z << std::endl;
        lout << "Spline order    = " << p.spline_order << std::endl;
        lout << "Alpha (Ewald)   = " << p.alpha << std::endl;
        lout << "Real-space rcut = " << p.rcut << std::endl;
        lout << "Accuracy target = " << p.accuracy_relative << std::endl;
        lout << "Volume          = " << p.volume << std::endl;
        lout << "Influence func  = " << p.influence_function.size() << " k-points" << std::endl;
        lout << "===================================" << std::endl;
        
        // Performance estimate
        size_t mesh_total = grid_size.nx * grid_size.ny * grid_size.nz;
        double fft_ops = mesh_total * std::log2(mesh_total) * 5; // Approximate FFT flops
        lout << "PME mesh points = " << mesh_total << std::endl;
        lout << "Est. FFT ops    = " << fft_ops / 1e6 << " Mflops" << std::endl;
        lout << "===================================" << std::endl;
      }

      // Set output cutoff
      *rcut_out = *rcut;
      *rcut_max = std::max(*rcut_max, *rcut);
    }
  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(coulombic_pme_init)
  {
    OperatorNodeFactory::instance()->register_factory("coulombic_pme_init", make_simple_operator<PMEInitOperator>);
  }

}
