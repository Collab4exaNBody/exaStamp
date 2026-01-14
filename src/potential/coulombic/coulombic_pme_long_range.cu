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
Particle Mesh Ewald (PME) Long-Range Electrostatics Operator
================================================================================

This operator computes long-range electrostatic interactions using the PME
algorithm, which is O(N log N) compared to O(N*K) for direct Ewald summation.

Algorithm steps:
1. Compute B-spline coefficients for each particle
2. Spread particle charges onto 3D mesh using B-spline interpolation
3. FFT the charge mesh to reciprocal space
4. Multiply by influence function (Green's function with B-spline correction)
5. Compute gradient in k-space (for forces)
6. Inverse FFT to get real-space forces
7. Interpolate forces back to particle positions

Note: Short-range interactions (erfc term) are computed in a separate operator.

================================================================================
*/

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>

#include <exaStamp/potential/coulombic/pme.h>
#include <exaStamp/potential/coulombic/pme_bspline.h>
#include <exaStamp/potential/coulombic/pme_kernels.h>
#include <exaStamp/unit_system.h>

#include <exanb/core/config.h>

#include <exanb/core/parallel_grid_algorithm.h>
#include <onika/cuda/cuda.h>
#include <exanb/core/xform.h>
#include <mpi.h>
#include <cufft.h>

namespace exaStamp
{
  using namespace exanb;
  using pme_constants::fpe0;
  using pme_constants::epsilonZero;
  using pme_constants::ONE_4PI_EPS0;

  //============================================================================
  // PME Long-Range Operator
  //============================================================================
  
  template<
    class GridT,
    class = AssertGridHasFields<GridT, field::_ep, field::_fx, field::_fy, field::_fz>
  >
  class PMELongRangePC : public OperatorNode
  {
    // I/O Slots
    ADD_SLOT(PMEParameters, pme_config, INPUT, OPTIONAL);
    ADD_SLOT(GridT, grid, INPUT_OUTPUT);
    ADD_SLOT(Domain, domain, INPUT, REQUIRED);
    ADD_SLOT(double, rcut_max, INPUT_OUTPUT, 0.0);
    ADD_SLOT(ParticleSpecies, species, INPUT, REQUIRED);
    ADD_SLOT(MPI_Comm, mpi, INPUT);
    ADD_SLOT(bool, trigger_thermo_state, INPUT, OPTIONAL);
    ADD_SLOT(double, potential_energy_shift, OUTPUT);
    
    // PME working storage (persistent across timesteps)
    //    ADD_SLOT(PMEMesh, pme_mesh, INPUT_OUTPUT, PMEMesh{});
    ADD_SLOT(ParticleBSplineData, pme_bspline, INPUT_OUTPUT, ParticleBSplineData{});

  public:
    inline void execute() override final
    {
      bool log_energy = false;
      if (trigger_thermo_state.has_value())
      {
        log_energy = *trigger_thermo_state;
      }
      
      ldbg << "------------------------------" << std::endl;
      ldbg << "PME Long-Range Electrostatics" << std::endl;
      ldbg << "------------------------------" << std::endl;

      if (!pme_config.has_value())
      {
        ldbg << "pme_config not set, skip PME computation" << std::endl;
        return;
      }

      *rcut_max = std::max(*rcut_max, pme_config->rcut);
      
      if (grid->number_of_cells() == 0) return;
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      ssize_t gl = grid->ghost_layers();

      const auto& params = *pme_config;
      const int nx = params.grid_size.nx;
      const int ny = params.grid_size.ny;
      const int nz = params.grid_size.nz;

      // // Initialize mesh if needed
      // if (pme_mesh->nx != nx || pme_mesh->ny != ny || pme_mesh->nz != nz)
      // {
      //   pme_mesh->resize(nx, ny, nz);
      //   pme_mesh->init_fft_plans();
      //   ldbg << "PME mesh initialized: " << nx << "x" << ny << "x" << nz << std::endl;
      // }

      cudaStream_t stream = global_cuda_ctx()->getThreadStream(0);

      // Count particles and gather positions/charges
      size_t n_particles = 0;
      const size_t n_cells = grid->number_of_cells();
      
      for (size_t c = 0; c < n_cells; ++c)
      {
        n_particles += cells[c].size();
      }
      
      ldbg << "Processing " << n_particles << " particles" << std::endl;

      if (n_particles == 0) return;

      // Allocate temporary arrays for gathered particle data
      onika::memory::CudaMMVector<double> pos_x(n_particles);
      onika::memory::CudaMMVector<double> pos_y(n_particles);
      onika::memory::CudaMMVector<double> pos_z(n_particles);
      onika::memory::CudaMMVector<double> charges(n_particles);
      onika::memory::CudaMMVector<double> forces_x(n_particles, 0.0);
      onika::memory::CudaMMVector<double> forces_y(n_particles, 0.0);
      onika::memory::CudaMMVector<double> forces_z(n_particles, 0.0);

      // Gather particle data
      // Note: In a real implementation, this would be done more efficiently
      // using the compute_cell_particles framework
      size_t pidx = 0;
      Mat3d xform = domain->xform();
      GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) nowait)
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          size_t n = cells[i].size();
          const uint8_t* __restrict__ types = cells[i].field_pointer_or_null(field::type); ONIKA_ASSUME_ALIGNED(types);
          const double* __restrict__ rx = cells[i][field::rx]; ONIKA_ASSUME_ALIGNED(rx);
          const double* __restrict__ ry = cells[i][field::ry]; ONIKA_ASSUME_ALIGNED(ry);
          const double* __restrict__ rz = cells[i][field::rz]; ONIKA_ASSUME_ALIGNED(rz);
          
          for(size_t j=0;j<n;j++)
            {
              const double q = 1.0;
              Vec3d r = { rx[j] , ry[j] , rz[j] };
              r = xform * r;
              pos_x[pidx] = r.x;
              pos_y[pidx] = r.y;
              pos_z[pidx] = r.z;
              charges[pidx] = (*species)[ types[j]].m_charge;
              ++pidx;
            }
        }
      GRID_OMP_FOR_END

      //========================================================================
      // PME Algorithm Steps
      //========================================================================
      
      // Step 1: Compute B-spline coefficients
      pme_compute_bsplines(
        pos_x.data(), pos_y.data(), pos_z.data(),
        n_particles, params, *pme_bspline, stream);
      
      // // Step 2: Spread charges onto mesh
      // pme_spread_charges(
      //   charges.data(), n_particles, params, *pme_bspline, *pme_mesh, stream);
      
      // // Step 3: MPI reduction of charge mesh (for parallel runs)
      // cudaStreamSynchronize(stream);
      // MPI_Allreduce(MPI_IN_PLACE, pme_mesh->charge_mesh.data(),
      //               pme_mesh->charge_mesh.size(), MPI_DOUBLE, MPI_SUM, *mpi);
      
      // // Step 4: Forward FFT and k-space solve
      // pme_solve_kspace(params, *pme_mesh, stream);
      
      // // Step 5: Compute k-space forces and inverse FFT
      // pme_compute_forces_kspace(params, *pme_mesh, stream);
      
      // // Step 6: Interpolate forces back to particles
      // pme_interpolate_forces(
      //   charges.data(), n_particles, params, *pme_bspline, *pme_mesh,
      //   forces_x.data(), forces_y.data(), forces_z.data(), stream);
      
      // cudaStreamSynchronize(stream);

      // // Scatter forces back to grid
      // pidx = 0;
      // for (size_t c = 0; c < n_cells; ++c)
      // {
      //   auto cell = cells[c];
      //   const size_t n_in_cell = cell.size();
        
      //   auto fx_field = cell[field::fx];
      //   auto fy_field = cell[field::fy];
      //   auto fz_field = cell[field::fz];
        
      //   for (size_t i = 0; i < n_in_cell; ++i)
      //   {
      //     fx_field[i] += forces_x[pidx];
      //     fy_field[i] += forces_y[pidx];
      //     fz_field[i] += forces_z[pidx];
      //     ++pidx;
      //   }
      // }

      // //========================================================================
      // // Energy computation (only when logging)
      // //========================================================================
      
      // if (log_energy)
      // {
      //   *potential_energy_shift = 0.0;
        
      //   // Reciprocal space energy is computed from k-space mesh
      //   // E_recip = 0.5 * Σ_k |ρ(k)|² * G(k)
      //   // This is done on CPU for simplicity; could be GPU kernel
        
      //   const int nz_complex = nz / 2 + 1;
      //   double e_recip = 0.0;
        
      //   // Note: This should be computed before kspace_solve modifies the mesh
      //   // In production, use a separate kernel or save |ρ(k)|²
        
      //   // Self-energy correction
      //   double e_self = 0.0;
      //   double sum_q2 = 0.0;
      //   for (size_t i = 0; i < n_particles; ++i)
      //   {
      //     double q = charges[i];
      //     sum_q2 += q * q;
      //   }
        
      //   // E_self = -α/√π * Σᵢ qᵢ²
      //   e_self = -params.alpha / std::sqrt(M_PI) * sum_q2 * ONE_4PI_EPS0;
        
      //   *potential_energy_shift = e_recip + e_self;
        
      //   ldbg << "PME Energy: reciprocal = " << e_recip 
      //        << ", self = " << e_self << std::endl;
      // }

      ldbg << "------------------------" << std::endl;
      ldbg << "End PME Long-Range" << std::endl;
      ldbg << "------------------------" << std::endl;
    }
  };

  template<class GridT>
  using PMELongRangePCTmpl = PMELongRangePC<GridT>;

  // Register factory
  ONIKA_AUTORUN_INIT(coulombic_pme_long_range)
  {
    OperatorNodeFactory::instance()->register_factory(
      "coulombic_pme_long_range",
      make_grid_variant_operator<PMELongRangePCTmpl>
    );
  }

}
