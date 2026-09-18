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

#include <memory>

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <onika/physics/units.h>
#include <onika/memory/allocator.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/source_term.h>

#include <mpi.h>
#include <iomanip>

#include "ttm_ti_deposit.h"
#include "ttm_langevin_rng.h"
#include "ttm_langevin_coupling.h"
#include "ttm_laplacian.h"

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_vx, field::_vy, field::_vz , field::_fx, field::_fy, field::_fz >
    >
  class IonicElectronicHeatTransfer : public OperatorNode
  {
    ADD_SLOT( MPI_Comm       , mpi          , INPUT , MPI_COMM_WORLD );
    
    ADD_SLOT( GridT          , grid         , INPUT , REQUIRED );
    ADD_SLOT( Domain         , domain       , INPUT , REQUIRED );
    ADD_SLOT( ParticleSpecies, species      , INPUT , REQUIRED );
    ADD_SLOT( double         , rcut_max     , INPUT_OUTPUT , 0.0 );  // neighborhood distance, in grid space
    
    ADD_SLOT( double         , dt           , INPUT , REQUIRED );
    ADD_SLOT( double         , physical_time, INPUT , REQUIRED );
    ADD_SLOT( long           , timestep     , INPUT , REQUIRED );

    ADD_SLOT( ScalarSourceTermInstance , te_source  , INPUT_OUTPUT , std::make_shared<ScalarSourceTerm>() );
    ADD_SLOT( ScalarSourceTermInstance , ti_source  , INPUT_OUTPUT , std::make_shared<ScalarSourceTerm>() );
    // electron-phonon (Langevin) coupling friction, LAMMPS fix-ttm's gamma_p convention:
    // a force/velocity friction coefficient, not mass-scaled like langevin_thermostat's gamma.
    ADD_SLOT( double         , gamma_p      , INPUT , 0.0 );
    // electronic stopping: friction boosted to (gamma_p+gamma_s) above the v_0 velocity threshold.
    ADD_SLOT( double         , gamma_s      , INPUT , 0.0 );
    ADD_SLOT( double         , v_0          , INPUT , 0.0 );
    // false (default): Gaussian noise, sqrt(2*kB*gamma_p*Te/dt) prefactor (this codebase's own
    // langevin_thermostat convention). true: LAMMPS fix-ttm's own uniform[-0.5,0.5) noise,
    // sqrt(24*kB*gamma_p*Te/dt) prefactor (same target variance, different RNG family).
    ADD_SLOT( bool           , lammps_noise , INPUT , false );
    // false (default): one explicit-Euler diffusion step per MD step, as before. true: LAMMPS
    // fix-ttm's own stability check (fix_ttm.cpp end_of_step) -- if a single step would exceed
    // the explicit-diffusion stability limit, subdivide it into several smaller inner steps.
    ADD_SLOT( bool           , substep_diffusion , INPUT , false );
    ADD_SLOT( double         , Ke           , INPUT , 1.0 );
    ADD_SLOT( double         , Ce           , INPUT , 1.0 );
    ADD_SLOT( double         , rho_e        , INPUT , 1.0 );
    ADD_SLOT( double         , splat_size   , INPUT , 1.0 );

    ADD_SLOT( bool           , copy_ti_te   , INPUT, false );
    ADD_SLOT( long           , grid_subdiv  , INPUT , 3 );
    ADD_SLOT( GridCellValues , grid_cell_values      , INPUT_OUTPUT );
    ADD_SLOT( double         , total_electronic_energy , INPUT_OUTPUT , 0.0 );
    // energy transferred from electrons to ions this MD step (LAMMPS fix-ttm's transfer_energy,
    // f_twotemp[2]): sum over the whole grid of the Langevin coupling work, integrated over dt.
    ADD_SLOT( double         , ion_transfer_energy     , INPUT_OUTPUT , 0.0 );

    // Persistent, GPU-visible (onika::memory::CudaMMVector) scratch storage for per-cell Ti,
    // Laplacian(Te), and the Langevin coupling energy sink -- reused across MD steps instead of
    // being freshly heap-allocated and zeroed on every execute() call.
    ADD_SLOT( onika::memory::CudaMMVector<double> , ttm_scratch_ti               , PRIVATE );
    ADD_SLOT( onika::memory::CudaMMVector<double> , ttm_scratch_lap_te           , PRIVATE );
    ADD_SLOT( onika::memory::CudaMMVector<double> , ttm_scratch_energy_transfer  , PRIVATE );
    // Se/Si source-term values, precomputed once per outer MD step (see execute(): they only
    // depend on cell center + physical_time, neither of which changes across inner substeps).
    ADD_SLOT( onika::memory::CudaMMVector<double> , ttm_scratch_se               , PRIVATE );
    ADD_SLOT( onika::memory::CudaMMVector<double> , ttm_scratch_si               , PRIVATE );
    // Per-species mass table, GPU-visible: sized to nSpecies (not a fixed MAX_PARTICLE_TYPES array
    // captured by value -- that blows past onika's GPU kernel-functor size cap, see ttm_ti_deposit.h).
    ADD_SLOT( onika::memory::CudaMMVector<double> , ttm_scratch_species_mass     , PRIVATE );

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      // compile time constant indicating if grid has type field
      using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
      static constexpr has_type_field_t has_type_field{};

      // compile time constant indicating if grid has an id field (Langevin coupling RNG keying, see
      // ttm_langevin_rng.h) -- same has_id_field_t fallback convention as langevin_thermostat.cpp.
      using has_id_field_t = typename GridT::CellParticles::template HasField < field::_id > ;
      static constexpr bool has_id_field = has_id_field_t::value;

      static constexpr double weight_sum_epsilon = 1.e-13; if constexpr (weight_sum_epsilon==0.0){}
      static constexpr double te_deviation_epsilon = 1.e-13;

      const double cell_size = domain->cell_size();
      const ssize_t subdiv = *grid_subdiv;
      const ssize_t n_subcells = subdiv * subdiv * subdiv;
      const double subcell_size = cell_size / subdiv;
      const double subcell_volume = subcell_size * subcell_size * subcell_size;

      // retreive field data accessor. create data field if needed
      if( ! grid_cell_values->has_field("te") )
      {
        grid_cell_values->add_field("te",subdiv,1);
      }
      assert( size_t(subdiv) == grid_cell_values->field("te").m_subdiv );
      assert( size_t(subdiv * subdiv * subdiv) == grid_cell_values->field("te").m_components );
      auto cell_te_data = grid_cell_values->field_data("te");
      
      // particle splatting size
      const double sp_size = *splat_size;
      // const double sp_volume = sp_size * sp_size * sp_size;

      if( sp_size > subcell_size )
      {
        lerr << "in " << pathname() << std::endl
             << "splat_size = "<<sp_size << " is larger than subcell size = "<<subcell_size<< std::endl
             << "Choose smaller splat_size or coarser grid_subdiv value" << std::endl;
        std::abort();
      }

      // Ghost region only has to cover the splat radius (splat_size/2). It must NOT depend on
      // subcell_size/cell_size: this runs during preinit_rcut_max, when domain->cell_size() is still
      // a placeholder (e.g. 17.6 ang) and rcut_max is never lowered afterwards, so a cell-size-based
      // value inflates the neighbor/ghost distance for the whole run (~15x slower eam_alloy_force).
      const double splat_reach = sp_size / 2.0;
      const double rcut_max_grid = (*rcut_max) / domain->xform_min_scale();
      if( splat_reach > rcut_max_grid )
      {
        ldbg << "in " << pathname() << std::endl
             << "splat_size/2.0 = "<<splat_reach<<" is larger than rcut_max/min_scale = "
             << *(rcut_max) <<'/'<<domain->xform_min_scale()<<" = "<<rcut_max_grid<<std::endl
             << "adjust rcut_max "<<rcut_max_grid<<" -> "<< splat_reach * domain->xform_min_scale() <<std::endl ;
        *rcut_max = std::max( *rcut_max , splat_reach * domain->xform_min_scale() );
      }

      if( grid->number_of_cells() == 0 )
      {
        return;
      }
  
//      static const double k = UnityConverterHelper::convert(onika::physics::boltzmann, "J/K");
      //ldbg << "cell_heat: dt="<<(*dt)<<std::endl;

      size_t nSpecies = species->size();
      if( nSpecies != 1 && !has_type_field )
      {
        lerr << pathname() << std::endl;
        lerr << "no type information, can't retreive masses" << std::endl;
        std::abort();
      }
      // Persistent, GPU-visible species-mass table (CudaMMVector, like the other ttm_scratch_*
      // buffers) instead of a std::vector, whose data() pointer isn't valid on-device.
      ttm_scratch_species_mass->assign( nSpecies , 0.0 );
      double* masses = ttm_scratch_species_mass->data();
      for(size_t i=0;i<nSpecies;i++) { masses[i] = species->at(i).m_mass; }

      auto cells = grid->cells();
      const ssize_t n_cells = grid->number_of_cells();
      const IJK dims = grid->dimension();
      const ssize_t gl = grid->ghost_layers();      

      // Persistent scratch storage for per cell Ti and per cell Laplacian(Te) -- .assign() reuses
      // existing capacity across calls (no realloc) instead of a fresh std::vector every step.
      ttm_scratch_ti->assign( n_cells * n_subcells , 0.0 );
      ttm_scratch_lap_te->assign( n_cells * n_subcells , 0.0 );
      double* Ti = ttm_scratch_ti->data();
      double* cell_L_Te = ttm_scratch_lap_te->data();

      const auto& te_source_func = * (*te_source);
      const auto& ti_source_func = * (*ti_source);

      // 1. computes per cell Ti (GPU-portable: exanb::compute_cell_particles + TtmTiDepositFunctor,
      // see ttm_ti_deposit.h -- replaces the old #pragma omp parallel / GRID_OMP_FOR_BEGIN block).
      {
        TtmTiDepositFunctor ti_deposit_func = {
          grid->origin(), grid->offset(), dims, subdiv,
          cell_size, subcell_size, subcell_volume, sp_size,
          masses,
          Ti
        };
        if constexpr ( has_type_field_t::value )
        {
          compute_cell_particles( *grid, true, ti_deposit_func,
            FieldSet<field::_rx,field::_ry,field::_rz,field::_vx,field::_vy,field::_vz,field::_type>{},
            parallel_execution_context() );
        }
        else
        {
          compute_cell_particles( *grid, true, ti_deposit_func,
            FieldSet<field::_rx,field::_ry,field::_rz,field::_vx,field::_vy,field::_vz>{},
            parallel_execution_context() );
        }
      }

      // if a simple copy Te <- Ti is requested, stop here
      if( *copy_ti_te )
      {
        ldbg << "copy_ti_te : dims=" <<dims<< std::endl;
        for(ssize_t i=0;i<n_cells;i++)
        {
          for(ssize_t j=0;j<n_subcells; j++)
          {
            cell_te_data.m_data_ptr[ i*cell_te_data.m_stride + j ] = Ti[i*n_subcells+j];
          }
        }
        return;
      }

      // Laplacian computation of Te will need some ghost
      if( gl <= 0 )
      {
        lerr << pathname() << std::endl;
        lerr << "No ghost layers, can't continue" << std::endl;
        std::abort();
      }
            
      // 2) Te<->Te and Te<->Ti heat transfer
      //const Vec3d grid_origin = grid->grid_bounds().bmin;
      const Mat3d xform = domain->xform();
      const double gamma_p_coupling = *gamma_p;
      const double Te_cond = *Ke;
      const double Ce_rho_e = (*Ce) * (*rho_e);
      const double delta_t = *dt;
      auto* Te = cell_te_data.m_data_ptr;

      ldbg << "Ke="<<Te_cond<<", Ce="<<(*Ce)<<", rho_e="<<(*rho_e)<<", gamma_p="<<gamma_p_coupling<< std::endl;

      // LAMMPS fix-ttm never touches Te (nor draws any coupling noise) until end_of_step() first
      // runs, which only happens after step 1 completes -- its own step-0 diagnostics reflect the
      // untouched initial condition (flangevin/net_energy_transfer start zero-initialized, and
      // setup()'s post_force_setup() only re-applies that still-zero force, drawing no noise).
      // Match that exactly here instead of running a full coupling+diffusion pass at timestep 0.
      if( *timestep <= 0 )
      {
        double sum_Te_initial = 0.0;
#       pragma omp parallel
        {
          GRID_OMP_FOR_BEGIN(dims,cell_i,cell_loc, schedule(static) reduction(+:sum_Te_initial) )
          {
            if( ! grid->is_ghost_cell(cell_loc) )
            {
              for(ssize_t sc=0;sc<n_subcells;sc++) { sum_Te_initial += Te[ cell_i*cell_te_data.m_stride + sc ]; }
            }
          }
          GRID_OMP_FOR_END
        }
        MPI_Allreduce(MPI_IN_PLACE,&sum_Te_initial,1,MPI_DOUBLE,MPI_SUM,*mpi);
        *total_electronic_energy = sum_Te_initial * subcell_volume * Ce_rho_e;
        *ion_transfer_energy = 0.0;
        return;
      }

      // Electron-phonon coupling: additive Langevin force on each particle, using the
      // (pre-diffusion-update) local Te, LAMMPS fix-ttm style: F = -gamma_p*v + noise*sqrt(2*kB*gamma_p*Te/dt).
      // Unlike the old deterministic g*(Te-Ti)-rescaling scheme, this can inject energy into
      // particles starting at rest and never divides by a near-zero ion temperature.
      // The work done on each particle (F.v) is deposited back as a Te sink, so the electron
      // bath cools down as it heats the lattice (energy-conserving, folded into dTe below).
      ttm_scratch_energy_transfer->assign( n_cells * n_subcells , 0.0 );
      double* cell_energy_transfer = ttm_scratch_energy_transfer->data();
      *ion_transfer_energy = 0.0;
      if( gamma_p_coupling != 0.0 )
      {
        const double kB = onika::physics::make_quantity( onika::physics::boltzmann, "J/K" ).convert();
        const bool use_lammps_noise = *lammps_noise;
        const double noise_variance_factor = use_lammps_noise ? 24.0 : 2.0; // uniform vs gaussian fluctuation-dissipation prefactor

        // Pass B, GPU-portable: exanb::compute_cell_particles + TtmLangevinCouplingFunctor, see
        // ttm_langevin_coupling.h -- replaces the old #pragma omp parallel / GRID_OMP_FOR_BEGIN block.
        TtmLangevinCouplingFunctor<decltype(cells),has_type_field_t::value,has_id_field> coupling_func = {
          cells,
          grid->origin(), grid->offset(), dims, gl, subdiv,
          cell_size, subcell_size, sp_size,
          gamma_p_coupling, *gamma_s, (*v_0)*(*v_0), kB, delta_t, noise_variance_factor, use_lammps_noise, uint64_t(*timestep),
          masses,
          Te, cell_te_data.m_stride,
          cell_energy_transfer
        };
        compute_cell_particles( *grid, true, coupling_func,
          FieldSet<field::_rx,field::_ry,field::_rz,field::_vx,field::_vy,field::_vz,field::_fx,field::_fy,field::_fz>{},
          parallel_execution_context() );

        // Pass C, confirmed already GPU-compatible as-is (staged plan's Stage 7): a plain
        // host-side reduction over cell_energy_transfer, which is already GPU-visible managed
        // memory (Stage 1) filled by Pass B's GPU kernel just above -- correct with an implicit
        // sync, nothing to port. LAMMPS fix-ttm's transfer_energy (f_twotemp[2]): total energy
        // transferred from electrons to ions over this whole MD step -- cell_energy_transfer
        // already holds it (fixed, computed once above), just sum it over this rank's own
        // (non-ghost) cells.
        double sum_ion_transfer = 0.0;
#       pragma omp parallel
        {
          GRID_OMP_FOR_BEGIN(dims,cell_i,cell_loc, schedule(static) reduction(+:sum_ion_transfer) )
          {
            if( ! grid->is_ghost_cell(cell_loc) )
            {
              for(ssize_t sc=0;sc<n_subcells;sc++) { sum_ion_transfer += cell_energy_transfer[ cell_i*n_subcells + sc ]; }
            }
          }
          GRID_OMP_FOR_END
        }
        MPI_Allreduce(MPI_IN_PLACE,&sum_ion_transfer,1,MPI_DOUBLE,MPI_SUM,*mpi);
        *ion_transfer_energy = sum_ion_transfer;
      }

      // Explicit-diffusion stability limit (LAMMPS fix_ttm.cpp end_of_step, isotropic grid so
      // dx=dy=dz=subcell_size): if a single step of size delta_t would violate it, subdivide
      // into several smaller inner steps of size inner_dt instead. The coupling sink
      // (cell_energy_transfer, computed once above from Te at the start of this MD step, exactly
      // like LAMMPS's post_force-once/end_of_step-substepped split) is reapplied unchanged at
      // every inner step -- its total contribution over the whole MD step is unaffected by how
      // many pieces the diffusion update is split into.
      long num_inner_timesteps = 1;
      double inner_dt = delta_t;
      if( *substep_diffusion )
      {
        const double diffusion_rate = Te_cond * ( 3.0 / (subcell_size*subcell_size) ); // sum_axes 1/dx^2, isotropic
        const double stability_criterion = 1.0 - 2.0*delta_t/Ce_rho_e*diffusion_rate;
        if( stability_criterion < 0.0 )
        {
          inner_dt = 0.5*Ce_rho_e / diffusion_rate;
          num_inner_timesteps = static_cast<long>( delta_t/inner_dt ) + 1;
          inner_dt = delta_t / double(num_inner_timesteps);
          ldbg << "substep_diffusion: stability_criterion="<<stability_criterion
               <<" -> num_inner_timesteps="<<num_inner_timesteps<<", inner_dt="<<inner_dt<< std::endl;
        }
      }

      // Precompute Se/Si once per outer MD step: (center, *physical_time) don't change across
      // inner substeps, so re-invoking the virtual ScalarSourceTerm calls every substep (as
      // before) was purely redundant work. This also removes the only per-subcell virtual
      // dispatch from the substep loop below.
      ttm_scratch_se->assign( n_cells * n_subcells , 0.0 );
      ttm_scratch_si->assign( n_cells * n_subcells , 0.0 );
      double* cell_Se = ttm_scratch_se->data();
      double* cell_Si = ttm_scratch_si->data();
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,cell_i,cell_loc, schedule(static) )
        {
          const Vec3d cell_origin = grid->cell_position( cell_loc );
          for(int ck=0;ck<subdiv;ck++)
          for(int cj=0;cj<subdiv;cj++)
          for(int ci=0;ci<subdiv;ci++)
          {
            IJK sc { ci, cj, ck };
            Vec3d scr { ci+0.5, cj+0.5, ck+0.5 };
            const size_t idx = cell_i*n_subcells + grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , sc );
            const Vec3d center = xform * ( cell_origin + scr*subcell_size );
            cell_Si[idx] = ti_source_func( center, *physical_time );
            cell_Se[idx] = te_source_func( center, *physical_time );
          }
        }
        GRID_OMP_FOR_END
      }

      for(long te_istep=0; te_istep<num_inner_timesteps; ++te_istep)
      {

      double sum_Te = 0.0;
      // norm_dTe/entropy_Te only feed the ldbg print below (sum_Te, unlike them, is also needed
      // unconditionally for the te_dev sanity check further down via old_sum_Te) -- skip their
      // per-cell computation and MPI payload in release builds.
#     ifndef NDEBUG
      double norm_dTe = 0.0;
      double entropy_Te = 0.0;
#     endif

      // Pass D (Laplacian), GPU-portable: onika::parallel::block_parallel_for + TtmLaplacianFunctor,
      // see ttm_laplacian.h -- replaces the old #pragma omp parallel / GRID_OMP_FOR_BEGIN stencil
      // loop. Also fixes a real bug the old loop had: it read Te via a bare n_subcells-based index
      // instead of cell_te_data.m_stride (the field's real per-cell stride, which only equals
      // n_subcells when "te" is the sole field on grid_cell_values -- every other Te access in this
      // file, e.g. Pass B/E and init_ttm.cpp, already used m_stride correctly). See ttm_laplacian.h.
      {
        TtmLaplacianFunctor laplacian_func = {
          dims, subdiv, subcell_size,
          Te, cell_te_data.m_stride,
          cell_L_Te
        };
        onika::parallel::block_parallel_for( n_cells * n_subcells, laplacian_func, parallel_execution_context() );
      }

      // Diagnostic reduction over the now GPU-filled cell_L_Te / Te, kept host-side post-kernel
      // (accepts an implicit sync -- same approach as Stage 3/7's diagnostics; managed memory makes
      // reading back from the host correct as-is once the kernel dispatch above has returned).
#     pragma omp parallel
      {
#       ifndef NDEBUG
        GRID_OMP_FOR_BEGIN(dims,cell_i,cell_loc, schedule(static) reduction(+:sum_Te,norm_dTe,entropy_Te) )
#       else
        GRID_OMP_FOR_BEGIN(dims,cell_i,cell_loc, schedule(static) reduction(+:sum_Te) )
#       endif
        {
          if( ! grid->is_ghost_cell(cell_loc) )
          {
            for(int ck=0;ck<subdiv;ck++)
            for(int cj=0;cj<subdiv;cj++)
            for(int ci=0;ci<subdiv;ci++)
            {
              IJK sc { ci, cj, ck };
              const size_t scindex = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , sc );
              const size_t idx_te = cell_i*cell_te_data.m_stride + scindex;
              sum_Te += Te[idx_te];
#             ifndef NDEBUG
              const size_t j = cell_i*n_subcells + scindex;
              const double dTe = ( (Te_cond*cell_L_Te[j]) / Ce_rho_e ) * inner_dt;
              norm_dTe += std::fabs( dTe );
              entropy_Te += Te[idx_te] * std::log(Te[idx_te]) * subcell_volume ;
#             endif
            }
          }
        }
        GRID_OMP_FOR_END
      }

      // check global absolute electronic energy variation (from dissipation)
#     ifndef NDEBUG
      {
        double tmp[3] = { sum_Te, norm_dTe, entropy_Te };
        MPI_Allreduce(MPI_IN_PLACE,tmp,3,MPI_DOUBLE,MPI_SUM,*mpi);
        sum_Te = tmp[0];
        norm_dTe = tmp[1];
        entropy_Te = tmp[2];
        ldbg <<"sum_Te="<<sum_Te <<" norm_dTe/sum_Te="<<norm_dTe/sum_Te<<" , entropy_Te="<<entropy_Te << std::endl;
      }
#     else
      {
        // sum_Te alone is still needed unconditionally (feeds old_sum_Te / te_dev below).
        MPI_Allreduce(MPI_IN_PLACE,&sum_Te,1,MPI_DOUBLE,MPI_SUM,*mpi);
      }
#     endif


      // 3. Compute Te dissipation, Te<->Ti transfer & source terms
      double old_sum_Te = sum_Te;
      double sum_dTe = 0.0;
      double sum_Se = 0.0;
      double sum_Si = 0.0;
      sum_Te = 0.0;
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,cell_i,cell_loc, schedule(static) reduction(+:sum_Te,sum_dTe,sum_Se,sum_Si) )
        {
          //const IJK cell_loc = loc + gl;
          //const size_t cell_i = grid_ijk_to_index( dims , cell_loc );

          for(int ck=0;ck<subdiv;ck++)
          for(int cj=0;cj<subdiv;cj++)
          for(int ci=0;ci<subdiv;ci++)
          {
            IJK sc { ci, cj, ck };
            const size_t scindex = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , sc );
            const size_t idx_ti = cell_i*n_subcells + scindex ;
            const size_t idx_te = cell_i*cell_te_data.m_stride + scindex;

            // source terms (precomputed once per outer MD step, above)
            const double Si = cell_Si[idx_ti];
            const double Se = cell_Se[idx_ti];

            // sum external contributions (source terms)
            sum_Se += Se;
            sum_Si += Si;

            // Te/Ti coupling: cell_energy_transfer holds a FIXED total energy for the whole
            // outer MD step (computed once above, before the diffusion update); divide by
            // delta_t (not inner_dt) to get the power density sink used at every inner step.
            const double coupling_sink = cell_energy_transfer[idx_ti] / subcell_volume / delta_t;

            // cell temperature increments
            const double dTe = ( Te_cond*cell_L_Te[idx_ti] - coupling_sink + Se ) / Ce_rho_e;
            Te[idx_te] += dTe * inner_dt;
            if( ! grid->is_ghost_cell(cell_loc) )
            {
              sum_Te += Te[idx_te];
              sum_dTe += dTe * inner_dt;
            }
          }
        }
        GRID_OMP_FOR_END
      }

      {
        double tmp[4] = { sum_Te, sum_dTe, sum_Se, sum_Si };
        MPI_Allreduce(MPI_IN_PLACE,tmp,4,MPI_DOUBLE,MPI_SUM,*mpi);
        sum_Te = tmp[0];
        sum_dTe = tmp[1];
        sum_Se = tmp[2];
        sum_Si = tmp[3];
      }

      if(sum_Te > 0.)
      {
        double te_dev = (sum_Te - old_sum_Te - sum_dTe) / sum_Te;
        ldbg << "Te dev="<< te_dev << " sum_dTe="<<sum_dTe << " sum_Se="<<sum_Se<<" sum_Si="<<sum_Si<< std::endl;
        if( te_dev > te_deviation_epsilon )
        {
          lerr << "Te deviation too big : "<<old_sum_Te<<" -> "<<sum_Te<<" , dev="<<te_dev<<std::endl;
        }
        // else { lout << "Te deviation Ok : "<<old_sum_Te<<" -> "<<sum_Te<<" , dev="<<te_dev<<std::endl; }
      }
      // total electron thermal energy, LAMMPS fix-ttm's e_energy = sum(Te*Ce*rho_e*cell_volume):
      // sum_Te above is a raw sum of temperatures, not energy -- must be weighted by the heat
      // capacity to become one. (Harmless to omit while Ce(Te) was hardcoded to 1.0, but wrong
      // now that Ce/rho_e are real physical inputs.)
      *total_electronic_energy = sum_Te * subcell_volume * Ce_rho_e;

      } // for( te_istep ... num_inner_timesteps )
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
Handles heat transfer between ionic and electronic temperatures. Electronic temperature (Te) is held by a rectilinear grid,
while ionic temperature (Ti) commes from particles kinetic energy.
1. compute per cell Ti (diagnostic / copy_ti_te only)
2. apply an additive Langevin force (friction gamma_p, boosted to gamma_p+gamma_s above velocity
   v_0, plus noise -- gaussian by default, or LAMMPS fix-ttm's own uniform noise if lammps_noise:
   true) to each particle using the local Te; the work done on particles is deposited back as a Te sink
3. solve the heat equation on the rectilinear grid for Te (conduction Ke/(Ce*rho_e), the coupling sink, and
   source terms); if substep_diffusion: true and a single MD step would exceed the explicit-diffusion
   stability limit (LAMMPS fix-ttm style), this step is subdivided into several smaller inner steps
)EOF";
    }

  };

  template<class GridT> using IonicElectronicHeatTransferTmpl = IonicElectronicHeatTransfer<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(ionic_electronic_heat_transfer)
  {
   OperatorNodeFactory::instance()->register_factory("ionic_eletronic_heat_transfer", make_grid_variant_operator< IonicElectronicHeatTransferTmpl > );
  }

}
