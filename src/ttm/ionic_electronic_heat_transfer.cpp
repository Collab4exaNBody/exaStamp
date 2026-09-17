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
#include <onika/parallel/random.h>

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

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      // compile time constant indicating if grid has type field
      using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
      static constexpr has_type_field_t has_type_field{};

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

      const double rcut_max_grid = (*rcut_max) / domain->xform_min_scale();
      if( subcell_size+sp_size/2.0 > rcut_max_grid )
      {
        ldbg << "in " << pathname() << std::endl
             << "subcell_size+sp_size/2.0 = "<<subcell_size<<'+'<<sp_size<<"/2.0 = "<<(subcell_size+sp_size/2.0) << std::endl
             << "is larger than" << std::endl
             << "rcut_max/min_scale = "<< *(rcut_max) <<'/'<<domain->xform_min_scale()<<" = "<<rcut_max_grid<<std::endl
             << "adjust rcut_max "<<rcut_max_grid<<" -> "<< (subcell_size+sp_size/2.0) * domain->xform_min_scale() <<std::endl ;
        *rcut_max = std::max( *rcut_max , (subcell_size+sp_size/2.0) * domain->xform_min_scale() );
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
      std::vector<double> masses(nSpecies);
      for(size_t i=0;i<nSpecies;i++) { masses[i] = species->at(i).m_mass; }

      auto cells = grid->cells();
      const ssize_t n_cells = grid->number_of_cells();
      const IJK dims = grid->dimension();
      const ssize_t gl = grid->ghost_layers();      

      // Temporary storage for per cell Ti and per cell Laplacian(Te)
      std::vector<double> tmp_storage_ti_LapTe( n_cells * n_subcells * 2 , 0.0 );
      double* Ti = tmp_storage_ti_LapTe.data();
      double* cell_L_Te = tmp_storage_ti_LapTe.data() + n_cells * n_subcells;

      const auto& te_source_func = * (*te_source);
      const auto& ti_source_func = * (*ti_source);

      // 1. computes per cell Ti
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,i,cell_loc, schedule(dynamic) )
        {
	      const Vec3d cell_origin = grid->cell_position( cell_loc );
	        
	      GridFieldSetPointerTuple< GridT, FieldSet<field::_rx,field::_ry,field::_rz,field::_vx,field::_vy,field::_vz> > ptrs;
	      cells[i].capture_pointers( ptrs );

          const auto* __restrict__ rx = ptrs[field::rx];
          const auto* __restrict__ ry = ptrs[field::ry];
          const auto* __restrict__ rz = ptrs[field::rz];

          const auto* __restrict__ vx = ptrs[field::vx];
          const auto* __restrict__ vy = ptrs[field::vy];
          const auto* __restrict__ vz = ptrs[field::vz];

          const auto* __restrict__ atom_type = cells[i].field_pointer_or_null(field::type);

          const unsigned int n = cells[i].size();
          for(unsigned int j=0;j<n;j++)
          {
            Vec3d r { rx[j] , ry[j] , rz[j] };
            Vec3d v { vx[j] , vy[j] , vz[j] };

            const double mass = get_mass( j, atom_type, masses.data(), has_type_field );
            const double v2 = norm2(v);
 
            IJK center_cell_loc;
            IJK center_subcell_loc;
            Vec3d rco = r - cell_origin;
            localize_subcell( rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc );
            center_cell_loc += cell_loc;

            [[maybe_unused]] double sum_w = 0.0;
            [[maybe_unused]] int nb_contribs = 0;
            [[maybe_unused]] int nb_neighbors = 0;

            for(int ck=-1;ck<=1;ck++)
            for(int cj=-1;cj<=1;cj++)
            for(int ci=-1;ci<=1;ci++)
            {
              ++ nb_neighbors;
              IJK nbh_cell_loc;
              IJK nbh_subcell_loc;
              gcv_subcell_neighbor( center_cell_loc, center_subcell_loc, subdiv, IJK{ci,cj,ck}, nbh_cell_loc, nbh_subcell_loc );
              if( grid->contains(nbh_cell_loc) )
              {
                ++ nb_contribs;
                ssize_t nbh_cell_i = grid_ijk_to_index( dims , nbh_cell_loc );
                ssize_t nbh_subcell_i = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , nbh_subcell_loc );
                assert( nbh_cell_i>=0 && nbh_cell_i<n_cells );
                assert( nbh_subcell_i>=0 && nbh_subcell_i<n_subcells );

                // compute weighted contribution of particle to sub cell
                Vec3d nbh_cell_origin = grid->cell_position(nbh_cell_loc);
                AABB subcell_box = { nbh_cell_origin + nbh_subcell_loc*subcell_size , nbh_cell_origin + (nbh_subcell_loc+1)*subcell_size };
                const double w = particle_smoothing(r, sp_size, subcell_box);
                
                sum_w += w;
                size_t scindex = nbh_cell_i * n_subcells + nbh_subcell_i;
                double Ti_contrib = ( mass * w * v2 ) / subcell_volume;
                
#               pragma omp atomic update
                Ti[ scindex ] += Ti_contrib;
              }

            }

            // particle's contribution accurately distributed
            assert( (nb_contribs<nb_neighbors) || std::fabs(sum_w-1.0) < weight_sum_epsilon );
          }

        }
        GRID_OMP_FOR_END
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
      std::vector<double> cell_energy_transfer( n_cells * n_subcells , 0.0 );
      *ion_transfer_energy = 0.0;
      if( gamma_p_coupling != 0.0 )
      {
        const double kB = onika::physics::make_quantity( onika::physics::boltzmann, "J/K" ).convert();
        const double gamma_s_coupling = *gamma_s;
        const double v_0_sq = (*v_0) * (*v_0);
        const bool use_lammps_noise = *lammps_noise;
        const double noise_variance_factor = use_lammps_noise ? 24.0 : 2.0; // uniform vs gaussian fluctuation-dissipation prefactor
#       pragma omp parallel
        {
          auto& re = onika::parallel::random_engine();
          std::normal_distribution<double> gauss_rand(0.,1.);
          std::uniform_real_distribution<double> uniform_rand(-0.5,0.5);
          auto noise = [&]() -> double { return use_lammps_noise ? uniform_rand(re) : gauss_rand(re); };
          GRID_OMP_FOR_BEGIN(dims,i,cell_loc, schedule(dynamic) )
          {
            const Vec3d cell_origin = grid->cell_position( cell_loc );
            const bool is_local_cell = ! grid->is_ghost_cell(cell_loc);

            GridFieldSetPointerTuple< GridT, FieldSet<field::_rx,field::_ry,field::_rz,field::_vx,field::_vy,field::_vz,field::_fx,field::_fy,field::_fz> > ptrs;
            cells[i].capture_pointers( ptrs );

            const auto* __restrict__ rx = ptrs[field::rx];
            const auto* __restrict__ ry = ptrs[field::ry];
            const auto* __restrict__ rz = ptrs[field::rz];
            const auto* __restrict__ vx = ptrs[field::vx];
            const auto* __restrict__ vy = ptrs[field::vy];
            const auto* __restrict__ vz = ptrs[field::vz];
            auto* __restrict__ fx = ptrs[field::fx];
            auto* __restrict__ fy = ptrs[field::fy];
            auto* __restrict__ fz = ptrs[field::fz];
            const auto* __restrict__ atom_type = cells[i].field_pointer_or_null(field::type);

            const unsigned int n = cells[i].size();
            for(unsigned int j=0;j<n;j++)
            {
              Vec3d r { rx[j] , ry[j] , rz[j] };
              Vec3d v { vx[j] , vy[j] , vz[j] };

              IJK center_cell_loc;
              IJK center_subcell_loc;
              Vec3d rco = r - cell_origin;
              localize_subcell( rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc );
              center_cell_loc += cell_loc;

              // gather local Te (weighted average, same splat kernel as the Ti deposit)
              // and remember per-neighbor weights to re-deposit the coupling work below.
              double nbh_w[27];
              size_t nbh_idx[27];
              int nbh_count = 0;
              double Te_local = 0.0;
              for(int ck=-1;ck<=1;ck++)
              for(int cj=-1;cj<=1;cj++)
              for(int ci=-1;ci<=1;ci++)
              {
                IJK nbh_cell_loc;
                IJK nbh_subcell_loc;
                gcv_subcell_neighbor( center_cell_loc, center_subcell_loc, subdiv, IJK{ci,cj,ck}, nbh_cell_loc, nbh_subcell_loc );
                if( grid->contains(nbh_cell_loc) )
                {
                  ssize_t nbh_cell_i = grid_ijk_to_index( dims , nbh_cell_loc );
                  ssize_t nbh_subcell_i = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , nbh_subcell_loc );
                  Vec3d nbh_cell_origin = grid->cell_position(nbh_cell_loc);
                  AABB subcell_box = { nbh_cell_origin + nbh_subcell_loc*subcell_size , nbh_cell_origin + (nbh_subcell_loc+1)*subcell_size };
                  const double w = particle_smoothing(r, sp_size, subcell_box);
                  const size_t scindex = nbh_cell_i * n_subcells + nbh_subcell_i;
                  Te_local += w * Te[ nbh_cell_i * cell_te_data.m_stride + nbh_subcell_i ];
                  nbh_w[nbh_count] = w;
                  nbh_idx[nbh_count] = scindex;
                  ++nbh_count;
                }
              }

              if( Te_local > 0.0 )
              {
                // electronic stopping: friction boosted above the v_0 velocity threshold (LAMMPS
                // fix-ttm convention) -- only the friction term is boosted, not the noise term.
                double friction = gamma_p_coupling;
                if( gamma_s_coupling != 0.0 && norm2(v) > v_0_sq ) { friction += gamma_s_coupling; }

                const double noise_amplitude = std::sqrt( noise_variance_factor * kB * gamma_p_coupling * Te_local / delta_t );
                const Vec3d f_langevin {
                  -friction * v.x + noise() * noise_amplitude ,
                  -friction * v.y + noise() * noise_amplitude ,
                  -friction * v.z + noise() * noise_amplitude
                };

                // energy actually injected into this particle over the whole MD step by holding
                // f_langevin constant over delta_t: 0.5*mass*((v+f*dt/mass)^2 - v^2), expanded.
                // The dot(f,v)*dt term alone (used previously) is zero-mean for an additive random
                // force and misses the dominant, always-positive fluctuation/self-heating term
                // 0.5*dt^2*|f|^2/mass -- exactly the missing piece that made total_electronic_energy
                // fail to show LAMMPS's real, smooth downward drift.
                const double mass = get_mass( j, atom_type, masses.data(), has_type_field );
                const double work = delta_t * dot(f_langevin,v) + 0.5 * delta_t * delta_t * dot(f_langevin,f_langevin) / mass;

                if( is_local_cell )
                {
                  fx[j] += f_langevin.x;
                  fy[j] += f_langevin.y;
                  fz[j] += f_langevin.z;
                }

                for(int k=0;k<nbh_count;k++)
                {
#                 pragma omp atomic update
                  cell_energy_transfer[ nbh_idx[k] ] += nbh_w[k] * work;
                }
              }
            }
          }
          GRID_OMP_FOR_END
        }

        // LAMMPS fix-ttm's transfer_energy (f_twotemp[2]): total energy transferred from
        // electrons to ions over this whole MD step -- cell_energy_transfer already holds it
        // (fixed, computed once above), just sum it over this rank's own (non-ghost) cells.
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

      // inspired from https://en.wikipedia.org/wiki/Discrete_Laplace_operator#Finite_differences
      static constexpr double Lap27Norm = 26.0;
      static constexpr double Lap27_compact[4] = { -88/Lap27Norm , 6/Lap27Norm, 3/Lap27Norm, 2/Lap27Norm };

#     ifndef NDEBUG
      static constexpr double Lap27 [3][3][3] = {
         { { 2, 3, 2 } ,
           { 3, 6, 3 } ,
           { 2, 3, 2 } } ,
         { { 3, 6, 3 } ,
           { 6,-88,6 } ,
           { 3, 6, 3 } } ,
         { { 2, 3, 2 } ,
           { 3, 6, 3 } ,
           { 2, 3, 2 } }
        };
#     endif

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

      for(long te_istep=0; te_istep<num_inner_timesteps; ++te_istep)
      {

      double sum_Te = 0.0;
      double norm_dTe = 0.0;
      double entropy_Te = 0.0;

      // compute Te laplace operator
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,cell_i,cell_loc, schedule(static) reduction(+:sum_Te,norm_dTe,entropy_Te) )
        {
          //const IJK cell_loc = loc + gl;
          //const size_t cell_i = grid_ijk_to_index( dims , cell_loc );

          for(int ck=0;ck<subdiv;ck++)
          for(int cj=0;cj<subdiv;cj++)
          for(int ci=0;ci<subdiv;ci++)
          {
            IJK sc { ci, cj, ck };
            size_t j = cell_i*n_subcells +  grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , sc );

            // Te discrete Laplacian operator
            double L_Te = 0.0;
            for(int nk=-1;nk<=1;nk++)
            for(int nj=-1;nj<=1;nj++)
            for(int ni=-1;ni<=1;ni++)
            {
              IJK nbh { ni, nj, nk };
              IJK nbh_cell_loc;
              IJK nbh_subcell_loc;
              gcv_subcell_neighbor( cell_loc, sc, subdiv, nbh, nbh_cell_loc, nbh_subcell_loc );
              if( grid->contains(nbh_cell_loc) )
              {
                ssize_t nbh_cell_i = grid_ijk_to_index( dims , nbh_cell_loc );
                ssize_t nbh_subcell_i = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , nbh_subcell_loc );
                assert( nbh_cell_i>=0 && nbh_cell_i<n_cells );
                assert( nbh_subcell_i>=0 && nbh_subcell_i<n_subcells );
                size_t nbh_j = nbh_cell_i*n_subcells + nbh_subcell_i;

                int lap_compact_index = std::abs(ni) + std::abs(nj) + std::abs(nk);
                assert( Lap27_compact[lap_compact_index] == Lap27[ni+1][nj+1][nk+1]/Lap27Norm );
//              L_Te += Te[nbh_j] * Lap27[ni+1][nj+1][nk+1] / Lap27Norm;
                L_Te += Te[nbh_j] * Lap27_compact[lap_compact_index];
              }
            }
            // normalize by h^2: raw stencil sum = h^2 * laplacian(Te) + O(h^4)
            cell_L_Te[j] = L_Te / (subcell_size*subcell_size);
            if( ! grid->is_ghost_cell(cell_loc) )
            {
              sum_Te += Te[j];
              const double dTe = ( (Te_cond*cell_L_Te[j]) / Ce_rho_e ) * inner_dt;
              norm_dTe += std::fabs( dTe );
              entropy_Te += Te[j] * std::log(Te[j]) * subcell_volume ;
            }
          }
        }
        GRID_OMP_FOR_END
      }

      // check global absolute electronic energy variation (from dissipation)
//#     ifndef NDEBUG
      {
        double tmp[3] = { sum_Te, norm_dTe, entropy_Te };
        MPI_Allreduce(MPI_IN_PLACE,tmp,3,MPI_DOUBLE,MPI_SUM,*mpi);
        sum_Te = tmp[0];
        norm_dTe = tmp[1];
        entropy_Te = tmp[2];
        ldbg <<"sum_Te="<<sum_Te <<" norm_dTe/sum_Te="<<norm_dTe/sum_Te<<" , entropy_Te="<<entropy_Te << std::endl;
      }
//#     endif


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
          const Vec3d cell_origin = grid->cell_position( cell_loc );

          for(int ck=0;ck<subdiv;ck++)
          for(int cj=0;cj<subdiv;cj++)
          for(int ci=0;ci<subdiv;ci++)
          {
            IJK sc { ci, cj, ck };
            Vec3d scr { ci+0.5, cj+0.5, ck+0.5 };
            const size_t scindex = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , sc );
            const size_t idx_ti = cell_i*n_subcells + scindex ;
            const size_t idx_te = cell_i*cell_te_data.m_stride + scindex;

            // source terms
            const Vec3d center = xform * ( cell_origin + scr*subcell_size );
            const double Si = ti_source_func( center, *physical_time );
            const double Se = te_source_func( center, *physical_time );

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


  private:
    static inline double get_mass( unsigned int j, const uint8_t* type_ptr, const double* masses, std::true_type )
    {
      return masses[type_ptr[j]];
    }

    static inline double get_mass( unsigned int j, const uint8_t* type_ptr, const double* masses, std::false_type )
    {
      return masses[0];
    }

    static inline void localize_subcell( const Vec3d& r, double cell_size, double sub_cellsize, ssize_t subdiv, IJK& cell_loc, IJK& subcell_loc )
    {
      cell_loc = make_ijk( r / cell_size );
      Vec3d ro = r - (cell_loc*cell_size);
      subcell_loc = vclamp( make_ijk(ro / sub_cellsize) , 0 , subdiv-1 );
    }

    // @return how much of this particle contributes to region cell_box.
    // sum of contributions for all disjoint cell_box paving the domain is guaranteed to be 1.0
    static inline double particle_smoothing(const Vec3d& r, double sp_size, const AABB& cell_box)
    {
      AABB contrib_box = { r - sp_size*0.5 , r + sp_size*0.5 };
      AABB sub_contrib_box = intersection( contrib_box , cell_box );
      double w = 0.0;
      if( ! is_empty(sub_contrib_box) ) { w = bounds_volume(sub_contrib_box) / (sp_size*sp_size*sp_size); }
      assert( w>=0. && w<=(1.0+1.e-13) );
      return w;
    }

  };

  template<class GridT> using IonicElectronicHeatTransferTmpl = IonicElectronicHeatTransfer<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(ionic_electronic_heat_transfer)
  {
   OperatorNodeFactory::instance()->register_factory("ionic_eletronic_heat_transfer", make_grid_variant_operator< IonicElectronicHeatTransferTmpl > );
  }

}
