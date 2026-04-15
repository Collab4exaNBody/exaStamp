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

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/log.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <onika/math/basic_types_stream.h>
#include <onika/physics/constants.h>
#include <onika/parallel/random.h>

#include <mpi.h>
#include <cstring>

#include <exanb/grid_cell_particles/particle_region.h>

namespace exaStamp
{
  // ================== Thermodynamic state compute operator ======================

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_vx, field::_vy, field::_vz >
    >
  class InitTemperatureNewNode : public OperatorNode
  {
    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;

    using has_field_id_t = typename GridT:: template HasField < field::_id >;
    static constexpr bool has_field_id = has_field_id_t::value;

    ADD_SLOT( MPI_Comm           , mpi                  , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( ParticleSpecies    , species              , INPUT , REQUIRED );
    ADD_SLOT( double             , T                    , INPUT , 300.0 );        // target temperature in Kelvins
    ADD_SLOT( bool               , override_velocities  , INPUT , true );         // overwrite velocities with random noise, then rescale to T
    ADD_SLOT( bool               , add_velocities       , INPUT , false );        // draw random delta-v at temperature T and ADD it to existing velocities
    ADD_SLOT( bool               , scale_velocities     , INPUT , false );        // rescale existing velocities so that their current temperature maps to T
    ADD_SLOT( bool               , deterministic_noise  , INPUT , false );        // use per-cell deterministic seed (single-threaded)
    ADD_SLOT( bool               , zero_linear_momentum , INPUT , true  );        // remove net linear momentum from selected particles after rescaling
    ADD_SLOT( std::string        , distribution         , INPUT , "gaussian" );   // velocity distribution: "gaussian" or "uniform"
    ADD_SLOT( GridT              , grid                 , INPUT_OUTPUT );
    ADD_SLOT( Domain             , domain               , INPUT );

    ADD_SLOT( ParticleRegions   , particle_regions , INPUT , OPTIONAL );
    ADD_SLOT( ParticleRegionCSG , region           , INPUT , OPTIONAL );
    ADD_SLOT( std::string       , atom_type        , INPUT , OPTIONAL );

  public:
    inline void execute () override final
    {
      // Conversion factor: (kB in SI / atomicMass in SI) gives [m²/s² / K]
      // We work in LAMMPS-style metal units where velocities are in Å/ps and
      // masses in g/mol, hence the factor 1e4 (see unit analysis).
      static constexpr double conv_temperature =
          1.e4 * onika::physics::atomicMass / onika::physics::boltzmann;

      // ---- validate mode flags ----
      const bool genvel  = *override_velocities;
      const bool addvel  = *add_velocities;
      const bool scalevel = *scale_velocities;
      const bool zero_p = *zero_linear_momentum;

      if( (genvel && addvel) || (genvel && scalevel) || (addvel && scalevel) )
      {
        fatal_error() << "init_temperature: 'override_velocities' and 'add_velocities' "
                         "are mutually exclusive (at most one may be true)" << std::endl;
      }

      // ---- validate distribution parameter ----
      const std::string& distrib = *distribution;
      if( distrib != "gaussian" && distrib != "uniform" )
      {
        fatal_error() << "init_temperature: unknown distribution '" << distrib
                      << "', expected 'gaussian' or 'uniform'" << std::endl;
      }

      // ---- build particle region selection ----
      ParticleRegionCSGShallowCopy prcsg;
      if( region.has_value() && !particle_regions.has_value() )
      {
        fatal_error() << "region is defined, but particle_regions has no value" << std::endl;
      }
      if( region.has_value() && region->m_nb_operands == 0 )
      {
        region->build_from_expression_string( particle_regions->data() , particle_regions->size() );
      }
      if( region.has_value() ) prcsg = *region;

      // ---- resolve atom type filter ----
      int target_type = -1;  // -1 means "all types"
      if( atom_type.has_value() )
      {
        for( unsigned int i = 0; i < species->size(); i++ )
        {
          if( species->at(i).name() == *atom_type ) { target_type = i; break; }
        }
        if( target_type == -1 )
        {
          lerr << "Warning: atom_type '" << *atom_type
               << "' not found, applying temperature to all atoms" << std::endl;
        }
      }

      // ---- cache per-species masses ----
      double atom_masses[ species->size() ];
      for( unsigned int i = 0; i < species->size(); i++ )
      {
        atom_masses[i] = species->at(i).m_mass;
      }

      MPI_Comm comm = *mpi;

      auto cells  = grid->cells();
      IJK  dims   = grid->dimension();
      const ssize_t gl    = grid->ghost_layers();
      const IJK    gstart = { gl, gl, gl };

      ldbg << "target temperature = "    << *T     << " K"
           << "  distribution = "        << distrib
           << "  override_velocities = " << genvel
           << "  add_velocities = "      << addvel
           << "  scale_velocities = "     << scalevel
           << "  zero_linear_momentum = "<< zero_p << std::endl;

      // ---- nthreads: limit to 1 when deterministic to guarantee reproducibility ----
      const int nthreads = *deterministic_noise ? 1 : omp_get_max_threads();

      // =========================================================
      // Branch: add_velocities
      //
      // Draw unit-variance random velocities, scale each component
      // so that <½ m dv²> = ½ kB T per degree of freedom, i.e.
      //   sigma = sqrt( T / (conv_temperature * m) )
      // then immediately add to the existing velocity.
      // No rescaling of the existing velocities takes place.
      // =========================================================
      if( addvel )
      {
#       pragma omp parallel num_threads(nthreads)
        {
          std::mt19937_64 det_re;
          std::mt19937_64& re = *deterministic_noise ? det_re : onika::parallel::random_engine();
          std::normal_distribution<double>       gauss_rand( 0., 1. );
          std::uniform_real_distribution<double> unif_rand ( -std::sqrt(3.), std::sqrt(3.) );
          // uniform on [-sqrt(3), sqrt(3)] has unit variance, matching the gaussian sigma

          GRID_OMP_FOR_BEGIN( dims - 2*gl, _, loc, schedule(dynamic) )
          {
            const auto i = grid_ijk_to_index( dims, loc + gstart );

            if( *deterministic_noise )
            {
              const auto domain_cell_idx =
                  grid_ijk_to_index( grid->offset() + loc + gstart, dims );
              det_re.seed( static_cast<uint64_t>(domain_cell_idx) * 1023ULL );
            }

            const auto* __restrict__ rx    = cells[i][field::rx];
            const auto* __restrict__ ry    = cells[i][field::ry];
            const auto* __restrict__ rz    = cells[i][field::rz];
            auto* __restrict__       vx    = cells[i][field::vx];
            auto* __restrict__       vy    = cells[i][field::vy];
            auto* __restrict__       vz    = cells[i][field::vz];
            const auto* __restrict__ atype = cells[i].field_pointer_or_null( field::type );
            const auto* __restrict__ ids   = cells[i].field_pointer_or_null( field::id   );

            const size_t cell_nb_particles = cells[i].size();

            for( size_t j = 0; j < cell_nb_particles; j++ )
            {
              int type = 0;
              if constexpr (has_type_field) { type = atype[j]; }
              if( target_type != -1 && type != target_type ) continue;

              uint64_t id = 0;
              if constexpr (has_field_id) { id = ids[j]; }
              const Vec3d r { rx[j], ry[j], rz[j] };
              if( !prcsg.contains(r, id) ) continue;

              // sigma such that <½ m sigma² xi²> = ½ kB T  =>  sigma = sqrt(T / (conv_temperature * m))
              // (xi is a unit-variance deviate)
              const double sigma = std::sqrt( (*T) / ( conv_temperature * atom_masses[type] ) );

              const double xi_x = ( distrib == "gaussian" ) ? gauss_rand(re) : unif_rand(re);
              const double xi_y = ( distrib == "gaussian" ) ? gauss_rand(re) : unif_rand(re);
              const double xi_z = ( distrib == "gaussian" ) ? gauss_rand(re) : unif_rand(re);

              vx[j] += sigma * xi_x;
              vy[j] += sigma * xi_y;
              vz[j] += sigma * xi_z;
            }
          }
          GRID_OMP_FOR_END
        }
        // add_velocities does not rescale and intentionally does not touch
        // zero_linear_momentum: the caller controls the existing momentum.
        return;
      }

      // =========================================================
      // Pass 1 – optionally draw new velocities, accumulate stats
      // =========================================================
      Vec3d  momentum        = {0., 0., 0.};
      Vec3d  kinetic_energy  = {0., 0., 0.};
      double total_mass      = 0.;
      size_t total_particles = 0;

#     pragma omp parallel num_threads(nthreads)
      {
        std::mt19937_64 det_re;
        std::mt19937_64& re = *deterministic_noise ? det_re : onika::parallel::random_engine();
        std::normal_distribution<double>       gauss_rand( 0., 1. );
        std::uniform_real_distribution<double> unif_rand ( -1., 1. );

        GRID_OMP_FOR_BEGIN( dims - 2*gl, _, loc,
            reduction(+: momentum, kinetic_energy, total_mass, total_particles)
            schedule(dynamic) )
        {
          const auto i = grid_ijk_to_index( dims, loc + gstart );

          if( *deterministic_noise )
          {
            const auto domain_cell_idx =
                grid_ijk_to_index( grid->offset() + loc + gstart, dims );
            det_re.seed( static_cast<uint64_t>(domain_cell_idx) * 1023ULL );
          }

          const auto* __restrict__ rx    = cells[i][field::rx];
          const auto* __restrict__ ry    = cells[i][field::ry];
          const auto* __restrict__ rz    = cells[i][field::rz];
          auto* __restrict__       vx    = cells[i][field::vx];
          auto* __restrict__       vy    = cells[i][field::vy];
          auto* __restrict__       vz    = cells[i][field::vz];
          const auto* __restrict__ atype = cells[i].field_pointer_or_null( field::type );
          const auto* __restrict__ ids   = cells[i].field_pointer_or_null( field::id   );

          Vec3d  local_momentum       = {0., 0., 0.};
          Vec3d  local_kinetic_energy = {0., 0., 0.};
          double local_mass           = 0.;
          size_t n                    = 0;

          const size_t cell_nb_particles = cells[i].size();

          for( size_t j = 0; j < cell_nb_particles; j++ )
          {
            int type = 0;
            if constexpr (has_type_field) { type = atype[j]; }
            if( target_type != -1 && type != target_type ) continue;

            uint64_t id = 0;
            if constexpr (has_field_id) { id = ids[j]; }
            const Vec3d r { rx[j], ry[j], rz[j] };
            if( !prcsg.contains(r, id) ) continue;

            if( genvel )
            {
              if( distrib == "gaussian" )
              {
                vx[j] = gauss_rand(re);
                vy[j] = gauss_rand(re);
                vz[j] = gauss_rand(re);
              }
              else  // "uniform"
              {
                vx[j] = unif_rand(re);
                vy[j] = unif_rand(re);
                vz[j] = unif_rand(re);
              }
            }

            const double pmass = atom_masses[type];
            const Vec3d  v { vx[j], vy[j], vz[j] };

            local_mass           += pmass;
            local_momentum       += v * pmass;
            local_kinetic_energy += v * v * pmass;  // factor ½ applied below
            ++n;
          }

          momentum        += local_momentum;
          kinetic_energy  += local_kinetic_energy;
          total_mass      += local_mass;
          total_particles += n;
        }
        GRID_OMP_FOR_END
      }

      // ---- apply ½ to KE sum ----
      kinetic_energy *= 0.5;

      // ---- MPI reduction ----
      {
        double tmp[8] = {
          momentum.x,       momentum.y,       momentum.z,
          kinetic_energy.x, kinetic_energy.y, kinetic_energy.z,
          total_mass,
          static_cast<double>(total_particles)
        };
        MPI_Allreduce( MPI_IN_PLACE, tmp, 8, MPI_DOUBLE, MPI_SUM, comm );
        momentum.x       = tmp[0];  momentum.y       = tmp[1];  momentum.z       = tmp[2];
        kinetic_energy.x = tmp[3];  kinetic_energy.y = tmp[4];  kinetic_energy.z = tmp[5];
        total_mass       = tmp[6];
        total_particles  = static_cast<size_t>( tmp[7] );
      }

      if( total_particles == 0 )
      {
        lerr << "Warning: init_temperature found 0 selected particles – nothing to do." << std::endl;
        return;
      }

      // ---- centre-of-mass velocity ----
      const Vec3d v_cm = momentum / total_mass;

      // ---- current temperature (CM-corrected) ----
      // KE' = KE - ½ M v_cm²  removes the CM contribution.
      // Equipartition: T = 2 * KE' / (3N * kB).
      // The factor 2 is explicit because kinetic_energy already has ½ folded in.
      const Vec3d   ke_corrected = kinetic_energy - Vec3d{ 0.5 * total_mass } * v_cm * v_cm;
      const double  temp_current =
          ( 2.0 * conv_temperature * ( ke_corrected.x + ke_corrected.y + ke_corrected.z ) )
          / ( 3.0 * static_cast<double>(total_particles) );

      if( temp_current <= 0. )
      {
        if( scalevel )
        {
          lerr << "Warning: init_temperature – scale_velocities requested but current "
                  "temperature is 0 K (all velocities are zero). Scaling will have no effect." << std::endl;
        }
        else
        {
          lerr << "Warning: init_temperature – computed temperature is non-positive ("
               << temp_current << " K). Check that velocities are not all zero." << std::endl;
        }
        return;
      }

      const double vel_scale = std::sqrt( (*T) / temp_current );

      ldbg << "current temperature = " << temp_current << " K"
           << "  velocity scale factor = " << vel_scale << std::endl;

      // =========================================================
      // Pass 2 – remove CM momentum (optional) and rescale to T
      // Note: zero_linear_momentum is skipped in scale_velocities
      // mode because a pure rescaling preserves the existing CM
      // velocity proportionally, which is the intended behaviour.
      // =========================================================
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN( dims - 2*gl, _, loc, schedule(dynamic) )
        {
          const auto i = grid_ijk_to_index( dims, loc + gstart );

          const auto* __restrict__ rx    = cells[i][field::rx];
          const auto* __restrict__ ry    = cells[i][field::ry];
          const auto* __restrict__ rz    = cells[i][field::rz];
          auto* __restrict__       vx    = cells[i][field::vx];
          auto* __restrict__       vy    = cells[i][field::vy];
          auto* __restrict__       vz    = cells[i][field::vz];
          const auto* __restrict__ atype = cells[i].field_pointer_or_null( field::type );
          const auto* __restrict__ ids   = cells[i].field_pointer_or_null( field::id   );

          const size_t n = cells[i].size();

#         pragma omp simd
          for( size_t j = 0; j < n; j++ )
          {
            int type = 0;
            if constexpr (has_type_field) { type = atype[j]; }
            if( target_type != -1 && type != target_type ) continue;

            uint64_t id = 0;
            if constexpr (has_field_id) { id = ids[j]; }
            const Vec3d r { rx[j], ry[j], rz[j] };
            if( !prcsg.contains(r, id) ) continue;

            if( zero_p && !scalevel )
            {
              vx[j] -= v_cm.x;
              vy[j] -= v_cm.y;
              vz[j] -= v_cm.z;
            }

            vx[j] *= vel_scale;
            vy[j] *= vel_scale;
            vz[j] *= vel_scale;
          }
        }
        GRID_OMP_FOR_END
      }
    }

   inline void yaml_initialize(const YAML::Node& node) override final
   {
      YAML::Node tmp;
      if( node.IsScalar() )
      {
        tmp["T"] = node;
      }
      else { tmp = node; }
      this->OperatorNode::yaml_initialize( tmp );
   }

  };

  template<class GridT> using InitTemperatureNewTmpl = InitTemperatureNewNode<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(init_temperature_new)
  {
   OperatorNodeFactory::instance()->register_factory( "init_temperature_new", make_grid_variant_operator< InitTemperatureNewTmpl > );
  }

}
