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
  class InitTemperatureNode : public OperatorNode
  {
    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;
 
    using has_field_id_t = typename GridT:: template HasField < field::_id >;
    static constexpr bool has_field_id = has_field_id_t::value;

    ADD_SLOT( MPI_Comm           , mpi                 , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( ParticleSpecies    , species             , INPUT , REQUIRED );
    ADD_SLOT( double             , temperature         , INPUT , 300.0 ); // in Kelvins
    ADD_SLOT( bool               , generate_velocity   , INPUT , false ); // set to true to erase velocity and replace it with gaussian noise
    ADD_SLOT( bool               , deterministic_noise , INPUT , false );    
    ADD_SLOT( GridT              , grid                , INPUT_OUTPUT );
    ADD_SLOT( Domain , domain    , INPUT );

    ADD_SLOT( ParticleRegions   , particle_regions , INPUT , OPTIONAL );
    ADD_SLOT( ParticleRegionCSG , region           , INPUT , OPTIONAL );
    ADD_SLOT( std::string       , atom_type        , INPUT , OPTIONAL );

  public:
    inline void execute () override final
    {
      static constexpr double conv_temperature = 1.e4 * onika::physics::atomicMass / onika::physics::boltzmann ;

      // build particle selection object 'prcsg' from region expression, if any specified
      ParticleRegionCSGShallowCopy prcsg;
      if( region.has_value() && !particle_regions.has_value() )
      {
        fatal_error() << "region is defined, but particle_regions has no value" << std::endl;
      }        
      if( region.has_value() && region->m_nb_operands==0 )
      {
        region->build_from_expression_string( particle_regions->data() , particle_regions->size() );
      }
      if( region.has_value() ) prcsg = *region;

      // look for target type
      int target_type = -1;
      if( atom_type.has_value() )
      {
        for(unsigned int i=0;i<species->size();i++)
        {
          if( species->at(i).name() == *atom_type ) target_type = i;
        }
        if( target_type == -1 )
        {
          lerr << "Warning: atom_type '"<<*atom_type<<"' not found, apply temperature to all atoms"<<std::endl;
        }
      }

      double atom_masses [ species->size() ];
      for(unsigned int i=0;i<species->size();i++)
      {
        atom_masses[i] = species->at(i).m_mass;
      }

      MPI_Comm comm = *mpi;

      const bool det_noise = *deterministic_noise;
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      size_t ghost_layers = grid->ghost_layers();
      IJK dims_no_ghost = dims - (2*ghost_layers);

      IJK gstart { ghost_layers, ghost_layers, ghost_layers };
      const auto dom_dims = domain->grid_dimension();
      const auto dom_start = grid->offset();
      
      Vec3d momentum;  // constructs itself with 0s
      Vec3d kinetic_energy;  // constructs itself with 0s
      double total_mass = 0.;
      size_t total_particles = 0;

      const bool genvel = *generate_velocity;
      ldbg << "target temperature is "<< *temperature << std::endl;

#     pragma omp parallel
      {
        std::mt19937_64 det_re;
        std::mt19937_64 & re = det_noise ? det_re : onika::parallel::random_engine() ;
        
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghosts, reduction(+:momentum,kinetic_energy,total_mass,total_particles) schedule(dynamic) )
        {
          
          IJK loc = loc_no_ghosts + ghost_layers;
          size_t cell_i = grid_ijk_to_index(dims,loc);
          const auto domain_cell_idx = grid_ijk_to_index( dom_dims , loc + gstart + dom_start );
          det_re.seed( domain_cell_idx * 1023 );

          const auto* __restrict__ rx = cells[cell_i][field::rx];
          const auto* __restrict__ ry = cells[cell_i][field::ry];
          const auto* __restrict__ rz = cells[cell_i][field::rz];

          auto* __restrict__ vx = cells[cell_i][field::vx];
          auto* __restrict__ vy = cells[cell_i][field::vy];
          auto* __restrict__ vz = cells[cell_i][field::vz];
          
          const auto* __restrict__ atype = cells[cell_i].field_pointer_or_null(field::type);
          const auto* __restrict__ ids = cells[cell_i].field_pointer_or_null(field::id);

          Vec3d local_momentum = {0.,0.,0.};
          Vec3d local_kinetic_ernergy = {0.,0.,0.};
          double local_mass = 0.;
          const size_t cell_nb_particles = cells[cell_i].size();
          size_t n = 0;

          std::normal_distribution<double> gaussian(0.0,1.0);
          
//#         pragma omp simd reduction(+:local_momentum,local_kinetic_ernergy,n)
          for(size_t j=0;j<cell_nb_particles;j++)
          {
            int type = 0;
            if constexpr (has_type_field) { type = atype[j]; }
            if( target_type==-1 || type==target_type )
            {
              uint64_t id = 0;
              if constexpr ( has_field_id ) { id = ids[j]; }
              const Vec3d r { rx[j], ry[j], rz[j] };
              if( prcsg.contains(r,id) )
              {
                if( genvel )
                {
                  vx[j] = gaussian(re);
                  vy[j] = gaussian(re);
                  vz[j] = gaussian(re);
                }
                const double pmass = atom_masses[type];
                const Vec3d v { vx[j], vy[j], vz[j] };
                local_mass += pmass;
                local_momentum += v * pmass;
                local_kinetic_ernergy += v * v * pmass; // * 0.5 later
                ++n;
              }
            }
          }

          momentum += local_momentum;
          kinetic_energy += local_kinetic_ernergy;
          total_mass += local_mass;
          total_particles += n;          
        }
        GRID_OMP_FOR_END
      }

      // normalization
      kinetic_energy *= 0.5; // later is here

      // reduce partial sums and share the result
      {
        // tmp size = 3*3 + 3 + 3 + 1 + 1 + 1 = 18
        double tmp[8] = {
          momentum.x, momentum.y, momentum.z,
          kinetic_energy.x, kinetic_energy.y, kinetic_energy.z,
          total_mass,
          static_cast<double>(total_particles) };
        MPI_Allreduce(MPI_IN_PLACE,tmp,8,MPI_DOUBLE,MPI_SUM,comm);
        momentum.x = tmp[0];
        momentum.y = tmp[1];
        momentum.z = tmp[2];
        kinetic_energy.x = tmp[3];
        kinetic_energy.y = tmp[4];
        kinetic_energy.z = tmp[5];
        total_mass = tmp[6];
        total_particles = tmp[7];
      }

      // temperature
      Vec3d temp = 2. * ( kinetic_energy - 0.5 * momentum * momentum / total_mass );
      double temp_scale = ( conv_temperature * ( temp.x + temp.y + temp.z ) / 3. ) / total_particles;
      temp_scale = (*temperature) / temp_scale;
      double vel_scale = std::sqrt( temp_scale );
      Vec3d momentum_shift = momentum / total_mass;

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghosts, schedule(dynamic) )
        {
          IJK loc = loc_no_ghosts + ghost_layers;
          size_t cell_i = grid_ijk_to_index(dims,loc);

          const auto* __restrict__ rx = cells[cell_i][field::rx];
          const auto* __restrict__ ry = cells[cell_i][field::ry];
          const auto* __restrict__ rz = cells[cell_i][field::rz];

          auto* __restrict__ vx = cells[cell_i][field::vx];
          auto* __restrict__ vy = cells[cell_i][field::vy];
          auto* __restrict__ vz = cells[cell_i][field::vz];
          
          const auto* __restrict__ atype = cells[cell_i].field_pointer_or_null(field::type);
          const auto* __restrict__ ids = cells[cell_i].field_pointer_or_null(field::id);

          size_t n = cells[cell_i].size();

#         pragma omp simd
          for(size_t j=0;j<n;j++)
          {
            int type = 0;
            if constexpr (has_type_field) { type = atype[j]; }
            if( target_type==-1 || type==target_type )
            {
              uint64_t id = 0;
              if constexpr ( has_field_id ) { id = ids[j]; }
              const Vec3d r { rx[j], ry[j], rz[j] };
              if( prcsg.contains(r,id) )
              {
                const double pmass = atom_masses[type];
                vx[j] -= momentum_shift.x * pmass;
                vy[j] -= momentum_shift.y * pmass;
                vz[j] -= momentum_shift.z * pmass;
                vx[j] *= vel_scale;
                vy[j] *= vel_scale;
                vz[j] *= vel_scale;
              }
            }
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
        tmp["temperature"] = node;
      }
      else { tmp = node; }
      this->OperatorNode::yaml_initialize( tmp );
   }

  };
    
  template<class GridT> using InitTemperatureTmpl = InitTemperatureNode<GridT>;
    
  // === register factories ===  
  ONIKA_AUTORUN_INIT(init_temperature)
  {
   OperatorNodeFactory::instance()->register_factory( "init_temperature", make_grid_variant_operator< InitTemperatureTmpl > );
  }

}

