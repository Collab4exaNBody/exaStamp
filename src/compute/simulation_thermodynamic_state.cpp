#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/log.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exaStamp/compute/thermodynamic_state.h>
#include <exanb/core/grid_fields.h>
#include <onika/math/basic_types.h>

#include <mpi.h>
#include <cstring>

#include <onika/omp/ompt_interface.h>

namespace exaStamp
{

  // =================== utility functions ==========================

  // get particle mass from its type. assume the type index is 0 if particle hasn't type field
  template<bool has_atom_type> static double get_particle_mass(const ParticleSpecies&, const uint8_t* __restrict__, size_t);
  template<> inline double get_particle_mass<true>( const ParticleSpecies& species , const uint8_t* __restrict__ atom_type, size_t particle)
  {
    uint8_t t = atom_type[particle];
    assert( t>=0 && t<species.size() );
    return species[ t ].m_mass;
  }
  template<> inline double get_particle_mass<false>( const ParticleSpecies& species , const uint8_t* __restrict__ , size_t )
  {
    assert( species.size() >= 1 );
    return species[0].m_mass;
  }

  // ================== Thermodynamic state compute operator ======================

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_ax, field::_ay, field::_az, field::_vx, field::_vy, field::_vz >
    >
  struct ThermodynamicStateNode : public OperatorNode
  {
    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;

    // compile time constant indicating if grid has type virial
    using has_virial_field_t = typename GridT::CellParticles::template HasField < field::_virial > ;
    static constexpr bool has_virial_field = has_virial_field_t::value;

    using has_ep_field_t = typename GridT::CellParticles::template HasField < field::_ep > ;
    static constexpr bool has_ep_field = has_ep_field_t::value;
  
    ADD_SLOT( MPI_Comm           , mpi                 , INPUT , MPI_COMM_WORLD);
    ADD_SLOT( GridT              , grid                , INPUT , REQUIRED);
    ADD_SLOT( Domain             , domain              , INPUT , REQUIRED);
    ADD_SLOT( ParticleSpecies    , species             , INPUT , REQUIRED);
    ADD_SLOT( double             , potential_energy_shift , INPUT , 0.0 );
    ADD_SLOT( ThermodynamicState , thermodynamic_state , OUTPUT );

    inline void execute () override final
    {
      MPI_Comm comm = *mpi;
      GridT& grid = *(this->grid);
      ParticleSpecies& species = *(this->species);
      ThermodynamicState& sim_info = *thermodynamic_state;

      auto cells = grid.cells();
      IJK dims = grid.dimension();
      size_t ghost_layers = grid.ghost_layers();
      IJK dims_no_ghost = dims - (2*ghost_layers);

      Mat3d virial; // constructs itself with 0s
      Mat3d kinetic_tensor; // constructs itself with 0s
      Vec3d momentum;  // constructs itself with 0s
      Vec3d kinetic_energy;  // constructs itself with 0s
      double potential_energy = 0.;
      double mass = 0.;
      size_t total_particles = 0;
      
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghosts, reduction(+:mass,momentum,potential_energy,kinetic_energy,virial,kinetic_tensor,total_particles) )
        {
          IJK loc = loc_no_ghosts + ghost_layers;
          size_t cell_i = grid_ijk_to_index(dims,loc);

          const double* __restrict__ vx = cells[cell_i][field::vx];
          const double* __restrict__ vy = cells[cell_i][field::vy];
          const double* __restrict__ vz = cells[cell_i][field::vz];
          const double* __restrict__ ep = cells[cell_i].field_pointer_or_null(field::ep);
          const uint8_t* __restrict__ atom_type = cells[cell_i].field_pointer_or_null(field::type);
          const Mat3d* __restrict__ vir = cells[cell_i].field_pointer_or_null(field::virial);

          Mat3d local_virial;
          Mat3d local_kinetic_tensor;	  
          Vec3d local_momentum = {0.,0.,0.};
          Vec3d local_kinetic_ernergy = {0.,0.,0.};
          double local_potential_energy = 0.;
          double local_mass = 0.;
          size_t n = cells[cell_i].size();

#         pragma omp simd reduction(+:local_mass,local_momentum,local_potential_energy,local_kinetic_ernergy,local_virial,local_kinetic_tensor)
          for(size_t j=0;j<n;j++)
          {
            double mass = get_particle_mass<has_type_field>( species, atom_type, j ) ;
            Vec3d v { vx[j], vy[j], vz[j] };
            local_mass += mass;
            local_momentum += v * mass;
            local_kinetic_ernergy += v * v * mass; // * 0.5 later
            if constexpr (has_ep_field) { local_potential_energy += ep[j]; }
	    if constexpr (has_virial_field) { local_virial += vir[j]; }
            local_kinetic_tensor += (tensor(v,v) * mass);
          }

          mass += local_mass;
          momentum += local_momentum;
          potential_energy += local_potential_energy;
          kinetic_energy += local_kinetic_ernergy;
          virial += local_virial;
          kinetic_tensor += local_kinetic_tensor;
          total_particles += n;
        }
        GRID_OMP_FOR_END
      }

      // normalization
      kinetic_energy *= 0.5; // later is here
      Mat3d ke_tensor = 0.5 * kinetic_tensor;

      // reduce partial sums and share the result
      {
        // tmp size = 3*3 + 3 + 3 + 1 + 1 + 1 = 18
        double tmp[27] = {
          virial.m11, virial.m12, virial.m13, virial.m21, virial.m22, virial.m23,  virial.m31, virial.m32, virial.m33,
          ke_tensor.m11, ke_tensor.m12, ke_tensor.m13, ke_tensor.m21, ke_tensor.m22, ke_tensor.m23,  ke_tensor.m31, ke_tensor.m32, ke_tensor.m33, 	  
          momentum.x, momentum.y, momentum.z,
          kinetic_energy.x, kinetic_energy.y, kinetic_energy.z,
          potential_energy,
          mass,
          static_cast<double>(total_particles) };
        assert( tmp[26] == total_particles );
        MPI_Allreduce(MPI_IN_PLACE,tmp,27,MPI_DOUBLE,MPI_SUM,comm);
        virial.m11 = tmp[0];
        virial.m12 = tmp[1];
        virial.m13 = tmp[2];
        virial.m21 = tmp[3];
        virial.m22 = tmp[4];
        virial.m23 = tmp[5];
        virial.m31 = tmp[6];
        virial.m32 = tmp[7];
        virial.m33 = tmp[8];
        ke_tensor.m11 = tmp[9];
        ke_tensor.m12 = tmp[10];
        ke_tensor.m13 = tmp[11];
        ke_tensor.m21 = tmp[12];
        ke_tensor.m22 = tmp[13];
        ke_tensor.m23 = tmp[14];
        ke_tensor.m31 = tmp[15];
        ke_tensor.m32 = tmp[16];
        ke_tensor.m33 = tmp[17];	
        momentum.x = tmp[18];
        momentum.y = tmp[19];
        momentum.z = tmp[20];
        kinetic_energy.x = tmp[21];
        kinetic_energy.y = tmp[22];
        kinetic_energy.z = tmp[23];
        potential_energy = tmp[24];
        mass = tmp[25];
        total_particles = tmp[26];
      }

      Vec3d temperature = 2. * ( kinetic_energy - 0.5 * momentum * momentum / mass );
      Mat3d ke_test     = 2. * (      ke_tensor - 0.5 * tensor(momentum,momentum) / mass);

      Vec3d virdiag = { virial.m11 , virial.m22, virial.m33 };
      Vec3d virdev  = { virial.m12 , virial.m13, virial.m23 };
      Vec3d ke_tensor_diag = { ke_test.m11 , ke_test.m22, ke_test.m33 };
      Vec3d ke_tensor_dev  = { ke_test.m12 , ke_test.m13, ke_test.m23 };
      
      // Volume
      double volume = 1.0;
      if( ! domain->xform_is_identity() )
      {
        Mat3d mat = domain->xform();
        Vec3d a { mat.m11, mat.m21, mat.m31 };
        Vec3d b { mat.m12, mat.m22, mat.m32 };
        Vec3d c { mat.m13, mat.m23, mat.m33 };
        volume = dot( cross(a,b) , c );
      }
      volume *= bounds_volume( domain->bounds() );
      
      Vec3d phydro = ( ke_tensor_diag + virdiag ) / volume;
      Vec3d pdev   = (  ke_tensor_dev +  virdev ) / volume;

      // write results to output
      sim_info.set_virial( virial );
      sim_info.set_ke_tensor( ke_test );
      sim_info.set_pressure( phydro );
      sim_info.set_deviator( pdev );
      sim_info.set_kinetic_energy( kinetic_energy );
      sim_info.set_temperature( temperature );
      sim_info.set_kinetic_momentum( momentum );
      sim_info.set_potential_energy( potential_energy + (*potential_energy_shift) );
      sim_info.set_internal_energy( 0. );
      sim_info.set_chemical_energy( 0. );
      sim_info.set_mass( mass );
      sim_info.set_volume( volume );
      sim_info.set_particle_count( total_particles );
    }
  };
    
  template<class SomeT> using ThermodynamicStateNodeTmpl = ThermodynamicStateNode<SomeT>;
    
  // === register factories ===  
  ONIKA_AUTORUN_INIT(simulation_thermodynamic_state)
  {
   OperatorNodeFactory::instance()->register_factory( "simulation_thermodynamic_state", make_grid_variant_operator< ThermodynamicStateNodeTmpl > );
  }

}

