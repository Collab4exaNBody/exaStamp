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

#include <exanb/core/quaternion_operators.h>
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
    class = AssertGridHasFields< GridT, field::_ep, field::_ax, field::_ay, field::_az, field::_vx, field::_vy, field::_vz, field::_angmom, field::_orient, field::_type >
    >
  struct ThermodynamicStateRigidmolNode : public OperatorNode
  {
    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;

    // compile time constant indicating if grid has type virial
    using has_virial_field_t = typename GridT::CellParticles::template HasField < field::_virial > ;
    static constexpr bool has_virial_field = has_virial_field_t::value;
  
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
      Vec3d momentum;  // constructs itself with 0s
      Vec3d kinetic_energy = {0., 0., 0.};  // constructs itself with 0s
      Vec3d rotational_energy = {0., 0., 0.};  // constructs itself with 0s
      Vec3d ndof = {0., 0., 0.}; //number of rotational degrees of freedom;
      double potential_energy = 0.;
      double masstotale = 0.;
      size_t total_particles = 0;
      
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghosts, reduction(+:masstotale,momentum,potential_energy,kinetic_energy,rotational_energy,ndof,virial,total_particles) )
        {
          IJK loc = loc_no_ghosts + ghost_layers;
          size_t cell_i = grid_ijk_to_index(dims,loc);

          const double* __restrict__ vx = cells[cell_i][field::vx];
          const double* __restrict__ vy = cells[cell_i][field::vy];
          const double* __restrict__ vz = cells[cell_i][field::vz];
          auto* __restrict__ angmom = cells[cell_i][field::angmom];
          auto* __restrict__ orient = cells[cell_i][field::orient];
          const double* __restrict__ ep = cells[cell_i][field::ep];
          const uint8_t* __restrict__ atom_type = cells[cell_i].field_pointer_or_null(field::type);
          const auto* __restrict__ type_atom = cells[cell_i][field::type];
          const Mat3d* __restrict__ vir = cells[cell_i].field_pointer_or_null(field::virial);

          Mat3d local_virial;
          Vec3d local_momentum = {0.,0.,0.};
          Vec3d local_kinetic_ernergy = {0.,0.,0.};
          Vec3d local_rotational_energy = {0.,0.,0.};
          double local_potential_energy = 0.;
          double local_mass = 0.;
          Vec3d local_ndof = {0.,0.,0.};
          size_t n = cells[cell_i].size();

#         pragma omp simd reduction(+:local_mass,local_momentum,local_potential_energy,local_kinetic_ernergy,local_rotational_energy,local_virial)
          for(size_t j=0;j<n;j++)
          {
            int t =type_atom[j];
            Vec3d minert = species[t].m_minert;
            double mass = get_particle_mass<has_type_field>( species, atom_type, j ) ;
            Vec3d v { vx[j], vy[j], vz[j] };
            local_mass += mass;
            local_momentum += v * mass;
            local_kinetic_ernergy += v * v * mass; // * 0.5 later

            //calcul du moment angulaire dans repere mobile
            Mat3d mat_lab_bf;
            mat_lab_bf.m11 = orient[j].w*orient[j].w + orient[j].x*orient[j].x - orient[j].y*orient[j].y - orient[j].z*orient[j].z;
            mat_lab_bf.m22 = orient[j].w*orient[j].w - orient[j].x*orient[j].x + orient[j].y*orient[j].y - orient[j].z*orient[j].z;
            mat_lab_bf.m33 = orient[j].w*orient[j].w - orient[j].x*orient[j].x - orient[j].y*orient[j].y + orient[j].z*orient[j].z;
            mat_lab_bf.m12 = 2.0 * (orient[j].x*orient[j].y + orient[j].w*orient[j].z ); 
            mat_lab_bf.m21 = 2.0 * (orient[j].x*orient[j].y - orient[j].w*orient[j].z );
            mat_lab_bf.m13 = 2.0 * (orient[j].x*orient[j].z - orient[j].w*orient[j].y );
            mat_lab_bf.m31 = 2.0 * (orient[j].x*orient[j].z + orient[j].w*orient[j].y );
            mat_lab_bf.m23 = 2.0 * (orient[j].y*orient[j].z + orient[j].w*orient[j].x );
            mat_lab_bf.m32 = 2.0 * (orient[j].y*orient[j].z - orient[j].w*orient[j].x );

            Vec3d angmom_m = mat_lab_bf * angmom[j]; 
            if (minert.x > 0.) {local_rotational_energy.x += angmom_m.x*angmom_m.x/minert.x; local_ndof.x += 1.;}
            if (minert.y > 0.) {local_rotational_energy.y += angmom_m.y*angmom_m.y/minert.y; local_ndof.y += 1.;}
            if (minert.z > 0.) {local_rotational_energy.z += angmom_m.z*angmom_m.z/minert.z; local_ndof.z += 1.;}
//         (omega=angmom/minert, Er=omega*omega*minert)
            local_potential_energy += ep[j];
            if constexpr (has_virial_field) { local_virial += vir[j]; }
          }

          masstotale += local_mass;
          momentum += local_momentum;
          potential_energy += local_potential_energy;
          kinetic_energy += local_kinetic_ernergy;
          rotational_energy += local_rotational_energy;
          virial += local_virial;
          total_particles += n;
          ndof += local_ndof;
        }
        GRID_OMP_FOR_END
      }

      // normalization
      kinetic_energy *= 0.5; // later is here
      rotational_energy *= 0.5;

      // reduce partial sums and share the result
      {
        // tmp size = 3*3 + 3 + 3 + 3 + 1 + 1 + 1 + 3 = 24
        double tmp[24] = {
          virial.m11, virial.m12, virial.m13, virial.m21, virial.m22, virial.m23,  virial.m31, virial.m32, virial.m33, 
          momentum.x, momentum.y, momentum.z,
          kinetic_energy.x, kinetic_energy.y, kinetic_energy.z,
          rotational_energy.x, rotational_energy.y, rotational_energy.z,
          potential_energy,
          masstotale,
          static_cast<double>(total_particles),
          ndof.x, ndof.y, ndof.z };

        assert( tmp[20] == total_particles );
        MPI_Allreduce(MPI_IN_PLACE,tmp,24,MPI_DOUBLE,MPI_SUM,comm);
        virial.m11 = tmp[0];
        virial.m12 = tmp[1];
        virial.m13 = tmp[2];
        virial.m21 = tmp[3];
        virial.m22 = tmp[4];
        virial.m23 = tmp[5];
        virial.m31 = tmp[6];
        virial.m32 = tmp[7];
        virial.m33 = tmp[8];
        momentum.x = tmp[9];
        momentum.y = tmp[10];
        momentum.z = tmp[11];
        kinetic_energy.x = tmp[12];
        kinetic_energy.y = tmp[13];
        kinetic_energy.z = tmp[14];
        rotational_energy.x = tmp[15];
        rotational_energy.y = tmp[16];
        rotational_energy.z = tmp[17];
        potential_energy = tmp[18];
        masstotale = tmp[19];
        total_particles = tmp[20];
        ndof.x = tmp[21];
        ndof.y = tmp[22];
        ndof.z = tmp[23];
      }

      double conv_temperature = 1.e4 * onika::physics::atomicMass / onika::physics::boltzmann ;
      Vec3d kinetic_temperature = 2. * ( kinetic_energy - 0.5 * momentum * momentum / masstotale) / total_particles ;

      Vec3d rotational_temperature = {0., 0., 0.};
      if (ndof.x > 0) { rotational_temperature.x = 2. * rotational_energy.x / ndof.x ;}
      if (ndof.y > 0) { rotational_temperature.y = 2. * rotational_energy.y / ndof.y ;}
      if (ndof.z > 0) { rotational_temperature.z = 2. * rotational_energy.z / ndof.z ;}

      Vec3d kinetic_energy_totale = kinetic_energy - 0.5 * momentum * momentum / masstotale + rotational_energy;

      Vec3d temperature = {0.,0.,0.};
      temperature.x = 2. * kinetic_energy_totale.x / ( total_particles + ndof.x);
      temperature.y = 2. * kinetic_energy_totale.y / ( total_particles + ndof.y);
      temperature.z = 2. * kinetic_energy_totale.z / ( total_particles + ndof.z);

      //lout<<"conv_temperature="<<conv_temperature<<" total_particles="<<total_particles<<" ndof="<<ndof<<std::endl;
      //lout<<" kinetic_temperature="<<kinetic_temperature*conv_temperature<<" rotational_energy="<<rotational_energy<<" rotational_temperature="<<rotational_temperature*conv_temperature<<" temperature="<<temperature*conv_temperature<<std::endl;

      Vec3d virdiag = { virial.m11 , virial.m22, virial.m33 };


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

      Vec3d pressure = ( temperature + virdiag ) / volume;

      // write results to output
      sim_info.set_virial( virial );
      sim_info.set_pressure( pressure );
      sim_info.set_kinetic_energy( kinetic_energy );
      sim_info.set_rotational_energy( rotational_energy );
      sim_info.set_temperature( temperature );
      sim_info.set_kinetic_temperature( kinetic_temperature );
      sim_info.set_rotational_temperature( rotational_temperature );
      sim_info.set_kinetic_momentum( momentum );
      sim_info.set_potential_energy( potential_energy + (*potential_energy_shift) );
      sim_info.set_internal_energy( 0. );
      sim_info.set_chemical_energy( 0. );
      sim_info.set_mass( masstotale );
      sim_info.set_ndof( ndof );
      sim_info.set_volume( volume );
      sim_info.set_particle_count( total_particles );
    }
  };
    
  template<class GridT> using ThermodynamicStateRigidmolNodeTmpl = ThermodynamicStateRigidmolNode<GridT>;
    
  // === register factories ===  
  ONIKA_AUTORUN_INIT(thermodynamic_state_rigidmol)
  {
   OperatorNodeFactory::instance()->register_factory( "simulation_thermodynamic_state_rigidmol", make_grid_variant_operator< ThermodynamicStateRigidmolNodeTmpl > );
  }

}

