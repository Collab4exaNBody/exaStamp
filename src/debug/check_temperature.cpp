#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/log.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/fields.h>
#include <exanb/core/basic_types_stream.h>

#include <mpi.h>
#include <cstring>

namespace exaStamp
{
  using namespace exanb;
  
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
    class = AssertGridHasFields< GridT, field::_vx, field::_vy, field::_vz >
    >
  struct CheckInitTemperatureNode : public OperatorNode
  {
    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;
 
    ADD_SLOT( MPI_Comm           , mpi                 , INPUT , REQUIRED);
    ADD_SLOT( ParticleSpecies    , species             , INPUT , REQUIRED);
    ADD_SLOT( double             , temperature         , INPUT , REQUIRED );
    ADD_SLOT( double             , tolerance           , INPUT , 0.01 ); // default 1% tolerance
    ADD_SLOT( GridT              , grid                , INPUT_OUTPUT );

    inline void execute () override final
    {
      MPI_Comm comm = *mpi;
      GridT& grid = *(this->grid);
      ParticleSpecies& species = *(this->species);

      auto cells = grid.cells();
      IJK dims = grid.dimension();
      size_t ghost_layers = grid.ghost_layers();
      IJK dims_no_ghost = dims - (2*ghost_layers);

      Vec3d momentum;  // constructs itself with 0s
      Vec3d kinetic_energy;  // constructs itself with 0s
      double total_mass = 0.;
      size_t total_particles = 0;

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghosts, reduction(+:momentum,kinetic_energy,total_mass,total_particles) schedule(dynamic) )
        {
          IJK loc = loc_no_ghosts + ghost_layers;
          size_t cell_i = grid_ijk_to_index(dims,loc);

          const auto* __restrict__ vx = cells[cell_i][field::vx];
          const auto* __restrict__ vy = cells[cell_i][field::vy];
          const auto* __restrict__ vz = cells[cell_i][field::vz];
          const uint8_t* __restrict__ atom_type = cells[cell_i].field_pointer_or_null(field::type);

          Vec3d local_momentum = {0.,0.,0.};
          Vec3d local_kinetic_ernergy = {0.,0.,0.};
          double local_mass = 0.;
          size_t n = cells[cell_i].size();

#         pragma omp simd reduction(+:local_momentum,local_kinetic_ernergy)
          for(size_t j=0;j<n;j++)
          {
            double pmass = get_particle_mass<has_type_field>( species, atom_type, j ) ;
            Vec3d v { vx[j], vy[j], vz[j] };
            local_mass += pmass;
            local_momentum += v * pmass;
            local_kinetic_ernergy += v * v * pmass; // * 0.5 later
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
      static constexpr double conv_temperature = 1.202712206318418424189076176845;
      Vec3d temp = 2. * ( kinetic_energy - 0.5 * momentum * momentum / total_mass );
      double temp_scale = ( conv_temperature * ( temp.x + temp.y + temp.z ) / 3. ) / total_particles;
  
      double min_temp = (*temperature)*(1.0-(*tolerance));
      double max_temp = (*temperature)*(1.0+(*tolerance));
      lout << "System temperature "<<temp_scale<<" K, target temperature range ["<<min_temp<<";"<<max_temp<<"]"<<std::endl;
      if( temp_scale < min_temp || temp_scale > max_temp )
      {
        std::abort();
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
    
  template<class GridT> using CheckTemperatureTmpl = CheckInitTemperatureNode<GridT>;
    
  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "check_temperature", make_grid_variant_operator< CheckTemperatureTmpl > );
  }

}

