#include <exanb/core/grid.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <exanb/core/cpp_utils.h>

#include <mpi.h>

#include <algorithm>
#include <limits>

namespace exaStamp
{

  /* This component computes minimum and maximum of particle charge */
  template<class GridT>
  class ChargeMinMax : public OperatorNode
  {      
    // compile time constants
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;

    using has_charge_field_t = typename GridT::CellParticles::template HasField < field::_charge > ;
    static constexpr bool has_charge_field = has_charge_field_t::value;

    // ========= I/O slots =======================
    ADD_SLOT( GridT, grid, INPUT, REQUIRED );
    ADD_SLOT( ParticleSpecies, species , INPUT , REQUIRED );
    ADD_SLOT( MPI_Comm  , mpi         , INPUT , MPI_COMM_WORLD );

  public:
    // Operator execution
    inline void execute () override final
    {      
      size_t n_cells = grid->number_of_cells();
      if( n_cells==0 )
      {
        return ;
      }
      
      if( !has_charge_field && !has_type_field )
      {
        lerr << "Cannot get particles' charge" << std::endl;
        return;
      }
      
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      
      double max_charge = std::numeric_limits<double>::lowest();
      double min_charge = std::numeric_limits<double>::max();
      
#     pragma omp parallel
      {        
        GRID_OMP_FOR_BEGIN(dims,i,loc, schedule(dynamic) reduction(min:min_charge) reduction(max:max_charge) )
        {
          const auto* __restrict__ particle_charge = cells[i].field_pointer_or_null(field::charge);
          const auto* __restrict__ particle_type = cells[i].field_pointer_or_null(field::charge);
          size_t n = cells[i].size();
          for(size_t j=0;j<n;j++)          
          {
            double charge = 0.;
            if( has_charge_field )
            {
              charge = particle_charge[ j ];
            }
            else if( has_type_field )
            {
              charge = (*species)[ particle_type[j] ].m_charge;
            }
            min_charge = std::min( min_charge , charge );
            max_charge = std::max( max_charge , charge );
          }
        }
        GRID_OMP_FOR_END

      }

      MPI_Allreduce(MPI_IN_PLACE,&min_charge,1,MPI_DOUBLE,MPI_MIN,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&max_charge,1,MPI_DOUBLE,MPI_MAX,*mpi);
      
      lout << "charges in range ["<<min_charge<<";"<<max_charge<<"]"<<std::endl;
    }

  std::vector< double > m_type_charge;
  
  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(charge_min_max)
  {  
    OperatorNodeFactory::instance()->register_factory( "charge_min_max" , make_grid_variant_operator< ChargeMinMax > );
  }

}


