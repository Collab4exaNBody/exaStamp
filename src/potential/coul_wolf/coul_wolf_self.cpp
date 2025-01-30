

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exaStamp/potential/coul_wolf/coul_wolf.h>

namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;

  struct CoulWolfSelfFunctor
  {
    const ParticleSpecie * __restrict__ m_species = nullptr;
    const CoulWolfParms m_params = {};
    ONIKA_HOST_DEVICE_FUNC inline void operator () (double& ep, unsigned int type) const
    {
      const double C = m_species[type].m_charge;
      double e_self = -(m_params.e_shift / 2.0 + m_params.alpha / sqrt(M_PI) ) * C * C * m_params.qqrd2e;

      double _ep=0.;
      _ep = UnityConverterHelper::convert(e_self, "eV");
      ep += _ep;
    }
  };

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep , field::_type >
    >
  
  class CoulWolfSelf : public OperatorNode
  {      
    // ========= I/O slots =======================
    ADD_SLOT( CoulWolfParms    , parameters        , INPUT , REQUIRED );
    ADD_SLOT( bool             , ghost             , INPUT , false );
    ADD_SLOT( ParticleSpecies  , species           , INPUT , REQUIRED );
    ADD_SLOT( GridT            , grid              , INPUT_OUTPUT );

    // ========= Internal types =======================

    // cell particles array type
    using CellParticles = typename GridT::CellParticles;

    // attributes processed during computation
    static inline constexpr FieldSet< field::_ep, field::_type > compute_fields = {};

  public:
    // Operator execution
    inline void execute () override final
    {
      //      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      
      size_t n_cells = grid->number_of_cells();

      // in this case, nothing to compute.
      // this is usefull case where compute_force is called at the very first to initialize rcut_max
      if( n_cells==0 )
      {
        return ;
      }
            
      CoulWolfSelfFunctor func = { species->data() , *parameters };
      compute_cell_particles( *grid , *ghost , func , compute_fields , parallel_execution_context() );
    }

  };

  template<class GridT> using CoulWolfSelfTmpl = CoulWolfSelf<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(coul_wolf_self)
  {  
    OperatorNodeFactory::instance()->register_factory( "coul_wolf_self" , make_grid_variant_operator<CoulWolfSelfTmpl> );
  }

}


