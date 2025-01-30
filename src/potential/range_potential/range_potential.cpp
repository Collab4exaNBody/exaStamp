

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>

#include "range_potential_ext.h"
#include "range_potential_force_op.h"

#include <exanb/particle_neighbors/chunk_neighbors.h>

// this allows for parallel compilation of templated operator for each available field set


namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
    >
  class RangePotential : public OperatorNode
  {      
    // ========= I/O slots =======================
    // Range Neighbors potential parameters
    ADD_SLOT( double                , r0min          , INPUT , 0.0 );
    ADD_SLOT( double                , r0max          , INPUT , 0.0 );
    ADD_SLOT( double                , r1min          , INPUT , 0.0 );
    ADD_SLOT( double                , r1max          , INPUT , 0.0 );
    ADD_SLOT( double                , r2min          , INPUT , 0.0 );
    ADD_SLOT( double                , r2max          , INPUT , 0.0 );

    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0 );    
    ADD_SLOT( exanb::GridChunkNeighbors    , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                  , ghost             , INPUT , false );
    ADD_SLOT( GridT                 , grid              , INPUT_OUTPUT );
    ADD_SLOT( Domain                , domain            , INPUT , REQUIRED );

    // ========= Internal types =======================

    // cell particles array type
    using CellParticles = typename GridT::CellParticles;

    // compile time constant indicating if grid has virial field
    static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz , field::_charge >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz, field::_charge, field::_virial >;
    using ComputeFields = std::conditional_t< has_virial_field , ComputeFieldsWithVirial , ComputeFieldsWithoutVirial >;

    // shortcut to the Compute buffer used (and passed to functor) by compute_cell_particle_pairs
    using ComputeBuffer = ComputePairBuffer2<false,false,RangeNeighborsExtraStorage,RangeNeighborsClassifier>;

  public:
    // Operator execution
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      double rcut = std::max( std::max(*r0max,*r1max) , *r2max );
      *rcut_max = std::max( *rcut_max , rcut );
      
      size_t n_cells = grid->number_of_cells();

      // in this case, nothing to compute.
      // this is usefull case where compute_force is called at the very first to initialize rcut_max
      if( n_cells==0 )
      {
        return ;
      }
      
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };

      const double sq_r0min = (*r0min) * (*r0min) ;
      const double sq_r0max = (*r0max) * (*r0max) ;
      const double sq_r1min = (*r1min) * (*r1min) ;
      const double sq_r1max = (*r1max) * (*r1max) ;
      const double sq_r2min = (*r2min) * (*r2min) ;
      const double sq_r2max = (*r2max) * (*r2max) ;

      ComputePairBufferFactory< ComputeBuffer , RangeNeighborsInitComputeBuffer > force_buf = { { {sq_r0min,sq_r0max,sq_r1min,sq_r1max,sq_r2min,sq_r2max} } };
      RangePotentialForceOp force_op { *r0min , *r0max , *r1min , *r1max , *r2min , *r2max };

      if( domain->xform_is_identity() )
      {
        NullXForm cp_xform;
        auto optional = make_compute_pair_optional_args( nbh_it, ComputePairNullWeightIterator{} , cp_xform, cp_locks );
        compute_cell_particle_pairs( *grid, rcut, *ghost, optional, force_buf, force_op , ComputeFields{} , parallel_execution_context() );
      }
      else
      {
        LinearXForm cp_xform { domain->xform() };
        auto optional = make_compute_pair_optional_args( nbh_it, ComputePairNullWeightIterator{} , cp_xform, cp_locks );
        compute_cell_particle_pairs( *grid, rcut, *ghost, optional, force_buf, force_op , ComputeFields{} , parallel_execution_context() );
      }

    }

  };

  template<class GridT> using RangePotentialTmpl = RangePotential<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(range_potential)
  {  
    OperatorNodeFactory::instance()->register_factory( "range_potential" , make_grid_variant_operator<RangePotentialTmpl> );
  }

}


