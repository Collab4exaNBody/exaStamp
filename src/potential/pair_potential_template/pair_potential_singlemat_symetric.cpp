#pragma xstamp_grid_variant

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <exanb/core/compact_grid_pair_weights.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_pair_optional_args.h>
#include <exanb/compute/compute_cell_particle_pairs.h>

#include <iostream>
#include <type_traits>

#include "pair_potential_template.h"

// this allows for parallel compilation of templated operator for each available field set


// defines operator class name and string from potential name
#define OPERATOR_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_compute_force_symetric)
#define OPERATOR_NAME_STR USTAMP_STR(OPERATOR_NAME)

// #define PAIR_POT_DEBUG_VERSION 1
#ifdef PAIR_POT_DEBUG_VERSION
#define DEBUG_ADDITIONAL_FIELDS ,field::id
#define DEBUG_ADDITIONAL_PARAMETERS int64_t id,
#define DEBUG_ADDITIONAL_CODE if(id==0) { std::cout<<"id="<<id<<" : nnbh="<<n<<" ep="<<ep<<" _ep="<<_ep<< std::endl; }
#else
#define DEBUG_ADDITIONAL_FIELDS /**/
#define DEBUG_ADDITIONAL_PARAMETERS /**/
#define DEBUG_ADDITIONAL_CODE /**/
#endif

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT,
    class = AssertGridHasFields< GridT ,field::_fx ,field::_fy ,field::_fz DEBUG_ADDITIONAL_FIELDS >
    >
  class OPERATOR_NAME : public OperatorNode
  {
//    using ParticleLock = decltype( ComputePairOptionalLocks<false>{}[0][0] );
    using CellParticles = typename GridT::CellParticles;

    // attributes processed during computation
    static inline constexpr DefaultFields position_fields {};

    // compile time constant indicating if grid has virial field
    static inline constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;
    static inline constexpr bool has_ep_field = GridHasField<GridT,field::_ep>::value;
    // compute fields
    using _ComputeFields = FieldSet< field::_fx ,field::_fy ,field::_fz DEBUG_ADDITIONAL_FIELDS >;
    using _ComputeFieldsVirial    = FieldSet< field::_fx ,field::_fy ,field::_fz ,field::_virial DEBUG_ADDITIONAL_FIELDS >;
    using _ComputeFieldsEp = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz DEBUG_ADDITIONAL_FIELDS >;
    using _ComputeFieldsEpVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_virial DEBUG_ADDITIONAL_FIELDS >;
    using ComputeFields = std::conditional_t< has_ep_field ,
                            std::conditional_t< has_virial_field , _ComputeFieldsEpVirial , _ComputeFieldsEp > ,
                            std::conditional_t< has_virial_field , _ComputeFieldsVirial , _ComputeFields > >;
    static constexpr ComputeFields compute_fields {};

    static inline constexpr DefaultFields center_ro_fields {};
    static inline constexpr auto center_rw_fields = compute_fields;
    static inline constexpr DefaultFields nbh_ro_fields {};
    static inline constexpr auto nbh_rw_fields = compute_fields;

    // task graph execution parameters
    static inline constexpr std::integral_constant<size_t,2> grid_grain {};
    static inline constexpr std::integral_constant<size_t,2> stencil_scale {};

    // helper functor, with templated call operator (with/without weights)
    struct SymetricForceOp
    {
      USTAMP_POTENTIAL_PARAMS p;
      //PairPotentialParameters
      PairPotentialMinimalParameters pair_params;
      double ecut;

      // version with virial computation
      template<bool UseWeights, class CellsAccessorT>
      inline void operator ()
        (
        size_t n,
        ComputePairBuffer2<UseWeights,true>& tab,
        double& _ep,
        double& _fx,
        double& _fy,
        double& _fz,
        CellsAccessorT cells
        ) const
      {
        FakeMat3d fake_virial;
        NullGridParticleLocks fake_locks;
        this-> operator () ( n,tab,_ep,_fx,_fy,_fz, fake_virial, cells, fake_locks , fake_locks[0][0] );
      }

      template<bool UseWeights, class CellsAccessorT>
      inline void operator ()
        (
        size_t n,
        ComputePairBuffer2<UseWeights,true>& tab,
        double& _fx,
        double& _fy,
        double& _fz,
        CellsAccessorT cells
        ) const
      {
        double _ep = 0.0;
        FakeMat3d fake_virial;
        NullGridParticleLocks fake_locks;
        this-> operator () ( n,tab,_ep,_fx,_fy,_fz, fake_virial, cells, fake_locks , fake_locks[0][0] );
      }

      template<bool UseWeights, class CellsAccessorT, class GridCellLocksT, class ParticleLockT>
      inline void operator ()
        (
        size_t n,
        ComputePairBuffer2<UseWeights,true>& tab,
        double& _ep,
        double& _fx,
        double& _fy,
        double& _fz,
        CellsAccessorT cells,
        GridCellLocksT& cell_locks,
        ParticleLockT& particle_lock
        ) const
      {
        FakeMat3d fake_virial;
        this-> operator () ( n,tab,_ep,_fx,_fy,_fz, fake_virial, cells, cell_locks , particle_lock );
      }

      template<bool UseWeights, class CellsAccessorT, class GridCellLocksT, class ParticleLockT>
      inline void operator ()
        (
        size_t n,
        ComputePairBuffer2<UseWeights,true>& tab,
        double& _fx,
        double& _fy,
        double& _fz,
        CellsAccessorT cells,
        GridCellLocksT& cell_locks,
        ParticleLockT& particle_lock
        ) const
      {
        double _ep = 0.0;
        FakeMat3d fake_virial;
        this-> operator () ( n,tab,_ep,_fx,_fy,_fz, fake_virial, cells, cell_locks , particle_lock );
      }

      template<bool UseWeights, class CellsAccessorT, class Mat3dT>
      inline void operator ()
        (
        size_t n,
        ComputePairBuffer2<UseWeights,true>& tab,
        double& _ep,
        double& _fx,
        double& _fy,
        double& _fz,
        Mat3dT& virial,
        CellsAccessorT cells
        ) const
      {
        NullGridParticleLocks fake_locks;
        this-> operator () ( n,tab,_ep,_fx,_fy,_fz, virial, cells, fake_locks , fake_locks[0][0] );
      }
      template<bool UseWeights, class CellsAccessorT, class Mat3dT>
      inline void operator ()
        (
        size_t n,
        ComputePairBuffer2<UseWeights,true>& tab,
        double& _fx,
        double& _fy,
        double& _fz,
        Mat3dT& virial,
        CellsAccessorT cells
        ) const
      {
        double _ep = 0.0;
        NullGridParticleLocks fake_locks;
        this-> operator () ( n,tab,_ep,_fx,_fy,_fz, virial, cells, fake_locks , fake_locks[0][0] );
      }
      
      template<bool UseWeights, class CellsAccessorT, class Mat3dT, class GridCellLocksT, class ParticleLockT >
      inline void operator ()
        (
        size_t n,
        ComputePairBuffer2<UseWeights,true>& tab,
        double& _fx,
        double& _fy,
        double& _fz,
        Mat3dT& virial,
        CellsAccessorT cells,
        GridCellLocksT locks,
        ParticleLockT& lock_a
        ) const
      {
        double _ep = 0.0;
        this-> operator () ( n,tab,_ep,_fx,_fy,_fz, virial, cells, locks , lock_a );
      }
      
      // version with virial computation
      template<bool UseWeights, class CellsAccessorT, class Mat3dT, class GridCellLocksT, class ParticleLockT >
      inline void operator ()
        (
        size_t n,
        ComputePairBuffer2<UseWeights,true>& tab,
        double& _ep,
        double& _fx,
        double& _fy,
        double& _fz,
        Mat3dT& virial,
        CellsAccessorT cells,
        GridCellLocksT locks,
        ParticleLockT& lock_a
        ) const
      {
        static constexpr bool compute_virial = std::is_same_v< Mat3dT , Mat3d >;

        double* __restrict__ ep_b = nullptr; 
        double* __restrict__ fx_b = nullptr; 
        double* __restrict__ fy_b = nullptr; 
        double* __restrict__ fz_b = nullptr;
        [[maybe_unused]] Mat3dT* __restrict__ vir_b = nullptr;
        
        size_t current_cell_b = std::numeric_limits<size_t>::max();

        double tmp_ep = 0.;
        double tmp_fx = 0.;
        double tmp_fy = 0.;
        double tmp_fz = 0.;
       
        [[maybe_unused]] Mat3dT tmp_vir; // initializes to all 0

        for(size_t i=0;i<n;i++)
        {
          size_t cell_b=0, p_b=0;
          tab.nbh.get(i, cell_b, p_b);
          if( cell_b != current_cell_b )
          {
            current_cell_b = cell_b;
            ep_b = cells[cell_b][field::ep];
            fx_b = cells[cell_b][field::fx];
            fy_b = cells[cell_b][field::fy];
            fz_b = cells[cell_b][field::fz];
            if constexpr (compute_virial) { vir_b = cells[cell_b][field::virial]; }
          }

          const double r = std::sqrt(tab.d2[i]); 

          // pair interaction weighting
          const auto weight = tab.nbh_data.get(i); //compute_pair_buffer_get_weight( tab , i );

          double e=0.0, de=0.0;

#         if USTAMP_POTENTIAL_HANDLE_FORCE_WEIGHTING
          USTAMP_POTENTIAL_COMPUTE( p, pair_params, r, e, de , weight );
#         else
          USTAMP_POTENTIAL_COMPUTE( p, pair_params, r, e, de );
          e *= weight;
          de *= weight;
#         endif
          e -= ecut * weight;
          de /= r;

          const double drx = tab.drx[i];
          const double dry = tab.dry[i];
          const double drz = tab.drz[i];
          const double fe_x = de * drx;
          const double fe_y = de * dry;
          const double fe_z = de * drz;
          const double energy = 0.5 * e;

          tmp_ep += energy;
          tmp_fx += fe_x ; 
          tmp_fy += fe_y ; 
          tmp_fz += fe_z ; 
          
          [[maybe_unused]] Mat3dT pair_vir;
          if constexpr (compute_virial)
          {
            pair_vir = tensor( Vec3d{fe_x,fe_y,fe_z}, Vec3d{drx,dry,drz} ) * -0.5;
            tmp_vir += pair_vir;
          }

          auto& lock_b = locks[cell_b][p_b];
          lock_b.lock();
          ep_b[p_b] += energy;
          fx_b[p_b] -= fe_x ; 
          fy_b[p_b] -= fe_y ; 
          fz_b[p_b] -= fe_z ; 
          if constexpr (compute_virial) { vir_b[p_b] += pair_vir; }
          lock_b.unlock();
        }

        lock_a.lock();
        _ep += tmp_ep;
        _fx += tmp_fx;
        _fy += tmp_fy;
        _fz += tmp_fz;
        if constexpr (compute_virial) { virial += tmp_vir; }
        lock_a.unlock();
      }
    };

    // ============ operator slots ====================    
    ADD_SLOT( USTAMP_POTENTIAL_PARAMS , parameters         , INPUT , REQUIRED );
    ADD_SLOT( double                  , rcut               , INPUT , REQUIRED );
    ADD_SLOT( double                  , rcut_max           , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors, chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( ParticleSpecies         , species            , INPUT , OPTIONAL );
    ADD_SLOT( std::string             , type               , INPUT , OPTIONAL );
    ADD_SLOT( GridT                   , grid               , INPUT_OUTPUT );
    ADD_SLOT( Domain                  , domain             , INPUT , REQUIRED );
    ADD_SLOT( CompactGridPairWeights  , compact_nbh_weight , INPUT , OPTIONAL );
    ADD_SLOT( bool                    , enable_pair_weights, INPUT , true );

    ADD_SLOT( GridParticleLocks      , particle_locks      , INPUT , OPTIONAL , DocString{"particle spin locks"} );

    // ===============================================    

  public:

    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      *rcut_max = std::max( *rcut , *rcut_max );
      size_t n_cells = grid->number_of_cells();
      if( n_cells == 0 ) { return ; }

      const bool has_weight = compact_nbh_weight.has_value() && (*enable_pair_weights);
      
      // find the specy definition to use for this potential
      size_t specy_index = 0; 
      if( ! species.has_value() )
      {
        lerr_stream() << "no species input, cannot continue" << std::endl;
        std::abort();
      }
      if( species->size() != 1 && ! type.has_value() )
      {
        lerr_stream() <<"Error: exactly 1 atom type expected (single material) but found "<<species->size()<<std::endl;
        std::abort();
      }
      if( type.has_value() )
      {
        std::string specy_name = *type;
        for(size_t s=0;s<species->size();s++)
        {
          if( species->at(s).m_name == specy_name ) { specy_index = s; }
        }
        ldbg_stream() << "specy_name="<<specy_name<< ", specy_index = "<<specy_index<<std::endl;
      }
      PairPotentialParameters pair_params { species->at(specy_index) , species->at(specy_index) };
      
      // compute energy cutoff at rcut distance
      const double ecut = energy_cutoff( *parameters, pair_params , *rcut );
      
      SymetricForceOp force_op {*parameters,pair_params,ecut};
      auto cp_force_buf = make_compute_pair_buffer< ComputePairBuffer2<false,true> >();
      exanb::GridChunkNeighborsLightWeightIt<true> nbh_it{ *chunk_neighbors };

      assert( particle_locks.has_value() );
      
      // specialize for different execution scenarios
      if( domain->xform_is_identity() )
      {
        NullXForm cp_xform = {};
        if( has_weight )
        {
          CompactPairWeightIterator cp_weight = { compact_nbh_weight->m_cell_weights.data() };
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, ComputePairOptionalLocks<true> { particle_locks->data() } );
          compute_cell_particle_pairs( *grid, *rcut, true, optional, cp_force_buf, force_op , compute_fields , parallel_execution_context() );
        }
        else
        {
          ComputePairNullWeightIterator cp_weight = {};
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, ComputePairOptionalLocks<true> { particle_locks->data() } );
          compute_cell_particle_pairs( *grid, *rcut, true, optional, cp_force_buf, force_op , compute_fields , parallel_execution_context() );
        }
      }
      else
      {
        LinearXForm cp_xform = { domain->xform() };
        if( has_weight )
        {
          CompactPairWeightIterator cp_weight = { compact_nbh_weight->m_cell_weights.data() };
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, ComputePairOptionalLocks<true> { particle_locks->data() } );
          compute_cell_particle_pairs( *grid, *rcut, true, optional, cp_force_buf, force_op , compute_fields , parallel_execution_context() );
        }
        else
        {
          ComputePairNullWeightIterator cp_weight = {};
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, ComputePairOptionalLocks<true> { particle_locks->data() } );
          compute_cell_particle_pairs( *grid, *rcut, true, optional, cp_force_buf, force_op , compute_fields , parallel_execution_context() );
        }
      }
    }

  private:

    // used only once to compute energy cutoff
    static inline double energy_cutoff(const USTAMP_POTENTIAL_PARAMS& p, const PairPotentialParameters& pair_params, double rcut)
    {
      double tmp_e = 0.0;
      if( rcut > 0.0 )
      {
        double tmp_de = 0.0;
        USTAMP_POTENTIAL_COMPUTE(p,pair_params,rcut,tmp_e,tmp_de);
      }
      return tmp_e;
    }

  };


  namespace TemplateHelper
  {
    template<class GridT> using OPERATOR_NAME = ::exaStamp::OPERATOR_NAME < GridT >;
  }
  
  // === register factories ===  
  ONIKA_AUTORUN_INIT(pair_potential_singlemat_symetric)
  {
    OperatorNodeFactory::instance()->register_factory( OPERATOR_NAME_STR , make_grid_variant_operator< TemplateHelper::OPERATOR_NAME > );
  }

}

#undef OPERATOR_NAME
#undef OPERATOR_NAME_STR

