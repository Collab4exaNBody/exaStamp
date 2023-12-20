#pragma xstamp_grid_variant

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/log.h>
#include <exanb/core/cpp_utils.h>

#include <exaStamp/potential/coul_wolf/coul_wolf.h>

#include <exanb/core/config.h> // for MAX_PARTICLE_NEIGHBORS constant
#include <exanb/particle_neighbors/chunk_neighbors.h>

// this allows for parallel compilation of templated operator for each available field set


namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz, field::_charge >
    >
  class CoulWolfPC : public OperatorNode
  {      
    // ========= I/O slots =======================
    ADD_SLOT( CoulWolfParms    , parameters        , INPUT , REQUIRED );
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

    
    // additional storage space added to compute buffer created by compute_pair_singlemat
    struct alignas(DEFAULT_ALIGNMENT) CoulWolfComputeBuffer
    {
      alignas(DEFAULT_ALIGNMENT) double charge[exanb::MAX_PARTICLE_NEIGHBORS];
    };
     
    // functor that populate compute buffer's extended storage for particle charges
    struct CopyParticleCharge
    {
      template<typename ComputeBufferT, typename FieldArraysT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2, const FieldArraysT * cells, size_t cell_b, size_t p_b, double weight) const noexcept
      {      
        assert( ssize_t(tab.count) < ssize_t(tab.MaxNeighbors) );
        tab.ext.charge[tab.count] = cells[cell_b][field::charge][p_b];
        exanb::DefaultComputePairBufferAppendFunc{} ( tab, dr, d2, cells, cell_b, p_b, weight );
      }
    };

    // shortcut to the Compute buffer used (and passed to functor) by compute_pair_singlemat
    using ComputeBuffer = ComputePairBuffer2<false,false,CoulWolfComputeBuffer,CopyParticleCharge>;

  public:
    // Operator execution
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      double rcut = parameters->rc;
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

      auto force_buf = make_compute_pair_buffer<ComputeBuffer>();
      ForceOp force_op { *parameters };
      
      if( domain->xform_is_identity() )
      {
        NullXForm cp_xform;
        auto optional = make_compute_pair_optional_args( nbh_it, ComputePairNullWeightIterator{} , cp_xform, cp_locks );
        compute_cell_particle_pairs( *grid, rcut, *ghost, optional, force_buf, force_op , ComputeFields{} , DefaultPositionFields{} , parallel_execution_context() );
      }
      else
      {
        LinearXForm cp_xform { domain->xform() };
        auto optional = make_compute_pair_optional_args( nbh_it, ComputePairNullWeightIterator{} , cp_xform, cp_locks );
        compute_cell_particle_pairs( *grid, rcut, *ghost, optional, force_buf, force_op , ComputeFields{} , DefaultPositionFields{} , parallel_execution_context() );
      }
      
    }

    // Wolf Summation Compute functor
    struct alignas(DEFAULT_ALIGNMENT) ForceOp 
    {
      // poetential parameters
      const CoulWolfParms m_params;

      inline void operator ()
        (        
        size_t n,
        const ComputeBuffer& tab,
        double& ep,
        double& fx,
        double& fy,
        double& fz,
        double charge,    // per particle (read only)
        CellParticles* unused
        ) const
      {
        FakeMat3d virial;
        this->operator() ( n,tab,ep, fx,fy,fz, charge, virial, unused );
      }

      template<class Mat3dT>
      inline void operator ()
        (        
        size_t n,
        const ComputeBuffer& tab,
        double& ep,
        double& fx,
        double& fy,
        double& fz,
        double charge,    // per particle (read only)
        Mat3dT& virial,
        CellParticles*
        ) const
      {

        // energy and force contributions to the particle
        double _ep = 0.;
        double _fx = 0.;
        double _fy = 0.;
        double _fz = 0.;

        Mat3dT _vir; // default constructor defines all elements to 0
        // assert( _vir.m11==0 && _vir.m12==0 && _vir.m13==0 && _vir.m21==0 && _vir.m22==0 && _vir.m23==0 && _vir.m31==0 && _vir.m32==0 && _vir.m33==0);

	double e_self = -(m_params.e_shift / 2.0 + m_params.alpha / sqrt(M_PI)) * charge * charge * m_params.qqrd2e;
	_ep = UnityConverterHelper::convert(e_self, "eV");
#       pragma omp simd reduction(+:_ep,_fx,_fy,_fz,_vir)
        for(size_t i=0;i<n;i++)
        {
          const double r = std::sqrt(tab.d2[i]);
          const double nbh_charge = tab.ext.charge[i];
          double e=0.0, de=0.0;
	  coul_wolf_compute_energy( m_params, charge, nbh_charge, r, e, de );	  
          de /= r;

          const double drx = tab.drx[i];
          const double dry = tab.dry[i];
          const double drz = tab.drz[i];
          const double fe_x = de * drx;
          const double fe_y = de * dry;
          const double fe_z = de * drz;
          _fx += fe_x;
          _fy += fe_y;
          _fz += fe_z;
          _ep += .5 * e;
          _vir += tensor( Vec3d{fe_x,fe_y,fe_z}, Vec3d{drx,dry,drz} ) * -0.5;
        }
        ep += _ep;
        fx += _fx; 
        fy += _fy; 
        fz += _fz; 
        virial += _vir;
      }
    };

  };

  template<class GridT> using CoulWolfPCTmpl = CoulWolfPC<GridT>;

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {  
    OperatorNodeFactory::instance()->register_factory( "coul_wolf_pc" , make_grid_variant_operator<CoulWolfPCTmpl> );
  }

}


