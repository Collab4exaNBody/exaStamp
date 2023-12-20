#pragma xstamp_cuda_enable

#pragma xstamp_grid_variant

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exanb/core/domain.h>

#include <onika/cuda/cuda.h>

#include <exanb/compute/compute_cell_particles.h>
#include <exanb/compute/compute_pair_optional_args.h>

#include <cmath>

namespace exaStamp
{
  using namespace exanb;
  using namespace onika;

  template<class XFormT>
  struct WallV2ComputeFunc
  {
    const Vec3d N;
    const double D;
    const double R;
    const long exposant = 12;
    const double epsilon = make_quantity(1.0e-19,"J").convert();
    XFormT xform;
    
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double rx, double ry, double rz, double& fx, double& fy, double& fz, double& ep ) const
    {
      Vec3d r { rx, ry, rz };
      r = xform.transformCoord(r);

      double d_sign = dot( r , N ) + D;
      double d = abs( d_sign );

      if( d <= R )
      {
        const double ratio = 1.0 - R / d;
        const double ratio_puis_exposant   = pow(ratio,exposant);

        double f_contrib =  -epsilon * exposant * ( R / (d*d) ) * ratio_puis_exposant / ratio ;
        Vec3d F = N * f_contrib;
        if (d_sign<0.) F = -F;

        // Energy
        ep += epsilon * ratio_puis_exposant;

        // Forces
        fx += F.x;
        fy += F.y;
        fz += F.z;
      }
    }

  };

}

namespace exanb
{
  template<class XFormT> struct ComputeCellParticlesTraits< exaStamp::WallV2ComputeFunc<XFormT> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{
  using namespace exanb;

  template<typename GridT
    , class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz, field::_ep >
    >
  class WallV2 : public OperatorNode
  {  
    ADD_SLOT( GridT  , grid    , INPUT_OUTPUT );
    ADD_SLOT( Vec3d  , normal  , INPUT , Vec3d{1.0,0.0,0.0} );
    ADD_SLOT( double , offset  , INPUT , 0.0 );
    ADD_SLOT( double , cutoff  , INPUT , REQUIRED );
    ADD_SLOT( Domain , domain  , INPUT , REQUIRED );
    ADD_SLOT( long   , exponent, INPUT , 12 );
    ADD_SLOT( double , epsilon , INPUT , make_quantity(1.0e-19,"J").convert() );

    static constexpr FieldSet<field::_rx,field::_ry,field::_rz, field::_fx, field::_fy, field::_fz, field::_ep > compute_field_set{};

  public:

    inline void execute () override final
    {
      if( ! domain->xform_is_identity() )
      {
        WallV2ComputeFunc<LinearXForm> func { *normal, - (*offset), *cutoff , *exponent, *epsilon, LinearXForm{ domain->xform() } };
        compute_cell_particles( *grid , false , func , compute_field_set , parallel_execution_context() );
      }
      else
      {
        WallV2ComputeFunc<NullXForm> func { *normal, - (*offset), *cutoff , *exponent, *epsilon, NullXForm{} };
        compute_cell_particles( *grid , false , func , compute_field_set , parallel_execution_context() );
      }
    }

  };
  
  template<class GridT> using WallV2Tmpl = WallV2<GridT>;
  
 // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "wall_v2", make_grid_variant_operator< WallV2Tmpl > );
  }

}

