#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <onika/memory/allocator.h>
#include <exaStamp/npt/npt.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types_stream.h>
#include <exanb/core/physics_constants.h>

#include <exanb/core/string_utils.h>
#include <exanb/core/print_utils.h>

#include <sstream>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_rx, field::_ry, field::_rz>
    >
  class RemapNPT : public OperatorNode
  {
    ADD_SLOT( Domain                  , domain     , INPUT_OUTPUT );
    ADD_SLOT( GridT                   , grid       , INPUT_OUTPUT );
    ADD_SLOT( NPTContext              , npt_ctx , INPUT_OUTPUT );
    ADD_SLOT( std::string             , file       , INPUT , "npt.dat" );
    ADD_SLOT( long                    , timestep            , INPUT, REQUIRED);
    ADD_SLOT( double                  , physical_time       , INPUT );
    ADD_SLOT( Mat3d                  , xform_npt       , OUTPUT );

    using PointerTuple = onika::soatl::FieldPointerTuple< GridT::CellParticles::Alignment , GridT::CellParticles::ChunkSize, field::_rx, field::_ry, field::_rz >;

  public:
    inline void execute () override final
    {
      Mat3d hprec = domain->xform() * diag_matrix( domain->extent() - domain->origin() );
      Mat3d h = hprec;
      
      double dto2 = npt_ctx->dto/2.0;
      double dto4 = npt_ctx->dto/4.0;
      double dto8 = npt_ctx->dto/8.0;

      if (npt_ctx->pstyle == "TRICLINIC") {
        if (npt_ctx->p_flag[4]) {
          const double expfac = exp(dto8*npt_ctx->omega_dot[0]);
          h.m13 *= expfac;
          h.m13 += dto4*(npt_ctx->omega_dot[5]*h.m23+npt_ctx->omega_dot[4]*h.m33);
          h.m13 *= expfac;
        }
        
        if (npt_ctx->p_flag[3]) {
          const double expfac = exp(dto4*npt_ctx->omega_dot[1]);
          h.m23 *= expfac;
          h.m23 += dto2*(npt_ctx->omega_dot[3]*h.m33);
          h.m23 *= expfac;
        }

        if (npt_ctx->p_flag[5]) {
          const double expfac = exp(dto4*npt_ctx->omega_dot[0]);
          h.m12 *= expfac;
          h.m12 += dto2*(npt_ctx->omega_dot[5]*h.m22);
          h.m12 *= expfac;
        }

        if (npt_ctx->p_flag[4]) {
          const double expfac = exp(dto8*npt_ctx->omega_dot[0]);
          h.m13 *= expfac;
          h.m13 += dto4*(npt_ctx->omega_dot[5]*h.m23+npt_ctx->omega_dot[4]*h.m33);
          h.m13 *= expfac;
        }
      }
      
      if (npt_ctx->p_flag[0]) {
      	const double expfac = exp(npt_ctx->dto*npt_ctx->omega_dot[0]);
        h.m11 *= expfac;
      }

      if (npt_ctx->p_flag[1]) {
      	const double expfac = exp(npt_ctx->dto*npt_ctx->omega_dot[1]);
        h.m22 *= expfac;
        if (npt_ctx->scalexy) h.m12 *= expfac;
      }

      if (npt_ctx->p_flag[2]) {
      	const double expfac = exp(npt_ctx->dto*npt_ctx->omega_dot[2]);
        h.m33 *= expfac;
        if (npt_ctx->scalexz) h.m13 *= expfac;
        if (npt_ctx->scaleyz) h.m23 *= expfac;        
      }

      if (npt_ctx->pstyle == "TRICLINIC") {

        if (npt_ctx->p_flag[4]) {
          const double expfac = exp(dto8*npt_ctx->omega_dot[0]);
          h.m13 *= expfac;
          h.m13 += dto4*(npt_ctx->omega_dot[5]*h.m23+npt_ctx->omega_dot[4]*h.m33);
          h.m13 *= expfac;
        }
      
        if (npt_ctx->p_flag[3]) {
          const double expfac = exp(dto4*npt_ctx->omega_dot[1]);
          h.m23 *= expfac;
          h.m23 += dto2*(npt_ctx->omega_dot[3]*h.m33);
          h.m23 *= expfac;
        }
      
        if (npt_ctx->p_flag[5]) {
          const double expfac = exp(dto4*npt_ctx->omega_dot[0]);
          h.m12 *= expfac;
          h.m12 += dto2*(npt_ctx->omega_dot[5]*h.m22);
          h.m12 *= expfac;
        }
      
        if (npt_ctx->p_flag[4]) {
          const double expfac = exp(dto8*npt_ctx->omega_dot[0]);
          h.m13 *= expfac;
          h.m13 += dto4*(npt_ctx->omega_dot[5]*h.m23+npt_ctx->omega_dot[4]*h.m33);
          h.m13 *= expfac;
      	}
      
      }
      
      *xform_npt = h * inverse( hprec );
      domain->set_xform( *xform_npt  * domain->xform() );
      
    }

  };


 template<class GridT> using RemapNPTTmpl = RemapNPT<GridT>;

 // === register factories ===  
  ONIKA_AUTORUN_INIT(remap_npt)
  {
   OperatorNodeFactory::instance()->register_factory( "remap_npt", make_grid_variant_operator< RemapNPTTmpl > );
  }

}

