#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <onika/memory/allocator.h>
#include <exaStamp/npt/npt.h>
#include <exanb/core/domain.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/physics_constants.h>

#include <exanb/core/string_utils.h>
#include <exanb/core/print_utils.h>

#include <sstream>

enum{NONE,XYZ,XY,YZ,XZ};
enum{ISO,ANISO,TRICLINIC};

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
      ldbg << std::fixed;
      ldbg << std::setprecision(10);
      
      //      const Mat3d hprec_inv = domain->inv_xform();     
      Mat3d hprec = domain->xform() * diag_matrix( domain->extent() - domain->origin() );
      Mat3d h = hprec;
      Mat3d hprec_inv = inverse (hprec);
      ldbg << "H beginning remap" << std::endl;
      ldbg << "h[0]=" << h.m11 << std::endl;
      ldbg << "h[1]=" << h.m22 << std::endl;
      ldbg << "h[2]=" << h.m33 << std::endl;      
      ldbg << "h[3]=" << h.m23 << std::endl;      
      ldbg << "h[4]=" << h.m13 << std::endl;      
      ldbg << "h[5]=" << h.m12 << std::endl;

      ldbg << "Omega exp vals" << std::endl;      
      for (int i = 0; i < 6; i++)
        ldbg << "omega["<<i<<"]=" << npt_ctx->omega_dot[i] << std::endl;
      
      int i;
      double oldlo,oldhi;
      double expfac;
      
      double dto2 = npt_ctx->dto/2.0;
      double dto4 = npt_ctx->dto/4.0;
      double dto8 = npt_ctx->dto/8.0;

      // off-diagonal components, first half
      // h[0] == domain_xform.m11
      // h[1] == domain_xform.m22
      // h[2] == domain_xform.m33      
      // h[3] == domain_xform.m23
      // h[4] == domain_xform.m13
      // h[5] == domain_xform.m12
      
      if (npt_ctx->p_flag[4]) {
        expfac = exp(dto8*npt_ctx->omega_dot[0]);
        ldbg << "expfac = " << expfac << std::endl;
        h.m13 *= expfac;
        h.m13 += dto4*(npt_ctx->omega_dot[5]*h.m23+npt_ctx->omega_dot[4]*h.m33);
        h.m13 *= expfac;
      }

      if (npt_ctx->p_flag[3]) {
        expfac = exp(dto4*npt_ctx->omega_dot[1]);
        ldbg << "expfac = " << expfac << std::endl;
        h.m23 *= expfac;
        h.m23 += dto2*(npt_ctx->omega_dot[3]*h.m33);
        h.m23 *= expfac;
      }

      if (npt_ctx->p_flag[5]) {
        expfac = exp(dto4*npt_ctx->omega_dot[0]);
        ldbg << "expfac = " << expfac << std::endl;        
        h.m12 *= expfac;
        h.m12 += dto2*(npt_ctx->omega_dot[5]*h.m22);
        h.m12 *= expfac;
      }

      if (npt_ctx->p_flag[4]) {
        expfac = exp(dto8*npt_ctx->omega_dot[0]);
        ldbg << "expfac = " << expfac << std::endl;
        h.m13 *= expfac;
        h.m13 += dto4*(npt_ctx->omega_dot[5]*h.m23+npt_ctx->omega_dot[4]*h.m33);
        h.m13 *= expfac;
      }

      // scale diagonal components
      // scale tilt factors with cell, if set

      if (npt_ctx->p_flag[0]) {
        //      	oldlo = domain->boxlo[0];
        //      	oldhi = domain->boxhi[0];
      	expfac = exp(npt_ctx->dto*npt_ctx->omega_dot[0]);
        ldbg << "expfac = " << expfac << std::endl;
        ldbg << "hm11 before = " << h.m11 << std::endl;        
        h.m11 *= expfac;
        ldbg << "hm11 after  = " << h.m11 << std::endl;                
        //      	domain->boxlo[0] = (oldlo-fixedpoint[0])*expfac + fixedpoint[0];
        //      	domain->boxhi[0] = (oldhi-fixedpoint[0])*expfac + fixedpoint[0];
      }

      if (npt_ctx->p_flag[1]) {
	//      	oldlo = domain->boxlo[1];
	//      	oldhi = domain->boxhi[1];
      	expfac = exp(npt_ctx->dto*npt_ctx->omega_dot[1]);
        ldbg << "expfac = " << expfac << std::endl;        
        h.m12 *= expfac;        
        h.m22 *= expfac;
	//      	domain->boxlo[1] = (oldlo-fixedpoint[1])*expfac + fixedpoint[1];
	//      	domain->boxhi[1] = (oldhi-fixedpoint[1])*expfac + fixedpoint[1];
        //      	if (npt_ctx->scalexy) h.m12 *= expfac;
      }

      if (npt_ctx->p_flag[2]) {
	//      	oldlo = domain->boxlo[2];
	//      	oldhi = domain->boxhi[2];
      	expfac = exp(npt_ctx->dto*npt_ctx->omega_dot[2]);
        ldbg << "expfac = " << expfac << std::endl;        
        h.m13 *= expfac;
        h.m23 *= expfac;
        h.m33 *= expfac;
	//      	domain->boxlo[2] = (oldlo-fixedpoint[2])*expfac + fixedpoint[2];
	//      	domain->boxhi[2] = (oldhi-fixedpoint[2])*expfac + fixedpoint[2];
      	// if (npt_ctx->scalexz) h.m13 *= expfac;
      	// if (npt_ctx->scaleyz) h.m23 *= expfac;
      }

      // off-diagonal components, second half

      //      if (npt_ctx->pstyle == TRICLINIC) {

      if (npt_ctx->p_flag[4]) {
        expfac = exp(dto8*npt_ctx->omega_dot[0]);
        ldbg << "expfac = " << expfac << std::endl;        
        h.m13 *= expfac;
        h.m13 += dto4*(npt_ctx->omega_dot[5]*h.m23+npt_ctx->omega_dot[4]*h.m33);
        h.m13 *= expfac;
      }
      
      if (npt_ctx->p_flag[3]) {
        expfac = exp(dto4*npt_ctx->omega_dot[1]);
        ldbg << "expfac = " << expfac << std::endl;        
        h.m23 *= expfac;
        h.m23 += dto2*(npt_ctx->omega_dot[3]*h.m33);
        h.m23 *= expfac;
      }
      
      if (npt_ctx->p_flag[5]) {
        expfac = exp(dto4*npt_ctx->omega_dot[0]);
        ldbg << "expfac = " << expfac << std::endl;        
        h.m12 *= expfac;
        h.m12 += dto2*(npt_ctx->omega_dot[5]*h.m22);
        h.m12 *= expfac;
      }
      
      if (npt_ctx->p_flag[4]) {
        expfac = exp(dto8*npt_ctx->omega_dot[0]);
        ldbg << "expfac = " << expfac << std::endl;        
        h.m13 *= expfac;
        h.m13 += dto4*(npt_ctx->omega_dot[5]*h.m23+npt_ctx->omega_dot[4]*h.m33);
        h.m13 *= expfac;
      	}
      
      //      }
      ldbg << "H after calcul" << std::endl;
      
      ldbg << "h[0]=" << h.m11 << std::endl;
      ldbg << "h[1]=" << h.m22 << std::endl;
      ldbg << "h[2]=" << h.m33 << std::endl;      
      ldbg << "h[3]=" << h.m23 << std::endl;      
      ldbg << "h[4]=" << h.m13 << std::endl;      
      ldbg << "h[5]=" << h.m12 << std::endl;

      // h = xform x hprec
      // xform = h * inverse (hprec);
      
      *xform_npt = h * inverse( hprec );
      // ldbg << "*xform_npt = " << *xform_npt << std::endl;
      // ldbg << "xform domain = " << domain->xform() << std::endl;
      // ldbg << "domain origin  = " << domain->origin() << std::endl;
      // ldbg << "domain extent  = " << domain->extent() << std::endl;
      // ldbg << "diag (ext-orig) = " << diag_matrix(domain->extent()-domain->origin()) << std::endl;
      // ldbg << "xform * (ext-ori)  = " << domain->xform() * diag_matrix(domain->extent()-domain->origin()) << std::endl;      
      domain->set_xform( *xform_npt  * domain->xform() );
      // ldbg << "xform domain = " << domain->xform() << std::endl;
      // ldbg << "domain origin  = " << domain->origin() << std::endl;
      // ldbg << "domain extent  = " << domain->extent() << std::endl;
      // ldbg << "diag (ext-orig) = " << diag_matrix(domain->extent()-domain->origin()) << std::endl;      
      ldbg << "xform * (ext-ori)  = " << domain->xform() * diag_matrix(domain->extent()-domain->origin()) << std::endl;      
      ldbg << "FIN REMAP" << std::endl;
      Mat3d hfinal = domain->xform() * diag_matrix( domain->extent() - domain->origin() );
      ldbg << "h[0]=" << hfinal.m11 << std::endl;
      ldbg << "h[1]=" << hfinal.m22 << std::endl;
      ldbg << "h[2]=" << hfinal.m33 << std::endl;      
      ldbg << "h[3]=" << hfinal.m23 << std::endl;      
      ldbg << "h[4]=" << hfinal.m13 << std::endl;      
      ldbg << "h[5]=" << hfinal.m12 << std::endl;
      ldbg <<"AAAAAAAAAAAAAAAAAAa" << std::endl;
//       GridT& grid = *(this->grid);      
//       auto cells = grid.cells();
//       IJK dims = grid.dimension();
//       ssize_t gl = grid.ghost_layers();

// #     pragma omp parallel
//       {
//         PointerTuple ptrs;   
//         GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
//         {
//           size_t i = grid_ijk_to_index( dims , loc + gl );
//           const unsigned int n = cells[i].size();
//           cells[i].capture_pointers( ptrs );

//           auto* __restrict__ rx = ptrs[ field::rx ];
//           auto* __restrict__ ry = ptrs[ field::ry ];
//           auto* __restrict__ rz = ptrs[ field::rz ];

// 	  for(unsigned int j=0;j<n;j++)
// 	    {
// 	      Vec3d lambda = hprec_inv * Vec3d{rx[j],ry[j],rz[j]};
// 	      Vec3d r = h * lambda;
// 	      rx[j] = r.x;
// 	      ry[j] = r.y;
// 	      rz[j] = r.z;
// 	    }
//         }
//         GRID_OMP_FOR_END
//       }
      //      std::abort();
    }

  };


 template<class GridT> using RemapNPTTmpl = RemapNPT<GridT>;

 // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "remap_npt", make_grid_variant_operator< RemapNPTTmpl > );
  }

}

