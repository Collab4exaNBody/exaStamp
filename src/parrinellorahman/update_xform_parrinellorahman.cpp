#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <onika/memory/allocator.h>
#include <exaStamp/parrinellorahman/parrinellorahman.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types_stream.h>
#include <onika/physics/constants.h>

#include <onika/string_utils.h>
#include <onika/print_utils.h>

#include <sstream>


namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_vx, field::_vy, field::_vz>
    >
  class UpdateXFormParrinelloRahman : public OperatorNode
  {
    ADD_SLOT( Domain                  , domain     , INPUT_OUTPUT );
    ADD_SLOT( GridT                   , grid       , INPUT_OUTPUT );
    ADD_SLOT( ParrinelloRahmanContext , parrinello_rahman_ctx , INPUT_OUTPUT );
    ADD_SLOT( std::string             , file       , INPUT , "parrinello_rahman.dat" );
    ADD_SLOT( long                    , timestep            , INPUT, REQUIRED);
    ADD_SLOT( double                  , physical_time       , INPUT );
    ADD_SLOT( bool                    , force_append_thermo , INPUT , false );

    using PointerTuple = onika::soatl::FieldPointerTuple< GridT::CellParticles::Alignment , GridT::CellParticles::ChunkSize, field::_vx, field::_vy, field::_vz >;

  public:
    inline void execute () override final
    {
      static const std::string header = "     Step     Time (ps)     La(ang)     Lb(ang)     Lc(ang)     alpha(degree)     beta(degree)     gamma(degree)     H11(ang)     H12(ang)     H13(ang)     H21(ang)     H22(ang)     H23(ang)     H31(ang)     H32(ang)     H33(ang)";

      std::ostringstream oss;

      if( *timestep == 0 )
      {
        oss << header;
        oss << '\n';
      }

      //ldbg << "UpdateXFormParrinelloRahman" << std::endl;
      
      const Mat3d oldXForm = domain->xform();     
      const Mat3d oldMat = oldXForm * diag_matrix( domain->bounds_size() ) ;
      //const Vec3d oldExt = { norm(column1(oldMat)) , norm(column2(oldMat)) , norm(column3(oldMat)) };
      ldbg << "oldXForm = "<< oldXForm << std::endl;
      //ldbg << "oldMat = "<< oldMat << std::endl;
      //ldbg << "oldExt = "<< oldExt << std::endl;
//#     ifndef NDEBUG
      parrinello_rahman_ctx->print( ldbg );
//#     endif
      
      const Mat3d newMat = parrinello_rahman_ctx->h;
      const Mat3d newXForm = multiply( newMat , diag_matrix( reciprocal( domain->bounds_size() ) ) );
      
      //const Vec3d newExt = { norm(column1(newMat)) , norm(column2(newMat)) , norm(column3(newMat)) };
      ldbg << "newXForm = "<< newXForm << std::endl;
      ldbg << "newMat = "<< newMat << std::endl;
      //      ldbg << "newExt = "<< newExt << std::endl;

      Vec3d a { newMat.m11, newMat.m21, newMat.m31 };
      Vec3d b { newMat.m12, newMat.m22, newMat.m32 };
      Vec3d c { newMat.m13, newMat.m23, newMat.m33 };
      double La = norm(a);
      double Lb = norm(b);
      double Lc = norm(c);

      double alpha = 180. * asin(norm(cross(a,b)) / (La * Lb)) / M_PI;
      double beta  = 180. * asin(norm(cross(b,c)) / (Lb * Lc)) / M_PI;
      double gamma = 180. * asin(norm(cross(c,a)) / (Lc * La)) / M_PI;

      oss << onika::format_string("%9ld % .6e % .10e  % .10e  % .10e  % .10e  % .10e  % .10e  % .10e  % .10e  % .10e  % .10e  % .10e  % .10e  % .10e  % .10e  % .10e \n",
			   *timestep,
			   *physical_time,
			   La,
			   Lb,
			   Lc,
			   alpha,
			   beta,
			   gamma,
			   newMat.m11,
			   newMat.m12,
			   newMat.m13,
			   newMat.m21,
			   newMat.m22,
			   newMat.m23,
			   newMat.m31,
			   newMat.m32,
			   newMat.m33);

      FileAppendWriteBuffer::instance().append_to_file( *file , oss.str(), *force_append_thermo);

      //const Mat3d newPM = parrinello_rahman_ctx->h * diag_matrix( reciprocal(newExt) );
      //ldbg << "newPM = "<< newPM << std::endl;

      auto cells = grid->cells();
      IJK dims = grid->dimension();
      ssize_t gl = grid->ghost_layers();

      //const Vec3d scaling = newExt / oldExt;
      const Mat3d delta_xform = newMat * inverse( oldMat ) ; //* diag_matrix( scaling );

      //ldbg << "scaling = " << scaling << std::endl;
      ldbg << "delta_xform     = " << delta_xform << std::endl;
      //ldbg << "newPM * scaling = " << newPM * diag_matrix( scaling ) << std::endl;

      // does not transform ghosts as ghost particles' positions will be updated right after
#     pragma omp parallel
      {
        PointerTuple ptrs;   
        GRID_OMP_FOR_BEGIN(dims -2*gl , i , loc)
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          int n = cells[i].size();
          cells[i].capture_pointers( ptrs );

          auto* __restrict__ vx = ptrs[ field::vx ]; ONIKA_ASSUME_ALIGNED(vx);
          auto* __restrict__ vy = ptrs[ field::vy ]; ONIKA_ASSUME_ALIGNED(vy);
          auto* __restrict__ vz = ptrs[ field::vz ]; ONIKA_ASSUME_ALIGNED(vz);

#         pragma omp simd
          for(int k=0;k<n;k++)
          {
            Vec3d v = {vx[k],vy[k],vz[k]};
            v = delta_xform * v;
            vx[k] = v.x;
            vy[k] = v.y;
            vz[k] = v.z;
          }
        }
        GRID_OMP_FOR_END
      }

      domain->set_xform( newXForm );
    }

  };


 template<class GridT> using UpdateXFormParrinelloRahmanTmpl = UpdateXFormParrinelloRahman<GridT>;

 // === register factories ===  
  ONIKA_AUTORUN_INIT(update_xform_parrinellorahman)
  {
   OperatorNodeFactory::instance()->register_factory( "update_xform_parrinellorahman", make_grid_variant_operator< UpdateXFormParrinelloRahmanTmpl > );
  }

}

