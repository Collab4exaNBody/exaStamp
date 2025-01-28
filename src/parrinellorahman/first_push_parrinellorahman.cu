

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <onika/memory/allocator.h>
#include <exaStamp/parrinellorahman/parrinellorahman.h>
#include <exaStamp/parrinellorahman/parrinellorahman_yaml.h>
#include <exaStamp/parrinellorahman/parrinellorahman_stream.h>
#include <exanb/core/domain.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/unityConverterHelper.h>
#include <exanb/core/quantity.h>
#include <exanb/core/physics_constants.h>
#include <exanb/core/unityConverterHelper.h>

#include <onika/soatl/field_pointer_tuple.h>
#include <memory>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <onika/math/basic_types_stream.h>


namespace exaStamp
{

  // get particle virial tensor. assume the virial is null if particle hasn't virial field
  template<bool has_virial>
  static inline Mat3d get_particle_virial(const Mat3d* virials, size_t p_i, std::integral_constant<bool,has_virial> )
  {
    if constexpr (has_virial) { return virials[p_i]; }
    return Mat3d();
  }
  
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_vx, field::_vy, field::_vz, field::_ax, field::_ay, field::_az, field::_id >
    >
  struct FirstPushParrinelloRahman : public OperatorNode
  {
    using PointerTuple = onika::soatl::FieldPointerTuple<
                          GridT::CellParticles::Alignment , GridT::CellParticles::ChunkSize , 
                          field::_rx, field::_ry, field::_rz,
                          field::_vx, field::_vy, field::_vz,
                          field::_ax, field::_ay, field::_az >;

    // compile time constant indicating if grid has type virial
    using has_virial_field_t = typename GridT::CellParticles::template HasField < field::_virial > ;
    static constexpr std::integral_constant<bool,has_virial_field_t::value> has_virial_field{};
    
    ADD_SLOT( MPI_Comm                , mpi        , INPUT , REQUIRED );    
    ADD_SLOT( GridT                   , grid       , INPUT_OUTPUT );
    ADD_SLOT( double                  , dt         , INPUT , REQUIRED );
    ADD_SLOT( double                  , dt_scale   , INPUT , REQUIRED );
    ADD_SLOT( ParrinelloRahmanContext , parrinello_rahman_ctx , INPUT_OUTPUT );
    ADD_SLOT( Domain                  , domain     , INPUT );    
    ADD_SLOT( ParticleSpecies         , species    , INPUT , REQUIRED );
    
    inline void execute () override final
    {
      static const double boltzmann_internal = UnityConverterHelper::convert(legacy_constant::boltzmann, "m^2*kg/s^2/K");
      static const double conv_gnvtv = UnityConverterHelper::convert(1.0,"m^2/s^2");
      static const double conv_time = UnityConverterHelper::convert(1.0,"1/s^2");
      ParrinelloRahmanContext& data = *parrinello_rahman_ctx;

      //ldbg << "ParrinelloRahman::firstPush begin" << std::endl;

      // compute extension
      Mat3d mat = domain->xform();
      //Vec3d mat_ext = { norm(column1(mat)) , norm(column2(mat)) , norm(column3(mat)) };
      //Vec3d mat_exti = reciprocal( mat_ext );
      Vec3d dom_ext = domain->bounds_size();
      Vec3d dom_exti = reciprocal( dom_ext );
      
      //Vec3d ext = mat_ext * dom_ext;
      //const Vec3d exti = reciprocal( ext );
      const Mat3d mcr = multiply( mat , diag_matrix(dom_ext) );
      const Mat3d mcri = multiply( diag_matrix(dom_exti), inverse(mat) );

      data.h = mcr; // multiply( mat ,  diag_matrix( ext ) ); // == mcr
      data.updateMembers();
      
      data.m_gammaNVTp = 0.;
      data.hp = make_zero_matrix(); 
      data.hpp = make_zero_matrix();
      
      const double raw_dt = *dt;
      const double scale = *dt_scale;

      auto cells = grid->cells();
      IJK dims = grid->dimension();
      ssize_t gl = grid->ghost_layers();
      // size_t nb_particles = grid->number_of_particles() - grid->number_of_ghost_particles() ;

      const double dt = raw_dt * scale;
      const double dt2 = dt*dt;

      //data.print( ldbg );

      size_t total_particles = 0;          
      Mat3d hpp_cell = {0.};
      double gammaNVTp_cell = 0.;

/*
      ldbg << "conv_gnvtv = " << conv_gnvtv << std::endl;
      ldbg << "conv_time = " << conv_time << std::endl;
      ldbg << "dt   = " << dt << std::endl;
      ldbg << "dt2  = " << dt2 << std::endl;
      ldbg << "mat      = " << mat << std::endl;
      ldbg << "mat_ext  = " << mat_ext << std::endl;
      ldbg << "mat_ext*dom_ext  = " << dom_ext*mat_ext << std::endl;
      ldbg << "dom_ext  = " << dom_ext << std::endl;
      ldbg << "ext  = " << ext << std::endl;
      ldbg << "exti = " << exti << std::endl;
      ldbg << "mcr  = " << mcr << std::endl ;
      ldbg << "mcri = " << mcri << std::endl ;
*/

#     pragma omp parallel
      {
        PointerTuple ptrs;   
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) reduction(+:hpp_cell,gammaNVTp_cell,total_particles) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          int n = cells[i].size();
          cells[i].capture_pointers( ptrs );
          
          auto* __restrict__ rx = ptrs[ field::rx ]; ONIKA_ASSUME_ALIGNED(rx);
          auto* __restrict__ ry = ptrs[ field::ry ]; ONIKA_ASSUME_ALIGNED(ry);
          auto* __restrict__ rz = ptrs[ field::rz ]; ONIKA_ASSUME_ALIGNED(rz);

          auto* __restrict__ vx = ptrs[ field::vx ]; ONIKA_ASSUME_ALIGNED(vx);
          auto* __restrict__ vy = ptrs[ field::vy ]; ONIKA_ASSUME_ALIGNED(vy);
          auto* __restrict__ vz = ptrs[ field::vz ]; ONIKA_ASSUME_ALIGNED(vz);
          
          const auto* __restrict__ ax = ptrs[ field::ax ]; ONIKA_ASSUME_ALIGNED(ax);
          const auto* __restrict__ ay = ptrs[ field::ay ]; ONIKA_ASSUME_ALIGNED(ay);
          const auto* __restrict__ az = ptrs[ field::az ]; ONIKA_ASSUME_ALIGNED(az);
          
          const uint8_t* __restrict__ types = cells[i].field_pointer_or_null(field::type); ONIKA_ASSUME_ALIGNED(types);
          const Mat3d* __restrict__ vir = cells[i].field_pointer_or_null(field::virial); ONIKA_ASSUME_ALIGNED(vir);

          const auto* __restrict__ ids = cells[i].field_pointer_or_null(field::id); ONIKA_ASSUME_ALIGNED(ids);

          Mat3d virk;
          double mk;          
          Mat3d hpp_k = {0.};
          double gammaNVTp_k = 0.;

#         pragma omp simd reduction(+:hpp_k,gammaNVTp_k)
          for(int k=0;k<n;k++)
          {
            mk = species->at(types[k]).m_mass; // FIXME: warning, cases with no type field will crash, replace with get_type method similar to get_virial
            virk = get_particle_virial( vir, k, has_virial_field);

            Vec3d r = dom_exti * Vec3d{ rx[k], ry[k], rz[k] };
            Vec3d v = mcri * Vec3d{ vx[k], vy[k], vz[k] };
            Vec3d f = Vec3d{ ax[k], ay[k], az[k] };

            const auto cnv_1 = ( data.Giht * f ) / conv_time;
            const auto cnv_2 = data.GiGp * v;
            const auto cnv_3 = ( data.m_gammaNVT * v ) / conv_gnvtv;
            const auto a = ( cnv_1 - cnv_2 - cnv_3 ) * conv_time;
            
            // const auto tmp1 = ext * (r + dt * v + 0.5 * dt2 * a);
            const auto tmp1_v2 = dom_ext * (r + dt * v + 0.5 * dt2 * a);
            const auto tmp2 = mcr * (v + 0.5 * dt * a);

            rx[k] = tmp1_v2.x;
            ry[k] = tmp1_v2.y;
            rz[k] = tmp1_v2.z;

            vx[k] = tmp2.x;
            vy[k] = tmp2.y;
            vz[k] = tmp2.z;

            hpp_k += (data.h * tensor(v, v)) * mk + transpose(data.hi * virk);
            gammaNVTp_k += dot(v, data.G * v) * mk;
          }

          hpp_cell += hpp_k;
          gammaNVTp_cell += gammaNVTp_k;
          total_particles += n;
        }
        GRID_OMP_FOR_END
      }

      // sum hpp and gammaNVTp over all processors
      
      //ldbg << "syncDataReduce" << std::endl;
      {
        double tmp[11] = { hpp_cell.m11, hpp_cell.m12, hpp_cell.m13, hpp_cell.m21, hpp_cell.m22, hpp_cell.m23, hpp_cell.m31, hpp_cell.m32, hpp_cell.m33,
                           gammaNVTp_cell, static_cast<double>(total_particles) };
        MPI_Allreduce(MPI_IN_PLACE,tmp,11,MPI_DOUBLE,MPI_SUM,*mpi);
        hpp_cell.m11 = tmp[0];
        hpp_cell.m12 = tmp[1];
        hpp_cell.m13 = tmp[2];
        hpp_cell.m21 = tmp[3];
        hpp_cell.m22 = tmp[4];
        hpp_cell.m23 = tmp[5];
        hpp_cell.m31 = tmp[6];
        hpp_cell.m32 = tmp[7];
        hpp_cell.m33 = tmp[8];
        gammaNVTp_cell = tmp[9];
        total_particles = static_cast<size_t>( tmp[10] );
      }

      //ldbg << "total particles = "<< total_particles<<std::endl;
      data.m_gammaNVTp = gammaNVTp_cell;
      data.hpp = hpp_cell;
      //data.print( ldbg );

      // WARNING: try to move it afterward for test
      // data.apply_mask();

      //std::cout << "after hp/hpp update" << std::endl;
      data.hpp = (data.hpp - data.hp * (data.m_gammaNVT * data.m_config.m_masseB) - comatrix(data.h) * data.m_config.m_Pext ) / data.m_config.m_masseB;  
      data.apply_mask(); // WARNING moved here for test // seems to work for NPT_iso_xy test case
      
      data.h = data.h + data.hp * dt + data.hpp * 0.5 * dt2;
      data.hp = data.hp + data.hpp * (0.5 * dt);
      data.m_gammaNVT += 0.5 * dt * data.m_gammaNVTp;

      //data.print( ldbg );

      const double oldtrace = trace_matrix(data.hpthp);
      //ldbg << "oldtrace = "<<oldtrace << std::endl;
      
      //ldbg << "updateMembers" << std::endl;
      // update scheme parameters
      data.updateMembers();
      //data.print( ldbg );

      // update gamma
				   
//      ldbg << "boltzmann =" << boltzmann_internal << std::endl;
//      ldbg << "Text =" << data.m_config.m_Text << std::endl;

      double tmp = data.m_config.m_masseB * oldtrace - 3. * total_particles * boltzmann_internal * data.m_config.m_Text;
//      ldbg << "tmp = " << tmp << std::endl;
      data.m_gammaNVTp = (data.m_gammaNVTp + tmp) / data.m_config.m_masseNVT;
      data.m_gammaNVT = data.m_gammaNVT + 0.5 * dt * data.m_gammaNVTp;

      // ldbg << "after gamma update" << std::endl;
      //data.print( ldbg );
      //ldbg << "--- end first push ---" << std::endl;
    }

  };

  template<class GridT> using FirstPushParrinelloRahmanTmpl = FirstPushParrinelloRahman<GridT>;
  
 // === register factories ===  
  ONIKA_AUTORUN_INIT(first_push_parrinellorahman)
  {
   OperatorNodeFactory::instance()->register_factory( "first_push_parrinellorahman", make_grid_variant_operator< FirstPushParrinelloRahmanTmpl > );
  }

}

