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
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <exaStamp/unit_system.h>
#include <exaStamp/compute/thermodynamic_state.h>
#include <onika/math/basic_types_stream.h>
#include <onika/value_streamer.h>

#include <onika/soatl/field_pointer_tuple.h>
#include <memory>
#include <mpi.h>

namespace exaStamp
{

  template<class T,size_t N> void my_array_tester( T array[N] )
  {
    std::cout << "array size = "<< N << std::endl;
  }

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_vx, field::_vy, field::_vz, field::_ax, field::_ay, field::_az, field::_virial>
    >
  class ConvergencePushParrinelloRahman : public OperatorNode
  {
    using PointerTuple = onika::soatl::FieldPointerTuple<
                          GridT::CellParticles::Alignment , GridT::CellParticles::ChunkSize , 
                          field::_vx, field::_vy, field::_vz,
                          field::_ax, field::_ay, field::_az,
                          field::_virial >;

    using PointerTuple2 = onika::soatl::FieldPointerTuple<
                          GridT::CellParticles::Alignment , GridT::CellParticles::ChunkSize , 
                          field::_vx, field::_vy, field::_vz,
                          field::_ax, field::_ay, field::_az >;

    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;

    ADD_SLOT( MPI_Comm           , mpi                 , INPUT , REQUIRED ); 
   
    ADD_SLOT( double             , dt                  , INPUT );
    ADD_SLOT( double             , dt_scale            , INPUT , 1.0 );

    ADD_SLOT( ParrinelloRahmanContext, parrinello_rahman_ctx , INPUT_OUTPUT );

    ADD_SLOT( Domain             , domain              , INPUT , REQUIRED );    

    ADD_SLOT( ParticleSpecies    , species             , INPUT , REQUIRED );
    ADD_SLOT( ThermodynamicState , thermodynamic_state , INPUT );

    ADD_SLOT( double             , epsilon             , INPUT , 1.0E-06 );
    ADD_SLOT( long               , max_iter            , INPUT , 100 );

    ADD_SLOT( GridT              , grid                , INPUT_OUTPUT );

  public:
    inline void execute () override final
    {
      using onika::ValueStreamer;

      static constexpr double conv_temperature = 1.e4 * onika::physics::atomicMass / onika::physics::boltzmann ;
      static constexpr double boltzmann_internal = EXASTAMP_CONST_QUANTITY( onika::physics::boltzmann * ( m^2 ) * kg / ( s^2 ) / K );
      static constexpr double conv_gammadt = EXASTAMP_CONST_QUANTITY( 1.0 / ( m^2 ) );
      static constexpr double conv_gammanvt = EXASTAMP_CONST_QUANTITY( 1.0 * s / ( m^2 ) );
      // ldbg << "conv_gammanvt = " << conv_gammanvt << std::endl;

      const double raw_dt = *(this->dt);
      const double scale = *dt_scale;
      const Mat3d inv_xform = domain->inv_xform();
      const ThermodynamicState& sim_info = *(this->thermodynamic_state);
      const double epsilon = *(this->epsilon);
      const size_t maxIter = *max_iter;
      ParrinelloRahmanContext& data = *parrinello_rahman_ctx;

      const Vec3d ext = domain->bounds_size();
      const Vec3d exti = reciprocal( ext );
      const Mat3d mcr = multiply( domain->xform() , diag_matrix(ext) );
      const Mat3d mcri = multiply( diag_matrix(exti), inv_xform );

      const auto gammaNVT_ref = data.m_gammaNVT;
      const auto hp_ref = data.hp;
      const auto hpp_ref = data.hpp;
      
      // double err = 0.;
      // bool hasConverged = false;
            
      const double dt = raw_dt * scale;
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      ssize_t gl = grid->ghost_layers();
      // size_t nb_particles = grid->number_of_particles() - grid->number_of_ghost_particles();
/*
      ldbg<<"ConvergencePushParrinelloRahman:"<<std::endl
          <<"dt   = " << raw_dt<<std::endl
          <<"dt_s = " << scale<<std::endl
          <<"ext  = " << ext << std::endl
          <<"exti = " << exti << std::endl
          <<"mcr  = " << mcr << std::endl
          <<"mcri = " << mcri << std::endl;
*/
      // backup original PR data to compute final velocity
      ParrinelloRahmanContext saved_prdata = data;

      // intialize temperature      
      double temperature_old = sim_info.temperature_scal() / sim_info.particle_count() * conv_temperature;

      // intial value over threshold
      size_t count = 0;
      double err = epsilon + 1.0;

      while ( err > epsilon && count < maxIter)
      {
        ++ count;
        const auto gammaNVT_old = data.m_gammaNVT;
        //const auto gammaNVTp_old = data.m_gammaNVTp;
        //const auto hp_old = data.hp;
/*
        ldbg << onika::default_stream_format
             << "iteration "<< count <<" / "<<maxIter<<std::endl
             << "temp_old="<<temperature_old<<" gammaNVT_old="<<gammaNVT_old<<" gammaNVTp_old="<<gammaNVTp_old<<std::endl
             << "hp_old="<<hp_old <<std::endl;
*/
        // secondPushPR begins
        // sum these two to drive convergence
        Mat3d hpp_cell; // = 0
        double gammaNVTp_cell = 0.;
        
        // sum these to compute temperature
        Vec3d momentum; // = 0
        Vec3d kinetic_energy; // = 0
        double total_mass = 0.;
        size_t total_particles = 0;

        // backup data used during last iteration to store final velocities
        saved_prdata = data;

        // equivalent to secondPushPR in exaStamp V1
#       pragma omp parallel
        {
          PointerTuple ptrs; 
          GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, reduction(+:hpp_cell,gammaNVTp_cell,momentum,kinetic_energy,total_mass,total_particles) schedule(dynamic))
          {
              size_t i = grid_ijk_to_index( dims , loc + gl );
              int n = cells[i].size();
              cells[i].capture_pointers( ptrs );

              const auto* __restrict__ vx = ptrs[ field::vx ]; ONIKA_ASSUME_ALIGNED(vx);
              const auto* __restrict__ vy = ptrs[ field::vy ]; ONIKA_ASSUME_ALIGNED(vy);
              const auto* __restrict__ vz = ptrs[ field::vz ]; ONIKA_ASSUME_ALIGNED(vz);
              
              const auto* __restrict__ ax = ptrs[ field::ax ]; ONIKA_ASSUME_ALIGNED(ax);
              const auto* __restrict__ ay = ptrs[ field::ay ]; ONIKA_ASSUME_ALIGNED(ay);
              const auto* __restrict__ az = ptrs[ field::az ]; ONIKA_ASSUME_ALIGNED(az);

              const auto* __restrict__ vir = ptrs[ field::virial ]; ONIKA_ASSUME_ALIGNED(vir);

              const auto* __restrict__ types = cells[i].field_pointer_or_null(field::type); ONIKA_ASSUME_ALIGNED(types);

              Mat3d hpp_k; // = 0
              double gammaNVTp_k = 0.;                        
              Vec3d local_momentum = {0.,0.,0.};
              Vec3d local_kinetic_energy = {0.,0.,0.};
              double local_mass = 0.;

#             pragma omp simd reduction(+:hpp_k,gammaNVTp_k,local_momentum,local_kinetic_energy,local_mass)
              for(int k=0;k<n;k++)
              {
                const Vec3d f = Vec3d{ax[k], ay[k], az[k]};           
                const Vec3d vSaved = mcri * Vec3d{vx[k], vy[k], vz[k]};
                const Vec3d p = (0.5 * dt) * (data.Giht * f) + vSaved;
                const Mat3d q = (0.5 * dt) * data.GiGp + (1. + 0.5 * dt * data.m_gammaNVT * conv_gammadt) * make_identity_matrix();
                const Vec3d v = inverse(q) * p;
                double mk = 0.0;
		if constexpr ( has_type_field ) mk = (*species)[ types[k] ].m_mass;
		else mk = (*species)[0].m_mass;
                const Mat3d hpp_part = mk * multiply(data.h, tensor(v, v)) + transpose(multiply(data.hi, vir[k]));
                hpp_k += hpp_part;
                const double gammaNVTp_part = mk * dot(v, data.G * v);
                gammaNVTp_k += gammaNVTp_part;

                // temperature computation part
                Vec3d temp_v = mcr * v;
                local_mass += mk;
                local_momentum += temp_v * mk;
                local_kinetic_energy += temp_v * temp_v * mk; // x0.5 later*
              }
              
              hpp_cell += hpp_k;
              gammaNVTp_cell += gammaNVTp_k;
              momentum += local_momentum;
              kinetic_energy += local_kinetic_energy;
              total_mass += local_mass;
              total_particles += n;          
          }
          GRID_OMP_FOR_END
        }
        // secondPushPR ends

        // normalization
        kinetic_energy *= 0.5; // *later is here

        // equivalent to syncDataReduce in exStamp V1
        // reduce partial sums and share the result
        {
          static constexpr size_t bufsize = 18; // 2x vec3, 3x scalars and 1x 3x3 matrix
          double tmp[bufsize];
          ValueStreamer<double>(tmp) << momentum << kinetic_energy << total_mass << total_particles << hpp_cell << gammaNVTp_cell;
          MPI_Allreduce(MPI_IN_PLACE,tmp,bufsize,MPI_DOUBLE,MPI_SUM,*mpi);
          ValueStreamer<double>(tmp) >> momentum >> kinetic_energy >> total_mass >> total_particles >> hpp_cell >> gammaNVTp_cell;
        }

        // buffers for error computation
        static constexpr size_t NB_ERRVALUES = 11; // 3x3 matrix + 2 scalars
        double old_values[NB_ERRVALUES]; ValueStreamer<double> olds(old_values);
        double new_values[NB_ERRVALUES]; ValueStreamer<double> news(new_values);

        data.hpp = hpp_cell;
        data.m_gammaNVTp = gammaNVTp_cell;

        data.apply_mask();

        /*
        ldbg << "--- end iteration computation ---" << std::endl
             << "hpp       = " << data.hpp << std::endl
             << "gammaNVTp = " << data.m_gammaNVTp << std::endl;
        */
        
        // final temperature computation
        Vec3d tempvec = 2. * ( kinetic_energy - 0.5 * momentum * momentum / total_mass );
        double temperature = ( conv_temperature * ( tempvec.x + tempvec.y + tempvec.z ) / 3. ) / total_particles;

        //ldbg << "temperature : "<<temperature_old<<" => "<<temperature<<std::endl;
        olds << temperature_old;
        news << temperature;

        // convergence step
        data.hpp = (data.hpp - data.m_config.m_Pext * comatrix(data.h) ) / data.m_config.m_masseB;

        // once again, there is a problem in this formula since gammaNVT has a unit ...
        // m_data.hp = (hp_ref + (0.5 * time) * m_data.hpp) / (1. + m_data.gammaNVT);
        data.hp = (hp_ref + (0.5 * dt) * data.hpp) / (1. + data.m_gammaNVT * conv_gammanvt );
        
        //ldbg << "hp          : "<<hp_old<<" => "<<data.hp<<std::endl;
        // olds << hp_old; news << data.hp;

        // update scheme parameters
        data.updateMembers();
        // data.print( ldbg );

        // update gamma (gammaNVTp is an energy), compute error
        double tmp = data.m_config.m_masseB * trace_matrix(data.hpthp) - 3. * total_particles * boltzmann_internal * data.m_config.m_Text;
        data.m_gammaNVTp = (data.m_gammaNVTp + tmp) / data.m_config.m_masseNVT;
        data.m_gammaNVT = gammaNVT_ref + 0.5 * dt * data.m_gammaNVTp;

        /*
        ldbg << "tmp       = "<< tmp<<std::endl
             << "gammaNVTp = "<<data.m_gammaNVTp<<std::endl
             << "gammaNVT  = "<<data.m_gammaNVT<<std::endl;
        */
        // data.print( ldbg );

        // compute residual error
        //ldbg << "gammaNVT    : "<<gammaNVT_old<<" => "<<data.m_gammaNVT<<std::endl;
        olds << gammaNVT_old;
        news << data.m_gammaNVT;

        ssize_t n_err_values = olds.buf - old_values;
        assert( (news.buf-new_values) == n_err_values );
        err = convergence_err_v1( old_values, new_values, n_err_values );
             
        // update temperature
        temperature_old = temperature;
      }
      // main convergence loop PR ends

      ldbg << count << " iterations : residual error = "<<err<< std::endl;

       // compute final velocities
#     pragma omp parallel
      {
        PointerTuple2 ptrs; 
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc)
        {
            size_t i = grid_ijk_to_index( dims , loc + gl );
            int n = cells[i].size();
            cells[i].capture_pointers( ptrs );
  
            auto* __restrict__ vx = ptrs[ field::vx ]; ONIKA_ASSUME_ALIGNED(vx);
            auto* __restrict__ vy = ptrs[ field::vy ]; ONIKA_ASSUME_ALIGNED(vy);
            auto* __restrict__ vz = ptrs[ field::vz ]; ONIKA_ASSUME_ALIGNED(vz);
            
            const auto* __restrict__ ax = ptrs[ field::ax ]; ONIKA_ASSUME_ALIGNED(ax);
            const auto* __restrict__ ay = ptrs[ field::ay ]; ONIKA_ASSUME_ALIGNED(ay);
            const auto* __restrict__ az = ptrs[ field::az ]; ONIKA_ASSUME_ALIGNED(az);
        
#           pragma omp simd
            for(int k=0;k<n;k++)
            {
              const Vec3d f = Vec3d{ax[k], ay[k], az[k]};           
              const Vec3d vSaved = mcri * Vec3d{vx[k], vy[k], vz[k]};
              const Vec3d p = (0.5 * dt) * (saved_prdata.Giht * f) + vSaved;
              const Mat3d q = (0.5 * dt) * saved_prdata.GiGp + (1. + 0.5 * dt * saved_prdata.m_gammaNVT * conv_gammadt) * make_identity_matrix();
              const Vec3d v = inverse(q) * p;
              const Vec3d tmp = mcr * v;
              vx[k] = tmp.x;
              vy[k] = tmp.y;
              vz[k] = tmp.z;
            }
        }
        GRID_OMP_FOR_END
      }

      // restore hpp value
      data.hpp = hpp_ref;
    }

  private:
  
    inline double convergence_err_v2( const double* old_values, const double* new_values, size_t n)
    {
      double err = 0.0;
      for(size_t i=0;i<n;i++)
      {
        double amp = std::max( std::abs(old_values[i]) , std::abs(new_values[i]) );
        if( amp > 0.0 )
        {
          err = std::max( err , std::abs(old_values[i]-new_values[i]) / amp );
        }
      }
      return err;
    }

    inline double convergence_err_v1( const double* old_values, const double* new_values, size_t n)
    {
      double err = 0.0;
      //size_t j=0;
      for(size_t i=0;i<n;i++)
      {
        if( old_values[i] != 0.0 )
        {
          double newerr = std::abs(1. - new_values[i]/old_values[i]);
          //if( newerr > err ) j=i;
          err = std::max( err , newerr );
        }
      }
      // ::exanb::ldbg << "err at i="<<j<<", old="<<old_values[j]<<", new="<<new_values[j]<< std::endl;
      return err;
    }


  };

  template<class GridT> using ConvergencePushParrinelloRahmanTmpl = ConvergencePushParrinelloRahman<GridT>;
  
 // === register factories ===  
  ONIKA_AUTORUN_INIT(convergence_push_parrinellorahman)
  {
   OperatorNodeFactory::instance()->register_factory( "convergence_push_parrinellorahman", make_grid_variant_operator< ConvergencePushParrinelloRahmanTmpl > );
  }

}

