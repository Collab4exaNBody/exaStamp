#include <chrono>
#include <ctime>
#include <mpi.h>
#include <string>
#include <numeric>

#include <onika/math/basic_types_yaml.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <onika/log.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/parallel/random.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <onika/physics/units.h>

#include <onika/math/quaternion_operators.h>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_vx, field::_vy, field::_vz, field::_angmom, /* field::_couple, */ field::_orient, field::_type >
    >
  class InitTemperatureRigidMol : public OperatorNode
  {
    ADD_SLOT( MPI_Comm        , mpi          , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( GridT           , grid         , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies , species      , INPUT_OUTPUT );
    ADD_SLOT( double          , T            , INPUT , REQUIRED );
  public:
    inline void execute () override final
    {
      // MPI Initialization
      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      auto cells = grid->cells();
      IJK dims = grid->dimension();
      size_t ghost_layers = grid->ghost_layers();
      IJK dims_no_ghost = dims - (2*ghost_layers);

      double sum_mass = 0.0;
      double sum_vx   = 0.0;
      double sum_vxc  = 0.0;
      double sum_vy   = 0.0;
      double sum_vyc  = 0.0;
      double sum_vz   = 0.0;
      double sum_vzc  = 0.0;
      double sum_wx   = 0.0;
      double sum_wxc  = 0.0;
      double sum_wy   = 0.0;
      double sum_wyc  = 0.0;
      double sum_wz   = 0.0;
      double sum_wzc  = 0.0;
      double sum_mix  = 0.0;
      double sum_miy  = 0.0;
      double sum_miz  = 0.0;
      double N        = 0.0;
      double nddl_x=0.;
      double nddl_y=0.;
      double nddl_z=0.;

      static const double k = UnityConverterHelper::convert(onika::physics::boltzmann, "J/K");
      const double T           = *(this->T);
      //double sum_nrj=0.0;

      // partie 1
#     pragma omp parallel
      {
        //creation graine pour distribution gaussienne
        auto& re = onika::parallel::random_engine();
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghosts, reduction(+:sum_mass,sum_vx,sum_vxc,sum_vy,sum_vyc,sum_vz,sum_vzc,sum_wx,sum_wxc,sum_wy,sum_wyc,sum_wz,sum_wzc,sum_mix,sum_miy,sum_miz,N,nddl_x,nddl_y,nddl_z) schedule(dynamic) )
        {
          IJK loc = loc_no_ghosts + ghost_layers;
          size_t cell_i = grid_ijk_to_index(dims,loc);

          auto* __restrict__ vx = cells[cell_i][field::vx];
          auto* __restrict__ vy = cells[cell_i][field::vy];
          auto* __restrict__ vz = cells[cell_i][field::vz];
//          const auto* __restrict__ rx = cells[cell_i][field::rx];
//          const auto* __restrict__ ry = cells[cell_i][field::ry];
//          const auto* __restrict__ rz = cells[cell_i][field::rz];
          auto* __restrict__ angmom = cells[cell_i][field::angmom];
          const auto* __restrict__ type_atom = cells[cell_i][field::type];

          size_t n = cells[cell_i].size();

          //boucle d'initialisation aleatoire des quantites
          for(size_t j=0;j<n;j++)
          {
            int t = type_atom[j];
            double mass = species->at(t).m_mass;
            Vec3d minert= species->at(t).m_minert;
            std::normal_distribution<double> f_rand_v(0.0 , std::sqrt(k * T / mass)) ;
//            /* if (minert.x > 0.) */ { std::normal_distribution<double> f_rand_wx(0.0, std::sqrt(k * T / minert.x ) ); }
//            /* if (minert.y > 0.) */ { std::normal_distribution<double> f_rand_wy(0.0, std::sqrt(k * T / minert.y ) ); }
//            /* if (minert.z > 0.) */ { std::normal_distribution<double> f_rand_wz(0.0, std::sqrt(k * T / minert.z ) ); }
            std::normal_distribution<double> f_rand_wx(0.0, std::sqrt(k * T / minert.x ) );
            std::normal_distribution<double> f_rand_wy(0.0, std::sqrt(k * T / minert.y ) );
            std::normal_distribution<double> f_rand_wz(0.0, std::sqrt(k * T / minert.z ) );
            angmom[j] = {0., 0., 0.};
            vx[j] = f_rand_v(re);
            vy[j] = f_rand_v(re);
            vz[j] = f_rand_v(re);
            if (minert.x > 0.) { angmom[j].x = f_rand_wx(re) * minert.x;nddl_x+=1.;}
            if (minert.y > 0.) { angmom[j].y = f_rand_wy(re) * minert.y;nddl_y+=1.;}
            if (minert.z > 0.) { angmom[j].z = f_rand_wz(re) * minert.z;nddl_z+=1.;}

            sum_mass += mass;
            sum_vx   += mass * vx[j];
            sum_vy   += mass * vy[j];
            sum_vz   += mass * vz[j];
            sum_wx   += angmom[j].x;
            sum_wy   += angmom[j].y;
            sum_wz   += angmom[j].z;
            sum_mix  += minert.x;
            sum_miy  += minert.y;
            sum_miz  += minert.z;
            N        += 1.0; // compteur temporaire

          }

        }
        GRID_OMP_FOR_END
      }

      // somme sur tous les processeurs
      MPI_Allreduce(MPI_IN_PLACE,&sum_mass,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_vx,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_vy,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_vz,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_wx,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_wy,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_wz,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_mix,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_miy,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_miz,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&N,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&nddl_x,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&nddl_y,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&nddl_z,1,MPI_DOUBLE,MPI_SUM,*mpi);

      // partie 2
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghosts, reduction(+:sum_vxc,sum_vyc,sum_vzc,sum_wxc,sum_wyc,sum_wzc) schedule(dynamic)  )
        {
          IJK loc = loc_no_ghosts + ghost_layers;
          size_t cell_i = grid_ijk_to_index(dims,loc);

         auto* __restrict__ vx = cells[cell_i][field::vx];
         auto* __restrict__ vy = cells[cell_i][field::vy];
         auto* __restrict__ vz = cells[cell_i][field::vz];
//         auto* __restrict__ rx = cells[cell_i][field::rx];
//         auto* __restrict__ ry = cells[cell_i][field::ry];
//         auto* __restrict__ rz = cells[cell_i][field::rz];
         auto* __restrict__ angmom = cells[cell_i][field::angmom];
//         auto* __restrict__ couple = cells[cell_i][field::couple];
//         const auto* __restrict__ orient = cells[cell_i][field::orient];
         const auto* __restrict__ type_atom = cells[cell_i][field::type];

         size_t n = cells[cell_i].size();

         //boucle d'annulation des grandeurs du CDM
         for(size_t j=0;j<n;j++)
         {
            const int t = type_atom[j];
            const Vec3d minert= species->at(t).m_minert;
            const double mass = species->at(t).m_mass;
            vx[j] = vx[j] - sum_vx / sum_mass;
            vy[j] = vy[j] - sum_vy / sum_mass;
            vz[j] = vz[j] - sum_vz / sum_mass;
            if (sum_mix > 0.) { angmom[j].x = angmom[j].x - sum_wx / sum_mix * minert.x;}
            if (sum_miy > 0.) { angmom[j].y = angmom[j].y - sum_wy / sum_miy * minert.y;}
            if (sum_miz > 0.) { angmom[j].z = angmom[j].z - sum_wz / sum_miz * minert.z;}
            sum_vxc  += mass * vx[j]*vx[j];
            sum_vyc  += mass * vy[j]*vy[j];
            sum_vzc  += mass * vz[j]*vz[j];
            if (minert.x > 0.) {sum_wxc  += angmom[j].x * angmom[j].x / minert.x;}
            if (minert.y > 0.) {sum_wyc  += angmom[j].y * angmom[j].y / minert.y;}
            if (minert.z > 0.) {sum_wzc  += angmom[j].z * angmom[j].z / minert.z;}
          }

        }
        GRID_OMP_FOR_END
      }

      MPI_Allreduce(MPI_IN_PLACE,&sum_vxc,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_vyc,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_vzc,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_wxc,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_wyc,1,MPI_DOUBLE,MPI_SUM,*mpi);
      MPI_Allreduce(MPI_IN_PLACE,&sum_wzc,1,MPI_DOUBLE,MPI_SUM,*mpi);


      {

        /* Modifications LS */
        Vec3d TTcal={0.,0.,0.};
        Vec3d EcTcal={0.,0.,0.};
        Vec3d TRcal={0.,0.,0.};
        Vec3d EcRcal={0.,0.,0.};
        /* Fin Modifications LS */

        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghosts, schedule(dynamic)  )
        {
          IJK loc = loc_no_ghosts + ghost_layers;
          size_t cell_i = grid_ijk_to_index(dims,loc);

         auto* __restrict__ vx = cells[cell_i][field::vx];
         auto* __restrict__ vy = cells[cell_i][field::vy];
         auto* __restrict__ vz = cells[cell_i][field::vz];
//         auto* __restrict__ rx = cells[cell_i][field::rx];
//         auto* __restrict__ ry = cells[cell_i][field::ry];
//         auto* __restrict__ rz = cells[cell_i][field::rz];
         auto* __restrict__ angmom = cells[cell_i][field::angmom];
//         auto* __restrict__ couple = cells[cell_i][field::couple];
         auto* __restrict__ orient = cells[cell_i][field::orient];
         const auto* __restrict__ type_atom = cells[cell_i][field::type];

         size_t n = cells[cell_i].size();

         //boucle de rescale pour avoir la bonne temperature
         if ( T>0.)
         {
           for(size_t j=0;j<n;j++)
           {
              const int t = type_atom[j];
              const Vec3d minert= species->at(t).m_minert;
              double mass = species->at(t).m_mass;

              vx[j] *= std::sqrt( N * k * T / sum_vxc);
              vy[j] *= std::sqrt( N * k * T / sum_vyc);
              vz[j] *= std::sqrt( N * k * T / sum_vzc);


              if (sum_wxc > 0.) { angmom[j].x *= std::sqrt( nddl_x * k * T / sum_wxc );}
              if (sum_wyc > 0.) { angmom[j].y *= std::sqrt( nddl_y * k * T / sum_wyc );}
              if (sum_wzc > 0.) { angmom[j].z *= std::sqrt( nddl_z * k * T / sum_wzc );}

             /* Modifications LS */
             if (minert.x > 0.) {EcRcal.x += 0.5*angmom[j].x*angmom[j].x/minert.x;}
             if (minert.y > 0.) {EcRcal.y += 0.5*angmom[j].y*angmom[j].y/minert.y;}
             if (minert.z > 0.) {EcRcal.z += 0.5*angmom[j].z*angmom[j].z/minert.z;}
             /* Fin Modifications LS */

              Mat3d mat_bf_lab;
              mat_bf_lab.m11 = orient[j].w*orient[j].w + orient[j].x*orient[j].x - orient[j].y*orient[j].y - orient[j].z*orient[j].z;
              mat_bf_lab.m22 = orient[j].w*orient[j].w - orient[j].x*orient[j].x + orient[j].y*orient[j].y - orient[j].z*orient[j].z;
              mat_bf_lab.m33 = orient[j].w*orient[j].w - orient[j].x*orient[j].x - orient[j].y*orient[j].y + orient[j].z*orient[j].z;
              mat_bf_lab.m21 = 2.0 * (orient[j].x*orient[j].y + orient[j].w*orient[j].z ); 
              mat_bf_lab.m12 = 2.0 * (orient[j].x*orient[j].y - orient[j].w*orient[j].z );
              mat_bf_lab.m31 = 2.0 * (orient[j].x*orient[j].z - orient[j].w*orient[j].y );
              mat_bf_lab.m13 = 2.0 * (orient[j].x*orient[j].z + orient[j].w*orient[j].y );
              mat_bf_lab.m32 = 2.0 * (orient[j].y*orient[j].z + orient[j].w*orient[j].x );
              mat_bf_lab.m23 = 2.0 * (orient[j].y*orient[j].z - orient[j].w*orient[j].x );



              //back to lab referential
              angmom[j] = mat_bf_lab*angmom[j];




             /* Modifications LS */
              EcTcal.x += 0.5 * mass *  vx[j] *vx[j];
              EcTcal.y += 0.5 * mass *  vy[j] *vy[j];
              EcTcal.z += 0.5 * mass *  vz[j] *vz[j];
             /* Fin Modifications LS */

            }
          }
        }
        GRID_OMP_FOR_END

        /* Modifications LS */
        MPI_Allreduce(MPI_IN_PLACE,&EcTcal.x,1,MPI_DOUBLE,MPI_SUM,*mpi);
        MPI_Allreduce(MPI_IN_PLACE,&EcTcal.y,1,MPI_DOUBLE,MPI_SUM,*mpi);
        MPI_Allreduce(MPI_IN_PLACE,&EcTcal.z,1,MPI_DOUBLE,MPI_SUM,*mpi);
        MPI_Allreduce(MPI_IN_PLACE,&EcRcal.x,1,MPI_DOUBLE,MPI_SUM,*mpi);
        MPI_Allreduce(MPI_IN_PLACE,&EcRcal.y,1,MPI_DOUBLE,MPI_SUM,*mpi);
        MPI_Allreduce(MPI_IN_PLACE,&EcRcal.z,1,MPI_DOUBLE,MPI_SUM,*mpi);

        TTcal = 2. * EcTcal /(N * k);
        if(nddl_x>0)
                TRcal.x = 2. * EcRcal.x /(nddl_x * k);
        if(nddl_y>0)
                TRcal.y = 2. * EcRcal.y /(nddl_y * k);
        if(nddl_z>0)
                TRcal.z = 2. * EcRcal.z /(nddl_z * k);

        lout<<"    Temperature cible : "<<T<<std::endl;
        lout<<"    Nombre d'atomes : "<<N<<std::endl;
        lout<<"    Nombre de degres de liberte rotationnels : "<<nddl_x<<" "<<nddl_y<<" "<<nddl_z<<std::endl;
        lout<<"    Temperature de translation : "<<TTcal<<std::endl;
        lout<<"    Temperature de rotation : "<<TRcal<<std::endl;

        
        /* Fin Modifications LS */

      }

    }
    
  };

  template<class GridT> using InitTemperatureRigidMolTmpl = InitTemperatureRigidMol<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(init_temperature_rigidmol)
  {
    OperatorNodeFactory::instance()->register_factory("init_temperature_rigidmol", make_grid_variant_operator< InitTemperatureRigidMolTmpl >);
  }

}
