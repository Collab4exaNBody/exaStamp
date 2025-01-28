



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

#include <onika/cuda/cuda.h>
#include <exanb/compute/compute_cell_particles.h>

//#include "quaternion_rotation.h"
#include <exanb/core/quaternion_operators.h>
#include <exanb/core/quaternion_to_matrix.h>

namespace exaStamp
{
  using namespace exanb;

  struct TorqueToQuaternionComputeFunc
  {
    const double dt = 0.0;
    const ParticleSpecie * __restrict__ species = nullptr;
    
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( const Vec3d& couple, const int t, Quaternion& orient, Vec3d& angmom ) const
    {
      Mat3d mat_lab_bf;
      mat_lab_bf.m11 = orient.w*orient.w + orient.x*orient.x - orient.y*orient.y - orient.z*orient.z;
      mat_lab_bf.m22 = orient.w*orient.w - orient.x*orient.x + orient.y*orient.y - orient.z*orient.z;
      mat_lab_bf.m33 = orient.w*orient.w - orient.x*orient.x - orient.y*orient.y + orient.z*orient.z;
      mat_lab_bf.m12 = 2.0 * (orient.x*orient.y + orient.w*orient.z ); 
      mat_lab_bf.m21 = 2.0 * (orient.x*orient.y - orient.w*orient.z );
      mat_lab_bf.m13 = 2.0 * (orient.x*orient.z - orient.w*orient.y );
      mat_lab_bf.m31 = 2.0 * (orient.x*orient.z + orient.w*orient.y );
      mat_lab_bf.m23 = 2.0 * (orient.y*orient.z + orient.w*orient.x );
      mat_lab_bf.m32 = 2.0 * (orient.y*orient.z - orient.w*orient.x );
      
      //recuperation moment inertiel
      Vec3d minert = species[t].m_minert;
      //calcul du moment angulaire dans le repere mobile
      Vec3d angmom_m = mat_lab_bf * angmom;
      //calcul du couple dans le repere mobile
      Vec3d couple_m = mat_lab_bf * couple;
      //calcul du moment angulaire dans repere fixe a t+dt/2
      angmom += couple*dt/2.0;
      //calcul du moment angulaire dans repere mobile a t+dt/2
      Vec3d omega_m={0.,0.,0.};
      if (minert.x>0.) {omega_m.x=angmom_m.x/minert.x;}
      if (minert.y>0.) {omega_m.y=angmom_m.y/minert.y;}
      if (minert.z>0.) {omega_m.z=angmom_m.z/minert.z;}
      Vec3d angmom_m_ddt = angmom_m + ( couple_m - cross( omega_m, angmom_m) )*dt/2.0;
      //calcul de la derivee du quaternion a t+dt/2
      Quaternion dquat_dt;
      Vec3d omega_m_ddt={0.,0.,0.};
      if (minert.x>0.) {omega_m_ddt.x=angmom_m_ddt.x/minert.x;}
      if (minert.y>0.) {omega_m_ddt.y=angmom_m_ddt.y/minert.y;}
      if (minert.z>0.) {omega_m_ddt.z=angmom_m_ddt.z/minert.z;}
      dquat_dt.w = 0.5 * (-orient.x*omega_m_ddt.x - orient.y*omega_m_ddt.y - orient.z*omega_m_ddt.z );
      dquat_dt.x = 0.5 * (orient.w*omega_m_ddt.x - orient.z*omega_m_ddt.y + orient.y*omega_m_ddt.z );
      dquat_dt.y = 0.5 * (orient.z*omega_m_ddt.x + orient.w*omega_m_ddt.y - orient.x*omega_m_ddt.z );
      dquat_dt.z = 0.5 * (-orient.y*omega_m_ddt.x + orient.x*omega_m_ddt.y + orient.w*omega_m_ddt.z );
      //calcul du quaternion a t+dt/2 
      Quaternion quat_ddt0 = normalize(orient + dquat_dt*dt/2.0);
      double conv = 1.0;
      static constexpr double conv_crit = 1.e-15;
      int niter = 0;
      static constexpr int niter_max = 30;
      while ( conv>conv_crit && niter<niter_max)
      {
        niter += 1;
        //calcul du moment angulaire dans le repere mobile
        mat_lab_bf.m11 = quat_ddt0.w*quat_ddt0.w + quat_ddt0.x*quat_ddt0.x - quat_ddt0.y*quat_ddt0.y - quat_ddt0.z*quat_ddt0.z;
        mat_lab_bf.m22 = quat_ddt0.w*quat_ddt0.w - quat_ddt0.x*quat_ddt0.x + quat_ddt0.y*quat_ddt0.y - quat_ddt0.z*quat_ddt0.z;
        mat_lab_bf.m33 = quat_ddt0.w*quat_ddt0.w - quat_ddt0.x*quat_ddt0.x - quat_ddt0.y*quat_ddt0.y + quat_ddt0.z*quat_ddt0.z;
        mat_lab_bf.m12 = 2.0 * (quat_ddt0.x*quat_ddt0.y + quat_ddt0.w*quat_ddt0.z ); 
        mat_lab_bf.m21 = 2.0 * (quat_ddt0.x*quat_ddt0.y - quat_ddt0.w*quat_ddt0.z );
        mat_lab_bf.m13 = 2.0 * (quat_ddt0.x*quat_ddt0.z - quat_ddt0.w*quat_ddt0.y );
        mat_lab_bf.m31 = 2.0 * (quat_ddt0.x*quat_ddt0.z + quat_ddt0.w*quat_ddt0.y );
        mat_lab_bf.m23 = 2.0 * (quat_ddt0.y*quat_ddt0.z + quat_ddt0.w*quat_ddt0.x );
        mat_lab_bf.m32 = 2.0 * (quat_ddt0.y*quat_ddt0.z - quat_ddt0.w*quat_ddt0.x );
        
        //calcul du moment angulaire dans le repere mobile a t+dt/2
        angmom_m_ddt = mat_lab_bf * angmom;
        //calcul de la vitesse angulaire dans le repere mobile a t+dt/2
        if (minert.x>0.) {omega_m_ddt.x=angmom_m_ddt.x/minert.x;}
        if (minert.y>0.) {omega_m_ddt.y=angmom_m_ddt.y/minert.y;}
        if (minert.z>0.) {omega_m_ddt.z=angmom_m_ddt.z/minert.z;}
        //calcul de la derivee du quaternion a t+dt/2
        dquat_dt.w = 0.5 * (-quat_ddt0.x*omega_m_ddt.x - quat_ddt0.y*omega_m_ddt.y - quat_ddt0.z*omega_m_ddt.z );
        dquat_dt.x = 0.5 * (quat_ddt0.w*omega_m_ddt.x - quat_ddt0.z*omega_m_ddt.y + quat_ddt0.y*omega_m_ddt.z );
        dquat_dt.y = 0.5 * (quat_ddt0.z*omega_m_ddt.x + quat_ddt0.w*omega_m_ddt.y - quat_ddt0.x*omega_m_ddt.z );
        dquat_dt.z = 0.5 * (-quat_ddt0.y*omega_m_ddt.x + quat_ddt0.x*omega_m_ddt.y + quat_ddt0.w*omega_m_ddt.z );
        //calcul du nouveau quaternion
        Quaternion quat_ddt1 = normalize(orient + dquat_dt*dt/2.0);
        //convergence
        conv = norm(quat_ddt1-quat_ddt0); //en realite on pourrait fusionner les deux lignes precedentes et n'avoir qu'un quaternion temporaire
        //stockage du nouveau quaternion obtenu
        quat_ddt0=quat_ddt1;
      }

      if( niter_max == niter && conv>conv_crit )
      {
        /* ABORT */
      }

      //calcul du quaternion a l'instant t+dt
      orient = normalize(orient + dquat_dt*dt);
    }
  };
}

namespace exanb
{
  template<> struct ComputeCellParticlesTraits< exaStamp::TorqueToQuaternionComputeFunc >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{
  using namespace exanb;

  inline std::ostream& operator << (  std::ostream& out , const exanb::Quaternion& q )
  {
    return out <<"("<< q.w<<","<<q.x<<","<<q.y<<","<<q.z<<")";
  }

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_angmom, field::_couple, field::_orient, field::_type >
    >
  class TorqueToQuaternionRigidMol : public OperatorNode
  {
    //ADD_SLOT( MPI_Comm        , mpi          , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( GridT           , grid         , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies , species      , INPUT_OUTPUT );
    ADD_SLOT( double          , dt           , INPUT , REQUIRED );
    
    static constexpr FieldSet< field::_couple, field::_type, field::_orient, field::_angmom > compute_field_set{};
    
  public:
    inline void execute () override final
    {
      compute_cell_particles( *grid , false , TorqueToQuaternionComputeFunc{*dt,species->data()} , compute_field_set , parallel_execution_context() );

#if 0
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      size_t ghost_layers = grid->ghost_layers();
      IJK dims_no_ghost = dims - (2*ghost_layers);
      const double dt           = *(this->dt);

      // partie 1
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghosts, schedule(dynamic) )
        {
          IJK loc = loc_no_ghosts + ghost_layers;
          size_t cell_i = grid_ijk_to_index(dims,loc);

          auto* __restrict__ angmom = cells[cell_i][field::angmom];
          auto* __restrict__ orient = cells[cell_i][field::orient];
          const auto* __restrict__ couple = cells[cell_i][field::couple];
          const auto* __restrict__ type_atom = cells[cell_i][field::type];

          size_t n = cells[cell_i].size();

          for(size_t j=0;j<n;j++)
          {
            int t = type_atom[j];
            Mat3d mat_lab_bf;
            mat_lab_bf.m11 = orient[j].w*orient[j].w + orient[j].x*orient[j].x - orient[j].y*orient[j].y - orient[j].z*orient[j].z;
            mat_lab_bf.m22 = orient[j].w*orient[j].w - orient[j].x*orient[j].x + orient[j].y*orient[j].y - orient[j].z*orient[j].z;
            mat_lab_bf.m33 = orient[j].w*orient[j].w - orient[j].x*orient[j].x - orient[j].y*orient[j].y + orient[j].z*orient[j].z;
            mat_lab_bf.m12 = 2.0 * (orient[j].x*orient[j].y + orient[j].w*orient[j].z ); 
            mat_lab_bf.m21 = 2.0 * (orient[j].x*orient[j].y - orient[j].w*orient[j].z );
            mat_lab_bf.m13 = 2.0 * (orient[j].x*orient[j].z - orient[j].w*orient[j].y );
            mat_lab_bf.m31 = 2.0 * (orient[j].x*orient[j].z + orient[j].w*orient[j].y );
            mat_lab_bf.m23 = 2.0 * (orient[j].y*orient[j].z + orient[j].w*orient[j].x );
            mat_lab_bf.m32 = 2.0 * (orient[j].y*orient[j].z - orient[j].w*orient[j].x );
            
            //recuperation moment inertiel
            Vec3d minert = species->at(t).m_minert;
            //calcul du moment angulaire dans le repere mobile
            Vec3d angmom_m = mat_lab_bf * angmom[j];
            //calcul du couple dans le repere mobile
            Vec3d couple_m = mat_lab_bf * couple[j];
            //calcul du moment angulaire dans repere fixe a t+dt/2
            angmom[j] += couple[j]*dt/2.0;
            //calcul du moment angulaire dans repere mobile a t+dt/2
            Vec3d omega_m={0.,0.,0.};
            if (minert.x>0.) {omega_m.x=angmom_m.x/minert.x;}
            if (minert.y>0.) {omega_m.y=angmom_m.y/minert.y;}
            if (minert.z>0.) {omega_m.z=angmom_m.z/minert.z;}
            Vec3d angmom_m_ddt = angmom_m + ( couple_m - cross( omega_m, angmom_m) )*dt/2.0;
            //calcul de la derivee du quaternion a t+dt/2
            Quaternion dquat_dt;
            Vec3d omega_m_ddt={0.,0.,0.};
            if (minert.x>0.) {omega_m_ddt.x=angmom_m_ddt.x/minert.x;}
            if (minert.y>0.) {omega_m_ddt.y=angmom_m_ddt.y/minert.y;}
            if (minert.z>0.) {omega_m_ddt.z=angmom_m_ddt.z/minert.z;}
            dquat_dt.w = 0.5 * (-orient[j].x*omega_m_ddt.x - orient[j].y*omega_m_ddt.y - orient[j].z*omega_m_ddt.z );
            dquat_dt.x = 0.5 * (orient[j].w*omega_m_ddt.x - orient[j].z*omega_m_ddt.y + orient[j].y*omega_m_ddt.z );
            dquat_dt.y = 0.5 * (orient[j].z*omega_m_ddt.x + orient[j].w*omega_m_ddt.y - orient[j].x*omega_m_ddt.z );
            dquat_dt.z = 0.5 * (-orient[j].y*omega_m_ddt.x + orient[j].x*omega_m_ddt.y + orient[j].w*omega_m_ddt.z );
            //calcul du quaternion a t+dt/2 
            Quaternion quat_ddt0 = normalize(orient[j] + dquat_dt*dt/2.0);
            double conv = 1.0;
            static constexpr double conv_crit = 1.e-15;
            int niter = 0;
            static constexpr int niter_max = 30;
            while ( conv>conv_crit && niter<niter_max)
            {
              niter += 1;
              //calcul du moment angulaire dans le repere mobile
              mat_lab_bf.m11 = quat_ddt0.w*quat_ddt0.w + quat_ddt0.x*quat_ddt0.x - quat_ddt0.y*quat_ddt0.y - quat_ddt0.z*quat_ddt0.z;
              mat_lab_bf.m22 = quat_ddt0.w*quat_ddt0.w - quat_ddt0.x*quat_ddt0.x + quat_ddt0.y*quat_ddt0.y - quat_ddt0.z*quat_ddt0.z;
              mat_lab_bf.m33 = quat_ddt0.w*quat_ddt0.w - quat_ddt0.x*quat_ddt0.x - quat_ddt0.y*quat_ddt0.y + quat_ddt0.z*quat_ddt0.z;
              mat_lab_bf.m12 = 2.0 * (quat_ddt0.x*quat_ddt0.y + quat_ddt0.w*quat_ddt0.z ); 
              mat_lab_bf.m21 = 2.0 * (quat_ddt0.x*quat_ddt0.y - quat_ddt0.w*quat_ddt0.z );
              mat_lab_bf.m13 = 2.0 * (quat_ddt0.x*quat_ddt0.z - quat_ddt0.w*quat_ddt0.y );
              mat_lab_bf.m31 = 2.0 * (quat_ddt0.x*quat_ddt0.z + quat_ddt0.w*quat_ddt0.y );
              mat_lab_bf.m23 = 2.0 * (quat_ddt0.y*quat_ddt0.z + quat_ddt0.w*quat_ddt0.x );
              mat_lab_bf.m32 = 2.0 * (quat_ddt0.y*quat_ddt0.z - quat_ddt0.w*quat_ddt0.x );
              //calcul du moment angulaire dans le repere mobile a t+dt/2
              angmom_m_ddt = mat_lab_bf * angmom[j];
              //calcul de la vitesse angulaire dans le repere mobile a t+dt/2
              if (minert.x>0.) {omega_m_ddt.x=angmom_m_ddt.x/minert.x;}
              if (minert.y>0.) {omega_m_ddt.y=angmom_m_ddt.y/minert.y;}
              if (minert.z>0.) {omega_m_ddt.z=angmom_m_ddt.z/minert.z;}
              //calcul de la derivee du quaternion a t+dt/2
              dquat_dt.w = 0.5 * (-quat_ddt0.x*omega_m_ddt.x - quat_ddt0.y*omega_m_ddt.y - quat_ddt0.z*omega_m_ddt.z );
              dquat_dt.x = 0.5 * (quat_ddt0.w*omega_m_ddt.x - quat_ddt0.z*omega_m_ddt.y + quat_ddt0.y*omega_m_ddt.z );
              dquat_dt.y = 0.5 * (quat_ddt0.z*omega_m_ddt.x + quat_ddt0.w*omega_m_ddt.y - quat_ddt0.x*omega_m_ddt.z );
              dquat_dt.z = 0.5 * (-quat_ddt0.y*omega_m_ddt.x + quat_ddt0.x*omega_m_ddt.y + quat_ddt0.w*omega_m_ddt.z );
              //calcul du nouveau quaternion
              Quaternion quat_ddt1 = normalize(orient[j] + dquat_dt*dt/2.0);
              //convergence
              conv = norm(quat_ddt1-quat_ddt0); //en realite on pourrait fusionner les deux lignes precedentes et n'avoir qu'un quaternion temporaire
              //stockage du nouveau quaternion obtenu
              quat_ddt0=quat_ddt1;
            }

	          if( niter_max == niter && conv>conv_crit )
            {
              std::cerr << "Convergence Quaternion "<<conv<<">"<<conv_crit<<" , quat["<<j<<"]="<<orient[j]<<" , couple="<<couple_m<<" , dquat_dt="<<dquat_dt<<std::endl;
              std::cerr << "Cell @"<<loc<<" , part #"<<j<<" , dims="<<dims<<", gl="<<ghost_layers<<" , offset="<<grid->offset()<<std::endl << std::flush;
              for(size_t k=0;k<j;k++)
              {
                std::cerr << "\t"<<k<<" : orient="<<orient[k]<<" , couple="<<couple[k]<<std::endl;                
              }
              std::cerr <<std::endl << std::flush;
              std::abort();
	          }

            //calcul du quaternion a l'instant t+dt
            orient[j] = normalize(orient[j] + dquat_dt*dt);
          }
        }
        GRID_OMP_FOR_END
      }
#endif

    }
    
  };

  // === register factories ===
  template<class GridT> using TorqueToQuaternionRigidMolTmpl = TorqueToQuaternionRigidMol<GridT>;

  ONIKA_AUTORUN_INIT(torque_to_quaternion)
  {
    OperatorNodeFactory::instance()->register_factory("torque_to_quaternion", make_grid_variant_operator< TorqueToQuaternionRigidMolTmpl >);
  }

}
