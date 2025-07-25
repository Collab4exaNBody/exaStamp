/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

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
#include <exaStamp/unit_system.h>

#include <onika/math/quaternion_operators.h>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_vx, field::_vy, field::_vz, field::_angmom /* , field::_couple */ , field::_orient, field::_type >
    >
  class LangevinRigidMol : public OperatorNode
  {
    ADD_SLOT( MPI_Comm        , mpi            , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( GridT           , grid           , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies , species        , INPUT_OUTPUT );
    ADD_SLOT( double          , dt             , INPUT );
    ADD_SLOT( double          , T              , INPUT , REQUIRED );
    ADD_SLOT( double          , friction       , INPUT , 1.0 );
    ADD_SLOT( double          , friction_ratio , INPUT , 1.e-2 );
    ADD_SLOT( double          , time_nve,        INPUT , 1e30 );    
    ADD_SLOT( int64_t         , timestep       , INPUT , REQUIRED );

  public:
    inline void execute () override final
    {
      static constexpr double k = EXASTAMP_CONST_QUANTITY( onika::physics::boltzmann * J / K );
    
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      size_t ghost_layers = grid->ghost_layers();
      IJK dims_no_ghost = dims - (2*ghost_layers);
      const double dt           = *(this->dt);
      const double T           = *(this->T);
      const double xsi       = *(this->friction);
      const double friction_ratio       = *(this->friction_ratio);
      const double xsir = xsi * friction_ratio;
      const double nvetime           = *(this->time_nve);
      const auto * __restrict__ species_ptr = species->data();

      double curtime           = dt * (*timestep);

      if(curtime < nvetime) {

        ldbg<<" ###### Langevin #######"<<std::endl<<std::flush;

      // partie 1
#     pragma omp parallel
      {
        //creation graine pour distribution gaussienne
        auto& re = onika::parallel::random_engine();
        std::normal_distribution<double> f_rand(0.0 , 1.0) ;
        GRID_OMP_FOR_BEGIN(dims_no_ghost,_,loc_no_ghosts, schedule(dynamic) )
        {
          IJK loc = loc_no_ghosts + ghost_layers;
          size_t cell_i = grid_ijk_to_index(dims,loc);

          auto* __restrict__ vx = cells[cell_i][field::vx];
          auto* __restrict__ vy = cells[cell_i][field::vy];
          auto* __restrict__ vz = cells[cell_i][field::vz];
          auto* __restrict__ angmom = cells[cell_i][field::angmom];
          auto* __restrict__ orient = cells[cell_i][field::orient];
//          const auto* __restrict__ couple = cells[cell_i][field::couple];
          const auto* __restrict__ type_atom = cells[cell_i][field::type];

          size_t n = cells[cell_i].size();

          for(size_t j=0;j<n;j++)
          {
            const int t = type_atom[j];
            //recuperation moment inertiel et de la masse
            const double mass = species_ptr[t].m_mass;
            //calcul du parametre alpha pour la vitesse du CM
            const double alphaCM = std::exp( - xsi * dt / mass ) ;
            //calcul du terme de dissipation pour le CM
            const double sigmaCM_langevin = std::sqrt( k * T * (1.-alphaCM*alphaCM) / mass ) ;

            //application du thermostat sur les vitesses de CM
            //creation vecteur aleatoire
            //distribution uniforme, on ajoute un test sur la norme du vecteur (http://corysimon.github.io/articles/uniformdistn-on-sphere/)
            Vec3d vec_alea_v={0.,0.,0.};
            while (norm(vec_alea_v)< 1.e-6)
            {
              vec_alea_v.x = f_rand(re);
              vec_alea_v.y = f_rand(re);
              vec_alea_v.z = f_rand(re);
            }
            vx[j] = alphaCM * vx[j] + sigmaCM_langevin * vec_alea_v.x;
            vy[j] = alphaCM * vy[j] + sigmaCM_langevin * vec_alea_v.y;
            vz[j] = alphaCM * vz[j] + sigmaCM_langevin * vec_alea_v.z;

            if( species_ptr[t].m_rigid_atom_count > 1 )
            {
              const Vec3d minert = species_ptr[t].m_minert;
              //calcul du parametre alpha pour le moment angulaire
              Vec3d alphaL = {0. , 0. , 0.};
              Vec3d sigmaL_langevin = {0. , 0. , 0.};
              //calcul du terme de dissipation pour le moment angulaire
              if (minert.x > 0.) 
              {
                alphaL.x = std::exp(-xsir * dt / minert.x);
                sigmaL_langevin.x = std::sqrt( k * T * (1. - alphaL.x * alphaL.x) / minert.x);
              }
              if (minert.y > 0.) 
              { 
                alphaL.y = std::exp(-xsir * dt / minert.y);
                sigmaL_langevin.y = std::sqrt( k * T * (1. - alphaL.y * alphaL.y) / minert.y);
              }
              if (minert.z > 0.)
              { 
                alphaL.z = std::exp(-xsir * dt / minert.z);
                sigmaL_langevin.z = std::sqrt( k * T * (1. - alphaL.z * alphaL.z) / minert.z);
              }
                          
              //creation vecteur aleatoire
              Vec3d vec_alea_L={0., 0., 0.};
              while (norm(vec_alea_L)<1.e-6)
              {
                vec_alea_L.x = f_rand(re);
                vec_alea_L.y = f_rand(re);
                vec_alea_L.z = f_rand(re);
              }
              //application du thermostat sur le moment angulaire
              Vec3d angmom_m;
              //calcul du moment angulaire dans le repere mobile
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
              angmom_m = mat_lab_bf * angmom[j];
              //application du thermostat
              angmom_m.x = alphaL.x * angmom_m.x + minert.x * sigmaL_langevin.x * vec_alea_L.x;
              angmom_m.y = alphaL.y * angmom_m.y + minert.y * sigmaL_langevin.y * vec_alea_L.y;
              angmom_m.z = alphaL.z * angmom_m.x + minert.z * sigmaL_langevin.z * vec_alea_L.z;
              //retour dans la base du labo
              Mat3d mat_bf_lab = transpose(mat_lab_bf);
              angmom[j] = mat_bf_lab * angmom_m;
            }
            
          }
        }
        GRID_OMP_FOR_END
      }
       
     } else {
        ldbg<<" ###### NVE #######"<<std::endl<<std::flush;

     }
    }
    
  };

  template<class GridT> using LangevinRigidMolTmpl = LangevinRigidMol<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(langevin_rigidmol)
  {
    OperatorNodeFactory::instance()->register_factory("langevin_rigidmol", make_grid_variant_operator< LangevinRigidMolTmpl >);
  }

}
