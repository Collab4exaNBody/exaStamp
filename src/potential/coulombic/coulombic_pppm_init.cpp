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

#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_stream.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/domain.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <exaStamp/potential/coulombic/pppm.h>
#include <mpi.h>
#include <vector>

namespace exaStamp
{
inline namespace coulombic_ewald
{
  using namespace exanb;

  class PPPMInitOperator : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( double            , accuracy_relative , INPUT , 1.0e-5 , DocString{"relative rms force accuracy (relative to the force between two unit charges at 1 ang)"} );
    ADD_SLOT( double            , g_ewald           , INPUT , 0.0 , DocString{"Ewald splitting parameter, 0 = automatic"} );
    ADD_SLOT( double            , radius            , INPUT , REQUIRED , DocString{"real space cutoff"} );
    ADD_SLOT( std::vector<long> , mesh              , INPUT , std::vector<long>{0,0,0} , DocString{"mesh points in each direction, 0 0 0 = automatic"} );
    ADD_SLOT( long              , order             , INPUT , 5 , DocString{"charge assignment order, 2 to 7"} );
    ADD_SLOT( Domain            , domain            , INPUT , OPTIONAL );
    ADD_SLOT( double            , sum_square_charge , INPUT );
    ADD_SLOT( double            , sum_charge        , INPUT );
    ADD_SLOT( uint64_t          , natoms            , INPUT , OPTIONAL );
    ADD_SLOT( MPI_Comm          , mpi               , INPUT );
    ADD_SLOT( PPPMParameters    , pppm_config       , INPUT_OUTPUT );
    ADD_SLOT( EwaldParameters   , ewald_config      , INPUT_OUTPUT , DocString{"real space parameters, for coulombic_ewald_short_range"} );
    ADD_SLOT( double            , rcut              , OUTPUT );
    ADD_SLOT( double            , rcut_max          , INPUT_OUTPUT , 0.0 );

  public:

    inline void execute () override final
    {
      int rank=0;
      MPI_Comm_rank(*mpi, &rank);

      if( mesh->size() != 3 )
      {
        fatal_error() << "coulombic_pppm_init : mesh must have 3 values" << std::endl;
      }
      const long mesh_user[3] = { (*mesh)[0] , (*mesh)[1] , (*mesh)[2] };

      if( domain.has_value() )
      {
        auto & p = *pppm_config;
        const Mat3d cell = ewald_cell_matrix( domain->xform() , domain->bounds_size() );

        // called before the system is set up (e.g. in init_parameters, to get rcut) : nothing to initialize yet
        const bool empty_domain = determinant( cell ) == 0.0;

        // all inputs are global (identical on all ranks) : every rank computes the same parameters, no broadcast.
        const bool need_init = ! empty_domain && (
             p.volume == 0.0
          || *g_ewald != p.g_ewald_user
          || *radius != p.radius
          || *accuracy_relative != p.accuracy_relative
          || *order != p.order
          || mesh_user[0] != p.mesh_user[0] || mesh_user[1] != p.mesh_user[1] || mesh_user[2] != p.mesh_user[2] );

        // cell change only (NPT, deformation) : keep mesh and g_ewald, update volume dependent quantities (LAMMPS PPPM::setup)
        const bool need_setup = ! empty_domain && ! need_init && ! ewald_same_cell( cell , p.cell );

        if( need_init )
        {
          if( ! ( domain->periodic_boundary_x() && domain->periodic_boundary_y() && domain->periodic_boundary_z() ) )
          {
            fatal_error() << "Domain must be entierly periodic, cannot initialize PPPM." << std::endl;
          }
          if( ! natoms.has_value() )
          {
            fatal_error() << "coulombic_pppm_init : natoms is required" << std::endl;
          }

          const bool first_init = ( p.volume == 0.0 );
          pppm_init_parameters( *g_ewald , *radius , *accuracy_relative , *order , mesh_user , cell , *natoms , *sum_square_charge , *sum_charge , p );

          if( rank == 0 && first_init )
          {
            lout << "====== PPPM configuration ======" << std::endl;
            lout << "cell    = "<< cell << std::endl;
            lout << "g_ewald = "<< p.g_ewald << std::endl;
            lout << "radius  = "<< p.radius << std::endl;
            lout << "mesh    = "<< p.nx <<" "<< p.ny <<" "<< p.nz << std::endl;
            lout << "order   = "<< p.order << std::endl;
            lout << "accuracy_relative  = "<< p.accuracy_relative << std::endl;
            lout << "estimated accuracy = "<< p.estimated_accuracy << " eV/ang (relative "<< p.estimated_accuracy / COULOMB_CONSTANT_EV_ANG <<")" << std::endl;
            lout << "qsum    = "<< p.qsum << std::endl;
            lout << "volume  = "<< p.volume << std::endl;
            lout << "================================" << std::endl;
          }
          else
          {
            ldbg << "PPPM re-initialized : cell="<<cell<<" g_ewald="<<p.g_ewald<<" mesh="<<p.nx<<"x"<<p.ny<<"x"<<p.nz<<std::endl;
          }
        }
        else if( need_setup )
        {
          pppm_setup( p , cell );
          ldbg << "PPPM setup for new cell "<<cell<<" volume="<<p.volume<<std::endl;
        }

        // real space part is the Ewald one, with PPPM's g_ewald
        if( p.volume > 0.0 )
        {
          auto & e = *ewald_config;
          e.g_ewald = p.g_ewald;
          e.radius = p.radius;
          e.accuracy_relative = p.accuracy_relative;
          e.qsum = p.qsum;
          e.volume = p.volume;
          e.cell = p.cell;
          e.gm_sr = ewald_constants::qqr2e;
          e.qqr2e = ewald_constants::qqr2e;
          e.bt_sr = 2. * p.g_ewald / std::sqrt(M_PI);
          e.nk = e.nknz = 0; // no k vectors : coulombic_ewald_long_range must not be used with this configuration
          e.Gdata.clear();
        }
      }

      *rcut = *radius;
      *rcut_max = std::max( *rcut_max , *radius );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(
Initializes PPPM long range coulomb (coulombic_pppm), same algorithm and parameter choice as LAMMPS kspace_style pppm
(ik differentiation). Works on orthogonal and triclinic periodic cells. Also fills ewald_config with the real space
parameters (g_ewald, radius) used by coulombic_ewald_short_range. Needs sum_square_charge, sum_charge and natoms
(sum_charges operator). When only the cell changes, mesh and g_ewald are kept and the influence function is rebuilt.
)EOF";
    }

  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(coulombic_pppm_init)
  {
    OperatorNodeFactory::instance()->register_factory( "coulombic_pppm_init" , make_simple_operator<PPPMInitOperator> );
  }

}
}
