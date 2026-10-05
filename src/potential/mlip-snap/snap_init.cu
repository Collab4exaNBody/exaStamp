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

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/log.h>
#include <onika/file_utils.h>

#include <md/snap/snap_params.h>
#include <md/snap/snap_read_lammps.h>
#include <md/snap/snap_config.h>
#include <md/snap/snap_context.h>
#include <md/snap/sna.h>

#include <algorithm>
#include <string>

// Early (init_parameters, before setup_system) SNAP context construction -- mirrors pod_init.cu /
// mtp_init.cu. Moved out of compute_descriptor_snap.cu (which used to build snap_ctx lazily inside
// its own execute(), too late for rcut_max to widen ghost/neighbor setup in time -- every .msp using
// it needed a manual `rcut_max: X ang` workaround under `global:`). Building it here once, before
// setup_system, sets rcut_max in time and removes that workaround, matching pod_init/mtp_init's
// existing pattern. compute_descriptor_snap.cu's OWN alternate constant-neighbor-count mode
// (nneigh_bispectrum/closest_bispectrum/neigh_margin) is unaffected -- those stay as standalone
// slots on that operator, never part of SnapContext.
namespace exaStamp
{
  using namespace exanb;

  class SnapInit : public OperatorNode
  {
    using RealT = double;
    using SnapContext = md::SnapXSContextRealT<RealT>;

    ADD_SLOT( md::SnapParms , parameters      , INPUT , REQUIRED , DocString{"LAMMPS-format SNAP parameter/coefficient files (param/coef), see compute_descriptor_snap/snap_force"} );
    ADD_SLOT( bool          , conv_coef_units , INPUT , false );
    ADD_SLOT( double        , rcut_max        , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( SnapContext   , snap_ctx        , OUTPUT );

  public:

    inline void execute() override final
    {
      ldbg << "Initializing SNAP potential" << std::endl;

      std::string lammps_param = onika::data_file_path( parameters->lammps_param );
      std::string lammps_coef  = onika::data_file_path( parameters->lammps_coef );
      SnapExt::snap_read_lammps( lammps_param, lammps_coef, snap_ctx->m_config, *conv_coef_units );

      const int nmat = snap_ctx->m_config.materials().size();
      snap_ctx->m_factor.assign( nmat, 1.0 );
      snap_ctx->m_radelem.assign( nmat, 0.0 );
      int cnt = 0;
      double max_radelem = 0.0;
      for ( const auto& mat : snap_ctx->m_config.materials() )
      {
        snap_ctx->m_factor[cnt]  = mat.weight();
        snap_ctx->m_radelem[cnt] = mat.radelem();
        max_radelem = std::max( max_radelem, mat.radelem() );
        cnt++;
      }

      // Real per-pair SNAP cutoff is (radelem[i]+radelem[j])*rcutfac (see BispectrumOpRealT's own
      // cutij, snap_bispectrum_op.h / real LAMMPS PairSNAP::init_one), NOT the bare rcutfac scale
      // factor -- rcutfac alone only equals the true cutoff when every material's radelem==0.5
      // (LAMMPS's convention for potentials with no per-element radii, true of every SNAP test
      // asset used so far -- Ta/WBe -- but not real per-element-radius potentials like InP, where
      // this silently produced a neighbor list too narrow to ever see a real neighbor). Worst-case
      // pair cutoff over all i,j is 2*max_i(radelem[i]) (max_i,j(r_i+r_j) = 2*max_i(r_i)).
      snap_ctx->m_rcut = snap_ctx->m_config.rcutfac() * 2.0 * max_radelem;
      *rcut_max = std::max( double(*rcut_max), double(snap_ctx->m_rcut) );
      ldbg << "SNAP cutoff radius: " << snap_ctx->m_rcut << std::endl;

      snap_ctx->sna = new SnapInternal::SNARealT<RealT>( new SnapInternal::Memory()
                          , snap_ctx->m_config.rfac0(), snap_ctx->m_config.twojmax(), snap_ctx->m_config.rmin0()
                          , snap_ctx->m_config.switchflag(), snap_ctx->m_config.bzeroflag(), snap_ctx->m_config.chemflag()
                          , snap_ctx->m_config.bnormflag(), snap_ctx->m_config.wselfallflag(), snap_ctx->m_config.nelements()
                          , snap_ctx->m_config.switchinnerflag() );
      snap_ctx->sna->init();
    }

  };

  ONIKA_AUTORUN_INIT(snap_init)
  {
    OperatorNodeFactory::instance()->register_factory("snap_init", make_simple_operator<SnapInit>);
  }

}
