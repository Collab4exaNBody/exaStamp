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

// Builds the SNAP context (parameter/coefficient files, per-material tables, SNA setup) once, in
// init_parameters, so that rcut_max is known before setup_system builds ghosts and neighbor lists.
namespace exaStamp
{
  using namespace exanb;

  class SnapInit : public OperatorNode
  {
    using RealT = double;
    using SnapContext = md::SnapXSContextRealT<RealT>;

    ADD_SLOT( md::SnapParms , parameters      , INPUT , REQUIRED , DocString{"SNAP parameter and coefficient files: { param: <file>, coef: <file> }"} );
    ADD_SLOT( bool          , conv_coef_units , INPUT , false , DocString{"Convert the coefficients from eV to internal energy units"} );
    ADD_SLOT( double        , rcut_max        , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( SnapContext   , snap_ctx        , OUTPUT );

  public:

    inline std::string documentation() const override final
    {
      return R"EOF(

Builds the SNAP context used by snap_force, snap_force_fp64 and compute_descriptor_snap from a SNAP
parameter file and coefficient file, and raises rcut_max to the largest pair cutoff
2*max(radelem)*rcutfac. Place it in init_parameters, after species.

Usage example:

init_parameters:
  - species
  - snap_init:
      parameters: { param: "W.snapparam", coef: "W.snapcoeff" }

)EOF";
    }

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

      // per-pair SNAP cutoff is (radelem[i]+radelem[j])*rcutfac (see snap_bispectrum_op.h), so the
      // largest one is 2*max(radelem)*rcutfac; rcutfac alone is only right when every radelem is 0.5
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
