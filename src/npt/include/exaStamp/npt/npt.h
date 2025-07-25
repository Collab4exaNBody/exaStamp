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

#pragma once

#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/core/domain.h>
#include <onika/print_utils.h>

namespace exaStamp
{

  struct NPTContext
  {
    using Mat3d = exanb::Mat3d;

    std::string algo;
    std::string pmode;    
    
    int dimension, which;
    double dtv, dtf, dthalf, dt4, dt8, dto;
    double boltz, nktv2p, tdof;
    double vol0;    // reference volume
    double t0;      // reference temperature

    // For thermostat purposes
    double t_start, t_stop;
    double t_current, t_target, ke_target;
    double t_freq;
    int tstat_flag;
    long start_it;
    
    // For barostat purposes
    std::string pstyle,pcouple;
    int p_flag[6];    // 1 if control P on this dim, 0 if not
    double p_start[6], p_stop[6];
    double p_period[6], p_freq[6], p_target[6];
    double omega[6], omega_dot[6];
    double omega_mass[6];
    double p_current[6];
    double drag;                     // drag factor specified by user (default = 0)
    double tdrag_factor;             // drag factor on thermostat
    double pdrag_factor;             // drag factor on barostat
    int pstat_flag;    // 1 if control P

    double p_temp;    // target temperature for barostat
    int p_temp_flag;

    double *eta, *eta_dot;    // chain thermostat for particles
    double *eta_dotdot;
    double *eta_mass;
    int mtchain;                 // length of chain
    int mtchain_default_flag;    // 1 = mtchain is default

    double *etap;    // chain thermostat for barostat
    double *etap_dot;
    double *etap_dotdot;
    double *etap_mass;
    int mpchain;    // length of chain

    bool mtk_flag;         // 0 if using Hoover barostat
    int pdim;             // number of barostatted dims
    double p_freq_max;    // maximum barostat frequency

    double p_hydro;    // hydrostatic target pressure

    int nc_tchain, nc_pchain;
    double factor_eta;
    double sigma[6];        // scaled target stress
    double fdev[6];         // deviatoric force on barostat
    int deviatoric_flag;    // 0 if target stress tensor is hydrostatic
    double h0_inv[6];       // h_inv of reference (zero strain) box
    int nreset_h0;          // interval for resetting h0

    double mtk_term1, mtk_term2;    // Martyna-Tobias-Klein corrections

    int eta_mass_flag;      // 1 if eta_mass updated, 0 if not.
    int omega_mass_flag;    // 1 if omega_mass updated, 0 if not.
    int etap_mass_flag;     // 1 if etap_mass updated, 0 if not.
    int dipole_flag;        // 1 if dipole is updated, 0 if not.
    int dlm_flag;           // 1 if using the DLM rotational integrator, 0 if not

    int scaleyz;     // 1 if yz scaled with lz
    int scalexz;     // 1 if xz scaled with lz
    int scalexy;     // 1 if xy scaled with ly

    double fixedpoint[3];    // location of dilation fixed-point

    inline void update_target_T_KE(long cur_it, long end_it)
    {
      double delta = cur_it - this->start_it;
      if (delta > 0) delta /= (end_it - this->start_it);
      this->t_target = this->t_start + delta * (this->t_stop-this->t_start);
      this->ke_target = this->tdof * this->boltz * this->t_target;      
    }

    inline void update_target_P(long cur_it, long end_it)
    {
      double delta = cur_it - this->start_it;
      if (delta > 0) delta /= (end_it - this->start_it);
      this->p_hydro = 0.0;
      for (int i = 0; i < 3; i++)
        if (this->p_flag[i]) {
          this->p_target[i] = this->p_start[i] + delta * (this->p_stop[i] - this->p_start[i]);
          this->p_hydro += this->p_target[i];
        }
      if (this->pdim > 0) this->p_hydro /= this->pdim;
      if (this->pstyle == "TRICLINIC")
        for (int i = 3; i < 6; i++)
          this->p_target[i] = this->p_start[i] + delta * (this->p_stop[i] - this->p_start[i]);
    }

    inline void update_sigma()
    {
      // generate upper-triangular half of
      // sigma = vol0*h0inv*(p_target-p_hydro)*h0inv^t
      // units of sigma are are PV/L^2 e.g. atm.A
      //
      // [ 0 5 4 ]   [ 0 5 4 ] [ 0 5 4 ] [ 0 - - ]
      // [ 5 1 3 ] = [ - 1 3 ] [ 5 1 3 ] [ 5 1 - ]
      // [ 4 3 2 ]   [ - - 2 ] [ 4 3 2 ] [ 4 3 2 ]
      ///      ldbg << "VOL0 = " << this->vol0 << std::endl;
      this->sigma[0] = this->vol0*(this->h0_inv[0]*((this->p_target[0]-this->p_hydro)*this->h0_inv[0] + this->p_target[5]*this->h0_inv[5]+this->p_target[4]*this->h0_inv[4]) + this->h0_inv[5]*(this->p_target[5]*this->h0_inv[0] + (this->p_target[1]-this->p_hydro)*this->h0_inv[5]+this->p_target[3]*this->h0_inv[4]) + this->h0_inv[4]*(this->p_target[4]*this->h0_inv[0]+this->p_target[3]*this->h0_inv[5] + (this->p_target[2]-this->p_hydro)*this->h0_inv[4]));
      this->sigma[1] = this->vol0*(this->h0_inv[1]*((this->p_target[1]-this->p_hydro)*this->h0_inv[1] + this->p_target[3]*this->h0_inv[3]) + this->h0_inv[3]*(this->p_target[3]*this->h0_inv[1] + (this->p_target[2]-this->p_hydro)*this->h0_inv[3]));
      this->sigma[2] = this->vol0*(this->h0_inv[2]*((this->p_target[2]-this->p_hydro)*this->h0_inv[2]));
      this->sigma[3] = this->vol0*(this->h0_inv[1]*(this->p_target[3]*this->h0_inv[2]) + this->h0_inv[3]*((this->p_target[2]-this->p_hydro)*this->h0_inv[2]));
      this->sigma[4] = this->vol0*(this->h0_inv[0]*(this->p_target[4]*this->h0_inv[2]) + this->h0_inv[5]*(this->p_target[3]*this->h0_inv[2]) + this->h0_inv[4]*((this->p_target[2]-this->p_hydro)*this->h0_inv[2]));
      this->sigma[5] = this->vol0*(this->h0_inv[0]*(this->p_target[5]*this->h0_inv[1]+this->p_target[4]*this->h0_inv[3]) + this->h0_inv[5]*((this->p_target[1]-this->p_hydro)*this->h0_inv[1]+this->p_target[3]*this->h0_inv[3]) + this->h0_inv[4]*(this->p_target[3]*this->h0_inv[1]+(this->p_target[2]-this->p_hydro)*this->h0_inv[3]));
    }

    inline void update_fdev(Mat3d h)
    {
      // ldbg << "h[0] = " << h.m11 << std::endl;
      // ldbg << "h[1] = " << h.m22 << std::endl;
      // ldbg << "h[2] = " << h.m33 << std::endl;
      // ldbg << "h[3] = " << h.m23 << std::endl;
      // ldbg << "h[4] = " << h.m13 << std::endl;
      // ldbg << "h[5] = " << h.m12 << std::endl;      

      // for (int i = 0; i < 6; i++) {
      //   ldbg << "sigma["<<i<<"] = " << this->sigma[i] << std::endl;
      // }      
      
      this->fdev[0] =
        h.m11*(this->sigma[0]*h.m11+this->sigma[5]*h.m12+this->sigma[4]*h.m13) +
        h.m12*(this->sigma[5]*h.m11+this->sigma[1]*h.m12+this->sigma[3]*h.m13) +
        h.m13*(this->sigma[4]*h.m11+this->sigma[3]*h.m12+this->sigma[2]*h.m13);
      this->fdev[1] =
        h.m22*(              this->sigma[1]*h.m22+this->sigma[3]*h.m23) +
        h.m23*(              this->sigma[3]*h.m22+this->sigma[2]*h.m23);
      this->fdev[2] =
        h.m33*(                            this->sigma[2]*h.m33);
      this->fdev[3] =
        h.m22*(                            this->sigma[3]*h.m33) +
        h.m23*(                            this->sigma[2]*h.m33);
      this->fdev[4] =
        h.m11*(                            this->sigma[4]*h.m33) +
        h.m12*(                            this->sigma[3]*h.m33) +
        h.m13*(                            this->sigma[2]*h.m33);
      this->fdev[5] =
        h.m11*(              this->sigma[5]*h.m22+this->sigma[4]*h.m23) +
        h.m12*(              this->sigma[1]*h.m22+this->sigma[3]*h.m23) +
        h.m13*(              this->sigma[3]*h.m22+this->sigma[2]*h.m23);
      
    }

  };

}

