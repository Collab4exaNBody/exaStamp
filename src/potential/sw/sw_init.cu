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

#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/unit_system.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/log.h>
#include <onika/file_utils.h>

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "sw_params.h"

namespace exaStamp
{

  using namespace exanb;

  class SWInit : public OperatorNode
  {
    ADD_SLOT( SWParams        , parameters , INPUT        , REQUIRED , DocString{"Stillinger-Weber parameter file and options"} );
    ADD_SLOT( double          , rcut_max   , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( ParticleSpecies , species    , INPUT        , REQUIRED );
    ADD_SLOT( SWContext       , sw_ctx     , OUTPUT );

    struct FileRecord
    {
      std::string e[3];
      double v[11];
    };

  public:

    inline std::string documentation() const override final
    {
      return R"EOF(

Reads a Stillinger-Weber parameter file and builds the per-species-triplet table used by sw_force,
then raises rcut_max to the SW cutoff. exaStamp species are matched to the file's elements by name,
every species triplet must have exactly one entry in the file.

File format : one entry per element triplet, 14 whitespace separated values (entries may span several
lines, '#' starts a comment) :
  elem1 elem2 elem3 epsilon(eV) sigma(ang) a lambda gamma costheta0 A B p q tol

Place it in init_parameters, after species.

Usage example:

init_parameters:
  - species
  - sw_init:
      parameters: { file: "Si.sw" }    # optional: threebody: false

)EOF";
    }

    inline void execute() override final
    {
      const std::string fname = onika::data_file_path( parameters->file );
      std::ifstream fin( fname );
      if( ! fin.good() ) { fatal_error() << "sw_init: cannot open parameter file '" << fname << "'" << std::endl; }

      // tokenize, stripping comments
      std::vector<std::string> tokens;
      std::string line;
      while( std::getline(fin,line) )
      {
        const auto hash = line.find('#');
        if( hash != std::string::npos ) line.erase(hash);
        std::istringstream iss(line);
        std::string tok;
        while( iss >> tok ) tokens.push_back(tok);
      }
      if( tokens.size() % 14 != 0 )
      {
        fatal_error() << "sw_init: '" << fname << "' contains " << tokens.size() << " values, not a multiple of 14" << std::endl;
      }

      std::vector<FileRecord> records( tokens.size() / 14 );
      for(size_t r=0; r<records.size(); r++)
      {
        for(int a=0;a<3;a++) records[r].e[a] = tokens[r*14+a];
        for(int a=0;a<11;a++)
        {
          try { records[r].v[a] = std::stod( tokens[r*14+3+a] ); }
          catch(...) { fatal_error() << "sw_init: bad numeric value '" << tokens[r*14+3+a] << "' in entry " << r << " of '" << fname << "'" << std::endl; }
        }
        for(int a=0;a<11;a++)
        {
          if( a != 5 && records[r].v[a] < 0.0 ) // all parameters but costheta0 must be non-negative
          {
            fatal_error() << "sw_init: illegal negative parameter in entry "
                          << records[r].e[0] << " " << records[r].e[1] << " " << records[r].e[2] << std::endl;
          }
        }
      }

      const auto & sp = *species;
      const int nt = sp.size();
      static constexpr double conv_energy = EXASTAMP_CONST_QUANTITY( 1. * eV );

      sw_ctx->ntypes = nt;
      sw_ctx->threebody = parameters->threebody;
      sw_ctx->triplets.assign( nt*nt*nt , SWTriplet{} );
      sw_ctx->rcut = 0.0;

      for(int i=0;i<nt;i++) for(int j=0;j<nt;j++) for(int k=0;k<nt;k++)
      {
        const std::string ni = sp[i].m_name, nj = sp[j].m_name, nk = sp[k].m_name;
        int found = -1;
        for(size_t r=0; r<records.size(); r++)
        {
          if( records[r].e[0]==ni && records[r].e[1]==nj && records[r].e[2]==nk )
          {
            if( found >= 0 ) { fatal_error() << "sw_init: duplicate entry for " << ni << " " << nj << " " << nk << std::endl; }
            found = r;
          }
        }
        if( found < 0 ) { fatal_error() << "sw_init: missing entry for " << ni << " " << nj << " " << nk << " in '" << fname << "'" << std::endl; }

        const double * v = records[found].v;
        SWTriplet & p = sw_ctx->triplets[ (i*nt + j)*nt + k ];
        p.epsilon  = v[0] * conv_energy;
        p.sigma    = v[1];
        p.littlea  = v[2];
        p.lambda   = parameters->threebody ? v[3] : 0.0;
        p.gamma    = v[4];
        p.costheta = v[5];
        p.biga     = v[6];
        p.bigb     = v[7];
        p.powerp   = v[8];
        p.powerq   = v[9];
        p.tol      = v[10];

        // derived constants
        p.cut = p.sigma * p.littlea;
        double rtmp = p.cut;
        if( p.tol > 0.0 )
        {
          if( p.tol > 0.01 ) p.tol = 0.01;
          if( p.gamma < 1.0 ) rtmp += p.gamma * p.sigma / std::log(p.tol);
          else                rtmp += p.sigma / std::log(p.tol);
        }
        p.cutsq = rtmp * rtmp;

        p.sigma_gamma     = p.sigma * p.gamma;
        p.lambda_epsilon  = p.lambda * p.epsilon;
        p.lambda_epsilon2 = 2.0 * p.lambda * p.epsilon;
        p.c1 = p.biga * p.epsilon * p.powerp * p.bigb * std::pow(p.sigma, p.powerp);
        p.c2 = p.biga * p.epsilon * p.powerq * std::pow(p.sigma, p.powerq);
        p.c3 = p.biga * p.epsilon * p.bigb * std::pow(p.sigma, p.powerp + 1.0);
        p.c4 = p.biga * p.epsilon * std::pow(p.sigma, p.powerq + 1.0);
        p.c5 = p.biga * p.epsilon * p.bigb * std::pow(p.sigma, p.powerp);
        p.c6 = p.biga * p.epsilon * std::pow(p.sigma, p.powerq);

        sw_ctx->rcut = std::max( sw_ctx->rcut , std::sqrt(p.cutsq) );

        ldbg << "SW " << ni << " " << nj << " " << nk << " : cut=" << p.cut << " rcut=" << std::sqrt(p.cutsq) << std::endl;
      }

      *rcut_max = std::max( *rcut_max , sw_ctx->rcut );
      lout << "Stillinger-Weber : " << nt << " species, cutoff = " << sw_ctx->rcut << " ang"
           << ( sw_ctx->threebody ? "" : ", three-body term disabled" ) << std::endl;
    }

  };

  ONIKA_AUTORUN_INIT(sw_init)
  {
    OperatorNodeFactory::instance()->register_factory("sw_init", make_simple_operator<SWInit>);
  }

}
