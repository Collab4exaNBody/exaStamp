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
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/algorithm.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/physics/units.h>
#include <exanb/core/particle_id_codec.h>

#include <memory>
#include <math.h>
#include <mpi.h>
#include <iostream>
#include <cmath>
//#include <yaml-cpp/yaml.h>
#include <functional>
#include <vector>
#include <string>

#include <onika/physics/units.h>
#include <onika/physics/units.h>
//#include "exanb/container_utils.h"
#include <onika/log.h>
#include <onika/file_utils.h>

#include "lchbop_inits.h"
//#include "parallel_build_dual_nbtbl.h"

using std::cout;
using std::cin;
using std::endl;



void Initialisation_LCHBOP(MPI_Comm,ParticleSpecies&,VLCHBOP&,NbCells&);
//void getnnbc(NbCells&);

//void printVLR(VLCHBOP&);
//void printGc(VLCHBOP&);
//void printGh(VLCHBOP&);
//void printH(VLCHBOP&);

namespace exaStamp
{
  class LCHBOP_inits : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( MPI_Comm        , mpi        , INPUT , MPI_COMM_WORLD , DocString{"MPI communicator"} );
    ADD_SLOT( bool            , NPT        , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies , species    , INPUT );
    ADD_SLOT( VLCHBOP         , parameters , OUTPUT );
    ADD_SLOT( NbCells         , nbclist    , INPUT_OUTPUT );
//    ADD_SLOT( bool            , flag_vmr   , INPUT_OUTPUT );
//    ADD_SLOT( double          , rcut_inc   , INPUT_OUTPUT );
    
    public:

    inline void execute () override final;
  };

  void LCHBOP_inits::execute ()
  {
//    bool flag_vmr = *(this->flag_vmr);
//    cout << "Initialisation LHCBOP: flag_vmr =" << std::boolalpha << flag_vmr << endl;
//    cin.get();

//    double rcut_inc= *(this->rcut_inc);
//    cout << "Initialisation LHCBOP: rcut_inc =" << rcut_inc << endl;
//    cin.get();

    MPI_Comm comm = *(this->mpi);
    ParticleSpecies& spcs = *(this->species);
    VLCHBOP& par = *(this->parameters);
    NbCells& nbc = *(this->nbclist);
//    getnnbc(nbc);

    Initialisation_LCHBOP(comm,spcs,par,nbc); 
  }

// === register factories ===  
  __attribute__((constructor)) static void register_factories()
  {
    OperatorNodeFactory::instance()->register_factory( "LCHBOP_inits",
    make_compatible_operator< LCHBOP_inits > );
  }

} // namespace

//=========================================================================//

void Initialisation_LCHBOP(MPI_Comm comm, ParticleSpecies& spcs,VLCHBOP& par,NbCells& nbc)
{
  FILE *ATPOS;
  ATPOS=fopen("TMPDAT/atpos.xyz","w");
  fclose(ATPOS);

  // MPI Initialization
  int rank=0;
  MPI_Comm_rank(comm, &rank);

  if (rank==0) lout << "======== LCHBOP ========" << endl;

  par.nspc=spcs.size();

  par.ispcc=0;
  if(spcs[0].name()=="H"){par.ispcc=1;}
  par.ispch=1-par.ispcc;

//  cout << "nspc=" << par.nspc << endl;
//  cout << "spcs[0]=" << spcs[0].m_name << endl;
//  cout << "par.ispcc=" << par.ispcc << endl;
//  cout << "par.ispch=" << par.ispch << endl;
//  abort();

  par.r2srm=0.;
  if(par.nspc==1)
  {
    if(par.ispcc==0){par.initVCC(comm);}
    else{par.initVHH(comm);}
  }  
  else  
  {
    par.initVCC(comm);
    if (rank==0) lout << "\tCC parameters are read" << endl;
    par.initVHH(comm);
    if (rank==0) lout << "\tHH parameters are read" << endl;
    par.initVCH(comm);
    if (rank==0) lout << "\tCH parameters are read" << endl;
  }
//  cout << "r2srm =" << par.r2srm << endl;

  nbc.getnbclst();

//  printVLR(par);
//  printGc(par);
//  printGh(par);
//  printH(par);
//  abort();

//  cout << "rvlr[isp][isp]=" << par.rvlr[isp][isp] << endl;
//  cout << "r1clr[isp][isp]=" << par.r1clr[isp][isp] << endl;
//  cout << "epsLJ[isp][isp]=" << par.epsLJ[isp][isp] << endl;
//  cout << "sigLJ[isp][isp]=" << par.sigLJ[isp][isp] << endl;
//  cout << "Clr[isp][isp]=" << par.Clr[isp][isp] << endl;
//  cout << "Dlr[isp][isp]=" << par.Dlr[isp][isp] << endl;
//  cout << "Elr[isp][isp]=" << par.Elr[isp][isp] << endl;

//  cout << "Initialisation finished " << endl;
  if (rank==0) lout << "========================" << endl << endl;

}

//=======================================================================//
