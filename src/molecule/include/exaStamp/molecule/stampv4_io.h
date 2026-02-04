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

#include <cstdint>

namespace exaStamp
{

  namespace stampv4
  {
    static inline constexpr unsigned int NMAXCONNECT = 5;
    static inline constexpr unsigned int NMAXLIAISON = 8000;

	  /* definition de la structure enregistrement des donnees moleculaires */
	  struct IOVersion
	  {
	    int32_t	version;
	  };

	  /* definition de la structure entete de l'enregistrement pour MPIIO */
	  template<class IndexT = int32_t> 
	  struct IOEntete
	  {
	    /* donnees structurales */
	    IndexT	Natomes;
	    double	long_a;
	    double	long_b;
	    double	long_c;
	    double	angle_a;
	    double	angle_b;
	    double	angle_g;
	    double	MatriceCR[3][3];
	    double	XCGeo,YCGeo,ZCGeo;

	    /* donnees temporelles */
	    IndexT	NumeroIterationAbsolu;
	    double	tempsPhysique;

	    /* donnees energetiques */
	    double	EnergieTotale;
	    double	EnergiePotentielle;
	    double	EnergieCinetique;
	    double	EnergieRotationnelle;

	    /* invariants */
	    double	invariant;

	    /* variables dynamiques de simulation */
	    double	dt_adaptatif;
	    double	LNVhug_T;
	    double	LNVhug_Tref;
	    double	NVT_gamma[3];
	    double	NVT_gammap[3];
	    double	NPT_qsi[3];
	    double	NPT_qsip[3];
	    double	NPH_omega[3];
	    double	NPH_omegap[3];
	    double	NPH_pi[3];

	    /* blocs presents dans la protection */
	    IndexT	bloc_molecules;
	    IndexT	bloc_dpd;
	    IndexT	bloc_graines;
	    IndexT	bloc_molrig;
	    IndexT	bloc_monomeres;
	    IndexT	bloc_polymerisation;
	    IndexT	bloc_posfiltre;
	    IndexT	bloc_posinit;
	    IndexT	bloc_numcelel;
	    IndexT	bloc_eeq;

	    /* espaces disponibles pour des donnees supplementaires, en phase de test, avant nouveau versionnage par exemple */
	    IndexT	bloc_Ijohndoe01;
	    IndexT	bloc_Ijohndoe02;
	    IndexT	bloc_Ijohndoe03;
	    IndexT	bloc_Ijohndoe04;
	    IndexT	bloc_Ijohndoe05;
	    IndexT	bloc_Ijohndoe06;
	    IndexT	bloc_Ijohndoe07;
	    double	bloc_Rjohndoe01;
	    double	bloc_Rjohndoe02;
	    double	bloc_Rjohndoe03;
	    double	bloc_Rjohndoe04;
	    double	bloc_Rjohndoe05;
	    double	bloc_Rjohndoe06;
	    double	bloc_Rjohndoe07;
	    double	bloc_Rjohndoe08;
	    double	bloc_Rjohndoe09;
	    double	bloc_Rjohndoe10;
	    
	    inline bool has_xsv2_extension() const { return bloc_Ijohndoe07 == 1002; }
	    inline void set_xsv2_extension() { bloc_Ijohndoe07 = 1002; }
	  };


	  /* definition de la structure enregistrement des donnees atomiques */
	  template<class IndexT=int32_t>
	  struct IOAtomes
	  {
	    char	typeAtome[16];
	    IndexT	numeroAtome;
	    double	charge;
	    double	Position[3];
	    double	Vitesse[3];
	  };


	  /* definition de la structure enregistrement des donnees moleculaires */
	  struct IOMolecules
	  {
	    char	typeAtomeFF[16];
	    int	typeMolecule;
	    int	numeroMolecule;
	    int	Connectivite[NMAXCONNECT];
	  };


	  /* definition de la structure enregistrement des graines - nombres aleatoires */
	  struct IOGraines
	  {
    	long	graine;
	  };


	  /* definition de la structure enregistrement des molecules rigides */
	  struct IOMolRigV4_1
	  {
	    double	quaternion[4];
	    double	momentangulaire[3];
	    double	orientation[3];
	  };

	  struct IOMolRigV4_2
	  {
	    double	quaternion[4];
	    double	momentangulaire[3];
	  };

	  /* definition de la structure enregistrement des donnees monomeres */
	  struct IOMonomeres
	  {
	    int	numeroMonomere;
	    int	reactivite;
	  };


	  /* definition de la structure enregistrement des donnees de polymerisation */
	  struct IOPolymerisation
	  {
	    int	fin;
	    double	temps;
	    int	iter;
	    int	phase;
	    int	maxLiaison;
	    int	nouvLiaisons[NMAXLIAISON];
	  };


	  /* definition de la structure enregistrement des positions/initiales/filtrees */
	  struct IOPosFiltre
	  {
	    double	PositionInitiale[3];
	    double	PositionFiltree[3];
	  };

	  /* definition de la structure enregistrement DPD */
	  struct IODpd
	  {
	    double	E_interne;
	    double	T_interne;
	    double	Reaction;
	    double	Reaction2;
	  };

    /* definition de la structure enregistrement numcelel */
    struct IONumcelel {

    int numeroCelEl;
    };

    /* definition de la structure enregistrement Eeq */
    struct IOEeq {

    double eeq;
    };

  }
  
}


