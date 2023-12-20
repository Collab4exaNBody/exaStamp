/// @file 
/// @brief

#ifndef __STAMP_LEGACY_IO_STRUCTURES_HPP_INCLUDED
#define __STAMP_LEGACY_IO_STRUCTURES_HPP_INCLUDED

#include <iostream>    
#include <fstream>
#include <sstream>
#include <map>
#include <stdlib.h> 
#include <string.h>

//#include "parallel/mympi.hpp"
#include <mpi.h>

// ExaStamp include for array class
//#include </cea/dsku/u-pelote/hal1/home/pmc/colombet/STAMP/exaSTAMP_project/exaSTAMP/trunk/include/utils/array/array.hpp>

// Header of the legacy Stamp V3 MPI-IO format
class LegacyHeaderIOStruct {
public:

  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
  double zmin = 0.0;
  double zmax = 0.0;
  double time = 0.0;
  double CPUtime = 0.0;
  double totalEnergy = 0.0;
  double potentialEnergy = 0.0;
  double internalEnergy = 0.0;
  double kineticEnergy = 0.0;
  double rotationalEnergy =0.0;
  int32_t particlesTotalNumber = 0;
  int32_t iterationNumber = 0;
  int32_t domainNumber = 0; 
  double johnDoe[5] = { 0., 0., 0., 0., 0. };

};
/* definition de la structure enregistrement des donnees moleculaires */
class EnregistrementExaStamp_version {
public:
	int	version;
};


// Header of the legacy Stamp V4.1 MPI-IO format
class LegacyHeaderIOStructV4_2 {
public:

	/* donnees structurales */
	long	Natomes;
	double	long_a;
	double	long_b;
	double	long_c;
	double	angle_a;
	double	angle_b;
	double	angle_g;
	double	MatriceCR[3][3];
	double	XCGeo,YCGeo,ZCGeo;

	/* donnees temporelles */
	long	NumeroIterationAbsolu;
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
	long	bloc_molecules;
	long	bloc_dpd;
	long	bloc_graines;
	long	bloc_molrig;
	long	bloc_monomeres;
	long	bloc_polymerisation;
	long	bloc_posfiltre;
	long	bloc_posinit;
	long	bloc_numcelel;
	long	bloc_eeq;

	/* espaces disponibles pour des donnees supplementaires, en phase de test, avant nouveau versionnage par exemple */
	long	bloc_Ijohndoe01;
	long	bloc_Ijohndoe02;
	long	bloc_Ijohndoe03;
	long	bloc_Ijohndoe04;
	long	bloc_Ijohndoe05;
	long	bloc_Ijohndoe06;
	long	bloc_Ijohndoe07;
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


using LegacyHeaderIOStructV4_1 = LegacyHeaderIOStructV4_2;

// Header of the legacy Stamp V4.1 MPI-IO format
class LegacyVersionIOStructV4_1 {
public:
	int Version;
};

class LegacyVersionIOStructV4_2 {
public:
	int Version;
};


// Particle Elementary type of the legacy Stamp V3 MPI-IO format
class LegacyParticleIOStruct {
public:

  double coordinates[3];
  double johnDoe[3];
  double velocity[3];
  double quaternion[4];
  double momentumAngular[3];
  double orientation[3];
  int32_t particleType;
  int32_t particleID;

};

// Particle Elementary type of the legacy Stamp V4.1 MPI-IO format
class LegacyParticleIOStructV4_2 {
public:

	char	particleType[16];
	long	particleID;
	double	charge;
	double	coordinates[3];
	double	velocity[3];
};

using LegacyParticleIOStructV4_1 = LegacyParticleIOStructV4_2;

struct LegacyParticleIOStructMolRigV4_2 {

	double	quaternion[4];
	double	momentangulaire[3];
};

using LegacyParticleIOStructMolRigV4_1 = LegacyParticleIOStructMolRigV4_2;

// All needed to read or write a legacy Stamp3 MPI IO Dump
// Some stuff are ugly (due to the legacy) sorry about that
class LegacySystemIOFile {


private :

  class PotentialCoupleData {
    
  public :
    std::string couple[2];
    std::string name;
    double epsilon, sigma, rcut,c, A0,n,m;
  };

  char *legacyDumpFileName;

  // FAtome file data
  int NumberOfParticleType;
  std::string* particleTypeArray;
  double* particleMassArray;

  int nbPotentialCouples ;
  int nbLJ, nbSUTTOCHEN, nbIG, nbLCHBOP; //Potentials "available" in Stamp3
  PotentialCoupleData* potentialCoupleDataArray;

  // DONNEE file data
  std::string* boundaryConditions ; // ConditionX, ConditionY, ConditionZ 
  double deltaT ;                   // Deltatemps 
  int outputRate;                   // ParaViewFrequence
  
  std::map<std::string,std::string> stamp3ToExaStamp ;
  std::map<std::string,int> exaStampPotentialForSwitch ;

  // For MPI IO
  MPI_Offset   currentOffset,HeaderOffset,HeaderOffsetV4_1,HeaderOffsetV4_2,ParticlesOffset,ParticlesOffsetV4_1,ParticlesOffsetV4_2,ParticlesOffsetMolRigV4_1,ParticlesOffsetMolRigV4_2,VersionOffsetV4_1,VersionOffsetV4_2;

  MPI_File     mpiioFileName;
  MPI_Status   status;
  MPI_Datatype HeaderMPIType, ParticlesMPIType;
  MPI_Datatype HeaderMPITypeV4_1,HeaderMPITypeV4_2, VersionNumberTypeV4_1, VersionNumberTypeV4_2,ParticlesMPITypeV4_1,ParticlesMPITypeV4_2,ParticlesMPITypeMolRigV4_1,ParticlesMPITypeMolRigV4_2;
  
public :
  LegacySystemIOFile() ;
  ~LegacySystemIOFile() ;

  // open StampV3 MPI IO Dump with read only status if rw = "r" ans write and create status if not
  void open(const char* filename, std::string rw ="r");
  void close() ;

  // read the data header of the dump
  void readHeader(LegacyHeaderIOStruct &entete);
  void readHeaderV4_1(LegacyHeaderIOStructV4_1 &entete);
  void readHeaderV4_2(LegacyHeaderIOStructV4_2 &entete);

  // read the one particle
  void readParticle(LegacyParticleIOStruct &uneParticule);

  // read un arry of particles of size count (must be allocated)
  void readArrayOfParticles(LegacyParticleIOStruct *particlesArray, int count);
  void readArrayOfParticlesV4_1(LegacyParticleIOStructV4_1 *particlesArrayV4_1, int count);
  void readArrayOfParticlesV4_2(LegacyParticleIOStructV4_2 *particlesArrayV4_2, int count);

  // write one particle
  void writeHeader(LegacyHeaderIOStruct &entete);
  void writeHeaderV4_1(LegacyHeaderIOStructV4_1 &entete);
  void writeHeaderV4_2(LegacyHeaderIOStructV4_2 &entete);

  // write version number
  void writeVersionNumber(LegacyVersionIOStructV4_1 &versionNumber);
  void writeVersionNumber(LegacyVersionIOStructV4_2 &versionNumber);

  // read version number
  void readVersionNumber(LegacyVersionIOStructV4_1 &versionNumber);
  void readVersionNumber(LegacyVersionIOStructV4_2 &versionNumber);

  // write an array of particles
  void writeArrayOfParticles(LegacyParticleIOStruct *particlesArray, int count);
  void writeArrayOfParticlesV4_1(LegacyParticleIOStructV4_1 *particlesArrayV4_1, int count);
  void writeArrayOfParticlesV4_2(LegacyParticleIOStructV4_2 *particlesArrayV4_2, int count);

  void writeArrayOfParticlesMolRigV4_1(LegacyParticleIOStructMolRigV4_1 *particlesArrayVMolRig4_1, int count);
  void writeArrayOfParticlesMolRigV4_2(LegacyParticleIOStructMolRigV4_2 *particlesArrayVMolRig4_2, int count);
 
  void readArrayOfParticlesMolRigV4_1(LegacyParticleIOStructMolRigV4_1 *particlesArrayVMolRig4_1, int count);

  // Read the DONNEE and FAtome files to know the simulation data stored in the MPIIO Dump
  void setSimulationData();
  
  int getNumberOfAtomicType();

  std::string getAtomicName(int type);

  double getAtomicMass(int type);

  int getNumberOfLJPotential();

  int getNumberOfIGPotential();

  int getNumberOfSUTTOCHENPotential();

  int getNumberOfLCHBOPPotential();

};



#endif // __PARTICLE_BUFFER_HPP_INCLUDED
