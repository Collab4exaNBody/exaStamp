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

#include <mpi.h>
#include <string>
#include <map>
#include <cstdlib>
#include <cstdint>

// ExaStamp include for array class
//#include </cea/dsku/u-pelote/hal1/home/pmc/colombet/STAMP/exaSTAMP_project/exaSTAMP/trunk/include/utils/array/array.h>

namespace exaStamp
{

  static constexpr union{ const char idstr[4]; int32_t x; } XsV2ExtensionMarker = { { 'x','s','v','2' } };

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

    bool operator == (const LegacyHeaderIOStruct& h) const ;
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



  // All needed to read or write a legacy Stamp3 MPI IO Dump
  // Some stuff are ugly (due to the legacy) sorry about that
  class LegacySystemIOFile {

  private :

    static constexpr int MAX_MPI_WRITE_RETRY = 10;

    class PotentialCoupleData {
      
    public :
      std::string couple[2];
      std::string name;
      double epsilon, sigma, rcut,c, A0,n,m;
    };

    std::string legacyDumpFileName;

    // FAtome file data
    int NumberOfParticleType;
    std::string* particleTypeArray;
    double* particleMassArray;

    int nbPotentialCouples ;
    int nbLJ, nbSUTTOCHEN, nbIG, nbLCHBOP; //Potentials "available" in Stamp3
    PotentialCoupleData* potentialCoupleDataArray;

    // DONNEE file data
    std::string boundaryConditions[3] ; // ConditionX, ConditionY, ConditionZ 
    double deltaT ;                   // Deltatemps 
    int outputRate;                   // ParaViewFrequence
    
    std::map<std::string,std::string> stamp3ToExaStamp ;
    std::map<std::string,int> exaStampPotentialForSwitch ;

    // For MPI IO
    MPI_Comm     m_comm = MPI_COMM_WORLD;
    MPI_Offset   currentOffset, HeaderOffset,ParticlesOffset;
    MPI_File     mpiioFileName;
    MPI_Status   status;
    MPI_Datatype HeaderMPIType, ParticlesMPIType;

  public :
    LegacySystemIOFile() ;
    ~LegacySystemIOFile() ;

    void setMpiComm(MPI_Comm comm);

    // open StampV3 MPI IO Dump with read only status if rw = "r" ans write and create status if not
    void open(const char* filename, std::string rw ="r");	
    void close() ;

    // increment offset to skip header data
    void skipHeader();

    // read the data header of the dump
    void readHeader(LegacyHeaderIOStruct &entete);

    // read the one particle
    void readParticle(LegacyParticleIOStruct &uneParticule);

    // read un arry of particles of size count (must be allocated)
    void readArrayOfParticles(LegacyParticleIOStruct *particlesArray, size_t count);

    // read a raw block
    void readExtendedBlock(void *buffer, size_t count);

    // write one particle
    void writeHeader(LegacyHeaderIOStruct &entete);

    // write an array of particles
    void writeArrayOfParticles(LegacyParticleIOStruct *particlesArray, size_t count);

    // write a raw block
    void writeExtendedBlock(void *buffer, size_t count);

    // Read the DONNEE and FAtome files to know the simulation data stored in the MPIIO Dump
    void setSimulationData();

    void incrementCurrentOffset(size_t count);
    
    int getNumberOfAtomicType();

    std::string getAtomicName(int type);

    double getAtomicMass(int type);

    int getNumberOfLJPotential();

    int getNumberOfIGPotential();

    int getNumberOfSUTTOCHENPotential();

    int getNumberOfLCHBOPPotential();

  };

}

