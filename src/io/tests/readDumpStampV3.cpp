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

#include <exaStamp/io/StampV3LegacyIOStructures.h>

#include <iostream>
#include <iomanip>

#include <set>

#define CHECK_POSITION_REPEATS 1

// Reading and writting Legacy MPI IO Stamp3 Dump
// Warning :
//      To read a stamp3 Dump, it's necesary to have the "DONNEE" file too !

using namespace exaStamp;

struct ParticleInfo
{
	double x=0.0;
        double y=0.0;
        double z=0.0;
	int64_t id=-1;
	inline bool operator < ( const ParticleInfo& p ) const
	{
		if( x < p.x ) return true;
		if( x > p.x ) return false;
		// x == p.x
		if( y < p.y ) return true;
		if( y > p.y ) return false;
		// y == p.y
		if( z < p.z ) return true;
		// z >= p.z
		return false;
	}
};

int main (int argc, char *argv[]) 
{
  if( argc < 2 )
  {
      std::cerr << "Usage: " << argv[0] << " <legacy-dump-file>" << std::endl;
      return -1;
  }

  bool header_only = false;
  if( argc >= 3 )
  {
    header_only = true;
  }
  
  const std::string inputFile( argv[1] );

  std::cout << "Reading " << inputFile << " ..." << std::endl;

  LegacyHeaderIOStruct entete ;
  LegacyParticleIOStruct *particlesArray;
   
  int rank, size,count;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  LegacySystemIOFile mpiioDumpFile ;

  mpiioDumpFile.open(inputFile.c_str());

  mpiioDumpFile.readHeader(entete);
  std::cout << "Iteration Number = " << entete.iterationNumber << std::endl;
  std::cout << "Total Particles  = "  << entete.particlesTotalNumber << std::endl;
  std::cout << "bounds           = (" << std::scientific << std::setprecision(10) << entete.xmin << ',' << entete.ymin << ',' << entete.zmin << ") - ("<< entete.xmax << ',' << entete.ymax << ',' << entete.zmax << ')' << std::endl;
  std::cout << "totalEnergy      = "  << entete.totalEnergy << std::endl;
  std::cout << "potentialEnergy  = "  << entete.potentialEnergy << std::endl;
  std::cout << "internalEnergy   = "  << entete.internalEnergy << std::endl;
  std::cout << "kineticEnergy    = "  << entete.kineticEnergy << std::endl;
  std::cout << "rotationalEnergy = "  << entete.rotationalEnergy << std::endl;
  std::cout << "Simulation Time  = "  << entete.time << std::endl;
  std::cout << " --------------------------------------------"<<std::endl << std::flush;
  
  if( header_only )
  {
    mpiioDumpFile.close();
    MPI_Finalize();
    return 0;
  }
  
  std::set<ParticleInfo> pset;

  count=entete.particlesTotalNumber;
  particlesArray = new LegacyParticleIOStruct[count];
  mpiioDumpFile.readArrayOfParticles(particlesArray, count);
  for (int i = 0;i<count; i++)
    {
      std::cout << "ID = " << particlesArray[i].particleID << std::endl;
      std::cout << "Type = " << particlesArray[i].particleType << std::endl;
      std::cout << "(x,y,z) = (" <<  particlesArray[i].coordinates[0] << " , " <<  particlesArray[i].coordinates[1] << " , " <<  particlesArray[i].coordinates[2] << ")"<<std::endl;
      std::cout << " --------------------------------------------"<<std::endl;
#     ifdef CHECK_POSITION_REPEATS
      ParticleInfo p = { particlesArray[i].coordinates[0] , particlesArray[i].coordinates[1], particlesArray[i].coordinates[2] , particlesArray[i].particleID };
      auto it = pset.find( p );
      if( it != pset.end() )
      {
	std::cout << "Particle #"<<it->id <<" (i="<<i<< ") is the same as "<<p.id<<" : coord="<<p.x<<","<<p.y<<","<<p.z<<  std::endl;
      }
#     endif
    } 
  delete[]  particlesArray;

  mpiioDumpFile.close();
  MPI_Finalize();

  return 0;
}
