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

// ---------------------------------------
// Reading an xyz atomes base
// and writting Legacy MPI IO Stamp3 Dump
// ---------------------------------------

#include <exaStamp/io/StampV3LegacyIOStructures.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace exaStamp;

int main (int argc, char *argv[]) 
{
  if( argc < 3 )
    {
        std::cerr << "Usage: " << argv[0] << " <xyz-input-file> <exastamp-mpio-output-file>" << std::endl;
        return -1;
    }

    const std::string inputFile( argv[1] );
    const std::string filename( argv[2] );

  std::cout << "Converting " << inputFile << " to " << filename << " ..." << std::endl;

  LegacyHeaderIOStruct entete ;
  LegacyParticleIOStruct uneParticule;
  LegacyParticleIOStruct *particlesArray;
   
  int rank, size, count, i;
  long PaticleID;
  double x,y,z;
  std::string UneLigne;
  char unType[100],unMot[100];
  char uneTabDeMots[5][100];
  double uneTableDeDouble[20];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::ifstream intialXyzFile;
  LegacySystemIOFile mpiioDumpFile ;
  
 
  intialXyzFile.open(inputFile.c_str());
  if ( !intialXyzFile.good() ) {
    std::cerr << "Impossible de lire dans le fichier !" << std::endl;
    exit(1);
  }


  getline(intialXyzFile,UneLigne);
  sscanf(UneLigne.c_str(),"%d",&count);
  std::cout << "Nombre de particules a traiter : " << count << std::endl;
  //std::cout << "--------------------" << std::endl;

  particlesArray = new LegacyParticleIOStruct[count];

  entete.iterationNumber = 0 ;
  entete.particlesTotalNumber = count;
  entete.time = 0 ;
  entete.domainNumber=1;

  getline(intialXyzFile,UneLigne);
  sscanf(UneLigne.c_str(),"%s energy=%lf virial=\"%lf %lf %lf %lf %lf %lf\" Lattice=\"%lf %lf %lf %lf %lf %lf %lf %lf %lf\" ",uneTabDeMots[0],&entete.potentialEnergy, &uneTableDeDouble[0], &uneTableDeDouble[1], &uneTableDeDouble[2], &uneTableDeDouble[3],&uneTableDeDouble[4],&uneTableDeDouble[5],&entete.xmax,&uneTableDeDouble[6],&uneTableDeDouble[7],&uneTableDeDouble[8], &entete.ymax,&uneTableDeDouble[9],&uneTableDeDouble[10],&uneTableDeDouble[11], &entete.zmax);

  std::cout << "max = " << entete.xmax << " , " << entete.ymax << " , " << entete.zmax << std::endl;

  // std::cout<< " Energie potentiel : " << entete.potentialEnergy << std::endl; // en A
  entete.totalEnergy=entete.potentialEnergy;
  entete.internalEnergy=entete.potentialEnergy;
  entete.kineticEnergy=0.0;
  entete.xmax=entete.xmax * 1.e-10; // transforme en metre
  entete.ymax=entete.ymax * 1.e-10; 
  entete.zmax=entete.zmax * 1.e-10; 

  entete.xmin=0.0;
  entete.ymin=0.0;
  entete.zmin=0.0;
  
  /*
  std::cout<< " xmax : " << entete.xmax << std::endl;
  std::cout<< " ymax : " << entete.ymax << std::endl;
  std::cout<< " zmax : " << entete.zmax << std::endl;
  std::cout << "--------------------" << std::endl << std::endl;
  */

  PaticleID=0;
  getline(intialXyzFile,UneLigne);
  while (!intialXyzFile.eof())
  {
    if( sscanf(UneLigne.c_str(),"%s %lf %lf %lf", unType, &x , &y, &z) == 4)
	  {
	    // std::cout<< "'" << unType << "' '" << x << "' '" << y << "' '" << z << "'" << std::endl;
	    particlesArray[PaticleID].coordinates[0] = x * 1.e-10;
	    particlesArray[PaticleID].coordinates[1] = y * 1.e-10;
	    particlesArray[PaticleID].coordinates[2] = z * 1.e-10;
	    
	    // Pour le moment un seul et meme type
	    particlesArray[PaticleID].particleType=0;
	    
	    // Pas de Vitesse
	    particlesArray[PaticleID].velocity[0]=0;
	    particlesArray[PaticleID].velocity[1]=0;
	    particlesArray[PaticleID].velocity[2]=0;

	    // Son unique ID
	    particlesArray[PaticleID].particleID= PaticleID+1;

	    // On passe a la suivante
	    PaticleID++;	  
	    getline(intialXyzFile,UneLigne);
	  }
    else
	  {
	    // erreur
	    std::cerr << "Erreur de lecture  : "<< UneLigne << std::endl;
	  }
  }

  intialXyzFile.close();

  std::cout << "Nombre de particules lues : " << PaticleID << std::endl;
  if (PaticleID != count)
    {
    std::cerr << "Pas de correspondance entre le nombre de particules lues et la taille indique dans le fichier !" << std::endl;
    exit(1);      
    }
  // std::cout << "--------------------" << std::endl << std::endl;

  std::cout << "Ecriture du fichier Dump pour ExaStamp " << std::endl;
  mpiioDumpFile.open(filename.c_str(), "w");
  mpiioDumpFile.writeHeader(entete);
  mpiioDumpFile.writeArrayOfParticles(particlesArray, count);
  mpiioDumpFile.close();
 
  MPI_Finalize();

  return 0;
}
