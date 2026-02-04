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


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
 
#include <cstdint>

#include <exaStamp/io/StampV3LegacyIOStructures.h>

using namespace exaStamp;

template<typename StreamT>
static inline void print_lammps_atoms(StreamT &flux, uint64_t *id, double *x, double *y, double *z, int number_of_atoms)
{
  for(int i=0; i<number_of_atoms; i++)
    flux << id[i] << " " << 1 << " " << x[i]*10 << " " << y[i]*10 << " " << z[i]*10 << std::endl;
};

template<typename StreamT>
static inline void print_lammps_vitesse(StreamT &flux, uint64_t *id, double *dx, double *dy, double *dz, int number_of_atoms)
{
  for(int i=0; i<number_of_atoms; i++)
    flux << id[i] << " " << dx[i]*10 << " " << dy[i]*10 << " " << dz[i]*10 << std::endl;
};


int main(int argc, char** argv)
{ 

    if( argc < 3 )
    {
	std::cerr << "Usage: " << argv[0] << " <xsp-input-file> <lamps-output-file>" << std::endl;
	return -1;
    }

    const std::string inputFile( argv[1] );
    const std::string filename( argv[2] );

    std::cout << "Converting " << inputFile << " to " << filename << " ..." << std::endl;

    int rank=0, sizeOfMpi=1;
    MPI_Init(&argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size(MPI_COMM_WORLD, &sizeOfMpi);


    { // open a block and close it before MPI_Finalize, so that all mpi resources are released before Fnialize

    LegacyHeaderIOStruct header ;  ///< header of a MpiIO file
    LegacySystemIOFile mpiioDumpFile ; ///< contains functions that we need to read MpiIO file
    LegacySystemIOFile mpiioDumpFile2 ; ///< contains functions that we need to read MpiIO file

    mpiioDumpFile.open(inputFile.c_str());
    mpiioDumpFile.readHeader(header);
    
    std::ofstream outputFile( filename.c_str() );
    
    // Use to sort the files in lexicographical order 
    std::string step1("a_");
    std::string step2("b_");
    std::string step3("c_");
    std::string step4("d_");
    
    //rescale
    double scale = 1E9;
    header.xmin *=scale;
    header.ymin *=scale;
    header.zmin *=scale;
    header.xmax *=scale;
    header.ymax *=scale;
    header.zmax *=scale;

    std::cout << "Iteration Number = " << header.iterationNumber << std::endl;
    std::cout << "Total Particles  = " << header.particlesTotalNumber << std::endl;
    std::cout << "Simulation Time  = " << header.time << std::endl;
    std::cout << " --------------------------------------------"<<std::endl;
    std::cout << "Tailles du domaines :"<< std::endl;
    std::cout << "Xmin x Xmax = " << header.xmin << " x " << header.xmax << std::endl;
    std::cout << "Ymin x Ymax = " << header.ymin << " x " << header.ymax << std::endl;
    std::cout << "Zmin x Zmax = " << header.zmin << " x " << header.zmax << std::endl;
    std::cout << " --------------------------------------------"<<std::endl;
    std::cout << "Energie totale         = " << header.totalEnergy << std::endl;
    std::cout << "Energie poetentiel     = " << header.potentialEnergy << std::endl;
    std::cout << "Energie interne        = " << header.internalEnergy << std::endl;
    std::cout << "Energie cinetique      = " << header.kineticEnergy << std::endl;
    std::cout << "Energie rotationnelle  = " << header.rotationalEnergy << std::endl;
    std::cout << " --------------------------------------------"<<std::endl;
    
    double *tmp_x, *tmp_y, *tmp_z, *tmp_dx, *tmp_dy, *tmp_dz;
    uint64_t *tmp_id;

    scale = 1E9;

    int N=5000;

    LegacyParticleIOStruct *particlesArray = new LegacyParticleIOStruct[N];
    tmp_x = new double[N];
    tmp_y = new double[N];
    tmp_z = new double[N];
    tmp_dx = new double[N];
    tmp_dy = new double[N];
    tmp_dz = new double[N];
    tmp_id = new uint64_t[N];

    int n = header.particlesTotalNumber;

    std::string header_str, header_bis_str;
    if(rank==0)
    {
      //std::string tmp("_header");
      //const std::string filenameheader = step1 + filename + tmp;
      //ofstream fluxheader( filenameheader.c_str() );
      std::ostringstream fluxheader;
    
      //tmp = "_header_bis";
      //const std::string filenameheaderbis = step3 + filename + tmp;
      //ofstream fluxheaderbis(filenameheaderbis.c_str(), std::ofstream::out);
      std::ostringstream fluxheaderbis;

      fluxheader << "LAMMPS data file" << std::endl << std::endl;
      fluxheader << n  << "  atoms" << std::endl<<std::endl;
      fluxheader <<  header.xmin*10 << " " <<  header.xmax*10 << " xlo xhi" <<std::endl;
      fluxheader <<  header.ymin*10 << " " <<  header.ymax*10 << " ylo yhi" <<std::endl;
      fluxheader <<  header.zmin*10 << " " << header.zmax*10 << " zlo zhi" <<std::endl;
      fluxheader << std::endl;
      fluxheader << "1 atom types "<< std::endl <<std::endl;
      fluxheader << "Masses" << std::endl<<std::endl << "1 118.71 "<< std::endl<<std::endl << "Atoms # atom" << std::endl<<std::endl;
      // 118.71 --> tin      

      //fluxheader.close();
    
      fluxheaderbis<< std::endl <<"Velocities" << std::endl;
      fluxheaderbis<< std::endl;
    
      //fluxheaderbis.close();
      header_str = fluxheader.str();
      header_bis_str = fluxheaderbis.str();
    }
    
    if(rank==0) {
      std::cout << " ATOMES " <<std::endl;
    }
    
    outputFile << header_str;
    
    //std::string strc_rank = std::to_string(rank);
    //std::string atoms("_atoms_");
    //std::string velocities("_velocities_");
    
    //const std::string filenameatoms = step2 +  filename + atoms + strc_rank;
    //const std::string filenamevelocities = step4 + filename + velocities + strc_rank;  

    std::string atoms_str;
    std::ostringstream fluxatoms; //(filenameatoms.c_str(), std::ofstream::out);
          
    // compute offset
    int blockPerMPI = header.particlesTotalNumber / sizeOfMpi;   
    int begin = rank*blockPerMPI;
    int end = (rank+1)*blockPerMPI;
    
    if(rank == sizeOfMpi-1) end = header.particlesTotalNumber;
        
    mpiioDumpFile.incrementCurrentOffset (begin);    
        
    for (int i = begin; i< end; i+=N)
    {
      if(rank==0)
        std::cout << double(i)/double(end-begin)*100 << " %" << "         \r";    
        
      int size = std::min (N, end-i);
      mpiioDumpFile.readArrayOfParticles(particlesArray, size);

      for(int it = 0; it < size ; it++)
      {
        tmp_id[it] = particlesArray[it].particleID;
        tmp_x[it] = particlesArray[it].coordinates[0]*scale;
        tmp_y[it] = particlesArray[it].coordinates[1]*scale;
        tmp_z[it] = particlesArray[it].coordinates[2]*scale;        
      }
      
      print_lammps_atoms(fluxatoms,tmp_id, tmp_x, tmp_y, tmp_z, size);
    }
    atoms_str = fluxatoms.str();
    
    
    std::vector<int> tailles( sizeOfMpi , 0 );
    int bufferSize = atoms_str.length();
    MPI_Gather(&bufferSize, 1, MPI_INT, tailles.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if(rank==0)
    {
      std::cout << "Write Atoms from P0" << std::endl;
      outputFile << atoms_str;
      std::vector<char> bufferrecv;
      for(int i=1;i<sizeOfMpi;i++)
      {
        bufferrecv.resize(tailles[i]);
        std::cout << "Recv "<<tailles[i]<<" Atoms from P" <<i<< std::endl;
        MPI_Recv(bufferrecv.data(),tailles[i],MPI_CHAR,i,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::string s(bufferrecv.data(),tailles[i]);
        std::cout << "Write Atoms from P" <<i<< std::endl;
        outputFile << s;
      }
    }
    else
    {
      std::cout << "Send "<<bufferSize<<" Atoms "<< std::endl;
      MPI_Send(atoms_str.data(), bufferSize, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
    
    //fluxatoms.close();
    mpiioDumpFile.close();

  
    if(rank==0) {
      outputFile << header_bis_str;
      std::cout << " VELOCITIES " <<std::endl;
    }
  
    mpiioDumpFile2.open(inputFile.c_str());
    mpiioDumpFile2.readHeader(header);
    
    mpiioDumpFile2.incrementCurrentOffset (begin); 

    std::string velocities_str;
    {
      std::ostringstream fluxvelocities; //(filenamevelocities.c_str(), std::ofstream::out);

      for (int i = begin ;i< end; i+=N)
      {

        if(rank==0)
          std::cout << double(i)/double(end-begin)*100 << " %" << "         \r";;
      
        int size = std::min (N, end-i);
        mpiioDumpFile2.readArrayOfParticles(particlesArray, size);

        for(int it = 0; it < size ; it++)
        {
          tmp_id[it] = particlesArray[it].particleID;
          tmp_x[it] = particlesArray[it].velocity[0]/1E3;
          tmp_y[it] = particlesArray[it].velocity[1]/1E3;
          tmp_z[it] = particlesArray[it].velocity[2]/1E3;  
        }
        
        print_lammps_vitesse(fluxvelocities,tmp_id, tmp_x, tmp_y, tmp_z, size);
      }
      velocities_str = fluxvelocities.str();
    }
    //fluxvelocities.close();
    mpiioDumpFile2.close();
        
    bufferSize = velocities_str.length();
    MPI_Gather(&bufferSize, 1, MPI_INT, tailles.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);    
    if(rank==0)
    {
      std::cout << "Write Velocities from P0" << std::endl;
      outputFile << velocities_str;
      std::vector<char> bufferrecv;
      for(int i=1;i<sizeOfMpi;i++)
      {
        bufferrecv.resize(tailles[i]);
        MPI_Recv(bufferrecv.data(),tailles[i],MPI_CHAR,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        std::string s(bufferrecv.data(),tailles[i]);
        std::cout << "Write Velocities from P" <<i<< std::endl;
        outputFile << s;
      }
    }
    else
    {
      MPI_Send(velocities_str.data(), bufferSize, MPI_CHAR, 0, rank, MPI_COMM_WORLD);
    }
    outputFile.close();

    } // release all resources

    MPI_Finalize();
    
    return 0;
}
