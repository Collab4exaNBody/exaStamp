//#include "StampV3LegacyIOStructures.h"
#include "saru/saru.h"
#include <cassert>
#include <iostream>
#include <string>

#include <exaStamp/io/StampV3LegacyIOStructures.h>
// Checking Legacy MPIIO Dump


/* ****************************************************************

  Simple example
  To run in the directory where MPI IO Dump and DONNEE are located

**************************************************************** */


using namespace exaStamp;

struct atom
{
  double x;
  double y;
  double z;
};

int main (int argc, char *argv[]) 
{
  LegacyHeaderIOStruct entete ;
  LegacyParticleIOStruct *particlesArray;
 
  int rank, size,count;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // ELEM 1 = rayon de la sphère
  // ELEM 2 = Rayon de coupure
  // ELEM 3 = taille de la lattice 
  // ELEM 4 = température 
  // ELEM 5 = type de lattice

  double rSphere;
  double rCut ;
  double Alattice;
  double initialTemperature;
  std::string Tlattice;
  
  std::cout << " nb args =  " << argc << std::endl;

	if( argc != 6 )
	{
		std::cerr << "Usage: " << argv[0] << " : rSphere (double), rCut (double), Alattice (double), initialTemperature (?K?)n, lattice type (bcc/fcc)  " << std::endl;
		return -1;
	}
	else
	{
	
	  rSphere = atof(argv[1]);
    rCut = atof(argv[2]);
    Alattice = atof(argv[3]);
    initialTemperature =  atof(argv[4]);
    Tlattice = std::string( argv[5] );
	  
	  if(rank == 0)
	  {
	    std::cout << " Parametres :\n"
	      << "          rSphere = " << rSphere << std::endl
	      << "          rCut = " << rCut << std::endl
	      << "          Alattice = " << Alattice << std::endl
	      << "          initialTemperature = " << initialTemperature << std::endl
	      << std::endl;
	  }
	}
	



  LegacySystemIOFile mpiioDumpFile ;
  LegacySystemIOFile newDumpFile ;
  
  
  newDumpFile.open("initial.MpiIO", "w");

   double SCALEmeters = 1E-9;

  // remplissage de l'entete
   entete.iterationNumber = 0;
   entete.time = 0;
   entete.xmax = rSphere+rCut;
   entete.ymax = rSphere+rCut;
   entete.zmax = rSphere+rCut;
   entete.xmin = -entete.xmax;
   entete.ymin = -entete.ymax;
   entete.zmin = -entete.zmax;
   
  
   double TailleDomain = entete.xmax*2; //(Cube centré en 0)
   assert(TailleDomain > Alattice);
   double nombreLatticeParDim = TailleDomain / Alattice;
   
  size_t nbL;
  atom* lattice;

  if(Tlattice == "bcc")
  {
    nbL = 2;
    lattice = new atom[nbL];

    lattice[0].x = 0.25;
    lattice[0].y = 0.25;
    lattice[0].z = 0.25; 

    lattice[1].x = 0.75;
    lattice[1].y = 0.75;
    lattice[1].z = 0.75; 
  } 
  else if(Tlattice == "fcc")
  {
    nbL = 4;
    lattice = new atom[nbL];

    lattice[0].x = 0.25;
    lattice[0].y = 0.25;
    lattice[0].z = 0.25; 

    lattice[1].x = 0.25;
    lattice[1].y = 0.75;
    lattice[1].z = 0.75;   

    lattice[2].x = 0.75;
    lattice[2].y = 0.25;
    lattice[2].z = 0.75; 
    
    lattice[3].x = 0.75;
    lattice[3].y = 0.75;
    lattice[3].z = 0.25; 
  }
  else
  {
    std::cout << " La lattice spécifiée n'est pas une structure bcc ou fcc " << std::endl;
    exit(0);
  }
    
    
   // sphère centrée en 0
   auto inSphere = [rSphere] (atom& elem) -> bool 
   {
     return (elem.x*elem.x +elem.y*elem.y+elem.z*elem.z) < rSphere * rSphere;
   };
   
   int nbAtom = 0;
   
   std::vector<atom> keepAtom;
   
   keepAtom.clear();
   
   //#pragma omp parallel for schedule(3) reduction(+:nbAtom)
   for(int z= 0; z<nombreLatticeParDim; z++)
    for(int y= 0; y<nombreLatticeParDim; y++)
      for(int x= 0; x<nombreLatticeParDim; x++)
      {
        atom check;
        for(int i = 0 ; i < nbL; i++)
        {
          check.x = (double(x) +lattice[i].x) * Alattice + entete.xmin ;
          check.y = (double(y) +lattice[i].y) * Alattice + entete.ymin ;
          check.z = (double(z) +lattice[i].z) * Alattice + entete.zmin ;
          
          
          if(inSphere(check))
            keepAtom.push_back(check);
        }
      }
      
    std::cout << " il y aura " << keepAtom.size() << " atomes " << std::endl;
      
    entete.particlesTotalNumber = keepAtom.size();   
    
    entete.xmax *= SCALEmeters;
    entete.ymax *= SCALEmeters;
    entete.zmax *= SCALEmeters;
    entete.xmin *= SCALEmeters;
    entete.ymin *= SCALEmeters;
    entete.zmin *= SCALEmeters; 
  
    newDumpFile.writeHeader(entete);
  
    
    size_t sizeOfBlock = 1000000;
    
    size_t nbBlock = 0;
    
    if(keepAtom.size() % sizeOfBlock == 0)
      nbBlock = keepAtom.size() / sizeOfBlock;
    else
      nbBlock = keepAtom.size() / sizeOfBlock + 1;
    
    
    int maxSize = std::min(keepAtom.size(), sizeOfBlock);
    particlesArray = new LegacyParticleIOStruct[keepAtom.size()]; 
    
    double sigma = std::sqrt( 1.380662e-23 * initialTemperature / 197.10101614E-27); // TIN
    
    
    for(int k = 0 ; k < nbBlock ; k++)
    {
      size_t shift = k * sizeOfBlock;
      
      assert(keepAtom.size()>shift);
      size_t nbElemeInThisBlock = std::min(sizeOfBlock, keepAtom.size()-shift);
      
      #pragma omp parallel for
      for(int i = 0 ; i < nbElemeInThisBlock ; i++)
      {
         particlesArray[i].coordinates[0] = keepAtom[shift+i].x*SCALEmeters;
         particlesArray[i].coordinates[1] = keepAtom[shift+i].y*SCALEmeters; 
         particlesArray[i].coordinates[2] = keepAtom[shift+i].z*SCALEmeters;   
         
         particlesArray[i].velocity[0] = sigma * Saru::randN(shift+i, -1);// - 800;
         particlesArray[i].velocity[1] = sigma * Saru::randN(shift+i, -2);
         particlesArray[i].velocity[2] = sigma * Saru::randN(shift+i, -3);
         
         for(int j=0; j<3; j++)  
         {
           particlesArray[i].johnDoe[j]=0;
           particlesArray[i].momentumAngular[j]=0;
           particlesArray[i].orientation[j]=0;     
           particlesArray[i].quaternion[j]=0;            
         } 
          
         particlesArray[i].quaternion[3]=0;
           
         particlesArray[i].particleType = 0;
         particlesArray[i].particleID = shift+i+1; // != 0                
      }
      
      newDumpFile.writeArrayOfParticles(particlesArray, nbElemeInThisBlock);
    }
  

    
    newDumpFile.close();
    delete[]  particlesArray;

  if(rank == 0)
  {
    std::cout << "Iteration Number = " << entete.iterationNumber << std::endl;
    std::cout << "Total Particles = " << entete.particlesTotalNumber << std::endl;
    std::cout << "Simulation Time = " << entete.time << std::endl;
    std::cout << " --------------------------------------------"<<std::endl;
    std::cout << " Taille du domaine : "<<std::endl;
    std::cout << " xmin x xmax " << entete.xmin << " x " << entete.xmax << std::endl;
    std::cout << " ymin x ymax " << entete.ymin << " x " << entete.ymax << std::endl;
    std::cout << " zmin x zmax " << entete.zmin << " x " << entete.zmax << std::endl;
    std::cout << " --------------------------------------------"<< std::endl;
  }


  MPI_Finalize();
  
  return 0;
}
