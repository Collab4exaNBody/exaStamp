/// @file
/// @brief

#include <exaStamp/io/StampV3LegacyIOStructures.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h> 
#include <string>
#include <cstring>
#include <cassert>
#include <thread>

namespace exaStamp
{

  bool LegacyHeaderIOStruct::operator == (const LegacyHeaderIOStruct& h) const
  {
    return memcmp(this,&h,sizeof(h)) == 0;
  }

  LegacySystemIOFile::LegacySystemIOFile(){

    currentOffset = 0 ;
    particleTypeArray = NULL;
    particleMassArray = NULL;
    
    HeaderOffset=sizeof(LegacyHeaderIOStruct);
    ParticlesOffset=sizeof(LegacyParticleIOStruct);

    //boundaryConditions = new std::string[3];

    stamp3ToExaStamp["Libre"] = "free"; 
    stamp3ToExaStamp["Periodique"] = "periodic"; 
    stamp3ToExaStamp["Mur"] = "fixed"; 

    exaStampPotentialForSwitch["LCHBOP"]=0;
    exaStampPotentialForSwitch["LJ"]=1;
    exaStampPotentialForSwitch["SuttonChen"]=2;
    exaStampPotentialForSwitch["IG"]=3;

    nbLJ=0;
    nbSUTTOCHEN=0;
    nbIG=0;
    nbLCHBOP=0;

    outputRate=0;


    //  Legacy MPI TYPE CONTINOUS using MPI_CHARACTER !!!
    MPI_Type_contiguous(sizeof(LegacyHeaderIOStruct), MPI_CHARACTER, &HeaderMPIType);
    MPI_Type_commit(&HeaderMPIType);

    MPI_Type_contiguous(sizeof(LegacyParticleIOStruct), MPI_CHARACTER, &ParticlesMPIType);
    MPI_Type_commit(&ParticlesMPIType);
    
  };



  LegacySystemIOFile:: ~LegacySystemIOFile()
  {
    close();

    if (particleTypeArray != nullptr)
    {
      delete [] particleTypeArray ;
      delete [] particleMassArray ;
    }

    MPI_Type_free(&HeaderMPIType);
    MPI_Type_free(&ParticlesMPIType);
  };

  void LegacySystemIOFile::setMpiComm(MPI_Comm comm)
  {
	  m_comm = comm;
  }

  void LegacySystemIOFile::open(const char* filename, std::string rw )
  {
    int rc = MPI_SUCCESS;
    legacyDumpFileName = filename;
    assert( * ( legacyDumpFileName.data() + legacyDumpFileName.size() ) == '\0' );
    
    if (rw == "r")
    {
      rc = MPI_File_open(m_comm,legacyDumpFileName.data(), MPI_MODE_RDONLY, MPI_INFO_NULL,&mpiioFileName);
    }
    else
    {
      rc = MPI_File_open(m_comm,legacyDumpFileName.data(),MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL,&mpiioFileName);
    }
    if( rc != MPI_SUCCESS )
    {
      std::cerr << "StampV3IO: unable to open file '" << filename << "'" << std::endl << std::flush ;
      std::abort();
    }
  }


    
  void LegacySystemIOFile::close()
  {
    if( ! legacyDumpFileName.empty() )
    {
      MPI_File_close(&mpiioFileName);
      legacyDumpFileName.clear();
    }
  }



  void LegacySystemIOFile::readHeader(LegacyHeaderIOStruct &entete)
  {
    MPI_File_read_at(mpiioFileName,currentOffset,&entete,1,HeaderMPIType,&status);
    currentOffset += HeaderOffset;
  }



  void LegacySystemIOFile::readParticle(LegacyParticleIOStruct &uneParticule)
  {
    MPI_File_read_at(mpiioFileName,currentOffset,&uneParticule,1,ParticlesMPIType,&status);
    currentOffset += ParticlesOffset;
  }



  void LegacySystemIOFile::readArrayOfParticles(LegacyParticleIOStruct *particlesArray, size_t count)
  {
    size_t n = static_cast<size_t>(count) * static_cast<size_t>(ParticlesOffset);
    size_t before_offset = currentOffset;

    MPI_File_read_at(mpiioFileName,currentOffset,particlesArray,count,ParticlesMPIType,&status);

    currentOffset += n;
    if( (currentOffset - before_offset) != n )
    {
      std::cerr<<"Internal error: Bad integer arithmetic (readArrayOfParticles)"<<std::endl;
      std::abort();
    }
  }

  void LegacySystemIOFile::readExtendedBlock(void *buffer, size_t count)
  {
    MPI_File_read_at(mpiioFileName,currentOffset,buffer,count,MPI_BYTE,&status);
    currentOffset += count;
  }


  void LegacySystemIOFile::writeHeader(LegacyHeaderIOStruct &entete)
  {
    int write_count = 0;
    MPI_File_write_at(mpiioFileName,currentOffset,&entete,1,HeaderMPIType,&status);
    MPI_Get_count(&status, HeaderMPIType, &write_count);
    if( write_count != 1 )
    {
      std::cerr<<"FATAL: StampV3IO: MPI_File_write_at failed (in writeHeader)" << std::endl << std::flush ;
      std::abort();
    }
    currentOffset+=HeaderOffset ;
  }

  void LegacySystemIOFile::skipHeader()
  {
    currentOffset += HeaderOffset ;
  }


  void LegacySystemIOFile::writeArrayOfParticles(LegacyParticleIOStruct *particlesArray, size_t count)
  {
    //using namespace std::chrono_literals;
    int retry = 0;
    int write_count = 0;
    do
    {
      if( retry > 0 )
      {
        std::this_thread::sleep_for( std::chrono::duration<double,std::milli>(retry*1000.0) );
      }
      MPI_File_write_at(mpiioFileName,currentOffset,particlesArray,count,ParticlesMPIType,&status);
      write_count = 0;
      MPI_Get_count(&status, ParticlesMPIType, &write_count);
      ++ retry;
    }
    while( write_count != ssize_t(count) && retry<=MAX_MPI_WRITE_RETRY );
    
    if( write_count != ssize_t(count) )
    {
      std::cerr<<"FATAL: StampV3IO: MPI_File_write_at failed (" << retry << " attempts)" << std::endl << std::flush ;
      std::abort();
    }
    
    size_t n = static_cast<size_t>(count) * static_cast<size_t>(ParticlesOffset);
    size_t before_offset = currentOffset;
    currentOffset += n;
    if( (currentOffset - before_offset) != n )
    {
      std::cerr<<"Internal error: Bad integer arithmetic (writeArrayOfParticles)"<<std::endl;
      std::abort();
    }
  }

  void LegacySystemIOFile::writeExtendedBlock(void *buffer, size_t count)
  {
    //using namespace std::chrono_literals;
    int retry = 0;
    int write_count = 0;
    do
    {
      if( retry > 0 )
      {
        std::this_thread::sleep_for( std::chrono::duration<double,std::milli>(retry*1000.0) );
      }
      MPI_File_write_at(mpiioFileName,currentOffset,buffer,count,MPI_BYTE,&status);
      write_count = 0;
      MPI_Get_count(&status, MPI_BYTE, &write_count);
      ++ retry;
    }
    while( write_count != ssize_t(count) && retry<=MAX_MPI_WRITE_RETRY );
    
    if( write_count != ssize_t(count) )
    {
      std::cerr<<"FATAL: StampV3IO: MPI_File_write_at failed (" << retry << " attempts)" << std::endl << std::flush ;
      std::abort();
    }
    
    currentOffset += count;
  }



  // Read the DONNEE and FAtome files to know the Type and Mass stored in the MPIIO Dump
  void LegacySystemIOFile::setSimulationData()
  {
    bool find = false ;
    int nbInfoFounded =0 ;
    std::ifstream legacyDataFile ;
    std::string currentLine;
    std::string instruction, fatomeFileName, dataOfInstruction;

    legacyDataFile.open("DONNEES");
    
    while(nbInfoFounded != 6 && !legacyDataFile.eof()){
      std::getline(legacyDataFile,currentLine);
      std::stringstream lineTooSplit(currentLine);
      lineTooSplit >> instruction;
      lineTooSplit >> dataOfInstruction;

      if (instruction == "FichierAtomes"){
        fatomeFileName = dataOfInstruction;
        nbInfoFounded++;
      }

      if (instruction == "ConditionX"){
        boundaryConditions[0]=stamp3ToExaStamp[dataOfInstruction];
        nbInfoFounded++;
      }   

      if (instruction == "ConditionY"){
        boundaryConditions[1]=stamp3ToExaStamp[dataOfInstruction];
        nbInfoFounded++;
      }   

      if (instruction == "ConditionZ"){
        boundaryConditions[2]=stamp3ToExaStamp[dataOfInstruction];
        nbInfoFounded++;
      }

      if (instruction == "Deltatemps"){
        deltaT=atof(dataOfInstruction.c_str());
        nbInfoFounded++;
      }

      if (instruction == "ParaViewFrequence"){
        outputRate=atoi(dataOfInstruction.c_str());
        nbInfoFounded++;
      }
    }

    legacyDataFile.close();
    legacyDataFile.open(fatomeFileName.c_str());
    if (legacyDataFile.fail()){
      std::cerr << "Cannot open file : "<< fatomeFileName <<" in function setSimulationData"<< std::endl;
      MPI_Abort(m_comm,-1);
      exit(-1); 
    }

    find=false ;
    
    int iter = 0;
    int iterPotential = 0;

    while(!find) {
      std::getline(legacyDataFile,currentLine);
      std::stringstream lineTooSplit(currentLine);
      lineTooSplit >> instruction;
      lineTooSplit >> dataOfInstruction;
      
      if (instruction == "nbAtomes"){
        NumberOfParticleType = atoi(dataOfInstruction.c_str());
        particleTypeArray = new std::string[NumberOfParticleType];
        particleMassArray = new double[NumberOfParticleType];

        nbPotentialCouples = (NumberOfParticleType *(NumberOfParticleType + 1))/2;
        potentialCoupleDataArray = new PotentialCoupleData[nbPotentialCouples];
      }
      
      if (NumberOfParticleType != 0 && instruction == "nom"){
        particleTypeArray[iter]=dataOfInstruction;
      }
      
      if (NumberOfParticleType != 0 && instruction == "masse"){
        particleMassArray[iter]= atof(dataOfInstruction.c_str());
        iter++;
      }
      
      if (NumberOfParticleType != 0 && instruction == "Potentiel"){
        potentialCoupleDataArray[iterPotential].couple[0]=particleTypeArray[atoi(dataOfInstruction.c_str())];
        lineTooSplit >> dataOfInstruction;
        potentialCoupleDataArray[iterPotential].couple[1]=particleTypeArray[atoi(dataOfInstruction.c_str())];
        
        // A potential is found
        lineTooSplit >> dataOfInstruction;
        switch (exaStampPotentialForSwitch[dataOfInstruction]){

        case 0 : // LCHBOP
	  std::cerr << "potential "<<dataOfInstruction<<" not implemented yet !"<< std::endl;
	  nbLCHBOP++;
	  break;
	  
        case 1 : // LJ
	  potentialCoupleDataArray[iterPotential].name="LJ";

	  for (int i=0;i<3;i++) {
	    lineTooSplit >> dataOfInstruction;
	    if (dataOfInstruction =="epsilon"){
	      lineTooSplit >> dataOfInstruction;
	      potentialCoupleDataArray[iterPotential].epsilon= atof(dataOfInstruction.c_str());
	    }
	    if (dataOfInstruction =="rc"){
	      lineTooSplit >> dataOfInstruction;
	      potentialCoupleDataArray[iterPotential].rcut= atof(dataOfInstruction.c_str());
	    }	
	    if (dataOfInstruction =="sigma"){
	      lineTooSplit >> dataOfInstruction;
	      potentialCoupleDataArray[iterPotential].sigma= atof(dataOfInstruction.c_str());
	    }
	  }
	  nbLJ++;
	  break;
	  
        case 2 : // SUTTOCHEN
	  {
	    potentialCoupleDataArray[iterPotential].name="SUTTOCHEN";
	    std::ifstream suttonChenDataFileName;
	    std::string lineOfFile, extracted, toBeAValue;
	    lineTooSplit >> dataOfInstruction;
	    
	    suttonChenDataFileName.open(dataOfInstruction.c_str());
	    if (suttonChenDataFileName.fail()){
	      std::cerr << "Cannot open file : "<< dataOfInstruction <<" in function setSimulationData"<< std::endl;
	      MPI_Abort(m_comm,-1);
	      exit(-1); 
	    }
	    std::getline(suttonChenDataFileName,lineOfFile);
	    
	    while (!suttonChenDataFileName.eof()) {
	      std::stringstream localLineTooSplit(lineOfFile);	  
	      localLineTooSplit >> extracted;
	      localLineTooSplit >> toBeAValue;
	      
	      if(extracted == "epsilonSuttonChen")
	        potentialCoupleDataArray[iterPotential].epsilon=atof(toBeAValue.c_str());
	      
	      if(extracted == "rc")
	        potentialCoupleDataArray[iterPotential].rcut=atof(toBeAValue.c_str());
	      
	      if(extracted == "a0SuttonChen")
	        potentialCoupleDataArray[iterPotential].A0=atof(toBeAValue.c_str());
	      
	      if(extracted == "cSuttonChen")
	        potentialCoupleDataArray[iterPotential].c=atof(toBeAValue.c_str());
	      
	      if(extracted == "nSuttonChen")
	        potentialCoupleDataArray[iterPotential].n=atof(toBeAValue.c_str());
	      
	      if(extracted == "mSuttonChen")
	        potentialCoupleDataArray[iterPotential].m=atof(toBeAValue.c_str());
	      
	      std::getline(suttonChenDataFileName,lineOfFile);
	    }
	    potentialCoupleDataArray[iterPotential].sigma=-1.0;
	    suttonChenDataFileName.close();
	    nbSUTTOCHEN++;
	    break;
	  }

        case 3 : // IG
	  potentialCoupleDataArray[iterPotential].name="IG";
	  nbLJ++;
	  break;

        default:
	  std::cerr << "potential "<<dataOfInstruction<<" not implemented yet !"<< std::endl;
	    break;
        }
        
        iterPotential++ ;
        if (iterPotential == nbPotentialCouples) find = true;      
      }
    }
    legacyDataFile.close();
    if (nbPotentialCouples != nbLJ+nbSUTTOCHEN+nbIG+nbLCHBOP){
      std::cerr << "An Atome couple have no potential ! please check de FAtome Stamp3 file"<< std::endl;
      MPI_Abort(m_comm,-1);
      exit(-1); 
    }
  }


  void LegacySystemIOFile::incrementCurrentOffset(size_t count)
  {
    size_t before_offset = currentOffset;
    size_t n = static_cast<size_t>(count) * static_cast<size_t>(ParticlesOffset);
    currentOffset += n;
    if( ( currentOffset - before_offset ) != n ) // may arise is MPI_Offset type is too short and count is to high
    {
      std::cerr << "Internal error: bad integer arithmetic (incrementCurrentOffset)" << std::endl;
      std::abort();
    }
  }


  int LegacySystemIOFile::getNumberOfAtomicType(){
    return NumberOfParticleType;
  }



  std::string LegacySystemIOFile::getAtomicName(int type){
    return  particleTypeArray[type];
  }


  double LegacySystemIOFile::getAtomicMass(int type){
    return particleMassArray[type];
  }

  int LegacySystemIOFile::getNumberOfLJPotential(){
    return nbLJ;
  }

  int LegacySystemIOFile::getNumberOfIGPotential(){
    return nbIG;
  }

  int LegacySystemIOFile::getNumberOfSUTTOCHENPotential(){
    return nbSUTTOCHEN;
  }

  int LegacySystemIOFile::getNumberOfLCHBOPPotential(){
    return nbLCHBOP;
  }

}

