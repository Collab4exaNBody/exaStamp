/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/
/// @file
/// @brief

#include "StampV3LegacyIOStructures.hpp"



LegacySystemIOFile::LegacySystemIOFile(){

  currentOffset = 0 ;
  legacyDumpFileName = NULL;
  particleTypeArray = NULL;
  particleMassArray = NULL;
  
  HeaderOffset=sizeof(LegacyHeaderIOStruct);
  HeaderOffsetV4_1=sizeof(LegacyHeaderIOStructV4_1);
  HeaderOffsetV4_2=sizeof(LegacyHeaderIOStructV4_2);

  VersionOffsetV4_1=sizeof(LegacyVersionIOStructV4_1);
  VersionOffsetV4_2=sizeof(LegacyVersionIOStructV4_2);

  ParticlesOffset=sizeof(LegacyParticleIOStruct);
  ParticlesOffsetV4_1=sizeof(LegacyParticleIOStructV4_1);
  ParticlesOffsetV4_2=sizeof(LegacyParticleIOStructV4_2);

  ParticlesOffsetMolRigV4_1=sizeof(LegacyParticleIOStructMolRigV4_1);
  ParticlesOffsetMolRigV4_2=sizeof(LegacyParticleIOStructMolRigV4_2);

  boundaryConditions = new std::string[3];

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

  MPI_Type_contiguous(sizeof(LegacyHeaderIOStructV4_1), MPI_CHARACTER, &HeaderMPITypeV4_1);
  MPI_Type_commit(&HeaderMPITypeV4_1);

  MPI_Type_contiguous(sizeof(LegacyHeaderIOStructV4_2), MPI_CHARACTER, &HeaderMPITypeV4_2);
  MPI_Type_commit(&HeaderMPITypeV4_2);

  MPI_Type_contiguous(sizeof(LegacyVersionIOStructV4_1), MPI_CHARACTER, &VersionNumberTypeV4_1);
  MPI_Type_commit(&VersionNumberTypeV4_1);

  MPI_Type_contiguous(sizeof(LegacyVersionIOStructV4_2), MPI_CHARACTER, &VersionNumberTypeV4_2);
  MPI_Type_commit(&VersionNumberTypeV4_2);

  MPI_Type_contiguous(sizeof(LegacyParticleIOStruct), MPI_CHARACTER, &ParticlesMPIType);
  MPI_Type_commit(&ParticlesMPIType);

  MPI_Type_contiguous(sizeof(LegacyParticleIOStructV4_1), MPI_CHARACTER, &ParticlesMPITypeV4_1);
  MPI_Type_commit(&ParticlesMPITypeV4_1); 

  MPI_Type_contiguous(sizeof(LegacyParticleIOStructV4_2), MPI_CHARACTER, &ParticlesMPITypeV4_2);
  MPI_Type_commit(&ParticlesMPITypeV4_2); 
  
  MPI_Type_contiguous(sizeof(LegacyParticleIOStructMolRigV4_1), MPI_CHARACTER, &ParticlesMPITypeMolRigV4_1);
  MPI_Type_commit(&ParticlesMPITypeMolRigV4_1); 

  MPI_Type_contiguous(sizeof(LegacyParticleIOStructMolRigV4_2), MPI_CHARACTER, &ParticlesMPITypeMolRigV4_2);
  MPI_Type_commit(&ParticlesMPITypeMolRigV4_2); 
};



LegacySystemIOFile:: ~LegacySystemIOFile(){

  if (legacyDumpFileName != NULL) delete legacyDumpFileName;

  if (particleTypeArray != NULL){
    delete[] particleTypeArray ;
    delete[] particleMassArray ;
  }

  //MPI_Type_free(&HeaderMPIType);
  //MPI_Type_free(&ParticlesMPIType);

};



void LegacySystemIOFile::open(const char* filename, std::string rw ){
  legacyDumpFileName=new char[strlen(filename)];
  strcpy(legacyDumpFileName,filename);
  currentOffset = 0 ;
  std::cout<<"Ouverture du fichier "<<legacyDumpFileName<<std::endl;
  if (rw == "r"){
    if(MPI_File_open(MPI_COMM_WORLD,legacyDumpFileName, MPI_MODE_RDONLY, MPI_INFO_NULL,&mpiioFileName) != MPI_SUCCESS){
	std::cout<<"\n##### Le protection n'existe pas #####\n"<<std::endl;
	MPI_Finalize();	
	exit(0);
    }
  		
  } else {
    if(MPI_File_open(MPI_COMM_WORLD,legacyDumpFileName,MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL,&mpiioFileName) != MPI_SUCCESS){
	std::cout<<"\n##### Probleme lors de la creation du fichier de protection #####\n"<<std::endl;
	MPI_Finalize();	
	exit(0);
    }
  }
}


  
void LegacySystemIOFile::close(){
  MPI_File_close(&mpiioFileName);
}



void LegacySystemIOFile::readHeader(LegacyHeaderIOStruct &entete){
  MPI_File_read_at(mpiioFileName,currentOffset,&entete,1,HeaderMPIType,&status);
  currentOffset+=HeaderOffset;
}





void LegacySystemIOFile::readParticle(LegacyParticleIOStruct &uneParticule){
  MPI_File_read_at(mpiioFileName,currentOffset,&uneParticule,1,ParticlesMPIType,&status);
  currentOffset+=ParticlesOffset;
}



void LegacySystemIOFile::readArrayOfParticles(LegacyParticleIOStruct *particlesArray, int count){
  //std::cout<<"On lit : "<<count<<" particule - offset :  "<<currentOffset<<std::endl; 	
  MPI_File_read_at(mpiioFileName,currentOffset,particlesArray,count,ParticlesMPIType,&status);
  //std::cout<<"Fin de la lecture"<<std::endl; 	
  currentOffset+=count*ParticlesOffset;
}

void LegacySystemIOFile::readArrayOfParticlesV4_1(LegacyParticleIOStructV4_1 *particlesArrayV4_1, int count){

  // Si currentOffset < VersionOffsetV4_1 + HeaderOffsetV4_1, cela implique que l'on n'est pas passer par la lecture du numero de version et de l'entete. 
  // On force l'offset à la bonne valeur initiale pour la lection des positions.
  if(currentOffset < VersionOffsetV4_1 + HeaderOffsetV4_1){
	MPI_Offset   currentOffsetAnt=currentOffset;
	currentOffset = VersionOffsetV4_1 + HeaderOffsetV4_1;	
	std::cout<<"**** Attention, l'offset de lecture des positions, initialement à "<<currentOffsetAnt<<", est forcé à "<<currentOffset<<" ****"<<std::endl;
  }
  //std::cout<<"\t\tOn lit : "<<count<<" particule - offset :  "<<currentOffset<<std::endl; 

  MPI_File_read_at(mpiioFileName,currentOffset,particlesArrayV4_1,count,ParticlesMPITypeV4_1,&status);
  currentOffset+=count*ParticlesOffsetV4_1;
  //std::cout<<"\t\tFin de la lecture - offet : "<<currentOffset<<std::endl; 	
}

void LegacySystemIOFile::writeHeader(LegacyHeaderIOStruct &entete){
	currentOffset = 0;
	//std::cout<<"On écrit l'entete avec un offset de  :  "<<currentOffset<<std::endl; 	
	MPI_File_write_at(mpiioFileName,currentOffset,&entete,1,HeaderMPIType,&status);

	currentOffset+=HeaderOffset ;
	//std::cout<<"offset après écriture :  "<<currentOffset<<std::endl; 	
}

void LegacySystemIOFile::writeVersionNumber(LegacyVersionIOStructV4_1 &versionNumber){
	currentOffset = 0;
	//std::cout<<"On écrit la version avec un offset de  :  "<<currentOffset<<std::endl; 	
	MPI_File_write_at(mpiioFileName,currentOffset,&versionNumber,1,VersionNumberTypeV4_1,&status);

	currentOffset = VersionOffsetV4_1 ;
	//std::cout<<"offset après écriture :  "<<currentOffset<<std::endl; 	
}

void LegacySystemIOFile::readVersionNumber(LegacyVersionIOStructV4_1 &versionNumber){
	currentOffset = 0;
	//std::cout<<"On lit la version avec un offset de  :  "<<currentOffset<<std::endl; 	
	MPI_File_read_at(mpiioFileName,currentOffset,&versionNumber,1,VersionNumberTypeV4_1,&status);

	currentOffset = VersionOffsetV4_1 ;
	//std::cout<<"offset après lecture :  "<<currentOffset<<std::endl; 	
}

void LegacySystemIOFile::readHeaderV4_1(LegacyHeaderIOStructV4_1 &entete){
	currentOffset = VersionOffsetV4_1 ;
	//std::cout<<"On lit l'entete avec un offset de  :  "<<currentOffset<<std::endl; 	
	MPI_File_read_at(mpiioFileName,currentOffset,&entete,1,HeaderMPITypeV4_1,&status);
	currentOffset+=HeaderOffsetV4_1;
	//std::cout<<"offset après lecture :  "<<currentOffset<<std::endl; 	
}


void LegacySystemIOFile::writeHeaderV4_1(LegacyHeaderIOStructV4_1 &entete){
	currentOffset = VersionOffsetV4_1 ;
	//std::cout<<"On écrit l'entete avec un offset de  :  "<<currentOffset<<std::endl; 	
	MPI_File_write_at(mpiioFileName,currentOffset,&entete,1,HeaderMPITypeV4_1,&status);

	currentOffset+=HeaderOffsetV4_1 ;
	//std::cout<<"offset après écriture :  "<<currentOffset<<std::endl; 	
}

void LegacySystemIOFile::readHeaderV4_2(LegacyHeaderIOStructV4_2 &entete){
	//std::cout<<"On lit l'entete avec un offset de  :  "<<currentOffset<<std::endl; 	
	MPI_File_read_at(mpiioFileName,currentOffset,&entete,1,HeaderMPITypeV4_2,&status);
	currentOffset+=HeaderOffsetV4_2;
	//std::cout<<"offset après lecture :  "<<currentOffset<<std::endl; 	
}

void LegacySystemIOFile::writeHeaderV4_2(LegacyHeaderIOStructV4_2 &entete){
	//std::cout<<"On écrit l'entete avec un offset de  :  "<<currentOffset<<std::endl; 	
	MPI_File_write_at(mpiioFileName,currentOffset,&entete,1,HeaderMPITypeV4_2,&status);

	currentOffset+=HeaderOffsetV4_2 ;
	//std::cout<<"offset après écriture :  "<<currentOffset<<std::endl; 	
}

void LegacySystemIOFile::writeArrayOfParticles(LegacyParticleIOStruct *particlesArray, int count){
  if(currentOffset == 0)
	currentOffset = HeaderOffset;
	//std::cout<<"On écrit "<<count<<" éléments d'un tableau de type LegacyParticleIOStruct avec un offset de  :  "<<currentOffset<<std::endl; 	
  
	MPI_File_write_at(mpiioFileName,currentOffset,particlesArray,count,ParticlesMPIType,&status);
  
	currentOffset+=count*ParticlesOffset; 
	//std::cout<<"offset après écriture :  "<<currentOffset<<std::endl; 	

}



void LegacySystemIOFile::writeArrayOfParticlesV4_1(LegacyParticleIOStructV4_1 *particlesArrayV4_1, int count){

 	 if(currentOffset == 0)
		currentOffset = VersionOffsetV4_1 + HeaderOffsetV4_1;

	//std::cout<<"On écrit "<<count<<" éléments d'un tableau de type LegacyParticleIOStructV4_1 avec un offset de  :  "<<currentOffset<<std::endl; 	
  
	MPI_File_write_at(mpiioFileName,currentOffset,particlesArrayV4_1,count,ParticlesMPITypeV4_1,&status);
  
	currentOffset+=count*ParticlesOffsetV4_1; 
	//std::cout<<"offset après écriture :  "<<currentOffset<<std::endl; 	

}


void LegacySystemIOFile::writeArrayOfParticlesV4_2(LegacyParticleIOStructV4_2 *particlesArrayV4_2, int count){

 	 if(currentOffset == 0)
		currentOffset = VersionOffsetV4_2 + HeaderOffsetV4_2;

	//std::cout<<"On écrit "<<count<<" éléments d'un tableau de type LegacyParticleIOStructV4_2 avec un offset de  :  "<<currentOffset<<std::endl; 	
  
	MPI_File_write_at(mpiioFileName,currentOffset,particlesArrayV4_2,count,ParticlesMPITypeV4_2,&status);
  
	currentOffset+=count*ParticlesOffsetV4_2; 
	//std::cout<<"offset après écriture :  "<<currentOffset<<std::endl; 	

}


void LegacySystemIOFile::writeArrayOfParticlesMolRigV4_1(LegacyParticleIOStructMolRigV4_1 *particlesArrayMolRigV4_1, int count){
	//std::cout<<"On écrit "<<count<<" éléments d'un tableau de type LegacyParticleIOStructMolRigV4_1 option molecule rigide avec un offset de  :  "<<currentOffset<<std::endl; 	
  
	MPI_File_write_at(mpiioFileName,currentOffset,particlesArrayMolRigV4_1,count,ParticlesMPITypeMolRigV4_1,&status);
  
	currentOffset+=count*ParticlesOffsetMolRigV4_1; 
	//std::cout<<"offset après écriture :  "<<currentOffset<<std::endl; 	

}

void LegacySystemIOFile::readArrayOfParticlesMolRigV4_1(LegacyParticleIOStructMolRigV4_1 *particlesArrayMolRigV4_1, int count){
	//std::cout<<"On écrit "<<count<<" éléments d'un tableau de type LegacyParticleIOStructMolRigV4_1 option molecule rigide avec un offset de  :  "<<currentOffset<<std::endl; 	
  
	MPI_File_read_at(mpiioFileName,currentOffset,particlesArrayMolRigV4_1,count,ParticlesMPITypeMolRigV4_1,&status);

	currentOffset+=count*ParticlesOffsetMolRigV4_1; 
	//std::cout<<"offset après écriture :  "<<currentOffset<<std::endl; 	

}

void LegacySystemIOFile::writeArrayOfParticlesMolRigV4_2(LegacyParticleIOStructMolRigV4_2 *particlesArrayMolRigV4_2, int count){
	//std::cout<<"On écrit "<<count<<" éléments d'un tableau de type LegacyParticleIOStructMolRigV4_2 option molecule rigide avec un offset de  :  "<<currentOffset<<std::endl; 	
  
	MPI_File_write_at(mpiioFileName,currentOffset,particlesArrayMolRigV4_2,count,ParticlesMPITypeMolRigV4_2,&status);
  
	currentOffset+=count*ParticlesOffsetMolRigV4_2; 
	//std::cout<<"offset après écriture :  "<<currentOffset<<std::endl; 	

}





// Read the DONNEE and FAtome files to know the Type and Mass stored in the MPIIO Dump
void LegacySystemIOFile::setSimulationData(){
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
    MPI_Abort(MPI_COMM_WORLD,-1);
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
	    MPI_Abort(MPI_COMM_WORLD,-1);
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
    MPI_Abort(MPI_COMM_WORLD,-1);
    exit(-1); 
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

