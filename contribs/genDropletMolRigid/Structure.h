

// BCC
NumeroStructure = 1;
strcpy(ListeStructures[NumeroStructure],"BCC");
NbAtomesMailleStructure[NumeroStructure] = 2;

atome = 0;
ParametresMaillesDirStructure[NumeroStructure][atome][0] = 0.;
ParametresMaillesDirStructure[NumeroStructure][atome][1] = 0.;
ParametresMaillesDirStructure[NumeroStructure][atome][2] = 0.;

atome = 1;
ParametresMaillesDirStructure[NumeroStructure][atome][0] = 0.5;
ParametresMaillesDirStructure[NumeroStructure][atome][1] = 0.5;
ParametresMaillesDirStructure[NumeroStructure][atome][2] = 0.5;

// CS
NumeroStructure = 2;
strcpy(ListeStructures[NumeroStructure],"CS");
NbAtomesMailleStructure[NumeroStructure] = 1;

atome = 0;
ParametresMaillesDirStructure[NumeroStructure][atome][0] = 0.5;
ParametresMaillesDirStructure[NumeroStructure][atome][1] = 0.5;
ParametresMaillesDirStructure[NumeroStructure][atome][2] = 0.5;

if(ECRITURE_BIBLI == OUI){
	for(int i=1;i<3;i++){
		std::cout << "Structure numero "<<i<<std::endl;	
		std::cout << "\tType :  "<<ListeStructures[i]<<std::endl;	
		std::cout << "\tNombre d'atomes :  "<<NbAtomesMailleStructure[i]<<std::endl;	
		std::cout << "\tPositions des atomes en coordonnées réduites  "<<std::endl;	
		for(int j=0;j<NbAtomesMailleStructure[i];j++){
			std::cout << "\tAtome  "<<j<<" : ";	
			for(int k=0;k<3;k++){
				std::cout << ParametresMaillesDirStructure[i][j][k]<<" ";	
			}
			std::cout<<std::endl;	
		}
	}
}

