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

