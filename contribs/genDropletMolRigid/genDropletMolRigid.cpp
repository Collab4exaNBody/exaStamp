#include "StampV3LegacyIOStructures.hpp"
#include <math.h>
#include <string>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <iostream>
#include <climits>

using namespace std;

const long DIMCHAR = 512;
const long PERIODIQUE = 1;
const long LIBRE = 0;
const long NONDEFINI = -1;
const long OUI = 1;
const long NON = 0;
long dimPaquet = 1000000;
const long FreqMessage = 10000000;
long TailleMaxFichier = 1; // en mega-octet

const double boltzman  = 1.380658e-23;
const double avogadro = 6.0221367e23;

const long NOMBRETYPESATOMES=1000;
const long NOMBRETYPESSTRUCTURES=200;
const long NOMBREATOMESSTRUCTURES=50;
const long NOMBRETYPESGEOMETRIES=50;
const long NOMBREPARAMETRESTYPEGEOMETRIE=50;

const long AtomeTest = 106000;

const double pi=acos(-1.);

long idum= 100000;
long idumVitesse= -88855;

double MasseAtomique[NOMBRETYPESATOMES];
char ListeAtomes[NOMBRETYPESATOMES][DIMCHAR];
long NumeroAtomique;

int NumeroGeometrie;
int NbParametresMailleGeometrie[NOMBRETYPESGEOMETRIES];
int ListeTypeGeometrie[NOMBRETYPESGEOMETRIES];

const int GEO_BULK = 0;
const int GEO_ALEATOIRE = 1;
const int GEO_BICOMPOSANT_SELON_X = 2;
const int GEO_RAINURE_SINUS_UN_MATERIAU = 3;
const int GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_SINUS = 4;
const int GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_TRIANGULAIRE = 5;
const int GEO_SPHERE_SEULE = 7;

char ListeGeometries[NOMBRETYPESGEOMETRIES][DIMCHAR];

int NumeroStructure;
int NbAtomesMailleStructure[NOMBRETYPESSTRUCTURES];
double ParametresMaillesDirStructure[NOMBRETYPESSTRUCTURES][NOMBREATOMESSTRUCTURES][3];
int atome;
int ECRITURE_BIBLI=OUI;

char ListeStructures[NOMBRETYPESSTRUCTURES][DIMCHAR];

int main (int argc, char *argv[]) 
{



char NomFichierDonnees[DIMCHAR];
char NomFichierXYZ[DIMCHAR];
char chaine[DIMCHAR];
char ligne[100][DIMCHAR];
char Version[DIMCHAR];
char NomAtome[NOMBRETYPESATOMES][16];
char MoleculeRigide[DIMCHAR];
char OrientationMoleculesRigides[DIMCHAR];
char FichierXYZ[DIMCHAR];
double ParametresMailles[NOMBRETYPESSTRUCTURES][10];
double convLong = 1.e-10;
double masse[NOMBRETYPESATOMES];
double masseTotale;
double CoorAtMaille[NOMBRETYPESATOMES][NOMBREATOMESSTRUCTURES][3];
double temperature[NOMBRETYPESATOMES];
double VitesseAdd[NOMBRETYPESATOMES][3];
double temperatureCalculee[NOMBRETYPESATOMES][3];
double VitesseTotale[NOMBRETYPESATOMES][3];
double VitessePaquet[3];
double beta[NOMBRETYPESATOMES];
double EcinetiqueGlobale;
double EcinetiqueAtome;
double VitesseCG[3];

float reel[100];
int STOP(int);

int F_GEO_ALEATOIRE(int,double *);
int F_GEO_BICOMPOSANT_SELON_X(int,double *,int,double);
int F_GEO_RAINURE_SINUS_UN_MATERIAU(int,double *,double,double);
int F_GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_SINUS(int,double *,int,double,double);
int F_GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_TRIANGULAIRE(int,double *,int,double,double);
int F_SPHERE_SEULE(int,double *,double,double,double);

double gasdev ( long * );
double ran1 ( long * );

int IdentificationAtome(char *);
int IdentificationStructure(char *);
int IdentificationGeometrie(char *);

long compteur;
long numligne;
long NbMailles[NOMBRETYPESATOMES][3];
long iOrientationMoleculesRigides;
long NbMaillesTotal;
long NbAtomesParMaille[NOMBREATOMESSTRUCTURES];
long NbAtomesTotal;
int rank,size;
long CL[3];
long erreur;
long NbParametresGeometrie;
long dimTab;
long TypeGeometrie;
long imasse;
long istructure;
long igeometrie;
long itype;
long autres_atomes;
long Nb_Type;
long compteurMessage;
long NnAtomesMaxi;;
 
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);


 

#include "Structure.h"
#include "Geometrie.h"
#include "MasseAtomique.h"

std::cout << "------------------------------------------ "<< std::endl;	  
std::cout << "Construction d'une configurations initiale "<< std::endl;
std::cout << "------------------------------------------\n\n "<< std::endl;	  
  
  
std::cout << "Nom du fichier de donnees : ConstructionConfigurationInitiale.don "<< std::endl;	   
  
ifstream donnees ("ConstructionConfigurationInitiale.don", ios::in);  

if(!donnees){
	std::cout << "Ouverture impossible du fichier ConstructionConfigurationInitiale.don"<<std::endl;	
	STOP(2);
}  
  
compteur=0;


//MPI_Offset jjj;
//std::cout<<sizeof(jjj)<<std::endl;  
//std::cout<<LONG_MAX<<std::endl;  
//std::cout<<INT_MAX<<std::endl;  
//STOP(888);

// Lecture du fichier de donnees
while(donnees >> ligne[compteur] ){
	std::cout<<ligne[compteur]<<std::endl;
	compteur++;
}  

donnees.close();

std::cout<<compteur<<" lignes lues dans le fichier ConstructionConfigurationInitiale.don"<<std::endl;  


// **************************************************************
// **********************  Atome principal **********************
// **************************************************************
itype = 0;

// Identification du type d'atomes
numligne=0;
strcpy(NomAtome[0],ligne[numligne]);
imasse = IdentificationAtome(ligne[numligne]);
masse[itype] = MasseAtomique[imasse];
masse[itype] /= 1000.;	

// Identification de la structure
numligne=1;
istructure = IdentificationStructure(ligne[numligne]);
NbAtomesParMaille[itype] = NbAtomesMailleStructure[istructure];
std::cout << "Structure numero "<<istructure<<" : "<<ListeStructures[istructure]<<std::endl;	

// Parametres geometriques de la structure
std::cout << "Parametres de maille : ";	
numligne=2;
for(int i=0;i<4;i++){
	sscanf(ligne[numligne+i], "%f",&reel[i] );
	ParametresMailles[itype][i] = (double)reel[i];
	std::cout << ParametresMailles[itype][i]<<"  ";	
}
std::cout <<std::endl;	

std::cout << "Coordonnees des atomes dans la maille : "<<std::endl;	


for(int i=0;i<NbAtomesParMaille[itype];i++){
	std::cout << "\tAtome "<<i<<" : ";	
	for(int j=0;j<3;j++){
		CoorAtMaille[itype][i][j] = ParametresMaillesDirStructure[istructure][i][j] * ParametresMailles[itype][j];
		std::cout << CoorAtMaille[itype][i][j]<<"  ";	
	}
	std::cout <<std::endl;	
}

// Nombre de motifs
numligne=6;
for(int i=0;i<3;i++){
	sscanf(ligne[numligne+i], "%d",&NbMailles[itype][i] );
}

NbMaillesTotal = NbMailles[itype][0] * NbMailles[itype][1] * NbMailles[itype][2];
NbAtomesTotal = NbMaillesTotal * NbAtomesParMaille[itype];

std::cout << "Nombre de mailles : "<<NbMailles[itype][0]<<" x "<<NbMailles[itype][1]<<" x" <<NbMailles[itype][2]<<" = "<<NbMaillesTotal<<", soit "<< NbAtomesTotal <<" atomes"<<std::endl;	     

// Conditions aux limites
numligne=9;

for(int i=0;i<3;i++){
	CL[i] = NONDEFINI;
	if(strcmp(ligne[numligne+i],"Periodique")==0)
		CL[i] = PERIODIQUE;
	if(strcmp(ligne[numligne+i],"Libre")==0)
		CL[i] = LIBRE;
		
	if(CL[i] == NONDEFINI){
		std::cout << "La condition aux limites "<<ligne[numligne+i]<<" n'est pas encore prévu"<<std::endl;	
		STOP(5);  		
	}
}
std::cout<<"Conditions aux limites : "<<ligne[numligne]<<" "<<ligne[numligne+1]<<" "<<ligne[numligne+2]<<std::endl;


// Temperature
numligne=12;
sscanf(ligne[numligne], "%f",&reel[0] );
temperature[itype] = (double)reel[0];
std::cout<<"Temperature : "<<temperature[itype]<<std::endl;  

beta[itype] = masse[itype] / (2.*boltzman*avogadro*temperature[itype]) ;
std::cout<<"beta : "<<beta[itype]<<std::endl;  

// Vitesse additionelle selon x
numligne=13;
sscanf(ligne[numligne], "%f",&reel[0] );
VitesseAdd[itype][0] = (double)reel[0];
std::cout<<"Vitesse additionelle selon X : "<<VitesseAdd[itype][0]<<std::endl;  

// Vitesse additionelle selon y
numligne=14;
sscanf(ligne[numligne], "%f",&reel[0] );
VitesseAdd[itype][1] = (double)reel[0];
std::cout<<"Vitesse additionelle selon Y : "<<VitesseAdd[itype][1]<<std::endl;  

// Vitesse additionelle selon z
numligne=15;
sscanf(ligne[numligne], "%f",&reel[0] );
VitesseAdd[itype][2] = (double)reel[0];
std::cout<<"Vitesse additionelle selon Z : "<<VitesseAdd[itype][2]<<std::endl;  

beta[itype] = masse[itype] / (2.*boltzman*avogadro*temperature[itype]) ;
std::cout<<"beta : "<<beta[itype]<<std::endl;  

// Type de geometrie
numligne=16;

igeometrie = IdentificationGeometrie(ligne[numligne]);
NbParametresGeometrie = NbParametresMailleGeometrie[igeometrie];
TypeGeometrie = ListeTypeGeometrie[igeometrie];;

std::cout<<"Type de geometrie : "<<ligne[numligne]<<" - Nombre de parametres : "<<NbParametresGeometrie<<std::endl;

// Ajout de parametres suplémentaires
if(   TypeGeometrie == GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_SINUS 
   || TypeGeometrie == GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_TRIANGULAIRE)
{
	NbParametresGeometrie++;
}
		
dimTab = NbParametresGeometrie;
if(NbParametresGeometrie==0)dimTab=1;
char ChaineParametresGeometrie[dimTab][DIMCHAR];
double ParametresGeometrie[dimTab];
int IntParametresGeometrie[dimTab];

// Version Stamp
numligne=17;
sscanf(ligne[numligne], "%s",Version);
std::cout<<"Version du format de la protection : "<<Version<<std::endl;

// Molecule rigide
numligne=18;
sscanf(ligne[numligne], "%s",MoleculeRigide);
std::cout<<"Option molécules rigides : "<<MoleculeRigide<<std::endl;

// Creation d'un fichier .xyz
numligne=19;
sscanf(ligne[numligne], "%s",FichierXYZ);
std::cout<<"Creation d'un fichier .xyz : "<<FichierXYZ<<std::endl;

// Non de l'éventuel fichier .xyz
numligne=20;
sscanf(ligne[numligne], "%s",NomFichierXYZ);
if(strcmp(FichierXYZ,"o")==0){
	std::cout<<"Nom du fichier .xyz : "<<NomFichierXYZ<<std::endl;
}
std::ofstream sotieXYZ(NomFichierXYZ,ios::out);
sotieXYZ << "NOMBREATOMES"<<std::endl;
sotieXYZ << "Commentaire"<<std::endl;

// Orientation des molecules rigides
numligne=21;
sscanf(ligne[numligne], "%s",OrientationMoleculesRigides);
if(strcmp(MoleculeRigide,"o")==0){
	std::cout<<"Orientation des molecules rigides : "<<OrientationMoleculesRigides<<std::endl;
}
iOrientationMoleculesRigides=NON;
if(strcmp(OrientationMoleculesRigides,"aleatoire")==0)
	iOrientationMoleculesRigides=OUI;
printf("iOrientationMoleculesRigides=%d\n",iOrientationMoleculesRigides);

// Lecture des parametres de la géométrie 
numligne=22;
for(int i=0;i<NbParametresGeometrie;i++){
	sscanf(ligne[numligne+i], "%s",ChaineParametresGeometrie[i] );
	std::cout << "\tParametre "<<i<<" : "<<ChaineParametresGeometrie[i]<<std::endl;		

}

// Chargement des parametres geometriques pour les cas mono-composant 

if(TypeGeometrie == GEO_ALEATOIRE){
	int j;
	j=0;
	sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
	ParametresGeometrie[j] = (double)reel[0];
}


if(TypeGeometrie == GEO_RAINURE_SINUS_UN_MATERIAU){
	int j;
	j=0;
	sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
	ParametresGeometrie[j] = (double)reel[0];

	j=1;
	sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
	ParametresGeometrie[j] = (double)reel[0];
}

if(TypeGeometrie == GEO_SPHERE_SEULE){
	int j;
	j=0;
	sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
	ParametresGeometrie[j] = (double)reel[0];

	j=1;
	sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
	ParametresGeometrie[j] = (double)reel[0];

	j=2;
	sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
	ParametresGeometrie[j] = (double)reel[0];
	

	j=3;
	sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
	ParametresGeometrie[j] = (double)reel[0];	
	
}

// **************************************************************
// **********************  Autres atomes   **********************
// **************************************************************
if(TypeGeometrie == GEO_BULK || TypeGeometrie == GEO_ALEATOIRE || TypeGeometrie == GEO_RAINURE_SINUS_UN_MATERIAU || TypeGeometrie == GEO_SPHERE_SEULE){
	autres_atomes = NON;
	Nb_Type = 1;
}else if(TypeGeometrie == GEO_BICOMPOSANT_SELON_X 
		|| TypeGeometrie == GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_SINUS 
		|| TypeGeometrie == GEO_RAINURE_SINUS_UN_MATERIAU
		|| TypeGeometrie == GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_TRIANGULAIRE)
{
	autres_atomes = OUI;
}else{
	std::cout << "Erreur de programmation"<<std::endl;		
	STOP(87);
}
if(autres_atomes == OUI){

	int j;
	itype = 1;	
	if(TypeGeometrie == GEO_BICOMPOSANT_SELON_X){
	
		Nb_Type = 2;

		j=0;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		j=3;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		j=4;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];
		
		j=5;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];
		
		j=6;
		sscanf(ChaineParametresGeometrie[j], "%d",&IntParametresGeometrie[j] );

		j=7;
		sscanf(ChaineParametresGeometrie[j], "%d",&IntParametresGeometrie[j] );

		j=8;
		sscanf(ChaineParametresGeometrie[j], "%d",&IntParametresGeometrie[j] );

		j=9;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		ParametresMailles[itype][0] = ParametresGeometrie[3]; 
		ParametresMailles[itype][1] = ParametresGeometrie[4]; 
		ParametresMailles[itype][2] = ParametresGeometrie[5]; 
		
		NbMailles[itype][0] = IntParametresGeometrie[6];
		NbMailles[itype][1] = IntParametresGeometrie[7];
		NbMailles[itype][2] = IntParametresGeometrie[8];

		// Identification de la structure
		j=2;
		istructure = IdentificationStructure(ChaineParametresGeometrie[j]);
		NbAtomesParMaille[itype] = NbAtomesMailleStructure[istructure];

		for(int i=0;i<NbAtomesParMaille[itype];i++){
			for(int j=0;j<3;j++){
				CoorAtMaille[itype][i][j] = ParametresMaillesDirStructure[istructure][i][j] * ParametresMailles[itype][j];
			}
		}
		
		// Identification du type d'atome
		j = 1;
		imasse = IdentificationAtome(ChaineParametresGeometrie[j]);
		masse[itype] = MasseAtomique[imasse];
		masse[itype] /= 1000.;
		strcpy(NomAtome[1],ChaineParametresGeometrie[j]);

		// Initialisation température
		temperature[itype] = ParametresGeometrie[9];
		if(temperature[itype] >0.)beta[itype] = masse[itype] / (2.*boltzman*avogadro*temperature[itype]) ;

	} else if(TypeGeometrie == GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_SINUS || GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_TRIANGULAIRE){
	
		Nb_Type = 2;
		
		// Parametre de maille selon x
		j=2;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		// Parametre de maille selon y
		j=3;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		// Parametre de maille selon z
		j=4;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		// Nombre de mailles selon x
		j=5;
		sscanf(ChaineParametresGeometrie[j], "%d",&IntParametresGeometrie[j] );

		// Nombre de mailles selon y		
		j=6;
		sscanf(ChaineParametresGeometrie[j], "%d",&IntParametresGeometrie[j] );

		// Nombre de mailles selon z
		j=7;
		sscanf(ChaineParametresGeometrie[j], "%d",&IntParametresGeometrie[j] );

		// Profondeur
		j=8;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		// Nombre de rainures
		j=9;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		// Temperature du milieu 2		
		j=10;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		// Taux de presence du milieu 1				
		j=11;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		// Taux de presence du milieu 2		
		j=12;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		// Vx additionnelle du milieu 1						
		j=13;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		// Vy additionnelle du milieu 1						
		j=14;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];
		
		// Vz additionnelle du milieu 1						
		j=15;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		// Vx additionnelle du milieu 2							
		j=16;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		// Vy additionnelle du milieu 2									
		j=17;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];

		// Vz additionnelle du milieu 2									
		j=18;
		sscanf(ChaineParametresGeometrie[j], "%f",&reel[0] );
		ParametresGeometrie[j] = (double)reel[0];
		
		ParametresMailles[itype][0] = ParametresGeometrie[2]; 
		ParametresMailles[itype][1] = ParametresGeometrie[3]; 
		ParametresMailles[itype][2] = ParametresGeometrie[4]; 

		NbMailles[itype][0] = IntParametresGeometrie[5];
		NbMailles[itype][1] = IntParametresGeometrie[6];
		NbMailles[itype][2] = IntParametresGeometrie[7];

		// Identification de la structure
		j=1;
		istructure = IdentificationStructure(ChaineParametresGeometrie[j]);
		NbAtomesParMaille[itype] = NbAtomesMailleStructure[istructure];

		for(int i=0;i<NbAtomesParMaille[itype];i++){
			for(int j=0;j<3;j++){
				CoorAtMaille[itype][i][j] = ParametresMaillesDirStructure[istructure][i][j] * ParametresMailles[itype][j];
			}
		}

		// Identification du type d'atome
		j = 0;
		strcpy(NomAtome[itype],ChaineParametresGeometrie[j]);
		imasse = IdentificationAtome(ChaineParametresGeometrie[j]);
		masse[itype] = MasseAtomique[imasse];
		masse[itype] /= 1000.;
		strcpy(NomAtome[1],ChaineParametresGeometrie[j]);
		
		// Initialisation température
		temperature[itype] = ParametresGeometrie[10];
		if(temperature[itype]>0.)beta[itype] = masse[itype] / (2.*boltzman*avogadro*temperature[itype]) ;
		
		// Initialisation des vitesses additionnelles
		VitesseAdd[0][0] = ParametresGeometrie[13];	
		VitesseAdd[0][1] = ParametresGeometrie[14];	
		VitesseAdd[0][2] = ParametresGeometrie[15];	
		
		VitesseAdd[itype][0] = ParametresGeometrie[16];	
		VitesseAdd[itype][1] = ParametresGeometrie[17];	
		VitesseAdd[itype][2] = ParametresGeometrie[18];	
						  				
	}else{
		std::cout << "La geometrie n'est pas dans la liste"<<std::endl;	
		STOP(20);
	}
}

for(int nType=0; nType<Nb_Type;nType++){
		std::cout<<scientific<<setprecision(8)<<"Masse atomique du type "<<nType<< " : "<<masse[nType]<<" ("<<masse[nType]/avogadro<<" kg)"<<std::endl;  
}

// Calcul du nombre de fichiers en cas de V4.2
long NbAtomesParFichiers_Positions;
if(strcmp(Version,"V4.2")==0){
	TailleMaxFichier *= 1e6;
	std::cout<<"Taille maxi des fichiers : "<<TailleMaxFichier<<std::endl;  

	NbAtomesParFichiers_Positions=TailleMaxFichier / sizeof(LegacyParticleIOStructV4_2);
	std::cout<<"Nombre d'atomes par fichiers _Positions_ : "<<NbAtomesParFichiers_Positions<<std::endl;  

	if(dimPaquet > NbAtomesParFichiers_Positions){
		dimPaquet=NbAtomesParFichiers_Positions;
		std::cout<<" Attention, adaptation de _dimpaquet_ : "<<dimPaquet<<std::endl;  
	}
}

// ************************************************************************
// **********************  Construction du systeme   **********************
// ************************************************************************


double pos[3];
double posMin[NOMBRETYPESSTRUCTURES][3];
double posMax[NOMBRETYPESSTRUCTURES][3];
double posMinGlobal[3];
double posMaxGlobal[3];
double Long[NOMBRETYPESSTRUCTURES][3];
double LongGlobal[3];
double posAtMin[NOMBRETYPESSTRUCTURES][3];
double posAtMax[NOMBRETYPESSTRUCTURES][3];
double posAtMinGlobal[3];
double posAtMaxGlobal[3];
int elegibilite;

LegacyHeaderIOStruct entete ;
LegacyHeaderIOStructV4_1 enteteV4_1 ;
LegacyParticleIOStruct *particlesArray;
LegacyParticleIOStructV4_1 *particlesArrayV4_1;
LegacyParticleIOStructV4_2 *particlesArrayV4_2;
LegacySystemIOFile Protection ;
LegacyParticleIOStructMolRigV4_1 *particlesArrayMolRigV4_1;

int modulo;
int div;
int jmax;
int numlec;
int NbPaquet;
long NbAtomesEcrits[NOMBRETYPESSTRUCTURES];
int compteurPaquet;

for(int i=0;i<3;i++){
	posMinGlobal[i] = 1.e30;
	posMaxGlobal[i] = -posMinGlobal[i];	
	LongGlobal[i] = 0.;	
}			

// Recherche des extrema - toutes les structures sont centrées en (0,0,0)
for(int nType=0; nType<Nb_Type;nType++){
	for(int i=0;i<3;i++){
		Long[nType][i] = NbMailles[nType][i] * ParametresMailles[nType][i];
		posMin[nType][i] = -0.5 * Long[nType][i];
		posMax[nType][i] =  0.5 * Long[nType][i];	
		
		if(posMin[nType][i] < posMinGlobal[i])
			posMinGlobal[i] = posMin[nType][i];

		if(posMax[nType][i] > posMaxGlobal[i])
			posMaxGlobal[i] = posMax[nType][i];

		if(Long[nType][i] > LongGlobal[i])
			LongGlobal[i] = Long[nType][i];

	}
		
	std::cout << "Structure "<<nType<<" : Dimensions du bulk : "<<Long[nType][0]<<" x "<<Long[nType][1]<<" x " <<Long[nType][2]<<std::endl;	     
	std::cout << "Structure "<<nType<<" : Extrema du bulk : "<<posMin[nType][0]<<" <x< "<<posMax[nType][0]<<", "<<posMin[nType][1]<<" <y< "<<posMax[nType][1]<<", "<<posMin[nType][2]<<" <z< "<<posMax[nType][2]<<std::endl;	     	     

}

std::cout <<"Dimensions finales du bulk : "<<LongGlobal[0]<<" x "<<LongGlobal[1]<<" x " <<LongGlobal[2]<<std::endl;	     
std::cout <<"Extrema finals du bulk : "<<posMinGlobal[0]<<" <x< "<<posMaxGlobal[0]<<", "<<posMinGlobal[1]<<" <y< "<<posMaxGlobal[1]<<", "<<posMinGlobal[2]<<" <z< "<<posMaxGlobal[2]<<std::endl;	     	     



// Calculs complementaires selon le type de geometrie
if(TypeGeometrie == GEO_RAINURE_SINUS_UN_MATERIAU){
	ParametresGeometrie[NbParametresGeometrie-1] = LongGlobal[2];
	NbParametresGeometrie++;
	
}

if(TypeGeometrie == GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_SINUS  || TypeGeometrie == GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_TRIANGULAIRE){
	ParametresGeometrie[NbParametresGeometrie-1] = LongGlobal[2];
}



if(TypeGeometrie == GEO_BULK){
	double masseTotale = NbAtomesParMaille[0] * masse[0] / avogadro;
	double VolumeTotal = ParametresMailles[0][0] * ParametresMailles[0][1] * ParametresMailles[0][2] * 1.e-30;
	std::cout << "Masse volumique "<<masseTotale/VolumeTotal<<std::endl;	
}



// Creation du fichier de protection
Protection.open("NomFichierProtection", "w");

if(strcmp(Version,"V3")==0){
	particlesArray = new LegacyParticleIOStruct[dimPaquet];	
}else if(strcmp(Version,"V4.1")==0){
	particlesArrayV4_1 = new LegacyParticleIOStructV4_1[dimPaquet];	
}else if(strcmp(Version,"V4.2")==0){
	particlesArrayV4_2 = new LegacyParticleIOStructV4_2[dimPaquet];	
}else{
	std::cout << " \nerreur - le format "<<Version<<" n'est pas prévu\n"<<std::endl;	
	STOP(124);
}

compteurPaquet = 0;
compteur = 0;
long NbAtomesEcritsGlobal = 0;
EcinetiqueGlobale = 0.;
masseTotale = 0.;
compteurMessage = 0;

for(int i=0;i<3;i++){
	posAtMinGlobal[i] = 1.e30;
	posAtMaxGlobal[i] = -posAtMinGlobal[i];	
	VitesseCG[i] = 0.;
}

/* Calcul du nombre d'atomes maximum */
NnAtomesMaxi = 0;
for(long nType=0; nType<Nb_Type;nType++){
	NnAtomesMaxi += (long)NbMailles[nType][0] * (long)NbMailles[nType][1] * (long)NbMailles[nType][2] * (long)NbAtomesParMaille[nType];
}
std::cout << "Nombre d'atomes maximum :  "<<NnAtomesMaxi<<" soit "<< NnAtomesMaxi/1000000<<" millions"  <<std::endl;	




for(long nType=0; nType<Nb_Type;nType++){

	NbAtomesEcrits[nType] = 0;

	for(int i=0;i<3;i++){
		posAtMin[nType][i] = 1.e30;
		posAtMax[nType][i] = -posAtMin[nType][i];	
		VitesseTotale[nType][i] = 0.;
		temperatureCalculee[nType][i] = 0.;
		VitessePaquet[i] = 0.;
	}		

	for(int i=0;i<NbMailles[nType][0];i++){
		for(int j=0;j<NbMailles[nType][1];j++){
			for(int k=0;k<NbMailles[nType][2];k++){
				for(int l=0;l<NbAtomesParMaille[nType];l++){


					// Calcul de la position de l'atome
					int m;
				
					m=0;
					pos[m] = CoorAtMaille[nType][l][m] + i * ParametresMailles[nType][m] + posMin[nType][m];

					m=1;
					pos[m] = CoorAtMaille[nType][l][m] + j * ParametresMailles[nType][m] + posMin[nType][m];
				
					m=2;
					pos[m] = CoorAtMaille[nType][l][m] + k * ParametresMailles[nType][m] + posMin[nType][m];
					
					// test de son élégibilité
					if(TypeGeometrie == GEO_BULK){		
						elegibilite = OUI;
					}else if(TypeGeometrie == GEO_ALEATOIRE){		
						elegibilite = (int)F_GEO_ALEATOIRE(NbParametresGeometrie,ParametresGeometrie);
					}else if(TypeGeometrie == GEO_BICOMPOSANT_SELON_X){		
						elegibilite = (int)F_GEO_BICOMPOSANT_SELON_X(NbParametresGeometrie,ParametresGeometrie,nType,pos[0]);
					}else if(TypeGeometrie == GEO_RAINURE_SINUS_UN_MATERIAU){	
						elegibilite = (int)F_GEO_RAINURE_SINUS_UN_MATERIAU(NbParametresGeometrie,ParametresGeometrie,pos[0],pos[2]);
					}else if(TypeGeometrie == GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_SINUS){	
						elegibilite = (int)F_GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_SINUS(NbParametresGeometrie,ParametresGeometrie,nType,pos[0],pos[2]);
					}else if(TypeGeometrie == GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_TRIANGULAIRE){	
						elegibilite = (int)F_GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_TRIANGULAIRE(NbParametresGeometrie,ParametresGeometrie,nType,pos[0],pos[2]);
					}else if(TypeGeometrie == GEO_SPHERE_SEULE){	
						elegibilite = (int)F_SPHERE_SEULE(NbParametresGeometrie,ParametresGeometrie,pos[0],pos[1],pos[2]);
					}else{	
						std::cout << "Ce type de geometrie ("<<TypeGeometrie<<") n'est pas prevues"<<std::endl;						
						STOP(10);				
					}

					if(elegibilite == OUI){
						compteur++;


						compteurMessage++;

						if(compteurMessage >= FreqMessage){
							std::cout<<"\t "<<compteur/1000000<<" millions d'atomes sur "<<NnAtomesMaxi<<" ("<<(double)compteur/(double)NnAtomesMaxi*100<<"%)"<<std::endl;			
							compteurMessage=0;
						}
				
						// Mise à jour des extrema
						for(int m=0;m<3;m++){
							if(posAtMin[nType][m] > pos[m])
								posAtMin[nType][m] = pos[m];
							if(posAtMax[nType][m] < pos[m])
								posAtMax[nType][m] = pos[m];
						}

						// Stockage dans la structure d'écriture
						if(strcmp(Version,"V3")==0){
							for(int m=0;m<3;m++){
								particlesArray[compteurPaquet].coordinates[m] = pos[m] * convLong;
								if(temperature[itype]>0.)particlesArray[compteurPaquet].velocity[m] = gasdev(&idumVitesse) * sqrt(1./(2.*beta[nType]));
								if(temperature[itype]==0.)particlesArray[compteurPaquet].velocity[m] = 0.;
								particlesArray[compteurPaquet].particleID = compteur;
								particlesArray[compteurPaquet].particleType = nType;
								VitessePaquet[m] += particlesArray[compteurPaquet].velocity[m];

								// Molecule rigide en V3
								if(strcmp(MoleculeRigide,"o")==0){

									if(iOrientationMoleculesRigides == OUI){
										double nQ=0.;
										for(int mm=0;mm<4;mm++){
											particlesArray[compteurPaquet].quaternion[mm] = ran1(&idum);
											nQ += particlesArray[compteurPaquet].quaternion[mm]*particlesArray[compteurPaquet].quaternion[mm];
										}
										nQ = sqrt(nQ);	
										for(int mm=0;mm<4;mm++){
											particlesArray[compteurPaquet].quaternion[mm] /= nQ;
										}		
									}else{
										particlesArray[compteurPaquet].quaternion[0] = 0.;
										for(int mm=1;mm<4;mm++){
											particlesArray[compteurPaquet].quaternion[mm] = 0.;
										}
									}
									for(int mm=0;mm<3;mm++){
										particlesArray[compteurPaquet].momentumAngular[mm] = 0.;
										particlesArray[compteurPaquet].orientation[mm] = 0.;
									}
								}
							}
						}else if(strcmp(Version,"V4.1")==0){
							for(int m=0;m<3;m++){
								particlesArrayV4_1[compteurPaquet].coordinates[m] = pos[m] / LongGlobal[m];

								if(fabs(particlesArrayV4_1[compteurPaquet].coordinates[m]) > 0.5){
									std::cout << "Erreur dans le calcul des positions"<<std::endl;	
									STOP(256);				
								}

								if(temperature[itype]>0.)particlesArrayV4_1[compteurPaquet].velocity[m] = gasdev(&idumVitesse) * sqrt(1./(2.*beta[nType]));
								if(temperature[itype]==0.)particlesArrayV4_1[compteurPaquet].velocity[m] = 0.;
								VitessePaquet[m] += particlesArrayV4_1[compteurPaquet].velocity[m];
							}
							particlesArrayV4_1[compteurPaquet].charge = 0.;
							particlesArrayV4_1[compteurPaquet].particleID = compteur;
							strcpy(particlesArrayV4_1[compteurPaquet].particleType,NomAtome[nType]);
						}else if(strcmp(Version,"V4.2")==0){
							for(int m=0;m<3;m++){
								particlesArrayV4_2[compteurPaquet].coordinates[m] = pos[m] / LongGlobal[m];

								if(fabs(particlesArrayV4_2[compteurPaquet].coordinates[m]) > 0.5){
									std::cout << "Erreur dans le calcul des positions"<<std::endl;	
									STOP(256);				
								}

								if(temperature[itype]>0.)particlesArrayV4_2[compteurPaquet].velocity[m] = gasdev(&idumVitesse) * sqrt(1./(2.*beta[nType]));
								if(temperature[itype]==0.)particlesArrayV4_2[compteurPaquet].velocity[m] = 0.;
								VitessePaquet[m] += particlesArrayV4_2[compteurPaquet].velocity[m];
							}
							particlesArrayV4_2[compteurPaquet].charge = 0.;
							particlesArrayV4_2[compteurPaquet].particleID = compteur;
							strcpy(particlesArrayV4_2[compteurPaquet].particleType,NomAtome[nType]);

						}else{
							std::cout << " \nerreur - le format "<<Version<<" n'est pas prévu\n"<<std::endl;	
							STOP(125);
						}

						compteurPaquet++;

						// Ecriture du paquet
						if(compteurPaquet == dimPaquet){

							// Annulation de la vitesse du centre de gravité
							for(int k=0;k<compteurPaquet;k++){
							
								if(strcmp(Version,"V3")==0){
									for(int m=0;m<3;m++){
										particlesArray[k].velocity[m] -= VitessePaquet[m] / compteurPaquet;
										VitesseTotale[nType][m] += particlesArray[k].velocity[m];
										EcinetiqueAtome = 0.5 * masse[nType] * particlesArray[k].velocity[m] * particlesArray[k].velocity[m] / avogadro;
									}
								}else if(strcmp(Version,"V4.1")==0){
									for(int m=0;m<3;m++){
										particlesArrayV4_1[k].velocity[m] -= VitessePaquet[m] / compteurPaquet;
										VitesseTotale[nType][m] += particlesArrayV4_1[k].velocity[m];
										EcinetiqueAtome = 0.5 * masse[nType] * particlesArrayV4_1[k].velocity[m] * particlesArrayV4_1[k].velocity[m] / avogadro;
									}
								}else if(strcmp(Version,"V4.2")==0){
									for(int m=0;m<3;m++){
										particlesArrayV4_2[k].velocity[m] -= VitessePaquet[m] / compteurPaquet;
										VitesseTotale[nType][m] += particlesArrayV4_2[k].velocity[m];
										EcinetiqueAtome = 0.5 * masse[nType] * particlesArrayV4_2[k].velocity[m] * particlesArrayV4_2[k].velocity[m] / avogadro;
									}
								}else{
									std::cout << " \nerreur - le format "<<Version<<" n'est pas prévu\n"<<std::endl;	
									STOP(126);
								}
								
								for(int m=0;m<3;m++){
									temperatureCalculee[nType][m] += EcinetiqueAtome;
									EcinetiqueGlobale += EcinetiqueAtome;
								}
							}
							
							// Vitesses additionnelles
							if(strcmp(Version,"V3")==0){
								for(int k=0;k<compteurPaquet;k++){
									for(int m=0;m<3;m++){
										particlesArray[k].velocity[m] += VitesseAdd[nType][m];
									}
								}
							
							}else if(strcmp(Version,"V4.1")==0){
								for(int k=0;k<compteurPaquet;k++){
									for(int m=0;m<3;m++){
										particlesArrayV4_1[k].velocity[m] += VitesseAdd[nType][m];
									}
								}
							}else{
								std::cout << " \nerreur - le format "<<Version<<" n'est pas prévu\n"<<std::endl;	
								STOP(127);
							}

							// Affichage atome test
							if(AtomeTest > 0){
								for(int jj=0;jj<compteurPaquet;jj++){
									if(strcmp(Version,"V3")==0){
										if(particlesArray[jj].particleID == AtomeTest){
											double NQ = 0.;
											for(int mm=0;mm<4;mm++){
												NQ += particlesArray[jj].quaternion[mm]*particlesArray[jj].quaternion[mm];
											}
											std::cout<< "###### Atome "<<AtomeTest<<" - Positions : "<<particlesArray[jj].coordinates[0]<<" "<<particlesArray[jj].coordinates[1]<<" "<<particlesArray[jj].coordinates[2]
											<<" - Vitesse : "<<particlesArray[jj].velocity[0]<<" "<<particlesArray[jj].velocity[1]<<" "<<particlesArray[jj].velocity[2]
											<<" - Quaternion : "<<particlesArray[jj].quaternion[0]<<" "<<particlesArray[jj].quaternion[1]<<" "<<particlesArray[jj].quaternion[2]<<" "<<particlesArray[jj].quaternion[3]
											<<" - Module du quaternion : "<<NQ
											<<std::endl;	

										}
									}else if(strcmp(Version,"V4.1")==0){
										if(particlesArrayV4_1[jj].particleID == AtomeTest){
											std::cout<< "###### Atome "<<particlesArrayV4_1[jj].particleID<<"  Nom : "<<particlesArrayV4_1[jj].particleType<<"  - Positions : "<<particlesArrayV4_1[jj].coordinates[0]<<" "<<particlesArrayV4_1[jj].coordinates[1]<<" "<<particlesArrayV4_1[jj].coordinates[2]<<"  - Vitesse : "<<particlesArrayV4_1[jj].velocity[0]<<" "<<particlesArrayV4_1[jj].velocity[1]<<" "<<particlesArrayV4_1[jj].velocity[2]<<std::endl;

										}
									}
								}
							}

							// Calcul de la vitesse du centre de gravité
							if(strcmp(Version,"V3")==0){
								for(int k=0;k<compteurPaquet;k++){
									for(int m=0;m<3;m++){
										VitesseCG[m] += masse[nType] * particlesArray[k].velocity[m];
									}
									masseTotale += masse[nType];
								}
							
							}else if(strcmp(Version,"V4.1")==0){
								for(int k=0;k<compteurPaquet;k++){
									for(int m=0;m<3;m++){
										VitesseCG[m] += masse[nType] * particlesArrayV4_1[k].velocity[m];
									}
									masseTotale += masse[nType];
								}
							}else{
								std::cout << " \nerreur - le format "<<Version<<" n'est pas prévu\n"<<std::endl;	
								STOP(140);
							}


							// Ecriture
							if(strcmp(Version,"V3")==0){
								Protection.writeArrayOfParticles(particlesArray, compteurPaquet);	
							}else if(strcmp(Version,"V4.1")==0){
								Protection.writeArrayOfParticlesV4_1(particlesArrayV4_1, compteurPaquet);	
							}else{
								std::cout << " \nerreur - le format "<<Version<<" n'est pas prévu\n"<<std::endl;	
								STOP(128);
							}															


							if(strcmp(FichierXYZ,"o")==0){

								if(strcmp(Version,"V3")==0){
									for(int k=0;k<compteurPaquet;k++){
										sotieXYZ<<NomAtome[0]<<" "<<particlesArray[k].coordinates[0]/convLong<<" "<<particlesArray[k].coordinates[1]/convLong<<" "<<particlesArray[k].coordinates[2]/convLong<<" "<<1.<<" "<<0.<<" "<<0.<<" "<<0.<<std::endl;
									}
								}else if(strcmp(Version,"V4.1")==0){
									for(int k=0;k<compteurPaquet;k++){
										sotieXYZ<<NomAtome[0]<<" "<<particlesArrayV4_1[k].coordinates[0]/convLong<<" "<<particlesArrayV4_1[k].coordinates[1]/convLong<<" "<<particlesArrayV4_1[k].coordinates[2]/convLong<<" "<<1.<<" "<<0.<<" "<<0.<<" "<<0.<<std::endl;
									}
								}
							}
						
							NbAtomesEcrits[nType] += compteurPaquet;	
							compteurPaquet = 0;


						}
					}								
				}
			}
		}
	}
	
	if(compteurPaquet != dimPaquet){
	
		// Annulation de la vitesse du centre de gravité
		
		if(strcmp(Version,"V3")==0){
			for(int k=0;k<compteurPaquet;k++){
				for(int m=0;m<3;m++){			
					particlesArray[k].velocity[m] -= VitessePaquet[m] / compteurPaquet;
					VitesseTotale[nType][m] += particlesArray[k].velocity[m];
					EcinetiqueAtome = 0.5 * masse[nType] * particlesArray[k].velocity[m] * particlesArray[k].velocity[m] / avogadro;
					temperatureCalculee[nType][m] += EcinetiqueAtome;
					EcinetiqueGlobale += EcinetiqueAtome;
				}
			}
		}else if(strcmp(Version,"V4.1")==0){
			for(int k=0;k<compteurPaquet;k++){
				for(int m=0;m<3;m++){			
					particlesArrayV4_1[k].velocity[m] -= VitessePaquet[m] / compteurPaquet;
					VitesseTotale[nType][m] += particlesArrayV4_1[k].velocity[m];
					EcinetiqueAtome = 0.5 * masse[nType] * particlesArrayV4_1[k].velocity[m] * particlesArrayV4_1[k].velocity[m] / avogadro;
					temperatureCalculee[nType][m] += EcinetiqueAtome;
					EcinetiqueGlobale += EcinetiqueAtome;
				}
			}
		}else{
			std::cout << " \nerreur - le format "<<Version<<" n'est pas prévu\n"<<std::endl;	
			STOP(129);
		}			

		for(int m=0;m<3;m++){
			VitessePaquet[m] = 0.;
		}
		
		// Vitesses additionnelles
		if(strcmp(Version,"V3")==0){
			for(int k=0;k<compteurPaquet;k++){
				for(int m=0;m<3;m++){
					particlesArray[k].velocity[m] += VitesseAdd[nType][m];
				}
			}
							
		}else if(strcmp(Version,"V4.1")==0){
			for(int k=0;k<compteurPaquet;k++){
				for(int m=0;m<3;m++){
					particlesArrayV4_1[k].velocity[m] += VitesseAdd[nType][m];
				}
			}
		}else{
			std::cout << " \nerreur - le format "<<Version<<" n'est pas prévu\n"<<std::endl;	
			STOP(130);
		}	


		// Affichage atome test
		if(AtomeTest > 0){
		for(int jj=0;jj<compteurPaquet;jj++){
			if(strcmp(Version,"V3")==0){
				if(particlesArray[jj].particleID == AtomeTest){
					std::cout<< "###### Atome "<<AtomeTest<<" Positions : "<<particlesArray[jj].coordinates[0]<<" "<<particlesArray[jj].coordinates[1]<<" "<<particlesArray[jj].coordinates[2]<<" Vitesse : "<<particlesArray[jj].velocity[0]<<" "<<particlesArray[jj].velocity[1]<<" "<<particlesArray[jj].velocity[2]<<std::endl;	
				}
			}else if(strcmp(Version,"V4.1")==0){
				if(particlesArrayV4_1[jj].particleID == AtomeTest){
					std::cout<< std::endl<<"###### Atome "<<particlesArrayV4_1[jj].particleID<<"  Nom : "<<particlesArrayV4_1[jj].particleType<<"  - Positions : "<<particlesArrayV4_1[jj].coordinates[0]<<" "<<particlesArrayV4_1[jj].coordinates[1]<<" "<<particlesArrayV4_1[jj].coordinates[2]<<"  - Vitesse : "<<particlesArrayV4_1[jj].velocity[0]<<" "<<particlesArrayV4_1[jj].velocity[1]<<" "<<particlesArrayV4_1[jj].velocity[2]<<"           "<<NomAtome[nType]<<std::endl<<std::endl;
				}
			}
		}
		}
		
		// Calcul de la vitesse du centre de gravité
		if(strcmp(Version,"V3")==0){
			for(int k=0;k<compteurPaquet;k++){
				for(int m=0;m<3;m++){
					VitesseCG[m] += masse[nType] * particlesArray[k].velocity[m];
				}
				masseTotale += masse[nType];
			}
							
		}else if(strcmp(Version,"V4.1")==0){
			for(int k=0;k<compteurPaquet;k++){
				for(int m=0;m<3;m++){
					VitesseCG[m] += masse[nType] * particlesArrayV4_1[k].velocity[m];
				}
				masseTotale += masse[nType];
			}
		}else{
			std::cout << " \nerreur - le format "<<Version<<" n'est pas prévu\n"<<std::endl;	
			STOP(141);
		}

		if(strcmp(Version,"V3")==0){
			Protection.writeArrayOfParticles(particlesArray, compteurPaquet);	
		}else if(strcmp(Version,"V4.1")==0){
			Protection.writeArrayOfParticlesV4_1(particlesArrayV4_1, compteurPaquet);	
		}else{
			std::cout << " \nerreur - le format "<<Version<<" n'est pas prévu\n"<<std::endl;	
			STOP(131);
		}							

		if(strcmp(FichierXYZ,"o")==0){

			if(strcmp(Version,"V3")==0){
				for(int k=0;k<compteurPaquet;k++){
					sotieXYZ<<NomAtome[0]<<" "<<particlesArray[k].coordinates[0]/convLong<<" "<<particlesArray[k].coordinates[1]/convLong<<" "<<particlesArray[k].coordinates[2]/convLong<<" "<<1.<<" "<<0.<<" "<<0.<<" "<<0.<<std::endl;
				}
			}else if(strcmp(Version,"V4.1")==0){
				for(int k=0;k<compteurPaquet;k++){
					sotieXYZ<<NomAtome[0]<<" "<<particlesArrayV4_1[k].coordinates[0]/convLong<<" "<<particlesArrayV4_1[k].coordinates[1]/convLong<<" "<<particlesArrayV4_1[k].coordinates[2]/convLong<<" "<<1.<<" "<<0.<<" "<<0.<<" "<<0.<<std::endl;
				}
			}
		}
							
		NbAtomesEcrits[nType] += compteurPaquet;	
		compteurPaquet = 0;
	}

	std::cout << "Type  "<<nType<<std::endl;
	std::cout << "Nombre d'atomes écrits pour ce type : "<<NbAtomesEcrits[nType]<<std::endl;
	if(NbAtomesEcrits[nType] > 0){	     
		std::cout << "Extrema des atomes : "<<posAtMin[nType][0]<<" <x< "<<posAtMax[nType][0]<<", "<<posAtMin[nType][1]<<" <y< "<<posAtMax[nType][1]<<", "<<posAtMin[nType][2]<<" <z< "<<posAtMax[nType][2]<<std::endl;
		std::cout << "Vitesse du centre de gravité : "<<VitesseTotale[nType][0]/NbAtomesEcrits[nType]<<" "<<VitesseTotale[nType][1]/NbAtomesEcrits[nType]<<" "<<VitesseTotale[nType][2]/NbAtomesEcrits[nType]<<std::endl;
		std::cout << "Energie cinetique /atome : "<<temperatureCalculee[nType][0]/NbAtomesEcrits[nType]<<" "<<temperatureCalculee[nType][1]/NbAtomesEcrits[nType]<<" "<<temperatureCalculee[nType][2]/NbAtomesEcrits[nType]<<std::endl;	
		std::cout << "Temperature : "<<2.*temperatureCalculee[nType][0]/(boltzman * NbAtomesEcrits[nType])<<" "<<2.*temperatureCalculee[nType][1]/(boltzman * NbAtomesEcrits[nType])<<" "<<2.*temperatureCalculee[nType][2]/(boltzman * NbAtomesEcrits[nType])<<std::endl;
	}	
	NbAtomesEcritsGlobal += NbAtomesEcrits[nType];
	
	for(int m=0;m<3;m++){
		if(posAtMinGlobal[m] > posAtMin[nType][m])
			posAtMinGlobal[m] = posAtMin[nType][m];
		if(posAtMaxGlobal[m] < posAtMax[nType][m])
			posAtMaxGlobal[m] = posAtMax[nType][m];
	}
						
}

std::cout << "Nombre final d'atomes écrits : "<<NbAtomesEcritsGlobal<<std::endl;	     
std::cout << "Extrema finals des atomes  (ang) : "<<posAtMinGlobal[0]<<" <x< "<<posAtMaxGlobal[0]<<", "<<posAtMinGlobal[1]<<" <y< "<<posAtMaxGlobal[1]<<", "<<posAtMinGlobal[2]<<" <z< "<<posAtMaxGlobal[2]<<std::endl;	
std::cout << "Energie cinétique totale : "<<EcinetiqueGlobale<<std::endl;	
std::cout << "Vitesse du centre de gravité: "<<VitesseCG[0]/masseTotale<<" "<<VitesseCG[1]/masseTotale<<" "<<VitesseCG[2]/masseTotale<<std::endl;	
std::ofstream EcritureNbAtomes("fort.1",ios::out);
EcritureNbAtomes << NbAtomesEcritsGlobal<< std::endl;	


if(strcmp(Version,"V3")==0){
	delete[]  particlesArray;
}else if(strcmp(Version,"V4.1")==0){
	delete[]  particlesArrayV4_1;
}else{
	std::cout << " \nerreur - le format "<<Version<<" n'est pas prévu\n"<<std::endl;	
	STOP(131);
}

// Molecule rigides en V4.1

if(strcmp(MoleculeRigide,"o")==0){
	std::cout << "Ecriture des quaternions "<<std::endl;	
	if(strcmp(Version,"V4.1")==0){

		int NbPaquet = NbAtomesEcritsGlobal / dimPaquet;
		int reste = NbAtomesEcritsGlobal - dimPaquet*NbPaquet;
		compteur = 0;
		particlesArrayMolRigV4_1 = new LegacyParticleIOStructMolRigV4_1[dimPaquet];		
		for(int i=0;i<NbPaquet;i++){
			for(int j=0;j<dimPaquet;j++){
				compteur++;
				if(iOrientationMoleculesRigides == OUI){

					double nQ=0.;
					for(int m=0;m<4;m++){
						particlesArrayMolRigV4_1[j].quaternion[m] = ran1(&idum);
						nQ += particlesArrayMolRigV4_1[j].quaternion[m]*particlesArrayMolRigV4_1[j].quaternion[m];
					}
					nQ = sqrt(nQ);	
					for(int m=0;m<4;m++){
						particlesArrayMolRigV4_1[j].quaternion[m] /= nQ;
					}		
				}else{

					particlesArrayMolRigV4_1[j].quaternion[0] = 1.;
					for(int m=1;m<4;m++){
						particlesArrayMolRigV4_1[j].quaternion[m] = 0.;
					}
				}
				for(int m=0;m<3;m++){
					particlesArrayMolRigV4_1[j].momentangulaire[m] = 0.;
					//particlesArrayMolRigV4_1[j].orientation[m] = 0.;
				}

				// Affichage atome test
				if(AtomeTest > 0){
				if(strcmp(Version,"V4.1")==0){
					if(compteur == AtomeTest){
							std::cout<< "###### Atome "<<compteur<<"  - Quaternion : "<<particlesArrayMolRigV4_1[j].quaternion[0]<<" "<<particlesArrayMolRigV4_1[j].quaternion[1]<<" "<<particlesArrayMolRigV4_1[j].quaternion[2]<<" "<<particlesArrayMolRigV4_1[j].quaternion[3]<<std::endl;
					}
				}
				}
			}
			Protection.writeArrayOfParticlesMolRigV4_1(particlesArrayMolRigV4_1, dimPaquet);	
		}
		delete[]  particlesArrayMolRigV4_1;

		if(reste > 0){
			particlesArrayMolRigV4_1 = new LegacyParticleIOStructMolRigV4_1[reste];		
			for(int j=0;j<reste;j++){
				compteur++;
				if(iOrientationMoleculesRigides == OUI){

					double nQ=0.;
					for(int m=0;m<4;m++){
						particlesArrayMolRigV4_1[j].quaternion[m] = ran1(&idum);
						nQ += particlesArrayMolRigV4_1[j].quaternion[m]*particlesArrayMolRigV4_1[j].quaternion[m];
					}
					nQ = sqrt(nQ);	
					for(int m=0;m<4;m++){
						particlesArrayMolRigV4_1[j].quaternion[m] /= nQ;
					}		
				}else{
					particlesArrayMolRigV4_1[j].quaternion[0] = 1.;
					for(int m=1;m<4;m++){
						particlesArrayMolRigV4_1[j].quaternion[m] = 0.;
					}
				}
				for(int m=0;m<3;m++){
					particlesArrayMolRigV4_1[j].momentangulaire[m] = 0.;
					//particlesArrayMolRigV4_1[j].momentangulaire[m] = 1.e-35;
					//if(ran1(&idum) > 0.5)particlesArrayMolRigV4_1[j].momentangulaire[m]*= -1.;
					//particlesArrayMolRigV4_1[j].orientation[m] = 0.;
				}

				// Affichage atome test
				if(AtomeTest > 0){
				if(strcmp(Version,"V4.1")==0){
					if(compteur == AtomeTest){
							std::cout<< "###### Atome "<<compteur<<"  - Quaternion : "<<particlesArrayMolRigV4_1[j].quaternion[0]<<" "<<particlesArrayMolRigV4_1[j].quaternion[1]<<" "<<particlesArrayMolRigV4_1[j].quaternion[2]<<" "<<particlesArrayMolRigV4_1[j].quaternion[3]<<std::endl;
					}
				}
				}
			}


			Protection.writeArrayOfParticlesMolRigV4_1(particlesArrayMolRigV4_1, reste);				
			delete[]  particlesArrayMolRigV4_1;		
		}

	}		
}

// En-tete

if(strcmp(Version,"V3")==0){
	entete.iterationNumber = 0;
	entete.particlesTotalNumber =  NbAtomesEcritsGlobal;
	entete.time = 0.0;
	
	if(CL[0] == LIBRE){
		entete.xmin = posAtMinGlobal[0] * convLong;
		entete.xmax = posAtMaxGlobal[0] * convLong;
	}

	if(CL[0] == PERIODIQUE){
		entete.xmin = posMinGlobal[0] * convLong;
		entete.xmax = posMaxGlobal[0] * convLong;
	}

	if(CL[1] == LIBRE){
		entete.ymin = posAtMinGlobal[1] * convLong;
		entete.ymax = posAtMaxGlobal[1] * convLong;
	}

	if(CL[1] == PERIODIQUE){
		entete.ymin = posMinGlobal[1] * convLong;
		entete.ymax = posMaxGlobal[1] * convLong;
	}

	if(CL[2] == LIBRE){
		entete.zmin = posAtMinGlobal[2] * convLong;
		entete.zmax = posAtMaxGlobal[2] * convLong;
	}

	if(CL[2] == PERIODIQUE){
		entete.zmin = posMinGlobal[2] * convLong;
		entete.zmax = posMaxGlobal[2] * convLong;
	}

	entete.totalEnergy = 0.;
	entete.potentialEnergy = 0.;
	entete.internalEnergy = 0.;
	entete.kineticEnergy = EcinetiqueGlobale;
	entete.rotationalEnergy = 0.;

	Protection.writeHeader(entete);
	Protection.close();	 
}else if(strcmp(Version,"V4.1")==0){



	enteteV4_1.Natomes = NbAtomesEcritsGlobal;
	
	enteteV4_1.long_a=LongGlobal[0] * convLong;
	enteteV4_1.long_b=LongGlobal[1] * convLong;
	enteteV4_1.long_c=LongGlobal[2] * convLong;

	enteteV4_1.angle_a=90.;
	enteteV4_1.angle_b=90.;
	enteteV4_1.angle_g=90.;
		
	enteteV4_1.MatriceCR[0][0] = 1.;
	enteteV4_1.MatriceCR[0][1] = 0.;
	enteteV4_1.MatriceCR[0][2] = 0.;
		
	enteteV4_1.MatriceCR[1][0] = 0.;
	enteteV4_1.MatriceCR[1][1] = 1.;
	enteteV4_1.MatriceCR[1][2] = 0.;
			
	enteteV4_1.MatriceCR[2][0] = 0.;
	enteteV4_1.MatriceCR[2][1] = 0.;
	enteteV4_1.MatriceCR[2][2] = 1.;	

	enteteV4_1.XCGeo = 0.;
	enteteV4_1.YCGeo = 0.;
	enteteV4_1.ZCGeo = 0.;

	//donnees temporelles
	enteteV4_1.NumeroIterationAbsolu =  0;
	enteteV4_1.tempsPhysique = 0;
		
	//donnees energetiques
	enteteV4_1.EnergieTotale = 0;
	enteteV4_1.EnergiePotentielle = 0;
	enteteV4_1.EnergieCinetique = EcinetiqueGlobale;
	enteteV4_1.EnergieRotationnelle = 0.;

	//invariants
	enteteV4_1.invariant = 5.;

	//parametres specifiques de la simulation
	enteteV4_1.dt_adaptatif		= 1e-15;
	enteteV4_1.LNVhug_T		= 0.;
	enteteV4_1.LNVhug_Tref		= 0.;
	for (int i=0;i<3;i++) {
		enteteV4_1.NVT_gamma[i]	= 0.;
		enteteV4_1.NVT_gammap[i] = 0.;

		enteteV4_1.NPT_qsi[i] = 0.;
		enteteV4_1.NPT_qsip[i] = 0.;

		enteteV4_1.NPH_omega[i] = 0.;
		enteteV4_1.NPH_omegap[i] = 0.;

		enteteV4_1.NPH_pi[i] = 0.;
	}

	enteteV4_1.bloc_molecules = 0.;
	enteteV4_1.bloc_dpd = 0.;
	enteteV4_1.bloc_graines = 0.;
	enteteV4_1.bloc_molrig = 0;

	if(strcmp(MoleculeRigide,"o")==0){
		enteteV4_1.bloc_molrig = 1;
	}

	enteteV4_1.bloc_monomeres = 0.;
	enteteV4_1.bloc_polymerisation = 0.;
	enteteV4_1.bloc_posfiltre = 0.;

	LegacyVersionIOStructV4_1 VersionData;
	VersionData.Version = 2; 

	enteteV4_1.bloc_Ijohndoe01 = 0;
	enteteV4_1.bloc_Ijohndoe02 = 0;
	enteteV4_1.bloc_Ijohndoe03 = 0;
	enteteV4_1.bloc_Ijohndoe04 = 0;
	enteteV4_1.bloc_Ijohndoe05 = 0;
	enteteV4_1.bloc_Ijohndoe06 = 0;
	enteteV4_1.bloc_Ijohndoe07 = 0;

	enteteV4_1.bloc_Rjohndoe01 = 0.;
	enteteV4_1.bloc_Rjohndoe02 = 0.;
	enteteV4_1.bloc_Rjohndoe03 = 0.;
	enteteV4_1.bloc_Rjohndoe04 = 0.;
	enteteV4_1.bloc_Rjohndoe05 = 0.;
	enteteV4_1.bloc_Rjohndoe06 = 0.;
	enteteV4_1.bloc_Rjohndoe07 = 0.;
	enteteV4_1.bloc_Rjohndoe08 = 0.;
	enteteV4_1.bloc_Rjohndoe09 = 0.;
	enteteV4_1.bloc_Rjohndoe10 = 0.;

	Protection.writeVersionNumber(VersionData);
	Protection.writeHeaderV4_1(enteteV4_1);
	Protection.close();	
			
}else{
	std::cout << " \nerreur - le format "<<Version<<" n'est pas prévu\n"<<std::endl;	
	STOP(123);
}



     
STOP(0);
return 0;
}

// #######################################################
int STOP(int erreur){

MPI_Finalize();	
if(erreur == 0){
	std::cout << " \n##### Fin normale du programme #######\n"<<std::endl;
}else{
	std::cout << " \n##### Fin du programme avec le code erreur "<<erreur<<" #######\n"<<std::endl;
}
exit(0);

}

// #######################################################

int IdentificationAtome(char NomAtome[512]){



int test=0;
int typeAtome=0;

for(int i=0;i<NOMBRETYPESATOMES && test == 0;i++){

	if(strcmp(NomAtome,ListeAtomes[i])==0){
		test = 1;
		typeAtome = i;
	}

}

if(test == 0){
	std::cout << "L'atome "<<NomAtome<<" n'est pas dans la liste"<<std::endl;
	STOP(99);
}

return typeAtome;

}

// #######################################################

int IdentificationStructure(char Structure[512]){



int test=0;
int typeStructure=0;

for(int i=0;i<NOMBRETYPESSTRUCTURES && test == 0;i++){

	if(strcmp(Structure,ListeStructures[i])==0){
		test = 1;
		typeStructure = i;
	}

}

if(test == 0){
	std::cout << "La structure "<<Structure<<" n'est pas dans la liste"<<std::endl;
	STOP(98);
}

return typeStructure;

}

// #######################################################

int IdentificationGeometrie(char Geometrie[512]){



int test=0;
int typeGeometrie=0;

for(int i=0;i<NOMBRETYPESGEOMETRIES && test == 0;i++){

	if(strcmp(Geometrie,ListeGeometries[i])==0){
		test = 1;
		typeGeometrie = i;
	}

}

if(test == 0){
	std::cout << "La geometrie "<<Geometrie<<" n'est pas dans la liste"<<std::endl;
	STOP(97);
}

return typeGeometrie;

}
// #######################################################
int F_GEO_ALEATOIRE(int NbParametresGeometrie,double Parametres[]){

double ran1 ( long * );
int elegibilite;
double nbalea;



nbalea= ran1(&idum);

if(nbalea > Parametres[0]){
	elegibilite = NON;
}else{
	elegibilite = OUI;
}

return elegibilite;
}

// #######################################################

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS (1.2E-07)
#define RNMX (1.0-EPS) 

double ran1 ( long *idum ) {
  /**
   * idum must be set as a negative integer to initialize and must never be 
   * modified elsewhere.
   *
   */
  int    j;
  long   k;
  double temp;
  static long iy=0;
  static long iv[NTAB];

  if ( (*idum <= 0) || (iy == 0) ) {		// initialize
    if (-(*idum) < 1 ) {			// be sure to prevent idum = 0
      *idum = 1 ;
    } else {
      *idum = -(*idum);
    }
    for(j=NTAB+7;j>=0;j--) {			// load the shuffle table (after 8 warm ups)
      k = *idum/IQ;
      *idum = IA*(*idum-k*IQ)-IR*k;
      if(*idum < 0) *idum += IM;
      if(j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = *idum / IQ;				// start here when not initializing
  *idum = IA * ( *idum - k * IQ ) - IR * k ;	// compute idum without overflows
  if( *idum < 0 ) *idum += IM ;			// 	schrage's method idum=(IA*idum)%IM
  j = iy / NDIV ;				// will be in the range 0 ... NTAB-1
  iy = iv[j] ;					// Output previously stored value and 
  iv[j] = *idum;				//	refill the shuffle table
  if ( (temp = AM * iy ) > RNMX ) {		// because users don't expect endpoint val
    return RNMX ;
  } else {
    return temp ;
  }
}

// #######################################################
double gasdev ( long *idum ) {
  static int    iset = 0 ;
  static double gset ;
  double        fac, rsq, v1, v2 ;
  if (iset == 0) {
    do {
      v1 = 2.0 * ran1(idum) - 1.0 ;
      v2 = 2.0 * ran1(idum) - 1.0 ;
      rsq = v1*v1 + v2*v2 ;
    } while ( rsq >= 1.0 || rsq == 0.0 ) ;
    fac = sqrt( -2.0 * log(rsq) / rsq ) ;
    gset = v1 * fac ;
    iset = 1 ;
    return (double) ( v2 * fac ) ;
  } else {
    iset = 0 ;
    return (double) gset ;
  }
}

// #######################################################
int F_GEO_BICOMPOSANT_SELON_X(int NbParametresGeometrie,double Parametres[],int itype,double x){

int elegibilite;

if(itype == 0){
	if(x > Parametres[0]){
		elegibilite = NON;
	}else{
		elegibilite = OUI;
	}
}

if(itype == 1){
	if(x <= Parametres[0]){
		elegibilite = NON;
	}else{
		elegibilite = OUI;
	}
}

return elegibilite;	
}

// #######################################################
int F_GEO_RAINURE_SINUS_UN_MATERIAU(int NbParametresGeometrie,double Parametres[],double x,double z){

int elegibilite = OUI;

double xLimite;

xLimite = 0.5 * Parametres[0] * (1. - cos(2.*Parametres[1]*pi*z/Parametres[NbParametresGeometrie-1]));

if(x > xLimite){
	elegibilite = NON;
}else{
	elegibilite = OUI;
}

std::cout << " A vérifier "<<std::endl;
STOP(999);

return elegibilite;	
}

// #######################################################
int F_GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_SINUS(int NbParametresGeometrie,double Parametres[],int type,double x,double z){

int elegibilite1 = OUI;
int elegibilite2 = OUI;
int elegibilite  = OUI;

double xLimite;
double nbalea1,nbalea2;

xLimite = 0.5 * Parametres[8] * (1. - cos(2.*Parametres[9]*pi*z/Parametres[NbParametresGeometrie-1]));

//printf("Profondeur : %g\n",Parametres[8]);
//printf("Largeur : %g\n",Parametres[NbParametresGeometrie-1]);
//printf("Nombre de rainures : %g\n",Parametres[9]);

nbalea1= ran1(&idum);
nbalea2= ran1(&idum);

if(type == 0){
	if(x > xLimite){
		elegibilite1 = NON;
	}else{
		elegibilite1 = OUI;
	}
	if(nbalea1 > Parametres[11]){
		elegibilite2 = NON;
	}else{
		elegibilite2 = OUI;
	}	
	
}else if(type == 1){
	if(x > xLimite){
		elegibilite1 = OUI;
	}else{
		elegibilite1 = NON;
	}
	if(nbalea2 > Parametres[12]){
		elegibilite2 = NON;
	}else{
		elegibilite2 = OUI;
	}
	
}else{
	STOP(76);
}


if(elegibilite1 == OUI && elegibilite2 == OUI){
	elegibilite = OUI;
}else{
	elegibilite = NON;
}


return elegibilite;	
}

// #######################################################
int F_GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_TRIANGULAIRE(int NbParametresGeometrie,double Parametres[],int type,double x,double z){

int elegibilite1 = OUI;
int elegibilite2 = OUI;
int elegibilite  = OUI;

double xLimite;
double pente;
double largeur1r,zmin,zorigine;
double nbalea1,nbalea2;

int numRainure,nbRainure;


nbRainure = (int)(Parametres[9] + 0.5);

largeur1r = Parametres[NbParametresGeometrie-1] / Parametres[9];
zmin = -0.5 * Parametres[NbParametresGeometrie-1];

numRainure = (z-zmin)/largeur1r;
if(numRainure > nbRainure-1 || numRainure < 0){
	STOP(888);
}

zorigine = zmin + 0.5 * (2*numRainure + 1 ) * largeur1r;
pente = Parametres[8] / (0.5 * largeur1r);
if(z < zorigine)
	pente = -pente;

xLimite = pente * (z - zorigine);
nbalea1= ran1(&idum);
nbalea2= ran1(&idum);

/*
printf("\n\nProfondeur : %g\n",Parametres[8]);
printf("Largeur totale : %g\n",Parametres[NbParametresGeometrie-1]);
printf("Largeur d'une rainure : %g\n",largeur1r);
printf("Nombre de rainures : %d\n",nbRainure);
printf("Pente : %g\n",pente);
printf("zmin : %g\n",zmin);
printf("x : %g\n",x);
printf("z : %g\n",z);
printf("Numero de la rainure : %d\n",numRainure);
printf("zorigine : %g\n",zorigine);
printf("xLimite : %g\n",xLimite);
if(numRainure==1)STOP(999);
*/



if(type == 0){
	if(x > xLimite){
		elegibilite1 = NON;
	}else{
		elegibilite1 = OUI;
	}
	if(nbalea1 > Parametres[11]){
		elegibilite2 = NON;
	}else{
		elegibilite2 = OUI;
	}	
	
}else if(type == 1){
	if(x > xLimite){
		elegibilite1 = OUI;
	}else{
		elegibilite1 = NON;
	}
	if(nbalea2 > Parametres[12]){
		elegibilite2 = NON;
	}else{
		elegibilite2 = OUI;
	}
	
}else{
	STOP(76);
}


if(elegibilite1 == OUI && elegibilite2 == OUI){
	elegibilite = OUI;
}else{
	elegibilite = NON;
}

return elegibilite;	
}

// #######################################################
int F_SPHERE_SEULE(int NbParametresGeometrie,double Parametres[],double x,double y,double z){

/*
	Parametres[0] : rayon
	Parametres[1] : x centre
	Parametres[2] : y centre
	Parametres[3] : z centre
*/


int elegibilite=NON;
double dx,dy,dz,distanceCarree;

dx = x - Parametres[1];
dy = y - Parametres[2];
dz = z - Parametres[3];

distanceCarree=dx*dx + dy*dy + dz*dz;

if(distanceCarree <= Parametres[0]*Parametres[0]){
	elegibilite=OUI;
}

return elegibilite;
}



