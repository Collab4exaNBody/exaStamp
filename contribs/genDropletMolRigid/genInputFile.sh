#!/bin/bash

####################################################################

cas=Sphere_R20_ConfIni
NettoyagePrealable=o
NomFichierProtection=ConfIni.MpiIO
Materiau=TIP4P
StructureCristallographic=CS
aStr=3.103787678
bStr=3.103787678
cStr=3.103787678
dStr=0.
NbMailleX=150
NbMailleY=150
NbMailleZ=150
CLX=Libre
CLY=Libre
CLZ=Libre
Temperature=1.
VitesseAddX=0.
VitesseAddY=0.
VitesseAddZ=0.
Version=V4.1
MoleculeRigide=o
FichierXYZ=n
NomFichierXYZ=test.xyz
OrientationMoleculeRigide=aleatoire
TypeGeometrie=sphere_seule

if [ $TypeGeometrie == 'sphere_seule' ]
then
	Rayon_sphere_seule=245
	Centre_sphere_seule_X=0
	Centre_sphere_seule_Y=0
	Centre_sphere_seule_Z=0
fi

####################################################################

module purge
module load mpi
pgm=ConstructionConfigurationInitiale

echo
echo '-------------------------------------------------'
echo Construction d\'une configurations initiale
echo '-------------------------------------------------'
echo

echo $Materiau > $pgm.don
echo $StructureCristallographic >> $pgm.don

echo $aStr >> $pgm.don
echo $bStr >> $pgm.don
echo $cStr >> $pgm.don
echo $dStr >> $pgm.don

echo $NbMailleX >> $pgm.don
echo $NbMailleY >> $pgm.don
echo $NbMailleZ >> $pgm.don

echo $CLX >> $pgm.don
echo $CLY >> $pgm.don
echo $CLZ >> $pgm.don

echo $Temperature  >> $pgm.don

echo $VitesseAddX >> $pgm.don
echo $VitesseAddY >> $pgm.don
echo $VitesseAddZ >> $pgm.don

echo $TypeGeometrie  >> $pgm.don
echo $Version  >> $pgm.don

echo $MoleculeRigide >> $pgm.don

echo $FichierXYZ >> $pgm.don
echo $NomFichierXYZ >> $pgm.don
echo $OrientationMoleculeRigide >> $pgm.don

if [ $TypeGeometrie == 'sphere_seule' ]
then

	echo $Rayon_sphere_seule  >>  $pgm.don	
	echo $Centre_sphere_seule_X  >>  $pgm.don	
	echo $Centre_sphere_seule_Y  >>  $pgm.don	
	echo $Centre_sphere_seule_Z  >>  $pgm.don	
fi


