/// @file
/// @brief Implementations for the particleWriterHercule
// ECom = to be commented or deleted by Estelle

#ifdef __use_lib_hercule

#include <ostream>

#include "referenceMap.hpp"

#include "io/particleOutputHercule.hpp"

#include "parallel/commManager.hpp"
#include "parallel/mympi.hpp"

#include "HIc.h"
#include "HIc_MPI.h"


/// @brief Destructor
///
///
ParticleWriterHercule::~ParticleWriterHercule() {
	if (!m_base.isNull()) {
		closeBase();
	}
	if (!m_api.isNull()) {
		delete m_api;
		m_api = HIc_Api();
	}
}


/// @brief [ECom] ParticleWriterHercule::write
void ParticleWriterHercule::write(int step, ParticleOutput* buff) {

	openBase();

	openStep(step, buff);
	HIc_Obj particleSet;

	writeMetaInfo(m_ctx,step);
	writePositions(m_ctx, particleSet);
	writeIndexes(m_ctx, particleSet);
	if (particles->writeType)
	  writeTypes(m_ctx, particleSet);
	if (particles->writeVelocity)
	  writeVelocities(m_ctx, particleSet);
	if (particles->writeEint)
	  writeEints(m_ctx, particleSet);
	if (particles->writeProgress)
	  writeProgresses(m_ctx, particleSet);
	closeStep();

}


/// @brief [ECom] ParticleWriterHercule::openBase
void ParticleWriterHercule::openBase() {
	if (!m_base.isNull())
		return;

	m_base = HIc_Base(m_api, "myBase", "parallele");

	m_base.setItemConf("mode", "write");
	m_base.setItemConf("write_dir", m_rep.c_str());
	m_base.setItemConf("read_dir", m_rep.c_str());
	m_base.setItemConf("bd_name", "HDep");
	m_base.setItemConf("multi_period", 1);
	m_base.open();
}


/// @brief [ECom] ParticleWriterHercule::closeBase
void ParticleWriterHercule::closeBase() {
	m_base.close();
	m_base = HIc_Base();
}


/// @brief [ECom] ParticleWriterHercule::openStep
void ParticleWriterHercule::openStep(uint _step, ParticleOutput* buff) {
	particles = buff;
	m_ctx = m_base.createCtxPar((double) _step, comm->getNumberOfNodes(),
			comm->getRank());

	m_ctx.setItemConf("mode", "write");

	m_ctx.open();

}


/// @brief [ECom] ParticleWriterHercule::closeStep
void ParticleWriterHercule::closeStep() {
	m_ctx.close();
	m_ctx = HIc_Ctx();
	particles = NULL;
}


/// @brief [ECom] ParticleWriterHercule::writeMetaInfo
void ParticleWriterHercule::writeMetaInfo(HIc_Ctx& _ctx,int step)
{
	HIc_Obj root = _ctx.getRoot();
	root.getAttr("contenu.topoEvol").setVal("ChangeToutLeTemps");
	const char* env_NCAS=getenv("NCAS");
	if (!env_NCAS){
		env_NCAS="etudeInconnu";
	}
	root.getAttr("etude.nom").setVal(env_NCAS);
	root.getAttr("contenu.numero_cycle").setVal(step);
}

/// @brief [ECom] ParticleWriterHercule::writePositions
void ParticleWriterHercule::writePositions(HIc_Ctx& _ctx,
		HIc_Obj& _obj_particleSet) {

	Array<vec3<double> >& r = particles->r;
	uint size = r.size();
	HIc_Obj root = _ctx.getRoot();

	_obj_particleSet = root.createWithSupport("particles", "Maillage-0D-NS");
	_obj_particleSet.setAttrVal("nbNoeuds", size);
	_obj_particleSet.setAttrVal("nbElements", size);

	//_obj_particleSet.createWithSupport("elements","Point[nbElements]");

	double *X = (double*) malloc(sizeof(double) * size);
	double *Y = (double*) malloc(sizeof(double) * size);
	double *Z = (double*) malloc(sizeof(double) * size);

	for (uint p = 0; p < size; ++p) {
		X[p] = r[p][0] * 10.0; // TODO OB: conversion en unite special A
		Y[p] = r[p][1] * 10.0;
		Z[p] = r[p][2] * 10.0;
	}

	_obj_particleSet.setAttrVal("noeuds.coord.X", X, size);
	_obj_particleSet.setAttrVal("noeuds.coord.Y", Y, size);
	_obj_particleSet.setAttrVal("noeuds.coord.Z", Z, size);

	free(X);
	free(Y);
	free(Z);

	/*if (writeCorners) {
	 // ecriture des Bounds global et de chaque sous-domaine
	 double minBounds[3];
	 double maxBounds[3];

	 _obj_particleSet.setAttrVal("global.minBounds",minBounds,3);
	 _obj_particleSet.setAttrVal("global.maxBounds",maxBounds,3);

	 vec3<double> limits[2];
	 limits[0] = Global::domainInfo.getMinBounds();
	 limits[1] = Global::domainInfo.getMaxBounds();

	 for (int i=0; i<2; ++i) {
	 	 for (int j=0; j<2; ++j) {
	 	 	 for (int k=0; k<2; ++k) {
	 	 	 	 vec3<double> R(limits[i].x, limits[j].y, limits[k].z);
	 	 	 	 //flux << R << std::endl;
	 	 	 }
	 	 }
	 }

	 _obj_particleSet.setAttrVal("global.minBounds",minBounds,3);
	 _obj_particleSet.setAttrVal("global.maxBounds",maxBounds,3);

	 }*/

}


/// @brief [ECom] ParticleWriterHercule::writeIndexes
void ParticleWriterHercule::writeIndexes(HIc_Ctx& _ctx,
		HIc_Obj& _obj_particleSet) {

	Array<uint>& id = particles->id;
	uint size = id.size();

	//if (comm->getRank() == 0)
	//std::cout << "ParticleWriterHercule::writeIndexes" << size <<std::endl;
	int_8 *idval = (int_8*) malloc(sizeof(int_8) * size);

	for (uint p = 0; p < size; ++p) {
		idval[p] = id[p];
	}

	_obj_particleSet.setAttrVal("elements.id", idval, size);
	HIc_Obj connNoeud = _obj_particleSet.createAttr("elements.connNoeud",
			"&Noeud[1]", _obj_particleSet.getAttr("noeuds"));
	int_8 *tab_connnNoeuds = (int_8*) malloc(sizeof(int_8) * size);
	for (uint p = 0; p < size; ++p) {
		tab_connnNoeuds[p] = p;
	}
	connNoeud.setVal(tab_connnNoeuds, size);

	free(tab_connnNoeuds);
	free(idval);
}


/// @brief [ECom] ParticleWriterHercule::writeTypes
void ParticleWriterHercule::writeTypes(HIc_Ctx& _ctx,
		HIc_Obj& _obj_particleSet) {

	Array<uint8_t>& type = particles->type;
	uint size = type.size();
	//if (comm->getRank() == 0)
	//std::cout << "ParticleWriterHercule::writeTypes " << size << std::endl;

	int_2 *typesval = (int_2*) malloc(sizeof(int_2) * size);

	//_obj_particleSet.setAttrVal("elements.real_element_type",typesval,size);
	uint8_t nTypes = Global::reference.getNumberOfTypes();
	std::vector < uint > nbZ(nTypes, 0);
	static bool first=true; // O.Bressand not beautifull but OK

	for (int i = 0; i < nTypes; i++) {
		TypeParticle * tp = Global::reference.find(i);
		TypeAtom *ta = dynamic_cast<TypeAtom*>(tp);
		if (ta) {
			if (comm->getRank() == 0 && first){
				std::cerr << "ParticleWriterHercule::writeTypes: Atom Type " << i << " name=" << tp->getName()
						<< " mass=" << tp->getMass() << " AtomicNumber=" << ta->getAtomicNumber()
						<< std::endl;
			}
			nbZ[i] = ta->getAtomicNumber();
		} else {
			if (comm->getRank() == 0 && first){
				std::cerr << "ParticleWriterHercule::writeTypes: Type " << i << " name=" << tp->getName() << " mass="
						<< tp->getMass() << std::endl;
			}
		}
	}
	first=false;

	for (uint p = 0; p < size; ++p) {
		typesval[p] = nbZ[(static_cast<int>(type[p]))];
	}

	HIc_Obj obj_nbZ = _obj_particleSet.getAttr("elements").createWithSupport(
			"nbZ", "GrandeurScalaireEntier");
	obj_nbZ.createAttr("val", "int_2").setVal(typesval, size);

	free(typesval);
}


/// @brief [ECom] ParticleWriterHercule::writeVelocities
void ParticleWriterHercule::writeVelocities(HIc_Ctx& _ctx,
		HIc_Obj& _obj_particleSet) {

	Array<vec3<double> >& v = particles->v;
	uint size = v.size();

	//	if (comm->getRank() == 0)
	//std::cout << "ParticleWriterHercule::writeVelocities " <<  size << std::endl;
	double *velocityval = (double*) malloc(sizeof(double) * size * 3);
	for (uint p = 0; p < size; ++p) {
		velocityval[p * 3 + 0] = v[p].x* 10.0; // TODO OB: conversion en unite special A
		velocityval[p * 3 + 1] = v[p].y* 10.0; // TODO OB: conversion en unite special A
		velocityval[p * 3 + 2] = v[p].z* 10.0; // TODO OB: conversion en unite special A
	}

	HIc_Obj nodes = _obj_particleSet.getAttr("noeuds");
	if (!nodes.isNull()) {
		HIc_Obj vitesse = nodes.createWithSupport("v", "Vitesse", nodes);
		vitesse.createAttr("val", "float_8[3]").setVal(velocityval, size * 3);
	} else {
		std::cout << "nodes.isNull()" << std::endl;
	}

	free(velocityval);
}


/// @brief [ECom] ParticleWriterHercule::writeEints
void ParticleWriterHercule::writeEints(HIc_Ctx& _ctx,
		HIc_Obj& _obj_particleSet) {

	Array<double>& ei = particles->ei;
	uint size = ei.size();

	//	if (comm->getRank() == 0)
	//std::cout << "ParticleWriterHercule::writeVelocities " <<  size << std::endl;
	double *eintval = (double*) malloc(sizeof(double) * size);
	for (uint p = 0; p < size; ++p) {
		eintval[p] = ei[p];
	}

	HIc_Obj nodes = _obj_particleSet.getAttr("noeuds");
	if (!nodes.isNull()) {
	  HIc_Obj eint = nodes.createWithSupport("ei", "GrandeurScalaire", nodes);
		eint.createAttr("val", "float_8").setVal(eintval, size);
	} else {
		std::cout << "nodes.isNull()" << std::endl;
	}

	free(eintval);

}


/// @brief [ECom] ParticleWriterHercule::writeProgresses
void ParticleWriterHercule::writeProgresses(HIc_Ctx& _ctx,
		HIc_Obj& _obj_particleSet) {

	Array<double>& progress = particles->progress;
	uint size = progress.size();

	//	if (comm->getRank() == 0)
	//std::cout << "ParticleWriterHercule::writeVelocities " <<  size << std::endl;
	double *progval = (double*) malloc(sizeof(double) * size);
	for (uint p = 0; p < size; ++p) {
		progval[p] = progress[p];
	}

	HIc_Obj nodes = _obj_particleSet.getAttr("noeuds");
	if (!nodes.isNull()) {
	  HIc_Obj prog = nodes.createWithSupport("progress", "GrandeurScalaire", nodes);
		prog.createAttr("val", "float_8").setVal(progval, size);
	} else {
		std::cout << "nodes.isNull()" << std::endl;
	}

	free(progval);

}

#endif
