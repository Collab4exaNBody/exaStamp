/// @file
/// @brief Implementations for the WriterHerculeDump
// ECom = to be commented or deleted by Estelle

#ifdef __use_lib_hercule

#include "io/WriterHerculeDump.hpp"

#include "parallel/commManager.hpp"
#include "parallel/mympi.hpp"

#include "HIc.h"
#include "HIc_MPI.h"


/// @brief [ECom] WriterHerculeDump::openBase
void WriterHerculeDump::openBase() {
	if (!m_base.isNull())
		return;

	m_base = HIc_Base(m_api, "myBase", "parallele");

	m_base.setItemConf("mode", "write");
	m_base.setItemConf("write_dir", m_rep.c_str());
	m_base.setItemConf("read_dir", m_rep.c_str());
	m_base.setItemConf("bd_name", "HProt");
	m_base.setItemConf("multi_period", 1);
	m_base.open();
}


/// @brief [ECom] WriterHerculeDump::closeBase
void WriterHerculeDump::closeBase() {
	m_base.close();
	m_base = HIc_Base();
}


/// @brief [ECom] WriterHerculeDump::openStep
void WriterHerculeDump::openStep(uint _step, int _numSSDom) {
	int numSSDom=_numSSDom;
	if (numSSDom == -1){
		numSSDom=comm->getRank();
	}
	m_ctx = m_base.createCtxPar((double) _step, comm->getNumberOfNodes(),
			numSSDom);

	m_ctx.setItemConf("mode", "write");

	m_ctx.open();

}


/// @brief [ECom] WriterHerculeDump::closeStep
void WriterHerculeDump::closeStep() {
	m_ctx.close();
	m_ctx = HIc_Ctx();
}


/// @brief [ECom] WriterHerculeDump::writeHeader
void WriterHerculeDump::writeHeader(Array<LegacyHeaderIOStruct>& _header){
	HIc_Obj root = m_ctx.getRoot();
	LegacyHeaderIOStruct& header=_header[0];

	root.createWithSupport("global_xmin","float_8").setVal(header.xmin);
	root.createWithSupport("global_xmax","float_8").setVal(header.xmax);
	root.createWithSupport("global_ymin","float_8").setVal(header.ymin);
	root.createWithSupport("global_ymax","float_8").setVal(header.ymax);
	root.createWithSupport("global_zmin","float_8").setVal(header.zmin);
	root.createWithSupport("global_zmax","float_8").setVal(header.zmax);

	root.createWithSupport("time","float_8").setVal(header.time);
	root.createWithSupport("CPUtime","float_8").setVal(header.CPUtime);
	root.createWithSupport("totalEnergy","float_8").setVal(header.totalEnergy);
	root.createWithSupport("potentialEnergy","float_8").setVal(header.potentialEnergy);
	root.createWithSupport("kineticEnergy","float_8").setVal(header.kineticEnergy);
	root.createWithSupport("internalEnergy","float_8").setVal(header.internalEnergy);
	root.createWithSupport("rotationalEnergy","float_8").setVal(header.rotationalEnergy);
	root.createWithSupport("particlesTotalNumber","int_8").setVal(header.particlesTotalNumber);
	root.createWithSupport("iterationNumber","int_8").setVal(header.iterationNumber);
	root.createWithSupport("domainLocalNumber","int_8").setVal(header.domainNumber);
	//root.createWithSupport("domainNumber","float_8[domainLocalNumber]").setVal(header.domainNumber,header.domainNumber);
}


/// @brief [ECom] WriterHerculeDump::writeParticles
void WriterHerculeDump::writeParticles(HerculeParticleIODumpStruct& _particles){
	// il serait bon de recuperer le nombre de cellule, leur positionnement et le nombre de particules dans chaque cellule.
	// cela permettrait de transmettre le paquet de cellule au seul processeur pouvant recevoir une partie ces particule de cette cellule.

	HIc_Obj root = m_ctx.getRoot();

	int_8 nbCell = _particles.nbCell;

	char tname_float_8[256];
	char tname_int_8[256];

	sprintf(tname_float_8, "float_8[%ld]", nbCell);
	sprintf(tname_int_8, "int_8[%ld]", nbCell);

	root.createWithSupport("nbCell","int_8").setVal(nbCell);
	root.createWithSupport("cells_xmin",tname_float_8).setVal(_particles.cells_xmin.data(),nbCell);
	root.createWithSupport("cells_ymin",tname_float_8).setVal(_particles.cells_xmin.data(),nbCell);
	root.createWithSupport("cells_zmin",tname_float_8).setVal(_particles.cells_xmin.data(),nbCell);
	root.createWithSupport("cells_xmax",tname_float_8).setVal(_particles.cells_xmin.data(),nbCell);
	root.createWithSupport("cells_ymax",tname_float_8).setVal(_particles.cells_xmin.data(),nbCell);
	root.createWithSupport("cells_zmax",tname_float_8).setVal(_particles.cells_xmin.data(),nbCell);
	root.createWithSupport("cells_nbParticles",tname_int_8).setVal(_particles.cells_nbParticles.data(),nbCell);

	int_8 nbParticle = _particles.nbParticle;

	sprintf(tname_float_8, "float_8[%ld]", nbParticle);
	sprintf(tname_int_8, "int_8[%ld]", nbParticle);

	root.createWithSupport("nbParticle","int_8").setVal(nbParticle);
	root.createWithSupport("particles_coordinatesX",tname_float_8).setVal(_particles.particles_coordinates[0].data(),nbParticle);
	root.createWithSupport("particles_coordinatesY",tname_float_8).setVal(_particles.particles_coordinates[1].data(),nbParticle);
	root.createWithSupport("particles_coordinatesZ",tname_float_8).setVal(_particles.particles_coordinates[2].data(),nbParticle);
	root.createWithSupport("particles_johnDoeX",tname_float_8).setVal(_particles.particles_johnDoe[0].data(),nbParticle);
	root.createWithSupport("particles_johnDoeY",tname_float_8).setVal(_particles.particles_johnDoe[1].data(),nbParticle);
	root.createWithSupport("particles_johnDoeZ",tname_float_8).setVal(_particles.particles_johnDoe[2].data(),nbParticle);
	root.createWithSupport("particles_velocityX",tname_float_8).setVal(_particles.particles_velocity[0].data(),nbParticle);
	root.createWithSupport("particles_velocityY",tname_float_8).setVal(_particles.particles_velocity[1].data(),nbParticle);
	root.createWithSupport("particles_velocityZ",tname_float_8).setVal(_particles.particles_velocity[2].data(),nbParticle);
	root.createWithSupport("particles_quaternion0",tname_float_8).setVal(_particles.particles_quaternion[0].data(),nbParticle);
	root.createWithSupport("particles_quaternion1",tname_float_8).setVal(_particles.particles_quaternion[1].data(),nbParticle);
	root.createWithSupport("particles_quaternion2",tname_float_8).setVal(_particles.particles_quaternion[2].data(),nbParticle);
	root.createWithSupport("particles_quaternion3",tname_float_8).setVal(_particles.particles_quaternion[3].data(),nbParticle);
	root.createWithSupport("particles_momentumAngular0",tname_float_8).setVal(_particles.particles_momentumAngular[0].data(),nbParticle);
	root.createWithSupport("particles_momentumAngular1",tname_float_8).setVal(_particles.particles_momentumAngular[1].data(),nbParticle);
	root.createWithSupport("particles_momentumAngular2",tname_float_8).setVal(_particles.particles_momentumAngular[2].data(),nbParticle);
	root.createWithSupport("particles_oriantationX",tname_float_8).setVal(_particles.particles_oriantation[0].data(),nbParticle);
	root.createWithSupport("particles_oriantationY",tname_float_8).setVal(_particles.particles_oriantation[0].data(),nbParticle);
	root.createWithSupport("particles_oriantationZ",tname_float_8).setVal(_particles.particles_oriantation[0].data(),nbParticle);
	root.createWithSupport("particles_particleType",tname_int_8).setVal(_particles.particles_type.data(),nbParticle*1);
	root.createWithSupport("particles_particleID",tname_int_8).setVal(_particles.particles_iD.data(),nbParticle*1);
}


/// @brief [ECom] WriterHerculeDump::writeParticles
void WriterHerculeDump::writeParticles(HerculeDPDEParticleIODumpStruct& _particles){
	// il serait bon de recuperer le nombre de cellule, leur positionnement et le nombre de particules dans chaque cellule.
	// cela permettrait de transmettre le paquet de cellule au seul processeur pouvant recevoir une partie ces particule de cette cellule.

	HIc_Obj root = m_ctx.getRoot();

	int_8 nbCell = _particles.nbCell;

	char tname_float_8[256];
	char tname_int_8[256];

	sprintf(tname_float_8, "float_8[%ld]", nbCell);
	sprintf(tname_int_8, "int_8[%ld]", nbCell);

	root.createWithSupport("nbCell","int_8").setVal(nbCell);
	root.createWithSupport("cells_xmin",tname_float_8).setVal(_particles.cells_xmin.data(),nbCell);
	root.createWithSupport("cells_ymin",tname_float_8).setVal(_particles.cells_xmin.data(),nbCell);
	root.createWithSupport("cells_zmin",tname_float_8).setVal(_particles.cells_xmin.data(),nbCell);
	root.createWithSupport("cells_xmax",tname_float_8).setVal(_particles.cells_xmin.data(),nbCell);
	root.createWithSupport("cells_ymax",tname_float_8).setVal(_particles.cells_xmin.data(),nbCell);
	root.createWithSupport("cells_zmax",tname_float_8).setVal(_particles.cells_xmin.data(),nbCell);
	root.createWithSupport("cells_nbParticles",tname_int_8).setVal(_particles.cells_nbParticles.data(),nbCell);

	int_8 nbParticle = _particles.nbParticle;

	sprintf(tname_float_8, "float_8[%ld]", nbParticle);
	sprintf(tname_int_8, "int_8[%ld]", nbParticle);

	root.createWithSupport("nbParticle","int_8").setVal(nbParticle);
	root.createWithSupport("particles_coordinatesX",tname_float_8).setVal(_particles.particles_coordinates[0].data(),nbParticle);
	root.createWithSupport("particles_coordinatesY",tname_float_8).setVal(_particles.particles_coordinates[1].data(),nbParticle);
	root.createWithSupport("particles_coordinatesZ",tname_float_8).setVal(_particles.particles_coordinates[2].data(),nbParticle);
	root.createWithSupport("particles_johnDoeX",tname_float_8).setVal(_particles.particles_johnDoe[0].data(),nbParticle);
	root.createWithSupport("particles_johnDoeY",tname_float_8).setVal(_particles.particles_johnDoe[1].data(),nbParticle);
	root.createWithSupport("particles_johnDoeZ",tname_float_8).setVal(_particles.particles_johnDoe[2].data(),nbParticle);
	root.createWithSupport("particles_velocityX",tname_float_8).setVal(_particles.particles_velocity[0].data(),nbParticle);
	root.createWithSupport("particles_velocityY",tname_float_8).setVal(_particles.particles_velocity[1].data(),nbParticle);
	root.createWithSupport("particles_velocityZ",tname_float_8).setVal(_particles.particles_velocity[2].data(),nbParticle);
	root.createWithSupport("particles_internalEnergy",tname_float_8).setVal(_particles.particles_internalEnergy.data(),nbParticle);
	root.createWithSupport("particles_internalTemperature",tname_float_8).setVal(_particles.particles_internalTemperature.data(),nbParticle);
	root.createWithSupport("particles_progress",tname_float_8).setVal(_particles.particles_progress.data(),nbParticle);
	root.createWithSupport("particles_quaternion0",tname_float_8).setVal(_particles.particles_quaternion[0].data(),nbParticle);
	root.createWithSupport("particles_quaternion1",tname_float_8).setVal(_particles.particles_quaternion[1].data(),nbParticle);
	root.createWithSupport("particles_quaternion2",tname_float_8).setVal(_particles.particles_quaternion[2].data(),nbParticle);
	root.createWithSupport("particles_quaternion3",tname_float_8).setVal(_particles.particles_quaternion[3].data(),nbParticle);
	root.createWithSupport("particles_momentumAngular0",tname_float_8).setVal(_particles.particles_momentumAngular[0].data(),nbParticle);
	root.createWithSupport("particles_momentumAngular1",tname_float_8).setVal(_particles.particles_momentumAngular[1].data(),nbParticle);
	root.createWithSupport("particles_momentumAngular2",tname_float_8).setVal(_particles.particles_momentumAngular[2].data(),nbParticle);
	root.createWithSupport("particles_oriantationX",tname_float_8).setVal(_particles.particles_oriantation[0].data(),nbParticle);
	root.createWithSupport("particles_oriantationY",tname_float_8).setVal(_particles.particles_oriantation[0].data(),nbParticle);
	root.createWithSupport("particles_oriantationZ",tname_float_8).setVal(_particles.particles_oriantation[0].data(),nbParticle);
	root.createWithSupport("particles_particleType",tname_int_8).setVal(_particles.particles_type.data(),nbParticle*1);
	root.createWithSupport("particles_particleID",tname_int_8).setVal(_particles.particles_iD.data(),nbParticle*1);
}


#endif
