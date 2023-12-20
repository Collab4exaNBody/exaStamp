/// @file
/// @brief Implementations for the ReaderHerculeDump
// ECom = to be commented or deleted by Estelle

#ifdef __use_lib_hercule

#include "io/ReaderHerculeDump.hpp"

#include "parallel/commManager.hpp"
#include "parallel/mympi.hpp"

#include "utils/stampUnits.hpp"

#include "HIc.h"
#include "HIc_MPI.h"


/// @brief [ECom] ReaderHerculeDump::openBase
void ReaderHerculeDump::openBase() {
	if (!m_base.isNull())
		return;

	if (isReader(-1)){
		m_base = HIc_Base(m_api, "myBase", "");

		m_base.setItemConf("mode", "read");
		m_base.setItemConf("read_dir", m_rep.c_str());
		m_base.setItemConf("bd_name", "HProt");
		m_base.setItemConf("access_mode", "sequential_access");
		m_base.setItemConf("multi_period", 1);
		m_base.open();
		vector<float_8> stepVector;
		m_base.getTimeList(stepVector);
		m_steps.alloc(stepVector.size());
		for(u_int i=0;i<stepVector.size();i++){
			m_steps[i]=stepVector[i];
		}
		vector<int_8> nbDomVector;
		m_base.getNbDomains(nbDomVector);
		m_nbDoms.alloc(nbDomVector.size());
		for(u_int i=0;i<nbDomVector.size();i++){
			m_nbDoms[i]=nbDomVector[i];
		}
		uint nbStep=m_steps.size();
		comm->broadcast(nbStep,masterIo(-1));
	}else{
		uint nbStep=0;
		comm->broadcast(nbStep,masterIo(-1));
		m_steps.alloc(nbStep);
		m_nbDoms.alloc(nbStep);
	}
		if (m_modeInitHercule != 0){
	comm->broadcast(m_steps,masterIo(-1));
	comm->broadcast(m_nbDoms,masterIo(-1));
		}
}


/// @brief [ECom] ReaderHerculeDump::nbStep
int ReaderHerculeDump::nbStep(void){
	return m_steps.size();
}


/// @brief [ECom] ReaderHerculeDump::steps
vector<int> ReaderHerculeDump::steps(){
	vector<int> ret(m_steps.size());
	for(u_int i=0;i<m_steps.size();i++){
		ret[i]=(int)m_steps[i];
	}
	return ret;
}


/// @brief [ECom] ReaderHerculeDump::nbSSDom
int ReaderHerculeDump::nbSSDom(int _step){
	for(u_int i=0;i<m_steps.size();i++){
		if ((int)m_steps[i] == _step){
			return m_nbDoms[i];
		}
	}
	std::cerr << "Can't find " << _step << std::endl;
	return 0;
}


/// @brief [ECom] ReaderHerculeDump::closeBase
void ReaderHerculeDump::closeBase() {
	if (isReader(-1)){
		m_base.close();
		m_base = HIc_Base();
	}
}


/// @brief [ECom] ReaderHerculeDump::readMetaInfo
void ReaderHerculeDump::readMetaInfo(){
	if (comm->getRank() == 0){
		// intialisation de la base en lecture du sous-domaine 0
		// lecture des infos de metaInfo
		// diffusion des infos
	}else{
		// reception des metaInfos
	}
}


/// @brief [ECom] ReaderHerculeDump::openStep
void ReaderHerculeDump::openStep(uint _step, int _numSDom) {
	if (_numSDom == -1){
		_numSDom=comm->getRank();
	}
	if (isReader(_numSDom)){
		m_ctx[_numSDom] = m_base.getCtxPar((double) _step, _numSDom);
		m_ctx[_numSDom].setItemConf("mode", "read");

		m_ctx[_numSDom].open();
	}

}


/// @brief [ECom] ReaderHerculeDump::closeStep
void ReaderHerculeDump::closeStep(int _numSDom) {
	if (isReader(_numSDom)){
		m_ctx[_numSDom].close();
		m_ctx[_numSDom] = HIc_Ctx();
	}
}


/// @brief [ECom] ReaderHerculeDump::readHeader
void ReaderHerculeDump::readHeader(int _numSDom, Array<LegacyHeaderIOStruct>& _header){
	HIc_Obj root;
	LegacyHeaderIOStruct& header=_header[0];

	if (isReader(_numSDom)){
		root = m_ctx[_numSDom].getRoot();

		root.get("global_xmin").getVal(header.xmin);
		root.get("global_xmax").getVal(header.xmax);
		root.get("global_ymin").getVal(header.ymin);
		root.get("global_ymax").getVal(header.ymax);
		root.get("global_zmin").getVal(header.zmin);
		root.get("global_zmax").getVal(header.zmax);

		root.get("time").getVal(header.time);
		root.get("CPUtime").getVal(header.CPUtime);
		root.get("totalEnergy").getVal(header.totalEnergy);
		root.get("potentialEnergy").getVal(header.potentialEnergy);
		root.get("kineticEnergy").getVal(header.kineticEnergy);
		root.get("internalEnergy").getVal(header.internalEnergy);
		root.get("rotationalEnergy").getVal(header.rotationalEnergy);
		root.get("particlesTotalNumber").getVal(header.particlesTotalNumber);
		root.get("iterationNumber").getVal(header.iterationNumber);
		root.get("domainLocalNumber").getVal(header.domainNumber);
		//root.get("domainNumber").getVal(header.domainNumber,header.domainNumber);
	}
	if (m_modeInitHercule != 0){
		comm->broadcast(header.xmin,masterIo(_numSDom));
		comm->broadcast(header.ymin,masterIo(_numSDom));
		comm->broadcast(header.zmin,masterIo(_numSDom));
		comm->broadcast(header.xmax,masterIo(_numSDom));
		comm->broadcast(header.ymax,masterIo(_numSDom));
		comm->broadcast(header.zmax,masterIo(_numSDom));

		comm->broadcast(header.time,masterIo(_numSDom));
		comm->broadcast(header.CPUtime,masterIo(_numSDom));
		comm->broadcast(header.totalEnergy,masterIo(_numSDom));
		comm->broadcast(header.potentialEnergy,masterIo(_numSDom));
		comm->broadcast(header.kineticEnergy,masterIo(_numSDom));
		comm->broadcast(header.internalEnergy,masterIo(_numSDom));
		comm->broadcast(header.rotationalEnergy,masterIo(_numSDom));
		comm->broadcast(header.particlesTotalNumber,masterIo(_numSDom));
		comm->broadcast(header.iterationNumber,masterIo(_numSDom));
		comm->broadcast(header.domainNumber,masterIo(_numSDom));
	}
}


/// @brief [ECom] ReaderHerculeDump::isReader
bool ReaderHerculeDump::isReader(int _numSDom){
	switch(m_modeInitHercule){
	default:
	case 0:
		return true;
	case 1:
	case 2:
		return comm->getRank() == 0;
	case 3:
		if (_numSDom == -1){
			return true;
		}
		return comm->getRank() == (_numSDom%comm->getNumberOfNodes());
	}
}


/// @brief [ECom] ReaderHerculeDump::masterIo
int ReaderHerculeDump::masterIo(int _numSDom){
	switch(m_modeInitHercule){
	case 0:
		return 0;
	default:
		if (_numSDom == -1){
			return -1;
		}
		return _numSDom % comm->getNumberOfNodes();
	}
}


/// @brief [ECom] ReaderHerculeDump::readParticles
void ReaderHerculeDump::readParticles(int _numSDom, HerculeParticleIODumpStruct& _particles){
	// il serait bon de recuperer le nombre de cellule, leur positionnement et le nombre de particules dans chaque cellule.
	// cela permettrait de transmettre le paquet de cellule au seul processeur pouvant recevoir une partie des particule de cette cellule.

	int_8 nbCell=0;
	int_8 nbParticle=0;
	HIc_Obj root;
	if (isReader(_numSDom) && m_modeInitHercule != 0){
		std::cerr << "le proc " << comm->getRank() << " lit le sous-domaine " << _numSDom << std::endl;
	}
	if (isReader(_numSDom)){
		root = m_ctx[_numSDom].getRoot();
		root.get("nbCell").getVal(nbCell);
		root.get("nbParticle").getVal(nbParticle);
	}
	if (m_modeInitHercule != 0){
		comm->broadcast(nbCell,masterIo(_numSDom));
		comm->broadcast(nbParticle,masterIo(_numSDom));
	}
	_particles.init(nbCell,nbParticle);

	const auto meterPerSecond = SI_Units_base::meter/SI_Units_base::second;

	if (nbCell){
		if (isReader(_numSDom)){
			root.get("cells_xmin").getVal(_particles.cells_xmin.data(),nbCell);
			root.get("cells_ymin").getVal(_particles.cells_ymin.data(),nbCell);
			root.get("cells_zmin").getVal(_particles.cells_zmin.data(),nbCell);
			root.get("cells_xmax").getVal(_particles.cells_xmax.data(),nbCell);
			root.get("cells_ymax").getVal(_particles.cells_ymax.data(),nbCell);
			root.get("cells_zmax").getVal(_particles.cells_zmax.data(),nbCell);
			root.get("cells_nbParticles").getVal(_particles.cells_nbParticles.data(),nbCell);
		}
		if (m_modeInitHercule != 0){
			//comm->broadcast(_particles.cells_xmin,masterIo(_numSDom));
			//comm->broadcast(_particles.cells_ymin,masterIo(_numSDom));
			//comm->broadcast(_particles.cells_zmin,masterIo(_numSDom));
			//comm->broadcast(_particles.cells_xmax,masterIo(_numSDom));
			//comm->broadcast(_particles.cells_ymax,masterIo(_numSDom));
			//comm->broadcast(_particles.cells_zmax,masterIo(_numSDom));
			//comm->broadcast(_particles.cells_nbParticles,masterIo(_numSDom));
		}
	}
	if (nbParticle){
		if (isReader(_numSDom)){
			//HIc_Obj particles_int = root.get("particles_int");
			//if (particles_int.isNull()){
			//}else{
			root.get("particles_coordinatesX").getVal(_particles.particles_coordinates[0].data(),nbParticle);
			root.get("particles_coordinatesY").getVal(_particles.particles_coordinates[1].data(),nbParticle);
			root.get("particles_coordinatesZ").getVal(_particles.particles_coordinates[2].data(),nbParticle);
			root.get("particles_johnDoeX").getVal(_particles.particles_johnDoe[0].data(),nbParticle);
			root.get("particles_johnDoeY").getVal(_particles.particles_johnDoe[1].data(),nbParticle);
			root.get("particles_johnDoeZ").getVal(_particles.particles_johnDoe[2].data(),nbParticle);
			root.get("particles_velocityX").getVal(_particles.particles_velocity[0].data(),nbParticle);
			root.get("particles_velocityY").getVal(_particles.particles_velocity[1].data(),nbParticle);
			root.get("particles_velocityZ").getVal(_particles.particles_velocity[2].data(),nbParticle);
			root.get("particles_quaternion0").getVal(_particles.particles_quaternion[0].data(),nbParticle);
			root.get("particles_quaternion1").getVal(_particles.particles_quaternion[1].data(),nbParticle);
			root.get("particles_quaternion2").getVal(_particles.particles_quaternion[2].data(),nbParticle);
			root.get("particles_quaternion3").getVal(_particles.particles_quaternion[3].data(),nbParticle);
			root.get("particles_momentumAngular0").getVal(_particles.particles_momentumAngular[0].data(),nbParticle);
			root.get("particles_momentumAngular1").getVal(_particles.particles_momentumAngular[1].data(),nbParticle);
			root.get("particles_momentumAngular2").getVal(_particles.particles_momentumAngular[2].data(),nbParticle);
			root.get("particles_oriantationX").getVal(_particles.particles_oriantation[0].data(),nbParticle);
			root.get("particles_oriantationY").getVal(_particles.particles_oriantation[0].data(),nbParticle);
			root.get("particles_oriantationZ").getVal(_particles.particles_oriantation[0].data(),nbParticle);
			root.get("particles_particleType").getVal(_particles.particles_type.data(),nbParticle*1);
			root.get("particles_particleID").getVal(_particles.particles_iD.data(),nbParticle*1);
			//}
		}
		if (m_modeInitHercule != 0){
			comm->broadcast(_particles.particles_coordinates[0],masterIo(_numSDom));
			comm->broadcast(_particles.particles_coordinates[1],masterIo(_numSDom));
			comm->broadcast(_particles.particles_coordinates[2],masterIo(_numSDom));
			comm->broadcast(_particles.particles_johnDoe[0],masterIo(_numSDom));
			comm->broadcast(_particles.particles_johnDoe[1],masterIo(_numSDom));
			comm->broadcast(_particles.particles_johnDoe[2],masterIo(_numSDom));
			comm->broadcast(_particles.particles_velocity[0],masterIo(_numSDom));
			comm->broadcast(_particles.particles_velocity[1],masterIo(_numSDom));
			comm->broadcast(_particles.particles_velocity[2],masterIo(_numSDom));
			comm->broadcast(_particles.particles_quaternion[0],masterIo(_numSDom));
			comm->broadcast(_particles.particles_quaternion[1],masterIo(_numSDom));
			comm->broadcast(_particles.particles_quaternion[2],masterIo(_numSDom));
			comm->broadcast(_particles.particles_quaternion[3],masterIo(_numSDom));
			comm->broadcast(_particles.particles_momentumAngular[0],masterIo(_numSDom));
			comm->broadcast(_particles.particles_momentumAngular[1],masterIo(_numSDom));
			comm->broadcast(_particles.particles_momentumAngular[2],masterIo(_numSDom));
			comm->broadcast(_particles.particles_oriantation[0],masterIo(_numSDom));
			comm->broadcast(_particles.particles_oriantation[1],masterIo(_numSDom));
			comm->broadcast(_particles.particles_oriantation[2],masterIo(_numSDom));
			comm->broadcast(_particles.particles_type,masterIo(_numSDom));
			comm->broadcast(_particles.particles_iD,masterIo(_numSDom));
		}
	}
#define HPROT_CONVERT_UNIT
	// pb coherence avec input
#ifdef HPROT_CONVERT_UNIT
	// changement d'unite
	for(int i=0;i<nbParticle;i++){
		_particles.particles_coordinates[0][i] = convert(_particles.particles_coordinates[0][i], SI_Units_base::meter, Stamp_Units::length);
		_particles.particles_coordinates[1][i] = convert(_particles.particles_coordinates[1][i], SI_Units_base::meter, Stamp_Units::length);
		_particles.particles_coordinates[2][i] = convert(_particles.particles_coordinates[2][i], SI_Units_base::meter, Stamp_Units::length);
		_particles.particles_velocity[0][i]= convert(_particles.particles_velocity[0][i], meterPerSecond, Stamp_Units::speed);
		_particles.particles_velocity[1][i]= convert(_particles.particles_velocity[1][i], meterPerSecond, Stamp_Units::speed);
		_particles.particles_velocity[2][i]= convert(_particles.particles_velocity[2][i], meterPerSecond, Stamp_Units::speed);
	}
#endif
}


/// @brief [ECom] ReaderHerculeDump::readParticles (for DPDE)
void ReaderHerculeDump::readParticles(int _numSDom, HerculeDPDEParticleIODumpStruct& _particles){
	// il serait bon de recuperer le nombre de cellule, leur positionnement et le nombre de particules dans chaque cellule.
	// cela permettrait de transmettre le paquet de cellule au seul processeur pouvant recevoir une partie des particule de cette cellule.

	int_8 nbCell=0;
	int_8 nbParticle=0;
	HIc_Obj root;
	if (isReader(_numSDom) && m_modeInitHercule != 0){
		std::cerr << "le proc " << comm->getRank() << " lit le sous-domaine " << _numSDom << std::endl;
	}
	if (isReader(_numSDom)){
		root = m_ctx[_numSDom].getRoot();
		root.get("nbCell").getVal(nbCell);
		root.get("nbParticle").getVal(nbParticle);
	}
	if (m_modeInitHercule != 0){
		comm->broadcast(nbCell,masterIo(_numSDom));
		comm->broadcast(nbParticle,masterIo(_numSDom));
	}
	_particles.init(nbCell,nbParticle);

	const auto meterPerSecond = SI_Units_base::meter/SI_Units_base::second;

	if (nbCell){
		if (isReader(_numSDom)){
			root.get("cells_xmin").getVal(_particles.cells_xmin.data(),nbCell);
			root.get("cells_ymin").getVal(_particles.cells_ymin.data(),nbCell);
			root.get("cells_zmin").getVal(_particles.cells_zmin.data(),nbCell);
			root.get("cells_xmax").getVal(_particles.cells_xmax.data(),nbCell);
			root.get("cells_ymax").getVal(_particles.cells_ymax.data(),nbCell);
			root.get("cells_zmax").getVal(_particles.cells_zmax.data(),nbCell);
			root.get("cells_nbParticles").getVal(_particles.cells_nbParticles.data(),nbCell);
		}
		if (m_modeInitHercule != 0){
			//comm->broadcast(_particles.cells_xmin,masterIo(_numSDom));
			//comm->broadcast(_particles.cells_ymin,masterIo(_numSDom));
			//comm->broadcast(_particles.cells_zmin,masterIo(_numSDom));
			//comm->broadcast(_particles.cells_xmax,masterIo(_numSDom));
			//comm->broadcast(_particles.cells_ymax,masterIo(_numSDom));
			//comm->broadcast(_particles.cells_zmax,masterIo(_numSDom));
			//comm->broadcast(_particles.cells_nbParticles,masterIo(_numSDom));
		}
	}
	if (nbParticle){
		if (isReader(_numSDom)){
			//HIc_Obj particles_int = root.get("particles_int");
			//if (particles_int.isNull()){
			//}else{
			root.get("particles_coordinatesX").getVal(_particles.particles_coordinates[0].data(),nbParticle);
			root.get("particles_coordinatesY").getVal(_particles.particles_coordinates[1].data(),nbParticle);
			root.get("particles_coordinatesZ").getVal(_particles.particles_coordinates[2].data(),nbParticle);
			root.get("particles_johnDoeX").getVal(_particles.particles_johnDoe[0].data(),nbParticle);
			root.get("particles_johnDoeY").getVal(_particles.particles_johnDoe[1].data(),nbParticle);
			root.get("particles_johnDoeZ").getVal(_particles.particles_johnDoe[2].data(),nbParticle);
			root.get("particles_velocityX").getVal(_particles.particles_velocity[0].data(),nbParticle);
			root.get("particles_velocityY").getVal(_particles.particles_velocity[1].data(),nbParticle);
			root.get("particles_velocityZ").getVal(_particles.particles_velocity[2].data(),nbParticle);
			root.get("particles_internalEnergy").getVal(_particles.particles_internalEnergy.data(),nbParticle);
			root.get("particles_internalTemperature").getVal(_particles.particles_internalTemperature.data(),nbParticle);
			root.get("particles_progress").getVal(_particles.particles_progress.data(),nbParticle);
			root.get("particles_quaternion0").getVal(_particles.particles_quaternion[0].data(),nbParticle);
			root.get("particles_quaternion1").getVal(_particles.particles_quaternion[1].data(),nbParticle);
			root.get("particles_quaternion2").getVal(_particles.particles_quaternion[2].data(),nbParticle);
			root.get("particles_quaternion3").getVal(_particles.particles_quaternion[3].data(),nbParticle);
			root.get("particles_momentumAngular0").getVal(_particles.particles_momentumAngular[0].data(),nbParticle);
			root.get("particles_momentumAngular1").getVal(_particles.particles_momentumAngular[1].data(),nbParticle);
			root.get("particles_momentumAngular2").getVal(_particles.particles_momentumAngular[2].data(),nbParticle);
			root.get("particles_oriantationX").getVal(_particles.particles_oriantation[0].data(),nbParticle);
			root.get("particles_oriantationY").getVal(_particles.particles_oriantation[0].data(),nbParticle);
			root.get("particles_oriantationZ").getVal(_particles.particles_oriantation[0].data(),nbParticle);
			root.get("particles_particleType").getVal(_particles.particles_type.data(),nbParticle*1);
			root.get("particles_particleID").getVal(_particles.particles_iD.data(),nbParticle*1);
			//}
		}
		if (m_modeInitHercule != 0){
			comm->broadcast(_particles.particles_coordinates[0],masterIo(_numSDom));
			comm->broadcast(_particles.particles_coordinates[1],masterIo(_numSDom));
			comm->broadcast(_particles.particles_coordinates[2],masterIo(_numSDom));
			comm->broadcast(_particles.particles_johnDoe[0],masterIo(_numSDom));
			comm->broadcast(_particles.particles_johnDoe[1],masterIo(_numSDom));
			comm->broadcast(_particles.particles_johnDoe[2],masterIo(_numSDom));
			comm->broadcast(_particles.particles_velocity[0],masterIo(_numSDom));
			comm->broadcast(_particles.particles_velocity[1],masterIo(_numSDom));
			comm->broadcast(_particles.particles_velocity[2],masterIo(_numSDom));
			comm->broadcast(_particles.particles_internalEnergy,masterIo(_numSDom));
			comm->broadcast(_particles.particles_internalTemperature,masterIo(_numSDom));
			comm->broadcast(_particles.particles_progress,masterIo(_numSDom));
			comm->broadcast(_particles.particles_quaternion[0],masterIo(_numSDom));
			comm->broadcast(_particles.particles_quaternion[1],masterIo(_numSDom));
			comm->broadcast(_particles.particles_quaternion[2],masterIo(_numSDom));
			comm->broadcast(_particles.particles_quaternion[3],masterIo(_numSDom));
			comm->broadcast(_particles.particles_momentumAngular[0],masterIo(_numSDom));
			comm->broadcast(_particles.particles_momentumAngular[1],masterIo(_numSDom));
			comm->broadcast(_particles.particles_momentumAngular[2],masterIo(_numSDom));
			comm->broadcast(_particles.particles_oriantation[0],masterIo(_numSDom));
			comm->broadcast(_particles.particles_oriantation[1],masterIo(_numSDom));
			comm->broadcast(_particles.particles_oriantation[2],masterIo(_numSDom));
			comm->broadcast(_particles.particles_type,masterIo(_numSDom));
			comm->broadcast(_particles.particles_iD,masterIo(_numSDom));
		}
	}
#define HPROT_CONVERT_UNIT
	// pb coherence avec input
#ifdef HPROT_CONVERT_UNIT
	// changement d'unite
	for(int i=0;i<nbParticle;i++){
		_particles.particles_coordinates[0][i] = convert(_particles.particles_coordinates[0][i], SI_Units_base::meter, Stamp_Units::length);
		_particles.particles_coordinates[1][i] = convert(_particles.particles_coordinates[1][i], SI_Units_base::meter, Stamp_Units::length);
		_particles.particles_coordinates[2][i] = convert(_particles.particles_coordinates[2][i], SI_Units_base::meter, Stamp_Units::length);
		_particles.particles_velocity[0][i] = convert(_particles.particles_velocity[0][i], meterPerSecond, Stamp_Units::speed);
		_particles.particles_velocity[1][i] = convert(_particles.particles_velocity[1][i], meterPerSecond, Stamp_Units::speed);
		_particles.particles_velocity[2][i] = convert(_particles.particles_velocity[2][i], meterPerSecond, Stamp_Units::speed);
		_particles.particles_internalEnergy[i] = convert(_particles.particles_internalEnergy[i], SI_Units_base::joule, Stamp_Units::energy);
	}
#endif
}


#endif
