/// @file 
/// @brief Implementation of node writing related methods


#include <iomanip>
#include <fstream>
#include <string>
#include <thread>
#include "omp.h"

#include "referenceMap.hpp"
//
#include "domain/domainInfo.hpp"
#include "domain/domainInterface.hpp"

#include "io/cellOutputVTK.hpp"
#include "io/particleOutputDAT.hpp"
#include "io/particleOutputHercule.hpp"
#include "io/particleOutputVTK.hpp"
#if __use_lib_hercule
#include "io/WriterHerculeDump.hpp"
#endif
#include "io/ParticleWriterLegacyDump.hpp"

#include "parallel/node.hpp"

#include "time/time.hpp"

#include "utils/stampUnits.hpp"


/// @brief Local function to print header in a flux
/// @param [in,out] f Print flux
static void writeHeader(std::ostream& f) {

  f << "# " << std::endl
    << "#       Step   Number of   LB State  Imbalance   Tot. Energy (eV/   Kin. Energy (eV/   "
    << "Pot. Energy (eV/   ";
  if (Global::reference.isDPD())
    f << "Int. Energy (eV/   ";
  if (Global::reference.isReactive())
    f << "Chem. Energy (eV/   ";
  f << "Tempera-  Pressure  Time (s/part  ";
      if(Global::reference.isTimers())
  f << " Atom/s/stp " << std::endl
    << "#             particles                         part)              part)              p"
    << "art)              ";
  f << "ture (K)   (Pa)      /stp/tsk)     " << std::endl << "# "
    << std::endl;

};



/// @brief Local function to print header in a flux
/// @param [in,out] f Print flux
static void writeHeader_amr(std::ostream& f) {

  f << "# " << std::endl
    << "#       Step   Number of   LB State  Imbalance   Tot. Energy (eV/   Kin. Energy (eV/   "
    << "Pot. Energy (eV/   ";
  f << "Tempera-  Time (s/part  ";
      if(Global::reference.isTimers())
  f << " Atom/s/stp    Mem Used" << std::endl
    << "#             particles                         part)              part)              part)              ";
  f << "ture (K)   /stp/tsk)                   (KB)  " << std::endl << "# "
    << std::endl;

};





/// @brief Local function to print the state of the system in a flux
/// @param [in,out] f Print flux
/// @param [in] step Current computation step
/// @param [in] totalNumberOfParticles Total number of particles
/// @param [in] state Load balancer state
/// @param [in] imbalance Imbalance rate
/// @param [in] energies Kinetic, potential, internal, chemical and total energies
/// @param [in] temperature Temperature
/// @param [in] pressure Pressure
/// @param [in] gt CPU time per iteration * nb of cores / nb particle per iteration
static void write_line(std::ostream& f, uint step,
		       uint64_t totalNumberOfParticles, LBS::State state, double imbalance,
		       const Array<double>& energies, double temperature, double pressure, double gt, int numberOfCores) {

  using namespace std;

  f << "  " << setw(10) << step << "  " << setw(10) << totalNumberOfParticles
    << "  " << setw(9) << LBS::str(state) << "  " << fixed
    << setw(9) << setprecision(3) << imbalance << "  " << scientific
    << setw(17) << setprecision(10) << energies[4] << "  " << scientific
    << setw(17) << setprecision(10) << energies[0] << "  " << scientific
    << setw(17) << setprecision(10) << energies[1] << "  ";
  if (Global::reference.isDPD())
    f << scientific << setw(17) << setprecision(10) << energies[2] << " ";
  if (Global::reference.isReactive())
    f << scientific << setw(17) << setprecision(10) << energies[3] << " ";
  f << fixed << setw(7) << setprecision(3) << temperature << "  " << scientific
    << setw(9) << setprecision(3) << pressure << " " << scientific
    << setw(13) << setprecision(3) << gt*numberOfCores << "  ";
       if(Global::reference.isTimers())
       {
         if(gt>0) f  << setw(11) <<  1./gt;     
         else f  << setw(11) <<  0.; 
       }
  
  f  << endl;


}


/// @brief Local function to print the state of the system in a flux
/// @param [in,out] f Print flux
/// @param [in] step Current computation step
/// @param [in] totalNumberOfParticles Total number of particles
/// @param [in] state Load balancer state
/// @param [in] imbalance Imbalance rate
/// @param [in] energies Kinetic, potential, internal, chemical and total energies
/// @param [in] temperature Temperature
/// @param [in] pressure Pressure
/// @param [in] gt CPU time per iteration * nb of cores / nb particle per iteration
static void write_line_amr(std::ostream& f, uint step,
		       uint64_t totalNumberOfParticles, LBS::State state, double imbalance,
		       const Array<double>& energies, double temperature, double gt, int numberOfCores, double memU) {

  using namespace std;

  f << "  " << setw(10) << step << "  " << setw(10) << totalNumberOfParticles
    << "  " << setw(9) << LBS::str(state) << "  " << fixed
    << setw(9) << setprecision(3) << imbalance << "  " << scientific
    << setw(17) << setprecision(10) << energies[4] << "  " << scientific
    << setw(17) << setprecision(10) << energies[0] << "  " << scientific
    << setw(17) << setprecision(10) << energies[1] << "  ";
  if (Global::reference.isDPD())
    f << scientific << setw(17) << setprecision(10) << energies[2] << " ";
  if (Global::reference.isReactive())
    f << scientific << setw(17) << setprecision(10) << energies[3] << " ";
  f << fixed << setw(7) << setprecision(3) << temperature << "  " << scientific
    << setw(13) << setprecision(3) << gt*numberOfCores << "  ";
       if(Global::reference.isTimers())
       {
         if(gt>0) f  << setw(11) <<  1./gt;     
         else f  << setw(11) <<  0.; 
       }
  f<< setw(11) << setprecision(3) << memU << "  ";
  
  f  << endl;


}


/// @brief Local function to print a time in a flux
/// @param [in,out] f Print flux
/// @param [in] str Type of time
/// @param [in] d Number of days
/// @param [in] h Number of hours
/// @param [in] m Number of minutes
/// @param [in] s Number of seconds
static void write_time(std::ostream& f, const std::string& str, uint64_t d,
		       uint64_t h, uint64_t m, double s) {

  using namespace std;

  f << str << setw(2) << d << " d. " << setw(2) << h << " h. " << setw(2) << m
    << " m. " << fixed << setw(5) << setprecision(2) << s << " s."
    << endl;

}


/// @brief Local function to print final considerations in a flux
/// @param [in,out] f Print flux
/// @param [in] totalNumberOfParticles Total number of particles
/// @param [in] energies Kinetic, potential, internal, chemical and total energies
/// @param [in] temperature Temperature
/// @param [in] timeZero Initial kinetic energy, potential energy, internal energy, chenmical energy, total energy, temperature and number of particles
/// @param [in] memUsage Memory usage
/// @param [in] dt CPU time per iteration
/// @param [in] gt CPU time * nb of cores / nb particle per iteration
/// @param [in] days Array of days for different times
/// @param [in] hrs Array of hours for different times
/// @param [in] min Array of minutes for different times
/// @param [in] sec Array of secondes for different times
static void write_end(std::ostream& f, uint64_t totalNumberOfParticles,
		      const Array<double>& energies, const double temperature, const double* timeZero,
		      double memUsage, double dt, double gt, const uint64_t* days,
		      const uint64_t* hrs, const uint64_t* min, const double* sec) {

  using namespace std;

  f << "# " << endl << "# --- Run stats ------------------" << endl << "# "
    << endl << "#  particle gain/loss    : "
    << static_cast<int>(timeZero[6] - totalNumberOfParticles) << endl
    << "#  abs(1 - Etot/Etot0)   : " << scientific << setprecision(3)
    << auxAbs(1 - energies[4] / timeZero[4]) << "  " << endl
    << "#  abs(1 - T/T0)         : " << scientific << setprecision(3)
    << auxAbs(1 - temperature / timeZero[5]) << "  " << endl;
  
  f << "# " << "  " << endl << "#  max res set size (MB) : " << scientific
    << setprecision(3) << memUsage * 1.0e-03 << "  " << endl
    << "#  avg. part. throughput : " << scientific << setprecision(3)
    << totalNumberOfParticles / dt << "  " << endl
    << "#  avg. eff. time        : " << scientific << setprecision(3)
    << gt << endl << "# " << "  " << endl;
}


/// @brief Print n end lines
/// @param [in,out] flux Print flux, default=std::cout
/// @param [in] n Number of end lines, default=1
void Node::writeEndLine(std::ostream& flux, int n) {

  if (isMaster)
    for (int i = 0; i < n; ++i)
      flux << std::endl;

}



static bool first_step = true;
static int compt = 0;

/// @brief Write energy log in flux and energy log file
/// @param [in,out] flux Print flux
/// @param [in] step Current step (-1 if end)
/// @param [in] dt CPU time per iteration
/// @param [in] state Load balancer state
// with step = time->CurrentStep() and dt =  cpuTime/iterations
void Node::writeEnergies(std::ostream& flux, int step, double dt, LBS::State state) {

  uint64_t totalNumberOfParticles = reduceNumberOfParticles();
  uint64_t numNodes = (uint64_t) numberOfNodes();

  // First call setup
  static bool firstCall = true;
  static std::ofstream out;
  static int numberOfCores = -1;

  // To store energies and temperature at init
  static double timeZero[7];

  // Reduce total energy, compute grain time and temperature

  Array<double> energies = reduceEnergy();
  if(isMaster) for (unsigned int i = 0; i < energies.size(); i++)
    energies[i] /= ((double) totalNumberOfParticles);

  double temperature = convert(energies[4], Stamp_Units::energy, SI_Units_base::joule);
  temperature *= 2. / (3. * Constant::boltzmann);

  energies[0] = convert(energies[0], Stamp_Units::energy, SI_Units_base::electronVolt);
  energies[1] = convert(energies[1], Stamp_Units::energy, SI_Units_base::electronVolt);
  energies[2] = convert(energies[2], Stamp_Units::energy, SI_Units_base::electronVolt);
  energies[3] = convert(energies[3], Stamp_Units::energy, SI_Units_base::electronVolt);
  energies[4] = energies[0] + energies[1] + energies[2] + energies[3];

  double pressure = reducePressure() / (3.*Global::domainInfo.getVolume());

  pressure = convert(pressure, Stamp_Units::pressure, SI_Units_base::pascal);
  
  double gt=0;

  double dtLocal(dt), dtGlobal(0.);

  commManager.reduce(dtLocal, dtGlobal, MPI_SUM);

  dtGlobal/=commManager.getNumberOfNodes();

  if(isMaster) gt = dtGlobal / totalNumberOfParticles;

  // First call : add tab entry, open file
  if (firstCall) {

    firstCall = false;

    Array<int> nLocalCores(1, setup.numberOfThreads);
    Array<int> nTotalCores(1, 0);
    commManager.reduce(nLocalCores, nTotalCores, MPI_SUM);

    if (isMaster) {

      numberOfCores = nTotalCores[0];

      std::string filename = "xstamp-energy.log";
      out.open(filename.c_str(), std::fstream::out | std::fstream::trunc);
      out << "# running on " << numNodes << " node(s) -- " << numberOfCores << " total cores " << std::endl;

      timeZero[0] = energies[0];
      timeZero[1] = energies[1];
      timeZero[2] = energies[2];
      timeZero[3] = energies[3];
      timeZero[4] = energies[4];
      timeZero[5] = temperature;
      timeZero[6] = totalNumberOfParticles;

      writeHeader_amr(flux);
      writeHeader_amr(out);

    }

  }

  uint64_t nRecvAtom(0), nSendAtom(0); 

  // Get total memory usage (kB)
  
  double memUsage=0;
  memUsage = metrics.memStats(commManager);

  // Get Imbalance rate
  double imbalance = 0;
  if (step>-1) 
    imbalance = imbalanceRate(true, step);
  else
  {

  uint64_t nRecvAtomLocal(0);
  uint64_t nSendAtomLocal(0);

  aferrgreg(nSendAtomLocal,nRecvAtomLocal);

  MPI_Reduce( &nRecvAtomLocal     , &nRecvAtom    , 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce( &nSendAtomLocal     , &nSendAtom    , 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

  }
  


  // Write all -- special version for final step
  if (isMaster) {

    if (step > -1) {
      write_line_amr(flux, step, totalNumberOfParticles, state, imbalance, energies, temperature, gt, numberOfCores, memUsage);
      write_line_amr(out, step, totalNumberOfParticles, state, imbalance,  energies, temperature, gt, numberOfCores, memUsage);
    }

    else {



      // convert total time in days/hours/min/sec
      uint64_t days[5], hrs[5], min[5];
      double sec[5];
      double tot = 0;
      for (uint i = 0; i < 5; ++i) {
	double tmp = 0;
	if (i < 4) {
	  tmp = metrics.total((Metrics::Time) i);
	  tot += tmp;
	} else {
	  tmp = tot;
	}
	days[i] = auxFloor<uint64_t>(tmp / 86400.);
	tmp -= days[i] * 86400.;
	hrs[i] = auxFloor<uint64_t>(tmp / 3600.);
	tmp -= hrs[i] * 3600.;
	min[i] = auxFloor<uint64_t>(tmp / 60.);
	tmp -= min[i] * 60.;
	sec[i] = tmp;
      }

      // write all on log and std::cout
      write_end(flux, totalNumberOfParticles, energies, temperature, timeZero,
		memUsage, dt, gt, days, hrs, min, sec);
      write_end(out, totalNumberOfParticles, energies, temperature, timeZero,
		memUsage, dt, gt, days, hrs, min, sec);



  	std::cout << "# The overall number of atoms transfered during the simulation is : " << nRecvAtom<< " (recv) and " << nRecvAtom << " (send)"<< std::endl;



    }
  }

}


/// @brief Print particles in the specified output type
/// @param [in] step Current step
/// @param [in] cpuTime CPU time (not used)
/// @param [in] physTime Physical time (not used)
/// @param [in] _rep Output directory (not used)
/// @param [in] type Output type (DAT, VTK or HERCULE_PARTICULE_DEP)
void Node::writeParticles(int step, double cpuTime, double physTime,
			  const std::string& _rep, OutputType type) {

  bool doEints = Global::reference.isDPD() || Global::reference.isSDPD();
  bool doProg = Global::reference.isReactive();

  switch (type) {

  case DAT: {
    ParticleOutput* buff = gatherParticles(step, true, true, doEints, doProg);

    uint tmp_rank = commManager.getRank();
    std::string tmp_name = "xstamp-particles_" + std::to_string(tmp_rank);
    ParticleWriterDAT dat(buff, true, true, doEints, doProg, tmp_name);
    //if (isMaster)
      dat.write(step);
    delete buff;
    break;
  }
  case VTK: {
    if (commManager.getRank() == 0){
      std::cout << "IO:Output with VTK ... " << step << std::flush;
    }
    ParticleOutput* buff = gatherParticles(step, true, true, doEints, doProg);
    ParticleWriterVTK vtk(buff, true, true, doEints, doProg, "xstamp-particles");
    if (isMaster)
      vtk.write(step);
    delete buff;
    if (commManager.getRank() == 0){
    	std::cout << std::endl;
    }
    break;
  }
  case HERCULE_PARTICULE_DEP: {
#ifdef __use_lib_hercule
    if (commManager.getRank() == 0){
    	std::cout << "IO:Output with Hercule ... " << step << std::flush;
    }
    ParticleOutput* buff = NULL;
    uint numDom = getNumberOfDomains();
    DomainInterface** domains = getDomains();
    if (numDom == 1) {
      buff = new ParticleOutput(domains[0]->getNumberOfParticles(), true,
				true, doEints, doProg);
      domains[0]->fillBuffer(buff);
    } else {
    	std::cerr << "IO:Output with Hercule . No support several domains."
	   << std::endl;
      break;
    }

    if (m_inputOutput->m_herculeDepW == NULL) {
      m_inputOutput->initHerculeDepW();
    }
    m_inputOutput->m_herculeDepW->write(step, buff);
    delete buff;
    if (commManager.getRank() == 0){
    	std::cout << std::endl;
    }
#else
    std::cerr << "Can't use Hercule output." << std::endl;
#endif 
    break;
  }
  default: // Oops
  	std::cerr << "No output format selected" << std::endl;
    break;

  }
}


/// @brief Print cells in the specified output type
/// @param [in] step Current step
/// @param [in] cpuTime CPU time (not used)
/// @param [in] physTime Physical time (not used)
/// @param [in] _rep Output directory (not used)
/// @param [in] type Output type (not used)
void Node::writeCells(int step, double cpuTime, double physTime,
			  const std::string& _rep, OutputType type) {

  if (commManager.getRank() == 0){
  	std::cout << "IO:Output with VTK ... " << step << std::flush;
  }
  CellOutput* buff = gatherCells(step);
  CellWriterVTK vtk(buff, "xstamp-cells");
  if (isMaster)
    vtk.write(step);
  delete buff;
  if (commManager.getRank() == 0){
  	std::cout << std::endl;
  }

}


/// @brief Write particles in dump output in the specified output type
/// @param [in] step Current step
/// @param [in] cpuTime CPU time
/// @param [in] physTime Physical time
/// @param [in] _rep Output directory (not used)
/// @param [in] type Dump output type (LEGACY_DUMP or HERCULE_DUMP)
void Node::dumpParticles(int step, double cpuTime, double physTime, const std::string& _rep, DumpType type) {

  switch (type) {

  case LEGACY_DUMP: {
    if (commManager.getRank() == 0){
    	std::cout << "IO:Legacy with MPIO ... " << step << std::flush;
    }
    if (Global::reference.isDPD() || Global::reference.isSDPD()) {
      Array<LegacyHeaderIOStruct> header;
      Array<LegacyDPDEParticleIOStruct> particlesArray;
      setLegacyDumpData(header, particlesArray, step, cpuTime, physTime);
      ParticleWriterLegacyDump stampv3Dump(&commManager,inputOutput());
      stampv3Dump.write(step, header, particlesArray);
    }
    else {
      Array<LegacyHeaderIOStruct> header;
      Array<LegacyParticleIOStruct> particlesArray;
      setLegacyDumpData(header, particlesArray, step, cpuTime, physTime);
      ParticleWriterLegacyDump stampv3Dump(&commManager,inputOutput());
      stampv3Dump.write(step, header, particlesArray);
    }
    if (commManager.getRank() == 0){
    	std::cout << std::endl;
    }
    break;
  }
  case HERCULE_DUMP: {
#ifdef __use_lib_hercule
    if (commManager.getRank() == 0){
    	std::cout << "IO:Dump with Hercule ... " << step << std::flush;
    }

    if (Global::reference.isDPD() || Global::reference.isSDPD()) {
      Array<LegacyHeaderIOStruct> header;
      HerculeDPDEParticleIODumpStruct particles;
      setLegacyDumpData(header, particles, step, cpuTime, physTime);
      
      if (m_inputOutput->m_herculeDumpW == NULL) {
	m_inputOutput->initHerculeDumpW();
      }
      m_inputOutput->m_herculeDumpW->write(header, particles, step, cpuTime, physTime);
      if (commManager.getRank() == 0){
	std::cout << std::endl;
      }
    }
    else {
      Array<LegacyHeaderIOStruct> header;
      HerculeParticleIODumpStruct particles;
      setLegacyDumpData(header, particles, step, cpuTime, physTime);
      
      if (m_inputOutput->m_herculeDumpW == NULL) {
	m_inputOutput->initHerculeDumpW();
      }
      m_inputOutput->m_herculeDumpW->write(header, particles, step, cpuTime, physTime);
      if (commManager.getRank() == 0){
	std::cout << std::endl;
      }
    }
#else
    std::cerr << "Can't use Hercule dump." << std::endl;
#endif

    break;
  }
  default: // Oops
  	std::cerr << "No output format selected" << std::endl;
    break;

  }

}


/// @brief Write dump output in the specified output type and flux
/// @param [in,out] flux Print flux
/// @param [in] step Current step
/// @param [in] dt CPU time per iteration (not used)
/// @param [in] physTime Physical time
/// @param [in] _rep Output directory
/// @param [in] type Dump output type (LEGACY_DUMP, HERCULE_DUMP or MOL_DUMP)
// with step = time->CurrentStep() and dt =  cpuTime/iterations, physTime = physical time
void Node::writeDump(std::ostream& flux, int step, double dt, double physTime,
		     const std::string& _rep, DumpType type) {
  dumpParticles(step, metrics.total(Metrics::COMPUTE), physTime, _rep, type);
}


/// @brief Write node general data in the specified flux as a check
///
/// Written data is :
/// number of nodes ans threads
/// I/O specifications
/// global geometry settings and domain decomposition
/// atom types from atomReference and potentialReference
/// @param [in,out] flux Print flux
// Function to print node content (debug)
void Node::print(std::ostream& flux) {

  using namespace std;

  if (isMaster) {
    flux << "GLOBAL " << endl
#if __use_orchestrator
         << "  Git version             : " << VERSION << endl
#endif
#if __restrict_1_thread_per_core
         << "  Hwloc restriction       : 1 thread per core" << endl
#elif __restrict_2_threads_per_core
         << "  Hwloc restriction       : 2 threads per core" << endl
#else
         << "  Hwloc restriction       : no" << endl
#endif
         << endl;
  }

  if (isMaster) {
    flux << "NODE " << endl << "  " << setw(22) << left << "Node(s)"
	 << " : " << setw(6) << numberOfNodes() << endl << "  "
	 << setw(22) << left << "Thread(s) per node" << " : " << setw(6)
	 << setup.numberOfThreads << endl;
  }

  writeEndLine(flux);

  if (isMaster) {
    flux << "IO " << endl
	 << "  Input type             : " << this->inputOutput()->initType()   << endl 
	 << "  Input dir              : " << this->inputOutput()->initDir()    << endl 
	 << "  Input step             : " << this->inputOutput()->initStep()   << endl 
	 << "  Output type            : " << this->inputOutput()->outputType() << endl 
	 << "  Output rate            : " << this->time->outputRate()          << endl
	 << "  Output dir             : " << this->inputOutput()->outputDir()  << endl 
	 << "  Dump type              : " << this->inputOutput()->dumpType()   << endl 
	 << "  Dump dir               : " << this->inputOutput()->dumpDir()    << endl 
	 << "  Dump rate              : " << this->time->dumpRate()            << endl
         << "  Number of steps        : " << this->time->numberOfSteps() << endl
	 << endl;
    Global::domainInfo.print(flux);
    this->printLocal(flux);
  }

  writeEndLine(flux);


  if(isMaster) {

    omp_sched_t kind;
    int chunk_size;
    omp_get_schedule(&kind,&chunk_size); 
 
    flux << "AMR Option"<<endl;
  
    if(kind==omp_sched_static)
      flux << "  Scheduling omp         : static , "<<chunk_size << endl;
    if(kind==omp_sched_dynamic)
      flux << "  Scheduling omp         : dynamic , "<<chunk_size << endl;
    if(kind==omp_sched_guided)
      flux << "  Scheduling omp         : guided , "<<chunk_size << endl;
    if(kind==omp_sched_auto)
      flux << "  Scheduling omp         : auto , "<<chunk_size << endl;

    flux << "  Dmax                   : " << Global::reference.Dmax << endl  
         << "  Criterion              : " << Global::reference.r_amrCriterion <<endl  ;

    if(!Global::reference.isBlockVerlet())
      flux << "  Neighbor method        : CellList + Verlet Lists" << endl ;
    else
      flux << "  Neighbor method        : CellList + Verlet Lists + Block Verlet" << endl ;

  }

  writeEndLine(flux);

  // Species data
  if (isMaster) {
    flux << "SPECIES " << endl;
    Global::reference.print(flux);
  }



  writeEndLine(flux, 2);

  flux << right;

}


/// @brief Print domains local info (like symmetrization) in specified flux
/// @param [in,out] flux Print flux
void NodeSingleDomain::printLocal(std::ostream& flux) {

  domain->printInfoBase(flux);

}


/// @brief Gather particles on all nodes
/// @param [in] step Step where the gathering is done (not used)
/// @param [in] doTypes Indicates if the types are written
/// @param [in] doVelocities Indicates if the velocities are written
/// @param [in] doEints Indicates if the internal energies are written
/// @param [in] doProg Indicate if the progress variables are written
/// @return A ParticleOutput containing all particles
ParticleOutput* NodeSingleDomain::gatherParticles(int step, bool doTypes, bool doVelocities,
						  bool doEints, bool doProg) {

  // Get the number of particles and nodes
  uint64_t totalNumberOfParticles = reduceNumberOfParticles();
  uint numNodes = (uint) numberOfNodes();

  // Set the recipent buffer on the master node
  ParticleOutput* bufferOut = nullptr;
  Array<int> counts;
  Array<int> disps;

  if (isMaster) {
    bufferOut = new ParticleOutput((uint) totalNumberOfParticles, doTypes,
				   doVelocities, doEints, doProg);
    counts = Array<int>(numNodes, 0);
    disps = Array<int>(numNodes, 0);
  } else {
    bufferOut = new ParticleOutput(1, doTypes, doVelocities, doEints, doProg);
  }

  // Fill local buffer (bufferIn)
  ParticleOutput* bufferIn = new ParticleOutput(domain->getNumberOfParticles(), doTypes, doVelocities, doEints,doProg);
  domain->fillBuffer(bufferIn);

  // Set the number of object to receive from each node and the position where they must be put
  Array<int> sendCounts(1, bufferIn->id.size());
  commManager.gather(sendCounts, counts);

  if (isMaster) {
    for (uint i = 1; i < numNodes; ++i)
      disps[i] = disps[i - 1] + counts[i - 1];
  }

  // Big gather
  bufferOut->gather(&commManager, bufferIn, counts, disps);
  delete bufferIn;

  return bufferOut;
  
  //return bufferIn;

}


/// @brief Gather cells on all nodes
/// @param [in] step Step where the gathering is done (not used)
/// @return A CellOutput containing all cells
CellOutput* NodeSingleDomain::gatherCells(int step) {

  // Get the number of cells and nodes
  uint64_t totalNumberOfCells = reduceNumberOfCells();
  uint numNodes = (uint) numberOfNodes();

  // Set the recipent buffer on the master node
  CellOutput* bufferOut = nullptr;
  Array<int> counts;
  Array<int> disps;

  if (isMaster) {
    bufferOut = new CellOutput((uint) totalNumberOfCells);
    counts = Array<int>(numNodes, 0);
    disps = Array<int>(numNodes, 0);
  } else {
    bufferOut = new CellOutput(1);
  }

  // Fill local buffer (bufferIn)
  CellOutput* bufferIn = new CellOutput(domain->getNumberOfRealCells());
  domain->fillBuffer(bufferIn);

  // Set the number of object to receive from each node and the position where they must be put
  Array<int> sendCounts(1, bufferIn->r.size());
  commManager.gather(sendCounts, counts);

  if (isMaster) {
    for (uint i = 1; i < numNodes; ++i)
      disps[i] = disps[i - 1] + counts[i - 1];
  }

  // Big gather
  bufferOut->gather(&commManager, bufferIn, counts, disps);
  delete bufferIn;

  return bufferOut;

}


/// @brief Create data for legacy dump writing (default version)
/// @param [out] header Created header
/// @param [out] particlesArray Created particles array
/// @param [in] step Step where the dump writing is done
/// @param [in] cpuTime CPU time
/// @param [in] physTime Physical time
void NodeSingleDomain::setLegacyDumpData(Array<LegacyHeaderIOStruct>& header,
					 Array<LegacyParticleIOStruct>& particlesArray, int step, double cpuTime,
					 double physTime) {

  uint64_t totalNumberOfParticles = reduceNumberOfParticles();

  const double c = setup.expandFactor - 1.0;

  //Filling Header
  header = Array<LegacyHeaderIOStruct>(1);
  header[0].iterationNumber = step;
  header[0].time = convert(physTime, Stamp_Units::time,
			   SI_Units_base::second);
  header[0].CPUtime = cpuTime;

  const vec3<double>& coordMin = Global::domainInfo.getMinBounds();
  const vec3<double>& coordMax = Global::domainInfo.getMaxBounds();
  const vec3<double>& coordExt = Global::domainInfo.getExtension();
  header[0].xmin = convert(coordMin.x, Stamp_Units::length,
			   SI_Units_base::meter);
  header[0].ymin = convert(coordMin.y, Stamp_Units::length,
			   SI_Units_base::meter);
  header[0].zmin = convert(coordMin.z, Stamp_Units::length,
			   SI_Units_base::meter);
  header[0].xmax = convert(coordMax.x, Stamp_Units::length,
			   SI_Units_base::meter);
  header[0].ymax = convert(coordMax.y, Stamp_Units::length,
			   SI_Units_base::meter);
  header[0].zmax = convert(coordMax.z, Stamp_Units::length,
			   SI_Units_base::meter);

  header[0].domainNumber = Global::domainInfo.getNumberOfDomains();
  header[0].particlesTotalNumber = (int) totalNumberOfParticles;

  const Array<double> energies = reduceEnergy();
  header[0].kineticEnergy = convert(energies[0], Stamp_Units::energy,
				    SI_Units_base::joule);
  header[0].potentialEnergy = convert(energies[1], Stamp_Units::energy,
				      SI_Units_base::joule);
  header[0].internalEnergy = convert(energies[2], Stamp_Units::energy,
				      SI_Units_base::joule);
  header[0].totalEnergy = convert(energies[0] + energies[1] + energies[2],
				  Stamp_Units::energy, SI_Units_base::joule);
  header[0].rotationalEnergy = 0.0;

  for (uint i = 0; i < 5; ++i)
    header[0].johnDoe[i] = 0.0;

  //
  if (!Global::domainInfo.noFreeBoundaryConditions()) {

    Array<int> localOutOfFreeBounds(6, 0);
    Array<int> globalOutOfFreeBounds(6, 0);
    domain->checkFreeBoundaries(localOutOfFreeBounds);
    commManager.allReduce(localOutOfFreeBounds, globalOutOfFreeBounds, MPI_SUM);

    if (globalOutOfFreeBounds[0] > 0) header[0].xmin -= convert(c * coordExt[0], Stamp_Units::length, SI_Units_base::meter);
    if (globalOutOfFreeBounds[1] > 0) header[0].ymin -= convert(c * coordExt[1], Stamp_Units::length, SI_Units_base::meter);
    if (globalOutOfFreeBounds[2] > 0) header[0].zmin -= convert(c * coordExt[2], Stamp_Units::length, SI_Units_base::meter);
    if (globalOutOfFreeBounds[3] > 0) header[0].xmax += convert(c * coordExt[0], Stamp_Units::length, SI_Units_base::meter);
    if (globalOutOfFreeBounds[4] > 0) header[0].ymax += convert(c * coordExt[1], Stamp_Units::length, SI_Units_base::meter);
    if (globalOutOfFreeBounds[5] > 0) header[0].zmax += convert(c * coordExt[2], Stamp_Units::length, SI_Units_base::meter);

  }

  //tout le monde : nb paticules pour mon domain, settaille, puis je remplie
  particlesArray = Array<LegacyParticleIOStruct>(domain->getNumberOfParticles());

  domain->fillBuffer(particlesArray);

}


/// @brief Create data for legacy dump writing (default version) for DPDE
/// @param [out] header Created header
/// @param [out] particlesArray Created particles array
/// @param [in] step Step where the dump writing is done
/// @param [in] cpuTime CPU time
/// @param [in] physTime Physical time
void NodeSingleDomain::setLegacyDumpData(Array<LegacyHeaderIOStruct>& header,
					 Array<LegacyDPDEParticleIOStruct>& particlesArray, int step, double cpuTime,
					 double physTime) {

  uint64_t totalNumberOfParticles = reduceNumberOfParticles();

  const double c = setup.expandFactor - 1.0;

  //Filling Header
  header = Array<LegacyHeaderIOStruct>(1);
  header[0].iterationNumber = step;
  header[0].time = convert(physTime, Stamp_Units::time,
			   SI_Units_base::second);
  header[0].CPUtime = cpuTime;

  const vec3<double>& coordMin = Global::domainInfo.getMinBounds();
  const vec3<double>& coordMax = Global::domainInfo.getMaxBounds();
  const vec3<double>& coordExt = Global::domainInfo.getExtension();
  header[0].xmin = convert(coordMin.x, Stamp_Units::length,
			   SI_Units_base::meter);
  header[0].ymin = convert(coordMin.y, Stamp_Units::length,
			   SI_Units_base::meter);
  header[0].zmin = convert(coordMin.z, Stamp_Units::length,
			   SI_Units_base::meter);
  header[0].xmax = convert(coordMax.x, Stamp_Units::length,
			   SI_Units_base::meter);
  header[0].ymax = convert(coordMax.y, Stamp_Units::length,
			   SI_Units_base::meter);
  header[0].zmax = convert(coordMax.z, Stamp_Units::length,
			   SI_Units_base::meter);

  header[0].domainNumber = Global::domainInfo.getNumberOfDomains();
  header[0].particlesTotalNumber = (int) totalNumberOfParticles;

  const Array<double> energies = reduceEnergy();
  header[0].kineticEnergy = convert(energies[0], Stamp_Units::energy,
				    SI_Units_base::joule);
  header[0].potentialEnergy = convert(energies[1], Stamp_Units::energy,
				      SI_Units_base::joule);
  header[0].internalEnergy = convert(energies[2], Stamp_Units::energy,
				      SI_Units_base::joule);
  header[0].totalEnergy = convert(energies[0] + energies[1] + energies[2],
				  Stamp_Units::energy, SI_Units_base::joule);
  header[0].rotationalEnergy = 0.0;

  for (uint i = 0; i < 5; ++i)
    header[0].johnDoe[i] = 0.0;

  //
  if (!Global::domainInfo.noFreeBoundaryConditions()) {

    Array<int> localOutOfFreeBounds(6, 0);
    Array<int> globalOutOfFreeBounds(6, 0);
    domain->checkFreeBoundaries(localOutOfFreeBounds);
    commManager.allReduce(localOutOfFreeBounds, globalOutOfFreeBounds, MPI_SUM);

    if (globalOutOfFreeBounds[0] > 0) header[0].xmin -= convert(c * coordExt[0], Stamp_Units::length, SI_Units_base::meter);
    if (globalOutOfFreeBounds[1] > 0) header[0].ymin -= convert(c * coordExt[1], Stamp_Units::length, SI_Units_base::meter);
    if (globalOutOfFreeBounds[2] > 0) header[0].zmin -= convert(c * coordExt[2], Stamp_Units::length, SI_Units_base::meter);
    if (globalOutOfFreeBounds[3] > 0) header[0].xmax += convert(c * coordExt[0], Stamp_Units::length, SI_Units_base::meter);
    if (globalOutOfFreeBounds[4] > 0) header[0].ymax += convert(c * coordExt[1], Stamp_Units::length, SI_Units_base::meter);
    if (globalOutOfFreeBounds[5] > 0) header[0].zmax += convert(c * coordExt[2], Stamp_Units::length, SI_Units_base::meter);

  }

  //tout le monde : nb particules pour mon domain, settaille, puis je remplis
  particlesArray = Array<LegacyDPDEParticleIOStruct>(domain->getNumberOfParticles());

  domain->fillBuffer(particlesArray);

}


/// @brief Create data for legacy dump writing (Hercule version)
/// @param [out] header Created header
/// @param [out] particles Created particles array
/// @param [in] step Step where the dump writing is done
/// @param [in] cpuTime CPU time
/// @param [in] physTime Physical time
void NodeSingleDomain::setLegacyDumpData(Array<LegacyHeaderIOStruct>& header,
					 HerculeParticleIODumpStruct& particles, int step, double cpuTime,
					 double physTime) {

  uint64_t totalNumberOfParticles = reduceNumberOfParticles();

  const double c = setup.expandFactor - 1.0;

  //Filling Header
  header = Array<LegacyHeaderIOStruct>(1);
  header[0].iterationNumber = step;
  header[0].time = convert(physTime, Stamp_Units::time,
			   SI_Units_base::second);
  header[0].CPUtime = cpuTime;

  const vec3<double>& coordMin = Global::domainInfo.getMinBounds();
  const vec3<double>& coordMax = Global::domainInfo.getMaxBounds();
  const vec3<double>& coordExt = Global::domainInfo.getExtension();
  header[0].xmin = coordMin.x;
  header[0].ymin = coordMin.y;
  header[0].zmin = coordMin.z;
  header[0].xmax = coordMax.x;
  header[0].ymax = coordMax.y;
  header[0].zmax = coordMax.z;
  header[0].domainNumber = Global::domainInfo.getNumberOfDomains();
  header[0].particlesTotalNumber = (int) totalNumberOfParticles;

  const Array<double> energies = reduceEnergy();
  header[0].kineticEnergy = energies[0];
  header[0].potentialEnergy = energies[1];
  header[0].internalEnergy = energies[2];
  header[0].totalEnergy = energies[0] + energies[1] + energies[2];
  header[0].rotationalEnergy = 0.0;

  for (uint i = 0; i < 5; ++i)
    header[0].johnDoe[i] = 0.0;

  //
  if (!Global::domainInfo.noFreeBoundaryConditions()) {

    Array<int> localOutOfFreeBounds(6, 0);
    Array<int> globalOutOfFreeBounds(6, 0);
    domain->checkFreeBoundaries(localOutOfFreeBounds);
    commManager.allReduce(localOutOfFreeBounds, globalOutOfFreeBounds,
			  MPI_SUM);

    if (globalOutOfFreeBounds[0] > 0) header[0].xmin -= c * coordExt[0];
    if (globalOutOfFreeBounds[1] > 0) header[0].ymin -= c * coordExt[1];
    if (globalOutOfFreeBounds[2] > 0) header[0].zmin -= c * coordExt[2];
    if (globalOutOfFreeBounds[3] > 0) header[0].xmax += c * coordExt[0];
    if (globalOutOfFreeBounds[4] > 0) header[0].ymax += c * coordExt[1];
    if (globalOutOfFreeBounds[5] > 0) header[0].zmax += c * coordExt[2];

  }

#define HPROT_CONVERT_UNIT
  // pb dans input
#ifdef HPROT_CONVERT_UNIT
  header[0].xmin = convert(header[0].xmin, Stamp_Units::length, SI_Units_base::meter);
  header[0].ymin = convert(header[0].ymin, Stamp_Units::length, SI_Units_base::meter);
  header[0].zmin = convert(header[0].zmin, Stamp_Units::length, SI_Units_base::meter);
  header[0].xmax = convert(header[0].xmax, Stamp_Units::length, SI_Units_base::meter);
  header[0].ymax = convert(header[0].ymax, Stamp_Units::length, SI_Units_base::meter);
  header[0].zmax = convert(header[0].zmax, Stamp_Units::length, SI_Units_base::meter);
  header[0].kineticEnergy = convert(header[0].kineticEnergy, Stamp_Units::energy, SI_Units_base::joule);
  header[0].potentialEnergy = convert(header[0].potentialEnergy, Stamp_Units::energy,   SI_Units_base::joule);
  header[0].totalEnergy = convert(header[0].totalEnergy, Stamp_Units::energy, SI_Units_base::joule);
#endif

  //tout le monde : nb paticules pour mon domain, settaille, puis je remplie
  //int nbCell=0;
  //vec3<int> &nbCellPerDim = Global::domainInfo.getNumberOfCellsPerDim();

  particles.init(0, domain->getNumberOfParticles());

  domain->fillBuffer(particles);

}


/// @brief Create data for legacy dump writing (Hercule version for DPDE)
/// @param [out] header Created header
/// @param [out] particles Created particles array
/// @param [in] step Step where the dump writing is done
/// @param [in] cpuTime CPU time
/// @param [in] physTime Physical time
void NodeSingleDomain::setLegacyDumpData(Array<LegacyHeaderIOStruct>& header,
					 HerculeDPDEParticleIODumpStruct& particles, int step, double cpuTime,
					 double physTime) {

  uint64_t totalNumberOfParticles = reduceNumberOfParticles();

  const double c = setup.expandFactor - 1.0;

    //Filling Header
  header = Array<LegacyHeaderIOStruct>(1);
  header[0].iterationNumber = step;
  header[0].time = convert(physTime, Stamp_Units::time,
			   SI_Units_base::second);
  header[0].CPUtime = cpuTime;

  const vec3<double>& coordMin = Global::domainInfo.getMinBounds();
  const vec3<double>& coordMax = Global::domainInfo.getMaxBounds();
  const vec3<double>& coordExt = Global::domainInfo.getExtension();
  header[0].xmin = coordMin.x;
  header[0].ymin = coordMin.y;
  header[0].zmin = coordMin.z;
  header[0].xmax = coordMax.x;
  header[0].ymax = coordMax.y;
  header[0].zmax = coordMax.z;
  header[0].domainNumber = Global::domainInfo.getNumberOfDomains();
  header[0].particlesTotalNumber = (int) totalNumberOfParticles;

  const Array<double> energies = reduceEnergy();
  header[0].kineticEnergy = energies[0];
  header[0].potentialEnergy = energies[1];
  header[0].internalEnergy = energies[2];
  header[0].totalEnergy = energies[0] + energies[1] + energies[2];
  header[0].rotationalEnergy = 0.0;

  for (uint i = 0; i < 5; ++i)
    header[0].johnDoe[i] = 0.0;

  //
  if (!Global::domainInfo.noFreeBoundaryConditions()) {

    Array<int> localOutOfFreeBounds(6, 0);
    Array<int> globalOutOfFreeBounds(6, 0);
    domain->checkFreeBoundaries(localOutOfFreeBounds);
    commManager.allReduce(localOutOfFreeBounds, globalOutOfFreeBounds,
			  MPI_SUM);

    if (globalOutOfFreeBounds[0] > 0) header[0].xmin -= c * coordExt[0];
    if (globalOutOfFreeBounds[1] > 0) header[0].ymin -= c * coordExt[1];
    if (globalOutOfFreeBounds[2] > 0) header[0].zmin -= c * coordExt[2];
    if (globalOutOfFreeBounds[3] > 0) header[0].xmax += c * coordExt[0];
    if (globalOutOfFreeBounds[4] > 0) header[0].ymax += c * coordExt[1];
    if (globalOutOfFreeBounds[5] > 0) header[0].zmax += c * coordExt[2];

  }

#define HPROT_CONVERT_UNIT
  // pb dans input
#ifdef HPROT_CONVERT_UNIT
  header[0].xmin = convert(header[0].xmin, Stamp_Units::length, SI_Units_base::meter);
  header[0].ymin = convert(header[0].ymin, Stamp_Units::length, SI_Units_base::meter);
  header[0].zmin = convert(header[0].zmin, Stamp_Units::length, SI_Units_base::meter);
  header[0].xmax = convert(header[0].xmax, Stamp_Units::length, SI_Units_base::meter);
  header[0].ymax = convert(header[0].ymax, Stamp_Units::length, SI_Units_base::meter);
  header[0].zmax = convert(header[0].zmax, Stamp_Units::length, SI_Units_base::meter);
  header[0].kineticEnergy = convert(header[0].kineticEnergy, Stamp_Units::energy, SI_Units_base::joule);
  header[0].potentialEnergy = convert(header[0].potentialEnergy, Stamp_Units::energy,   SI_Units_base::joule);
  header[0].totalEnergy = convert(header[0].totalEnergy, Stamp_Units::energy, SI_Units_base::joule);
#endif

  //tout le monde : nb paticules pour mon domain, settaille, puis je remplie
  //int nbCell=0;
  //vec3<int> &nbCellPerDim = Global::domainInfo.getNumberOfCellsPerDim();

  particles.init(0, domain->getNumberOfParticles());

  domain->fillBuffer(particles);

}





