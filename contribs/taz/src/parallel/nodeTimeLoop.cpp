/// @file
/// @brief Implementation of node time related methods


#include <chrono>
#include <fstream>
#include <thread>

#include "domain/domainInterface.hpp"

#include "parallel/node.hpp"

#include "time/time.hpp"
#include "time/timeIntegration.hpp"




#ifdef __use_lib_ccc_user
extern "C"{
#include "ccc_user.h"
}
#endif


extern bool updateVerletLists; ///< Need for initialisation

/*!
  \def gettid()
  get syscall
*/
#define gettid() syscall(SYS_gettid)


/// @brief Update cells on the domain
/// @warning Used only one time in the whole program despite multiple tests
/// @param [in,out] upToDate Indicates if cells update is needed
inline void NodeSingleDomain::updateCells(bool& upToDate) {
  if (!upToDate) {

    domain->updateCells();
    upToDate = true;
    updateVerletLists=true;
  }
}


/// @brief Do one step of the calculation
///
///
void NodeSingleDomain::oneStep() {


    metrics.tic(Metrics::COMPUTE);

    // integration of motion equations
    scheme->oneStep(domain, time->delta());

    metrics.toc(Metrics::COMPUTE);
    metrics.store(Metrics::COMPUTE);

}


/// @brief Load balancing at step 0
///
/// Initialize load balancer
/// @param [out] state Load balancer state
/// @param [out] upToDate Indicates if cells update is needed
void NodeSingleDomain::first_balance(LBS::State& state, bool& upToDate) {

  metrics.tic(Metrics::BALANCE);

  // balance cells, get state, update query functions
  if (time->balanceRate()!=-1) {

    loadBalancer.setMethod(loadBalancer.method(), LBS::PARTITION);

    state = domain->balance(&loadBalancer);
    loadBalancer.setCallbackQueryFunctions(this->getNumberOfDomains(), this->getDomains());
    upToDate = true;

  }

  metrics.toc(Metrics::BALANCE);
  metrics.store(Metrics::BALANCE);
  
}


/// @brief Load balancing
/// @param [out] state Load balancer state
/// @param [out] upToDate Indicates if cells update is needed
void NodeSingleDomain::balance(LBS::State& state, bool& upToDate) {

  metrics.tic(Metrics::BALANCE);

  // balance cells, get state, update querry functions    
  if ((time->balance()) && time->currentStep()!=time->numberOfSteps()) {

    domain->updateCellsInside();

    state = domain->balance(&loadBalancer);
    loadBalancer.setCallbackQueryFunctions(this->getNumberOfDomains(), this->getDomains());
    upToDate = true;
  }

  metrics.toc(Metrics::BALANCE);
  metrics.store(Metrics::BALANCE);
  
}

/// @brief Function to end the iteration
///
/// This function ends the iteration and especially communicates with
/// the orchestrator in case there are data analytics to perform.
/// It replaces and generalizes the writeIO function.
/// Logs and dumps are still performed by the simulation master threads but outputs are handled thanks to the orchestrator.
/// @todo For the moment, only the particles are copied into the shared structure and the cells are handled by the simulation
/// @param [in,out] state Load balancer state
/// @param [in,out] upToDate Indicates if cells update is needed
/// @param [in,out] tmpTime Compute time since last log writing
void NodeSingleDomain::endIteration(LBS::State& state, bool& upToDate, double& tmpTime)
{



  metrics.tic(Metrics::IO); // measure time of IO

  // Update cells, write log
  if (time->log()) {
    updateCells(upToDate);
    writeEnergies(std::cout, time->currentStep(), tmpTime/time->logRate(), state);
    state   = LBS::VOID;
    tmpTime = 0.;
  }
  
  // check simulation computed values against a reference DB, if enabled
  verifySimulationValues(timeManager()->currentStep());
  
  // Update cells, write dump
  if (time->dump()) {

    domain->updateCellsInside();
    domain->initForces(false,0.);
    writeDump(std::cout, time->currentStep(), tmpTime/time->logRate(), time->elapsed(), m_inputOutput->dumpDir(), m_inputOutput->dumpType());
  }
  

  metrics.toc(Metrics::IO);
  metrics.store(Metrics::IO);
  
}



/// @brief Check remaining time and particle leak at free boundaries
/// @param [in,out] upToDate Indicates if cells update is needed
/// @param [in] tmpTime Time since last log writing
void NodeSingleDomain::goAheadChecking(bool& upToDate, double& tmpTime) {

#ifdef __use_lib_ccc_user

  if (time->log()) {

    // at least 600 seconds are needed to properly stop the calculation
    constexpr double minTimeNeededToStop = 600.0;       

    Array<int> abort(1, 0) ;
    double remainingTime = 99999.99;

    if (isMaster) {                           

      ccc_tremain(&remainingTime);
      std::ifstream test("stop");

      abort[0] = (remainingTime<minTimeNeededToStop || test.good()) ? 1:0 ;

      test.close();

    }
    
    commManager.broadcast(abort);

    if (abort[0]) {

      time->abort();

      if (time->dumpRate() != -1) {

        domain->updateCellsInside();

	writeDump(std::cout, time->currentStep(), tmpTime/time->logRate(), time->elapsed(), m_inputOutput->dumpDir(), m_inputOutput->dumpType());
      }
      
    }

  }
    
#endif

  if (!setup.expandFreeLimits) return;

  Array<int> outOfFreeBounds_local (6, 0);
  domain->checkFreeBoundaries(outOfFreeBounds_local);
  commManager.allReduce(outOfFreeBounds_local, m_outOfFreeBounds, MPI_SUM);

  bool abort = false;
  for (uint8_t i=0; i<6; ++i) {
    if (m_outOfFreeBounds[i]>0) {
      abort = true;
      break;
    }
  }

  if (abort) {

    time->abort();

    if (isMaster) {
      std::cout<< std::endl
	       << "[WARNING] at iteration " << time->currentStep() << " : particle(s) out of 'free' boundary " << "("
	       <<   "xmin " << m_outOfFreeBounds[0] << ", xmax " << m_outOfFreeBounds[3]
	       << ", ymin " << m_outOfFreeBounds[1] << ", ymax " << m_outOfFreeBounds[4]
	       << ", zmin " << m_outOfFreeBounds[2] << ", zmax " << m_outOfFreeBounds[5]
	       << "). " << std::endl << "          " 
	       << "ExaSTAMP will stop and will produce a dump if the input file has been set properly. "
	       << std::endl << std::endl;
    }

    if (time->dumpRate() != -1) {
    
        domain->updateCellsInside();
        writeDump(std::cout, time->currentStep(), tmpTime/time->logRate(), time->elapsed(), m_inputOutput->dumpDir(), m_inputOutput->dumpType());
    }
    
  }

}


/// @brief Check whether walls should be stopped
///
void NodeSingleDomain::wallChecking() {

  if (time->stopWalls()) {
    domain->stopWalls();

    if (isMaster) {
      std::cout<< std::endl
	       << "[WARNING] at iteration " << time->currentStep() << " wall movement has been stopped " 
	       << std::endl << std::endl;
    }
  }
  
}

/// @brief Check whether Verlet list should be updated
///  @param [in,out] collect Indicates if the verlet list will be updated
void NodeSingleDomain::verletChecking(bool& collect)  {

  if(collect==false) {
    collect = domain->checkVerlet();
    int tmp = collect;
    int tmp2;
    MPI_Allreduce( &tmp, &tmp2, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    collect = tmp2;
  }
}


void NodeSingleDomain::collectTimer()  {

  double max_global_init;
  double local_init = metrics.total(Metrics::INIT);

  double max_global_compute;
  double local_compute = metrics.total(Metrics::COMPUTE);

  double max_global_io;
  double local_io = metrics.total(Metrics::IO);

  double max_global_balance;
  double local_balance = metrics.total(Metrics::BALANCE);

  double max_global_timeloop;
  double local_timeloop = metrics.total(Metrics::TIMELOOP);

  MPI_Reduce( &local_init     , &max_global_init    , 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce( &local_compute  , &max_global_compute , 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce( &local_io       , &max_global_io      , 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce( &local_balance  , &max_global_balance , 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce( &local_timeloop , &max_global_timeloop, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


  if(isMaster) {

  std::cout << std::fixed << std::setw(5) << std::setprecision(2) ;
  std::cout << "# Timers "<<std::endl<<std::endl
            << "#" << std::endl
            << "#    init    : " << max_global_init     << " s " << std::endl
            << "#    compute : " << max_global_compute  << " s " << std::endl
            << "#    ioMaster: " << local_io       << " s " << std::endl
            << "#    balance : " << max_global_balance  << " s " << std::endl
            << "#    total   : " << max_global_timeloop << " s " << std::endl<<std::endl;

  }

  if(Global::reference.isTimers())
  {
     metrics.data.write_end(max_global_compute, isMaster);
  }
}

/// @brief Do all the work (except initialisation)
///
///
void NodeSingleDomain::doComputeWork() {
  
  double tmpTime   = 0.;
  
  bool   upToDate ;
  if(Global::reference.isVerlet()) upToDate = true;
  else upToDate = false;

  auto   firstStep = time->currentStep();
  LBS::State state = LBS::VOID;

  
    if(Global::reference.isTimers())
    {
      ptM = getMetrics();
      ptM->data.reserve(time->final()-time->start()+1);
    }


  std::chrono::steady_clock::time_point startTime, endTime;
  
  // Print check
  print(std::cout);


  metrics.toc(Metrics::INIT);
  metrics.store(Metrics::INIT);


  // Balance at step 0
  first_balance(state, upToDate);

  // Write log at step 0
  if (time->logRate() != -1)
    {
    //  metrics.tic(Metrics::IO);
      writeEnergies(std::cout, time->currentStep(), 0., state);
    //  metrics.toc(Metrics::IO);
    //  metrics.store(Metrics::IO);
    }
  state = LBS::VOID;

 
  if (time->outputRate()!= -1){
    metrics.tic(Metrics::IO);
    writeParticles(time->currentStep(), 0., 0., m_inputOutput->outputDir(), m_inputOutput->outputType());
    if(m_inputOutput->outputCells()) writeCells(time->currentStep(), 0., 0., m_inputOutput->outputDir(), m_inputOutput->outputType());
    metrics.toc(Metrics::IO);
    metrics.store(Metrics::IO);
  }
  
  // check batch remaining time 
  goAheadChecking(upToDate, tmpTime);


  double deb, end;
  getManager().barrier();
  deb = MPI_Wtime();

  // // Big time loop
  while( !time->isFinished() ) {
  
    if(Global::reference.isTimers())
    {
      ptM = getMetrics();
      ptM->data.resize_one();
    }

    // increment time 
    ++(*time);

    metrics.tic(Metrics::TIMELOOP);

    // define if the verlet lists must be updated
    if(Global::reference.isVerlet() )
      verletChecking(updateVerletLists);
      

    // one step, save time    
    oneStep();
    tmpTime += metrics.elapsed(Metrics::COMPUTE);

    // wall checking
    wallChecking();
    
    // load balancing
    balance(state, upToDate);

    // end of iteration
    endIteration(state, upToDate, tmpTime);

    // check batch remaining time 
    goAheadChecking(upToDate, tmpTime);

    metrics.toc(Metrics::TIMELOOP);
    metrics.store(Metrics::TIMELOOP);
  }
  
  getManager().barrier();
  end = MPI_Wtime();  

  if(isMaster)
  {
    std::cout << "# Le vrai temps total de la simulation est de "<< end-deb << " secondes. " << std::endl;  
  }
  
  
#if __use_orchestrator
  // wait for last analytics to complete
  metrics.tic(Metrics::SIMULATIONSLEEP);
  while(!this->m_shared->dataRefreshable)
    {
      std::this_thread::sleep_for(std::chrono::microseconds(10));
    }
  metrics.toc(Metrics::SIMULATIONSLEEP);
  metrics.store(Metrics::SIMULATIONSLEEP);
#endif /*__use_orchestrator */
  
#if __use_orchestrator
  std::this_thread::sleep_for(std::chrono::microseconds(20)); // test for last analytics
  this->m_shared->stopOrchestrator = true;
#endif /* __use_orchestrator */

#if __debug_tasks 
  endTime = std::chrono::steady_clock::now();
  double totalTime = std::chrono::duration_cast< std::chrono::microseconds > (endTime - startTime).count()*1.0e-06;
#endif /* __debug_tasks */
  
  // Write final log and timers
  writeEnergies(std::cout, -1, metrics.total(Metrics::COMPUTE)/(time->numberOfSteps()-firstStep), LBS::VOID);
  writeEndLine(std::cout);

#if __debug_tasks
  writeTasksInfo(totalTime);
#endif /* __debug_tasks */

  // AMR
  collectTimer();

  // End AMR
}
