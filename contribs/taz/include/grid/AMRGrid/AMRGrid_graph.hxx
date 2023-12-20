#pragma once

#include <functional>

#if _OPENMP >= 201511

	#define DEPS_EVAL(X) tmp_depGraph[pTask[X]]

	#define DEPS1 DEPS_EVAL(0)
	#define DEPS2 DEPS1, DEPS_EVAL(1)
	#define DEPS3 DEPS2, DEPS_EVAL(2) 
	#define DEPS4 DEPS3, DEPS_EVAL(3) 
	#define DEPS5 DEPS4, DEPS_EVAL(4) 
	#define DEPS6 DEPS5, DEPS_EVAL(5) 
	#define DEPS7 DEPS6, DEPS_EVAL(6) 
	#define DEPS8 DEPS7, DEPS_EVAL(7) 
	#define DEPS9 DEPS8, DEPS_EVAL(8) 
	#define DEPS10 DEPS9, DEPS_EVAL(9) 
	#define DEPS11 DEPS10, DEPS_EVAL(10) 
	#define DEPS12 DEPS11, DEPS_EVAL(11) 
	#define DEPS13 DEPS12, DEPS_EVAL(12) 
	#define DEPS14 DEPS13, DEPS_EVAL(13) 
	#define DEPS15 DEPS14, DEPS_EVAL(14) 
	#define DEPS16 DEPS15, DEPS_EVAL(15) 
	#define DEPS17 DEPS16, DEPS_EVAL(16) 
	#define DEPS18 DEPS17, DEPS_EVAL(17) 
	#define DEPS19 DEPS18, DEPS_EVAL(18)  
	#define DEPS20 DEPS19, DEPS_EVAL(19) 
	#define DEPS21 DEPS20, DEPS_EVAL(20) 
	#define DEPS22 DEPS21, DEPS_EVAL(21) 
	#define DEPS23 DEPS22, DEPS_EVAL(22) 
	#define DEPS24 DEPS23, DEPS_EVAL(23) 
	#define DEPS25 DEPS24, DEPS_EVAL(24) 
	#define DEPS26 DEPS25, DEPS_EVAL(25) 


	#define EVAL_N(X) DEPS ## X
	#define EVAL_PRAGMA(x) _Pragma (#x)
	#define DO_PRAGMA(x) EVAL_PRAGMA(x)
	#define OMP_TASK_DEPEND(N,OUT) DO_PRAGMA(omp task depend(in: EVAL_N(N)) depend(out:tmp_depGraph[OUT]))

	#define CASE(N, ...) case N:\
		{\
			OMP_TASK_DEPEND(N,uTask)\
			potential(idOctree);\
		} break;

#endif

TMPLSG template<typename F>
inline void TMPL_AMRGrid::OpenMPTask(std::size_t nbDep, std::size_t uTask, int idOctree, std::size_t * pTask, F potential) {

	// nbDep = number of dependencies for octree idOctree
	assert(nbDep <=26);
	// version openmp 4.5
#if _OPENMP >= 201511

	int * tmp_depGraph = depGraph;

	switch(nbDep)
	{
		case 0: 
		{
			#pragma omp task depend(out: tmp_depGraph[uTask])
			potential(idOctree);
		} break;
		CASE(1);
		CASE(2);
		CASE(3);
		CASE(4);
		CASE(5);
		CASE(6);
		CASE(7);
		CASE(8);
		CASE(9);
		CASE(10);
		CASE(11);
		CASE(12);
		CASE(13);
		CASE(14);
		CASE(15);
		CASE(16);
		CASE(17);
		CASE(18);
		CASE(19);
		CASE(20);
		CASE(21);
		CASE(22);
		CASE(23);
		CASE(24);
		CASE(25);
		CASE(26);
		default: std::cout << " problème dans le graphe de dépendance " <<std::endl;
	}
#endif
 }

TMPLSG template<typename F> void TMPL_AMRGrid::launchGraph(F potential)
{
    std::cout << "launchGraph, nwaves="<<cWave.size()<<std::endl;
	#pragma omp parallel
	#pragma omp master  /* the master thread launches tasks*/
	for(size_t i = 0; i < cWave.size() ; i++)  /* cWave is pre-sorted according to the wave number  */
		OpenMPTask(cWave.numberOfDependencies(i), cWave.valueOfDependency(i), cWave.getIndex(i), cWave.arrayOfDependencies(i), potential);
}


TMPLSG template<typename F> void TMPL_AMRGrid::launchGraphEdge(F potential)
{
	#pragma omp parallel
	#pragma omp master  /* the master thread launches tasks*/
	for(size_t i = 0; i < cWaveEdge.size() ; i++)  /* cWave is pre-sorted according to the wave number  */
		OpenMPTask(cWaveEdge.numberOfDependencies(i), cWaveEdge.valueOfDependency(i), cWaveEdge.getIndex(i), cWaveEdge.arrayOfDependencies(i), potential);
}

TMPLSG void TMPL_AMRGrid::BuildWaveMethod()
{
	// OpenMp
	depGraph = new int[ (numberOfCellsPerDim.x+4) * (numberOfCellsPerDim.y+4) * (numberOfCellsPerDim.z+4)];
}
