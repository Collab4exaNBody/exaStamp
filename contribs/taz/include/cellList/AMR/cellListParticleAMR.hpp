#ifndef __CELL_LIST_AMR_HPP
#define __CELL_LIST_AMR_HPP


// Je copie tout les includes de cellListParticles
#include <algorithm>
#include <list>
#include <ostream>
#include <tuple>
#include <vector>

#include "cellList/cellListBase.hpp"

#include "io/particleOutput.hpp"
#include "io/cellOutputVTK.hpp"

#include "parallel/types/MPI_exchangeEAM.hpp"
#include "parallel/types/MPI_exchangeGhostV3.hpp"

#include "potential/allPotentials.hpp"
//
#include "simd/neighbors.hpp"
#include "simd/push.hpp"
#include "simdAMR/push.hpp"


#include "utils/neighbor.hpp"
#include "utils/stampUnits.hpp"
#include "io/particleInSitu.hpp"

typedef std::tuple<uint, uint16_t, uint8_t> nbr_id;


/// @brief Macros for enable if
///
/// Usage : @n
/// <tt>/<class Pot_t, __ENABLE_IF(Pot_t::has_simd_opt)/> function</tt>
/// will be used if the potential Pot_t has an SIMD optimized version,
/// otherwise <tt>function</tt> will be used
#define __ENABLE_IF(_bool_) class = typename std::enable_if<(_bool_)>::type


class Octree;

// @brief A class that handle the particles in a cell
/// @tparam elem_chunk Base chunk for the arrays
class leafCell : public CellListBase {

	typedef leafCell self_type;

	public :

	//info
	uint level;
	vec3<int> position;
	int isleaf;

	self_type *daughter[8];
	self_type* mother;
	Octree* granny;

	std::vector<self_type *>  neighborCells;
	uint size;
	uint shift;
  mat3<double> pressureTensor; ///< Local pressure tensor 

	public :

	/* constructors */
	leafCell() :  size(0), CellListBase(), 
		level(0), position(0,0,0),
		isleaf(1), mother(nullptr), 
		granny(nullptr), shift(0),
		neighborCells(0), pressureTensor() 
	{
		for(uint8_t i=0;i<8;i++) 
		  daughter[i]=nullptr;
	}


	leafCell(const leafCell &a) {}

	/* reference parameters is used to fill numberOfParticles in the CellList */
	leafCell(const uint& i) : CellListBase(),
		level(0), position(0,0,0),
		isleaf(1), mother(nullptr), 
		granny(nullptr), size(0), shift(0),
		neighborCells(0) , pressureTensor() 
	{
		for(uint8_t i=0;i<8;i++) 
			daughter[i]=nullptr;
	}

	leafCell(const Octree& mortherCell) {};


	leafCell(uint begin, uint numberOfElements, leafCell *motherCell, uint8_t posChild) :
	shift(begin), size(numberOfElements), 
	CellListBase(),
	neighborCells(0) 
	{
		// define infos;
		isleaf   = true;
		level    = motherCell->level +1; 

		/* poschild =(x={0,1}y={0,1}z={0,1})_bit */
		/* pos_daughter =  pos_mother + vec3_convert posChild */ 
		position = motherCell->position*2 + vec3<int>( (posChild>>2) & 1 , (posChild>>1) & 1 , posChild & 1 );

		/* define family information */
		granny     = motherCell->granny;
		mother     = motherCell;

		for(uint8_t i=0;i<8;i++) 
		  daughter[i]=nullptr;

	}
	       
	~leafCell() {}

	template<class Cell> inline void addNeighborCell (Cell * cell);
	inline void clearNeighborCells();
	inline int getLevel ();
	inline uint getSize ();
	inline vec3<int> getPosition ();
	inline vec3<int>& getRefPosition ();
	inline self_type** getNeighborCells ();
	inline uint getNeighborCellsSize ();
	inline self_type* getMotherCell ();
	inline self_type* getDaughterCell (size_t i);
  	inline void resetPtrDaughterCell (size_t i);
	inline int getIsLeaf ();
	inline void setIsLeaf (int i);
	inline void avoidDoublon(vec3<int> & numberOfCellsPerDim, int TotalNumberOfCells);
	inline void setInfo(vec3<int> pos, int l);
	inline vec3<int> getOffset(); 
	inline bool getIsGhost(); 
	inline bool getIsEdge(); 
  	size_t getNumberOfParticles();

	virtual inline uint64_t* getId();
	virtual inline uint8_t* getType();
	const inline uint32_t* getMorton() const;
	virtual inline double* getPotentialEnergy();
	virtual inline double* getForceX();
	virtual inline double* getForceY();
	virtual inline double* getForceZ();
	virtual inline double* getVelocityX();
	virtual inline double* getVelocityY();
	virtual inline double* getVelocityZ();
	virtual inline double* getPositionX();
	virtual inline double* getPositionY();
	virtual inline double* getPositionZ();  

	virtual int id_capacity() const { return 0; } ; // warning
 
	virtual inline double getPotentialEnergy(const int i) const;
	virtual inline uint64_t  getId(const int i) const;
	virtual inline uint8_t  getType(const int i) const;
	virtual inline double getForceX(const int i) const;
	virtual inline double getForceY(const int i) const;
	virtual inline double getForceZ(const int i) const;
	virtual inline double getVelocityX(const int i) const;
	virtual inline double getVelocityY(const int i) const;
	virtual inline double getVelocityZ(const int i) const;
	virtual inline double getPositionX(const int i) const;
	virtual inline double getPositionY(const int i) const;
	virtual inline double getPositionZ(const int i) const;

	inline vec3<double> getRealPosition(vec3<double>& sizeOfCell);
	inline void checkLeafLocals();
	template <bool fill_emb> uint       fillForceBufferVerlet(const uint i, const uint8_t typeIndex);  

	uint fillForceBufferVerletMEAM( const uint i, const uint8_t typeIndex); 

	uint resizeForceBuffer(uint size); 
	     
	// Mutator

	inline void incForceX(const uint i, const double value);
	inline void incForceY(const uint i, const double value);
	inline void incForceZ(const uint i, const double value);
	inline void incPontentialEnergy(const uint i, const double value);

	// impl
	void adjust(const size_t refinement_max, size_t *nAtomPerCell) ;
	virtual void createChild( const size_t refinement_max, size_t *nAtomPerCell);
	
	void defineNewCell( size_t begin, size_t numberOfElements, self_type* motherCell, size_t posChild ) ;
	void updateOldCell( size_t begin, size_t numberOfElements) ;

	template<typename F> void refine(const size_t refinement_max,  size_t* nAtomPerCell, F criterion);

	// compute forces
	#define __TMPL_CellList_computeForce_baseArgs const uint8_t typeIndexA, const uint8_t typeIndexB, CList* cells, const uint8_t* ghostLayer, const Correcter& correcter, const bool symmetrize

	template <class CList> void computeForcePairVerlet(PairPotential* pot, __TMPL_CellList_computeForce_baseArgs);
	template <class CList> void computeForceRhoVerlet  (EAMPotential* pot, __TMPL_CellList_computeForce_baseArgs);
	template <class CList> void computeForceFinalVerlet(EAMPotential* pot, __TMPL_CellList_computeForce_baseArgs);


	template <class Pot_t, class CList, __ENABLE_IF(Pot_t::has_simd_opt)> void computeForcePairVerlet (Pot_t* pot, __TMPL_CellList_computeForce_baseArgs);
	template <class Pot_t, class CList, __ENABLE_IF(Pot_t::has_simd_opt)> void computeForcePairVerlet_mutex (Pot_t* pot, __TMPL_CellList_computeForce_baseArgs); // mutex veriosn
	template <class Pot_t, class CList, __ENABLE_IF(Pot_t::has_simd_opt)> void computeForceRhoVerlet  (Pot_t* pot, __TMPL_CellList_computeForce_baseArgs);
	template <class Pot_t, class CList, __ENABLE_IF(Pot_t::has_simd_opt)> void computeForceFinalVerlet(Pot_t* pot, __TMPL_CellList_computeForce_baseArgs);
	template <bool mutex, class Pot_t, class CList> void computeForcesAndPotentialMEAMVerlet(Pot_t* pot,const uint8_t typeIndexA, const uint8_t typeIndexB, int cell, CList* cells, const uint8_t* ghostLayer, const Correcter& correcter);
	template <class Pot_t> void computeForceRhoVerlet_2  (Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, const bool symmetrize);
	template <class Pot_t> void computeForceRhoVerlet_without_buffer  (Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, const bool symmetrize);
	template <class Pot_t> void computeForceRhoVerlet_without_buffer_per_block  (Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB);
	template <class Pot_t> void computeForceFinalVerlet_2(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, const bool symmetrize);
	template <class Pot_t> void computeForceFinalVerlet_without_buffer(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, const bool symmetrize);
	template <class Pot_t> void computeForceFinalVerlet_without_buffer_per_block(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB);
	template <class Pot_t> void computeForcePairVerlet (Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB);
	template <class Pot_t> void computeForceEmb (Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB);
	template <class Pot_t> void computeForcePair_perso (Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB);       // without mutex
	template <class Pot_t> void computeForcePair_buffer (Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB);   // without mutex
	template <class Pot_t> void computeForcePairVerlet_mutex (Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB); // mutex version
	
	template <class Pot_t> void computeForcePair_perso_per_block (Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB);
	template <class Pot_t> void computeForcePair_perso_per_block_mutex (Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB);

	void makeNeighborsAMR(const size_t cellIndex, self_type& cell, const bool* mask, const bool symmetrize, const bool perBlock);
	void makeNeighborsAMR_self(const size_t myIndex, const bool symmetrize, const bool perBlock);

	template <bool mutex, class Pot_t> void computeForcesAndPotentialMEAMVerlet_AMR(Pot_t* pot,const uint8_t typeIndexA, const uint8_t typeIndexB, int cell, const uint8_t* ghostLayer);

	#undef __TMPL_CellList_computeForce_baseArgs

	inline void getProjAMR(const vec3<int>& coords, const vec3<int>& dcell, bool mask[]);
	void writeForce( const uint i, const uint8_t typeIndex, const bool symmetrize, uint nbrSize);
	void writeRho(const uint i, const uint8_t typeIndex, const bool symmetrize, uint nbrSize);

	template<bool mutex> void writeForceMEAM(const uint cell, const uint i, const uint8_t typeIndex, const uint8_t* ghostLayer, const uint size);

	/// @brief Shortcut for the arguments of a MEAM force computation
	#define __TMPL_CellList_Write_baseArgs Pot_t* pot, const uint cell, const uint i, const uint8_t typeIndex, const uint8_t* ghostLayer,const uint size
	/// @brief Shortcut for the arguments of a MEAM force computation (rho0)
	#define __Write_rho0_MEAM double rho0,vec3<double> &rhod0
	/// @brief Shortcut for the arguments of a MEAM force computation (rho1)
	#define __Write_rho1_MEAM double rho1,vec3<double> &rhod1
	/// @brief Shortcut for the arguments of a MEAM force computation (rho2)
	#define __Write_rho2_MEAM double rho2,vec3<double> &rhod2
	/// @brief Shortcut for the arguments of a MEAM force computation (rho3)
	#define __Write_rho3_MEAM double rho3,vec3<double> &rhod3
	
	template <bool mutex, class Pot_t> void writeRho_0123_MEAM(__TMPL_CellList_Write_baseArgs,__Write_rho0_MEAM,__Write_rho1_MEAM,__Write_rho2_MEAM,__Write_rho3_MEAM);

	#undef __TMPL_CellList_Write_baseArgs
	#undef __Write_rho0_MEAM
	#undef __Write_rho1_MEAM
	#undef __Write_rho2_MEAM
	#undef __Write_rho3_MEAM

	/* Manage particles */
	/* The storage is included in the root cell or Octree */ 
	template <bool resize_velocity, bool resize_fe, class CList> inline void ghostCopyAMR(CList& to, const Correcter& correcter);

	template <class CList> void embCopyAMR(CList& to)  const;
	template <class Q> void ghostCopyAMR(const int destDomain, const vec3<int>& destCell, std::vector< std::tuple<int, Q> >& leaving) const;
	inline double& embAMR(size_t i);
	inline double& rhoAMR(size_t i);
	void m_eamStorageReset();
	void embCopyAMR(const int destDomain, const vec3<int>& destCell, std::vector< std::tuple<int, ExchangeEAM> >& leaving) const ;

	inline void checkMtx();

	/* assert functions */

	template <class CList> bool disjoint(CList& to);
	bool atomDoublon();
	bool keepRightNumberOfAtoms();
	
	template <class DumpStruct> 
	void fillBufferAMR(Array<DumpStruct>& buffer, const uint start);
	
	void neighborListAddNeighbor(const uint16_t index, const nbr_id& nbr);
  	void neighborListGetNeighbors(const uint16_t index, const uint8_t type, nbr_id*& start, uint& size);
	leafCell** neighborListGetLeafCellNeighbors(const uint16_t index);
	void sortNeighbors();
	

};

#include "cellList/AMR/CellListParticleAMR_Base.hpp"

#include "cellList/AMR/leafCell/cellListParticleAMR.hxx"
#include "cellList/AMR/leafCell/cellListParticleAMR_accessors.hxx"
#include "cellList/AMR/leafCell/cellListParticleAMR_refine.hxx"
#include "cellList/AMR/leafCell/cellListParticleAMR_compute.hxx"
#include "cellList/AMR/leafCell/cellListParticleAMR_neighbors.hxx"
#include "cellList/AMR/leafCell/cellListParticleAMR_force.hxx"
#include "cellList/AMR/leafCell/cellListParticleAMR_force_lj.hxx"
#include "cellList/AMR/leafCell/cellListParticleAMR_force_eam.hxx"
#include "cellList/AMR/leafCell/cellListParticleAMR_force_meam.hxx"
#include "cellList/AMR/leafCell/cellListParticleAMR_force_buffer.hxx"
#include "cellList/AMR/leafCell/cellListParticleAMR_communication.hxx"
#include "cellList/AMR/leafCell/cellListParticleAMR_assert.hxx"


#include "cellList/AMR/rootCell/cellListParticleAMR_Base.hxx"
#include "cellList/AMR/rootCell/cellListParticleAMR_Base_neighbor.hxx"
#include "cellList/AMR/rootCell/cellListParticleAMR_Base_sort.hxx"
#include "cellList/AMR/rootCell/cellListParticleAMR_Base_refine.hxx"
#include "cellList/AMR/rootCell/cellListParticleAMR_Base_communication.hxx"
#include "cellList/AMR/rootCell/cellListParticleAMR_Base_accessor.hxx"
#include "cellList/AMR/rootCell/cellListParticleAMR_Base_force.hxx"
#include "cellList/AMR/rootCell/cellListParticleAMR_Base_buffer.hxx"
#include "cellList/AMR/rootCell/cellListParticleAMR_Base_balance.hxx"

#undef leafCell 

#endif

