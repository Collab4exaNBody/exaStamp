#pragma once

class Octree : public leafCell {

	typedef Octree self_type;

	public :

	vec3<int> offset;
	bool m_isGhost;

	std::vector<uint64_t> id; ///< IDs of the particles
	std::vector<uint8_t> ti; ///< Types of the particles
	std::vector<uint32_t> morton; ///< Types of the particles  

	std::vector<double> ep; ///< Potential energies of the particles

	std::vector<double> rx; ///< X components of the positions of the particles
	std::vector<double> ry; ///< Y components of the positions of the particles
	std::vector<double> rz; ///< Z components of the positions of the particles

	std::vector<double> vx; ///< X components of the velocities of the particles
	std::vector<double> vy; ///< Y components of the velocities of the particles
	std::vector<double> vz; ///< Z components of the velocities of the particles

	std::vector<double> fx; ///< X components of the forces on the particles
	std::vector<double> fy; ///< Y components of the forces on the particles
	std::vector<double> fz; ///< Z components of the forces on the particles

	Octree() : 
	 leafCell(this->size),
		offset(0,0,0), m_isGhost(false),
		id(0), ti(0), morton(0), ep(0), 
		rx(0), ry(0), rz(0), 
		vx(0), vy(0), vz(0), 
		fx(0), fy(0), fz(0) 
	{ 
		this->granny = this; 
		this->size = 0; 
		this->shift=0;
	}
	       
	      
	Octree(const Octree &a) {}

	~Octree() {}

	//accessor
	inline bool getIsGhost();
	inline bool getIsReal();
	
	virtual inline uint64_t* getId() override;
	virtual inline uint8_t* getType() override;
	virtual inline double* getPotentialEnergy() override;
	virtual inline double* getForceX() override;
	virtual inline double* getForceY() override;
	virtual inline double* getForceZ() override;
	virtual inline double* getVelocityX() override;
	virtual inline double* getVelocityY() override;
	virtual inline double* getVelocityZ() override;
	virtual inline double* getPositionX() override;
	virtual inline double* getPositionY() override;
	virtual inline double* getPositionZ() override;  

	virtual int id_capacity() const override;
	virtual inline double getPotentialEnergy(const int i) const override;
	virtual inline uint64_t  getId(const int i) const override;
	virtual inline uint8_t  getType(const int i) const override;
	virtual inline double getForceX(const int i) const override;
	virtual inline double getForceY(const int i) const override;
	virtual inline double getForceZ(const int i) const override;
	virtual inline double getVelocityX(const int i) const override;
	virtual inline double getVelocityY(const int i) const override;
	virtual inline double getVelocityZ(const int i) const override;
	virtual inline double getPositionX(const int i) const override;
	virtual inline double getPositionY(const int i) const override;
	virtual inline double getPositionZ(const int i) const override;

	inline void setOffset(vec3<int>& o, bool ghost);
	inline vec3<int> getOffset();
	inline void updateGranny();
	
	//communication

	inline  void ghostVelCopy(self_type& to) const ;
	inline  void ghostVelCopy(const int destDomain, const vec3<int>& destCell, std::vector< std::tuple<int, ExchangeGhostV3> >& leaving) const;


  void fillVerlet();
  bool checkVerlet(double rVerlet);

	void computeForceEmb(EAMPotential* pot, const uint8_t typeIndexA, const uint8_t typeIndexB);

	inline void resize(const uint& size);
	inline void resizeInfos(const uint& size);
	inline void resizePositions(const uint& size);
	inline void resizeVelocities(const uint& size);
	inline void resizeForces(const uint& size);

	inline  void remove(const uint& index);  
	inline  void remove_last();
	inline  void clear();

	void add(const MPI__ParticleBase& q);
	void add(const MPI__Particle& q);
	void add(const MPI__Mesoparticle& q);

	//refine
	inline void adjust(const size_t refinement_max);
	template<typename F> inline void refine(const size_t refinement_max, F criterion);
	virtual inline void createChild(const size_t refinement_max, size_t *nAtomPerCell) override;
	inline void fill_indx_morton( const vec3<double>& inv_sizeOfCell, const vec3<double> position);


	// sort 
	inline void sort(const size_t refinement_max , const vec3<double>& sizeOfCell, const vec3<double> &inf);

	void __debug_print(std::ostream& flux) const;
	std::tuple<uint, uint, double, double> __debug_diag() const;

	inline void ghostAdjustAMR(const Correcter& correcter);

	void checkEAMData() {
		this->m_eamStorage.check(this->size);
	}

	inline double computePotentialEnergyAMR() const;
	inline void setEmbAMR(const ExchangeEAM& exchange);  

	int8_t computeNeighborMov(const uint i, const vec3<double>& minBounds, const vec3<double>& maxBounds) const;
  void internalReorganization(std::list< std::tuple<int, MPI__Particle> >& list, const Array<int>& nbrs, const vec3<double>& minBounds, const vec3<double>& maxBounds);
	void getExchange(MPI__Particle& p, const uint i) const;

  // io buffer
  void fillBuffer(ParticleOutput* buffer, const uint start);
  template <class DumpStruct> void fillBuffer(DumpStruct& buffer, const uint start);
  void fillBuffer(ParticleInSitu* buffer, const uint start);
  void fillGhostBuffer(ParticleInSitu* buffer, const uint start);
  void fillBuffer(CellOutput* buffer, const uint start);
  
  // balance
  double computeMemory() const;
  inline double computeWorkloadAMR() const;
};

