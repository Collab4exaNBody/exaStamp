TMPLSG inline void TMPL_AMRGrid::refineCells() {

	assert(doRefine >=0);

	if(doRefine==0)
	{
		/* Apply refinement criteria to determine if a cell have to be refined */
		Refine();
		doRefine = refineStepMax; 
	}  
	else
	{
		/* Atoms are moved, they have to be re sorted                               */
		/* then pointers are adjusted on the storage included into their root cell. */
		Adjust();
		std::cout << " adjust " << std::endl;
		doRefine--;
	}
}


template <class cell> inline void dLeafCells(cell* C)
{
  for(size_t i=0; i<8; i++)
    if(C->getDaughterCell(i)!=nullptr)
      dLeafCells(C->getDaughterCell(i));
  delete C;  
}

TMPLSG inline void TMPL_AMRGrid::destructorLeafCells() {

	const auto& cells = this->getTraversal(TraversalManager::REAL);

	#pragma omp parallel for schedule(runtime)
	for(size_t i=0; i<cells.size(); ++i)
	{ 
		size_t idx = cells[i];

    if(particles[idx].getDaughterCell(0) != nullptr) // never refined ?
      for(size_t child=0; child<8; child++)
      {
        dLeafCells(particles[idx].getDaughterCell(child));
        particles[idx].resetPtrDaughterCell(child);
      }

		particles[idx]. setIsLeaf(1);
	}


}


TMPLSG inline void TMPL_AMRGrid::Refine() {

	const int crit = g_amrCriterionValue;
	auto criterionFunction = [crit](Cell *C)->bool {return C->size > crit ;};


	UpdateTraversalAMRoctrees();

	#pragma omp parallel for schedule(runtime)
	for(int i=0; i<octrees.size(); ++i)
	{ 
		const vec3<double> cellLength = Global::domainInfo.getCellLength();
		const vec3<double> minBounds  = Global::domainInfo.getMinBounds();

		// redefine root cell as being a leaf cell
		octrees[i]->setIsLeaf(1);

		if(octrees[i]->getNumberOfParticles() != 0)
		{
			octrees[i]->sort(levelMax, cellLength, minBounds);
			octrees[i]->refine(levelMax,criterionFunction);
			assert(octrees[i]->keepRightNumberOfAtoms());
		}
	}

	/* Update the traversals containing octree and leaf cell information */
	UpdateTraversalAMR();

	if (Global::reference.isMEAM())
		UpdateTraversalAMR_MEAM();

	/* For each leaf cell, the list of its neighbouring leaf cells is updated. */
	fillneighborCells_List();
}

TMPLSG inline void TMPL_AMRGrid::Adjust() {

	static const vec3<double>& cellLength = Global::domainInfo.getCellLength();
	static const vec3<double>& minBounds  = Global::domainInfo.getMinBounds();


	UpdateTraversalAMRoctrees();

	#pragma omp parallel for schedule(runtime)
	for(int i=0; i<octrees.size(); ++i)
	{
		if(octrees[i]->getNumberOfParticles() != 0) 
		{
			octrees[i]->sort(levelMax, cellLength, minBounds);
			octrees[i]->adjust(levelMax);
			assert(octrees[i]->keepRightNumberOfAtoms());

		}
	}
}


template<class leaf>
inline void findLeafCell(leaf* C, std::vector<leaf*> & toFill)
{
	if(C->getIsLeaf() == 1)
	{
	  if(C->size > 0)
		  toFill.push_back(C) ;
  	}
	else 
		for(uint8_t i = 0 ; i<8 ; ++i)
			findLeafCell(C->getDaughterCell(i), toFill) ;
}

TMPLSG void TMPL_AMRGrid::UpdateTraversalAMR(){

	leafcells.clear();

	const auto& cells = this->getTraversal(TraversalManager::REAL);

	for (int i = 0; i < cells.size(); ++i)
	{
	  if(particles[cells[i]].size>0)
	  {
		if(particles[cells[i]].getIsLeaf() == 1) 
		{
		  if(particles[cells[i]].size >0)
			  leafcells.push_back( static_cast<Cell*>(&particles[cells[i]]));
	  	}
		else
			for(uint8_t j = 0 ; j<=7 ; ++j)
			  findLeafCell(particles[cells[i]].getDaughterCell(j), leafcells);
	  }
	} 
}



TMPLSG void TMPL_AMRGrid::UpdateTraversalAMRoctrees(){

	octrees.clear();

	const auto& cells = this->getTraversal(TraversalManager::REAL);

	for (int i = 0; i < cells.size(); ++i)
	{
	  if(particles[cells[i]].size>0)
	  {
		  octrees.push_back(&particles[cells[i]]);
	  }
	} 
}

TMPLSG inline void TMPL_AMRGrid::UpdateTraversalAMR_MEAM(){

	octreesMEAM.clear();
	leafcellsMEAM.clear();

	const auto& cells = this->getTraversal(TraversalManager::ALL);

	for (int i = 0; i < cells.size(); ++i)
	{
	  assert(cells[i]>=0);
	  assert(cells[i]<particles.size());
	  
	  if(particles[cells[i]].getNumberOfParticles() > 0 || particles[cells[i]].getIsGhost() )
	  {
		octreesMEAM.push_back(cells[i]);

	  
		if(particles[cells[i]].getIsLeaf() == 1) 
		{
			  leafcellsMEAM.push_back( static_cast<Cell*>(&particles[cells[i]]));
		}
		else
			for(uint8_t j = 0 ; j<8 ; ++j)
			  findLeafCell(particles[cells[i]].getDaughterCell(j), leafcellsMEAM);
	  }
	} 
}
